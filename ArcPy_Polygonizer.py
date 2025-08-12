# ArcGIS Pro / ArcPy 3.x
# Polygonizer that RE-WRITES the same outputs each run (truncate + append).
# Builds intersections from the street network, creates plus-shaped intersection polygons
# with legs, segments interconnects by level, trims overlaps, and writes stable outputs.

import arcpy
from collections import defaultdict
import math
import os

# =========================
# INPUTS
# =========================
GDB = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb"
CTN_FC = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb\CTN_Test"

SR = arcpy.SpatialReference(2277)  # StatePlane Texas Central (US Feet)

# Buffer widths by street level (feet)
BUFFER_WIDTH_BY_LEVEL = {1: 40, 2: 40, 3: 50, 4: 75}
# Target interconnect segment length by level (feet)
TARGET_LEN_BY_LEVEL   = {1: 150, 2: 150, 3: 250, 4: 250}

# Intersection logic
GOAL_LEG_LEN = 200.0   # max leg length (ft)
TOUCH_TOL    = 5.0     # ft (snap/search tolerance)
HALF_GAP_EPS = 0.05    # ft (anti-overlap epsilon)

# =========================
# PERSISTENT OUTPUTS (rewritten each run)
# =========================
OUT_INTERSECT_PTS      = os.path.join(GDB, "Polygonizer_BuiltIntersections")   # POINT
OUT_TRUE_INTERSECT_PTS = os.path.join(GDB, "Polygonizer_TrueIntersections")    # POINT
OUT_INTERSECT_POLYS    = os.path.join(GDB, "Polygonizer_IntersectionPolys")    # POLYGON
OUT_INTERCONNECT_POLYS = os.path.join(GDB, "Polygonizer_InterconnectPolys")    # POLYGON
OUT_FINAL_POLYS        = os.path.join(GDB, "Polygonizer_FinalPolys")           # POLYGON

# =========================
# TEMPS (created & deleted each run)
# =========================
TMP_SJ                  = os.path.join(GDB, "_pg_tmp_sj")
TMP_RAW_INTERSECT       = os.path.join(GDB, "_pg_tmp_raw_intersections")       # may be MultiPoint
TMP_TRUE_INTERSECT      = os.path.join(GDB, "_pg_tmp_true_intersections")      # may be MultiPoint
TMP_INTER_TRIM          = os.path.join(GDB, "_pg_tmp_inter_trim")
TMP_FINAL_MERGE         = os.path.join(GDB, "_pg_tmp_final_merge")

# =========================
# ENV
# =========================
arcpy.env.workspace = GDB
arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = SR

# =========================
# HELPERS
# =========================
def round_half_up(x): return int(math.floor(x + 0.5))
def clamp(v, lo, hi): return max(lo, min(hi, v))

def street_level_width(level):
    try: return BUFFER_WIDTH_BY_LEVEL.get(int(level), 40)
    except: return 40

def street_level_target_len(level):
    try: return TARGET_LEN_BY_LEVEL.get(int(level), 150)
    except: return 150

def is_empty(g):
    if g is None: return True
    try: return g.isEmpty
    except: pass
    try:
        if hasattr(g, "area"):   return (g.area == 0) and (getattr(g, "partCount", 1) == 0)
        if hasattr(g, "length"): return (g.length == 0) or (getattr(g, "pointCount", 1) == 0)
        if hasattr(g, "pointCount"): return g.pointCount == 0
    except: return False
    return False

def rm(*paths):
    for p in paths:
        if p and arcpy.Exists(p):
            arcpy.Delete_management(p)

def ensure_fc(path, geom_type, sr, fields=None):
    """Ensure a reusable FC: if exists -> TRUNCATE + add missing fields; else create."""
    if arcpy.Exists(path):
        if arcpy.Describe(path).shapeType.lower() != geom_type.lower():
            arcpy.Delete_management(path)
            arcpy.CreateFeatureclass_management(os.path.dirname(path), os.path.basename(path), geom_type, spatial_reference=sr)
        if fields:
            existing = {f.name for f in arcpy.ListFields(path)}
            for (nm, ftype, prec, scale, length) in fields:
                if nm not in existing:
                    arcpy.AddField_management(path, nm, ftype, prec, scale, length)
        arcpy.management.TruncateTable(path)
    else:
        arcpy.CreateFeatureclass_management(os.path.dirname(path), os.path.basename(path), geom_type, spatial_reference=sr)
        if fields:
            for (nm, ftype, prec, scale, length) in fields:
                arcpy.AddField_management(path, nm, ftype, prec, scale, length)
    return path

def as_point_geometry(g):
    if g is None: return None
    try:
        if g.type.lower() == "point": return g
    except: pass
    try:
        pt = getattr(g, "firstPoint", None)
        if pt: return arcpy.PointGeometry(pt, g.spatialReference)
    except: pass
    try:
        cent = g.centroid
        return arcpy.PointGeometry(cent, g.spatialReference)
    except: return None

def explode_points_to_point_fc(src_fc, dst_fc):
    """Write every point inside Point/MultiPoint features from src_fc into dst_fc (POINT)."""
    icur = arcpy.da.InsertCursor(dst_fc, ["SHAPE@"])
    with arcpy.da.SearchCursor(src_fc, ["SHAPE@"]) as sc:
        for (g,) in sc:
            if is_empty(g): 
                continue
            # PointGeometry straight in
            try:
                if g.type.lower() == "point":
                    icur.insertRow([g])
                    continue
            except:
                pass
            # MultiPoint / anything else: try Array access
            sr = g.spatialReference
            try:
                arr = g.getPart(0)  # arcpy.Array of points (for MultiPoint)
                n = getattr(arr, "count", None)
                if n is not None:
                    for i in range(n):
                        pt = arr.getObject(i)
                        if pt:
                            icur.insertRow([arcpy.PointGeometry(pt, sr)])
                    continue
            except:
                pass
            # Generic iteration fallback
            try:
                for pt in g:
                    icur.insertRow([arcpy.PointGeometry(pt, sr)])
                continue
            except:
                pass
            # Last resort: centroid
            icur.insertRow([arcpy.PointGeometry(g.centroid, sr)])
    del icur

# Cleanup temps at start
rm(TMP_SJ, TMP_RAW_INTERSECT, TMP_TRUE_INTERSECT, TMP_INTER_TRIM, TMP_FINAL_MERGE, "streets_lyr")

# =========================
# 1) READ STREETS
# =========================
streets = []
with arcpy.da.SearchCursor(CTN_FC, ["OID@", "SHAPE@", "STREET_LEVEL", "STREET_NAME"]) as cur:
    for oid, shp, lvl, name in cur:
        if is_empty(shp): continue
        streets.append({"oid": oid, "geom": shp, "level": lvl, "name": name, "length": shp.length})
if not streets:
    raise RuntimeError("No street features found in CTN_Test.")

# Layer for quick selects
arcpy.MakeFeatureLayer_management(CTN_FC, "streets_lyr")

# =========================
# 2) BUILD RAW INTERSECTION POINTS (TEMP) & WRITE TO PERSISTENT OUT_INTERSECT_PTS (explode MultiPoint -> Point)
# =========================
arcpy.Intersect_analysis([CTN_FC], TMP_RAW_INTERSECT, "ALL", f"{TOUCH_TOL} Feet", "POINT")
try: arcpy.DeleteIdentical_management(TMP_RAW_INTERSECT, ["SHAPE"])
except: pass
try: arcpy.management.Integrate([TMP_RAW_INTERSECT, CTN_FC], f"{TOUCH_TOL} Feet")
except: pass

ensure_fc(OUT_INTERSECT_PTS, "POINT", SR)
explode_points_to_point_fc(TMP_RAW_INTERSECT, OUT_INTERSECT_PTS)
rm(TMP_RAW_INTERSECT)

# =========================
# 3) KEEP TRUE INTERSECTIONS (>=2 incident lines) -> TEMP -> PERSISTENT (also explode to POINT)
# =========================
arcpy.SpatialJoin_analysis(
    target_features=OUT_INTERSECT_PTS,
    join_features=CTN_FC,
    out_feature_class=TMP_SJ,
    join_operation="JOIN_ONE_TO_ONE",
    join_type="KEEP_COMMON",
    match_option="INTERSECT",
    search_radius=f"{TOUCH_TOL} Feet"
)

def get_join_count_field(fc):
    names = [f.name for f in arcpy.ListFields(fc)]
    for n in names:
        if n.lower().startswith("join_count"): return n
    for n in names:
        if "join" in n.lower() and "count" in n.lower(): return n
    raise RuntimeError(f"Could not find Join_Count field on {fc}. Fields: {names}")

JOIN_COUNT_FIELD = get_join_count_field(TMP_SJ)

arcpy.FeatureClassToFeatureClass_conversion(
    TMP_SJ, GDB, os.path.basename(TMP_TRUE_INTERSECT),
    where_clause=f"{arcpy.AddFieldDelimiters(TMP_SJ, JOIN_COUNT_FIELD)} >= 2"
)
rm(TMP_SJ)

ensure_fc(OUT_TRUE_INTERSECT_PTS, "POINT", SR)
explode_points_to_point_fc(TMP_TRUE_INTERSECT, OUT_TRUE_INTERSECT_PTS)
rm(TMP_TRUE_INTERSECT)

# Pack true nodes
x_oid_field = arcpy.Describe(OUT_TRUE_INTERSECT_PTS).OIDFieldName
true_nodes = []
with arcpy.da.SearchCursor(OUT_TRUE_INTERSECT_PTS, [x_oid_field, "SHAPE@"]) as cur:
    for oid, g in cur:
        if is_empty(g): continue
        true_nodes.append({"oid": oid, "geom": g})
if not true_nodes:
    raise RuntimeError("No true intersections found (>= 2 incident streets).")

# Map street -> [(int_id, measure, pgeom)]
street_to_intersections = defaultdict(list)
for node in true_nodes:
    pgeom = as_point_geometry(node["geom"])
    if pgeom is None: continue
    pbuf = pgeom.buffer(TOUCH_TOL)
    arcpy.SelectLayerByLocation_management("streets_lyr", "INTERSECT", pbuf, selection_type="NEW_SELECTION")
    with arcpy.da.SearchCursor("streets_lyr", ["OID@", "SHAPE@"]) as sc:
        for soid, sgeom in sc:
            try:
                if sgeom.distanceTo(pgeom) <= TOUCH_TOL:
                    m = sgeom.measureOnLine(pgeom, False)
                    street_to_intersections[soid].append((node["oid"], m, pgeom))
            except:
                try:
                    m_on = sgeom.measureOnLine(pgeom, False)
                    proj_pt = sgeom.positionAlongLine(m_on, False)
                    m = sgeom.measureOnLine(proj_pt, False)
                    street_to_intersections[soid].append((node["oid"], m, pgeom))
                except:
                    continue
arcpy.SelectLayerByAttribute_management("streets_lyr", "CLEAR_SELECTION")

# =========================
# 4) LEG LENGTHS PER STREET/INTERSECTION (NO OVERLAP)
# =========================
computed_leg_len = defaultdict(dict)
street_lookup = {s["oid"]: s for s in streets}

for s in streets:
    soid = s["oid"]
    if soid not in street_to_intersections: continue
    entries = sorted(street_to_intersections[soid], key=lambda x: x[1])
    measures = [m for (_pid, m, _pt) in entries]
    for i, (pid, m, _pt) in enumerate(entries):
        nn_gap = s["length"] if len(measures) == 1 else min(abs(m - measures[j]) for j in range(len(measures)) if j != i)
        leg_len = min(GOAL_LEG_LEN, 0.5 * nn_gap - HALF_GAP_EPS, 0.5 * s["length"])
        computed_leg_len[soid][pid] = max(0.0, leg_len)

def extract_leg(street_geom, meet_point, leg_len):
    m = street_geom.measureOnLine(meet_point, False)
    L = street_geom.length
    if m <= (L - m):
        start_m, end_m = 0.0, clamp(leg_len, 0.0, L)
    else:
        start_m, end_m = clamp(L - leg_len, 0.0, L), L
    return street_geom.segmentAlongLine(start_m, end_m, False)

def end_is_intersection(street_oid, at_start=True):
    if street_oid not in street_to_intersections: return False
    L = street_lookup[street_oid]["length"]
    for (_pid, m, _pt) in street_to_intersections[street_oid]:
        if at_start and m <= TOUCH_TOL: return True
        if (not at_start) and abs(L - m) <= TOUCH_TOL: return True
    return False

# =========================
# 5) INTERSECTION POLYGONS (PERSISTENT; TRUNCATE+INSERT)
# =========================
ensure_fc(OUT_INTERSECT_POLYS, "POLYGON", SR, fields=[
    ("X_ID", "LONG", None, None, None),
    ("LEGS", "SHORT", None, None, None),
    ("NOTE", "TEXT", None, None, 64),
])
ins_int = arcpy.da.InsertCursor(OUT_INTERSECT_POLYS, ["SHAPE@", "X_ID", "LEGS", "NOTE"])

for node in true_nodes:
    center = as_point_geometry(node["geom"])
    if center is None: continue

    incidents = [soid for soid, lst in street_to_intersections.items()
                 if any(pid == node["oid"] for pid, _m, _pt in lst)]
    if not incidents: continue

    # Build leg buffers (flat ends)
    leg_buffers = []
    for soid in incidents:
        srec   = street_lookup[soid]
        width  = street_level_width(srec["level"])
        leglen = computed_leg_len[soid].get(node["oid"], 0.0)
        if leglen <= 0.0: continue
        leg_line = extract_leg(srec["geom"], center, leglen)
        tmp = arcpy.Buffer_analysis(leg_line, os.path.join("in_memory", "_leg"), f"{width} Feet", "FULL", "FLAT", "NONE")
        with arcpy.da.SearchCursor(tmp, ["SHAPE@"]) as bc:
            for (g,) in bc:
                if not is_empty(g): leg_buffers.append(g)
        arcpy.Delete_management(tmp)

    # If no legs (edge case), still put a small center fill
    if not leg_buffers:
        max_w = max(street_level_width(street_lookup[soid]["level"]) for soid in incidents)
        center_fill = center.buffer(max_w)
        ins_int.insertRow([center_fill, node["oid"], 0, "CenterOnly"])
        continue

    # Center fill and geometry‑only union (no temp FCs, avoids transactions)
    max_w = max(street_level_width(street_lookup[soid]["level"]) for soid in incidents)
    merged_geom = leg_buffers[0]
    for g in leg_buffers[1:]:
        merged_geom = merged_geom.union(g)
    merged_geom = merged_geom.union(center.buffer(max_w))

    if not is_empty(merged_geom):
        ins_int.insertRow([merged_geom, node["oid"], len(incidents), "OK"])

del ins_int

# =========================
# 6) INTERCONNECT POLYGONS (PERSISTENT; TRUNCATE+INSERT)
# =========================
ensure_fc(OUT_INTERCONNECT_POLYS, "POLYGON", SR, fields=[
    ("STREET_OID", "LONG", None, None, None),
    ("STREET_NAME", "TEXT", None, None, 128),
    ("STREET_LEVEL", "SHORT", None, None, None),
    ("PART_ID", "LONG", None, None, None),
])
ins_ic = arcpy.da.InsertCursor(OUT_INTERCONNECT_POLYS, ["SHAPE@", "STREET_OID", "STREET_NAME", "STREET_LEVEL", "PART_ID"])

def end_measures_for(street_oid):
    srec = street_lookup[street_oid]; L = srec["length"]
    has_start = has_end = False; start_len = end_len = 0.0
    if street_oid in street_to_intersections:
        near_start = [(pid, m) for (pid, m, _pt) in street_to_intersections[street_oid] if m <= TOUCH_TOL]
        near_end   = [(pid, m) for (pid, m, _pt) in street_to_intersections[street_oid] if abs(L - m) <= TOUCH_TOL]
        if near_start:
            has_start = True; pid = min(near_start, key=lambda x: x[1])[0]
            start_len = computed_leg_len[street_oid].get(pid, 0.0)
        if near_end:
            has_end = True; pid = max(near_end, key=lambda x: x[1])[0]
            end_len = computed_leg_len[street_oid].get(pid, 0.0)
    return has_start, start_len, has_end, end_len

for s in streets:
    soid, geom, lvl, name, L = s["oid"], s["geom"], s["level"], s["name"], s["length"]
    width = street_level_width(lvl); target_len = street_level_target_len(lvl)

    has_start_x, start_leg, has_end_x, end_leg = end_measures_for(soid)
    inter_len = L - (start_leg + end_leg)
    if inter_len <= 0.01: continue

    n_parts = max(1, round_half_up(inter_len / target_len))
    part_len = inter_len / n_parts

    for i in range(n_parts):
        a = start_leg + i * part_len
        b = start_leg + (i + 1) * part_len
        seg = geom.segmentAlongLine(a, b, False)

        tmp = arcpy.Buffer_analysis(seg, os.path.join("in_memory", "_seg"), f"{width} Feet", "FULL", "FLAT", "NONE")
        seg_poly = None
        with arcpy.da.SearchCursor(tmp, ["SHAPE@"]) as tc:
            for (tg,) in tc:
                seg_poly = tg; break
        arcpy.Delete_management(tmp)

        add_caps = []
        if i == 0 and not has_start_x:
            p = geom.positionAlongLine(a, False); add_caps.append(p.buffer(width))
        if i == n_parts - 1 and not has_end_x:
            p = geom.positionAlongLine(b, False); add_caps.append(p.buffer(width))
        # ✅ union caps individually (no lists to Geometry.union)
        for cap in add_caps:
            seg_poly = seg_poly.union(cap)

        if not is_empty(seg_poly):
            ins_ic.insertRow([seg_poly, soid, name, lvl, i + 1])

del ins_ic

# =========================
# 7) TRIM + FINAL MERGE (PERSISTENT; TRUNCATE+APPEND)
# =========================
# Dissolve intersections into single geometry in memory
x_geos = []
with arcpy.da.SearchCursor(OUT_INTERSECT_POLYS, ["SHAPE@"]) as xc:
    for (g,) in xc:
        if not is_empty(g): x_geos.append(g)
x_union = None
if x_geos:
    x_union = x_geos[0]
    for g in x_geos[1:]: x_union = x_union.union(g)

# Build trimmed interconnects into temp FC, then append into OUT_FINAL_POLYS
rm(TMP_INTER_TRIM, TMP_FINAL_MERGE)
arcpy.CreateFeatureclass_management(GDB, os.path.basename(TMP_INTER_TRIM), "POLYGON", spatial_reference=SR)
for nm, ftype, prec, scale, length in [
    ("STREET_OID", "LONG", None, None, None),
    ("STREET_NAME", "TEXT", None, None, 128),
    ("STREET_LEVEL", "SHORT", None, None, None),
    ("PART_ID", "LONG", None, None, None),
]:
    arcpy.AddField_management(TMP_INTER_TRIM, nm, ftype, prec, scale, length)

ins_trim = arcpy.da.InsertCursor(TMP_INTER_TRIM, ["SHAPE@", "STREET_OID", "STREET_NAME", "STREET_LEVEL", "PART_ID"])
with arcpy.da.SearchCursor(OUT_INTERCONNECT_POLYS, ["SHAPE@", "STREET_OID", "STREET_NAME", "STREET_LEVEL", "PART_ID"]) as sc:
    for g, oid, nm, lv, pid in sc:
        if is_empty(g): continue
        if x_union:
            diff = g.difference(x_union)
            if not is_empty(diff): ins_trim.insertRow([diff, oid, nm, lv, pid])
        else:
            ins_trim.insertRow([g, oid, nm, lv, pid])
del ins_trim

arcpy.Merge_management([OUT_INTERSECT_POLYS, TMP_INTER_TRIM], TMP_FINAL_MERGE)

ensure_fc(OUT_FINAL_POLYS, "POLYGON", SR)
arcpy.Append_management(TMP_FINAL_MERGE, OUT_FINAL_POLYS, "NO_TEST")

# cleanup temps
rm(TMP_INTER_TRIM, TMP_FINAL_MERGE, "streets_lyr")

print("✅ ArcPy Polygonizer completed (outputs rewritten).")
print(f"CTN_FC: {CTN_FC}")
print(f"- Built intersections (rewritten): {OUT_INTERSECT_PTS}")
print(f"- True intersections (rewritten): {OUT_TRUE_INTERSECT_PTS}")
print(f"- Intersection polygons (rewritten): {OUT_INTERSECT_POLYS}")
print(f"- Interconnect polygons (rewritten): {OUT_INTERCONNECT_POLYS}")
print(f"- Final polygons (rewritten): {OUT_FINAL_POLYS}")
