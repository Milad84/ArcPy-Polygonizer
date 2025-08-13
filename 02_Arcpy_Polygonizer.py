# ArcGIS Pro / ArcPy 3.x
# ArcPy Polygonizer (resumable + cleaned subtraction + street_levels + detailed logging)

import arcpy
from collections import defaultdict
import math
import os
import time

# =========================
# CONFIG
# =========================
GDB   = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb"
CTN_FC = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb\CTN_Test"
SR    = arcpy.SpatialReference(2277)  # NAD 1983 StatePlane Texas Central (US Feet)

# Geometry rules
BUFFER_WIDTH_BY_LEVEL = {1: 40, 2: 40, 3: 50, 4: 75}     # half-width feet
TARGET_LEN_BY_LEVEL   = {1: 150, 2: 150, 3: 250, 4: 250} # segmentation feet
GOAL_LEG_LEN = 200.0
TOUCH_TOL    = 5.0
HALF_GAP_EPS = 0.05
CLEAN_BUF    = 0.1      # gap-killer out–in buffer (ft)

# Force-rebuild flags (set to True to recompute even if output exists)
FORCE_RAW_INTERSECTION_PTS   = False
FORCE_TRUE_INTERSECTION_PTS  = False
FORCE_INTERSECTION_POLYS     = False
FORCE_INTERCONNECT_POLYS     = False
FORCE_TRIM_INTERCONNECTS     = False
FORCE_FINAL_POLYS            = False
FORCE_ENRICH_ATTRIBUTES      = False   # you can rerun attribute write without rebuilding geometries

# Stable outputs
OUT_INTERSECT_PTS      = os.path.join(GDB, "Polygonizer_BuiltIntersections")   # POINT
OUT_TRUE_INTERSECT_PTS = os.path.join(GDB, "Polygonizer_TrueIntersections")    # POINT
OUT_INTERSECT_POLYS    = os.path.join(GDB, "Polygonizer_IntersectionPolys")    # POLYGON
OUT_INTERCONNECT_POLYS = os.path.join(GDB, "Polygonizer_InterconnectPolys")    # POLYGON
OUT_FINAL_POLYS        = os.path.join(GDB, "Polygonizer_FinalPolys")           # POLYGON

# Temporary datasets in GDB (safe for reuse)
TMP_SJ             = os.path.join(GDB, "_pg_tmp_sj")
TMP_RAW_INTERSECT  = os.path.join(GDB, "_pg_tmp_raw_intersections")
TMP_TRUE_INTERSECT = os.path.join(GDB, "_pg_tmp_true_intersections")
TMP_INTER_TRIM     = os.path.join(GDB, "_pg_tmp_inter_trim")
TMP_FINAL_MERGE    = os.path.join(GDB, "_pg_tmp_final_merge")
TMP_SJ_POLY_CTN    = os.path.join(GDB, "_pg_tmp_sj_poly_ctn")
TMP_SJ_POLY_XPTS   = os.path.join(GDB, "_pg_tmp_sj_poly_xpts")
TMP_X_UNION_DISS   = os.path.join(GDB, "_pg_tmp_x_union_diss")
TMP_X_UNION_BUFOUT = os.path.join(GDB, "_pg_tmp_x_union_bufout")
TMP_X_UNION_BUFIN  = os.path.join(GDB, "_pg_tmp_x_union_bufin")

# Environment
arcpy.env.workspace = GDB
arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = SR

# =========================
# LOGGING & HELPERS
# =========================
_start = time.time()
def log(msg):
    stamp = f"[{time.strftime('%H:%M:%S')}] "
    try:
        arcpy.AddMessage(stamp + msg)
    except:
        pass
    print(stamp + msg)

def has_features(fc):
    return arcpy.Exists(fc) and int(arcpy.GetCount_management(fc)[0]) > 0

def field_exists(fc, name):
    return name in [f.name for f in arcpy.ListFields(fc)]

def field_type(fc, name):
    for f in arcpy.ListFields(fc):
        if f.name == name:
            return f.type
    return None

def ensure_text_field(fc, name, length):
    if not field_exists(fc, name):
        arcpy.AddField_management(fc, name, "TEXT", field_length=length)
    elif field_type(fc, name) != "String":
        arcpy.DeleteField_management(fc, [name])
        arcpy.AddField_management(fc, name, "TEXT", field_length=length)

def rm(*paths):
    for p in paths:
        if p and arcpy.Exists(p):
            arcpy.Delete_management(p)

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
    """Copy any point/multipoint features into a point FC."""
    if not arcpy.Exists(dst_fc):
        arcpy.CreateFeatureclass_management(os.path.dirname(dst_fc), os.path.basename(dst_fc), "POINT", spatial_reference=SR)
    icur = arcpy.da.InsertCursor(dst_fc, ["SHAPE@"])
    n = 0
    with arcpy.da.SearchCursor(src_fc, ["SHAPE@"]) as sc:
        for (g,) in sc:
            if is_empty(g): 
                continue
            try:
                if g.type.lower() == "point":
                    icur.insertRow([g]); n += 1; continue
            except: pass
            sr = g.spatialReference
            # try treat as multipoint/polyline/polygon ring
            try:
                arr = g.getPart(0)
                cnt = getattr(arr, "count", None)
                if cnt is not None:
                    for i in range(cnt):
                        pt = arr.getObject(i)
                        if pt:
                            icur.insertRow([arcpy.PointGeometry(pt, sr)]); n += 1
                    continue
            except: pass
            try:
                for pt in g:
                    icur.insertRow([arcpy.PointGeometry(pt, sr)]); n += 1
                continue
            except: pass
            icur.insertRow([arcpy.PointGeometry(g.centroid, sr)]); n += 1
    del icur
    return n

def get_join_count_field(fc):
    names = [f.name for f in arcpy.ListFields(fc)]
    for n in names:
        if n.lower().startswith("join_count"): return n
    for n in names:
        if "join" in n.lower() and "count" in n.lower(): return n
    raise RuntimeError(f"Could not find Join_Count field on {fc}. Fields: {names}")

def find_join_field(join_fc, target_fc, preferred_names):
    """Find the field from the JOIN (CTN) even if Spatial Join suffixed it (e.g., *_1)."""
    tgt_fields = {f.name for f in arcpy.ListFields(target_fc)}
    jf = [f.name for f in arcpy.ListFields(join_fc)]
    # exact and not on target
    for base in preferred_names:
        if base in jf and base not in tgt_fields:
            return base
    # suffixed
    uprefs = [b.upper() for b in preferred_names]
    for name in jf:
        up = name.upper()
        for b in uprefs:
            if up.startswith(b + "_"):
                return name
    # substring fallback
    for name in jf:
        up = name.upper()
        for b in uprefs:
            if b in up:
                return name
    raise RuntimeError(f"Could not locate join field for any of {preferred_names} in {join_fc}")

def round_half_up(x): return int(math.floor(x + 0.5))
def clamp(v, lo, hi): return max(lo, min(hi, v))

def street_level_width(level):
    try: return BUFFER_WIDTH_BY_LEVEL.get(int(level), 40)
    except: return 40

def street_level_target_len(level):
    try: return TARGET_LEN_BY_LEVEL.get(int(level), 150)
    except: return 150

# =========================
# 1) LOAD STREETS
# =========================
log("Checking CTN and loading streets...")
ctn_fields = [f.name for f in arcpy.ListFields(CTN_FC)]
use_full_name = "FULL_STREET_NAME" in ctn_fields
name_field_to_read = "FULL_STREET_NAME" if use_full_name else "STREET_NAME"

streets = []
with arcpy.da.SearchCursor(CTN_FC, ["OID@", "SHAPE@", "STREET_LEVEL", name_field_to_read]) as cur:
    for oid, shp, lvl, nm in cur:
        if is_empty(shp): continue
        streets.append({"oid": oid, "geom": shp, "level": lvl, "name": nm, "length": shp.length})
log(f"Loaded {len(streets)} street features.")
if not streets:
    raise RuntimeError("No street features found in CTN_Test.")
arcpy.MakeFeatureLayer_management(CTN_FC, "streets_lyr")

# =========================
# 2) RAW INTERSECTION CANDIDATES (POINTS)
# =========================
if not has_features(OUT_INTERSECT_PTS) or FORCE_RAW_INTERSECTION_PTS:
    log("Computing raw intersection candidates (line-line)…")
    rm(TMP_RAW_INTERSECT)
    arcpy.Intersect_analysis([CTN_FC], TMP_RAW_INTERSECT, "ALL", f"{TOUCH_TOL} Feet", "POINT")
    try: arcpy.DeleteIdentical_management(TMP_RAW_INTERSECT, ["SHAPE"])
    except: pass
    try: arcpy.management.Integrate([TMP_RAW_INTERSECT, CTN_FC], f"{TOUCH_TOL} Feet")
    except: pass

    rm(OUT_INTERSECT_PTS)
    arcpy.CreateFeatureclass_management(os.path.dirname(OUT_INTERSECT_PTS), os.path.basename(OUT_INTERSECT_PTS),"POINT",spatial_reference=SR)
    npts = explode_points_to_point_fc(TMP_RAW_INTERSECT, OUT_INTERSECT_PTS)
    rm(TMP_RAW_INTERSECT)
    log(f"Raw candidate points: {npts}")
else:
    log("Raw intersection points already exist—skipping.")

# =========================
# 3) TRUE INTERSECTION POINTS (Join_Count >= 2)
# =========================
if not has_features(OUT_TRUE_INTERSECT_PTS) or FORCE_TRUE_INTERSECTION_PTS:
    log("Selecting true intersections (>=2 incident streets)…")
    rm(TMP_SJ, TMP_TRUE_INTERSECT, OUT_TRUE_INTERSECT_PTS)
    arcpy.SpatialJoin_analysis(
        target_features=OUT_INTERSECT_PTS,
        join_features=CTN_FC,
        out_feature_class=TMP_SJ,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_COMMON",
        match_option="INTERSECT",
        search_radius=f"{TOUCH_TOL} Feet"
    )
    JOIN_COUNT_FIELD = get_join_count_field(TMP_SJ)
    arcpy.FeatureClassToFeatureClass_conversion(
        TMP_SJ, GDB, os.path.basename(TMP_TRUE_INTERSECT),
        where_clause=f"{arcpy.AddFieldDelimiters(TMP_SJ, JOIN_COUNT_FIELD)} >= 2"
    )
    rm(TMP_SJ)

    arcpy.CreateFeatureclass_management(os.path.dirname(OUT_TRUE_INTERSECT_PTS), os.path.basename(OUT_TRUE_INTERSECT_PTS),"POINT",spatial_reference=SR)
    ntrue = explode_points_to_point_fc(TMP_TRUE_INTERSECT, OUT_TRUE_INTERSECT_PTS)
    rm(TMP_TRUE_INTERSECT)
    log(f"True intersections: {ntrue}")
else:
    log("True intersection points already exist—skipping.")

# Build in-memory node list (always; it’s fast)
x_oid_field = arcpy.Describe(OUT_TRUE_INTERSECT_PTS).OIDFieldName
true_nodes = []
with arcpy.da.SearchCursor(OUT_TRUE_INTERSECT_PTS, [x_oid_field, "SHAPE@"]) as cur:
    for oid, g in cur:
        if is_empty(g): continue
        true_nodes.append({"oid": oid, "geom": g})
if not true_nodes:
    raise RuntimeError("No true intersections found (>= 2 incident streets).")

# =========================
# 4) MAP NODES TO STREETS & LEG LENGTHS
# =========================
log("Mapping intersections to incident streets…")
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
log("Incident lists built.")

log("Computing non-overlapping leg lengths per incident street…")
computed_leg_len = defaultdict(dict)
street_lookup = {s["oid"]: s for s in streets}
for s in streets:
    soid = s["oid"]
    if soid not in street_to_intersections: continue
    entries  = sorted(street_to_intersections[soid], key=lambda x: x[1])
    measures = [m for (_pid, m, _pt) in entries]
    for i, (pid, m, _pt) in enumerate(entries):
        nn_gap = s["length"] if len(measures) == 1 else min(abs(m - measures[j]) for j in range(len(measures)) if j != i)
        leg_len = min(GOAL_LEG_LEN, 0.5 * nn_gap - HALF_GAP_EPS, 0.5 * s["length"])
        computed_leg_len[soid][pid] = max(0.0, leg_len)
log("Leg lengths computed.")

def extract_leg(street_geom, meet_point, leg_len):
    m = street_geom.measureOnLine(meet_point, False)
    L = street_geom.length
    if m <= (L - m):
        start_m, end_m = 0.0, clamp(leg_len, 0.0, L)
    else:
        start_m, end_m = clamp(L - leg_len, 0.0, L), L
    return street_geom.segmentAlongLine(start_m, end_m, False)

# =========================
# 5) INTERSECTION POLYGONS
# =========================
if not has_features(OUT_INTERSECT_POLYS) or FORCE_INTERSECTION_POLYS:
    log("Building intersection polygons (plus shapes with legs)…")
    rm(OUT_INTERSECT_POLYS)
    arcpy.CreateFeatureclass_management(GDB, os.path.basename(OUT_INTERSECT_POLYS), "POLYGON", spatial_reference=SR)
    for fld in [("X_ID","LONG"),("LEGS","SHORT"),("NOTE","TEXT",64),("location_name","TEXT",250),("street_levels","TEXT",50)]:
        if len(fld)==2: arcpy.AddField_management(OUT_INTERSECT_POLYS,fld[0],fld[1])
        else:           arcpy.AddField_management(OUT_INTERSECT_POLYS,fld[0],fld[1], field_length=fld[2])
    ins_int = arcpy.da.InsertCursor(OUT_INTERSECT_POLYS, ["SHAPE@", "X_ID", "LEGS", "NOTE", "location_name", "street_levels"])
    n_inter_polys = 0
    for node in true_nodes:
        center = as_point_geometry(node["geom"])
        if center is None: continue

        incidents = [soid for soid, lst in street_to_intersections.items()
                     if any(pid == node["oid"] for pid, _m, _pt in lst)]
        if not incidents: continue

        name_set, level_set = set(), set()
        for soid in incidents:
            srec = street_lookup[soid]
            if srec["name"]:
                name_set.add(str(srec["name"]).strip())
            if srec["level"] is not None:
                try: level_set.add(int(srec["level"]))
                except: level_set.add(str(srec["level"]))

        leg_buffers = []
        for soid in incidents:
            srec   = street_lookup[soid]
            width  = street_level_width(srec["level"])
            leglen = computed_leg_len[soid].get(node["oid"], 0.0)
            if leglen <= 0.0: continue
            leg_line = extract_leg(srec["geom"], center, leglen)
            tmp = arcpy.Buffer_analysis(leg_line, os.path.join("in_memory", "_leg"), f"{width} Feet",
                                        line_side="FULL", line_end_type="FLAT", dissolve_option="NONE")
            with arcpy.da.SearchCursor(tmp, ["SHAPE@"]) as bc:
                for (g,) in bc:
                    if not is_empty(g): leg_buffers.append(g)
            arcpy.Delete_management(tmp)

        max_w = max(street_level_width(street_lookup[soid]["level"]) for soid in incidents)
        if not leg_buffers:
            merged_geom = center.buffer(max_w)
        else:
            merged_geom = leg_buffers[0]
            for g in leg_buffers[1:]:
                merged_geom = merged_geom.union(g)
            merged_geom = merged_geom.union(center.buffer(max_w))

        loc_txt = ", ".join(sorted(name_set))[:250] if name_set else None
        try:
            levels_sorted = sorted(level_set, key=lambda v: int(v))
        except:
            levels_sorted = sorted(map(str, level_set))
        lvl_txt = ", ".join(str(v) for v in levels_sorted) if level_set else None

        if not is_empty(merged_geom):
            ins_int.insertRow([merged_geom, node["oid"], len(incidents), "OK", loc_txt, lvl_txt])
            n_inter_polys += 1
    del ins_int
    log(f"Intersection polygons created: {n_inter_polys}")
else:
    log("Intersection polygons already exist—skipping.")

# =========================
# 6) INTERCONNECT POLYGONS
# =========================
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

if not has_features(OUT_INTERCONNECT_POLYS) or FORCE_INTERCONNECT_POLYS:
    log("Building interconnect polygons (segmented buffers)…")
    rm(OUT_INTERCONNECT_POLYS)
    arcpy.CreateFeatureclass_management(GDB, os.path.basename(OUT_INTERCONNECT_POLYS), "POLYGON", spatial_reference=SR)
    for fld in [("STREET_OID","LONG"),("STREET_NAME","TEXT",128),("STREET_LEVEL","SHORT"),("PART_ID","LONG")]:
        if len(fld)==2: arcpy.AddField_management(OUT_INTERCONNECT_POLYS,fld[0],fld[1])
        else:           arcpy.AddField_management(OUT_INTERCONNECT_POLYS,fld[0],fld[1], field_length=fld[2])

    ins_ic = arcpy.da.InsertCursor(OUT_INTERCONNECT_POLYS, ["SHAPE@", "STREET_OID", "STREET_NAME", "STREET_LEVEL", "PART_ID"])
    n_interconnect_polys = 0
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

            tmp = arcpy.Buffer_analysis(seg, os.path.join("in_memory", "_seg"), f"{width} Feet",
                                        line_side="FULL", line_end_type="FLAT", dissolve_option="NONE")
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
            for cap in add_caps:
                seg_poly = seg_poly.union(cap)

            if not is_empty(seg_poly):
                ins_ic.insertRow([seg_poly, soid, name, lvl, i + 1])
                n_interconnect_polys += 1
    del ins_ic
    log(f"Interconnect polygons created: {n_interconnect_polys}")
else:
    log("Interconnect polygons already exist—skipping.")

# =========================
# 7) CLEAN SUBTRACTION (DISSOLVE + OUT–IN BUFFER) & MERGE
# =========================
if (not has_features(TMP_INTER_TRIM) or not has_features(TMP_FINAL_MERGE) or FORCE_TRIM_INTERCONNECTS):
    log("Dissolving intersection polygons…")
    rm(TMP_X_UNION_DISS, TMP_X_UNION_BUFOUT, TMP_X_UNION_BUFIN)
    arcpy.Dissolve_management(OUT_INTERSECT_POLYS, TMP_X_UNION_DISS)

    log("Cleaning dissolved surface with out–in buffer…")
    arcpy.Buffer_analysis(TMP_X_UNION_DISS, TMP_X_UNION_BUFOUT, f"{CLEAN_BUF} Feet")
    arcpy.Buffer_analysis(TMP_X_UNION_BUFOUT, TMP_X_UNION_BUFIN, f"-{CLEAN_BUF} Feet")

    # Read cleaned surface geometry
    x_clean_geom = None
    with arcpy.da.SearchCursor(TMP_X_UNION_BUFIN, ["SHAPE@"]) as xc:
        for (g,) in xc:
            if not is_empty(g):
                x_clean_geom = g
                break

    rm(TMP_INTER_TRIM, TMP_FINAL_MERGE)
    if x_clean_geom is None:
        log("Warning: cleaned surface empty; copying interconnects as-is.")
        arcpy.CopyFeatures_management(OUT_INTERCONNECT_POLYS, TMP_INTER_TRIM)
    else:
        log("Subtracting cleaned intersection surface from interconnect polygons…")
        arcpy.CreateFeatureclass_management(GDB, os.path.basename(TMP_INTER_TRIM), "POLYGON", spatial_reference=SR)
        for fld in [("STREET_OID","LONG"),("STREET_NAME","TEXT",128),("STREET_LEVEL","SHORT"),("PART_ID","LONG")]:
            if len(fld)==2: arcpy.AddField_management(TMP_INTER_TRIM,fld[0],fld[1])
            else:           arcpy.AddField_management(TMP_INTER_TRIM,fld[0],fld[1], field_length=fld[2])

        ins_trim = arcpy.da.InsertCursor(TMP_INTER_TRIM, ["SHAPE@", "STREET_OID", "STREET_NAME", "STREET_LEVEL", "PART_ID"])
        n_trimmed = 0
        with arcpy.da.SearchCursor(OUT_INTERCONNECT_POLYS, ["SHAPE@", "STREET_OID", "STREET_NAME", "STREET_LEVEL", "PART_ID"]) as sc:
            for g, oid, nm, lv, pid in sc:
                if is_empty(g): continue
                diff = g.difference(x_clean_geom)
                if not is_empty(diff):
                    ins_trim.insertRow([diff, oid, nm, lv, pid]); n_trimmed += 1
        del ins_trim
        log(f"Interconnect polygons trimmed: {n_trimmed}")

    # Merge intersection polys + trimmed interconnects
    arcpy.Merge_management([OUT_INTERSECT_POLYS, TMP_INTER_TRIM], TMP_FINAL_MERGE)
else:
    log("Trimmed interconnects and merge already exist—skipping.")

# =========================
# 8) FINAL POLYS (rewrite once; or skip if exists)
# =========================
if not has_features(OUT_FINAL_POLYS) or FORCE_FINAL_POLYS:
    log("Writing final polygon layer…")
    rm(OUT_FINAL_POLYS)
    arcpy.CreateFeatureclass_management(GDB, os.path.basename(OUT_FINAL_POLYS), "POLYGON", spatial_reference=SR)
    ensure_text_field(OUT_FINAL_POLYS, "location_name", 250)
    ensure_text_field(OUT_FINAL_POLYS, "street_levels", 50)
    if not field_exists(OUT_FINAL_POLYS, "is_intersection"):
        arcpy.AddField_management(OUT_FINAL_POLYS, "is_intersection", "SHORT")
    arcpy.Append_management(TMP_FINAL_MERGE, OUT_FINAL_POLYS, "NO_TEST")
else:
    log("Final polygons already exist—skipping rebuild; will just update attributes.")

# =========================
# 9) ATTRIBUTE ENRICHMENT (coalesce; SINGLE 'street_levels')
# =========================
if not has_features(OUT_FINAL_POLYS):
    raise RuntimeError("Final polygons missing; cannot enrich attributes.")

if FORCE_ENRICH_ATTRIBUTES or True:
    log("Enriching attributes (names/levels; coalescing with prefilled intersections)…")
    # Spatial Join polygons <- CTN (1:M)
    rm(TMP_SJ_POLY_CTN)
    arcpy.SpatialJoin_analysis(
        target_features=OUT_FINAL_POLYS,
        join_features=CTN_FC,
        out_feature_class=TMP_SJ_POLY_CTN,
        join_operation="JOIN_ONE_TO_MANY",
        join_type="KEEP_COMMON",
        match_option="INTERSECT"
    )
    poly_oid_name = arcpy.Describe(OUT_FINAL_POLYS).OIDFieldName
    target_id_field = "TARGET_FID" if field_exists(TMP_SJ_POLY_CTN, "TARGET_FID") else poly_oid_name

    # Detect CTN join fields
    try:
        JOIN_FULL_NAME = find_join_field(TMP_SJ_POLY_CTN, OUT_FINAL_POLYS, ["FULL_STREET_NAME"])
    except:
        JOIN_FULL_NAME = find_join_field(TMP_SJ_POLY_CTN, OUT_FINAL_POLYS, ["STREET_NAME"])
    JOIN_LEVEL = find_join_field(TMP_SJ_POLY_CTN, OUT_FINAL_POLYS, ["STREET_LEVEL"])

    names_by_poly  = defaultdict(set)
    levels_by_poly = defaultdict(set)
    with arcpy.da.SearchCursor(TMP_SJ_POLY_CTN, [target_id_field, JOIN_FULL_NAME, JOIN_LEVEL]) as sc:
        for pid, sname, slevel in sc:
            if pid is None: continue
            if sname:
                names_by_poly[pid].add(str(sname).strip())
            if slevel is not None:
                try: levels_by_poly[pid].add(int(slevel))
                except: levels_by_poly[pid].add(str(slevel))

    # is_intersection via Spatial Join to true intersection points (1:1)
    rm(TMP_SJ_POLY_XPTS)
    arcpy.SpatialJoin_analysis(
        target_features=OUT_FINAL_POLYS,
        join_features=OUT_TRUE_INTERSECT_PTS,
        out_feature_class=TMP_SJ_POLY_XPTS,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        match_option="INTERSECT"
    )
    join_count_field2 = get_join_count_field(TMP_SJ_POLY_XPTS)
    is_intersection_by_poly = {}
    with arcpy.da.SearchCursor(TMP_SJ_POLY_XPTS, [poly_oid_name, join_count_field2]) as sc2:
        for oid_val, cnt in sc2:
            try:
                is_intersection_by_poly[oid_val] = 1 if (cnt and int(cnt) >= 1) else 0
            except:
                is_intersection_by_poly[oid_val] = 0

    # Ensure final fields exist
    ensure_text_field(OUT_FINAL_POLYS, "location_name", 250)
    ensure_text_field(OUT_FINAL_POLYS, "street_levels", 50)
    if not field_exists(OUT_FINAL_POLYS, "is_intersection"):
        arcpy.AddField_management(OUT_FINAL_POLYS, "is_intersection", "SHORT")

    # Coalesce & write back
    count_upd = 0
    with arcpy.da.UpdateCursor(OUT_FINAL_POLYS, [poly_oid_name, "location_name", "street_levels", "is_intersection"]) as uc:
        for oid_val, locname_existing, levels_existing, isx in uc:
            # Names
            names_join = sorted(n for n in names_by_poly.get(oid_val, set()) if n)
            loc_final = ", ".join(names_join)[:250] if names_join else locname_existing

            # Levels (unique list -> TEXT)
            lvls_join = list(levels_by_poly.get(oid_val, set()))
            if lvls_join:
                try:  lvls_join = sorted(lvls_join, key=lambda v: int(v))
                except: lvls_join = sorted(map(str, lvls_join))
                levels_final = ", ".join(str(v) for v in lvls_join)
            else:
                levels_final = levels_existing  # keep prefilled (intersections)

            ix = is_intersection_by_poly.get(oid_val, isx if isx is not None else 0)
            uc.updateRow([oid_val, loc_final, levels_final, ix])
            count_upd += 1
    log(f"Attributes updated on {count_upd} polygons.")
    rm(TMP_SJ_POLY_CTN, TMP_SJ_POLY_XPTS)

# =========================
# DONE
# =========================
elapsed = time.time() - _start
log(f"✅ Completed. Elapsed: {elapsed:0.1f} sec")
log(f"Outputs:")
log(f"  • Intersections (points):   {OUT_TRUE_INTERSECT_PTS}")
log(f"  • Intersection polygons:    {OUT_INTERSECT_POLYS}")
log(f"  • Interconnect polygons:    {OUT_INTERCONNECT_POLYS}")
log(f"  • Final polygons:           {OUT_FINAL_POLYS}")
log("Fields on final: location_name (TEXT), street_levels (TEXT), is_intersection (SHORT)")
