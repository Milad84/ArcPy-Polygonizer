

<img width="823" height="626" alt="image" src="https://github.com/user-attachments/assets/492d7507-c232-4257-9ef0-dfc506c6d133" />

# ArcPy Polygonizer (ArcGIS Pro / ArcPy 3.x)

Generate clean, attribute-rich street area polygons from a centerline network — fast, reproducible, and resumable.
This implementation uses pairwise geoprocessing tools, a robust in-GDB near table, and FLAT-ended line buffers. The final layer contains a single street_levels field and an is_intersection flag.

# Highlights

Resumable: each stage can be skipped unless forced to rebuild.

Faster: Pairwise* tools + parallelProcessingFactor="100%".

Robust intersections: GenerateNearTable written to the GDB (not in_memory).

FLAT buffer ends for legs and interconnects (matches baseline cartography).

Clean attributes: location_name, single street_levels, is_intersection.

# Requirements

ArcGIS Pro 3.x with ArcPy.

Feet-based projected CRS (default: EPSG 2277 – NAD 1983 StatePlane Texas Central (US Feet)).

# Tools used:

PairwiseIntersect, PairwiseDissolve, PairwiseBuffer

GenerateNearTable

SpatialJoin, Merge, Buffer, Append, Dissolve

If a pairwise tool isn’t available in your environment, swap to the classic equivalent (potentially slower).

# Input Assumptions

CTN_FC is a polyline centerline layer with:

STREET_LEVEL (integer or convertible text).

A name field: prefers FULL_STREET_NAME, otherwise STREET_NAME.

Centerlines are reasonably snapped (default matching tolerance TOUCH_TOL = 5 ft).

#Outputs

All outputs are written to the target geodatabase.

Dataset	Type	Description
Polygonizer_BuiltIntersections	Point	Raw self line–line candidates from PairwiseIntersect
Polygonizer_TrueIntersections	Point	Points with ≥ 2 incident streets within TOUCH_TOL
Polygonizer_IntersectionPolys	Polygon	Node cores + FLAT-ended leg buffers (unioned)
Polygonizer_InterconnectPolys	Polygon	Segmented FLAT-ended buffers between intersections
Polygonizer_FinalPolys	Polygon	Final merged surface with attributes

# Final attribute schema

location_name (TEXT, 250): comma-joined unique street names intersecting the polygon.

street_levels (TEXT, 50): sorted, comma-joined unique levels.

is_intersection (SHORT): 1 if the polygon intersects ≥1 true-intersection point; else 0.

The writer guarantees exactly one street_levels field on the final output.

# Configuration

Edit these values at the top of the script:

```python
# Paths & spatial reference
GDB   = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb"
CTN_FC = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb\CTN_Test"
SR    = arcpy.SpatialReference(2277)  # projected in US Feet

# Geometry rules
BUFFER_WIDTH_BY_LEVEL = {1: 40, 2: 40, 3: 50, 4: 75}     # half-width (ft)
TARGET_LEN_BY_LEVEL   = {1: 150, 2: 150, 3: 250, 4: 250} # interconnect segment length (ft)
GOAL_LEG_LEN = 200.0
TOUCH_TOL    = 5.0
HALF_GAP_EPS = 0.05
CLEAN_BUF    = 0.1  # tiny out–in buffer used only on the dissolved intersection surface

# Rebuild controls (set False to reuse existing outputs)
FORCE_RAW_INTERSECTION_PTS   = True
FORCE_TRUE_INTERSECTION_PTS  = True
FORCE_INTERSECTION_POLYS     = True
FORCE_INTERCONNECT_POLYS     = True
FORCE_TRIM_INTERCONNECTS     = True
FORCE_FINAL_POLYS            = True
FORCE_ENRICH_ATTRIBUTES      = True
```

# Tips

If matching seems inconsistent, try TOUCH_TOL = 8–10.

Customize widths & segment lengths per STREET_LEVEL.

CLEAN_BUF = 0.1 ft gently removes micro-slivers after dissolve.

Quick Start

Open ArcGIS Pro 3.x → Python Window / Notebook.

Paste the script (or run as a Script Tool).

Set GDB, CTN_FC, and (if needed) SR.

Run with defaults for a full rebuild.

To only refresh attributes later:

FORCE_RAW_INTERSECTION_PTS   = False
FORCE_TRUE_INTERSECTION_PTS  = False
FORCE_INTERSECTION_POLYS     = False
FORCE_INTERCONNECT_POLYS     = False
FORCE_TRIM_INTERCONNECTS     = False
FORCE_FINAL_POLYS            = False
FORCE_ENRICH_ATTRIBUTES      = True

How It Works (Pipeline)
1) Load Streets

Reads OID@, SHAPE@, STREET_LEVEL, and a name field (FULL_STREET_NAME or STREET_NAME), caching length and level per feature.

2) Raw Intersection Candidates

PairwiseIntersect([CTN_FC], …, "POINT") finds self intersections.
Deduplicate with DeleteIdentical. Write points to Polygonizer_BuiltIntersections.

3) True Intersection Points

GenerateNearTable (raw points ↔ CTN within TOUCH_TOL) is written to the GDB and verified.
Count incident streets by IN_FID; keep points with count ≥ 2 → Polygonizer_TrueIntersections.

4) Map Nodes to Streets & Leg Lengths

Using the near table, associate each node to incident streets and compute non-overlapping leg lengths:

leg_len = min(
  GOAL_LEG_LEN,
  0.5 * nearest_gap - HALF_GAP_EPS,
  0.5 * street_length
)


extract_leg() keeps the half that runs from the node toward the nearer end.

5) Intersection Polygons (FLAT Legs)

For each node:

Build a circular core using the max half-width among incident street levels.

Buffer each incident leg with Buffer_analysis(..., line_end_type="FLAT") at that level’s half-width.

Union the core with all leg buffers and write X_ID, LEGS, location_name, street_levels.

6) Interconnect Polygons (FLAT Segments)

For each street:

Detect nodes at ends; carve off end legs.

Split the remaining center segment into n_parts using TARGET_LEN_BY_LEVEL[level].

Buffer each segment with FLAT end caps.

If a part touches an open end (no node), add a round cap by unioning a point buffer (cul-de-sac behavior).

Write STREET_OID, STREET_NAME, STREET_LEVEL, PART_ID.

7) Clean Subtraction & Merge

PairwiseDissolve all intersection polygons to a single surface.

Out–in PairwiseBuffer (+CLEAN_BUF then -CLEAN_BUF) removes micro artifacts.

For each interconnect polygon, subtract the cleaned surface via geometry.difference().

Merge intersections + trimmed interconnects into a staging FC.

8) Final Layer Writing

Create Polygonizer_FinalPolys with the target schema and append the staging FC.

Defensively delete any suffixed duplicates of street_levels (rare).

9) Attribute Enrichment

SpatialJoin (1:M) FinalPolys ← CTN: collect unique names & levels per polygon.

SpatialJoin (1:1) FinalPolys ← TrueIntersections: compute is_intersection.

Update:

location_name = sorted, comma-joined names (≤ 250 chars)

street_levels = sorted, comma-joined unique levels

is_intersection = 0/1

Performance Choices

Pairwise tools where possible: PairwiseIntersect, PairwiseDissolve, PairwiseBuffer.

Near table is persisted in the GDB (stable, lock-resilient).

FLAT buffer caps for legs/interconnects match the intended cartography.

arcpy.env.parallelProcessingFactor = "100%" enables multi-core where supported.

Resumability (Force Flags)

Common modes:

# Full rebuild (first run)
FORCE_RAW_INTERSECTION_PTS   = True
FORCE_TRUE_INTERSECTION_PTS  = True
FORCE_INTERSECTION_POLYS     = True
FORCE_INTERCONNECT_POLYS     = True
FORCE_TRIM_INTERCONNECTS     = True
FORCE_FINAL_POLYS            = True
FORCE_ENRICH_ATTRIBUTES      = True

# Geometry frozen, refresh attributes only
FORCE_RAW_INTERSECTION_PTS   = False
FORCE_TRUE_INTERSECTION_PTS  = False
FORCE_INTERSECTION_POLYS     = False
FORCE_INTERCONNECT_POLYS     = False
FORCE_TRIM_INTERCONNECTS     = False
FORCE_FINAL_POLYS            = False
FORCE_ENRICH_ATTRIBUTES      = True

Troubleshooting

PairwiseDissolve not found → Use arcpy.analysis.PairwiseDissolve (Analysis toolbox), not arcpy.management.

Near table cursor error (cannot open ... _pg_tmp_near) → The script writes the near table to the GDB and checks existence; inspect geoprocessing messages printed right after GenerateNearTable.

MakeFeatureLayer crash → No field_info is passed; a simple view avoids brittle mappings.

Duplicate street_levels → Final writer creates the field once, appends, then deletes any suffixed duplicates defensively.

Unexpected shapes → Ensure US Feet and arcpy.env.outputCoordinateSystem = SR. Legs & interconnects use FLAT caps; cul-de-sacs are rounded by design.

Design Notes & Constraints

Leg allocator prevents overlap by capping at:

½ of the nearest inter-node gap minus HALF_GAP_EPS

½ of the street length

GOAL_LEG_LEN

Interconnect segmentation yields consistent part sizes by STREET_LEVEL.

Cleaning uses a tiny out–in buffer on the intersection dissolve surface only, minimizing drift.

Temp datasets are prefixed _pg_tmp_* and can be deleted safely between runs.

License

Add a license of your choice (e.g., MIT) in a LICENSE file.

Code Reference

This README documents the last working script titled:

ArcPy Polygonizer — SPEED-OPTIMIZED + robust near-table (in-GDB),
pairwise tools, FLAT-ended buffers (match baseline geometry), single 'street_levels'


Place the script in your repo (e.g., polygonizer.py) and edit the Configuration block before running.

You said:
ThI need the headings to be formatted 
ChatGPT said:
Thought for 13s

Got it — here’s a copy-paste-ready README.md with proper Markdown headings (#, ##, ###) and nothing fancy. Drop this straight into the GitHub editor.

ArcPy Polygonizer (ArcGIS Pro / ArcPy 3.x)

Generate clean, attribute-rich street area polygons from a centerline network — fast, reproducible, and resumable.
This implementation uses pairwise geoprocessing tools, a robust in-GDB near table, and FLAT-ended line buffers. The final layer contains a single street_levels field and an is_intersection flag.

Highlights

Resumable: each stage can be skipped unless forced to rebuild.

Faster: Pairwise* tools + parallelProcessingFactor="100%".

Robust intersections: GenerateNearTable written to the GDB (not in_memory).

FLAT buffer ends for legs and interconnects (matches baseline cartography).

Clean attributes: location_name, single street_levels, is_intersection.

Requirements

ArcGIS Pro 3.x with ArcPy.

Feet-based projected CRS (default: EPSG 2277 – NAD 1983 StatePlane Texas Central (US Feet)).

Tools used:

PairwiseIntersect, PairwiseDissolve, PairwiseBuffer

GenerateNearTable

SpatialJoin, Merge, Buffer, Append, Dissolve

If a pairwise tool isn’t available in your environment, swap to the classic equivalent (potentially slower).

Input Assumptions

CTN_FC is a polyline centerline layer with:

STREET_LEVEL (integer or convertible text).

A name field: prefers FULL_STREET_NAME, otherwise STREET_NAME.

Centerlines are reasonably snapped (default TOUCH_TOL = 5 ft).

Outputs

All outputs are written to the target geodatabase.

Dataset	Type	Description
Polygonizer_BuiltIntersections	Point	Raw self line–line candidates from PairwiseIntersect
Polygonizer_TrueIntersections	Point	Points with ≥ 2 incident streets within TOUCH_TOL
Polygonizer_IntersectionPolys	Polygon	Node cores + FLAT-ended leg buffers (unioned)
Polygonizer_InterconnectPolys	Polygon	Segmented FLAT-ended buffers between intersections
Polygonizer_FinalPolys	Polygon	Final merged surface with attributes

Final attribute schema

location_name (TEXT, 250): comma-joined unique street names intersecting the polygon.

street_levels (TEXT, 50): sorted, comma-joined unique levels.

is_intersection (SHORT): 1 if the polygon intersects ≥1 true-intersection point; else `

You said:
You are still using bullet points...Why do you do this to me?
ChatGPT said:
Thought for 26s
ArcPy Polygonizer (ArcGIS Pro / ArcPy 3.x)

This project generates clean, attribute-rich street-area polygons from a centerline network. It produces intersection polygons around multi-leg nodes using FLAT-ended line buffers, interconnect polygons along streets between intersections segmented by street level and target length, subtracts cleaned intersection surfaces from interconnects, merges the results, and enriches the final polygons with names, levels, and an is_intersection flag. The implementation uses pairwise geoprocessing tools where possible, stores the near table in the geodatabase for robustness, and guarantees a single street_levels field on the final output. This README documents the last working script and intentionally excludes any dent-fix post-processing.

Requirements

The workflow requires ArcGIS Pro 3.x with ArcPy available and a feet-based projected coordinate system. By default the script uses EPSG 2277 (NAD 1983 StatePlane Texas Central, US Feet) and sets arcpy.env.outputCoordinateSystem accordingly. The tools used include PairwiseIntersect, PairwiseDissolve, PairwiseBuffer, GenerateNearTable, SpatialJoin, Buffer, Merge, Dissolve, Append, and standard ArcPy cursors. If a pairwise tool is missing in your environment, you may substitute the classic equivalent at the cost of performance. The script enables parallelism by setting arcpy.env.parallelProcessingFactor = "100%".

Input assumptions

The input CTN_FC is a polyline centerline feature class. It contains a street level field (STREET_LEVEL, integer or text convertible to integer) that drives buffer half-widths and interconnect segmentation lengths, and a name field. The script prefers FULL_STREET_NAME and falls back to STREET_NAME when the full name is not present. Centerlines should be reasonably snapped for a default node-to-street tolerance of five feet to work reliably.

Outputs

All datasets are written into the configured file geodatabase.

Dataset	Type	Purpose
Polygonizer_BuiltIntersections	Point	Raw self line–line candidates from PairwiseIntersect written as points.
Polygonizer_TrueIntersections	Point	Points that have at least two incident streets within the tolerance, derived from a near table.
Polygonizer_IntersectionPolys	Polygon	Union of a circular node core and FLAT-ended leg buffers sized by street level.
Polygonizer_InterconnectPolys	Polygon	Segmented FLAT-ended buffers along streets between intersections, with optional round caps for cul-de-sacs.
Polygonizer_FinalPolys	Polygon	Merged and cleaned polygons with attributes location_name, street_levels, and is_intersection.

The final attribute schema consists of location_name (TEXT, length 250) that stores comma-joined unique street names intersecting the polygon, street_levels (TEXT, length 50) that stores sorted, comma-joined unique levels, and is_intersection (SHORT) that is set to one if a polygon intersects at least one true-intersection point and zero otherwise. The writer guarantees exactly one street_levels column on the final output by creating the field up front, appending, and deleting any suffixed duplicates defensively.

Configuration

The configuration block at the top of the script defines the file geodatabase path, the centerline feature class, the spatial reference, and the geometry rules governing buffer widths by level, segmentation targets by level, leg length caps, and tolerances. It also holds force-rebuild flags for each stage so you can rerun only what you need.

# Paths & spatial reference
GDB   = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb"
CTN_FC = r"G:\ATD\ACTIVE TRANS\Vision Zero\GIS\New Vision Zero Polygons_2025\ArcPyPolygonizer.gdb\CTN_Test"
SR    = arcpy.SpatialReference(2277)  # US Feet

# Geometry rules
BUFFER_WIDTH_BY_LEVEL = {1: 40, 2: 40, 3: 50, 4: 75}     # half-width (ft)
TARGET_LEN_BY_LEVEL   = {1: 150, 2: 150, 3: 250, 4: 250} # segment length (ft)
GOAL_LEG_LEN = 200.0
TOUCH_TOL    = 5.0
HALF_GAP_EPS = 0.05
CLEAN_BUF    = 0.1  # tiny out–in buffer on the dissolved intersection surface

# Rebuild controls
FORCE_RAW_INTERSECTION_PTS   = True
FORCE_TRUE_INTERSECTION_PTS  = True
FORCE_INTERSECTION_POLYS     = True
FORCE_INTERCONNECT_POLYS     = True
FORCE_TRIM_INTERCONNECTS     = True
FORCE_FINAL_POLYS            = True
FORCE_ENRICH_ATTRIBUTES      = True


If matching appears inconsistent in your data, increase TOUCH_TOL to a value between eight and ten feet. Update BUFFER_WIDTH_BY_LEVEL and TARGET_LEN_BY_LEVEL to reflect your roadway hierarchy. The CLEAN_BUF value of one tenth of a foot is sufficient to remove microscopic artifacts after dissolving intersections while avoiding measurable geometry drift.

Quick start

Open ArcGIS Pro 3.x and a Python Window or Notebook. Paste the script or run it as a Script Tool. Adjust GDB, CTN_FC, and SR if needed. Run the script with the default force flags for a full build the first time. For subsequent runs when you only need to refresh attributes while preserving geometry, set all geometry-stage force flags to False and leave FORCE_ENRICH_ATTRIBUTES set to True.

FORCE_RAW_INTERSECTION_PTS   = False
FORCE_TRUE_INTERSECTION_PTS  = False
FORCE_INTERSECTION_POLYS     = False
FORCE_INTERCONNECT_POLYS     = False
FORCE_TRIM_INTERCONNECTS     = False
FORCE_FINAL_POLYS            = False
FORCE_ENRICH_ATTRIBUTES      = True

Processing pipeline
Step 1 — Load streets

The script reads OID@, SHAPE@, STREET_LEVEL, and either FULL_STREET_NAME or STREET_NAME, discards empty geometry, caches per-street length, and builds a lookup by object ID. It also prepares a feature layer for spatial selections if needed later.

Step 2 — Build raw intersection candidates

The pipeline uses arcpy.analysis.PairwiseIntersect([CTN_FC], ..., "POINT") to compute self line–line intersections. After removing duplicates with DeleteIdentical, it writes the results as a simple point feature class named Polygonizer_BuiltIntersections.

Step 3 — Derive true intersection points

The script calls arcpy.analysis.GenerateNearTable between the raw candidate points and the centerline feature class using the configured tolerance. The near table is written into the geodatabase rather than in_memory for stability. It counts the number of incident streets per candidate point and keeps only those with counts greater than or equal to two, which are copied into Polygonizer_TrueIntersections.

Step 4 — Map nodes to streets and compute leg lengths

Using the near table pairs, each true intersection node is associated with its incident streets. For each node–street pairing the script computes a chainage measure along the street and derives a non-overlapping leg length using the minimum of the goal leg length, half of the nearest inter-node gap minus a small epsilon, and half of the street length. A helper extracts the shorter side segment from the node toward the nearer end.

Step 5 — Build intersection polygons with FLAT legs

For each true intersection node the script constructs a circular core using the maximum buffer half-width among the incident street levels. It extracts each leg segment, buffers it at the level’s half-width with line_end_type="FLAT", and unions those legs with the node core. It records the node identifier, the number of incident legs, a comma-joined location_name of unique incident street names, and a comma-joined street_levels of unique incident levels.

Step 6 — Build interconnect polygons with FLAT segments

For every street the script determines whether nodes lie at its ends and trims the computed end legs away from the center. It then splits the remaining center length into parts based on the level’s target length. Each center segment is buffered with FLAT end caps. If a part lies at an open end with no node, the script adds a round cap by unioning a point buffer to represent a cul-de-sac. It writes STREET_OID, STREET_NAME, STREET_LEVEL, and PART_ID for each interconnect polygon.

Step 7 — Clean subtraction and merge

All intersection polygons are dissolved with PairwiseDissolve into one surface. A minimal out–in cleanup using PairwiseBuffer expands by CLEAN_BUF and then contracts by the same amount to remove microscopic artifacts. For each interconnect polygon the script computes a geometry difference against this cleaned intersection surface. The cleaned interconnects are merged with the intersection polygons into a staging feature class.

Step 8 — Write the final polygons

The script creates Polygonizer_FinalPolys with the target schema, ensures that location_name, street_levels, and is_intersection exist with the correct types, and appends the staging features into it. If the append creates any suffixed variants of street_levels, those are deleted immediately so that only one street_levels field remains.

Step 9 — Enrich attributes

A one-to-many spatial join from the final polygons to the centerlines gathers unique street names and levels per polygon. A one-to-one spatial join from the final polygons to the true intersection points provides a join count that is converted into the boolean-like is_intersection. The script writes back a sorted, comma-joined name string limited to 250 characters, a sorted, comma-joined level string, and the is_intersection value for every polygon.

Performance choices

Pairwise tools are used where available because they scale better and support parallel processing. The near table is persisted in the geodatabase rather than held in memory because that avoids transient schema locks and intermittent cursor failures. FLAT end caps are used on leg and interconnect buffers to match the intended cartography. The environment’s parallel processing factor is set to one hundred percent to utilize available cores for supported tools.

Resumability and force flags

Every stage is guarded by a force flag so the pipeline can be resumed or partially rebuilt. A full rebuild sets all force flags to True. Attribute-only refreshes set all geometry stages to False and FORCE_ENRICH_ATTRIBUTES to True, leaving geometry intact while recomputing names, levels, and the intersection flag. You can selectively toggle flags to rebuild only steps six through eight when you need to change interconnect segmentation while preserving earlier outputs.

Troubleshooting

If PairwiseDissolve appears missing, call it from the Analysis toolbox as arcpy.analysis.PairwiseDissolve instead of arcpy.management. If a cursor cannot open the near table, the script already logs geoprocessing messages from GenerateNearTable; inspect those messages and verify that _pg_tmp_near was created in the geodatabase. If MakeFeatureLayer fails, note that this script avoids passing field_info to keep joins robust; ensure the centerline feature class is accessible and not exclusively locked. If you ever observe duplicate street_levels fields, confirm that the final writer ran to completion; it creates the field once and deletes any suffixed duplicates after appending.

Design notes and constraints

Leg allocation prevents overlap by capping the leg length at half of the nearest inter-node gap minus a small epsilon, half of the street length, and the configured maximum leg length. Interconnect segmentation yields consistent part sizes per level and adds rounded caps only when a center segment terminates at an open end without a true intersection. Cleaning is applied solely to the dissolved intersection surface with a very small out–in buffer so the subtraction from interconnects is crisp with negligible drift. Temporary datasets are clearly prefixed with _pg_tmp_ and can be safely deleted between runs.

Code reference

This README documents the script titled “ArcPy Polygonizer — SPEED-OPTIMIZED + robust near-table (in-GDB), pairwise tools, FLAT-ended buffers (match baseline geometry), single street_levels” and describes the exact steps implemented in that file. Place the script in your repository (for example as polygonizer.py) and edit the configuration block before running in ArcGIS Pro.
