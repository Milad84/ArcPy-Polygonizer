# ArcPy Polygonizer (ArcGIS Pro / ArcPy 3.x)

Generate clean, attribute-rich street area polygons from a centerline network — fast, reproducible, and resumable.
This implementation utilizes a series of tools already available in ArcPro, including but not limited to pairwise geoprocessing tools, a robust in-GDB near table, and FLAT-ended line buffers, to name a few. The final layer can inherit values from the underlying network of streets. These Polygons are used to track the distribution of crashes in a segmented fashion to identify places with a higher frequency of fatal and serious injury crashes.

<img width="829" height="697" alt="image" src="https://github.com/user-attachments/assets/c1c2187a-cee5-4eb7-91c0-1495936f0394" />

# Background
* This script is inspired by the work done by my colleague Frank Hereford (https://github.com/frankhereford) at the Data and Technology Services of the Austin transportation and Public Works Department. The QGIS plug-in is available at https://github.com/frankhereford/qgis-polygonizer/tree/main. 
* The Vision Zero team, as part of the Austin Transportation and Public Works Department, operates a database of crashes in and around the Austin city limits. A point position tracks these crashes, and to conduct an analysis of the city's roadways about crash location and crash severity, a set of polygons was created based on the street network. These polygons are of approximately equal area and have been formed in such a way that intersections are built into a single shape, and interconnecting roads between intersections are divided into equal segments. By aggregating crashes into these polygons, comparisons become possible between intersections and road segments.

# Goal
* Having the possibility of creating polygons inside the ArcPro environment.
  
# Highlights

* Resumable: each stage can be skipped unless forced to rebuild.
*  Faster: Pairwise* tools + parallelProcessingFactor="100%".
*  Robust intersections: GenerateNearTable written to the GDB (not in-memory).
*  FLAT buffer ends for legs and interconnects (matches baseline cartography).
*  Clean attributes: location_name, single street_levels, is_intersection.

# Requirements

*  ArcGIS Pro 3.x with ArcPy.
*  Feet-based projected CRS (default: EPSG 2277 – NAD 1983 StatePlane Texas Central (US Feet)).

# Tools used:

* PairwiseIntersect 
* PairwiseDissolve
* PairwiseBuffer
* GenerateNearTable
* SpatialJoin
* Merge
* Buffer
* Append
* Dissolve

If a pairwise tool isn’t available in your environment, swap to the classic equivalent (potentially slower).

# Input Assumptions

* CTN_FC is a polyline centerline layer with:
* STREET_LEVEL (integer or convertible text).
* A name field: prefers FULL_STREET_NAME, otherwise STREET_NAME.
* Centerlines are reasonably snapped (default matching tolerance TOUCH_TOL = 5 ft).

# Outputs
All datasets are written into the configured file geodatabase.
| Dataset                          | Type    | Purpose                                                                                                     |
| -------------------------------- | ------- | ----------------------------------------------------------------------------------------------------------- |
| `Polygonizer_BuiltIntersections` | Point   | Raw self line–line candidates from PairwiseIntersect written as points.                                     |
| `Polygonizer_TrueIntersections`  | Point   | Points that have at least two incident streets within the tolerance, derived from a near table.             |
| `Polygonizer_IntersectionPolys`  | Polygon | Union of a circular node core and FLAT-ended leg buffers sized by street level.                             |
| `Polygonizer_InterconnectPolys`  | Polygon | Segmented FLAT-ended buffers along streets between intersections, with optional round caps for cul-de-sacs. |
| `Polygonizer_FinalPolys`         | Polygon | Merged and cleaned polygons with attributes `location_name`, `street_levels`, and `is_intersection`.        |


# Final attribute schema

* location_name (TEXT, 250): comma-joined unique street names intersecting the polygon.
* street_levels (TEXT, 50): sorted, comma-joined unique levels.
* is_intersection (SHORT): 1 if the polygon intersects ≥1 true-intersection point; else 0.
* The writer guarantees exactly one street_levels field on the final output.

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

* If matching seems inconsistent, try TOUCH_TOL = 8–10.
* Customize widths & segment lengths per STREET_LEVEL.
* CLEAN_BUF = 0.1 ft gently removes micro-slivers after dissolve.

# Quick Start

* Open ArcGIS Pro 3.x → Python Window / Notebook.
* Paste the script (or run as a Script Tool).
* Set GDB, CTN_FC, and (if needed) SR.
* Run with defaults for a full rebuild.

## To only refresh attributes later:

```python
FORCE_RAW_INTERSECTION_PTS   = False
FORCE_TRUE_INTERSECTION_PTS  = False
FORCE_INTERSECTION_POLYS     = False
FORCE_INTERCONNECT_POLYS     = False
FORCE_TRIM_INTERCONNECTS     = False
FORCE_FINAL_POLYS            = False
FORCE_ENRICH_ATTRIBUTES      = True
```
# How It Works (Pipeline)

## 1) Load Streets

Reads OID@, SHAPE@, STREET_LEVEL, and a name field (FULL_STREET_NAME or STREET_NAME), caching length and level per feature.

## 2) Raw Intersection Candidates

PairwiseIntersect([CTN_FC], …, "POINT") finds self intersections.
Deduplicate with DeleteIdentical. Write points to Polygonizer_BuiltIntersections.

## 3) True Intersection Points

GenerateNearTable (raw points ↔ CTN within TOUCH_TOL) is written to the GDB and verified.
Count incident streets by IN_FID; keep points with count ≥ 2 → Polygonizer_TrueIntersections.

## 4) Map Nodes to Streets & Leg Lengths

Using the near table, associate each node to incident streets and compute non-overlapping leg lengths:
```python 
leg_len = min(
  GOAL_LEG_LEN,
  0.5 * nearest_gap - HALF_GAP_EPS,
  0.5 * street_length
)
```

extract_leg() keeps the half that runs from the node toward the nearer end.

## 5) Intersection Polygons (FLAT Legs)

For each node:

* Build a circular core using the max half-width among incident street levels.
* Buffer each incident leg with Buffer_analysis(..., line_end_type="FLAT") at that level’s half-width.
*Union the core with all leg buffers and write X_ID, LEGS, location_name, street_levels.

## 6) Interconnect Polygons (FLAT Segments)

### For each street:

* Detect nodes at ends; carve off end legs.
* Split the remaining center segment into n_parts using TARGET_LEN_BY_LEVEL[level].
* Buffer each segment with FLAT end caps.
* If a part touches an open end (no node), add a round cap by unioning a point buffer (cul-de-sac behavior).
* Write STREET_OID, STREET_NAME, STREET_LEVEL, PART_ID.

## 7) Clean Subtraction & Merge

* PairwiseDissolve all intersection polygons to a single surface.
* Out–in PairwiseBuffer (+CLEAN_BUF then -CLEAN_BUF) removes micro artifacts.
* For each interconnect polygon, subtract the cleaned surface via geometry.difference().
* Merge intersections + trimmed interconnects into a staging FC.

## 8) Final Layer Writing

* Create Polygonizer_FinalPolys with the target schema and append the staging FC.
* Defensively delete any suffixed duplicates of street_levels (rare).

## 9) Attribute Enrichment

* SpatialJoin (1:M) FinalPolys ← CTN: collect unique names & levels per polygon.
* SpatialJoin (1:1) FinalPolys ← TrueIntersections: compute is_intersection.
  
## Update:

* location_name = sorted, comma-joined names (≤ 250 chars)
* street_levels = sorted, comma-joined unique levels
* is_intersection = 0/1

# Performance Choices

Pairwise tools where possible: PairwiseIntersect, PairwiseDissolve, PairwiseBuffer.
Near table is persisted in the GDB (stable, lock-resilient).
FLAT buffer caps for legs/interconnects match the intended cartography.
arcpy.env.parallelProcessingFactor = "100%" enables multi-core where supported.

# Resumability (Force Flags)

Common modes:
```python
# Full rebuild (first run)
FORCE_RAW_INTERSECTION_PTS   = True
FORCE_TRUE_INTERSECTION_PTS  = True
FORCE_INTERSECTION_POLYS     = True
FORCE_INTERCONNECT_POLYS     = True
FORCE_TRIM_INTERCONNECTS     = True
FORCE_FINAL_POLYS            = True
FORCE_ENRICH_ATTRIBUTES      = True
```
```python
# Geometry frozen, refresh attributes only
FORCE_RAW_INTERSECTION_PTS   = False
FORCE_TRUE_INTERSECTION_PTS  = False
FORCE_INTERSECTION_POLYS     = False
FORCE_INTERCONNECT_POLYS     = False
FORCE_TRIM_INTERCONNECTS     = False
FORCE_FINAL_POLYS            = False
FORCE_ENRICH_ATTRIBUTES      = True
```
# Troubleshooting

* PairwiseDissolve not found → Use arcpy.analysis.PairwiseDissolve (Analysis toolbox), not arcpy.management.
* Near table cursor error (cannot open ... _pg_tmp_near) → The script writes the near table to the GDB and checks existence; inspect geoprocessing messages printed right after GenerateNearTable.
* MakeFeatureLayer crash → No field_info is passed; a simple view avoids brittle mappings.
* Duplicate street_levels → Final writer creates the field once, appends, then deletes any suffixed duplicates defensively.
* Unexpected shapes → Ensure US Feet and arcpy.env.outputCoordinateSystem = SR. Legs & interconnects use FLAT caps; cul-de-sacs are rounded by design.

# Design Notes & Constraints

### Leg allocator prevents overlap by capping at:
* ½ of the nearest inter-node gap minus HALF_GAP_EPS
* ½ of the street length
### GOAL_LEG_LEN
* Interconnect segmentation yields consistent part sizes by STREET_LEVEL.
* Cleaning uses a tiny out–in buffer on the intersection dissolve surface only, minimizing drift.
* Temp datasets are prefixed _pg_tmp_* and can be deleted safely between runs.

