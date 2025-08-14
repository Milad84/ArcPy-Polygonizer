# ArcPy-Polygonizer
This script creates polygons across a selected network of streets keeping an eye on the intersection and by managing the size and shape of the polygons over the network of streets.
<img width="823" height="626" alt="image" src="https://github.com/user-attachments/assets/492d7507-c232-4257-9ef0-dfc506c6d133" />
ArcPy Polygonizer (ArcGIS Pro / ArcPy 3.x)

Generate clean, attribute-rich street area polygons from a centerline network — fast, reproducible, and resumable.
This implementation uses pairwise geoprocessing tools, a robust in-GDB near table, and FLAT-ended line buffers. The final layer contains a single street_levels field and an is_intersection flag.

Table of Contents

Highlights

Requirements

Input Assumptions

Outputs

Configuration

Quick Start

How It Works (Pipeline)

1) Load Streets

2) Raw Intersection Candidates

3) True Intersection Points

4) Map Nodes to Streets & Leg Lengths

5) Intersection Polygons (FLAT Legs)

6) Interconnect Polygons (FLAT Segments)

7) Clean Subtraction & Merge

8) Final Layer Writing

9) Attribute Enrichment

Performance Choices

Resumability (Force Flags)

Troubleshooting

Design Notes & Constraints

License

Full Script (Reference)

Highlights

🔁 Resumable: each stage can be skipped unless forced to rebuild.

🚀 Fast: uses Pairwise* tools where available and runs with parallelProcessingFactor="100%".

🧭 Robust topology: true intersections derived via GenerateNearTable stored in the GDB (not in_memory).

🧱 FLAT caps on line buffers (legs and interconnects) to match baseline cartography.

🧾 Clean attributes: location_name, single street_levels, and is_intersection.

Requirements

ArcGIS Pro 3.x with ArcPy available.

Projected, feet-based coordinate system (default: EPSG 2277).

Access to Analysis tools used here:

PairwiseIntersect, PairwiseDissolve, PairwiseBuffer

GenerateNearTable

SpatialJoin, Merge, Dissolve, Buffer, Append

If a pairwise tool isn’t licensed, you can swap to the classic equivalent at a potential performance cost.

Input Assumptions

CTN_FC is a polyline centerline layer with:

STREET_LEVEL (integer or convertible text) used for widths & segmentation lengths.

A name field: prefers FULL_STREET_NAME, otherwise STREET_NAME.

Centerlines are reasonably snapped (default matching tolerance TOUCH_TOL = 5 ft).

Outputs

All outputs are written to the configured geodatabase.

Dataset	Type	Description
Polygonizer_BuiltIntersections	Point	Raw self-intersection candidates from PairwiseIntersect
Polygonizer_TrueIntersections	Point	Points with ≥2 incident streets within TOUCH_TOL
Polygonizer_IntersectionPolys	Polygon	Node cores + FLAT-ended leg buffers (unioned)
Polygonizer_InterconnectPolys	Polygon	Segmented FLAT-ended buffers between intersections
Polygonizer_FinalPolys	Polygon	Final merged surface with attributes

Final fields

location_name (TEXT, 250): comma-joined unique street names intersecting the polygon.

street_levels (TEXT, 50): sorted, comma-joined unique levels intersecting the polygon.

is_intersection (SHORT): 1 if polygon intersects ≥1 true-intersection point; otherwise 0.

The writer guarantees exactly one street_levels column.
