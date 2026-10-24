# Changelog

This log records material project changes. Entries made before local source
control is initialized may be reconstructed from local files and timestamps.

## 2026-10-24 — Selectable U/V/W multi-body stack axis

Purpose: allow the user to join continuous bodies along any one logical axis,
while rejecting layouts that vary on more than the selected axis.

Changes:

- Added `body_stack_axis = u|v|w` to the public run configuration; `-1`
  remains available programmatically for unambiguous inference.
- Replaced the W-only body origin with a three-component logical grid origin.
- Required counts on both non-stack axes to be identical across all bodies;
  differences outside the selected axis are reported as configuration errors.
- Changed CSV writing, ASAF updating, IRIT thickness lookup, geometry tagging,
  and tile ownership naming to use canonical global flat-grid indices.
- Kept `NeighbourMatchingB_V2()` body-independent over the global grid.

Verification:

- Visual C++ 2022 compilation succeeded with only existing external-header
  warnings.
- Synthetic U-, V-, and W-stacked two-body layouts each produced complete,
  unique global index coverage; a non-selected-axis difference was detected.
- `git diff --check` reported no whitespace errors.
- No complete multi-body optimization run was performed.

## 2026-10-24 — Global-W multi-body neighbour matching

Purpose: treat continuous macro bodies as consecutive W ranges of one logical
3D lattice so thickness matching uses the same flattened structure everywhere.

Changes:

- Added global U/V/W lattice counts to the workflow state.
- Required every continuous body to use identical U and V tile counts; W is
  allowed to differ and global W is their sum.
- Recorded each body's cumulative logical W origin and verified that its
  existing flat `TileOffset` equals the corresponding global-grid offset.
- Replaced body-local `NeighbourMatchingB()` with
  `NeighbourMatchingB_V2()`, which traverses and sorts one flattened global
  `(U,V,W)` grid and therefore averages members across consecutive body W
  interfaces.

Verification:

- Visual C++ 2022 compilation succeeded. Only existing IRIT/Windows header
  warnings were emitted.
- A synthetic two-body layout (`U=2`, `V=1`, body W counts `1,2`) verified the
  cumulative offsets, six-tile global flattening, and shared logical coordinate
  of the first body's W+ and second body's W- interface members.
- `git diff --check` reported no whitespace errors.
- No complete multi-body optimization run was performed.

## 2026-10-24 — Persistent finite-element topology

Purpose: retain the optimization model's node, element, and element-to-subtile
topology across iterations while allowing geometry coordinates, solver fields,
and thickness design variables to change.

Changes:

- Extended `OptElementStruct` with the eight persistent C3D8 node IDs.
- Changed repeated INP parsing to update node coordinates in place and verify
  that node IDs, element IDs, element connectivity, and counts remain stable.
- Removed per-iteration destruction of the node and element arrays.
- Bound element pointers to subtiles once and preserved those associations
  while clearing only solver-derived tile/subtile values.
- Preserved the first exported mesh-ownership sequence and verified every later
  geometry generation against it.
- Created Git baseline commit `e94f17c` before the refactor so the change has a
  reviewable before/after history.

Verification:

- Compiled `test.c` successfully with Visual C++ 2022 using the IRIT Windows
  defines and local NLopt headers. Only existing IRIT/Windows header warnings
  were emitted.
- Scanned `C:\IRIT\ASAF\RUN\TEMP_SET\optimized_geo.inp`: 340,948 sequential
  nodes and 229,376 sequential C3D8 elements were found; every element row had
  one element ID and eight node IDs.
- A complete multi-iteration optimization run was not performed.

## 2026-10-24 — Local Git history initialized

Purpose: provide exact, recoverable, line-by-line history without depending on
Visual Studio integration or a GitHub account.

Changes:

- Initialized a standalone Git repository at `C:\irit\irit\test`.
- Added `.gitignore` rules for compiler, Visual Studio, CalculiX, ANSYS, and
  temporary generated files.
- Kept source files, reference examples, historical comparison copies,
  configuration, and project documentation eligible for tracking.
- Used repository-local author identity `Asaf Tzur <asaf.tzur@local>` for the
  initial offline snapshot; this placeholder should be changed before pushing
  to a hosted service if a different public commit identity is desired.

Verification:

- Git repository root and tracked-file status checked from the command line.
- No GitHub remote was created; the repository remains only on this computer.

## 2026-10-24 — Durable project documentation established

Purpose: preserve technical context, change history, design decisions, and
thesis-relevant evidence on the local computer for use in future work.

Changes:

- Added root `AGENTS.md` with mandatory documentation and recording rules.
- Added the `docs/project-history` documentation set.
- Clarified the meaning of mesh ownership and documented that it was existing
  behavior, not a newly invented system.
- Documented the current absence of cross-body continuity enforcement.
- Corrected `test_c_repairs_2026-10-23.html` to describe the October 23
  multi-body change as consolidation/refactoring.

Verification:

- Compared the current `test.c` with timestamped copies under `t1`.
- Inspected body loading, tile offsets, geometry collection, mesh ownership,
  solver-result mapping, and neighbour matching in the current source.

## 2026-10-23 — Configuration and workflow refactor

Status: reconstructed from `t1/test - Copy.c`, the current files, and local
timestamps; there was no Git history proving the exact edit sequence.

Changes:

- Moved public configuration structures into `opt_config.h`.
- Added the callable `IritOptRun()` entry point.
- Added default-path resolution from `root_folder` while retaining selected
  explicit path overrides.
- Added configuration copying and cleanup for per-body count arrays.
- Added validation for the macro-shape `.itd` filename and skin thickness.
- Added workflow finish-reason recording and distinguished iteration failure
  from ordinary continuation.
- Consolidated the existing multi-body and mesh-ownership implementation into
  `test.c`; this was primarily relocation/refactoring, not a new ownership
  concept.
- Replaced generic auxiliary comments with more precise descriptions of the
  associated code.

Verification evidence:

- Current `test.c` modification time: 2026-10-23 23:54:46 local time.
- Comparison against `t1/test - Copy.c` showed 735 insertions and 210 deletions.
- No full runtime multi-body continuity test is recorded.

## 2026-10-21 — Skin-boundary repair

Purpose: apply a configurable exterior skin without out-of-bounds boundary
array access, unintended internal sheets, or the known approximate-FFD failure.

Changes and evidence are recorded in `SKIN_BOUNDARY_CHANGELOG.md`.

Important results:

- Boundary thickness uses exactly four corner values.
- Skin is applied only to selected exterior V faces for lattice type 0.
- Exact local deformation is used because approximate/global FFD produced
  negative-Jacobian trivariates on the curved Winglet map.
- Compilation was reported successful; the record says a freshly rebuilt
  executable was still required for runtime verification.

