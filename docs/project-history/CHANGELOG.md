# Changelog

This log records material project changes. Entries made before local source
control is initialized may be reconstructed from local files and timestamps.

## 2026-11-07 — Simplify configuration cleanup and make field growth atomic

Purpose: remove duplicated ownership cleanup and prevent partially resized
parallel scalar-field arrays.

Changes:

- Replaced the duplicated lattice-count-array cleanup loops with one
  `OptFreeConfigArrays()` helper used by caller and workflow configurations.
- Changed `ScalarFieldEnsureCapacity()` to allocate both replacement arrays
  before changing the field, copy existing values only after both allocations
  succeed, and then release the old arrays.

Verification:

- Visual C++ 2022 compiled the `test.c` translation unit successfully through
  the Release x64 `ClCompile` target.
- `git diff --check -- test.c` reported no whitespace errors.

## 2026-11-07 — Consolidate all configuration into one structure

Purpose: prepare the optimization entry point for a future IRIT-script wrapper
that supplies every setting through one C configuration object.

Changes:

- Replaced the nested `OptRunInputStruct` and `OptConfigInputStruct` types with
  one flat `OptConfigStruct` in `opt_config.h`.
- Moved `UseTargetVolume` and `GeometryOnly` into `OptConfigStruct`, removing
  the former global flag and separate `IritOptRun()` mode argument.
- Changed workflow state to own one deep-copied `OptConfigStruct`, including
  independent copies of the dynamic per-body lattice-count arrays.
- Updated parsing, defaults, path resolution, validation, and all workflow
  field access to use the same configuration object.

Verification:

- Visual C++ 2022 compiled the `test.c` translation unit successfully through
  the `ClCompile` MSBuild target in Release x64 configuration.
- `git diff --check -- test.c opt_config.h` reported no whitespace errors.

## 2026-11-07 — Repair skinned-tile degree unification

Purpose: remove linker failures caused by calls from `test.c` to private
`static` helpers in IRIT's executable layer.

Changes:

- Reimplemented the local uniform-degree adapter with the public
  `IritTrivTVsUniformCommonDegree` API.
- Preserved tile and trivariate attributes on the newly allocated geometry.
- Replaced and freed the original tile after successful degree unification.
- Kept the call site behind the user's `if (FALSE)` development toggle, so the
  behavior remains disabled until explicitly enabled.

Verification:

- Source search confirmed that `test.c` no longer references
  `IritMicroTileGetMaxOrder` or `IritMicroTileRaiseUniformDegree`.
- A full Visual Studio build was attempted but could not complete because
  MSBuild's `FileTracker` received an access-denied error in the environment.

## 2026-11-04 — Reduce the tracked repository to active project files

Purpose: keep the synchronized GitHub working tree focused on the active
optimization source, configuration, and durable project memory.

Changes:

- Added `No_Git/` to `.gitignore` for local comparison and archive material.
- Removed old IRIT test bundles, historical source copies, the imported cleanup
  package, and the obsolete HTML repair summary from the tracked tree.
- Retained `AGENTS.md` and `docs/project-history/` because they provide the
  cross-chat instructions and technical history required by the active
  project workflow.

Verification:

- The removal commit was scoped explicitly; active source/configuration files
  and unrelated uncommitted edits were not staged.
- Removed files remain recoverable from Git history and local copies under
  `No_Git/` where available.

## 2026-11-04 — Import ANSYS winglet INP support-face conversion

Purpose: retain the newer ANSYS-compatible INP generation from `test1.c`
without restoring its cancelled passive-skin optimization mapping.

Changes:

- Added compile-time selection of the supported macro body and its parametric
  V-min or V-max boundary.
- Exported `SUPPORT` as an element-face surface and converted its C3D8 faces
  through the persistent FEM connectivity into a unique `SUPPORT_NODES` set.
- Applied the ANSYS `*BOUNDARY` card to `SUPPORT_NODES` while preserving the
  `PRESSURE` element surface used by `*DSLOAD`.
- Kept the strict `test.c` element-to-subtile binder. The unverified
  `SET_n`-range inference and passive-skin branch from `test1.c` were not
  imported.

Verification:

- Source search confirmed that `SUPPORT_NODES` and the selectable support-side
  path are present in `test.c`.
- Source search confirmed that `ReadGeneratedTileElementRanges`, the
  `Skin-aware topology mapping` diagnostics, and the passive-skin branch remain
  absent from `test.c`.
- `git diff --check -- test.c` reported no whitespace errors after cleanup.
- Compilation was not run because `cl.exe` was not available on the active
  command-runner PATH.

## 2026-11-02 — Remove unsuccessful skin-element exclusion attempt

Purpose: return the optimizer's FEM ownership mapping to a clean baseline
before redesigning how passive skin elements are distinguished from lattice
elements.

Changes:

- Removed the skin-aware `SET_n` range parser from `test.c`.
- Removed the branch that assumed each generated tile set stored all lattice
  elements first and all skin elements afterward.
- Restored strict one-to-one FEM binding: the element count must equal
  `tile count * subtile count * mesh sample` before elements are assigned to
  lattice design variables.
- Preserved skin geometry generation, INP export, material assignment, and
  solver participation. No replacement skin-exclusion method was introduced.

Behavioral consequence:

- A skinned model with additional FEM elements now stops with the existing
  mesh/order mismatch instead of silently classifying elements by an
  unverified export-order assumption.

Verification:

- Confirmed by source search that the `ReadGeneratedTileElementRanges` helper,
  skin-aware mapping messages, passive-skin branch, and its cleanup label are
  absent from `test.c`.
- `git diff --check` reached the working tree and reported only two existing
  trailing-whitespace lines in unrelated portions of the broader uncommitted
  `test.c` changes; this rollback introduced no reported whitespace error.
- Compilation was not run because `cl.exe` was not available on the active
  command-runner PATH.

## 2026-10-28 — Configurable microstructure FFD mode

Purpose: allow the run configuration to select the existing approximate/global
FFD microstructure path without recompiling `test.c`.

Changes:

- Added `use_microstructure_ffd = true|false` to `opt_config.txt` and
  `OptRunInputStruct`.
- `GenerateMicroStructures()` now reads the setting instead of using a local
  hard-coded `UseMicroStructFFD` constant.
- Defaulted the setting to `true`, preserving the currently selected behavior.
- Kept `UseFFDTileAggregation` internal and unchanged.

Verification:

- Visual C++ 2022 compiled `test.c` successfully with `/W4`; only the existing
  nameless-union warnings from IRIT's `mvar_lib.h` were emitted.

## 2026-10-28 — Unified-lattice multi-body IRIT generation

Purpose: support macro-shape files containing multiple trivariates without
making the optimization lattice, design update, or neighbour model body-aware.

Changes:

- Load every macro-shape trivariate in ITD order and require one configured
  U/V/W count value per trivariate.
- Validate that counts differ on at most one stack axis; infer that axis when
  possible and accept `body_stack_axis = u|v|w` for explicit or ambiguous
  equal-count layouts.
- Resolve the count lists into one unified logical U/V/W lattice. Existing
  tile/subtile allocation, design updates, CSV data, and neighbour-group
  construction continue to operate on that single lattice.
- Call the IRIT tiler separately for each trivariate. The callback translates
  body-local tile coordinates through a temporary grid origin and reads the
  corresponding thickness from the unified list.
- Combine all generated body objects before the single ITD/INP export.
- Bind the exported FEM elements in per-body generation order back into their
  canonical global tiles, preserving the persistent element-ID mapping.
- Kept the manual four-element skin-boundary assignment that replaced the
  ineffective `memcpy` during interactive validation.

Verification:

- Visual C++ 2022 compiled `test.c` successfully with `/W4`, `__WINNT__`, and
  the configured IRIT/NLopt include paths. Only existing nameless-union
  warnings from `mvar_lib.h` were emitted.
- `git diff --check` reported no whitespace errors.
- Synthetic two-body U-, V-, and W-stack layouts produced complete unique
  global-index coverage (84/84, 72/72, and 54/54 tiles respectively).
- No complete multi-body geometry, solver, or optimization run was performed.

## 2026-10-27 — Import persistent-topology cleanup from external ChatGPT work

Purpose: replace the tracked implementation with the newer cleanup package
prepared on another machine while retaining its source bundle and recording
the behavioral tradeoffs in the base repository.

Changes:

- Replaced `test.c`, `opt_config.h`, and `opt_config.txt` with the versions from
  `persistent_topology_cleanup/`.
- Introduced explicit FEM and lattice model containers, permanent per-subtile
  element-ID ownership, and persistent neighbour groups built once per run.
- Changed INP processing so topology is parsed and bound after the first real
  geometry export and is not reread on later iterations.
- Rejected `remove_thin_tiles=true`, because removing geometry would invalidate
  the fixed finite-element topology.
- Retained the supplied cleanup package as import provenance.
- Intentionally replaced the later multi-body stack-axis implementation: the
  imported version parses multiple U/V/W values but executes only the first
  count set. Cross-body neighbour matching is therefore no longer active.

Verification:

- Pending compilation and Git whitespace checks in this session.

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

Problem observed:

- An early attempt wrote `BPrm->Bndry[4]`, although the boundary array has
  exactly four entries indexed 0 through 3. IRIT creates a face skin only when
  all four corner thicknesses are positive.
- Setting all four values to `0.05` while using approximate/global FFD mapping
  produced many negative-Jacobian trivariates on the curved Winglet map and
  stopped before processing the second body.

Repair recorded at the time:

- `PreProcessTile1FaceParam` accepted
  `const CagdRType SkinThickness[4]` and copied exactly four values into the
  boundary structure.
- Lattice type 0 used four-value uniform-skin and no-skin arrays.
- Skin was assigned only to exterior V-min and V-max faces; interior V faces
  and all U/W faces used zero thickness to avoid internal solid sheets.
- Exact local deformation was selected because the approximate/global FFD
  path inverted skin trivariates on the Winglet deformation map.

Verification recorded at the time:

- The old negative-Jacobian failure was reproduced.
- Compilation of the repaired source was reported successful.
- A freshly rebuilt executable was still required for runtime verification;
  no completed runtime validation was recorded.

Later evolution: the skin thickness and FFD choice subsequently became
configuration-driven, and the implementation changed from the originally
recorded `memcpy` form. The durable invariant remains that there must be no
`Bndry[4]` access and every face boundary uses exactly four corner values.

