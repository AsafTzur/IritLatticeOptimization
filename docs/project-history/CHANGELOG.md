# Changelog

This log records material project changes. Entries made before local source
control is initialized may be reconstructed from local files and timestamps.

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

