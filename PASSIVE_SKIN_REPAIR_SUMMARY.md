# Passive Skin Ownership Repair

Implemented and committed on a new branch:

- Branch: `feature/passive-skin-element-ownership`
- Commit: `be61317` — `Separate passive skin from lattice optimization`

The repair now:

- Keeps skin and lattice elements in one connected INP simulation.
- Tags geometry explicitly as lattice or passive skin.
- Binds only lattice elements to tile/subtile design variables.
- Excludes skin SED, volume, compliance, sensitivity, and thickness from design updates.
- Keeps skin stiffness active in the simulation.
- Accounts for fixed skin volume in the MMA volume target.
- Validates ownership remains unchanged between iterations.

Verification:

- `test.c` and the parser compiled successfully with Visual C++ 2022 x64 `/W4`.
- `git diff --check` passed.
- A complete skinned solver/optimization run has not yet been performed.
- The checked-in Visual Studio project could not run because its referenced
  `..\IRIT_SM\GlobalRules.props` is missing at this repository location.

The branch is committed locally but not pushed.

## Change size

Commit `be61317` added 422 lines and removed 51 lines across six files, for a
net increase of 371 lines.

In `test.c` specifically, it added 340 lines and removed 50 lines, for a net
increase of 290 lines. The remaining 82 added lines and one removed line are
project-history documentation.

## Functions added

Six new function definitions were added to `test.c`:

1. `ComputeOwnedObjVol()` — computes geometry volume for one explicit owner
   type, allowing lattice and passive-skin volumes to be separated.
2. `IsFullFaceSkinTV()` — identifies a full-face skin trivariate in the
   undeformed unit tile.
3. `TagTileMeshOwnership()` — assigns owner type, tile index, and subtile index
   metadata to each synthesized tile trivariate.
4. `CountMeshPatchOwners()` — counts trivariate patches in final INP export
   traversal order.
5. `FillMeshPatchOwners()` — builds the explicit patch-owner table from the
   final scaled geometry.
6. `BuildMeshPatchOwners()` — allocates the table and validates persistent
   lattice/skin ownership on later iterations.

The existing `BindElementsToSubTiles()` function was substantially modified.
It now consumes every exported FEM patch, binds lattice patches to subtiles,
marks skin elements as passive, and verifies that every expected lattice patch
was bound exactly once.
