# Project Context

Last reviewed: 2026-10-27

## Purpose

The project develops an IRIT-oriented iterative optimization workflow for
heterogeneous lattice structures. It generates lattice geometry inside one or
more trivariate macro bodies, exports a finite-element model, runs a structural
solver, maps solver results back to lattice design variables, and updates local
thicknesses over repeated iterations.

## Principal files

- `test.c`: optimization driver, geometry generation, solver integration,
  result mapping, design updates, and utility code.
- `opt_config.h`: public configuration structures and `IritOptRun()` entry
  point declaration.
- `opt_config.txt`: editable run, material, optimization, and path settings.

The `t1` directory contains older comparison copies. It is not a substitute
for version control, and its files must not be assumed to be authoritative
without checking their timestamps and contents.

## Runtime workflow

1. Parse and validate the unified configuration.
2. Resolve default paths relative to `root_folder`.
3. Load every trivariate macro body from the configured ITD file.
4. Assign per-body lattice counts and non-overlapping global tile offsets.
5. Allocate the persistent global tile/subtile design structure.
6. Generate each body's microstructure and combine the results for export.
7. Record the exported geometry-patch-to-tile/subtile ownership map.
8. Export ITD and INP geometry; initialize FE topology on the first export and
   verify the same topology on later exports while updating node coordinates.
9. Run ANSYS or CalculiX.
10. Read element SED/volume results through persistent element-to-subtile
    associations.
11. Update thickness design variables and match neighbours on the combined
    global grid, including interfaces on the selected U, V, or W stack axis.
12. Record results and repeat until convergence or the iteration limit.

## Terminology

### Body

One trivariate macro-shape map loaded from the ITD input. Each body can have
its own U, V, and W tile counts.

### Tile and subtile

A tile is one lattice cell. A subtile is a generated component or design
variable within that tile, such as one face arm or diagonal rod.

### Mesh ownership

Mesh ownership is bookkeeping, not filesystem ownership and not an additional
mesh. It maps an exported geometry patch—and therefore its ordered finite
elements—back to the global tile and subtile that generated it. The optimizer
uses this map to accumulate solver quantities into the correct design variable.

### Continuity

Three different continuity requirements must remain distinct:

- Design continuity: paired interface members use compatible thicknesses.
- Geometric continuity: generated solids meet without gaps or overlaps.
- Analysis continuity: the finite-element mesh is conforming or otherwise
  coupled across the interface.

Matching thickness values alone does not prove geometric or analysis
continuity.

## Current multi-body limitation

Status update, 2026-10-27: the active imported cleanup parses multiple U/V/W
count values but executes only the first set. The selected-axis multi-body
implementation described below remains available in Git history through commit
`4ceb75b`, but it is not active in the current source.

Status update, 2026-10-28: superseded. The active implementation resolves all
body count sets into one logical lattice and calls IRIT separately for every
macro trivariate. Optimization data remains body-independent. Counts must
match on the two non-stack axes. This restores logical cross-interface
neighbour matching but still does not establish geometric or FE continuity.

Bodies are placed consecutively on one configured logical U, V, or W stack
axis; counts on the other axes must match. Design-thickness matching crosses
those logical interfaces. There is still no verification of compatible
parameter orientation, geometry weld, node merge, or solver coupling generated
by `test.c`. Geometric and analysis continuity are therefore not guaranteed.

## Required future cross-body interface workflow

1. Identify touching body-face pairs.
2. Verify positional coincidence of the macro-body interfaces.
3. Determine parameter orientation, including reversed U/V directions.
4. Match compatible interface grids or define a nonmatching correspondence.
5. Compare connection points in world coordinates.
6. Build persistent cross-body adjacency records.
7. Share or constrain paired interface thickness variables.
8. Join geometry and create a conforming mesh, or define an explicit solver
   coupling when a conforming mesh is not possible.
9. Validate positional continuity, requested tangent continuity, watertightness,
   and solver connectivity.

## Verification status

- The source was inspected to establish the current workflow and limitations.
- The October 23 history was reconstructed from timestamped local copies
  because no Git repository existed at the time.
- No new multi-body interface-continuity implementation or runtime validation
  has been completed as of this review.

