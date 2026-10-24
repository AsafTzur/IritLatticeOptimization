# Thesis Notes

This file collects technically supportable material for a future thesis. Draft
claims must remain tied to code, experiment logs, solver output, or figures.

## System contribution

The current system integrates parametric lattice generation, heterogeneous
local sizing, finite-element analysis, and iterative design updates. Its
multi-body extension assigns each macro body an independent trivariate map and
lattice resolution while maintaining global tile indices for combined analysis.

## Mesh-to-design result transfer

Generated geometry is tagged with global tile and subtile identifiers. During
export, an ordered ownership table is constructed. After analysis, element
results are assigned back through this table and aggregated into the persistent
tile/subtile design structure. This mechanism connects solver-space quantities
such as strain-energy density and element volume to parametric design variables.

Potential thesis figure:

`macro body -> tile -> subtile -> exported patch -> FE elements -> solver result -> thickness update`

### Persistent-topology invariant

The optimization now treats node IDs, element IDs, C3D8 connectivity, exported
patch order, and element-to-subtile ownership as persistent topology. Later INP
files update node coordinates and are rejected if their topology differs. This
is appropriate only while lattice layout, mesh sampling, and export ordering
remain fixed. The implementation was compile-checked and its C3D8 parsing
assumption was checked against one existing 229,376-element INP; a complete
multi-iteration runtime experiment is still required.

## Multi-body method currently implemented

- Load multiple trivariate macro bodies from one ITD input.
- Validate that each body has configured U/V/W lattice counts.
- Assign a global tile-offset range to every body.
- Generate lattice geometry independently in each body's parameterization.
- Combine generated body geometry into one export.
- Preserve global body/tile/subtile identity during result mapping.

## Known limitation: body-interface continuity

The revised neighbour-matching algorithm represents continuous bodies as
consecutive W ranges of one integer lattice when their U/V counts match. It can
therefore equalize thicknesses across those logical W interfaces. Compatible
parameter orientation and physical coincidence of the macro-body interfaces
are assumed rather than geometrically verified.

This leaves three unproven properties:

1. Geometric contact or watertightness between generated body lattices.
2. Tangent or derivative continuity at the interface.
3. Finite-element load transfer through a conforming or coupled interface.

This limitation is a suitable thesis motivation for validating the logical
correspondence in world space and adding a geometric/solver coupling method.

## Proposed interface-continuity research task

Develop a correspondence method that maps boundary lattice connection points
from each body's parametric domain into common world coordinates. Pair points
within a scale-aware tolerance, account for reversed parameter directions, and
represent matched members with shared or constrained thickness variables.
Geometry and analysis continuity must then be established by node merging,
watertight construction, or an explicit solver coupling.

## Skin-boundary observation

For the curved Winglet case, approximate/global FFD deformation produced
negative-Jacobian skin trivariates. Exact local deformation was selected for
the repaired path. This is an engineering observation, not yet a general proof;
future experiments should record geometry, deformation settings, Jacobian
statistics, runtime, and failure counts.

## Evidence still needed

- Reproducible build command and compiler/linker configuration.
- Baseline single-body and multi-body input datasets.
- Solver versions and complete run configurations.
- Interface gap and overlap measurements.
- Mesh-connectivity checks at body interfaces.
- Convergence plots and objective/volume histories.
- Comparison of exact deformation and approximate/global FFD paths.
- Screenshots or rendered figures tied to named experiment runs.

## Writing rule

Do not convert an implementation intention into a thesis claim. Use wording
such as “implemented,” “verified,” “observed,” or “proposed” only when the
corresponding evidence level is recorded.

