# Design Decisions

## D-011 — Expose one flat optimization configuration object

Date: 2026-11-07

Decision: every externally supplied run, lattice, solver, material, path, and
execution-mode setting is stored directly in `OptConfigStruct`; `IritOptRun()`
accepts only a pointer to that object.

Reason: the planned IRIT-script interconnection needs one unambiguous API
payload. Nested run/environment types, a process-global volume-search flag,
and a separate geometry-only argument created multiple configuration channels
and made scripted invocation state-dependent.

Ownership: callers own their configuration and count arrays. The workflow
deep-copies the complete structure and its dynamic count arrays, so runtime
updates and cleanup do not alter caller-owned input.

## D-001 — Store durable project records beside the source

Date: 2026-10-24

Decision: Keep root instructions in `AGENTS.md` and detailed records in
`docs/project-history` under `C:\irit\irit\test`.

Reason: Future Codex work launched from the project directory can discover the
instructions and read the same records as the source. The records remain local
and can later be committed to Git without requiring GitHub.

## D-002 — Treat mesh ownership as result-mapping metadata

Date: 2026-10-24

Decision: Use “mesh ownership” only for the mapping from ordered exported
geometry/finite elements to global tile and subtile indices.

Reason: The name can otherwise be mistaken for mesh creation, memory ownership,
or geometric joining. The mechanism is used when assigning solver results to
the design structure and does not enforce body-interface continuity.

## D-003 — Do not claim that current multi-body support is continuous

Date: 2026-10-24

Decision: Describe the present feature as multi-body loading, indexing,
generation, export, and result mapping—not cross-body geometric continuity.

Reason: `NeighbourMatchingB()` processes each body independently. No current
step creates cross-body adjacency, welds geometry, merges mesh nodes, or adds
solver coupling.

Status update, 2026-10-24: the body-local design-matching portion was
superseded by D-006. The absence of a geometry weld, merged mesh nodes, or
solver coupling remains true.

## D-004 — Keep continuity requirements separate

Date: 2026-10-24

Decision: Evaluate design-variable continuity, geometric continuity, and
finite-element connectivity independently.

Reason: Equal thickness parameters do not prove that solids meet, and solids
that appear to meet do not prove that the solver mesh transfers load across the
interface.

## Open decision — Interface continuity target

Before implementation, choose the required geometric target:

- C0/positional continuity: bodies meet without a gap.
- G1/tangent continuity: tangent planes align visually/geometrically.
- C1 parametric continuity: first derivatives also match under a compatible
  parameterization.

The finite-element connection strategy must also be selected: conforming merged
nodes versus an explicit tied/contact/coupling formulation.

## D-005 — Persist and verify finite-element topology

Date: 2026-10-24

Decision: Allocate node and element containers on the first generated INP,
retain element connectivity and element-to-subtile pointers for the entire
optimization, and treat later INP files as coordinate updates that must pass a
strict topology comparison.

Reason: Thickness changes deform the generated node coordinates but, with a
fixed lattice layout and mesh sampling, should not change node IDs, element
IDs, C3D8 connectivity, patch order, or ownership. Preserving this invariant
removes redundant topology reconstruction and prevents silent result-mapping
errors. A detected topology change aborts processing rather than leaving stale
pointers active.

## D-006 — Model continuous bodies as consecutive global W ranges

Date: 2026-10-24

Decision: Require identical U/V lattice counts for all continuous bodies,
allow W counts to differ, and interpret body order as consecutive ranges in a
single global W direction. Run neighbour matching on the resulting flattened
global grid without body boundaries.

Reason: With U-fast storage, concatenating body blocks of size
`U * V * BodyW` is exactly the flattening of one grid sized
`U * V * sum(BodyW)`. This preserves the established tile structure and CSV
order while allowing logical interface members to be matched. It assumes
compatible parameter orientation and does not by itself establish geometric
or finite-element continuity.

Status update, 2026-10-24: superseded by D-007, which generalizes the same
logical-grid model to a user-selected U, V, or W stack axis.

## D-007 — Allow exactly one selected multi-body stack axis

Date: 2026-10-24

Decision: Let the user select U, V, or W through `body_stack_axis`. Bodies are
placed consecutively on that logical axis in ITD order. Counts on the other two
axes must be identical; any difference there is an error. If the public API
leaves the axis at `-1`, infer it only when exactly one count list differs.

Reason: A canonical global 3D structure and flattening should be independent of
how IRIT divides the macro geometry. Explicit body-local-to-global indexing is
required for U and V stacking because body-contiguous storage is equivalent to
global U-fast flattening only for W stacking.

Status update, 2026-10-27: superseded in the active implementation by D-008.
The imported persistent-topology cleanup executes only the first configured
U/V/W count set, so the multi-body stack-axis behavior remains in Git history
but is not present in the current source.

## D-008 — Adopt fixed-topology single-body cleanup

Date: 2026-10-27

Decision: Adopt the externally prepared persistent-topology cleanup as the
active implementation. Store FEM topology in an explicit model, store stable
element IDs on subtiles, build neighbour groups once, and reject thin-tile
removal during an optimization run.

Reason: This makes the fixed-topology invariant explicit and removes repeated
topology parsing, pointer rebinding, and geometric neighbour discovery. The
tradeoff accepted by importing this version is that only the first configured
U/V/W count set is executed; prior multi-body stack-axis support is preserved
in Git history for possible later integration.

Status update, 2026-10-28: the single-body execution restriction is superseded
by D-009. The explicit persistent FEM/lattice containers and fixed-topology
requirements remain active.

## D-009 — Keep multi-body handling at the IRIT generation boundary

Date: 2026-10-28

Decision: Represent all configured macro bodies as one rectangular logical
U/V/W lattice. Keep tiles, subtiles, neighbour groups, design updates, solver
results, and CSV data body-independent. Call IRIT once per macro trivariate and
translate each callback's local tile coordinate into the unified lattice.

Reason: Body partitioning is a property of the macro-shape parameterization,
not of the optimization variables. Restricting body-specific logic to loading,
validation, IRIT dispatch, and export-order FEM binding lets the established
single-lattice optimization pipeline remain unchanged.

Limitation: matching counts on the two non-stack axes provides a logical
correspondence only. It does not verify coincident faces, parameter direction,
watertight geometry, merged mesh nodes, or solver coupling.

