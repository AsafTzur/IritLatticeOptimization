# Design Decisions

## D-014 — Carry explicit passive-skin ownership through one connected INP

Date: 2026-11-13

Decision: keep lattice and skin geometry in the same connected INP model, tag
each generated trivariate as lattice or passive skin before deformation, and
bind only explicitly lattice-owned FEM patches to optimization subtiles.

Reason: a separate independently meshed skin INP would still require node
merging or solver coupling, while splitting by element order would repeat the
unverified assumption rejected in D-010. Geometry-carried ownership preserves
one structural model and allows the skin to contribute stiffness without
receiving design updates.

Validation rule: every exported patch must have exactly one recognized owner;
all expected lattice subtiles must bind once; all remaining recognized patches
must be passive skin; and the mapping must remain identical across iterations.

## D-013 — Treat NLopt MMA as an optional build feature

Date: 2026-11-11

Decision: compile the MMA implementation only when `HAVE_NLOPT` is defined.
The Visual Studio test project controls this through `UseNLopt`, which defaults
to `false` for portability.

Reason: geometry generation and optimization algorithms 0, 2, and 3 do not
depend on NLopt and should remain buildable without its headers and library.
Algorithm 1 remains available in explicitly NLopt-enabled builds.

## D-012 — Keep IritOptRun independent of the text parser

Date: 2026-11-11

Decision: text parsing is implemented in `optimization_config_parser.c`, while
the optimization engine accepts only a completed `OptConfigStruct` through
`IritOptRun(const OptConfigStruct *)`.

Reason: the text file is the current input adapter, not part of the optimizer's
core contract. A future IRIT-language wrapper can initialize and populate the
same structure, validate it, and call `IritOptRun()` without invoking or
depending on the text parser.

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

## D-010 — Reject export-order inference for passive skin elements

Date: 2026-11-02

Decision: Remove the attempted skin-exclusion logic that parsed generated
`SET_n` ranges and treated the first fixed number of elements in every tile
range as lattice elements, with the remainder treated as passive skin.

Reason: the approach depends on an unverified ordering convention and does not
carry explicit geometry-to-FEM identity for skin versus optimizable lattice
members. A wrong classification can feed skin SED, stress, or volume into a
lattice design variable without an obvious failure. Until a replacement uses
explicit, validated ownership metadata, extra skin elements must trigger the
strict mesh/order mismatch.

Scope: skin remains part of geometry, the INP model, and structural analysis.
This decision concerns only classification for optimization and does not
define the replacement method.

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

