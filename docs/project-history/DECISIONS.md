# Design Decisions

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

