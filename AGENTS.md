# Project Instructions

This directory contains Asaf Tzur's IRIT lattice-optimization research code.

Before modifying code, read:

1. `docs/project-history/PROJECT_CONTEXT.md`
2. `docs/project-history/CHANGELOG.md`
3. `docs/project-history/DECISIONS.md`
4. `docs/project-history/THESIS_NOTES.md`

Also read any topic-specific record linked from those files, such as
`SKIN_BOUNDARY_CHANGELOG.md`.

## Cross-chat continuity

Treat the Markdown records in `docs/project-history/` as the persistent project
memory shared across separate Codex chats. A new chat must not assume that it
remembers earlier conversations.

At the start of every project-related chat—not only before modifying code—read
the four core project-history files listed above. Then read
`docs/project-history/SESSION_HANDOFF.md` when it exists, followed by any
topic-specific records relevant to the request. Use these records to recover
the current state before answering questions, diagnosing behavior, planning
work, or making changes.

Maintain `docs/project-history/SESSION_HANDOFF.md` as a concise description of
the latest working state. Create it when substantial project work first needs
to be handed to a future chat. Update it during or at the end of substantial
work when omitting the information could cause a future chat to repeat work,
lose an unresolved issue, or make a materially different decision. Record:

- the user's current objective and agreed scope;
- important findings and decisions, with links to the durable history record
  where applicable;
- files changed and the reason for each material change;
- commands or tests actually run and their observed results;
- unresolved problems, limitations, and concrete next steps.

Do not use the handoff as a full chat transcript or duplicate the permanent
history files. Keep durable architecture decisions in `DECISIONS.md`, material
source-code changes in `CHANGELOG.md`, thesis-relevant information in
`THESIS_NOTES.md`, and stable project structure and terminology in
`PROJECT_CONTEXT.md`. Never claim to remember an earlier chat unless its
content is present in the active conversation or recorded in the workspace.
If a required record cannot be read or updated, state that limitation clearly
and do not imply that cross-chat continuity is complete.

## Change-recording rules

- Preserve the user's existing code and distinguish a new feature from a
  refactor, relocation, or validation of existing behavior.
- For every material source-code change, add a dated entry to
  `docs/project-history/CHANGELOG.md` describing the purpose, files changed,
  behavior affected, and verification performed.
- Record important architectural choices and rejected alternatives in
  `docs/project-history/DECISIONS.md`.
- Add thesis-relevant methods, observations, limitations, and experimental
  evidence to `docs/project-history/THESIS_NOTES.md`.
- Update `docs/project-history/PROJECT_CONTEXT.md` when architecture,
  terminology, dependencies, or workflow stages change.
- Never claim a test passed unless its command and result were observed.
- Label reconstructed history as reconstructed when no version-control record
  proves the exact change.
- Do not store passwords, API keys, license keys, or other secrets in these
  documentation files.

## Current critical limitation

Multi-body generation processes each macro body separately. The present
neighbour-thickness matching operates independently inside each body and does
not guarantee geometric or finite-element continuity across body interfaces.

