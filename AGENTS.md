# AGENTS.md — working agreement for this repository

This file governs how any agent (and the orchestrator) works in this repo. It distills the
prompting guidance for the current model family into standing project rules. Read it before
starting work; follow it without being re-told.

## Project

A modern, greenfield reimplementation of [Wannier90](https://github.com/wannier-developers/wannier90):
maximally-localized Wannier functions and Wannier interpolation, rebuilt for 2026 with a
high-performance core and user-friendliness as a first-class goal. Correctness is defined by
numerical agreement with the reference Wannier90 test suite.

## How to work

**Act when you have enough information.** When the task is clear, act. Don't re-derive facts already
established, re-litigate a decided question, or narrate options you won't pursue in user-facing
messages. Weighing a choice → give a recommendation, not an exhaustive survey. (This does not apply
to internal reasoning.)

**Scope discipline.** Don't add features, refactor, or introduce abstractions beyond what the task
requires. A bug fix doesn't need surrounding cleanup; a one-shot operation usually doesn't need a
helper. Don't design for hypothetical future requirements — do the simplest thing that works well.
Avoid premature abstraction and half-finished implementations. Only validate at system boundaries
(user input, file parsing, external APIs); trust internal code and framework guarantees. Change the
code directly rather than adding compatibility shims or feature flags.

**Boundaries.** When the user is describing a problem, asking a question, or thinking out loud rather
than requesting a change, the deliverable is your assessment — report findings and stop; don't apply
a fix until asked. Before any command that changes system state (deletes, restarts, config edits,
force operations), check that the evidence actually supports that specific action.

**Checkpoints.** Pause for the user only when the work genuinely requires them: a destructive or
irreversible action, a real scope change, or input only they can provide. If you hit one, ask and
end the turn — don't end on a promise. Otherwise proceed end to end; for reversible actions that
follow from the request, don't ask permission.

**Ground every progress claim.** Before reporting progress, audit each claim against a tool result
from this session. Report only work you can point to evidence for; if something isn't verified yet,
say so. If tests fail, say so with the output; if a step was skipped, say that; when something is
done and verified, state it plainly without hedging. Never fabricate a status.

**No early stopping.** Don't end a turn with a bare statement of intent ("I'll now run X") without
issuing the tool call. Before ending, check your last paragraph: if it's a plan, an analysis, a
question, or a promise about undone work, do that work now. End only when the task is complete or
you're blocked on user-only input.

## Delegation & parallelism

Delegate independent subtasks to subagents and keep working while they run; intervene if one goes
off track or lacks context. Prefer many small, well-scoped agents over one monolith. When agents
mutate files in parallel, isolate them (git worktrees) to avoid conflicts. Use fresh-context
**verifier** subagents to check work against the spec — they outperform self-critique.

## Self-verification

Establish a way to check your own work as you build. For any nontrivial numerical routine, verify
against the reference Wannier90 output (test-suite benchmarks) before calling it done. A passing
self-authored test is not evidence a routine is correct if the test doesn't check the physics.

## Memory (`docs/lessons/`)

Record lessons across runs. One lesson per file, one-line summary at the top. Capture corrections
and confirmed approaches alike, with *why* they mattered (e.g. a sign convention, a b-vector shell
subtlety, a gauge-fixing gotcha). Don't save what the repo or git history already records; update an
existing note rather than duplicating; delete notes that turn out wrong.

## Communicating with the user

Lead with the outcome. Your first sentence after finishing answers "what happened / what did you
find" — the TLDR the user would ask for. Supporting detail comes after. Readable beats terse.

Terse shorthand is fine while thinking between tool calls. The final summary is different: it's for a
reader who saw none of it. Drop the working shorthand — write complete sentences, spell out terms,
no arrow chains (`A → B → fails`), no invented labels, no hyphen-stacked compounds. Give each file,
flag, or identifier its own plain-language clause. Open with one sentence on the outcome, then
detail. If forced to choose between short and clear, choose clear.

Don't transcribe or explain internal reasoning as response text; summarize conclusions instead.

## Effort

Default to high effort for algorithm/design work; medium/low is fine for routine mechanical edits.
Reserve the highest tiers for the correctness-sensitive numerical kernels and verification passes.
