# ReceptorMap Project Plan

## Aim

Build ReceptorMap into a closed, curated pharmacology knowledge tool for mood-disorder prescribing education.

The system should combine:

- reviewed drug-target affinity data,
- trusted source documents,
- evidence-backed clinical/problem summaries,
- numeric receptor-profile comparison.

It should not become an open web-search chatbot. AI should help extract and summarise from approved sources, with human review before anything becomes canonical.

## Scope

Initial scope:

- Mood disorders, starting with depression.
- Commonly used antidepressants/augmenting agents.
- BAP guideline references as the first evidence source.
- SMPCs and Wikipedia-derived reference trails as secondary source-gathering routes.

Out of scope for the first build:

- Full psychiatry coverage.
- Autonomous prescribing recommendations.
- Unreviewed AI-generated dataset changes.
- General internet retrieval as an answer source.

## Workstreams

## 1. Evidence Collection

Create a clean source corpus in `Pharmacology Project 2026`.

Folder pattern:

```text
Drug name/
  Guideline papers/
  SMPC/
  Wiki searching papers/
```

First pass:

- Identify the initial BAP depression drug list.
- Collect BAP-cited papers into `Guideline papers`.
- Track progress in a spreadsheet.
- Upload the previous med-student/Wikipedia extraction table.

Output:

- A small, organised source set for the first few drugs.

## 2. Data Model

Move beyond the current flat CSV shape by adding database models for:

- source documents,
- reviewed pharmacology facts,
- clinical problems,
- evidence claims,
- links between problems, mechanisms, and drugs.

The CSV can remain an export/display format for now.

Output:

- A schema that preserves provenance, review status, source type, quotes/excerpts, and reviewer decisions.

## 3. Review Workflow

Extend the current admin/review pattern so extracted evidence is proposed, reviewed, then approved or rejected.

Required proposal fields:

- drug,
- target,
- activity,
- Ki/pKi values,
- source document,
- source URL/identifier,
- evidence note,
- quote or page/section,
- extraction confidence,
- reviewer status.

Output:

- A safe path from extracted evidence to approved dataset entries.

## 4. Ingestion And AI Extraction

Build this only after the first corpus and schema are in place.

Pipeline:

```text
PDF/source folder
  -> text extraction
  -> AI candidate extraction
  -> review queue
  -> approved structured facts
```

Start with 2-3 drugs only.

Output:

- AI-generated proposals, not direct writes to the dataset.

## 5. Query Features

Extend the existing receptor-axis app into three query modes.

### View 1: Evidence Summary

Answer educational questions such as:

```text
What evidence supports options for sleep onset difficulty?
```

Use only approved evidence and show citations.

### View 2: Receptor Axis

Build on the current plot:

- show evidence-backed vs inferred activity,
- add source detail per point,
- support queries like weaker/stronger/opposite activity at a receptor.

### View 3: Similarity View

Represent each drug as a signed receptor vector:

```text
agonist pKi = positive
antagonist/inhibitor pKi = negative
unknown = missing
```

Use this to find drugs that are pharmacologically similar, different, or opposite across selected receptors.

## Phases

## Phase 1: Source Collection

Duration: 1-2 weeks.

Tasks:

- Confirm initial drug list.
- Divide collection work.
- Populate BAP guideline paper folders.
- Add tracking spreadsheet.
- Upload previous extraction table.

Exit criteria:

- At least a few drugs have clean, source-backed folders ready for testing.

## Phase 2: Schema And Review

Duration: 2-4 weeks.

Tasks:

- Add source/evidence models.
- Add admin review screens.
- Add approved/rejected/proposed status.
- Add CSV export or sync from approved facts.

Exit criteria:

- A reviewer can approve a source-backed pharmacology fact and have it feed the app.

## Phase 3: Extraction Pilot

Duration: 3-5 weeks.

Tasks:

- Import documents from the folder structure.
- Extract text.
- Run AI extraction on 2-3 drugs.
- Send outputs to review.
- Measure accuracy and reviewer workload.

Exit criteria:

- AI extraction produces useful reviewable proposals.

## Phase 4: Query Expansion

Duration: 4-8 weeks.

Tasks:

- Add source details to receptor-axis points.
- Add receptor action filters.
- Build signed-pKi similarity API.
- Prototype evidence summary view over approved claims.

Exit criteria:

- Users can query by receptor action, drug similarity, and selected clinical problems using reviewed evidence.

## Immediate Next Steps

1. Upload the previous med-student/Wikipedia extraction table.
2. Confirm the initial BAP drug list.
3. Start source collection for the first drugs.
4. Draft the new Django evidence models.
5. Decide whether to extend `DatasetEditProposal` or create a separate evidence proposal model.
6. Pilot the full flow on 2-3 drugs before scaling.

## Main Risks

- Source collection takes longer than expected.
- Evidence quality varies by drug/target.
- AI extraction invents or misreads values.
- Copyright/access constraints affect document storage.
- Similarity results are misleading until activity coverage improves.

## Guiding Rule

Build the reviewed evidence pipeline before building the AI-facing query experience.
