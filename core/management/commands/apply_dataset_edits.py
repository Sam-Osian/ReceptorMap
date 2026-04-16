from __future__ import annotations

from pathlib import Path
import math

import pandas as pd
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone

from core.models import DatasetEditProposal


CSV_COLUMNS = [
    "drug",
    "drug_class",
    "target",
    "activity",
    "activity_recoded",
    "ki_lower_nm",
    "ki_upper_nm",
    "modelled_ki_nm",
    "ki_is_range",
    "log10_ki_nm",
    "pKi",
]


def _dataset_path() -> Path:
    return settings.BASE_DIR / "data" / "antidepressants_binding_affinities.csv"


def _derive_activity_recoded(proposal: DatasetEditProposal) -> str:
    activity_text = (proposal.activity or "").strip()
    lowered = activity_text.lower()
    if not lowered:
        raise CommandError(f"Proposal {proposal.id}: `activity` is required for add/update.")

    if lowered in {"agonist", "full agonist", "partial agonist"}:
        return "Agonist"
    if lowered in {"antagonist", "inverse agonist", "inhibitor", "inhibition", "blocker"}:
        return "Antagonist"
    if "agonist" in lowered and "antagon" not in lowered and "inverse" not in lowered:
        return "Agonist"
    if any(token in lowered for token in ("antagon", "inhib", "block", "inverse agon")):
        return "Antagonist"

    raise CommandError(
        f"Proposal {proposal.id}: could not auto-derive `activity_recoded` from activity '{activity_text}'. "
        "Use a recognised activity term."
    )


def _proposal_row_values(proposal: DatasetEditProposal) -> dict:
    ki_lower_nm = proposal.ki_lower_nm
    ki_upper_nm = proposal.ki_upper_nm
    has_lower = ki_lower_nm is not None
    has_upper = ki_upper_nm is not None

    if not has_lower and not has_upper:
        raise CommandError(
            f"Proposal {proposal.id}: provide `ki_lower_nm` and/or `ki_upper_nm` for add/update."
        )
    if has_lower and ki_lower_nm <= 0:
        raise CommandError(f"Proposal {proposal.id}: `ki_lower_nm` must be > 0 when provided.")
    if has_upper and ki_upper_nm <= 0:
        raise CommandError(f"Proposal {proposal.id}: `ki_upper_nm` must be > 0 when provided.")

    if has_lower and has_upper and ki_lower_nm > ki_upper_nm:
        ki_lower_nm, ki_upper_nm = ki_upper_nm, ki_lower_nm

    ki_is_range = bool(ki_lower_nm is not None and ki_upper_nm is not None)
    modelled_ki_nm = math.sqrt(ki_lower_nm * ki_upper_nm) if ki_is_range else (ki_lower_nm or ki_upper_nm)
    log10_ki_nm = math.log10(modelled_ki_nm)
    pki = 9 - log10_ki_nm

    # Persist generated fields on the proposal record so audit entries are self-contained.
    proposal.ki_lower_nm = ki_lower_nm
    proposal.ki_upper_nm = ki_upper_nm
    proposal.ki_is_range = ki_is_range
    proposal.modelled_ki_nm = modelled_ki_nm
    proposal.pKi = pki
    proposal.activity_recoded = _derive_activity_recoded(proposal)

    return {
        "drug": proposal.drug,
        "drug_class": proposal.drug_class,
        "target": proposal.target,
        "activity": proposal.activity,
        "activity_recoded": proposal.activity_recoded,
        "ki_lower_nm": ki_lower_nm,
        "ki_upper_nm": ki_upper_nm,
        "modelled_ki_nm": modelled_ki_nm,
        "ki_is_range": ki_is_range,
        "log10_ki_nm": log10_ki_nm,
        "pKi": pki,
    }


def _validate_add_or_update(proposal: DatasetEditProposal) -> None:
    required_text = {
        "drug": proposal.drug,
        "target": proposal.target,
        "drug_class": proposal.drug_class,
        "activity": proposal.activity,
    }
    for field_name, value in required_text.items():
        if not (value or "").strip():
            raise CommandError(f"Proposal {proposal.id}: `{field_name}` is required for add/update.")

    if proposal.ki_lower_nm is None and proposal.ki_upper_nm is None:
        raise CommandError(
            f"Proposal {proposal.id}: provide `ki_lower_nm` and/or `ki_upper_nm` for add/update."
        )


class Command(BaseCommand):
    help = "Apply approved dataset edit proposals to the receptor CSV dataset."

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true", help="Validate and preview changes without writing to CSV.")
        parser.add_argument("--limit", type=int, default=0, help="Apply at most N approved proposals (0 = all).")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        limit = max(0, options["limit"])
        dataset_path = _dataset_path()
        if not dataset_path.exists():
            raise CommandError(f"Dataset file not found: {dataset_path}")

        df = pd.read_csv(dataset_path)
        missing_columns = [col for col in CSV_COLUMNS if col not in df.columns]
        if missing_columns:
            raise CommandError(f"Dataset missing required columns: {', '.join(missing_columns)}")
        extra_columns = [col for col in df.columns if col not in CSV_COLUMNS]

        proposals_qs = DatasetEditProposal.objects.filter(
            status=DatasetEditProposal.Status.APPROVED,
            applied_at__isnull=True,
        ).order_by("created_at", "id")
        if limit:
            proposals = list(proposals_qs[:limit])
        else:
            proposals = list(proposals_qs)

        if not proposals:
            self.stdout.write(self.style.WARNING("No approved unapplied proposals found."))
            return

        applied = 0
        skipped = 0
        changed_df = False

        with transaction.atomic():
            for proposal in proposals:
                mask = (df["drug"] == proposal.drug) & (df["target"] == proposal.target)
                exists = bool(mask.any())

                try:
                    row_values = None
                    if proposal.action in (
                        DatasetEditProposal.Action.ADD,
                        DatasetEditProposal.Action.UPDATE,
                    ):
                        _validate_add_or_update(proposal)
                        row_values = _proposal_row_values(proposal)

                    if proposal.action == DatasetEditProposal.Action.ADD:
                        if exists:
                            raise CommandError(
                                f"Proposal {proposal.id}: cannot add, row already exists for ({proposal.drug}, {proposal.target})."
                            )
                        row_with_optional_columns = {column: "" for column in extra_columns}
                        row_with_optional_columns.update(row_values)
                        df = pd.concat([df, pd.DataFrame([row_with_optional_columns])], ignore_index=True)
                        changed_df = True

                    elif proposal.action == DatasetEditProposal.Action.UPDATE:
                        if not exists:
                            raise CommandError(
                                f"Proposal {proposal.id}: cannot update, row not found for ({proposal.drug}, {proposal.target})."
                            )
                        current_drug_class = str(df.loc[mask, "drug_class"].iloc[0]).strip()
                        proposed_drug_class = str(proposal.drug_class or "").strip()
                        drug_class_changed = current_drug_class != proposed_drug_class

                        for col, value in row_values.items():
                            df.loc[mask, col] = value

                        # Business rule: in update workflow, a changed drug class is drug-wide.
                        if drug_class_changed:
                            drug_mask = df["drug"] == proposal.drug
                            if not bool(drug_mask.any()):
                                raise CommandError(
                                    f"Proposal {proposal.id}: cannot apply drug-wide class update; "
                                    f"no rows found for drug '{proposal.drug}'."
                                )
                            df.loc[drug_mask, "drug_class"] = proposal.drug_class
                        changed_df = True

                    elif proposal.action == DatasetEditProposal.Action.DELETE:
                        if not exists:
                            raise CommandError(
                                f"Proposal {proposal.id}: cannot delete, row not found for ({proposal.drug}, {proposal.target})."
                            )
                        df = df.loc[~mask].copy()
                        changed_df = True
                    else:
                        raise CommandError(f"Proposal {proposal.id}: unknown action `{proposal.action}`.")

                    proposal.processing_error = ""
                    if not dry_run:
                        proposal.applied_at = timezone.now()
                        proposal.save(
                            update_fields=[
                                "ki_lower_nm",
                                "ki_upper_nm",
                                "ki_is_range",
                                "modelled_ki_nm",
                                "pKi",
                                "activity_recoded",
                                "processing_error",
                                "applied_at",
                            ]
                        )
                    applied += 1
                except CommandError as exc:
                    proposal.processing_error = str(exc)
                    if not dry_run:
                        proposal.save(update_fields=["processing_error"])
                    skipped += 1
                    self.stderr.write(self.style.ERROR(str(exc)))

            if dry_run:
                self.stdout.write(self.style.WARNING("Dry run only; no file or DB changes were committed."))
                transaction.set_rollback(True)
                self.stdout.write(
                    self.style.SUCCESS(f"Validated {len(proposals)} proposal(s): {applied} pass, {skipped} fail.")
                )
                return

            if changed_df:
                output_columns = CSV_COLUMNS + [col for col in df.columns if col not in CSV_COLUMNS]
                df = df[output_columns]
                df.to_csv(dataset_path, index=False)

        self.stdout.write(
            self.style.SUCCESS(
                f"Processed {len(proposals)} proposal(s): {applied} applied, {skipped} failed. "
                f"Dataset path: {dataset_path}"
            )
        )
