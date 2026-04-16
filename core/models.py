from __future__ import annotations

from django.conf import settings
from django.db import models


class DatasetEditProposal(models.Model):
    class Action(models.TextChoices):
        ADD = "add", "Add"
        UPDATE = "update", "Update"
        DELETE = "delete", "Delete"

    class Status(models.TextChoices):
        PENDING = "pending", "Pending"
        APPROVED = "approved", "Approved"
        REJECTED = "rejected", "Rejected"

    action = models.CharField(max_length=16, choices=Action.choices, default=Action.UPDATE)
    status = models.CharField(max_length=16, choices=Status.choices, default=Status.PENDING)

    # Row identifier in the CSV (drug + target is currently unique in this dataset).
    drug = models.CharField(max_length=255)
    target = models.CharField(max_length=255)

    # Row values for add/update actions.
    drug_class = models.CharField(max_length=255, blank=True)
    activity = models.CharField(max_length=255, blank=True)
    activity_recoded = models.CharField(max_length=255, blank=True)
    ki_lower_nm = models.FloatField(null=True, blank=True)
    ki_upper_nm = models.FloatField(null=True, blank=True)
    modelled_ki_nm = models.FloatField(null=True, blank=True)
    ki_is_range = models.BooleanField(default=False)
    pKi = models.FloatField(null=True, blank=True)

    citation = models.TextField(blank=True)
    notes = models.TextField(blank=True)
    rejection_reason = models.TextField(blank=True)
    processing_error = models.TextField(blank=True)

    submitted_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="submitted_dataset_edits",
    )
    reviewed_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="reviewed_dataset_edits",
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    reviewed_at = models.DateTimeField(null=True, blank=True)
    applied_at = models.DateTimeField(null=True, blank=True)

    class Meta:
        ordering = ["-created_at"]
        permissions = [
            ("can_propose_dataset_edit", "Can create/update dataset edit proposals"),
            ("can_approve_dataset_edit", "Can approve/reject dataset edit proposals"),
        ]

    def __str__(self) -> str:
        return f"{self.get_action_display()} {self.drug} / {self.target} [{self.get_status_display()}]"


class CsvEditRecord(models.Model):
    """Immutable audit record for every field change made via the direct CSV editor."""

    session_id = models.CharField(max_length=64, db_index=True)
    drug = models.CharField(max_length=255)
    target = models.CharField(max_length=255)
    field_name = models.CharField(max_length=64)
    old_value = models.TextField(blank=True)
    new_value = models.TextField(blank=True)
    edited_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="csv_edits",
    )
    edited_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-edited_at"]
        permissions = [
            ("can_edit_csv", "Can directly edit the CSV dataset via the spreadsheet editor"),
        ]

    def __str__(self) -> str:
        return f"{self.drug}/{self.target}/{self.field_name}: {self.old_value!r} → {self.new_value!r}"
