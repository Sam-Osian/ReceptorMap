from __future__ import annotations

import csv
from functools import lru_cache
import json
from pathlib import Path
from urllib.parse import quote_plus

from django import forms
from django.conf import settings
from django.contrib import admin
from django.contrib import messages
from django.core.paginator import Paginator
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.template.response import TemplateResponse
from django.urls import path, reverse
from django.utils.html import format_html
from django.utils import timezone

from .models import DatasetEditProposal


DATASET_COLUMNS = [
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
PROPOSAL_VALUE_FIELDS = (
    "drug_class",
    "activity",
    "ki_lower_nm",
    "ki_upper_nm",
)
UPDATE_COMPARISON_FIELDS = PROPOSAL_VALUE_FIELDS
NEW_CHOICE_VALUE = "__new__"


def _dataset_path() -> Path:
    return settings.BASE_DIR / "data" / "antidepressants_binding_affinities.csv"


def _to_float(value: object) -> float | None:
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _to_bool(value: object) -> bool:
    text = str(value).strip().lower()
    return text in {"1", "true", "yes"}


@lru_cache(maxsize=4)
def _load_dataset_rows_cached(dataset_mtime_ns: int) -> list[dict]:
    # `dataset_mtime_ns` is used as a cache key to auto-refresh when CSV changes.
    _ = dataset_mtime_ns
    rows: list[dict] = []
    with _dataset_path().open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row = {
                "drug": (raw.get("drug") or "").strip(),
                "drug_class": (raw.get("drug_class") or "").strip(),
                "target": (raw.get("target") or "").strip(),
                "activity": (raw.get("activity") or "").strip(),
                "activity_recoded": (raw.get("activity_recoded") or "").strip(),
                "ki_lower_nm": _to_float(raw.get("ki_lower_nm")),
                "ki_upper_nm": _to_float(raw.get("ki_upper_nm")),
                "modelled_ki_nm": _to_float(raw.get("modelled_ki_nm")),
                "ki_is_range": _to_bool(raw.get("ki_is_range")),
                "log10_ki_nm": _to_float(raw.get("log10_ki_nm")),
                "pKi": _to_float(raw.get("pKi")),
            }
            rows.append(row)
    rows.sort(key=lambda item: (item["drug"].lower(), item["target"].lower()))
    return rows


def _load_dataset_rows() -> list[dict]:
    return _load_dataset_rows_cached(_dataset_path().stat().st_mtime_ns)


def _encode_row_key(drug: str, target: str) -> str:
    return f"{quote_plus(drug)}::{quote_plus(target)}"


def _build_dataset_row_map() -> dict[str, dict]:
    return {_encode_row_key(row["drug"], row["target"]): row for row in _load_dataset_rows()}


def _choice_values(column: str) -> list[str]:
    values = {
        str(row.get(column)).strip()
        for row in _load_dataset_rows()
        if row.get(column) not in {None, ""}
    }
    return sorted(values, key=lambda item: item.lower())


def _format_row_snapshot(row: dict) -> str:
    lines = []
    for column in DATASET_COLUMNS:
        value = row.get(column)
        text = "" if value is None else str(value)
        lines.append(f"{column}: {text}")
    return "\n".join(lines)


def _display_value(value: object) -> str:
    if value in {None, ""}:
        return "-"
    return str(value)


def _proposal_change_summary(proposal: DatasetEditProposal, row_map: dict[str, dict]) -> list[str]:
    action = proposal.action
    if action == DatasetEditProposal.Action.ADD:
        return [
            f"New row: {proposal.drug} | {proposal.target}",
            f"Drug class: {_display_value(proposal.drug_class)}",
            f"Activity: {_display_value(proposal.activity)}",
            f"Ki lower bound (nM): {_display_value(proposal.ki_lower_nm)}",
            f"Ki upper bound (nM): {_display_value(proposal.ki_upper_nm)}",
        ]

    if action == DatasetEditProposal.Action.DELETE:
        key = _encode_row_key(proposal.drug or "", proposal.target or "")
        current = row_map.get(key)
        if not current:
            return ["Delete selected, but current dataset row was not found."]
        return [
            f"Delete row: {current['drug']} | {current['target']}",
            f"Drug class: {_display_value(current.get('drug_class'))}",
            f"Activity: {_display_value(current.get('activity'))}",
            f"Ki lower bound (nM): {_display_value(current.get('ki_lower_nm'))}",
            f"Ki upper bound (nM): {_display_value(current.get('ki_upper_nm'))}",
        ]

    if action == DatasetEditProposal.Action.UPDATE:
        key = _encode_row_key(proposal.drug or "", proposal.target or "")
        current = row_map.get(key)
        if not current:
            return ["Current dataset row not found. Cannot compute field-level differences."]

        changes: list[str] = []
        comparisons = (
            ("drug_class", "Drug class"),
            ("activity", "Activity"),
            ("ki_lower_nm", "Ki lower bound (nM)"),
            ("ki_upper_nm", "Ki upper bound (nM)"),
        )
        for field_name, label in comparisons:
            before = current.get(field_name)
            after = getattr(proposal, field_name)
            if before != after:
                changes.append(f"{label}: {_display_value(before)} -> {_display_value(after)}")
                if field_name == "drug_class":
                    drug_rows = sum(1 for row in row_map.values() if row.get("drug") == proposal.drug)
                    changes.append(
                        f"Scope note: Drug class updates are applied to all rows for {proposal.drug} "
                        f"({drug_rows} row(s))."
                    )
        if not changes:
            return ["No differences from current row."]
        return changes

    return ["No summary available."]


class DatasetEditProposalForm(forms.ModelForm):
    existing_row = forms.ChoiceField(
        required=False,
        label="Existing row",
        help_text="For update/delete, select an existing dataset row to avoid typing case-sensitive keys.",
    )
    current_row_snapshot = forms.CharField(
        required=False,
        label="Current row values",
        widget=forms.Textarea(attrs={"readonly": "readonly", "rows": 12}),
        help_text="Read-only snapshot of the currently selected row.",
    )
    delete_reason = forms.CharField(
        required=False,
        label="Reason for deletion",
        widget=forms.Textarea(attrs={"rows": 3}),
        help_text="Required for delete proposals.",
    )

    drug_select = forms.ChoiceField(required=False, label="Drug")
    drug_new = forms.CharField(required=False, label="Add new drug")
    target_select = forms.ChoiceField(required=False, label="Target")
    target_new = forms.CharField(required=False, label="Add new target")
    drug_class_select = forms.ChoiceField(required=False, label="Drug class")
    drug_class_new = forms.CharField(required=False, label="Add new drug class")
    activity_select = forms.ChoiceField(required=False, label="Activity")
    activity_new = forms.CharField(required=False, label="Add new activity")

    class Meta:
        model = DatasetEditProposal
        fields = "__all__"

    class Media:
        js = ("core/admin/dataset_edit_proposal_form.js",)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._row_map = _build_dataset_row_map()

        choices = [("", "---------")]
        for key, row in self._row_map.items():
            choices.append((key, f"{row['drug']} | {row['target']}"))
        self.fields["existing_row"].choices = choices
        self.fields["existing_row"].widget.attrs["data-row-map"] = json.dumps(
            self._row_map,
            separators=(",", ":"),
            ensure_ascii=True,
        )

        self._configure_choice_field("drug_select", _choice_values("drug"), "Drug")
        self._configure_choice_field("target_select", _choice_values("target"), "Target")
        self._configure_choice_field("drug_class_select", _choice_values("drug_class"), "Drug class")
        self._configure_choice_field("activity_select", _choice_values("activity"), "Activity")

        self.fields["drug_new"].widget.attrs["placeholder"] = "Enter new drug"
        self.fields["target_new"].widget.attrs["placeholder"] = "Enter new target"
        self.fields["drug_class_new"].widget.attrs["placeholder"] = "Enter new drug class"
        self.fields["activity_new"].widget.attrs["placeholder"] = "Enter new activity"
        self.fields["ki_lower_nm"].label = "Ki lower bound (nM)"
        self.fields["ki_upper_nm"].label = "Ki upper bound (nM)"
        self.fields["drug_class_select"].help_text = (
            "Update workflow: changing drug class updates all rows for the selected drug."
        )
        self.fields["ki_lower_nm"].help_text = (
            "If you do not have a Ki range, enter one Ki value in either lower or upper."
        )
        self.fields["ki_upper_nm"].help_text = (
            "If you do not have a Ki range, enter one Ki value in either lower or upper."
        )

        selected_key = self._selected_row_key()
        if selected_key and selected_key in self._row_map:
            selected_row = self._row_map[selected_key]
            self.fields["existing_row"].initial = selected_key
            self.fields["current_row_snapshot"].initial = _format_row_snapshot(selected_row)

            if not self.is_bound and self._current_action() in {
                DatasetEditProposal.Action.UPDATE,
                DatasetEditProposal.Action.DELETE,
            }:
                for field_name in PROPOSAL_VALUE_FIELDS:
                    if self.initial.get(field_name) in {None, ""}:
                        self.initial[field_name] = selected_row.get(field_name)

        if not self.is_bound:
            self._set_select_or_new_initial("drug_select", "drug_new", self.initial.get("drug"))
            self._set_select_or_new_initial("target_select", "target_new", self.initial.get("target"))
            self._set_select_or_new_initial("drug_class_select", "drug_class_new", self.initial.get("drug_class"))
            self._set_select_or_new_initial("activity_select", "activity_new", self.initial.get("activity"))

        self._apply_help_icons()

    def _add_help_icon(self, field_name: str, tooltip_text: str) -> None:
        if field_name not in self.fields:
            return
        base_label = self.fields[field_name].label or field_name.replace("_", " ").title()
        self.fields[field_name].label = format_html(
            '{} <span title="{}" style="display:inline-flex;align-items:center;justify-content:center;'
            'width:1rem;height:1rem;margin-left:0.25rem;border-radius:999px;'
            'border:1px solid #8ea0c8;color:#2f4a87;font-weight:700;font-size:0.72rem;cursor:help;">?</span>',
            base_label,
            tooltip_text,
        )

    def _apply_help_icons(self) -> None:
        self._add_help_icon("action", "Select Add, Update, or Delete to set the proposal workflow.")
        self._add_help_icon("existing_row", "Choose the dataset row to update or delete.")
        self._add_help_icon("drug_select", "Pick an existing drug name, or choose Add new.")
        self._add_help_icon("target_select", "Pick an existing target name, or choose Add new.")
        self._add_help_icon("drug_class_select", "Pick an existing drug class, or choose Add new.")
        self._add_help_icon("activity_select", "Pick an existing activity label, or choose Add new.")
        self._add_help_icon("ki_lower_nm", "Lower Ki bound in nanomolar (nM).")
        self._add_help_icon("ki_upper_nm", "Upper Ki bound in nanomolar (nM).")
        self._add_help_icon("citation", "Source supporting the proposed change.")
        self._add_help_icon("notes", "Optional context for add/update proposals.")
        self._add_help_icon("delete_reason", "Why this row should be removed.")

    def _configure_choice_field(self, field_name: str, values: list[str], noun: str) -> None:
        choices = [(NEW_CHOICE_VALUE, f"Add new {noun.lower()}..."), ("", "---------")]
        choices.extend((value, value) for value in values)
        self.fields[field_name].choices = choices

    def _set_select_or_new_initial(self, select_field: str, new_field: str, value: object) -> None:
        text = (str(value).strip() if value not in {None, ""} else "")
        if not text:
            return
        choice_values = {choice[0] for choice in self.fields[select_field].choices}
        if text in choice_values:
            self.initial[select_field] = text
        else:
            self.initial[select_field] = NEW_CHOICE_VALUE
            self.initial[new_field] = text

    def _resolve_select_or_new(self, cleaned_data: dict, select_field: str, new_field: str, label: str) -> str:
        selected = (cleaned_data.get(select_field) or "").strip()
        typed = (cleaned_data.get(new_field) or "").strip()
        if selected == NEW_CHOICE_VALUE:
            if not typed:
                self.add_error(new_field, f"Enter a value for {label.lower()}.")
                return ""
            return typed
        if selected:
            return selected
        if typed:
            return typed
        self.add_error(select_field, f"Select {label.lower()} or choose 'Add new'.")
        return ""

    def _selected_row_key(self) -> str:
        if self.is_bound:
            return (self.data.get(self.add_prefix("existing_row")) or "").strip()
        if self.instance and self.instance.pk and self.instance.drug and self.instance.target:
            return _encode_row_key(self.instance.drug, self.instance.target)
        return ""

    def _current_action(self) -> str:
        if self.is_bound:
            return (self.data.get(self.add_prefix("action")) or "").strip()
        if self.initial.get("action"):
            return str(self.initial["action"])
        if self.instance and self.instance.pk:
            return str(self.instance.action or "")
        return ""

    def clean(self):
        cleaned_data = super().clean()
        action = cleaned_data.get("action")
        selected_key = (cleaned_data.get("existing_row") or "").strip()

        selected_row = None
        if selected_key:
            selected_row = self._row_map.get(selected_key)
            if selected_row is None:
                self.add_error("existing_row", "Selected row was not found in the current dataset.")
            else:
                cleaned_data["drug"] = selected_row["drug"]
                cleaned_data["target"] = selected_row["target"]
                cleaned_data["current_row_snapshot"] = _format_row_snapshot(selected_row)

        if action == DatasetEditProposal.Action.ADD:
            cleaned_data["drug"] = self._resolve_select_or_new(cleaned_data, "drug_select", "drug_new", "Drug")
            cleaned_data["target"] = self._resolve_select_or_new(cleaned_data, "target_select", "target_new", "Target")
            cleaned_data["drug_class"] = self._resolve_select_or_new(
                cleaned_data, "drug_class_select", "drug_class_new", "Drug class"
            )
            cleaned_data["activity"] = self._resolve_select_or_new(
                cleaned_data, "activity_select", "activity_new", "Activity"
            )
            cleaned_data["activity_recoded"] = ""

        drug = (cleaned_data.get("drug") or "").strip()
        target = (cleaned_data.get("target") or "").strip()

        if action in {DatasetEditProposal.Action.UPDATE, DatasetEditProposal.Action.DELETE}:
            if not selected_key:
                raise forms.ValidationError("For update/delete, select an existing dataset row.")
            if not drug or not target:
                raise forms.ValidationError("For update/delete, select an existing dataset row.")

            current_key = _encode_row_key(drug, target)
            current_row = self._row_map.get(current_key)
            if current_row is None:
                raise forms.ValidationError(
                    "The selected drug/target does not exist in the current CSV dataset. "
                    "Use Add instead or pick an existing row."
                )
            cleaned_data["current_row_snapshot"] = _format_row_snapshot(current_row)

            if action == DatasetEditProposal.Action.UPDATE:
                cleaned_data["drug_class"] = self._resolve_select_or_new(
                    cleaned_data, "drug_class_select", "drug_class_new", "Drug class"
                )
                cleaned_data["activity"] = self._resolve_select_or_new(
                    cleaned_data, "activity_select", "activity_new", "Activity"
                )
                cleaned_data["activity_recoded"] = ""

                if cleaned_data.get("ki_lower_nm") in {None, ""}:
                    cleaned_data["ki_lower_nm"] = current_row.get("ki_lower_nm")
                if cleaned_data.get("ki_upper_nm") in {None, ""}:
                    cleaned_data["ki_upper_nm"] = current_row.get("ki_upper_nm")

                changed = False
                for field_name in UPDATE_COMPARISON_FIELDS:
                    if cleaned_data.get(field_name) != current_row.get(field_name):
                        changed = True
                        break
                if not changed:
                    raise forms.ValidationError(
                        "No changes detected for update. Edit at least one value before saving."
                    )
            else:
                reason = (cleaned_data.get("delete_reason") or "").strip()
                if not reason:
                    self.add_error("delete_reason", "Provide a reason for deletion.")
                cleaned_data["notes"] = reason

        if action == DatasetEditProposal.Action.ADD:
            if selected_key:
                raise forms.ValidationError("For add, do not select an existing row.")
        if action == DatasetEditProposal.Action.ADD and drug and target:
            if _encode_row_key(drug, target) in self._row_map:
                raise forms.ValidationError(
                    "This drug/target row already exists in the dataset. Use Update instead."
                )

        citation = (cleaned_data.get("citation") or "").strip()
        if action in {DatasetEditProposal.Action.ADD, DatasetEditProposal.Action.DELETE} and not citation:
            self.add_error("citation", "Citation is required for add and delete proposals.")

        return cleaned_data


@admin.register(DatasetEditProposal)
class DatasetEditProposalAdmin(admin.ModelAdmin):
    form = DatasetEditProposalForm
    change_list_template = "admin/core/dataset_edit_proposal/change_list.html"

    list_display = (
        "id",
        "action",
        "status",
        "drug",
        "target",
        "submitted_by",
        "reviewed_by",
        "created_at",
        "applied_at",
    )
    list_filter = ("action", "status", "activity_recoded", "ki_is_range", "created_at", "applied_at")
    search_fields = ("drug", "target", "drug_class", "activity", "citation", "notes", "processing_error")
    ordering = ("-created_at",)
    readonly_fields = ("submitted_by", "reviewed_by", "created_at", "updated_at", "reviewed_at", "applied_at")

    fieldsets = (
        (
            "Workflow",
            {
                "description": "Choose workflow first. The form sections below adapt automatically.",
                "fields": ("action", "submitted_by", "reviewed_by", "reviewed_at", "applied_at"),
            },
        ),
        (
            "Row Selection",
            {
                "classes": ("workflow-row-selection",),
                "fields": ("existing_row", "current_row_snapshot"),
            },
        ),
        (
            "Row Values (for add/update)",
            {
                "classes": ("workflow-row-values",),
                "fields": (
                    "drug_select",
                    "drug_new",
                    "target_select",
                    "target_new",
                    "drug_class_select",
                    "drug_class_new",
                    "activity_select",
                    "activity_new",
                    "ki_lower_nm",
                    "ki_upper_nm",
                )
            },
        ),
        ("Evidence", {"classes": ("workflow-evidence",), "fields": ("citation",)}),
        ("Research Context", {"classes": ("workflow-research-context",), "fields": ("notes",)}),
        ("Deletion", {"classes": ("workflow-delete-details",), "fields": ("delete_reason",)}),
        ("Timestamps", {"fields": ("created_at", "updated_at")}),
    )

    def get_urls(self):
        opts = self.model._meta
        custom_urls = [
            path(
                "review-queue/",
                self.admin_site.admin_view(self.review_queue_view),
                name=f"{opts.app_label}_{opts.model_name}_review_queue",
            ),
            path(
                "current-dataset/",
                self.admin_site.admin_view(self.current_dataset_view),
                name=f"{opts.app_label}_{opts.model_name}_current_dataset",
            )
        ]
        return custom_urls + super().get_urls()

    def changelist_view(self, request, extra_context=None):
        opts = self.model._meta
        extra_context = extra_context or {}
        extra_context["current_dataset_url"] = reverse(
            f"admin:{opts.app_label}_{opts.model_name}_current_dataset"
        )
        if self._can_approve(request):
            extra_context["review_queue_url"] = reverse(
                f"admin:{opts.app_label}_{opts.model_name}_review_queue"
            )
        return super().changelist_view(request, extra_context=extra_context)

    def _can_approve(self, request) -> bool:
        return request.user.is_superuser or request.user.has_perm("core.can_approve_dataset_edit")

    def _can_propose(self, request) -> bool:
        return request.user.is_superuser or request.user.has_perm("core.can_propose_dataset_edit")

    def _can_access(self, request) -> bool:
        return self._can_propose(request) or self._can_approve(request)

    def review_queue_view(self, request):
        if not self._can_approve(request):
            raise PermissionDenied

        opts = self.model._meta
        queue_url = reverse(f"admin:{opts.app_label}_{opts.model_name}_review_queue")

        if request.method == "POST":
            proposal_id = (request.POST.get("proposal_id") or "").strip()
            review_action = (request.POST.get("review_action") or "").strip()
            proposal = get_object_or_404(DatasetEditProposal, pk=proposal_id)

            if proposal.status != DatasetEditProposal.Status.PENDING:
                self.message_user(
                    request,
                    f"Proposal #{proposal.id} is no longer pending and cannot be reviewed.",
                    level=messages.WARNING,
                )
                return HttpResponseRedirect(queue_url)

            if review_action == "approve":
                proposal.status = DatasetEditProposal.Status.APPROVED
                proposal.rejection_reason = ""
                proposal.reviewed_by = request.user
                proposal.reviewed_at = timezone.now()
                proposal.save(update_fields=["status", "rejection_reason", "reviewed_by", "reviewed_at"])
                self.message_user(request, f"Approved proposal #{proposal.id}.")
                return HttpResponseRedirect(queue_url)

            if review_action == "reject":
                reason = (request.POST.get("reject_reason") or "").strip()
                if not reason:
                    self.message_user(
                        request,
                        "A rejection reason is required.",
                        level=messages.ERROR,
                    )
                    return HttpResponseRedirect(queue_url)
                proposal.status = DatasetEditProposal.Status.REJECTED
                proposal.rejection_reason = reason
                proposal.reviewed_by = request.user
                proposal.reviewed_at = timezone.now()
                proposal.save(update_fields=["status", "rejection_reason", "reviewed_by", "reviewed_at"])
                self.message_user(request, f"Rejected proposal #{proposal.id}.")
                return HttpResponseRedirect(queue_url)

            self.message_user(request, "Invalid review action.", level=messages.ERROR)
            return HttpResponseRedirect(queue_url)

        query = (request.GET.get("q") or "").strip()
        pending_qs = DatasetEditProposal.objects.filter(status=DatasetEditProposal.Status.PENDING).select_related(
            "submitted_by"
        )
        if query:
            pending_qs = pending_qs.filter(
                Q(drug__icontains=query)
                | Q(target__icontains=query)
                | Q(drug_class__icontains=query)
                | Q(activity__icontains=query)
                | Q(citation__icontains=query)
            )
        pending_qs = pending_qs.order_by("created_at", "id")

        paginator = Paginator(pending_qs, per_page=50)
        page_obj = paginator.get_page(request.GET.get("page"))
        row_map = _build_dataset_row_map()
        for proposal in page_obj.object_list:
            proposal.review_change_summary = _proposal_change_summary(proposal, row_map)

        context = {
            **self.admin_site.each_context(request),
            "opts": opts,
            "title": "Review queue",
            "search_query": query,
            "page_obj": page_obj,
            "is_paginated": page_obj.has_other_pages(),
            "changelist_url": reverse(f"admin:{opts.app_label}_{opts.model_name}_changelist"),
            "queue_url": queue_url,
        }
        return TemplateResponse(request, "admin/core/dataset_edit_proposal/review_queue.html", context)

    def current_dataset_view(self, request):
        if not self.has_view_permission(request):
            raise PermissionDenied

        query = (request.GET.get("q") or "").strip().lower()
        rows = _load_dataset_rows()
        if query:
            rows = [
                row
                for row in rows
                if query in row["drug"].lower()
                or query in row["target"].lower()
                or query in row["activity"].lower()
                or query in row["activity_recoded"].lower()
            ]

        paginator = Paginator(rows, per_page=50)
        page_obj = paginator.get_page(request.GET.get("page"))

        table_rows = []
        for row in page_obj.object_list:
            table_rows.append(["" if row.get(column) is None else row.get(column) for column in DATASET_COLUMNS])

        opts = self.model._meta
        context = {
            **self.admin_site.each_context(request),
            "opts": opts,
            "title": "Current Dataset (Read-only)",
            "columns": DATASET_COLUMNS,
            "rows": table_rows,
            "search_query": request.GET.get("q", ""),
            "total_rows": paginator.count,
            "page_obj": page_obj,
            "is_paginated": page_obj.has_other_pages(),
            "changelist_url": reverse(f"admin:{opts.app_label}_{opts.model_name}_changelist"),
        }
        return TemplateResponse(request, "admin/core/dataset_edit_proposal/current_dataset.html", context)

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        if self._can_approve(request):
            return qs
        if self._can_propose(request):
            return qs.filter(submitted_by=request.user)
        return qs.none()

    def has_module_permission(self, request):
        return self._can_access(request)

    def has_view_permission(self, request, obj=None):
        if not self._can_access(request):
            return False
        if self._can_approve(request) or obj is None:
            return True
        return obj.submitted_by_id == request.user.id

    def has_add_permission(self, request):
        return self._can_propose(request)

    def has_change_permission(self, request, obj=None):
        if self._can_approve(request):
            return True
        if not self._can_propose(request):
            return False
        if obj is None:
            return True
        return obj.submitted_by_id == request.user.id and obj.status == DatasetEditProposal.Status.PENDING

    def has_delete_permission(self, request, obj=None):
        if self._can_approve(request):
            return True
        if self._can_propose(request):
            if obj is None:
                return True
            return obj.submitted_by_id == request.user.id and obj.status == DatasetEditProposal.Status.PENDING
        return False

    def get_readonly_fields(self, request, obj=None):
        fields = list(super().get_readonly_fields(request, obj))
        if not self._can_approve(request):
            fields.extend(["reviewed_by", "reviewed_at"])
        return fields

    def save_model(self, request, obj, form, change):
        previous_status = None
        if change and obj.pk:
            previous = DatasetEditProposal.objects.filter(pk=obj.pk).only("status").first()
            previous_status = previous.status if previous else None

        if not obj.submitted_by_id:
            obj.submitted_by = request.user
        if not obj.pk:
            obj.status = DatasetEditProposal.Status.PENDING

        can_approve = self._can_approve(request)
        if can_approve and obj.status in {DatasetEditProposal.Status.APPROVED, DatasetEditProposal.Status.REJECTED}:
            if obj.status != previous_status:
                obj.reviewed_by = request.user
                obj.reviewed_at = timezone.now()
        elif not can_approve and change and previous_status is not None:
            obj.status = previous_status

        super().save_model(request, obj, form, change)
