from __future__ import annotations

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


def create_dataset_editors_group(apps, schema_editor):
    Group = apps.get_model("auth", "Group")
    Permission = apps.get_model("auth", "Permission")
    ContentType = apps.get_model("contenttypes", "ContentType")

    editors_group, _ = Group.objects.get_or_create(name="Dataset Editors")

    content_type = ContentType.objects.filter(
        app_label="core", model="csveditrecord"
    ).first()
    if content_type is None:
        return

    perm = Permission.objects.filter(
        content_type=content_type,
        codename="can_edit_csv",
    ).first()
    if perm:
        editors_group.permissions.add(perm)


def reverse_dataset_editors_group(apps, schema_editor):
    Group = apps.get_model("auth", "Group")
    Group.objects.filter(name="Dataset Editors").delete()


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ("core", "0003_create_dataset_edit_groups"),
    ]

    operations = [
        migrations.CreateModel(
            name="CsvEditRecord",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("session_id", models.CharField(db_index=True, max_length=64)),
                ("drug", models.CharField(max_length=255)),
                ("target", models.CharField(max_length=255)),
                ("field_name", models.CharField(max_length=64)),
                ("old_value", models.TextField(blank=True)),
                ("new_value", models.TextField(blank=True)),
                (
                    "edited_by",
                    models.ForeignKey(
                        blank=True,
                        null=True,
                        on_delete=django.db.models.deletion.SET_NULL,
                        related_name="csv_edits",
                        to=settings.AUTH_USER_MODEL,
                    ),
                ),
                ("edited_at", models.DateTimeField(auto_now_add=True)),
            ],
            options={
                "ordering": ["-edited_at"],
                "permissions": [
                    (
                        "can_edit_csv",
                        "Can directly edit the CSV dataset via the spreadsheet editor",
                    )
                ],
            },
        ),
        migrations.RunPython(
            create_dataset_editors_group,
            reverse_dataset_editors_group,
        ),
    ]
