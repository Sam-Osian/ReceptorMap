from django.apps import AppConfig
from django.db.models.signals import post_migrate


class CoreConfig(AppConfig):
    name = "core"

    def ready(self):
        from django.contrib.auth import get_user_model

        from .permissions import ensure_dataset_edit_groups

        post_migrate.connect(
            ensure_dataset_edit_groups,
            dispatch_uid="core.ensure_dataset_edit_groups",
        )

        # UI wording only: show "Admin access" instead of Django's "staff status".
        user_model = get_user_model()
        try:
            is_staff_field = user_model._meta.get_field("is_staff")
            is_staff_field.verbose_name = "Admin access"
        except Exception:
            pass
