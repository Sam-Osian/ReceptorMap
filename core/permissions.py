from django.contrib.auth.models import Group, Permission


def ensure_dataset_edit_groups(**kwargs):
    proposers_group, _ = Group.objects.get_or_create(name="Dataset Proposers")
    approvers_group, _ = Group.objects.get_or_create(name="Dataset Approvers")

    propose_perm = Permission.objects.filter(
        codename="can_propose_dataset_edit",
        content_type__app_label="core",
        content_type__model="dataseteditproposal",
    ).first()
    approve_perm = Permission.objects.filter(
        codename="can_approve_dataset_edit",
        content_type__app_label="core",
        content_type__model="dataseteditproposal",
    ).first()

    if propose_perm is not None:
        proposers_group.permissions.add(propose_perm)
    if approve_perm is not None:
        approvers_group.permissions.add(approve_perm)
