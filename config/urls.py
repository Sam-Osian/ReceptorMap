from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.urls import path
from core.views import (
    about,
    axis_data,
    dataset,
    dataset_download,
    editor_audit_api,
    editor_save_api,
    editor_view,
    home,
)

urlpatterns = [
    path("", home, name="home"),
    path("about/", about, name="about"),
    path("dataset/", dataset, name="dataset"),
    path("api/axis-data/", axis_data, name="axis_data"),
    path("api/dataset/", dataset_download, name="dataset_download"),
    path("admin/", admin.site.urls),

    # ── Dataset editor ───────────────────────────────────────────────────────
    path("editor/", editor_view, name="editor"),
    path(
        "editor/login/",
        auth_views.LoginView.as_view(
            template_name="core/editor_login.html",
            redirect_authenticated_user=True,
            next_page="/editor/",
        ),
        name="editor_login",
    ),
    path(
        "editor/logout/",
        auth_views.LogoutView.as_view(next_page="/editor/login/"),
        name="editor_logout",
    ),
    path("editor/api/save/", editor_save_api, name="editor_save"),
    path("editor/api/audit/", editor_audit_api, name="editor_audit"),
]
