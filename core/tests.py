from django.test import TestCase
from django.urls import reverse

from core.views import _load_affinity_data


class HomeViewTests(TestCase):
    def test_home_renders(self) -> None:
        response = self.client.get(reverse("home"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Receptor Axis (Antagonist vs Agonist)")

    def test_home_accepts_drug_filter(self) -> None:
        df = _load_affinity_data()
        sample_drug = str(df["drug"].iloc[0])
        response = self.client.get(reverse("home"), {"drug": sample_drug})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, f"Top Modelled pKi Targets for {sample_drug}")
