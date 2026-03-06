from django.test import TestCase
from django.urls import reverse

from core.views import _load_affinity_data


class HomeViewTests(TestCase):
    def test_home_renders(self) -> None:
        response = self.client.get(reverse("home"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Axis Plot")

    def test_home_accepts_drug_and_receptor_filter(self) -> None:
        df = _load_affinity_data()
        sample_drug = str(df["drug"].iloc[0])
        sample_receptor = str(df.loc[df["drug"] == sample_drug, "target"].iloc[0])
        response = self.client.get(reverse("home"), {"drug": sample_drug, "receptor": sample_receptor})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, f"Drug: {sample_drug}")
        self.assertContains(response, f"Receptor: {sample_receptor}")

    def test_axis_data_endpoint_returns_json(self) -> None:
        df = _load_affinity_data()
        sample_drug = str(df["drug"].iloc[0])
        sample_receptor = str(df.loc[df["drug"] == sample_drug, "target"].iloc[0])
        response = self.client.get(reverse("axis_data"), {"drug": sample_drug, "receptor": sample_receptor})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "application/json")
        payload = response.json()
        self.assertEqual(payload["selected_drug"], sample_drug)
        self.assertEqual(payload["selected_receptor"], sample_receptor)
        self.assertIn("comparison_rows", payload)
