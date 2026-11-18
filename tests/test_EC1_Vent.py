"""Tests de non-régression pour le module vent (EN 1991-1-4 + AN).

Couvre notamment les correctifs :
- zonage DOM (vitesse cyclonique via le code département INSEE),
- formule 4.14-NA du coefficient orographique C0(z) en Cas 1,
- restauration de b_bat/d_bat après l'itération vent 90°,
- terrasses à rives vives (hp = 0),
- géométrie des zones à 90° (fig. 7.7(c) et 7.8(b)),
- erreur explicite au-delà de z = 200 m.
"""

import math as mt

import forallpeople as si
import pytest

from ourocode.eurocode.ec1.vent import (
    Murs_verticaux,
    Toiture_1_pant,
    Toiture_2_pants,
    Toiture_terrasse_acrotere,
    Vent,
)

si.environment("structural")

METRO = {"terrain": "II", "oro": "Aucun", "code_INSEE": 1001, "alt": 200, "h_bat": 10, "d_bat": 20, "b_bat": 30, "alpha_toit": 0}
DOM_GLP = {"terrain": "II", "oro": "Aucun", "code_INSEE": 97101, "alt": 10, "h_bat": 10, "d_bat": 20, "b_bat": 30, "alpha_toit": 0}


class TestZonageVent:
    def test_zone_metropole(self):
        # INSEE 01001 (Ain) : zone 2 -> Vb,0 = 24 m/s
        assert Vent(z=10, **METRO).zone_vent == 2
        assert Vent(z=10, **METRO).Vb_0.value == pytest.approx(24)

    def test_zone_dom_guadeloupe(self):
        v = Vent(z=10, **DOM_GLP)
        assert v.zone_vent == "Guadeloupe"
        assert v.Vb_0.value == pytest.approx(36)

    def test_zone_dom_mayotte(self):
        v = Vent(z=10, **{**DOM_GLP, "code_INSEE": 98503})
        assert v.zone_vent == "Mayotte"
        assert v.Vb_0.value == pytest.approx(30)

    def test_commune_sans_zonage(self):
        v = Vent(z=10, **{**DOM_GLP, "code_INSEE": 97501})
        with pytest.raises(ValueError, match="zonage vent"):
            _ = v.zone_vent


class TestPressionDynamique:
    def test_qp_zone2_terrain_II(self):
        """qp(z=10) = (1 + 7·Iv)·½·ρ·vm², vérifié à ~828 Pa (zones 2, cat. II)."""
        v = Vent(z=10, **METRO)
        assert v.Vb_0.value == pytest.approx(24)
        assert v.Qp_z[1].value == pytest.approx(827.6, rel=0.01)

    def test_z_sup_200_leve_erreur(self):
        v = Vent(z=250, **{**METRO, "h_bat": 250})
        with pytest.raises(ValueError):
            _ = v.Qp_z


class TestCoefficientOro:
    def test_c0_cas1_formule_414NA(self):
        """C0(z) = 1 + 0,004·ΔAC·exp(−0,014·max(0, z−10)) ; max à 1,15."""
        v = Vent(z=30, **{**METRO, "alt": 500, "h_bat": 30})
        v.oro = "Cas 1"
        v.dict_alti_oro = {k: 450 * si.m for k in ["An1", "An2", "Ae1", "Ae2", "As1", "As2", "Ao1", "Ao2"]}
        v._calc_delta_AC()
        attendu = 1 + 0.004 * 40 * mt.exp(-0.014 * 20)
        assert v.Co_z[1] == pytest.approx(attendu, rel=1e-4)

    def test_c0_z_inferieur_10(self):
        """Sous z = 10 m, l'exposant est borné à 0 : C0 = 1 + 0,004·ΔAC."""
        v = Vent(z=8, **{**METRO, "alt": 500, "h_bat": 8})
        v.oro = "Cas 1"
        v.dict_alti_oro = {k: 450 * si.m for k in ["An1", "An2", "Ae1", "Ae2", "As1", "As2", "Ao1", "Ao2"]}
        v._calc_delta_AC()
        assert v.Co_z[1] == pytest.approx(1 + 0.004 * 40, rel=1e-4)


class TestRestitutionBatiment:
    """La permutation b_bat/d_bat pour la direction 90° ne doit pas persister."""

    @pytest.mark.parametrize("cls", [Toiture_1_pant, Toiture_2_pants])
    def test_bd_restores_after_90(self, cls):
        t = cls(load_area=10, z=8, **{**METRO, "h_bat": 8, "alpha_toit": 15})
        assert t.b_bat.value == pytest.approx(30)
        assert t.d_bat.value == pytest.approx(20)

    def test_bd_restores_murs_verticaux(self):
        m = Murs_verticaux(load_area=10, z=10, **METRO)
        assert m.b_bat.value == pytest.approx(30)
        assert m.d_bat.value == pytest.approx(20)


class TestTerrasse:
    def test_rives_vives_hp_zero(self):
        t = Toiture_terrasse_acrotere(load_area=10, hp=0, h=8, z=8, **{**METRO, "h_bat": 8})
        assert t.type_terrasse == "Rives vives"
        assert t.get_Cpe().loc["CPE 10", "F"] == pytest.approx(-1.8)

    def test_acroteres(self):
        t = Toiture_terrasse_acrotere(load_area=10, hp=0.4, h=8, z=8, **{**METRO, "h_bat": 8})
        assert t.type_terrasse == "Acrotères"
        cpe = t.get_Cpe()
        assert cpe[cpe["Cpe"] == "CPE 10"]["F"].iloc[0] == pytest.approx(-1.4)

    def test_petit_rapport_hp_h_borne_basse(self):
        """hp/h < 0,025 retombe sur la ligne tabulée la plus basse au lieu de lever une erreur."""
        t = Toiture_terrasse_acrotere(load_area=10, hp=0.1, h=8, z=8, **{**METRO, "h_bat": 8})
        cpe = t.get_Cpe()
        assert cpe[cpe["Cpe"] == "CPE 10"]["F"].iloc[0] == pytest.approx(-1.6)


class TestGeometrie90:
    def test_toiture_1_pant(self):
        """Fig. 7.7(c) : G = b − e/2, H = b en longueur."""
        t = Toiture_1_pant(load_area=10, z=8, **{**METRO, "h_bat": 8, "d_bat": 20, "b_bat": 30, "alpha_toit": 15})
        geo = t.get_geo("90°")
        assert geo["G"]["Longueur"].value == pytest.approx(20 - 16 / 2)
        assert geo["G"]["Surface"].value == pytest.approx(
            16 / 10 * (20 - 16 / 2) / mt.cos(mt.radians(15))
        )
        assert geo["H"]["Longueur"].value == pytest.approx(20)
        assert geo["H"]["Surface"].value == pytest.approx(
            (16 / 2 - 16 / 10) * 20 / mt.cos(mt.radians(15))
        )

    def test_toiture_2_pants(self):
        """Fig. 7.8(b) : G par versant = b/2 − e/2 ; Surface = Longueur × Largeur / cos(α)."""
        t = Toiture_2_pants(load_area=10, z=8, **{**METRO, "h_bat": 8, "alpha_toit": 15})
        g = t.get_geo("90°")["G"]
        assert g["Longueur"].value == pytest.approx(20 / 2 - 16 / 2)
        assert g["Surface"].value == pytest.approx(g["Longueur"].value * g["Largeur"].value / mt.cos(mt.radians(15)))


class TestAccesseurs:
    def test_get_geo_sans_argument(self):
        t = Toiture_2_pants(load_area=10, z=8, **{**METRO, "h_bat": 8, "alpha_toit": 15})
        assert list(t.get_geo().keys()) == ["F", "G", "H", "I", "J"]

    def test_get_cpe_sans_argument(self):
        t = Toiture_1_pant(load_area=10, z=8, **{**METRO, "h_bat": 8, "alpha_toit": 15})
        assert "F-" in t.get_Cpe().columns
