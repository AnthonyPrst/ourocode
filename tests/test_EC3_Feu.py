"""Tests unitaires pour ec3/feu.py : _TemperatureGaz, _CoeffFeu, Feu_acier."""
import pytest
import math
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.ec3.feu import _TemperatureGaz, _CoeffFeu, Feu_acier
from ourocode.eurocode.ec3.element_droit.plat import Plat


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def plat():
    return Plat(t=10, h=200, b=100, classe_acier="S235", classe_transv=1)


@pytest.fixture
def feu_acier(plat):
    return Feu_acier._from_parent_class(
        plat,
        am_v=100,
        time=30,
        theta_a0=20.0,
        epsilon_m=0.7,
        k_sh=1.0,
    )


# ---------------------------------------------------------------------------
# Test_TemperatureGaz
# ---------------------------------------------------------------------------

class Test_TemperatureGaz:

    def test_iso_834_t0(self):
        """A t=0s (0 min) la formule donne 20 + 345*log10(1) = 20°C."""
        assert _TemperatureGaz.iso_834(0) == 20.0

    def test_iso_834_t60s(self):
        """A 60s (1 min) : 20 + 345*log10(9) ≈ 349°C."""
        expected = 20.0 + 345.0 * math.log10(8.0 * 1.0 + 1.0)
        assert abs(_TemperatureGaz.iso_834(60) - expected) < 1e-6

    def test_iso_834_t3600s(self):
        """A 3600s (60 min) : valeur positive et > 800°C."""
        val = _TemperatureGaz.iso_834(3600)
        assert val > 800

    def test_courbe_iso_longueur_30min(self):
        """courbe_iso(30) doit retourner une liste avec 30*60/5 + 1 = 361 points."""
        courbe = _TemperatureGaz.courbe_iso(30)
        assert len(courbe) == 30 * 60 // 5 + 1

    def test_courbe_iso_premier_point(self):
        """Le premier point de la courbe à 30 min doit être 20°C (t=0)."""
        courbe = _TemperatureGaz.courbe_iso(30)
        assert courbe[0] == 20.0

    def test_courbe_iso_croissant(self):
        """La courbe de température doit être strictement croissante."""
        courbe = _TemperatureGaz.courbe_iso(10)
        assert all(courbe[i] < courbe[i + 1] for i in range(len(courbe) - 1))


# ---------------------------------------------------------------------------
# Test_CoeffFeu
# ---------------------------------------------------------------------------

class Test_CoeffFeu:

    def test_ky_theta_20_est_1(self):
        """A 20°C ky,θ = 1.0."""
        coef = _CoeffFeu(theta=20)
        assert coef.ky_theta == 1.0

    def test_ky_theta_1200_est_0(self):
        """A 1200°C ky,θ = 0.0."""
        coef = _CoeffFeu(theta=1200)
        assert coef.ky_theta == 0.0

    def test_ky_theta_interpolation_550(self):
        """A 550°C : interpolation linéaire entre (500, 0.78) et (600, 0.47)."""
        coef = _CoeffFeu(theta=550)
        expected = 0.78 + (550 - 500) / (600 - 500) * (0.47 - 0.78)
        assert abs(coef.ky_theta - expected) < 1e-9

    def test_kb_theta_bornes(self):
        """k_b,θ = 1 à 20°C et 0 à 1000°C."""
        assert _CoeffFeu(theta=20).kb_theta == 1.0
        assert _CoeffFeu(theta=1000).kb_theta == 0.0

    def test_kb_theta_interpolation_250(self):
        """k_b,θ à 250°C interpolé entre (200, 0.935) et (300, 0.903)."""
        coef = _CoeffFeu(theta=250)
        expected = 0.935 + (250 - 200) / (300 - 200) * (0.903 - 0.935)
        assert abs(coef.kb_theta - expected) < 1e-9

    def test_kp_theta_20(self):
        """k_p,θ = 1 à 20°C."""
        assert _CoeffFeu(theta=20).kp_theta == 1.0

    def test_kE_theta_20(self):
        """k_E,θ = 1 à 20°C."""
        assert _CoeffFeu(theta=20).kE_theta == 1.0

    def test_kw_theta_20(self):
        """k_w,θ = 1 à 20°C."""
        assert _CoeffFeu(theta=20).kw_theta == 1.0

    def test_synthese_contient_toutes_les_cles(self):
        """synthese() doit contenir exactement les 6 clés attendues."""
        coef = _CoeffFeu(theta=400)
        s = coef.synthese()
        assert set(s.keys()) == {
            "theta_C", "k_y_theta", "k_p_theta", "k_E_theta", "k_b_theta", "k_w_theta"
        }

    def test_synthese_valeurs_cohérentes(self):
        """Les valeurs de synthese() doivent correspondre aux propriétés."""
        coef = _CoeffFeu(theta=600)
        s = coef.synthese()
        assert s["theta_C"] == 600
        assert s["k_y_theta"] == coef.ky_theta
        assert s["k_b_theta"] == coef.kb_theta

    def test_interp_en_dessous_borne_inférieure(self):
        """Température < 20°C retourne la valeur en 20°C."""
        coef = _CoeffFeu(theta=0)
        assert coef.ky_theta == 1.0

    def test_interp_au_dessus_borne_superieure(self):
        """Température > 1200°C retourne la valeur en 1200°C."""
        coef = _CoeffFeu(theta=1500)
        assert coef.ky_theta == 0.0


# ---------------------------------------------------------------------------
# Test_Feu_acier
# ---------------------------------------------------------------------------

class Test_Feu_acier:

    def test_init(self, feu_acier):
        """Attributs initialisés correctement."""
        assert feu_acier.time == 30
        assert feu_acier.am_v == 100
        assert feu_acier.epsilon_m == 0.7
        assert feu_acier.k_sh == 1.0

    def test_fire_data_non_vide(self, feu_acier):
        """fire_data doit être un DataFrame non vide."""
        import pandas as pd
        assert isinstance(feu_acier.fire_data, pd.DataFrame)
        assert not feu_acier.fire_data.empty

    def test_fire_data_colonnes(self, feu_acier):
        """fire_data doit contenir toutes les colonnes attendues."""
        cols = feu_acier.fire_data.columns
        for col in ["t (s)", "θg (°C)", "θa (°C)", "ky,θ", "kp,θ", "kE,θ", "kb,θ", "kw,θ"]:
            assert col in cols

    def test_theta_a_croissant(self, feu_acier):
        """La température de l'acier doit augmenter au cours du temps."""
        theta_debut = feu_acier.fire_data["θa (°C)"].iloc[0]
        theta_fin = feu_acier.fire_data["θa (°C)"].iloc[-1]
        assert theta_fin > theta_debut

    def test_get_fire_data_retourne_dataframe(self, feu_acier):
        """get_fire_data() retourne le DataFrame complet."""
        import pandas as pd
        df = feu_acier.get_fire_data()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == feu_acier.time * 60 // 5 + 1

    def test_get_specific_time_data_valide(self, feu_acier):
        """Accès à t=15 min retourne une série avec les bonnes clés."""
        row = feu_acier.get_specific_time_data(15)
        assert "θa (°C)" in row.index
        assert row["θa (°C)"] > 20

    def test_get_specific_time_invalide_raises(self, feu_acier):
        """Temps hors-bornes lève ValueError."""
        with pytest.raises(ValueError):
            feu_acier.get_specific_time_data(60)

    def test_get_specific_time_negatif_raises(self, feu_acier):
        """Temps négatif lève ValueError."""
        with pytest.raises(ValueError):
            feu_acier.get_specific_time_data(-1)

    def test_c_a_branche_moins_600(self):
        """_c_a < 600°C : formule polynomiale."""
        val = Feu_acier._c_a(300)
        assert val > 425

    def test_c_a_branche_600_735(self):
        """_c_a entre 600 et 735°C : pic de capacité thermique."""
        val = Feu_acier._c_a(700)
        assert val > 600

    def test_c_a_branche_735_900(self):
        """_c_a entre 735 et 900°C."""
        val = Feu_acier._c_a(800)
        assert val > 545

    def test_c_a_branche_900_1200(self):
        """_c_a entre 900 et 1200°C : constante à 650 J/kgK."""
        assert Feu_acier._c_a(1000) == 650

    def test_c_a_en_dessous_20_raises(self):
        """Température < 20°C lève ValueError."""
        with pytest.raises(ValueError):
            Feu_acier._c_a(10)

    def test_c_a_au_dessus_1200_raises(self):
        """Température > 1200°C lève ValueError."""
        with pytest.raises(ValueError):
            Feu_acier._c_a(1300)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
