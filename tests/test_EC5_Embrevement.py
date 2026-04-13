"""Tests unitaires pour ec5/assemblage/embrevement.py : Embrevement."""
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.ec5.element_droit.barre import Barre
from ourocode.eurocode.ec5.assemblage.embrevement import Embrevement


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def barre_1():
    """Pièce 1 : barre GL24h 100x200 (avec embrèvement en bout)."""
    return Barre(
        b=100, h=200,
        section="Rectangulaire",
        classe="GL24h",
        cs=2,
        effet_systeme=False,
    )


@pytest.fixture
def barre_2():
    """Pièce 2 : barre GL24h 120x220 (pièce entaillée)."""
    return Barre(
        b=120, h=220,
        section="Rectangulaire",
        classe="GL24h",
        cs=2,
        effet_systeme=False,
    )


@pytest.fixture
def embrevement_auto(barre_1, barre_2):
    """Embrèvement alpha=45° Bissectrice, profondeur calculée automatiquement."""
    return Embrevement(
        beam_1=barre_1,
        beam_2=barre_2,
        alpha=45,
        type_embrevement="Bissectrice",
        prof_embrevement=0,
        l_talon=500,
        l1=10000,
    )


@pytest.fixture
def embrevement_manuel(barre_1, barre_2):
    """Embrèvement alpha=45° Bissectrice, profondeur = 30mm imposée."""
    return Embrevement(
        beam_1=barre_1,
        beam_2=barre_2,
        alpha=45,
        type_embrevement="Bissectrice",
        prof_embrevement=30,
        l_talon=500,
        l1=10000,
    )


# ---------------------------------------------------------------------------
# Test_Embrevement_init
# ---------------------------------------------------------------------------

class Test_Embrevement_init:

    def test_alpha_stocke(self, embrevement_auto):
        assert embrevement_auto.alpha == 45

    def test_type_embrevement_stocke(self, embrevement_auto):
        assert embrevement_auto.type_embrevement == "Bissectrice"

    def test_profondeur_auto_non_nulle(self, embrevement_auto):
        """Profondeur calculée auto (alpha=45 -> tv <= h2/4)."""
        assert embrevement_auto.prof_embrevement > 0

    def test_profondeur_auto_valeur(self, embrevement_auto, barre_2):
        """Pour alpha=45°, profondeur = h2/4 = 220/4 = 55mm."""
        import forallpeople as si
        expected = barre_2.h * si.mm / 4
        assert abs(float(embrevement_auto.prof_embrevement) - float(expected)) < 1e-3

    def test_profondeur_manuel_imposee(self, embrevement_manuel):
        """Profondeur imposée = 30mm."""
        assert embrevement_manuel.prof_embrevement == 30 * si.mm

    def test_l_talon_stocke(self, embrevement_auto):
        assert embrevement_auto.l_talon == 500 * si.mm

    def test_beam_1_stocke(self, embrevement_auto, barre_1):
        assert embrevement_auto.beam_1 is barre_1

    def test_beam_2_stocke(self, embrevement_auto, barre_2):
        assert embrevement_auto.beam_2 is barre_2


# ---------------------------------------------------------------------------
# Test_Embrevement_geometrie
# ---------------------------------------------------------------------------

class Test_Embrevement_geometrie:
    """Tests des méthodes géométriques d'Embrevement."""

    def test_type_embrevement_equerre_p2_init(self, barre_1, barre_2):
        """Embrevement 'Equerre à la pièce 2' s'instancie sans erreur."""
        emb = Embrevement(
            beam_1=barre_1,
            beam_2=barre_2,
            alpha=45,
            type_embrevement="Equerre à la pièce 2",
            prof_embrevement=30,
            l_talon=500,
            l1=10000,
        )
        assert emb.type_embrevement == "Equerre à la pièce 2"
        assert emb.prof_embrevement == 30 * si.mm

    def test_type_embrevement_equerre_p1_init(self, barre_1, barre_2):
        """Embrevement 'Equerre à la pièce 1' s'instancie sans erreur."""
        emb = Embrevement(
            beam_1=barre_1,
            beam_2=barre_2,
            alpha=45,
            type_embrevement="Equerre à la pièce 1",
            prof_embrevement=30,
            l_talon=500,
            l1=10000,
        )
        assert emb.type_embrevement == "Equerre à la pièce 1"

    def test_profondeur_auto_ne_depasse_pas_h_sur_4(self, embrevement_auto, barre_2):
        """Profondeur auto <= h2/4."""
        assert embrevement_auto.prof_embrevement <= barre_2.h / 4

    def test_alpha_45_stocke(self, embrevement_auto):
        assert embrevement_auto.alpha == 45

    def test_beams_stockes(self, embrevement_auto, barre_1, barre_2):
        assert embrevement_auto.beam_1 is barre_1
        assert embrevement_auto.beam_2 is barre_2

    def test_l_talon_stocke(self, embrevement_auto):
        assert embrevement_auto.l_talon == 500 * si.mm

    def test_l1_stocke(self, embrevement_auto):
        assert embrevement_auto.l1 == 10000 * si.mm

    def test_profondeur_manuel_30(self, embrevement_manuel):
        assert embrevement_manuel.prof_embrevement == 30 * si.mm

    def test_alpha_15(self, barre_1, barre_2):
        """Alpha=15° s'instancie correctement."""
        emb = Embrevement(
            beam_1=barre_1, beam_2=barre_2, alpha=15,
            type_embrevement="Bissectrice", prof_embrevement=20,
            l_talon=400, l1=8000,
        )
        assert emb.alpha == 15

    def test_alpha_30(self, barre_1, barre_2):
        """Alpha=30° s'instancie correctement."""
        emb = Embrevement(
            beam_1=barre_1, beam_2=barre_2, alpha=30,
            type_embrevement="Bissectrice", prof_embrevement=25,
            l_talon=400, l1=8000,
        )
        assert emb.alpha == 30


if __name__ == "__main__":
    pytest.main(["-v", __file__])
