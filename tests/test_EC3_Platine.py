"""Tests unitaires pour ec3/assemblage/platine.py :
Platine_assise_compression_beton et Platine_assise_compression_bois."""
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.ec3.assemblage.platine import (
    Platine_assise_compression_beton,
    Platine_assise_compression_bois,
)
from ourocode.eurocode.ec3.assemblage.tige import Tige
from ourocode.eurocode.ec5.element_droit.barre import Barre


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tige_base():
    """Tige S235 16mm pour servir de base aux platines."""
    return Tige(
        d=16, d0=18,
        qualite="8.8",
        t=15, h=150, b=150,
        classe_acier="S235", classe_transv=1,
    )


@pytest.fixture
def barre_bois():
    """Barre bois C24 100x200 pour platine bois."""
    return Barre(
        b=100, h=200,
        section="Rectangulaire",
        classe="C24",
        cs=1,
        effet_systeme=False,
    )


@pytest.fixture
def platine_beton(tige_base):
    return Platine_assise_compression_beton._from_parent_class(
        tige_base,
        gorge=5,
        classe_beton="C25/30",
    )


@pytest.fixture
def platine_bois(tige_base, barre_bois):
    return Platine_assise_compression_bois._from_parent_class(
        tige_base,
        wood_beam=barre_bois,
        gorge=5,
        alpha_wood=0,
    )


# ---------------------------------------------------------------------------
# Test_Platine_assise_beton
# ---------------------------------------------------------------------------

class Test_Platine_assise_beton:

    def test_init_gorge(self, platine_beton):
        assert platine_beton.gorge == 5 * si.mm

    def test_init_classe_beton(self, platine_beton):
        assert platine_beton.classe_beton == "C25/30"

    def test_init_fck_positif(self, platine_beton):
        assert platine_beton.fck.value > 0

    def test_init_fck_valeur(self, platine_beton):
        """fck pour C25/30 = 25 MPa."""
        assert abs(platine_beton.fck.value - 25e6) < 1

    def test_f_jd_retourne_tuple(self, platine_beton):
        result = platine_beton.f_jd
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_f_jd_positif(self, platine_beton):
        _, fjd = platine_beton.f_jd
        assert fjd.value > 0

    def test_c_retourne_tuple(self, platine_beton):
        result = platine_beton.c()
        assert isinstance(result, tuple)

    def test_c_positif(self, platine_beton):
        _, C = platine_beton.c()
        assert C.value > 0

    def test_taux_compression_positif(self, platine_beton):
        _, taux = platine_beton.taux_compression(N_c_Ed=50, Aef=5000)
        assert float(taux) > 0

    def test_taux_compression_synthese(self, platine_beton):
        """L'appel à taux_compression alimente _synthese_taux_df."""
        platine_beton.taux_compression(N_c_Ed=50, Aef=5000)
        df = platine_beton.synthese_taux_travail()
        assert len(df) >= 1


# ---------------------------------------------------------------------------
# Test_Platine_assise_bois
# ---------------------------------------------------------------------------

class Test_Platine_assise_bois:

    def test_init_wood_beam(self, platine_bois, barre_bois):
        assert platine_bois.wood_beam is barre_bois

    def test_init_gorge(self, platine_bois):
        assert platine_bois.gorge == 5 * si.mm

    def test_init_alpha_wood(self, platine_bois):
        assert platine_bois.alpha_wood == 0

    def test_f_jd_bois_positif(self, platine_bois):
        _, fjd = platine_bois._f_jd("Court terme", "Fondamentales")
        assert float(fjd) > 0

    def test_c_bois_positif(self, platine_bois):
        _, C = platine_bois.c("Court terme", "Fondamentales")
        assert float(C) > 0

    def test_taux_compression_bois_positif(self, platine_bois):
        platine_bois._f_jd("Court terme", "Fondamentales")
        _, taux = platine_bois.taux_compression(N_c_Ed=30, Aef=3000)
        assert float(taux) > 0


if __name__ == "__main__":
    pytest.main(["-v", __file__])
