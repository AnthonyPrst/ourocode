"""Tests unitaires pour ec5/blc.py : Poutre_simple_decroissance."""
import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.ec5.blc import Poutre_simple_decroissance


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def poutre():
    """Poutre à simple décroissance GL24h 100x300, alpha=5°."""
    return Poutre_simple_decroissance(
        alpha=5,
        b=100,
        h=300,
        section="Rectangulaire",
        classe="GL24h",
        cs=2,
        effet_systeme=False,
        lo_rel_y=4000,
        lo_rel_z=4000,
        coeflef_y=1.0,
        coeflef_z=1.0,
        pos="Charge sur fibre tendue",
    )


@pytest.fixture
def poutre_alpha_30():
    """Poutre à angle 30° pour tester une valeur différente."""
    return Poutre_simple_decroissance(
        alpha=30,
        b=100,
        h=300,
        section="Rectangulaire",
        classe="GL24h",
        cs=2,
        effet_systeme=False,
        lo_rel_y=4000,
        lo_rel_z=4000,
        coeflef_y=1.0,
        coeflef_z=1.0,
        pos="Charge sur fibre tendue",
    )


# ---------------------------------------------------------------------------
# Test_Poutre_simple_decroissance
# ---------------------------------------------------------------------------

class Test_Poutre_simple_decroissance:

    def test_init_alpha(self, poutre):
        """L'attribut alpha doit être stocké correctement."""
        assert poutre.alpha == 5

    def test_init_alpha_30(self, poutre_alpha_30):
        assert poutre_alpha_30.alpha == 30

    def test_herite_de_flexion(self, poutre):
        """Poutre_simple_decroissance doit hériter de Flexion."""
        from ourocode.eurocode.ec5.element_droit.flexion import Flexion
        assert isinstance(poutre, Flexion)

    def test_f_type_d_fvk_positif(self, poutre):
        """_f_type_d('fvk') doit retourner une résistance positive."""
        _, fvd = poutre._f_type_d("fvk", "Court terme", "Fondamentales")
        assert fvd.value > 0

    def test_f_type_d_fm0k_positif(self, poutre):
        """f_m_d() doit retourner une résistance positive."""
        _, fmd = poutre.f_m_d("Court terme", "Fondamentales")
        assert fmd.value > 0

    def test_f_type_d_ft90k_positif(self, poutre):
        """_f_type_d('ft90k') doit retourner une résistance positive."""
        _, ft90d = poutre._f_type_d("ft90k", "Court terme", "Fondamentales")
        assert ft90d.value > 0

    def test_f_type_d_fc90k_positif(self, poutre):
        """_f_type_d('fc90k') doit retourner une résistance positive."""
        _, fc90d = poutre._f_type_d("fc90k", "Court terme", "Fondamentales")
        assert fc90d.value > 0

    def test_attributs_section(self, poutre):
        """b et h doivent être stockés en mm."""
        assert poutre.b == 100 * si.mm
        assert poutre.h == 300 * si.mm

    def test_sigma_m_d(self, poutre):
        """sigma_m_d doit retourner un dict avec les axes y et z."""
        _, res = poutre.sigma_m_d(My=10, Mz=5)
        assert "y" in res and "z" in res
        assert res["y"].value > 0
        assert res["z"].value > 0


if __name__ == "__main__":
    pytest.main(["-v", __file__])
