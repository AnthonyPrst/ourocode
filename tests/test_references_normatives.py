#! env\Scripts\python.exe
# Encoding in UTF-8
"""Tests de référence normative (golden values) — EN 1995-1-1.

Chaque cas est calculé à la main à partir des tableaux normatifs pour servir
de référence de non-régression indépendante du code (contrairement aux tests
internes qui peuvent masquer un bug de formule).

Cas de référence commun : barre bois massif C24, b=100 mm, h=300 mm,
classe de service 1, Hi = Hf (pas de retrait → b_calcul = b, h_calcul = h).

    A   = 100 x 300            = 30 000 mm²
    Iy  = 100 x 300³ / 12      = 225 000 000 mm⁴
    Wy  = 2·Iy/h               = 1 500 000 mm³
    kmod (massif, CS1, moyen terme) = 0.8     (Tableau 3.1)
    gamma_M (massif, fondamentales) = 1.3     (Tableau 2.3 AN)
    K_h = 1 pour h = 300 mm >= 150 mm
    K_m = 0.7 pour section rectangulaire massif (§6.1.6(2))
    K_cr = 0.67 pour massif h > 150 mm en CS1 (§6.1.7(2) AN)
"""

import sys
import pytest
import forallpeople as si

si.environment("structural")
sys.path.insert(1, "./")

from ourocode.eurocode.ec5 import element_droit as EC5
from ourocode.eurocode.core.objet import Objet


COMBI = {"loadtype": "Moyen terme", "typecombi": "Fondamentales"}


@pytest.fixture
def barre_c24():
    return EC5.Barre(b=100, h=300, section="Rectangulaire", classe="C24", cs=1)


class Test_Flexion_C24:
    """Flexion simple C24 — σ = M/W, f_m,d = f_m,k · kmod/γM, taux §6.1.6."""

    def test_sigma_m_d_golden(self, barre_c24):
        """σ_m,y,d = 10 kN·m / 1.5e6 mm³ = 6.6667 MPa."""
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        flexion.sigma_m_d(My=10, Mz=0)
        assert flexion.sigma_m_rd["y"].value == pytest.approx(6.666667e6, rel=1e-4)
        assert flexion.sigma_m_rd["z"].value == 0

    def test_f_m_d_golden(self, barre_c24):
        """f_m,d = 24 MPa x 0.8 / 1.3 = 14.769 MPa."""
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        flexion.f_m_d(**COMBI)
        assert flexion.f_type_rd.value == pytest.approx(14.76923e6, rel=1e-4)

    def test_taux_m_d_golden(self, barre_c24):
        """equ6.11 = σ/f = 6.6667/14.769 = 0.4514 ; equ6.12 = K_m·σ/f = 0.3160."""
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        flexion.sigma_m_d(My=10, Mz=0)
        flexion.f_m_d(**COMBI)
        flexion.taux_m_d()
        assert flexion.taux_m_rd["equ6.11"] == pytest.approx(0.45139, rel=1e-3)
        assert flexion.taux_m_rd["equ6.12"] == pytest.approx(0.31597, rel=1e-3)


class Test_Traction_C24:
    """Traction axiale C24 — σ = N/A, f_t,0,d = f_t,0,k · kmod/γM, équa. 6.1."""

    def test_sigma_t_golden(self, barre_c24):
        """σ_t,0,d = 30 kN / 30 000 mm² = 1.0 MPa... = 1 MPa attendu : 30e3/30e3."""
        traction = EC5.Traction._from_parent_class(barre_c24)
        traction.sigma_t_0_d(30)
        assert traction.sigma_t_0_rd.value == pytest.approx(1.0e6, rel=1e-4)

    def test_f_t_0_d_golden(self, barre_c24):
        """f_t,0,d = 14.5 MPa x 0.8 / 1.3 = 8.923 MPa."""
        traction = EC5.Traction._from_parent_class(barre_c24)
        traction.f_t_0_d(**COMBI)
        assert traction.f_type_rd.value == pytest.approx(8.92308e6, rel=1e-4)

    def test_taux_t_golden(self, barre_c24):
        """equ6.1 = 1.0 / 8.923 = 0.1121."""
        traction = EC5.Traction._from_parent_class(barre_c24)
        traction.sigma_t_0_d(30)
        traction.f_t_0_d(**COMBI)
        traction.taux_t_0_d()
        assert traction.taux_t_0_rd["equ6.1"] == pytest.approx(0.11207, rel=1e-3)


class Test_Cisaillement_C24:
    """Cisaillement C24 — τ = 1.5·V/(K_cr·b·h), f_v,d = f_v,k · kmod/γM."""

    def test_tau_d_golden(self, barre_c24):
        """τ_d = 1.5 x 10 kN / (0.67 x 30 000 mm²) = 0.7463 MPa."""
        cisaillement = EC5.Cisaillement._from_parent_class(barre_c24)
        cisaillement.tau_d(10)
        assert cisaillement.tau_rd.value == pytest.approx(0.74627e6, rel=1e-3)


class Test_EffetSysteme:
    """Régression : k_sys doit valoir 1 par défaut et 1.1 seulement si activé.

    Historiquement, le défaut tuple ("False", "True") était truthy et
    activait l'effet de système à tort.
    """

    def test_k_sys_default(self, barre_c24):
        assert barre_c24.effet_systeme is False
        assert barre_c24.k_sys == 1

    def test_k_sys_true_string(self, barre_c24):
        barre = EC5.Barre(b=100, h=300, classe="C24", cs=1, effet_systeme="True")
        assert barre.k_sys == 1.1

    def test_k_sys_false_string(self, barre_c24):
        barre = EC5.Barre(b=100, h=300, classe="C24", cs=1, effet_systeme="False")
        assert barre.k_sys == 1


class Test_Flag:
    """Normalisation des pseudo-enums booléens via Objet._flag."""

    def test_bool_passthrough(self):
        assert Objet._flag(True) is True
        assert Objet._flag(False) is False

    def test_strings(self):
        assert Objet._flag("True") is True
        assert Objet._flag("False") is False
        assert Objet._flag("true") is True

    def test_tuple_default_first_element(self):
        assert Objet._flag(("False", "True")) is False
        assert Objet._flag(("True", "False")) is True

    def test_none_and_numbers(self):
        assert Objet._flag(None) is False
        assert Objet._flag(1) is True
        assert Objet._flag(0) is False

    def test_unknown_string_raises(self):
        with pytest.raises(ValueError):
            Objet._flag("peut-être")


class Test_CheckChoice:
    """Validation des choix pseudo-enum via Objet._check_choice."""

    def test_valid_value(self):
        assert Objet._check_choice("C24", ["C14", "C24"], "classe") == "C24"

    def test_container_default_takes_first(self):
        assert Objet._check_choice(("C14", "C24"), ["C14", "C24"], "classe") == "C14"

    def test_dict_default_takes_first_key(self):
        assert Objet._check_choice({"A": 1, "B": 2}, ("A", "B"), "cat") == "A"

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            Objet._check_choice("C99", ["C14", "C24"], "classe")


class Test_CheckDimension:
    """Validation des dimensions via Objet._check_dimension."""

    def test_positive_ok(self):
        Objet._check_dimension(100, "b")

    def test_zero_strict_raises(self):
        with pytest.raises(ValueError):
            Objet._check_dimension(0, "b")

    def test_zero_non_strict_ok(self):
        Objet._check_dimension(0, "t", strict=False)

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            Objet._check_dimension(-5, "t", strict=False)

    def test_non_numeric_raises(self):
        with pytest.raises(ValueError):
            Objet._check_dimension("cent", "b")


class Test_FromParentClassFiltering:
    """_from_parent_class ne doit propager que les paramètres des constructeurs.

    Régression : auparavant tout le __dict__ du parent était injecté dans le
    constructeur fils (états calculés, caches, dataframes, attributs ajoutés
    dynamiquement), ce qui pouvait écraser des paramètres ou propager un état
    mutable non pertinent.
    """

    def test_constructor_params_propagated(self, barre_c24):
        """Les paramètres du constructeur Barre remontent dans la classe fille."""
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        assert flexion.b.value == barre_c24.b.value
        assert flexion.h.value == barre_c24.h.value
        assert flexion.classe == "C24"
        assert flexion.cs == barre_c24.cs

    def test_physical_units_preserved(self, barre_c24):
        """La valeur physique ne doit pas être multipliée deux fois par l'unité."""
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        # .value est en unité SI de base (m) : 100 mm = 0.1 m.
        # Si la valeur était re-multipliée par si.mm on aurait 1e-4 m.
        assert flexion.b.value == pytest.approx(0.1)
        assert flexion.h.value == pytest.approx(0.3)

    def test_dynamic_attribute_not_propagated(self, barre_c24):
        """Un attribut ajouté dynamiquement sur le parent n'est pas un paramètre
        du constructeur : il ne doit pas être propagé."""
        barre_c24.attribut_interne_calcule = {"fake": 1}
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        assert not hasattr(flexion, "attribut_interne_calcule")

    def test_private_attribute_not_propagated(self, barre_c24):
        barre_c24._cache_prive = "secret"
        flexion = EC5.Flexion._from_parent_class(barre_c24, lo_rel_y=3000, lo_rel_z=3000)
        assert not hasattr(flexion, "_cache_prive")

    def test_kwargs_override_inherited(self, barre_c24):
        """Les kwargs explicites écrasent les attributs hérités du parent."""
        flexion = EC5.Flexion._from_parent_class(barre_c24, b=50, lo_rel_y=3000, lo_rel_z=3000)
        assert flexion.b.value == pytest.approx(0.05)
