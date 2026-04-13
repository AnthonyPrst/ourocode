"""Tests unitaires pour Projet, Objet et les mixins (Serialization, MathUtils, Synthese, DataLoader)."""
import pytest
import math
import sys
import os
import json
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.core.projet import Projet
from ourocode.eurocode.core.objet import Objet
from ourocode.eurocode.mixins.serialization import SerializationMixin
from ourocode.eurocode.mixins.math_utils import MathUtilsMixin
from ourocode.eurocode.mixins.synthese import SyntheseMixin
from ourocode.eurocode.ec3.element_droit.plat import Plat


# ---------------------------------------------------------------------------
# Fixtures communes
# ---------------------------------------------------------------------------

@pytest.fixture
def projet_vide():
    return Projet()


@pytest.fixture
def projet_rempli():
    return Projet(
        ingenieur="Dupont",
        num_project="2024-001",
        name="Bâtiment A",
        adresse="1 rue de la Paix",
        code_INSEE=75056,
        pays="France",
        alt=300,
    )


@pytest.fixture
def plat_simple():
    return Plat(t=10, h=100, b=200, classe_acier="S235", classe_transv=1)


# ---------------------------------------------------------------------------
# Test_Projet
# ---------------------------------------------------------------------------

class Test_Projet:
    def test_init_defaults(self, projet_vide):
        p = projet_vide
        assert p.ingenieur is None
        assert p.num_project is None
        assert p.name is None
        assert p.adresse is None
        assert p.code_INSEE is None
        assert p.pays == "France"
        assert p.alt == 0 * si.m

    def test_init_avec_valeurs(self, projet_rempli):
        p = projet_rempli
        assert p.ingenieur == "Dupont"
        assert p.num_project == "2024-001"
        assert p.name == "Bâtiment A"
        assert p.code_INSEE == 75056
        assert p.alt == 300 * si.m

    def test_kwargs_dynamiques(self):
        p = Projet(commentaire="test", version=2)
        assert p.commentaire == "test"
        assert p.version == 2

    def test_from_dict(self):
        p = Projet._from_dict({"ingenieur": "Martin", "alt": 100})
        assert p.ingenieur == "Martin"
        assert p.alt == 100 * si.m

    def test_from_parent_class_simple(self, projet_rempli):
        p2 = Projet._from_parent_class(projet_rempli, name="Bâtiment B")
        assert p2.ingenieur == "Dupont"
        assert p2.name == "Bâtiment B"
        assert p2.code_INSEE == 75056

    def test_from_parent_class_list(self):
        p1 = Projet(ingenieur="Alice", alt=100)
        p2 = Projet(ingenieur="Bob", name="Projet X", alt=200)
        merged = Projet._from_parent_class([p1, p2])
        # Le dernier objet écrase le premier -> ingenieur de p2
        assert merged.ingenieur == "Bob"
        assert merged.name == "Projet X"
        assert merged.alt == 200 * si.m

    def test_classe_acier_invalide_raises(self):
        with pytest.raises(ValueError):
            Plat(t=10, h=100, b=200, classe_acier="S999", classe_transv=1)

    def test_classe_transv_invalide_raises(self):
        with pytest.raises(ValueError):
            Plat(t=10, h=100, b=200, classe_acier="S235", classe_transv=4)


# ---------------------------------------------------------------------------
# Test_Objet_get_value
# ---------------------------------------------------------------------------

class Test_Objet_get_value:
    def test_get_value_dict_key(self, projet_vide):
        data = {"a": 1, "b": 2}
        assert projet_vide.get_value(data, key="a", get_keys=False) == 1

    def test_get_value_list_index(self, projet_vide):
        data = [10, 20, 30]
        assert projet_vide.get_value(data, index=1) == 20

    def test_get_value_get_keys(self, projet_vide):
        data = {"x": 1, "y": 2}
        keys = projet_vide.get_value(data, get_keys=True)
        assert set(keys) == {"x", "y"}

    def test_get_value_json_string(self, projet_vide):
        data = "{'alpha': 42, 'beta': 99}"
        assert projet_vide.get_value(data, key="alpha") == 42


# ---------------------------------------------------------------------------
# Test_SerializationMixin
# ---------------------------------------------------------------------------

class Test_SerializationMixin:
    def test_physical_to_dict_simple(self, projet_vide):
        val = 500 * si.mm
        result = projet_vide._physical_to_dict(val)
        assert isinstance(result, dict)
        assert "_physical_value" in result
        assert "_physical_unit" in result

    def test_physical_to_dict_roundtrip_mm(self, projet_vide):
        val = 250 * si.mm
        d = projet_vide._physical_to_dict(val)
        restored = projet_vide._dict_to_physical(d)
        assert abs(float(restored) - float(val)) < 1e-9

    def test_physical_to_dict_roundtrip_kN(self, projet_vide):
        val = 12.5 * si.kN
        d = projet_vide._physical_to_dict(val)
        restored = projet_vide._dict_to_physical(d)
        assert abs(float(restored) - float(val)) < 1e-6

    def test_physical_to_dict_roundtrip_MPa(self, projet_vide):
        val = 235 * si.MPa
        d = projet_vide._physical_to_dict(val)
        restored = projet_vide._dict_to_physical(d)
        assert abs(float(restored) - float(val)) < 1e-3

    def test_physical_to_dict_nested(self, projet_vide):
        data = {"force": 10 * si.kN, "liste": [5 * si.mm, 3 * si.mm]}
        result = projet_vide._physical_to_dict(data)
        assert isinstance(result["force"], dict)
        assert isinstance(result["liste"], list)
        assert isinstance(result["liste"][0], dict)

    def test_resolve_unit_simple(self):
        m = SerializationMixin._resolve_unit_expr("m")
        mm = SerializationMixin._resolve_unit_expr("mm")
        kN = SerializationMixin._resolve_unit_expr("kN")
        assert float(1 * m) == float(1 * si.m)
        assert float(1 * mm) == float(1 * si.mm)
        assert float(1 * kN) == float(1 * si.kN)

    def test_resolve_unit_power(self):
        mm2 = SerializationMixin._resolve_unit_expr("mm*mm")
        assert float(1 * mm2) == float(1 * si.mm**2)

    def test_resolve_unit_product(self):
        kNm = SerializationMixin._resolve_unit_expr("kN*m")
        assert float(1 * kNm) == float(1 * si.kN * si.m)

    def test_resolve_unit_unknown_raises(self):
        with pytest.raises(ValueError):
            SerializationMixin._resolve_unit_expr("furlongs")

    def test_save_load_json_roundtrip(self, projet_vide, tmp_path):
        data = {"masse": 50 * si.kN, "longueur": 3000 * si.mm}
        path = str(tmp_path / "test_data.json")
        projet_vide.save_data(data, type_data="JSON", path=path)
        loaded = projet_vide.load_data(type_data="JSON", path=path)
        assert abs(float(loaded["masse"]) - float(data["masse"])) < 1e-6
        assert abs(float(loaded["longueur"]) - float(data["longueur"])) < 1e-9


# ---------------------------------------------------------------------------
# Test_MathUtilsMixin
# ---------------------------------------------------------------------------

class Test_MathUtilsMixin:
    def test_abs_value(self, projet_vide):
        assert projet_vide.abs_value(-5.3) == 5.3
        assert projet_vide.abs_value(3.7) == 3.7

    def test_max_et_min(self, projet_vide):
        assert projet_vide.max(3, 7) == 7
        assert projet_vide.min(3, 7) == 3

    def test_extract_numbers_list(self, projet_vide):
        nums = projet_vide._extract_numbers([1.0, -2.0, 3.0], absolute=False)
        assert set(nums) == {1.0, -2.0, 3.0}

    def test_extract_numbers_absolute_false(self, projet_vide):
        nums = projet_vide._extract_numbers([-5.0, 3.0], absolute=False)
        assert -5.0 in nums

    def test_extract_numbers_nested_dict(self, projet_vide):
        data = {"a": 1.0, "b": {"c": 2.0}}
        nums = projet_vide._extract_numbers(data)
        assert 1.0 in nums and 2.0 in nums

    def test_extract_numbers_physical(self, projet_vide):
        nums = projet_vide._extract_numbers([10 * si.kN, 5 * si.kN])
        assert len(nums) == 2
        assert all(n > 0 for n in nums)

    def test_max_list(self, projet_vide):
        assert projet_vide.max_list([1.0, -5.0, 3.0]) == 5.0

    def test_min_list(self, projet_vide):
        assert projet_vide.min_list([1.0, -5.0, 3.0]) == 1.0

    def test_operation_between_values(self, projet_vide):
        assert projet_vide.operation_between_values(4, 2, "+") == 6
        assert projet_vide.operation_between_values(4, 2, "-") == 2
        assert projet_vide.operation_between_values(4, 2, "x") == 8
        assert projet_vide.operation_between_values(4, 2, "/") == 2

    def test_operation_invalide_raises(self, projet_vide):
        with pytest.raises(ValueError):
            projet_vide.operation_between_values(4, 2, "**")

    def test_get_trigonometric_sin_90(self, projet_vide):
        assert abs(projet_vide.get_trigonometric_value(90, "SIN") - 1.0) < 1e-10

    def test_get_trigonometric_cos_0(self, projet_vide):
        assert abs(projet_vide.get_trigonometric_value(0, "COS") - 1.0) < 1e-10

    def test_get_trigonometric_tan_45(self, projet_vide):
        assert abs(projet_vide.get_trigonometric_value(45, "TAN") - 1.0) < 1e-10

    def test_get_trigonometric_invalide_raises(self, projet_vide):
        with pytest.raises(ValueError):
            projet_vide.get_trigonometric_value(45, "SINH")

    def test_convert_unit_physical_m_to_mm(self):
        result = MathUtilsMixin._convert_unit_physical(1.5, si.m, si.mm)
        assert abs(result - 1500.0) < 1e-9

    def test_convert_unit_physical_N_to_kN(self):
        result = MathUtilsMixin._convert_unit_physical(5000, si.N, si.kN)
        assert abs(result - 5.0) < 1e-9

    def test_convert_unit_physical_m2_to_mm2(self):
        result = MathUtilsMixin._convert_unit_physical(1.0, si.m**2, si.mm**2)
        assert abs(result - 1e6) < 1e-6


# ---------------------------------------------------------------------------
# Test_SyntheseMixin
# ---------------------------------------------------------------------------

class Test_SyntheseMixin:
    @pytest.fixture
    def obj(self):
        return Projet()

    def test_add_synthese_dict_format(self, obj):
        lignes = [{"effort": "Traction", "SN": 0.5, "SI": 0.3}]
        df = obj._add_synthese_taux_travail(lignes)
        assert len(df) == 1
        assert df.iloc[0, 1] == 0.5

    def test_add_synthese_tuple_format(self, obj):
        lignes = [("Compression", 0.8, 0.6)]
        df = obj._add_synthese_taux_travail(lignes)
        assert df.iloc[0]["Situation normale"] == 0.8

    def test_add_synthese_cumul(self, obj):
        obj._add_synthese_taux_travail([("Traction", 0.4, 0.2)])
        df = obj._add_synthese_taux_travail([("Cisaillement", 0.7, 0.5)])
        assert len(df) == 2

    def test_synthese_taux_travail_max(self, obj):
        obj._add_synthese_taux_travail([
            ("Traction", 0.4, 0.2),
            ("Flexion", 0.9, 0.7),
        ])
        df_final = obj.synthese_taux_travail()
        max_row = df_final[df_final.iloc[:, 0] == "Max"].iloc[0]
        assert max_row["Situation normale"] == 0.9
        assert max_row["Situation d'incendie"] == 0.7

    def test_synthese_type_invalide_raises(self, obj):
        with pytest.raises(TypeError):
            obj._add_synthese_taux_travail(["invalide"])


# ---------------------------------------------------------------------------
# Test_DataLoaderMixin
# ---------------------------------------------------------------------------

class Test_DataLoaderMixin:
    def test_data_from_csv_retourne_dataframe(self, projet_vide):
        df = projet_vide._data_from_csv("caracteristique_meca_bois.csv")
        assert not df.empty
        assert "C24" in df.index or any("C24" in str(i) for i in df.index)

    def test_data_from_csv_cache(self, projet_vide):
        df1 = projet_vide._data_from_csv("caracteristique_meca_acier.csv")
        df2 = projet_vide._data_from_csv("caracteristique_meca_acier.csv")
        assert df1.equals(df2)

    def test_data_from_csv_index_correct(self, projet_vide):
        df = projet_vide._data_from_csv("caracteristique_meca_acier.csv")
        assert "S235" in df.index


if __name__ == "__main__":
    pytest.main(["-v", __file__])
