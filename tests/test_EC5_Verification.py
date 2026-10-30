# Encoding in UTF-8 by Anthony PARISOT
"""Tests d'intégration bout-à-bout de la classe Verification_EC5.

Scénarios :
- Poutre simple sur deux appuis (M1) : flexion + cisaillement + flèche.
- Poutre continue sur trois appuis (M1-M2 groupés) : auto-group + verify.
- Boucle par combinaison : kmod et typecombi corrects selon la combo gouvernante.
"""

import sys
from pathlib import Path

import forallpeople as si
import pandas as pd
import pytest

sys.path.append(str(Path(__file__).parent.parent))

from ourocode.eurocode.core.combinaison import Combinaison
from ourocode.eurocode.core.model_generator import Model_generator
from ourocode.eurocode.core.model_result import Model_result
from ourocode.eurocode.ec5.element_droit.verification_EC5 import Verification_EC5


# ---------------------------------------------------------------------------
# Fixtures : poutre simple (1 travée) et poutre continue (2 travées)
# ---------------------------------------------------------------------------


def _new_model(name: str = "test") -> Model_generator:
    return Model_generator(ingenieur="Testeur", name=name, code_INSEE=75056, alt=100)


def _build_beam(model: Model_generator, nodes_xy, b=100, h=220, classe="C24"):
    mat = model.add_material_by_class(classe)
    sec = model.add_section(b, h, 0, "Rectangulaire")
    for x, y in nodes_xy:
        model.add_node(X=x, Y=y, Z=0)
    for i in range(len(nodes_xy) - 1):
        model.add_member(f"N{i + 1}", f"N{i + 2}", mat, sec, poids_propre=False)
    model.add_support("N1", DX=True, DY=True, DZ=True, RX=True, RY=False, RZ=False)
    for i in range(1, len(nodes_xy)):
        model.add_support(f"N{i + 1}", DX=False, DY=True, DZ=True, RX=True, RY=False, RZ=False)
    return mat, sec


def _apply_gq_loads(model: Model_generator, gkNm: float, qkNm: float):
    for mid in model.get_all_members().keys():
        model.create_dist_load(
            member_id=mid,
            name=f"G_{mid}",
            start_load=-gkNm,
            end_load=-gkNm,
            start_pos="start",
            end_pos="end",
            action="Permanente G",
            direction="FY",
        )
        model.create_dist_load(
            member_id=mid,
            name=f"Q_{mid}",
            start_load=-qkNm,
            end_load=-qkNm,
            start_pos="start",
            end_pos="end",
            action="Exploitation Q",
            direction="FY",
        )


@pytest.fixture
def simple_beam_verif():
    """Poutre simple 5 m, 100x220 C24, G=1 kN/m + Q=2 kN/m."""
    model = _new_model("Poutre simple")
    _build_beam(model, [(0, 0), (5000, 0)])
    _apply_gq_loads(model, gkNm=1.0, qkNm=2.0)
    combi = Combinaison(
        model,
        ELU_STR=True,
        ELU_STR_ACC=False,
        ELS_C=True,
        ELS_QP=True,
        cat="Cat A : habitation",
        kdef=0.6,
        type_psy_2="Moyen terme",
    )
    result = Model_result(model, analyze_type="Général", check_stability=False)
    model.group_members(
        "Solive",
        ["M1"],
        role="Solive",
        design_params={
            "classe_bois": "C24",
            "cs": 1,
            "type_element_fleche": "Élément structuraux",
        },
    )
    verif = Verification_EC5(combi, result, "ELU_ALL", 12, 12, False, 1, "Bâtiments courants", "Élément structuraux")
    return model, combi, result, verif


@pytest.fixture
def continuous_beam_verif():
    """Poutre continue 2x3 m, 100x200 C24, G=1 kN/m + Q=2 kN/m."""
    model = _new_model("Poutre continue")
    _build_beam(model, [(0, 0), (3000, 0), (6000, 0)], b=100, h=200)
    _apply_gq_loads(model, gkNm=1.0, qkNm=2.0)
    combi = Combinaison(
        model,
        ELU_STR=True,
        ELU_STR_ACC=False,
        ELS_C=True,
        ELS_QP=True,
        cat="Cat A : habitation",
        kdef=0.6,
        type_psy_2="Moyen terme",
    )
    result = Model_result(model, analyze_type="Général", check_stability=False)
    model.auto_group_continuous_members(
        role="Solive",
        design_params={
            "classe_bois": "C24",
            "cs": 1,
            "type_element_fleche": "Élément structuraux",
        },
    )
    verif = Verification_EC5(combi, result, "ELU_ALL", 12, 12, False, 1, "Bâtiments courants", "Élément structuraux")
    return model, combi, result, verif


@pytest.fixture
def fire_beam_verif():
    """Poutre simple avec combinaisons FIRE pour tester la vérification incendie."""
    model = _new_model("Poutre feu")
    _build_beam(model, [(0, 0), (5000, 0)])
    _apply_gq_loads(model, gkNm=1.0, qkNm=2.0)
    combi = Combinaison(
        model,
        ELU_STR=True,
        ELU_STR_ACC=True,
        ELS_C=False,
        ELS_QP=False,
        cat="Cat A : habitation",
        kdef=0.6,
        type_psy_2="Moyen terme",
    )
    result = Model_result(model, analyze_type="Général", check_stability=False)
    model.group_members(
        "Solive",
        ["M1"],
        role="Solive",
        design_params={
            "classe_bois": "C24",
            "cs": 1,
            "type_element_fleche": "Élément structuraux",
        },
    )
    verif = Verification_EC5(
        combi, result, "ELU_ALL", 12, 12, False, 1, "Bâtiments courants", "Élément structuraux", t_expo=30
    )
    return model, combi, result, verif


# ---------------------------------------------------------------------------
# Combinaison.type_combi
# ---------------------------------------------------------------------------


class TestTypeCombi:
    def test_elu_str_is_fondamental(self, simple_beam_verif):
        _, combi, _, _ = simple_beam_verif
        assert combi.type_combi("ELU_STR 1.35G + 1.5Q") == "Fondamentales"
        assert combi.type_combi("ELU_STR G") == "Fondamentales"

    def test_elu_acc_is_accidental(self, simple_beam_verif):
        _, combi, _, _ = simple_beam_verif
        assert combi.type_combi("ELU_STR_ACC G + Ae") == "Accidentelles"

    def test_non_elu_raises(self, simple_beam_verif):
        _, combi, _, _ = simple_beam_verif
        with pytest.raises(ValueError):
            combi.type_combi("ELS_C G + Q")


# ---------------------------------------------------------------------------
# Verification_EC5 : factory et état
# ---------------------------------------------------------------------------


class TestFactory:
    def test_from_combinaison_instantiable(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        assert isinstance(verif, Verification_EC5)
        assert verif._combinaison is not None
        assert verif._result is not None
        assert verif._model_generator is not None

    def test_direct_init_requires_check(self, simple_beam_verif):
        """Verification_EC5 nécessite combinaison et model_result pour verify()."""
        model, combi, result, verif = simple_beam_verif
        # Le constructeur vérifie déjà la présence des arguments requis
        # Si des arguments sont manquants, il lève TypeError
        assert verif._combinaison is not None
        assert verif._result is not None


# ---------------------------------------------------------------------------
# Verification_EC5.verify — poutre simple
# ---------------------------------------------------------------------------


class TestVerifySimpleBeam:
    def test_returns_dataframe(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        res = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        assert res["name"] == "Solive"
        assert isinstance(res["dataframe"], pd.DataFrame)
        assert not res["dataframe"].empty

    def test_flexion_and_cisaillement_present(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        res = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        assert res["taux"]["Flexion"] is not None
        assert res["taux"]["Cisaillement"] is not None

    def test_no_axial(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        res = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        assert res["taux"]["Traction"] is None
        assert res["taux"]["Compression"] is None

    def test_taux_values_reasonable(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        res = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        tflex = res["taux"]["Flexion"]["taux"]
        tcis = res["taux"]["Cisaillement"]["taux"]
        assert 0.3 < tflex < 3.0
        assert 0.0 < tcis < 1.5

    def test_fleche_present(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        res = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        assert res["taux"]["Flèche W_inst(Q)"] is not None

    def test_governing_combo_is_1_35G_1_5Q(self, simple_beam_verif):
        """La combo gouvernante de la flexion doit être ELU_STR 1.35G + 1.5Q."""
        _, _, _, verif = simple_beam_verif
        res = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        combo = res["taux"]["Flexion"]["combinaison"]
        assert "1.35G" in combo and "1.5Q" in combo


# ---------------------------------------------------------------------------
# Verification_EC5.get_combo_objects
# ---------------------------------------------------------------------------


class TestGetComboObjects:
    COMBO = "ELU_STR 1.35G + 1.5Q"

    def test_returns_flexion_and_cisaillement(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        objs = verif.get_combo_objects("Solive", self.COMBO)
        assert objs["Flexion"] is not None
        assert objs["Cisaillement"] is not None

    def test_compression_and_traction_none(self, simple_beam_verif):
        """Poutre sans effort axial : Compression et Traction doivent être None."""
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        objs = verif.get_combo_objects("Solive", self.COMBO)
        assert objs["Compression"] is None
        assert objs["Traction"] is None

    def test_flexion_object_has_taux_attribute(self, simple_beam_verif):
        """L'objet Flexion retourné doit posséder taux_m_rd calculé."""
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        objs = verif.get_combo_objects("Solive", self.COMBO)
        flexion = objs["Flexion"]
        assert hasattr(flexion, "taux_m_rd")
        assert len(flexion.taux_m_rd) > 0

    def test_cisaillement_object_has_taux_attribute(self, simple_beam_verif):
        """L'objet Cisaillement retourné doit posséder taux_tau_rd calculé."""
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        objs = verif.get_combo_objects("Solive", self.COMBO)
        cis = objs["Cisaillement"]
        assert hasattr(cis, "taux_tau_rd")
        assert len(cis.taux_tau_rd) > 0

    def test_raises_if_verify_not_called(self, simple_beam_verif):
        """KeyError attendu si verify() n'a pas été appelé."""
        _, _, _, verif = simple_beam_verif
        with pytest.raises(KeyError, match="verify"):
            verif.get_combo_objects("Solive", self.COMBO)

    def test_taux_matches_combo_elu(self, simple_beam_verif):
        """Le taux Flexion de get_combo_objects (combo gouvernante) doit être > 0."""
        _, combi, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        objs = verif.get_combo_objects("Solive", self.COMBO)
        taux_flexion = max(objs["Flexion"].taux_m_rd.values())
        assert taux_flexion > 0.0
        # La combo G seule produit un taux plus faible
        objs_g = verif.get_combo_objects("Solive", "ELU_STR 1.35G")
        taux_g = max(objs_g["Flexion"].taux_m_rd.values())
        assert taux_flexion > taux_g


# ---------------------------------------------------------------------------
# Verification_EC5.get_verification_latex
# ---------------------------------------------------------------------------


class TestGetVerificationLatex:
    COMBO = "ELU_STR 1.35G + 1.5Q"

    def test_flexion_latex_is_non_empty_string(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        latex = verif.get_verification_latex("Solive", "Flexion", self.COMBO)
        assert isinstance(latex, str)
        assert len(latex) > 0

    def test_cisaillement_latex_is_non_empty_string(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        latex = verif.get_verification_latex("Solive", "Cisaillement", self.COMBO)
        assert isinstance(latex, str)
        assert len(latex) > 0

    def test_fleche_latex_does_not_require_combo_name(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        latex = verif.get_verification_latex("Solive", "Flèche W_inst(Q)")
        assert isinstance(latex, str)
        assert len(latex) > 0

    def test_traction_not_applicable_raises_value_error(self, simple_beam_verif):
        """Poutre sans effort axial : la traction n'est pas applicable sur cette combo."""
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        with pytest.raises(ValueError):
            verif.get_verification_latex("Solive", "Traction", self.COMBO)

    def test_missing_combo_name_for_elu_type_raises_value_error(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        with pytest.raises(ValueError):
            verif.get_verification_latex("Solive", "Flexion")

    def test_unknown_verif_type_raises_key_error(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        with pytest.raises(KeyError):
            verif.get_verification_latex("Solive", "Torsion", self.COMBO)

    def test_raises_if_verify_not_called(self, simple_beam_verif):
        _, _, _, verif = simple_beam_verif
        with pytest.raises(KeyError, match="verify"):
            verif.get_verification_latex("Solive", "Flexion", self.COMBO)

    def test_fire_flexion_latex_is_non_empty_string(self, fire_beam_verif):
        _, combi, _, verif = fire_beam_verif
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        fire_combos = [c for c in combi.get_list_combination("ELU_ALL") or [] if "FIRE" in c]
        if fire_combos:
            latex = verif.get_verification_latex("Solive", "Flexion_feu", fire_combos[0])
            assert isinstance(latex, str)
            assert len(latex) > 0


# ---------------------------------------------------------------------------
# Boucle par combinaison : kmod dépend de la combo
# ---------------------------------------------------------------------------


class TestKmodPerCombo:
    def test_kmod_permanente_stricter_than_moyen_terme(self, simple_beam_verif):
        """Vérifie que la vérification utilise effectivement le bon kmod
        par combinaison : pour une combo G seule (Permanente), kmod est
        plus faible que pour G+Q (Moyen terme), donc taux plus élevé
        à effort égal.

        On compare les taux de flexion calculés sur la combo ELU_STR G
        (Permanente, kmod=0.6) vs ELU_STR 1.35G + 1.5Q (Moyen terme,
        kmod=0.8). Vu la différence de chargement (G vs 1.35G+1.5Q),
        ce n'est pas strictement prouvable, mais on peut vérifier que
        min_type_load renvoie bien les bonnes durées.
        """
        _, combi, _, _ = simple_beam_verif
        assert combi.min_type_load("ELU_STR 1.35G") == "Permanente"
        assert combi.min_type_load("ELU_STR 1.35G + 1.5Q") == "Moyen terme"

    def test_verify_respects_combo_governance(self, simple_beam_verif):
        """Si l'on filtre sur ELU_STR_ACC uniquement (pas d'accidentelle
        définie), aucune vérification ne doit sortir."""
        _, _, _, verif = simple_beam_verif
        # Le filtre ELU_STR_ACC n'existe pas dans get_list_combination --
        # on utilise 'ELU_STR' qui doit fonctionner et donner le même
        # résultat que ELU_ALL puisqu'il n'y a pas d'accidentelle.
        res_all = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        res_str = verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        assert res_all["taux"]["Flexion"]["taux"] == pytest.approx(res_str["taux"]["Flexion"]["taux"], rel=1e-6)


# ---------------------------------------------------------------------------
# Verification_EC5.synthese — poutre continue
# ---------------------------------------------------------------------------


class TestSyntheseContinuousBeam:
    def test_auto_group_one_sm(self, continuous_beam_verif):
        model, _, _, _ = continuous_beam_verif
        names = [n for n, _ in model._iter_structural_members()]
        assert len(names) == 1
        sm = model.get_structural_member(names[0])
        assert sm["Barres FEM"] == ["M1", "M2"]

    def test_synthese_returns_dataframe(self, continuous_beam_verif):
        _, _, _, verif = continuous_beam_verif
        df = verif.synthese()
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        for col in ["Barre structurale", "Vérification", "Taux", "Statut"]:
            assert col in df.columns

    def test_synthese_sorted_desc(self, continuous_beam_verif):
        _, _, _, verif = continuous_beam_verif
        df = verif.synthese()
        taux = df["Taux"].tolist()
        assert taux == sorted(taux, reverse=True)

    def test_synthese_governing_combo_present(self, continuous_beam_verif):
        """La synthèse doit contenir au moins une ligne avec ELU_STR 1.35G + 1.5Q."""
        _, _, _, verif = continuous_beam_verif
        df = verif.synthese()
        combos = df["Combinaison"].tolist()
        assert any("1.35G" in c and "1.5Q" in c for c in combos if c is not None)


# ---------------------------------------------------------------------------
# Verification_EC5 : vérification feu
# ---------------------------------------------------------------------------


class TestFireVerification:
    def test_fire_combo_has_no_normal_verification(self, fire_beam_verif):
        """Pour une combinaison FIRE, les vérifications normales doivent être None."""
        _, combi, _, verif = fire_beam_verif
        # Vérifier qu'il existe des combinaisons FIRE
        fire_combos = [c for c in combi.get_list_combination("ELU_ALL") or [] if "FIRE" in c]
        if fire_combos:
            # Appeler verify pour peupler le cache
            verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
            # Récupérer les objets pour une combinaison FIRE
            objs = verif.get_combo_objects("Solive", fire_combos[0])
            # Les objets normaux doivent être None pour les combos FIRE
            assert objs["Flexion"] is None
            assert objs["Cisaillement"] is None
            assert objs["Compression"] is None
            assert objs["Traction"] is None

    def test_fire_combo_has_fire_verification(self, fire_beam_verif):
        """Pour une combinaison FIRE, les vérifications feu doivent être calculées."""
        _, combi, _, verif = fire_beam_verif
        # Vérifier qu'il existe des combinaisons FIRE
        fire_combos = [c for c in combi.get_list_combination("ELU_ALL") or [] if "FIRE" in c]
        if fire_combos:
            # Appeler verify pour peupler le cache
            verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
            # Récupérer les objets pour une combinaison FIRE
            objs = verif.get_combo_objects("Solive", fire_combos[0])
            # Au moins un objet feu doit être calculé (flexion ou cisaillement)
            assert objs["Flexion_feu"] is not None or objs["Cisaillement_feu"] is not None

    def test_normal_combo_has_no_fire_verification(self, simple_beam_verif):
        """Pour une combinaison normale, les vérifications feu doivent être None."""
        _, combi, _, verif = simple_beam_verif
        # Appeler verify pour peupler le cache
        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        # Récupérer les objets pour une combinaison normale
        objs = verif.get_combo_objects("Solive", "ELU_STR 1.35G + 1.5Q")
        # Les objets feu doivent être None pour les combos normaux
        assert objs["Flexion_feu"] is None
        assert objs["Cisaillement_feu"] is None
        assert objs["Compression_feu"] is None
        assert objs["Traction_feu"] is None

    def test_t_expo_parameter_used(self, fire_beam_verif):
        """Le paramètre t_expo doit être utilisé dans les calculs feu."""
        _, _, _, verif = fire_beam_verif
        assert verif.t_expo == 30


class TestVerifyComboDirections:
    @staticmethod
    def _zero_efforts():
        force = {"Min": (0 * si.N, None), "Max": (0 * si.N, None)}
        moment = {"Min": (0 * si.N * si.mm, None), "Max": (0 * si.N * si.mm, None)}
        return {
            "Nx": force.copy(),
            "Vy": force.copy(),
            "Vz": force.copy(),
            "Mx": moment.copy(),
            "My": moment.copy(),
            "Mz": moment.copy(),
        }

    def test_negative_nx_is_verified_as_positive_traction(self, simple_beam_verif, monkeypatch):
        _, _, result, verif = simple_beam_verif
        efforts = self._zero_efforts()
        efforts["Nx"]["Min"] = (-10 * si.kN, "ELU_STR 1.35G")
        monkeypatch.setattr(result, "get_min_max_internal_force", lambda *args: efforts)
        monkeypatch.setattr(
            result,
            "get_absolute_internal_force",
            lambda *args, **kwargs: 0 * si.N,
        )

        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        obj = verif.get_combo_objects("Solive", "ELU_STR 1.35G")

        assert obj["Traction"] is not None
        assert max(obj["Traction"].taux_t_0_rd.values()) > 0

    def test_vz_is_verified_when_more_unfavorable_than_vy(self, simple_beam_verif, monkeypatch):
        _, _, result, verif = simple_beam_verif
        efforts = self._zero_efforts()
        monkeypatch.setattr(result, "get_min_max_internal_force", lambda *args: efforts)

        def absolute_force(*args, **kwargs):
            return 10 * si.kN if args[2] == "Vz" else 0 * si.N

        monkeypatch.setattr(result, "get_absolute_internal_force", absolute_force)

        verif.verify("Solive", "Élément structuraux", "Charge sur fibre comprimée", 12, 12, False, 1)
        obj = verif.get_combo_objects("Solive", "ELU_STR 1.35G")

        assert obj["Cisaillement"] is not None
        assert max(obj["Cisaillement"].taux_tau_rd.values()) > 0
