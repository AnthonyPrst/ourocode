# Encoding in UTF-8 by Anthony PARISOT
# Tests pour la détection et le regroupement de barres structurales
# dans Model_generator.
import sys
from pathlib import Path

import pytest

sys.path.append(str(Path(__file__).parent.parent))

from ourocode.eurocode.core.model_generator import Model_generator


def _make_model() -> Model_generator:
    return Model_generator(
        ingenieur="Testeur", name="Projet Test", code_INSEE=75056, alt=100
    )


def _make_colinear_chain(nb_members: int, material: str = "C24", section=None):
    """Crée une chaîne colinéaire de ``nb_members`` barres alignées sur l'axe X.

    Returns:
        tuple[Model_generator, str, list[str]]: (model, material_id, member_ids)
    """
    model = _make_model()
    mat = model.add_material_by_class(material)
    sec = section if section else model.add_section(100, 200, 0, "Rectangulaire")
    # Création des noeuds N1..N{nb+1} alignés sur X
    for i in range(nb_members + 1):
        model.add_node(X=1000 * i, Y=0, Z=0)
    member_ids = []
    for i in range(nb_members):
        mid = model.add_member(
            f"N{i + 1}", f"N{i + 2}", mat, sec, poids_propre=False
        )
        member_ids.append(mid)
    return model, mat, member_ids


# ---------------------------------------------------------------------------
# detect_continuous_members
# ---------------------------------------------------------------------------

class TestDetectContinuousMembers:
    def test_single_member(self):
        model, _, ids = _make_colinear_chain(1)
        chains = model.detect_continuous_members()
        assert chains == [ids]

    def test_continuous_chain_of_three(self):
        model, _, ids = _make_colinear_chain(3)
        chains = model.detect_continuous_members()
        assert len(chains) == 1
        assert chains[0] == ids

    def test_break_on_different_section(self):
        model = _make_model()
        mat = model.add_material_by_class("C24")
        sec1 = model.add_section(100, 200, 0, "Rectangulaire")
        sec2 = model.add_section(150, 250, 0, "Rectangulaire")
        for i in range(3):
            model.add_node(X=1000 * i, Y=0, Z=0)
        m1 = model.add_member("N1", "N2", mat, sec1, poids_propre=False)
        m2 = model.add_member("N2", "N3", mat, sec2, poids_propre=False)
        chains = model.detect_continuous_members()
        # Les deux barres ne doivent pas fusionner (section différente)
        assert sorted(map(tuple, chains)) == sorted([(m1,), (m2,)])

    def test_break_on_different_material(self):
        model = _make_model()
        mat1 = model.add_material_by_class("C24")
        mat2 = model.add_material_by_class("GL24h")
        sec = model.add_section(100, 200, 0, "Rectangulaire")
        for i in range(3):
            model.add_node(X=1000 * i, Y=0, Z=0)
        model.add_member("N1", "N2", mat1, sec, poids_propre=False)
        model.add_member("N2", "N3", mat2, sec, poids_propre=False)
        chains = model.detect_continuous_members()
        assert all(len(c) == 1 for c in chains)
        assert len(chains) == 2

    def test_break_on_hinge_release(self):
        model, _, ids = _make_colinear_chain(3)
        # Rotule en rotation z à la fin de M1 → rompt la continuité M1/M2
        model.add_release(
            ids[0], position="end",
            u=False, v=False, w=False,
            teta_x=False, teta_y=False, teta_z=True,
        )
        chains = model.detect_continuous_members()
        chains_sorted = sorted(map(tuple, chains), key=lambda c: c[0])
        assert chains_sorted == [(ids[0],), (ids[1], ids[2])]

    def test_t_junction_perpendicular_branch_does_not_break_chain(self):
        """Une contrefiche perpendiculaire fixée en milieu d'arbalétrier ne doit pas
        rompre la continuité de l'arbalétrier (les 2 barres de l'arbalétrier sont
        colinéaires — la non-colinéarité de la contrefiche l'exclut naturellement)."""
        model, mat, ids = _make_colinear_chain(2)
        sec = list(model.get_all_sections().keys())[0]
        # Contrefiche perpendiculaire fixée au noeud intermédiaire N2
        model.add_node(X=1000, Y=1000, Z=0)  # N4
        contrefiche = model.add_member("N2", "N4", mat, sec, poids_propre=False)
        chains = model.detect_continuous_members()
        # L'arbalétrier (ids[0]+ids[1]) forme une chaîne continue ;
        # la contrefiche reste isolée.
        assert len(chains) == 2
        arba_chain = next(c for c in chains if len(c) == 2)
        assert set(arba_chain) == set(ids)
        contrefiche_chain = next(c for c in chains if len(c) == 1)
        assert contrefiche_chain == [contrefiche]

    def test_t_junction_with_release_breaks_chain(self):
        """Si la contrefiche a un relâchement en rotation à son noeud de fixation
        sur l'arbalétrier, la continuité de l'arbalétrier reste intacte (le
        relâchement ne concerne que la contrefiche)."""
        model, mat, ids = _make_colinear_chain(2)
        sec = list(model.get_all_sections().keys())[0]
        model.add_node(X=1000, Y=1000, Z=0)  # N4
        contrefiche = model.add_member("N2", "N4", mat, sec, poids_propre=False)
        # Relâchement en rotation au pied de la contrefiche (noeud N2)
        model.add_release(
            contrefiche, position="start",
            u=False, v=False, w=False,
            teta_x=False, teta_y=False, teta_z=True,
        )
        chains = model.detect_continuous_members()
        # L'arbalétrier reste continu ; la contrefiche isolée.
        assert len(chains) == 2
        arba_chain = next(c for c in chains if len(c) == 2)
        assert set(arba_chain) == set(ids)

    def test_break_on_non_colinearity(self):
        model = _make_model()
        mat = model.add_material_by_class("C24")
        sec = model.add_section(100, 200, 0, "Rectangulaire")
        model.add_node(X=0, Y=0, Z=0)         # N1
        model.add_node(X=1000, Y=0, Z=0)      # N2
        model.add_node(X=1000, Y=1000, Z=0)   # N3 -> angle 90°
        m1 = model.add_member("N1", "N2", mat, sec, poids_propre=False)
        m2 = model.add_member("N2", "N3", mat, sec, poids_propre=False)
        chains = model.detect_continuous_members()
        assert sorted(map(tuple, chains)) == sorted([(m1,), (m2,)])

    def test_angle_tolerance(self):
        """Un léger désalignement est accepté si dans la tolérance."""
        model = _make_model()
        mat = model.add_material_by_class("C24")
        sec = model.add_section(100, 200, 0, "Rectangulaire")
        model.add_node(X=0, Y=0, Z=0)         # N1
        model.add_node(X=1000, Y=0, Z=0)      # N2
        # Angle de ~0.5° (dy=10 sur dx=1000 -> atan ≈ 0.57°)
        model.add_node(X=2000, Y=10, Z=0)     # N3
        ids = [
            model.add_member("N1", "N2", mat, sec, poids_propre=False),
            model.add_member("N2", "N3", mat, sec, poids_propre=False),
        ]
        # Tolérance 1° : doivent fusionner
        assert model.detect_continuous_members(angle_tol_deg=1.0) == [ids]
        # Tolérance 0.1° : ne doivent pas fusionner
        chains_strict = model.detect_continuous_members(angle_tol_deg=0.1)
        assert all(len(c) == 1 for c in chains_strict)

    def test_reversed_member_orientation(self):
        """Deux barres partageant un noeud via leurs extrémités "start" doivent aussi fusionner."""
        model = _make_model()
        mat = model.add_material_by_class("C24")
        sec = model.add_section(100, 200, 0, "Rectangulaire")
        model.add_node(X=-1000, Y=0, Z=0)  # N1
        model.add_node(X=0, Y=0, Z=0)      # N2 (partagé, deux "start")
        model.add_node(X=1000, Y=0, Z=0)   # N3
        # Les deux barres partent de N2 -> même extrémité "start"
        m1 = model.add_member("N2", "N1", mat, sec, poids_propre=False)
        m2 = model.add_member("N2", "N3", mat, sec, poids_propre=False)
        chains = model.detect_continuous_members()
        assert len(chains) == 1
        assert set(chains[0]) == {m1, m2}


# ---------------------------------------------------------------------------
# group_members / auto_group_continuous_members / iter / get / del
# ---------------------------------------------------------------------------

class TestGroupMembers:
    def test_group_and_get(self):
        model, _, ids = _make_colinear_chain(3)
        model.group_members(
            "Solive_1", ids, role="Solive",
            design_params={"classe_bois": "C24", "cs": 1, "Hi": 12, "Hf": 12},
        )
        sm = model.get_structural_member("Solive_1")
        assert sm["Barres FEM"] == ids
        assert sm["Rôle"] == "Solive"
        assert sm["Design"]["classe_bois"] == "C24"
        # Vérification des types d'appuis automatiques
        assert "type_appuis_y" in sm["Design"]
        assert "type_appuis_z" in sm["Design"]
        # Sans appuis définis, les rotations sont libres des deux côtés
        assert sm["Design"]["type_appuis_y"] == "Rotule - Rotule"
        assert sm["Design"]["type_appuis_z"] == "Rotule - Rotule"
        # Vérification des paramètres de flexion automatiques (déversement)
        assert "lo_rel_y" in sm["Design"]
        assert "lo_rel_z" in sm["Design"]
        assert "coeflef_y" in sm["Design"]
        assert "coeflef_z" in sm["Design"]
        # Sans charges définies, coeflef = 1.0 (moment constant)
        assert sm["Design"]["coeflef_y"] == 1.0
        assert sm["Design"]["coeflef_z"] == 1.0
        # Longueur totale = somme des 3 barres de 1000 mm
        assert abs(sm["Longueur totale"].value - 3.0) < 1e-6  # en mètres (si.m)

    def test_empty_member_ids_raises(self):
        model = _make_model()
        with pytest.raises(ValueError, match="ne peut pas être vide"):
            model.group_members("X", [])

    def test_unknown_member_raises(self):
        model, _, ids = _make_colinear_chain(1)
        with pytest.raises(ValueError, match="n'existe pas"):
            model.group_members("X", ids + ["M999"])

    def test_duplicate_name_raises(self):
        model, _, ids = _make_colinear_chain(1)
        model.group_members("X", ids)
        with pytest.raises(ValueError, match="existe déjà"):
            model.group_members("X", ids)

    def test_del_structural_member(self):
        model, _, ids = _make_colinear_chain(1)
        model.group_members("X", ids, role="Poteau")
        data = model.del_structural_member("X")
        assert data["Rôle"] == "Poteau"
        assert "X" not in model.get_all_structural_members()

    def test_iter_structural_members(self):
        model, _, ids = _make_colinear_chain(2)
        model.group_members("A", [ids[0]])
        model.group_members("B", [ids[1]])
        names = [name for name, _ in model._iter_structural_members()]
        assert names == ["A", "B"]

    def test_auto_group_creates_all(self):
        model, _, ids = _make_colinear_chain(3)
        created = model.auto_group_continuous_members(
            role="Solive", design_params={"classe_bois": "C24", "cs": 1},
        )
        assert created == ["SM1"]
        sm = model.get_structural_member("SM1")
        assert sm["Barres FEM"] == ids
        assert sm["Rôle"] == "Solive"

    def test_auto_group_only_continuous(self):
        """Une barre isolée doit être ignorée avec only_continuous=True."""
        # Chaîne continue M1-M2 + barre isolée M3 (autre section -> rupture)
        model = _make_model()
        mat = model.add_material_by_class("C24")
        sec1 = model.add_section(100, 200, 0, "Rectangulaire")
        sec2 = model.add_section(150, 250, 0, "Rectangulaire")
        for i in range(3):
            model.add_node(X=1000 * i, Y=0, Z=0)
        model.add_member("N1", "N2", mat, sec1, poids_propre=False)
        model.add_member("N2", "N3", mat, sec2, poids_propre=False)
        created = model.auto_group_continuous_members(only_continuous=True)
        assert created == []  # aucune chaîne de longueur ≥ 2

    def test_lo_rel_with_lateral_bracing(self):
        """Appui DZ intermédiaire en N2 → lo_rel_y = 1000 mm (demi-portée)."""
        model, _, ids = _make_colinear_chain(2)  # N1-N2-N3, chaque barre = 1000 mm
        # Appui simple en N1 et N3 (bloquant DY et DZ → appuis d'about)
        model.add_support("N1", DX=False, DY=True, DZ=True, RX=False, RY=False, RZ=False)
        model.add_support("N3", DX=False, DY=True, DZ=True, RX=False, RY=False, RZ=False)
        # Appui latéral intermédiaire en N2 : bloque DZ uniquement (contreventement Z)
        model.add_support("N2", DX=False, DY=False, DZ=True, RX=False, RY=False, RZ=False)
        model.group_members("Poutre", ids)
        sm = model.get_structural_member("Poutre")
        # lo_rel_y (déversement selon Y) bloqué par DZ : segments de 1000 mm
        assert sm["Design"]["lo_rel_y"] == pytest.approx(1000.0)
        # lo_rel_z (déversement selon Z) bloqué par DY : N2 ne bloque pas DY
        # → longueur totale = 2000 mm
        assert sm["Design"]["lo_rel_z"] == pytest.approx(2000.0)
