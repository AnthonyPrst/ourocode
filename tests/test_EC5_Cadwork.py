# Encoding in UTF-8 by Anthony PARISOT
"""Tests de la passerelle Cadwork -> MOB (export JSON + MOB.load_cadwork_data)."""

import json
import math
import sys
import types
from pathlib import Path

import forallpeople as si
import numpy as np
import pytest

si.environment("structural")
from ourocode.eurocode.core.batiment import Batiment
from ourocode.eurocode.ec5 import cvt as EC5
from ourocode.eurocode.ec5.assemblage import Pointe
from ourocode.eurocode.ec5.cadwork import ourocode_export_mob as cw
from ourocode.eurocode.ec5.element_droit.barre import Barre


def _rot(p, angle):
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    return (p[0] * c - p[1] * s, p[0] * s + p[1] * c, p[2])


def box(x0, x1, y0, y1, z0, z1, angle=0):
    pts = [(x, y, z) for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)]
    return [_rot(p, angle) for p in pts]


def model(angle=0):
    """Bâtiment 10 x 8 m, parois de 200 mm d'axe x=0, x=10, y=0, y=8, OSB 12 mm."""
    H, t, e = 2700, 100, 12
    xl = {"x": _rot((1, 0, 0), angle), "y": _rot((0, 1, 0), angle)}
    walls = [
        {"id": 1, "name": "Sud", "storey": "RDC", "xl": xl["x"], "vertices": box(-t, 10000 + t, -t, t, 0, H, angle)},
        {
            "id": 2,
            "name": "Nord",
            "storey": "RDC",
            "xl": xl["x"],
            "vertices": box(-t, 10000 + t, 8000 - t, 8000 + t, 0, H, angle),
        },
        {"id": 3, "name": "Ouest", "storey": "RDC", "xl": xl["y"], "vertices": box(-t, t, -t, 8000 + t, 0, H, angle)},
        {
            "id": 4,
            "name": "Est",
            "storey": "RDC",
            "xl": xl["y"],
            "vertices": box(10000 - t, 10000 + t, -t, 8000 + t, 0, H, angle),
        },
    ]
    openings = [{"id": 50, "name": "Fenêtre", "vertices": box(4000, 5200, -t, t, 900, 2100, angle)}]

    def osb(i, x0, x1, y0, y1, z0=0, z1=H, name="OSB 12"):
        return {
            "id": i,
            "name": name,
            "material": "OSB/3",
            "storey": "RDC",
            "vertices": box(x0, x1, y0, y1, z0, z1, angle),
        }

    panels = []
    # Sud (extérieur y < -100) : 3 x 1250 à gauche de la fenêtre, 4 à droite, 1 au-dessus de la fenêtre
    for k in range(3):
        panels.append(osb(100 + k, 1250 * k, 1250 * (k + 1), -t - e, -t))
    for k, (a, b) in enumerate([(5200, 6450), (6450, 7700), (7700, 8950), (8950, 10000)]):
        panels.append(osb(110 + k, a, b, -t - e, -t))
    panels.append(osb(120, 4000, 5200, -t - e, -t, 2100, H))
    # Isolant : filtré par les mots-clés
    panels.append(
        {"id": 130, "name": "Laine de bois", "material": "Fibre", "vertices": box(0, 1250, -t, -t + 60, 0, H, angle)}
    )
    # Nord : 8 x 1250 face intérieure (y < 7900)
    for k in range(8):
        panels.append(osb(200 + k, 1250 * k, 1250 * (k + 1), 8000 - t - e, 8000 - t))
    # Ouest : extérieur x < -100, 6 x 1250 + 500
    for k in range(6):
        panels.append(osb(300 + k, -t - e, -t, 1250 * k, 1250 * (k + 1)))
    panels.append(osb(306, -t - e, -t, 7500, 8000))
    # Est : panneaux empilés 2 x 1350 de haut, extérieur x > 10100
    for k in range(4):
        panels.append(osb(400 + 2 * k, 10000 + t, 10000 + t + e, 2000 * k, 2000 * (k + 1), 0, 1350))
        panels.append(osb(401 + 2 * k, 10000 + t, 10000 + t + e, 2000 * k, 2000 * (k + 1), 1350, H))
    return walls, panels, openings


@pytest.fixture
def export():
    return cw.build_mob_export(*model())


@pytest.fixture
def mob():
    batiment = Batiment(h_bat=3, d_bat=10.0, b_bat=20.0, alpha_toit=10)
    panel = Barre(b=12, h=1250, section="Rectangulaire", classe="OSB/3 11-18 mm", cs=1, effet_systeme=False)
    montant = Barre(b=45, h=145, section="Rectangulaire", classe="C24", cs=1, effet_systeme=False)
    kw = {
        "d": 2.5,
        "dh": 5,
        "l": 55,
        "qualite": "6.8",
        "n": 1,
        "type_organe": "Pointe circulaire lisse",
        "percage": False,
        "nfile": 1,
        "nCis": 1,
    }
    conn = Pointe(beam_1=panel, beam_2=montant, **kw)
    conn2 = Pointe(beam_1=panel, beam_2=montant, **kw)
    return EC5.MOB._from_parent_class(batiment, connecteurs_int=conn, connecteurs_ext=conn2)


def _by_name(export):
    return {s["name"]: s for s in export["sys_walls"]}


class TestCadworkExport:
    def test_systems_and_geometry(self, export):
        sys_walls = _by_name(export)
        assert set(sys_walls) == {"Sud", "Nord", "Ouest", "Est"}
        sud = sys_walls["Sud"]
        assert sud["etage"] == "RDC" and sud["h_sys_MOB"] == 2.7 and sud["l_sys_MOB"] == 10.2
        # Vue depuis l'extérieur (sud), les abscisses vont vers l'est
        assert sud["alpha_sw"] == 0
        assert sys_walls["Nord"]["alpha_sw"] == 180
        assert sys_walls["Est"]["alpha_sw"] == 90
        assert sys_walls["Ouest"]["alpha_sw"] == 270
        assert sud["Yg"] == 0 and sys_walls["Est"]["Xg"] == 10

    def test_openings_split_walls(self, export):
        sud = _by_name(export)["Sud"]
        assert len(sud["walls"]) == 2
        left, right = sud["walls"]
        assert (left["position"], left["l_MOB"]) == (0.1, 3.75)
        assert (right["position"], right["l_MOB"]) == (5.3, 4.8)
        assert [(p["number"], p["b_panel"]) for p in left["panels"]] == [(3, 1250)]
        assert sorted((p["number"], p["b_panel"]) for p in right["panels"]) == [(1, 1050), (3, 1250)]
        assert {p["position_panel"] for w in sud["walls"] for p in w["panels"]} == {"Extérieur"}
        reasons = {i["cadwork_id"]: i["reason"] for i in export["ignored"]}
        assert reasons[120] == "chevauche une ouverture"
        assert 130 not in reasons  # filtré par mot-clé, jamais candidat

    def test_sides_and_stack(self, export):
        sys_walls = _by_name(export)
        nord = sys_walls["Nord"]["walls"][0]["panels"]
        assert nord[0]["position_panel"] == "Intérieur" and nord[0]["number"] == 8
        est = sys_walls["Est"]["walls"][0]["panels"]
        assert len(est) == 8
        tops = [p for p in est if p["on_top_of"]]
        assert len(tops) == 4
        keys = {p["key"]: p for p in est}
        assert all(keys[p["on_top_of"]]["on_top_of"] is None for p in tops)
        assert all(p["h_panel"] == 1350 for p in est)

    def test_rotation_invariant(self, export):
        rotated = cw.build_mob_export(*model(angle=30))
        for a, b in zip(export["sys_walls"], rotated["sys_walls"]):
            assert a["l_sys_MOB"] == b["l_sys_MOB"]
            assert [w["l_MOB"] for w in a["walls"]] == [w["l_MOB"] for w in b["walls"]]
            assert math.isclose((a["alpha_sw"] + 30) % 360, b["alpha_sw"], abs_tol=0.01)
            xg, yg, _ = _rot((a["Xg"], a["Yg"], 0), 30)
            assert math.isclose(xg, b["Xg"], abs_tol=2e-3) and math.isclose(yg, b["Yg"], abs_tol=2e-3)

    def test_no_keyword_filter_rejects_thick_panel(self):
        data = cw.build_mob_export(*model(), options={"panel_keywords": []})
        reasons = {i["cadwork_id"]: i["reason"] for i in data["ignored"]}
        assert reasons[130].startswith("non parallèle")


class TestLoadCadworkData:
    def test_regenerates_mob(self, export, mob, tmp_path):
        path = tmp_path / "cadwork.json"
        path.write_text(json.dumps(export), encoding="utf-8")
        mapping = mob.load_cadwork_data(str(path), couturage_ext_rive=100, couturage_ext_inter=200)
        assert mapping == {"S1": "Sud", "S2": "Nord", "S3": "Ouest", "S4": "Est"}
        walls = mob.data_walls
        assert list(walls) == ["S1_W1", "S1_W2", "S2_W1", "S3_W1", "S4_W1"]
        assert walls["S1_W1"]["Couturage ext. rive"] == 100 * si.mm
        assert "Couturage int. rive" not in walls["S1_W1"]
        assert "Couturage ext. rive" not in walls["S2_W1"] and walls["S2_W1"]["Couturage int. rive"] == 150 * si.mm
        panels = mob.data_panels
        assert panels["S1_W1_P1"]["Nombre"] == 3 and panels["S1_W1_P1"]["Largeur"] == 1250 * si.mm
        assert len(mob._list_on_top_of) == 4
        assert panels["S4_W1_P2"]["On top of"] == "S4_W1_P1"

        Ks = mob.get_Kser_wall()
        assert all(Ks["Kser système de mur"].apply(lambda k: k.value) > 0)
        res = mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=15, XFw=5, YFw=4, etage="RDC")
        F = sum(
            res.loc[s, "Fv,Ed,i ELU"].value * math.cos(math.radians(mob.data["RDC"][s]["alpha_sw"])) for s in res.index
        )
        assert np.isclose(abs(F), 15e3)

        save = tmp_path / "mob.json"
        mob.save_walls_data(str(save))
        reloaded = EC5.MOB._from_parent_class(mob)
        reloaded.load_walls_data(str(save))
        assert list(reloaded.data_panels) == list(panels)

    def test_matches_manual_definition(self, export, mob, tmp_path):
        path = tmp_path / "cadwork.json"
        path.write_text(json.dumps(export), encoding="utf-8")
        mob.load_cadwork_data(str(path))
        k_import = mob.get_Kser_wall().loc["S1", "Kser système de mur"].value

        manual = EC5.MOB._from_parent_class(
            mob, connecteurs_int=mob._connecteurs["Intérieur"], connecteurs_ext=mob._connecteurs["Extérieur"]
        )
        s = manual.add_sys_wall(h_etage=2.7, h_sys_MOB=2.7, l_sys_MOB=10.2, etage="RDC")
        w1 = manual.add_wall_to_sys(s, position=0.1, l_MOB=3.75, couturage_int_rive=0, couturage_int_inter=0)
        manual.add_panel_to_wall(w1, h_panel=2700, b_panel=1250, number=3, position_panel="Extérieur")
        w2 = manual.add_wall_to_sys(s, position=5.3, l_MOB=4.8, couturage_int_rive=0, couturage_int_inter=0)
        manual.add_panel_to_wall(w2, h_panel=2700, b_panel=1250, number=3, position_panel="Extérieur")
        manual.add_panel_to_wall(w2, h_panel=2700, b_panel=1050, number=1, position_panel="Extérieur")
        k_manual = manual.get_Kser_wall().loc["S1", "Kser système de mur"].value
        assert np.isclose(k_import, k_manual)

    def test_stacked_panels_kser_in_series(self, export, mob, tmp_path):
        path = tmp_path / "cadwork.json"
        path.write_text(json.dumps(export), encoding="utf-8")
        mob.load_cadwork_data(str(path))
        mob.get_Kser_wall()
        k_wall = mob.data_Kser_walls.loc["S4_W1", "Kser mur"].value
        kp = mob.data_Kser_panels["Kser,p unit."].apply(lambda k: k.value)
        panels = {n: p for n, p in mob.data_panels.items() if n.startswith("S4_W1_")}
        tops = {n: p["On top of"] for n, p in panels.items() if p.get("On top of")}
        assert tops
        expected = sum(panels[b]["Nombre"] * kp[t] * kp[b] / (kp[t] + kp[b]) for t, b in tops.items())
        expected += sum(p["Nombre"] * kp[n] for n, p in panels.items() if n not in tops and n not in tops.values())
        assert np.isclose(k_wall, expected)

    def test_missing_connector(self, export, tmp_path):
        batiment = Batiment(h_bat=3, d_bat=10.0, b_bat=20.0, alpha_toit=10)
        panel = Barre(b=12, h=1250, section="Rectangulaire", classe="OSB/3 11-18 mm", cs=1, effet_systeme=False)
        montant = Barre(b=45, h=145, section="Rectangulaire", classe="C24", cs=1, effet_systeme=False)
        conn = Pointe(
            d=2.5,
            dh=5,
            l=55,
            qualite="6.8",
            n=1,
            type_organe="Pointe circulaire lisse",
            percage=False,
            beam_1=panel,
            beam_2=montant,
            nfile=1,
            nCis=1,
        )
        mob = EC5.MOB._from_parent_class(batiment, connecteurs_ext=conn)
        path = tmp_path / "cadwork.json"
        path.write_text(json.dumps(export), encoding="utf-8")
        with pytest.raises(ValueError, match="Intérieur"):
            mob.load_cadwork_data(str(path))
        assert mob.data == {}

    def test_wrong_format(self, mob, tmp_path):
        path = tmp_path / "other.json"
        path.write_text(json.dumps({"RDC": {}}), encoding="utf-8")
        with pytest.raises(ValueError, match="Cadwork"):
            mob.load_cadwork_data(str(path))


class TestCadworkApi:
    """Simule les contrôleurs Cadwork pour vérifier la lecture du modèle et l'écriture du JSON."""

    @pytest.fixture
    def fake_cadwork(self, monkeypatch):
        walls, panels, openings = model()
        elements = {e["id"]: ("wall", e) for e in walls}
        elements.update({e["id"]: ("panel", e) for e in panels})
        elements.update({e["id"]: ("opening", e) for e in openings})
        messages = []
        mods = {
            "cadwork": {},
            "element_controller": {
                "get_active_identifiable_element_ids": list,
                "get_visible_identifiable_element_ids": lambda: list(elements),
            },
            "attribute_controller": {
                "get_name": lambda i: elements[i][1]["name"],
                "is_framed_wall": lambda i: elements[i][0] == "wall",
                "is_opening": lambda i: elements[i][0] == "opening",
                "is_panel": lambda i: elements[i][0] == "panel",
                "get_element_material_name": lambda i: elements[i][1].get("material", ""),
                "get_group": lambda i: "",
                "get_subgroup": lambda i: "",
                "get_comment": lambda i: "",
            },
            "geometry_controller": {
                "get_element_vertices": lambda i: elements[i][1]["vertices"],
                "get_xl": lambda i: elements[i][1]["xl"],
            },
            "bim_controller": {
                "get_storey": lambda i: "RDC",
                "get_building": lambda i: "B",
                "get_all_buildings": lambda: ["B"],
                "get_all_storeys": lambda b: ["RDC", "R+1"],
                "get_storey_height": lambda b, s: {"RDC": 0.0, "R+1": 2950.0}[s],
            },
            "utility_controller": {
                "get_3d_file_path": lambda: "C:/projet.3d",
                "get_new_user_file_from_dialog": lambda f: "",
                "print_message": lambda msg, r, c: messages.append(msg),
            },
        }
        for name, funcs in mods.items():
            monkeypatch.setitem(sys.modules, name, types.SimpleNamespace(**funcs))
        return messages

    def test_export_from_cadwork(self, fake_cadwork, tmp_path):
        path = cw.export_from_cadwork(str(tmp_path / "export"))
        assert path.endswith("export.json")
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        assert data["format"] == cw.FORMAT and data["source"] == "C:/projet.3d"
        assert len(data["sys_walls"]) == 4 and data["sys_walls"][0]["h_etage"] == 2.95
        assert "4 parois, 5 murs" in fake_cadwork[0]

    def test_elements_hors_mob_ignores(self, fake_cadwork, tmp_path, monkeypatch):
        """Les éléments hors MOB ne sont pas interrogés (leur géométrie peut être illisible)."""
        reference = json.loads(Path(cw.export_from_cadwork(str(tmp_path / "ref.json"))).read_text(encoding="utf-8"))
        gc, ac, ec = (sys.modules[n] for n in ("geometry_controller", "attribute_controller", "element_controller"))
        ids = ec.get_visible_identifiable_element_ids()

        def vertices(i, f=gc.get_element_vertices):
            if i == 999:
                raise TypeError("Unable to convert function return value to a Python type!")
            return f(i)

        monkeypatch.setattr(ec, "get_visible_identifiable_element_ids", lambda: [*ids, 999])
        originals = {n: getattr(ac, n) for n in ("is_framed_wall", "is_opening", "is_panel")}
        for name, func in originals.items():
            monkeypatch.setattr(ac, name, lambda i, f=func: i != 999 and f(i))
        monkeypatch.setattr(gc, "get_element_vertices", vertices)
        data = json.loads(Path(cw.export_from_cadwork(str(tmp_path / "out.json"))).read_text(encoding="utf-8"))
        assert data["sys_walls"] == reference["sys_walls"]

    def test_dialog_cancelled(self, fake_cadwork):
        assert cw.export_from_cadwork() == ""


class TestCadworkEdgeCases:
    H = 2700

    def _facade_and_refend(self):
        walls = [
            {
                "id": 1,
                "name": "Facade",
                "storey": "RDC",
                "xl": (1, 0, 0),
                "vertices": box(0, 10000, -100, 100, 0, self.H),
            },
            {
                "id": 2,
                "name": "Refend",
                "storey": "RDC",
                "xl": (0, 1, 0),
                "vertices": box(4900, 5100, 100, 4000, 0, self.H),
            },
        ]
        return walls

    def test_t_junction_panel_kept_on_facade(self):
        walls = self._facade_and_refend()
        # Panneau intérieur de façade centré au droit du refend : centre aussi contenu dans le refend
        panels = [{"id": 10, "name": "OSB", "vertices": box(4400, 5600, 100, 112, 0, self.H)}]
        export = cw.build_mob_export(walls, panels)
        facade = _by_name(export)["Facade"]
        assert not export["ignored"]
        assert facade["walls"][0]["panels"][0]["b_panel"] == 1200

    def test_overlapping_layers_not_stacked(self):
        walls = self._facade_and_refend()[:1]
        panels = [
            {"id": 10, "name": "OSB", "vertices": box(0, 1250, -112, -100, 0, self.H)},
            {"id": 11, "name": "OSB", "vertices": box(0, 1250, -124, -112, 0, self.H)},
        ]
        export = cw.build_mob_export(walls, panels)
        out = export["sys_walls"][0]["walls"][0]["panels"]
        assert all(p["on_top_of"] is None for p in out)

    def test_storey_heights_per_building(self):
        levels = {"A": {"RDC": 0.0, "R+1": 3000.0}, "B": {"RDC": 500.0, "R+1": 3300.0}}
        bc = types.SimpleNamespace(
            get_all_buildings=lambda: list(levels),
            get_all_storeys=lambda b: list(levels[b]),
            get_storey_height=lambda b, s: levels[b][s],
        )
        heights = cw._storey_heights(bc, {"A / RDC", "B / RDC"}, multi=True)
        assert heights == {"A / RDC": 3.0, "B / RDC": 2.8}

    def test_layers_of_different_thickness_not_grouped(self):
        walls = self._facade_and_refend()[:1]
        panels = [
            {"id": 10, "name": "OSB", "vertices": box(0, 1250, -112, -100, 0, self.H)},
            {"id": 11, "name": "OSB", "vertices": box(0, 1250, -127, -112, 0, self.H)},
        ]
        out = cw.build_mob_export(walls, panels)["sys_walls"][0]["walls"][0]["panels"]
        assert sorted(p["thickness"] for p in out) == [12, 15]

    def test_unassigned_storey_keeps_elevation_fallback(self):
        assert cw._storey_label("B", "", True) == ""
        assert cw._storey_label("B", "RDC", True) == "B / RDC"

    @pytest.mark.parametrize(("default", "expected"), [(None, "Niveau 0.00 m"), ("RDC", "RDC")])
    def test_unassigned_storeys_split_by_building(self, default, expected):
        walls = self._facade_and_refend()
        for w, building in zip(walls, "AB", strict=True):
            w.update(storey="", building=building)
        panels = [
            {"id": 10, "name": "OSB", "vertices": box(0, 1250, -112, -100, 0, self.H)},
            {"id": 11, "name": "OSB", "vertices": box(4888, 4900, 200, 1450, 0, self.H)},
        ]
        export = cw.build_mob_export(walls, panels, options={"default_storey": default})
        assert {s["etage"] for s in export["sys_walls"]} == {f"A / {expected}", f"B / {expected}"}
