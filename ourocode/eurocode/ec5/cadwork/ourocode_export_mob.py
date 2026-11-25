# Encoding in UTF-8 by Anthony PARISOT
"""Export des murs ossature bois (MOB) de Cadwork 3D vers ourocode.

Ce module est **autonome** (bibliothèque standard uniquement) : il s'exécute dans
l'interpréteur Python de Cadwork 3D, sans ourocode installé. Il produit un fichier JSON
relu côté ourocode par ``MOB.load_cadwork_data``, qui régénère les appels
``add_sys_wall`` / ``add_wall_to_sys`` / ``add_panel_to_wall``.

Installation dans Cadwork (plugin API) :
    copier ce fichier dans ``<userprofil>/3d/API.x64/ourocode_export_mob/ourocode_export_mob.py``
    puis lancer le plugin « ourocode_export_mob » depuis la barre des plugins.

Règles de conversion (voir ``build_mob_export``) :
    - chaque « paroi ossature » (``attribute_controller.is_framed_wall``) devient un système de mur ;
    - les ouvertures (``is_opening``) situées dans la paroi la découpent en murs (trumeaux) ;
    - les panneaux (``is_panel``) situés dans la paroi, parallèles à son plan, filtrés par mots-clés,
      deviennent des panneaux de contreventement ; ceux qui chevauchent une ouverture sont ignorés ;
    - la face extérieure est celle opposée au centre de gravité en plan des parois de l'étage ;
    - deux panneaux superposés dans un même trumeau sont liés par ``on_top_of``.

Le contenu de la sélection active est exporté ; à défaut, tous les éléments visibles.
"""

from __future__ import annotations

import itertools
import json
import math

FORMAT = "ourocode.MOB.cadwork"
VERSION = 1

# Options par défaut, modifiables en tête de script ou via build_mob_export(options=...)
DEFAULT_OPTIONS = {
    # Mots-clés (insensibles à la casse) recherchés dans le nom, le matériau, le groupe, le sous-groupe
    # et le commentaire des panneaux. Liste vide = tous les panneaux de la paroi sont retenus.
    "panel_keywords": ["OSB", "contrevent", "CTB", "contreplaqu"],
    # Épaisseur maximale (mm) d'un panneau de contreventement (élimine isolants, panneaux à plat...).
    "max_panel_thickness": 40.0,
    # Tolérance géométrique (mm) d'appartenance d'un élément à une paroi.
    "tolerance": 20.0,
    # Longueur minimale (mm) d'un trumeau conservé.
    "min_wall_length": 100.0,
    # Nom de l'étage si l'élément n'est rattaché à aucun étage BIM.
    "default_storey": None,
}


# ----------------------------------------------------------------------
# Géométrie (pure Python, testable hors Cadwork)
# ----------------------------------------------------------------------


def _dot2(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _extent(vertices, origin, u, n):
    """Étendues (min, max) d'un nuage de points selon u, n (plan) et z."""
    s = [_dot2((v[0] - origin[0], v[1] - origin[1]), u) for v in vertices]
    t = [_dot2((v[0] - origin[0], v[1] - origin[1]), n) for v in vertices]
    z = [v[2] for v in vertices]
    return (min(s), max(s)), (min(t), max(t)), (min(z), max(z))


def _horizontal_axis(xl, vertices):
    """Direction horizontale unitaire de la paroi (XL projeté, sinon plus grande étendue en plan)."""
    if xl is not None:
        norm = math.hypot(xl[0], xl[1])
        if norm > 1e-6:
            return (xl[0] / norm, xl[1] / norm)
    best, best_len = (1.0, 0.0), 0.0
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            d = math.hypot(b[0] - a[0], b[1] - a[1])
            if d > best_len:
                best, best_len = ((b[0] - a[0]) / d, (b[1] - a[1]) / d), d
    return best


def _text_of(element):
    keys = ("name", "material", "group", "subgroup", "comment")
    return " ".join(str(element.get(k) or "") for k in keys).lower()


def _subtract_intervals(base, holes):
    """Retire les intervalles holes de l'intervalle base. Retourne la liste des intervalles restants."""
    segments = [base]
    for h0, h1 in sorted(holes):
        new = []
        for a, b in segments:
            if h1 <= a or h0 >= b:
                new.append((a, b))
                continue
            if h0 > a:
                new.append((a, h0))
            if h1 < b:
                new.append((h1, b))
        segments = new
    return segments


def _wall_frame(wall):
    """Repère de la paroi : origine, axe u, normale n, étendues et épaisseur (mm)."""
    verts = wall["vertices"]
    a = _horizontal_axis(wall.get("xl"), verts)
    n = (-a[1], a[0])
    p0 = (verts[0][0], verts[0][1])
    (s0, s1), (t0, t1), (z0, z1) = _extent(verts, p0, a, n)
    return {"p0": p0, "u": a, "n": n, "s": (s0, s1), "t": (t0, t1), "z": (z0, z1), "vertices": verts}


def _round(x, nd=3):
    return round(float(x), nd)


def build_mob_export(
    walls: list, panels: list, openings: list = (), options: dict | None = None, storeys: dict | None = None
) -> dict:
    """Construit le dictionnaire d'export MOB à partir de la géométrie brute Cadwork.

    Args:
        walls (list): Parois ossature : dicts ``{"id", "name", "storey", "vertices", "xl"}``,
            ``vertices`` = liste de points (x, y, z) globaux en mm, ``xl`` = axe de longueur.
        panels (list): Panneaux candidats : dicts ``{"id", "name", "material", "group", "subgroup",
            "comment", "storey", "vertices"}``.
        openings (list): Ouvertures : dicts ``{"id", "name", "vertices"}``.
        options (dict, optional): Surcharge de DEFAULT_OPTIONS.
        storeys (dict, optional): ``{nom_étage: h_etage en m}`` (hauteur d'étage, plancher à plancher).

    Returns:
        dict: Données JSON-sérialisables (longueurs de plan en m, dimensions de panneaux en mm).
    """
    opt = dict(DEFAULT_OPTIONS)
    opt.update(options or {})
    tol = opt["tolerance"]
    keywords = [k.lower() for k in opt["panel_keywords"]]
    storeys = storeys or {}
    warnings, ignored = [], []

    frames = {}
    for wall in walls:
        frame = _wall_frame(wall)
        frame["storey"] = wall.get("storey") or opt["default_storey"] or f"Niveau {frame['z'][0] / 1000:.2f} m"
        frames[wall["id"]] = frame

    # Face extérieure : opposée au centre de gravité en plan des parois du même étage
    centroids = {}
    for frame in frames.values():
        mid = _mid_plan(frame)
        centroids.setdefault(frame["storey"], []).append(mid)
    for frame in frames.values():
        pts = centroids[frame["storey"]]
        g = (sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts))
        mid = _mid_plan(frame)
        sign = 1 if _dot2((mid[0] - g[0], mid[1] - g[1]), frame["n"]) >= 0 else -1
        n_ext = (sign * frame["n"][0], sign * frame["n"][1])
        # Vue depuis l'extérieur, l'axe des abscisses va de gauche à droite : u = (-n_ext) x z
        u = (-n_ext[1], n_ext[0])
        (s0, s1), (t0, t1), (z0, z1) = _extent(frame["vertices"], frame["p0"], u, n_ext)
        frame.update({"u": u, "n": n_ext, "s": (s0, s1), "t": (t0, t1), "z": (z0, z1)})

    def locate(element):
        """Parois contenant le centre de l'élément (la plus proche du plan médian en premier)."""
        candidates = []
        for wid, f in frames.items():
            (s0, s1), (t0, t1), (z0, z1) = _extent(element["vertices"], f["p0"], f["u"], f["n"])
            cs, ct, cz = (s0 + s1) / 2, (t0 + t1) / 2, (z0 + z1) / 2
            if not (f["s"][0] - tol <= cs <= f["s"][1] + tol and f["z"][0] - tol <= cz <= f["z"][1] + tol):
                continue
            if not (f["t"][0] - tol <= ct <= f["t"][1] + tol):
                continue
            dist = abs(ct - (f["t"][0] + f["t"][1]) / 2)
            candidates.append((dist, wid, (s0, s1), (t0, t1), (z0, z1)))
        return sorted(candidates, key=lambda c: c[0])

    holes = {wid: [] for wid in frames}
    for op in openings:
        locs = locate(op)
        if locs:
            loc = locs[0]
            f = frames[loc[1]]
            holes[loc[1]].append((loc[2][0] - f["s"][0], loc[2][1] - f["s"][0]))

    wall_panels = {wid: [] for wid in frames}
    for p in panels:
        if keywords and not any(k in _text_of(p) for k in keywords):
            continue
        locs = locate(p)
        if not locs:
            ignored.append({"cadwork_id": p["id"], "name": p.get("name"), "reason": "hors paroi ossature"})
            continue
        # Aux jonctions en T, le centre d'un panneau peut aussi tomber dans la paroi perpendiculaire
        for _, wid, (s0, s1), (t0, t1), (z0, z1) in locs:
            thickness, b, h = t1 - t0, s1 - s0, z1 - z0
            if thickness <= opt["max_panel_thickness"] and thickness < min(b, h):
                break
        else:
            ignored.append(
                {"cadwork_id": p["id"], "name": p.get("name"), "reason": "non parallèle à la paroi ou trop épais"}
            )
            continue
        f = frames[wid]
        side = "Extérieur" if (t0 + t1) / 2 >= (f["t"][0] + f["t"][1]) / 2 else "Intérieur"
        wall_panels[wid].append(
            {
                "id": p["id"],
                "name": p.get("name"),
                "material": p.get("material"),
                "s": (s0 - f["s"][0], s1 - f["s"][0]),
                "z": (z0, z1),
                "b": b,
                "h": h,
                "thickness": thickness,
                "side": side,
            }
        )

    sys_walls = []
    wall_by_id = {w["id"]: w for w in walls}
    for wid, f in frames.items():
        L = f["s"][1] - f["s"][0]
        H = f["z"][1] - f["z"][0]
        segments = [(a, b) for a, b in _subtract_intervals((0.0, L), holes[wid]) if b - a >= opt["min_wall_length"]]
        sub_walls, weights = [], []
        remaining = list(wall_panels[wid])
        for a, b in segments:
            inside = [p for p in remaining if p["s"][0] >= a - tol and p["s"][1] <= b + tol]
            remaining = [p for p in remaining if p not in inside]
            if not inside:
                continue
            pos0 = min(p["s"][0] for p in inside)
            pos1 = max(p["s"][1] for p in inside)
            panels_out = _group_panels(inside, warnings, wall_by_id[wid].get("name"), tol)
            for p in inside:
                weights.append((p["b"], (p["s"][0] + p["s"][1]) / 2))
            sub_walls.append(
                {"position": _round(pos0 / 1000), "l_MOB": _round((pos1 - pos0) / 1000), "panels": panels_out}
            )
        for p in remaining:
            ignored.append({"cadwork_id": p["id"], "name": p["name"], "reason": "chevauche une ouverture"})
        if not sub_walls:
            warnings.append(f"Paroi {wall_by_id[wid].get('name') or wid} : aucun panneau de contreventement retenu.")
            continue
        xc = sum(w * c for w, c in weights) / sum(w for w, _ in weights)
        origin = (f["p0"][0] + f["s"][0] * f["u"][0], f["p0"][1] + f["s"][0] * f["u"][1])
        t_mid = (f["t"][0] + f["t"][1]) / 2
        origin = (origin[0] + t_mid * f["n"][0], origin[1] + t_mid * f["n"][1])
        Xg, Yg = origin[0] + xc * f["u"][0], origin[1] + xc * f["u"][1]
        sys_walls.append(
            {
                "cadwork_id": wid,
                "name": wall_by_id[wid].get("name"),
                "etage": f["storey"],
                "h_etage": _round(storeys.get(f["storey"]) or H / 1000),
                "h_sys_MOB": _round(H / 1000),
                "l_sys_MOB": _round(L / 1000),
                "Xg": _round(Xg / 1000),
                "Yg": _round(Yg / 1000),
                "alpha_sw": _round(math.degrees(math.atan2(f["u"][1], f["u"][0])) % 360, 2),
                "walls": sub_walls,
            }
        )
    return {
        "format": FORMAT,
        "version": VERSION,
        "units": {"plan": "m", "panneaux": "mm"},
        "options": opt,
        "sys_walls": sys_walls,
        "ignored": ignored,
        "warnings": warnings,
    }


def _mid_plan(frame):
    s = (frame["s"][0] + frame["s"][1]) / 2
    t = (frame["t"][0] + frame["t"][1]) / 2
    return (
        frame["p0"][0] + s * frame["u"][0] + t * frame["n"][0],
        frame["p0"][1] + s * frame["u"][1] + t * frame["n"][1],
    )


def _group_panels(panels, warnings, wall_label, tol):
    """Regroupe les panneaux d'un trumeau : piles verticales (on_top_of) et panneaux identiques (number)."""
    out = []
    for side in ("Extérieur", "Intérieur"):
        columns = []
        for p in sorted((p for p in panels if p["side"] == side), key=lambda p: (p["s"][0], p["z"][0])):
            for col in columns:
                ref = col[0]
                overlap = min(ref["s"][1], p["s"][1]) - max(ref["s"][0], p["s"][0])
                stacked = all(p["z"][0] >= q["z"][1] - tol or p["z"][1] <= q["z"][0] + tol for q in col)
                if overlap > 0.5 * min(ref["b"], p["b"]) and stacked:
                    col.append(p)
                    break
            else:
                columns.append([p])
        singles = {}
        for col in columns:
            col.sort(key=lambda p: p["z"][0])
            if len(col) == 1:
                p = col[0]
                key = (round(p["h"]), round(p["b"]))
                singles.setdefault(key, []).append(p)
                continue
            if len(col) > 2:
                warnings.append(
                    f"Paroi {wall_label} : pile de {len(col)} panneaux, seuls les deux premiers sont liés (on_top_of)."
                )
            below = None
            for i, p in enumerate(col):
                entry = _panel_entry([p], side)
                if i == 1:
                    entry["on_top_of"] = below
                out.append(entry)
                below = entry["key"]
        for group in singles.values():
            out.append(_panel_entry(group, side))
    return out


def _panel_entry(group, side):
    p = group[0]
    return {
        "key": f"cw{p['id']}",
        "cadwork_ids": [g["id"] for g in group],
        "position_panel": side,
        "h_panel": _round(p["h"], 1),
        "b_panel": _round(p["b"], 1),
        "number": len(group),
        "thickness": _round(p["thickness"], 1),
        "material": p.get("material"),
        "on_top_of": None,
    }


# ----------------------------------------------------------------------
# Lecture du modèle Cadwork (exécuté dans Cadwork 3D)
# ----------------------------------------------------------------------


def _pt(p):
    return (float(p[0]), float(p[1]), float(p[2]))


def _storey_label(building, storey, multi):
    """Nom d'étage unique : préfixé par le bâtiment BIM si le modèle en contient plusieurs."""
    return f"{building} / {storey}" if multi and building else storey


def _storey_heights(bc, storey_names, multi=False):
    """Hauteurs d'étage (m) déduites des altitudes BIM, par bâtiment (licence IFC requise, sinon dict vide)."""
    try:
        heights = {}
        for building in bc.get_all_buildings():
            levels = sorted((bc.get_storey_height(building, s), s) for s in bc.get_all_storeys(building))
            for (z, name), (z_up, _) in itertools.pairwise(levels):
                label = _storey_label(building, name, multi)
                if label in storey_names:
                    heights[label] = (z_up - z) / 1000
        return heights
    except Exception:  # noqa: BLE001 - API Cadwork : licence IFC absente, etc.
        return {}


def collect_from_cadwork(options: dict | None = None) -> dict:
    """Lit la sélection active (ou les éléments visibles) dans Cadwork et construit l'export MOB."""
    import attribute_controller as ac
    import element_controller as ec
    import geometry_controller as gc

    try:
        import bim_controller as bc
    except ImportError:
        bc = None

    ids = ec.get_active_identifiable_element_ids() or ec.get_visible_identifiable_element_ids()
    walls, panels, openings = [], [], []
    for eid in ids:
        base = {"id": int(eid), "name": ac.get_name(eid), "vertices": [_pt(v) for v in gc.get_element_vertices(eid)]}
        if not base["vertices"]:
            continue
        base["storey"] = bc.get_storey(eid) if bc else ""
        base["building"] = bc.get_building(eid) if bc else ""
        if ac.is_framed_wall(eid):
            base["xl"] = _pt(gc.get_xl(eid))
            walls.append(base)
        elif ac.is_opening(eid):
            openings.append(base)
        elif ac.is_panel(eid):
            base.update(
                {
                    "material": ac.get_element_material_name(eid),
                    "group": ac.get_group(eid),
                    "subgroup": ac.get_subgroup(eid),
                    "comment": ac.get_comment(eid),
                }
            )
            panels.append(base)
    multi = len({w["building"] for w in walls}) > 1
    for w in walls:
        w["storey"] = _storey_label(w["building"], w["storey"], multi)
    storeys = _storey_heights(bc, {w["storey"] for w in walls}, multi) if bc else {}
    return build_mob_export(walls, panels, openings, options=options, storeys=storeys)


def export_from_cadwork(path: str | None = None, options: dict | None = None) -> str:
    """Exporte les MOB du modèle Cadwork ouvert dans un fichier JSON lisible par MOB.load_cadwork_data.

    Args:
        path (str, optional): Fichier de destination. Si None, une boîte de dialogue Cadwork s'ouvre.
        options (dict, optional): Surcharge de DEFAULT_OPTIONS.

    Returns:
        str: Chemin du fichier écrit (chaîne vide si annulé).
    """
    import utility_controller as uc

    data = collect_from_cadwork(options)
    data["source"] = uc.get_3d_file_path()
    if not path:
        path = uc.get_new_user_file_from_dialog("JSON (*.json)")
    if not path:
        return ""
    if not path.lower().endswith(".json"):
        path += ".json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4, ensure_ascii=False)
    n_walls = sum(len(s["walls"]) for s in data["sys_walls"])
    n_panels = sum(len(w["panels"]) for s in data["sys_walls"] for w in s["walls"])
    msg = f"ourocode : {len(data['sys_walls'])} parois, {n_walls} murs, {n_panels} panneaux exportés -> {path}"
    if data["warnings"]:
        msg += f" ({len(data['warnings'])} avertissement(s))"
    uc.print_message(msg, 1, 1)
    return path


if __name__ == "__main__":
    export_from_cadwork()
