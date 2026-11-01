#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
from math import sqrt, radians, cos, sin, floor
from copy import deepcopy
import pandas as pd
import numpy as np
import forallpeople as si

si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.core.batiment import Batiment
from ourocode.eurocode.ec5.assemblage import Assemblage


class MOB(Batiment):
    CRITERE_ELS = 200
    RATIO_KU_KSER = 2 / 3  # Ku = 2/3 Kser, EN 1995-1-1 §2.2.2(1)b
    POSITION_PANEL = ("Extérieur", "Intérieur")
    COEFF_DOUBLE_PANEL = {
        "Non": 1,
        "Panneau et organe identiques": 1,
        "Panneau différent & glissement identique": 0.75,
        "Panneau différent & glissement différent": 0.5,
    }

    def __init__(self, connecteurs_int: object = None, connecteurs_ext: object = None, **kwargs):
        """Initialise un système de murs à ossature bois (MOB) selon l'EN 1995-1-1 §9.2.4 méthode A.

        Un système MOB est composé de plusieurs murs de contreventement (interrompus par portes/fenêtres),
        eux-mêmes constitués de plusieurs panneaux. Le workflow recommandé est :

        1. Instancier MOB avec les connecteurs.
        2. Ajouter les systèmes de mur via add_sys_wall.
        3. Ajouter les murs internes via add_wall_to_sys.
        4. Ajouter les panneaux via add_panel_to_wall.
        5. (Optionnel) Ajouter des murs béton (refends) via add_concrete_wall.
        6. Calculer les efforts via calculate_loads_and_deformations ou en fournissant Fv_Ed directement dans add_sys_wall.
        7. Récupérer les taux de vérification via taux_walls et taux_ancrage, la torsion via synthese_torsion.
        8. Visualiser la vue en plan via show_plan.

        Cette classe ne calcule pas les montants au flambement ni la compression perpendiculaire
        sur la lisse basse : utiliser Barre et les classes de vérification EC5_Element_droit.

        Args:
            connecteurs_int (Agrafe | Pointe, optional): Connecteur intérieur (panneau face intérieure)
                issu du module EC5_Assemblage. Instancier avec n_file=1 et n=1. Defaults to None.
            connecteurs_ext (Agrafe | Pointe, optional): Connecteur extérieur (panneau face extérieure)
                issu du module EC5_Assemblage. Instancier avec n_file=1 et n=1. Defaults to None.
            **kwargs: Arguments transmis à la classe parent Batiment.
        """
        super().__init__(**kwargs)
        self._connecteurs = {"Intérieur": connecteurs_int, "Extérieur": connecteurs_ext}
        self._data_MOB = {}
        self._data_connecteurs = {}
        self._data_concrete = {}
        self._list_on_top_of = []
        self._setup_MOB()

    def _setup_MOB(self):
        """Initialise les matériaux, sections et valeurs de calcul pour les murs à ossature bois."""
        for position, connecteur in self._connecteurs.items():
            if connecteur is None:
                continue
            effet_corde = False
            if connecteur.type_organe in [
                "Tirefond",
                "Pointe carrée lisse",
                "Pointe circulaire lisse",
                "Autres pointes",
            ]:
                connecteur.Fax_Rk()
                effet_corde = True
            Fvrk = connecteur.FvRk(effet_corde=effet_corde)
            Kser = connecteur.Kser()
            self._data_connecteurs[position] = {"Ff,Rk": Fvrk[1], "Kser": Kser[1]}

    def _get_all_sys_walls_names(self):
        """Retourne la liste de tous les noms des systèmes de mur."""
        data = []
        for etage, sys in self.data.items():
            for sys_name in sys.keys():
                data.append(sys_name)
        return data

    # ------------------------------------------------------------------
    # Outils internes : étages, unités, solveur diaphragme rigide
    # ------------------------------------------------------------------

    @staticmethod
    def _val(x) -> float:
        """Retourne la valeur SI (float) d'une grandeur Physical, ou la valeur telle quelle."""
        return float(getattr(x, "value", x))

    def _all_floors(self) -> list:
        """Retourne la liste ordonnée des étages contenant des MOB ou des murs béton."""
        floors = list(self.data.keys())
        for wall in self.data_concrete_walls.values():
            if wall["Etage"] not in floors:
                floors.append(wall["Etage"])
        return floors

    def _resolve_etage(self, etage=None) -> str:
        """Résout l'étage demandé (un tuple/liste = choix par défaut -> premier étage existant)."""
        if etage is None or isinstance(etage, (tuple, list)):
            floors = self._all_floors()
            for e in etage or []:
                if e in floors:
                    return e
            return floors[0] if floors else None
        return etage

    @staticmethod
    def _solve_diaphragm(elements: list, Fx: float, Fy: float, XF: float, YF: float, e_acc: float = 0) -> dict:
        """Résout l'équilibre d'un diaphragme rigide sur appuis élastiques (murs de contreventement).

        Chaque élément i est défini par (x_i, y_i, alpha_i, K_i) en unités SI (m, rad, N/m) :
        il ne reprend un effort que dans son plan (direction (cos alpha, sin alpha)).
        Inconnues : translations (u, v) et rotation theta (anti-horaire) de l'origine du repère.
        Déplacement dans le plan de l'élément : delta_i = u cos(a) + v sin(a) + theta * d_i,
        avec d_i = x_i sin(a) - y_i cos(a) le bras de levier signé de l'élément.

        Args:
            elements (list): Liste de tuples (x, y, alpha, K) en SI.
            Fx, Fy (float): Composantes de la force résultante en N.
            XF, YF (float): Point d'application de la force en m.
            e_acc (float): Excentricité accidentelle en m (enveloppe des cas +e_acc et -e_acc).

        Returns:
            dict: delta (m), F (N) par élément, centre de raideur xO, yO, raideurs Kx, Ky, K_theta (au CR),
                rayons de torsion r_x, r_y, excentricités ex, ey, translations u, v et rotation theta du cas de base.
        """
        arr = np.asarray(elements, dtype=float).reshape(-1, 4)
        x, y, al, K = arr.T
        c, s = np.cos(al), np.sin(al)
        d = x * s - y * c
        A, B, C = (K * c * c).sum(), (K * c * s).sum(), (K * s * s).sum()
        D, E, G = (K * c * d).sum(), (K * s * d).sum(), (K * d * d).sum()
        det = A * C - B**2
        if len(arr) < 3 or det <= 1e-9 * (A + C) ** 2:
            raise ValueError(
                "Contreventement instable : il faut au moins 3 murs non tous parallèles "
                "et non concourants pour reprendre les translations et la torsion."
            )
        xO = (A * E - B * D) / det
        yO = (B * E - C * D) / det
        d_cr = d - (xO * s - yO * c)
        K_theta = (K * d_cr**2).sum()
        M3 = np.array([[A, B, D], [B, C, E], [D, E, G]])
        if abs(np.linalg.det(M3)) <= 1e-12 * (A + C) * G:
            raise ValueError("Contreventement instable en torsion : tous les murs concourent en un point.")

        F_res = np.hypot(Fx, Fy)
        cases = []
        for sign in [0] if not e_acc else [1, -1]:
            Mz = XF * Fy - YF * Fx + sign * F_res * e_acc
            u, v, theta = np.linalg.solve(M3, [Fx, Fy, Mz])
            cases.append((u, v, theta, u * c + v * s + theta * d))
        deltas = np.array([case[3] for case in cases])
        idx = np.abs(deltas).argmax(axis=0)
        delta = deltas[idx, np.arange(len(arr))]
        u, v, theta, _ = cases[0]
        return {
            "delta": delta,
            "F": K * delta,
            "xO": xO,
            "yO": yO,
            "Kx": A,
            "Ky": C,
            "K_theta": K_theta,
            "r_x": sqrt(K_theta / C) if C else float("inf"),
            "r_y": sqrt(K_theta / A) if A else float("inf"),
            "ex": XF - xO,
            "ey": YF - yO,
            "u": u,
            "v": v,
            "theta": theta,
        }

    def _init_loads_frames(self):
        """Crée (ou complète) les tableaux d'efforts en conservant les résultats des autres étages."""
        frames = {
            "data_sys_walls_loads": (self._get_all_sys_walls_names(), ["Étage", "delta,v,Ed,i ELS", "Fv,Ed,i ELU"]),
            "data_walls_loads": (list(self.data_walls.keys()), ["Étage", "Fv,Ed,i ELU", "Ft,Ed,i ELU (soulèvement)"]),
            "data_panels_loads": (list(self.data_panels.keys()), ["Étage", "Fv,Ed,i ELU"]),
            "data_concrete_loads": (
                list(self.data_concrete_walls.keys()),
                ["Étage", "delta,v,Ed,i ELS", "Fv,Ed,i ELU", "M,base,Ed ELU"],
            ),
        }
        for attr, (index, columns) in frames.items():
            df = getattr(self, attr, None)
            if df is None:
                df = pd.DataFrame(index=index, columns=columns, dtype=object)
            else:
                df = df.reindex(index).astype(object)
            setattr(self, attr, df)
        if getattr(self, "data_torsion", None) is None:
            self.data_torsion = pd.DataFrame(
                columns=["xO", "yO", "ex", "ey", "Kx", "Ky", "K_theta", "r_x", "r_y", "theta ELS", "xO ELU", "yO ELU"],
                dtype=object,
            )
            self.data_torsion.index.name = "Étage"

    def calculate_loads_and_deformations(
        self,
        Fw_Ed_ELS: si.kN,
        Fw_Ed_ELU: si.kN,
        XFw: si.m = 0,
        YFw: si.m = 0,
        alpha_Fw: float = 0,
        etage: str = Batiment.ETAGE,
        e_acc: si.m = 0,
    ) -> "pd.DataFrame":
        """Calcule les efforts de contreventement et les déformations par système de mur (diaphragme rigide).

        Répartit l'effort global de contreventement entre les systèmes de mur MOB et les murs béton
        (voir add_concrete_wall) en fonction de leurs raideurs dans leur plan et de leurs positions
        par rapport au centre de raideur (centre de torsion), selon la méthode AQCEN MOB généralisée
        (résolution matricielle u, v, theta). Nécessite que le diaphragme (plancher ou toiture) soit
        considéré parfaitement rigide (faibles ouvertures et trémies). Chaque étage est traité
        indépendamment (pas de couplage vertical des refends).

        Deux répartitions sont réalisées :

        - ELS : raideurs Kser des MOB et raideur des murs béton -> déplacements delta,v,Ed,i ELS ;
        - ELU : raideurs Ku = 2/3 Kser des MOB (EN 1995-1-1 §2.2.2(1)b, RATIO_KU_KSER) et raideur
          inchangée des murs béton -> efforts Fv,Ed,i ELU sous Fw_Ed_ELU.

        Pour un contreventement uniquement bois, la répartition ELU est proportionnelle à la
        répartition ELS ; en présence de murs béton, ceux-ci reprennent une part plus importante à l'ELU.

        Args:
            Fw_Ed_ELS (float): Force résultante de contreventement à l'ELS en kN, au niveau de l'étage.
            Fw_Ed_ELU (float): Force résultante de contreventement à l'ELU en kN, au niveau de l'étage.
            XFw (float): Coordonnée X de la force résultante depuis l'origine du repère cartésien
                du calculateur, en m. Defaults to 0.
            YFw (float): Coordonnée Y de la force résultante depuis l'origine du repère cartésien
                du calculateur, en m. Defaults to 0.
            alpha_Fw (float): Angle d'application de la force résultante par rapport au repère cartésien,
                en degrés dans le sens trigonométrique (sens inverse des aiguilles d'une montre). Defaults to 0.
            etage (str): Nom de l'étage où s'applique la force résultante. Defaults to premier étage.
            e_acc (float): Excentricité accidentelle en m (ex. 0,05 x dimension du bâtiment perpendiculaire
                à l'action, EN 1998-1 §4.3.2). Si non nulle, les cas +e_acc et -e_acc sont enveloppés
                (effort maximal en valeur absolue par mur). Defaults to 0.

        Returns:
            pd.DataFrame: Tableau des efforts et déplacements par système de mur avec les colonnes
                ["Étage", "delta,v,Ed,i ELS", "Fv,Ed,i ELU"]. Les résultats des murs béton sont dans
                data_concrete_loads et les caractéristiques de torsion dans data_torsion.
        """
        if not Fw_Ed_ELS:
            raise ValueError("Fw_Ed_ELS doit être non nul.")
        etage = self._resolve_etage(etage)
        F_ELS = Fw_Ed_ELS * 1e3
        F_ELU = Fw_Ed_ELU * 1e3
        alpha_Fw_rad = radians(alpha_Fw)

        Ks = self.get_Kser_wall()
        self._init_loads_frames()

        # Éléments résistants de l'étage : (type, nom, x, y, alpha, K ELS, K ELU) en SI
        elements = []
        for sys_name, sys_wall in self.data.get(etage, {}).items():
            K = self._val(Ks.loc[sys_name, "Kser système de mur"])
            if K > 0:
                elements.append(
                    (
                        "MOB",
                        sys_name,
                        self._val(sys_wall.get("Xg", 0)),
                        self._val(sys_wall.get("Yg", 0)),
                        radians(sys_wall.get("alpha_sw", 0)),
                        K,
                        K * self.RATIO_KU_KSER,
                    )
                )
        for name, wall in self.data_concrete_walls.items():
            if wall["Etage"] == etage:
                elements.append(
                    (
                        "Béton",
                        name,
                        self._val(wall["X"]),
                        self._val(wall["Y"]),
                        radians(wall["alpha"]),
                        self._val(wall["Kser"]),
                        self._val(wall["Kser"]),
                    )
                )
        if not elements:
            raise ValueError(f"Aucun élément de contreventement raide trouvé pour l'étage {etage}.")

        res = self._solve_diaphragm(
            [e[2:6] for e in elements],
            F_ELS * cos(alpha_Fw_rad),
            F_ELS * sin(alpha_Fw_rad),
            XFw,
            YFw,
            e_acc,
        )
        res_ELU = self._solve_diaphragm(
            [(*e[2:5], e[6]) for e in elements],
            F_ELU * cos(alpha_Fw_rad),
            F_ELU * sin(alpha_Fw_rad),
            XFw,
            YFw,
            e_acc,
        )
        self.data_torsion.loc[etage] = [
            res["xO"] * si.m,
            res["yO"] * si.m,
            res["ex"] * si.m,
            res["ey"] * si.m,
            res["Kx"] * si.N / si.m,
            res["Ky"] * si.N / si.m,
            res["K_theta"] * si.N * si.m,
            res["r_x"] * si.m,
            res["r_y"] * si.m,
            res["theta"],
            res_ELU["xO"] * si.m,
            res_ELU["yO"] * si.m,
        ]
        self._last_load = {
            "Etage": etage,
            "Fw_Ed_ELS": Fw_Ed_ELS,
            "Fw_Ed_ELU": Fw_Ed_ELU,
            "XFw": XFw,
            "YFw": YFw,
            "alpha_Fw": alpha_Fw,
            "e_acc": e_acc,
        }

        for (kind, name, *_), delta, F in zip(elements, res["delta"], res_ELU["F"]):
            Fv_ELU = F * si.N
            if kind == "Béton":
                h = self.data_concrete_walls[name]["Hauteur"]
                self.data_concrete_loads.loc[name] = [etage, delta * si.m, Fv_ELU, Fv_ELU * h]
                continue

            self.data_sys_walls_loads.loc[name] = [etage, delta * si.m, Fv_ELU]
            Ksw = Ks.loc[name, "Kser système de mur"]
            for wall_name, wall in self.data[etage][name]["Murs"].items():
                Kw = self.data_Kser_walls.loc[wall_name, "Kser mur"]
                Fv_wall = Fv_ELU * Kw / Ksw
                self.data_walls_loads.loc[wall_name, "Étage"] = etage
                self.data_walls_loads.loc[wall_name, "Fv,Ed,i ELU"] = Fv_wall
                # traction dans le montant d'extrémité
                self.data_walls_loads.loc[wall_name, "Ft,Ed,i ELU (soulèvement)"] = self._FtEd(wall_name)

                for panel_name, panel in wall["Panneaux"].items():
                    if not panel["Contreventement"]:
                        continue
                    Kp = self.data_Kser_panels.loc[panel_name, "Kser,p unit."]
                    self.data_panels_loads.loc[panel_name, "Étage"] = etage
                    self.data_panels_loads.loc[panel_name, "Fv,Ed,i ELU"] = Fv_wall * Kp / Kw

        return self.data_sys_walls_loads

    def synthese_torsion(self) -> "pd.DataFrame":
        """Retourne la synthèse de torsion par étage (après calculate_loads_and_deformations).

        Colonnes : centre de raideur (xO, yO), excentricités de la force (ex, ey), raideurs de translation
        globales (Kx, Ky), raideur de torsion au centre de raideur (K_theta), rayons de torsion
        (r_x = sqrt(K_theta/Ky), r_y = sqrt(K_theta/Kx), cf. EN 1998-1 §4.2.3.2) et rotation ELS du diaphragme,
        calculés avec les raideurs ELS (Kser). xO ELU, yO ELU : centre de raideur avec Ku = 2/3 Kser pour
        les MOB (différent du CR ELS uniquement en présence de murs béton).

        Returns:
            pd.DataFrame: Tableau data_torsion indexé par étage.
        """
        if getattr(self, "data_torsion", None) is None:
            raise ValueError("Lancer calculate_loads_and_deformations avant synthese_torsion.")
        return self.data_torsion

    # ------------------------------------------------------------------
    # Murs béton (refends)
    # ------------------------------------------------------------------

    @property
    def data_concrete_walls(self) -> dict:
        """Retourne le dictionnaire des murs béton (refends) participant au contreventement."""
        if not hasattr(self, "_data_concrete"):
            self._data_concrete = {}
        return self._data_concrete

    def add_concrete_wall(
        self,
        l_wall: si.m,
        t_wall: si.m,
        h_wall: si.m,
        X: si.m = 0,
        Y: si.m = 0,
        alpha_w: float = 0,
        etage: str = Batiment.ETAGE,
        classe_beton: str = "C25/30",
        Ecm: si.MPa = None,
        coeff_fissuration: float = 0.5,
        nu: float = 0.2,
        K_override: si.kN / si.mm = None,
    ) -> str:
        """Ajoute un mur béton (refend, voile) participant au contreventement de l'étage.

        Les murs béton sont nommés automatiquement B1, B2, … Seule la raideur dans le plan du mur
        est prise en compte (raideurs hors plan et de torsion propre négligées).
        Par défaut la raideur est celle d'une console encastrée en pied (flexion + cisaillement) :

            K = 1 / (h³ / (3·E·I) + 1,2·h / (G·A)),  I = t·l³/12,  A = t·l,  G = E / (2(1+ν)),
            E = coeff_fissuration · Ecm

        Pour un voile tenu en tête par un plancher rigide sur plusieurs niveaux ou un voile percé,
        renseigner K_override (raideur issue d'un modèle dédié).

        Args:
            l_wall (float): Longueur du mur dans son plan en m.
            t_wall (float): Épaisseur du mur en m.
            h_wall (float): Hauteur du mur (de l'encastrement au diaphragme) en m.
            X (float): Coordonnée X du milieu du mur depuis l'origine du repère du calculateur, en m.
            Y (float): Coordonnée Y du milieu du mur depuis l'origine du repère du calculateur, en m.
            alpha_w (float): Orientation du mur dans le sens trigonométrique, en degrés. Defaults to 0.
            etage (str): Nom de l'étage auquel appartient le mur. Defaults to premier étage.
            classe_beton (str): Classe de résistance du béton (ex. "C25/30"), utilisée pour Ecm.
            Ecm (float, optional): Module sécant du béton en MPa. Écrase la valeur de la classe. Defaults to None.
            coeff_fissuration (float): Coefficient de réduction de raideur pour fissuration
                (0,5 usuel, EN 1998-1 §4.3.1(7) ; 1,0 si béton non fissuré). Defaults to 0.5.
            nu (float): Coefficient de Poisson du béton. Defaults to 0.2.
            K_override (float, optional): Raideur dans le plan en kN/mm. Si fournie, remplace le calcul.

        Returns:
            str: Nom attribué au mur béton (ex. : "B1").
        """
        etage = etage[0] if isinstance(etage, (tuple, list)) else etage
        if Ecm is None:
            Ecm = float(self._data_from_csv("caracteristique_meca_beton.csv").loc[classe_beton, "Ecm"])
        if K_override:
            K = K_override * 1e6
        else:
            if l_wall <= 0 or t_wall <= 0 or h_wall <= 0:
                raise ValueError("l_wall, t_wall et h_wall doivent être strictement positifs.")
            E = coeff_fissuration * Ecm * 1e6
            G = E / (2 * (1 + nu))
            I = t_wall * l_wall**3 / 12
            A = t_wall * l_wall
            K = 1 / (h_wall**3 / (3 * E * I) + 1.2 * h_wall / (G * A))

        name = f"B{len(self.data_concrete_walls) + 1}"
        self.data_concrete_walls[name] = {
            "Etage": etage,
            "Longueur": l_wall * si.m,
            "Epaisseur": t_wall * si.m,
            "Hauteur": h_wall * si.m,
            "X": X * si.m,
            "Y": Y * si.m,
            "alpha": alpha_w,
            "Classe béton": classe_beton,
            "Ecm": Ecm * si.MPa,
            "Coeff. fissuration": coeff_fissuration,
            "Kser": K * si.N / si.m,
        }
        return name

    # ------------------------------------------------------------------
    # Visualisation
    # ------------------------------------------------------------------

    def _sys_wall_geometry(self, sys_name: str, sys_wall: dict) -> dict:
        """Calcule la géométrie globale d'un système de mur (origine, direction, murs).

        L'origine locale (coin gauche) est déduite de (Xg, Yg), considéré comme le centre de raideur
        du système : x_c,local = Σ(Kw · (position + L/2)) / ΣKw (pondération par les longueurs à défaut).
        """
        alpha = radians(sys_wall.get("alpha_sw", 0))
        ux, uy = cos(alpha), sin(alpha)
        weights, centers, segments = [], [], {}
        Kw_df = getattr(self, "data_Kser_walls", None)
        for wall_name, wall in sys_wall["Murs"].items():
            x0 = self._val(wall["Position en X"])
            L = self._val(wall["Longueur"])
            segments[wall_name] = (x0, x0 + L)
            Kw = Kw_df.loc[wall_name, "Kser mur"] if Kw_df is not None and wall_name in Kw_df.index else 0
            weights.append(self._val(Kw))
            centers.append(x0 + L / 2)
        if not weights or sum(weights) <= 0:
            weights = [b - a for a, b in segments.values()]
        xc = float(np.average(centers, weights=weights)) if centers and sum(weights) > 0 else 0
        Xg, Yg = self._val(sys_wall.get("Xg", 0)), self._val(sys_wall.get("Yg", 0))
        ox, oy = Xg - xc * ux, Yg - xc * uy

        def to_global(t):
            return ox + t * ux, oy + t * uy

        return {
            "origin": (ox, oy),
            "dir": (ux, uy),
            "center": (Xg, Yg),
            "length": self._val(sys_wall["Longueur du système de mur"]),
            "walls": {name: (to_global(a), to_global(b)) for name, (a, b) in segments.items()},
            "to_global": to_global,
        }

    def show_plan(
        self,
        etage: str = None,
        show_forces: bool = True,
        show_labels: bool = True,
        ax=None,
        filepath: str = None,
        show: bool = True,
    ):
        """Affiche la vue en plan du contreventement d'un étage (matplotlib).

        Représente les systèmes de mur MOB (murs pleins en trait épais, ouvertures en pointillés),
        les murs béton (hachurés), le centre de raideur / torsion (CR), la force résultante appliquée
        et son excentricité, ainsi que les efforts ELU repris par chaque élément (flèches dans leur plan).
        Si taux_walls a été appelé, les murs MOB sont colorés selon leur taux de cisaillement
        (vert < 70 %, orange < 100 %, rouge ≥ 100 %).

        Args:
            etage (str, optional): Étage à représenter. Defaults to dernier étage calculé, sinon premier étage.
            show_forces (bool): Affiche la force appliquée et les efforts repris. Defaults to True.
            show_labels (bool): Affiche les noms et valeurs. Defaults to True.
            ax (matplotlib.axes.Axes, optional): Axe existant sur lequel tracer. Defaults to None.
            filepath (str, optional): Chemin d'export de l'image (PNG, PDF, SVG…). Defaults to None.
            show (bool): Appelle plt.show(). Defaults to True.

        Returns:
            tuple: (fig, ax) matplotlib.
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon

        last = getattr(self, "_last_load", None)
        if etage is None and last:
            etage = last["Etage"]
        etage = self._resolve_etage(etage)
        if etage is None:
            raise ValueError("Aucun mur à représenter.")
        if getattr(self, "data_Kser_walls", None) is None and self.data:
            self.get_Kser_wall()

        fig, ax = plt.subplots(figsize=(10, 8)) if ax is None else (ax.figure, ax)
        taux = getattr(self, "taux_FvRd_walls", None)
        sys_loads = getattr(self, "data_sys_walls_loads", None)
        concrete_loads = getattr(self, "data_concrete_loads", None)
        loaded = bool(show_forces and last and last["Etage"] == etage)

        def color_taux(t):
            return "tab:green" if t < 70 else ("tab:orange" if t < 100 else "tab:red")

        def load_of(df, name):
            if df is None or name not in df.index or df.loc[name, "Étage"] != etage:
                return None
            F = df.loc[name, "Fv,Ed,i ELU"]
            return None if F is None or (isinstance(F, float) and np.isnan(F)) else self._val(F)

        points, arrows = [], []
        for sys_name, sys_wall in self.data.get(etage, {}).items():
            geo = self._sys_wall_geometry(sys_name, sys_wall)
            p0, p1 = geo["to_global"](0), geo["to_global"](geo["length"])
            points += [p0, p1]
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], ls=":", color="grey", lw=1.5, zorder=1)
            for wall_name, (a, b) in geo["walls"].items():
                color = "saddlebrown"
                if taux is not None and wall_name in taux.index:
                    color = color_taux(abs(taux.loc[wall_name, "Taux cisaillement (%)"]))
                ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=6, solid_capstyle="butt", zorder=2)
                if show_labels:
                    ax.annotate(
                        wall_name,
                        ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2),
                        xytext=(0, -12),
                        textcoords="offset points",
                        ha="center",
                        fontsize=7,
                        color=color,
                    )
            ax.plot(*geo["center"], marker="o", ms=4, color="saddlebrown", zorder=3)
            if show_labels:
                ax.annotate(
                    sys_name, geo["center"], xytext=(4, 6), textcoords="offset points", fontsize=9, weight="bold"
                )
            F = load_of(sys_loads, sys_name) if loaded else None
            if F is not None:
                arrows.append((geo["center"], geo["dir"], F))

        for name, wall in self.data_concrete_walls.items():
            if wall["Etage"] != etage:
                continue
            al = radians(wall["alpha"])
            ux, uy, nx, ny = cos(al), sin(al), -sin(al), cos(al)
            X, Y = self._val(wall["X"]), self._val(wall["Y"])
            L, t = self._val(wall["Longueur"]) / 2, self._val(wall["Epaisseur"]) / 2
            corners = [
                (X + i * L * ux + j * t * nx, Y + i * L * uy + j * t * ny)
                for i, j in ((-1, -1), (1, -1), (1, 1), (-1, 1))
            ]
            points += corners
            ax.add_patch(
                Polygon(corners, closed=True, facecolor="lightgrey", edgecolor="dimgrey", hatch="///", zorder=2)
            )
            if show_labels:
                ax.annotate(
                    name, (X, Y), xytext=(4, 6), textcoords="offset points", fontsize=9, weight="bold", color="dimgrey"
                )
            F = load_of(concrete_loads, name) if loaded else None
            if F is not None:
                arrows.append(((X, Y), (ux, uy), F))

        torsion = getattr(self, "data_torsion", None)
        cr = cr_ELU = None
        if torsion is not None and etage in torsion.index:
            cr = (self._val(torsion.loc[etage, "xO"]), self._val(torsion.loc[etage, "yO"]))
            cr_ELU = (self._val(torsion.loc[etage, "xO ELU"]), self._val(torsion.loc[etage, "yO ELU"]))
            if np.hypot(cr_ELU[0] - cr[0], cr_ELU[1] - cr[1]) < 0.01:
                cr_ELU = None
            points += [p for p in (cr, cr_ELU) if p]

        if loaded:
            points.append((last["XFw"], last["YFw"]))
        pts = np.array(points) if points else np.zeros((1, 2))
        extent = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]), 1.0)

        if arrows:
            F_max = max(abs(F) for *_, F in arrows) or 1
            for (x, y), (ux, uy), F in arrows:
                k = 0.15 * extent * F / F_max
                ax.annotate(
                    "",
                    xy=(x + k * ux, y + k * uy),
                    xytext=(x, y),
                    arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1.5),
                    zorder=4,
                )
                if show_labels:
                    ax.annotate(
                        f"{F / 1e3:.1f} kN",
                        (x + k * ux, y + k * uy),
                        xytext=(3, 3),
                        textcoords="offset points",
                        fontsize=7,
                        color="tab:blue",
                    )

        if cr is not None:
            ax.plot(*cr, marker="X", ms=12, color="tab:purple", zorder=5, label="Centre de raideur / torsion")
            if show_labels:
                ax.annotate("CR", cr, xytext=(6, -12), textcoords="offset points", color="tab:purple", weight="bold")
        if cr_ELU is not None:
            ax.plot(*cr_ELU, marker="x", ms=9, color="tab:pink", zorder=5, label="Centre de raideur ELU (Ku)")

        if loaded:
            XF, YF, a = last["XFw"], last["YFw"], radians(last["alpha_Fw"])
            L = 0.25 * extent
            ax.annotate(
                "",
                xy=(XF + L * cos(a), YF + L * sin(a)),
                xytext=(XF, YF),
                arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=2.5, mutation_scale=20),
                zorder=6,
            )
            ax.plot(XF, YF, marker="o", color="tab:red", zorder=6, label="Force résultante Fw")
            if show_labels:
                ax.annotate(
                    f"Fw,ELU = {last['Fw_Ed_ELU']:.1f} kN\nFw,ELS = {last['Fw_Ed_ELS']:.1f} kN",
                    (XF + L * cos(a), YF + L * sin(a)),
                    xytext=(5, 5),
                    textcoords="offset points",
                    color="tab:red",
                    fontsize=8,
                )
            if cr is not None:
                ax.plot([cr[0], XF], [cr[1], YF], ls="--", color="tab:purple", lw=1, zorder=4, label="Excentricité")

        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(True, ls=":", alpha=0.6)
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        ax.set_title(f"Contreventement - {etage}")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="best", fontsize=8)
        if filepath:
            fig.savefig(filepath, dpi=150, bbox_inches="tight")
        if show:
            plt.show()
        return fig, ax

    @property
    def data(self):
        return self._data_MOB

    @data.setter
    def data(self, data: dict):
        self._data_MOB = data

    @property
    def data_sys_walls(self) -> dict:
        """Retourne un dictionnaire des systèmes de mur sans leur sous-clé 'Murs'.

        Returns:
            dict: Clés = noms des systèmes de mur, valeurs = dict de leurs propriétés
                (étage, hauteurs, longueur, position, efforts éventuels).
        """
        sys_wall_data = {}
        _data = deepcopy(self._data_MOB)
        for etage, sys in _data.items():
            for sys_name, sys_wall in sys.items():
                sys_wall.pop("Murs", None)
                sys_wall_data[sys_name] = sys_wall
        return sys_wall_data

    @property
    def data_walls(self) -> dict:
        """Retourne un dictionnaire des murs sans leur sous-clé 'Panneaux'.

        Returns:
            dict: Clés = noms des murs (ex. "S1_W1"), valeurs = dict de leurs propriétés
                (position, longueur, couturage, connecteurs).
        """
        wall_data = {}
        _data = deepcopy(self._data_MOB)
        for etage, sys in _data.items():
            for sys_name, sys_wall in sys.items():
                for wall_name, wall in sys_wall["Murs"].items():
                    wall.pop("Panneaux", None)
                    wall_data[wall_name] = wall
        return wall_data

    @property
    def data_panels(self) -> dict:
        """Retourne un dictionnaire de tous les panneaux de tous les murs et systèmes.

        Returns:
            dict: Clés = noms des panneaux (ex. "S1_W1_P1"), valeurs = dict de leurs
                propriétés (type, position, dimensions, nombre, contreventement, panneau double).
        """
        panel_data = {}
        for etage, sys in self._data_MOB.items():
            for sys_name, sys_wall in sys.items():
                for wall_name, wall in sys_wall["Murs"].items():
                    for panel_name, panel in wall["Panneaux"].items():
                        panel_data[panel_name] = panel
        return panel_data

    def add_sys_wall(
        self,
        h_etage: si.m,
        h_sys_MOB: si.m,
        l_sys_MOB: si.m,
        etage: str = Batiment.ETAGE,
        Fv_Ed_ELS: si.kN = None,
        Fv_Ed_ELU: si.kN = None,
        Xg: si.m = 0,
        Yg: si.m = 0,
        alpha_sw: float = 0,
    ) -> str:
        """Ajoute un système de mur à ossature bois (contreventement complet) au niveau d'un étage.

        Les systèmes de mur sont nommés automatiquement S1, S2, …
        Le système de mur est défini comme un repère cartésien, coin inférieur gauche
        en (0, 0), la face externe définit le sens positif des abscisses.

        Args:
            h_etage (float): Hauteur de l'étage en m.
            h_sys_MOB (float): Hauteur du système de mur en m (hauteur libre entre lisse basse et lisse haute).
            l_sys_MOB (float): Longueur totale du système de mur en m (somme des murs et des ouvertures).
            etage (str): Nom de l'étage auquel appartient le système de mur.
            Fv_Ed_ELS (float, optional): Effort de contreventement de calcul à l'ELS en kN.
                À renseigner uniquement pour écraser la valeur calculée par calculate_loads_and_deformations.
                Defaults to None.
            Fv_Ed_ELU (float, optional): Effort de contreventement de calcul à l'ELU en kN.
                À renseigner uniquement pour écraser la valeur calculée par calculate_loads_and_deformations.
                Defaults to None.
            Xg (float, optional): Coordonnée X du centre de rigidité du système de mur depuis l'origine
                du repère cartésien du calculateur, en m. Defaults to 0.
            Yg (float, optional): Coordonnée Y du centre de rigidité du système de mur depuis l'origine
                du repère cartésien du calculateur, en m. Defaults to 0.
            alpha_sw (float, optional): Angle d'orientation du système de mur dans le sens trigonométrique
                (anti-horaire) depuis l'origine, en degrés. Defaults to 0.

        Returns:
            str: Nom attribué au système de mur (ex. : "S1", "S2").
        """
        # Nommer les systèmes de mur de façon continue indépendamment des étages
        total_sys = sum(len(sys) for sys in self.data.values()) if self.data else 0
        name = f"S{total_sys + 1}"

        sys_wall = {
            "Etage": etage,
            "Hauteur étage": h_etage * si.m,
            "Hauteur du système de mur": h_sys_MOB * si.m,
            "Longueur du système de mur": l_sys_MOB * si.m,
            "Murs": {},
        }
        if Xg != 0 or Yg != 0 or alpha_sw != 0:
            sys_wall["Xg"] = Xg * si.m
            sys_wall["Yg"] = Yg * si.m
            sys_wall["alpha_sw"] = alpha_sw

        if Fv_Ed_ELU is not None:
            sys_wall["Fv_Ed_ELU"] = Fv_Ed_ELU * si.kN
        if Fv_Ed_ELS is not None:
            sys_wall["Fv_Ed_ELS"] = Fv_Ed_ELS * si.kN

        if self.data.get(etage) is None:
            self.data[etage] = {}
        self.data[etage][name] = sys_wall
        return name

    def add_wall_to_sys(
        self,
        sys_name: str,
        position: si.m,
        l_MOB: si.m,
        couturage_ext_rive: si.mm = 150,
        couturage_ext_inter: si.mm = 300,
        couturage_int_rive: si.mm = 0,
        couturage_int_inter: si.mm = 0,
        Ff_Rk_ext: si.N = None,
        Ff_Rk_int: si.N = None,
        Kser_ext: si.N / si.mm = None,
        Kser_int: si.N / si.mm = None,
    ) -> str:
        """Ajoute un mur contreventant interne à un système de mur à ossature bois.

        Les murs sont nommés automatiquement {sys_name}_W1, {sys_name}_W2, …
        Le système de mur est défini comme un repère cartésien, coin inférieur gauche
        en (0, 0), la face externe définit le sens positif des abscisses.

        Args:
            sys_name (str): Nom du système de mur auquel le mur est rattaché.
            position (float): Coordonnée X du coin gauche du mur dans le système de mur, en m.
            l_MOB (float): Longueur du mur en m (sans les ouvertures).
            couturage_ext_rive (float): Entraxe des organes en rive sur la face extérieure du panneau, en mm.
                Defaults to 150 mm.
            couturage_ext_inter (float): Entraxe des organes intermédiaires sur la face extérieure, en mm.
                Defaults to 300 mm.
            couturage_int_rive (float): Entraxe des organes en rive sur la face intérieure du panneau, en mm.
                0 si pas de panneau intérieur. Defaults to 0.
            couturage_int_inter (float): Entraxe des organes intermédiaires sur la face intérieure, en mm.
                Defaults to 0.
            Ff_Rk_ext (float, optional): Capacité résistante caractéristique de cisaillement unitaire du connecteur
                extérieur en N (valeur brute, sans Kmod ni gamma_M). Écrase la valeur issue de connecteurs_ext.
                Defaults to None.
            Ff_Rk_int (float, optional): Capacité résistante caractéristique de cisaillement unitaire du connecteur
                intérieur en N (valeur brute, sans Kmod ni gamma_M). Écrase la valeur issue de connecteurs_int.
                Defaults to None.
            Kser_ext (float, optional): Raideur ELS du connecteur extérieur en N/mm.
                Écrase la valeur issue de connecteurs_ext. Defaults to None.
            Kser_int (float, optional): Raideur ELS du connecteur intérieur en N/mm.
                Écrase la valeur issue de connecteurs_int. Defaults to None.

        Returns:
            str: Nom attribué au mur (ex. : "S1_W1", "S1_W2").

        Raises:
            ValueError: Si sys_name ne correspond à aucun système de mur enregistré.
        """
        etage_found = None
        for etage, sys in self.data.items():
            if sys_name in sys.keys():
                etage_found = etage
                break
        if etage_found is None:
            raise ValueError(f"Le système de mur {sys_name} n'existe pas")

        name = f"{sys_name}_W{len(self.data[etage_found][sys_name]['Murs']) + 1}"
        wall = {"Position en X": position * si.m, "Longueur": l_MOB * si.m, "Panneaux": {}}

        if couturage_ext_rive:
            wall["Couturage ext. rive"] = couturage_ext_rive * si.mm
            wall["Couturage ext. intermédiaire"] = couturage_ext_inter * si.mm
            if Ff_Rk_ext is None or Ff_Rk_ext == 0 * si.N:
                wall["Ff,Rk ext"] = self._data_connecteurs["Extérieur"]["Ff,Rk"]
            else:
                wall["Ff,Rk ext"] = Ff_Rk_ext * si.N
            if Kser_ext is None or Kser_ext == 0 * si.N / si.mm:
                wall["Kser ext"] = self._data_connecteurs["Extérieur"]["Kser"]
            else:
                wall["Kser ext"] = Kser_ext * si.N / si.mm

        if couturage_int_rive:
            wall["Couturage int. rive"] = couturage_int_rive * si.mm
            wall["Couturage int. intermédiaire"] = couturage_int_inter * si.mm
            if Ff_Rk_int is None or Ff_Rk_int == 0 * si.N:
                wall["Ff,Rk int"] = self._data_connecteurs["Intérieur"]["Ff,Rk"]
            else:
                wall["Ff,Rk int"] = Ff_Rk_int * si.N
            if Kser_int is None or Kser_int == 0 * si.N / si.mm:
                wall["Kser int"] = self._data_connecteurs["Intérieur"]["Kser"]
            else:
                wall["Kser int"] = Kser_int * si.N / si.mm

        self.data[etage_found][sys_name]["Murs"][name] = wall
        return name

    def add_panel_to_wall(
        self,
        wall_name: str,
        h_panel: si.mm,
        b_panel: si.mm,
        number: int = 1,
        position_panel: str = POSITION_PANEL,
        on_top_of: str = None,
    ) -> str:
        """Ajoute un panneau de contreventement à un mur interne.

        Les panneaux sont nommés automatiquement {wall_name}_P1, {wall_name}_P2, …
        Un panneau est automatiquement marqué non contreventant si son rapport
        hauteur/largeur (ou largeur/hauteur) dépasse 4 (EN 1995-1-1 §9.2.4.2).
        Les panneaux empilés verticalement (on_top_of) sont traités en série de raideur.

        Args:
            wall_name (str): Nom du mur auquel le panneau est rattaché.
            h_panel (float): Hauteur du panneau en mm (généralement la dimension verticale).
            b_panel (float): Largeur du panneau en mm (généralement la dimension horizontale).
            number (int): Nombre de panneaux identiques côte à côte. Defaults to 1.
            position_panel (str): Face du panneau par rapport à l'ossature.
                "Extérieur" : face extérieure du mur.
                "Intérieur" : face intérieure du mur.
            on_top_of (str, optional): Nom du panneau inférieur si ce panneau est posé
                en empilé vertical (panneaux en H/2). Defaults to None.

        Returns:
            str: Nom attribué au panneau (ex. : "S1_W1_P1").

        Raises:
            ValueError: Si le système de mur parent est introuvable dans data.
        """
        sys_name = wall_name.split("_")[0]
        etage_found = None
        for etage, sys in self.data.items():
            if sys_name in sys.keys():
                etage_found = etage
                break
        if etage_found is None:
            raise ValueError(f"Le système de mur {sys_name} n'existe pas")

        name = f"{wall_name}_P{len(self.data[etage_found][sys_name]['Murs'][wall_name]['Panneaux']) + 1}"
        wall_data = self.data[etage_found][sys_name]["Murs"][wall_name]

        panel = {
            "Type": self._connecteurs[position_panel].beam_1.classe,
            "Position": position_panel,
            "Panneau double": "Non",
            "Contreventement": True,
            "Nombre": number,
            "Epaisseur": self._connecteurs[position_panel].beam_1.b,
            "Hauteur": h_panel * si.mm,
            "Largeur": b_panel * si.mm,
        }

        # On test si le panneau est contreventant
        if h_panel > b_panel:
            if h_panel / 4 > b_panel:
                panel["Contreventement"] = False
        else:
            if b_panel / 4 > h_panel:
                panel["Contreventement"] = False

        if on_top_of:
            panel["On top of"] = on_top_of
            self._list_on_top_of.append(name)

        # Si le mur est composé potentiellement de deux panneaux, vérifier que les 2 côtés existent réellement
        if (
            self._connecteurs["Intérieur"]
            and self._connecteurs["Extérieur"]
            and ("Kser ext" in wall_data)
            and ("Kser int" in wall_data)
        ):
            # Si les Kser sont différents
            if wall_data["Kser ext"] != wall_data["Kser int"]:
                panel["Panneau double"] = "Panneau différent & glissement différent"
            # Si les classes sont identiques
            elif self._connecteurs["Intérieur"].beam_1.classe == self._connecteurs["Extérieur"].beam_1.classe:
                panel["Panneau double"] = "Panneau différent & glissement identique"
                # Si les largeurs sont identiques
                if self._connecteurs["Intérieur"].beam_1.b == self._connecteurs["Extérieur"].beam_1.b:
                    if (
                        ("Ff,Rk ext" in wall_data)
                        and ("Ff,Rk int" in wall_data)
                        and (wall_data["Ff,Rk ext"] != wall_data["Ff,Rk int"])
                    ):
                        panel["Panneau double"] = "Panneau et organe identiques"
            else:
                panel["Panneau double"] = "Panneau différent & glissement identique"
        self.data[etage_found][sys_name]["Murs"][wall_name]["Panneaux"][name] = panel
        return name

    def save_walls_data(self, path: str = None) -> None:
        """Sauvegarde la structure complète des murs MOB dans un fichier JSON.

        Sérialise le dictionnaire data (systèmes de mur, murs, panneaux) via le
        mécanisme Batiment.save_data. Ce fichier peut ensuite être rechargé avec load_walls_data.

        Args:
            path (str, optional): Chemin absolu du fichier de destination (.json).
                Si None, une boîte de dialogue s'ouvre pour sélectionner le fichier.
                Defaults to None.
        """
        data = {"_MOB": self.data, "_concrete": self.data_concrete_walls} if self.data_concrete_walls else self.data
        super().save_data(data, type_data="JSON", path=path)

    def load_walls_data(self, path: str = None) -> dict:
        """Charge la structure des murs MOB depuis un fichier JSON produit par save_walls_data.

        Reconstruit le dictionnaire data interne (_data_MOB) à partir du fichier JSON.

        Args:
            path (str, optional): Chemin absolu du fichier source (.json).
                Si None, une boîte de dialogue s'ouvre pour sélectionner le fichier.
                Defaults to None.

        Returns:
            dict: Dictionnaire data structuré par étage → système de mur → mur → panneau.
        """
        data = super().load_data(type_data="JSON", path=path)
        if "_concrete" in data:
            self._data_concrete = data["_concrete"]
            data = data["_MOB"]
        self.data = data
        self._list_on_top_of = [name for name, panel in self.data_panels.items() if panel.get("On top of")]
        return self.data

    def _get_Kser_panel(self, b: float, h: float, k_ser_r: float, k_ser_i: float):
        """
        Calcul la raideur du panneau.

        Args:
            b (float): largeur du panneau.
            h (float): hauteur du panneau.
            k_ser_r (float): raideur du panneau.
            k_ser_i (float): raideur du panneau.
        """
        if b == 0 * si.mm or h == 0 * si.mm:
            raise ValueError("Dimensions de panneau invalides: b et h doivent être non nulles.")
        if k_ser_r == 0 * si.N / si.mm**2:
            raise ValueError("Raideur linéique de rive nulle (k_ser_r): vérifier Kser et les couturages de rive.")
        K_p = (2 / (k_ser_r * (b + h / 3) + k_ser_i * h / 6) + 2 * h**2 / (k_ser_r * b**2 * (h + b / 3))) ** -1
        return K_p

    def _get_Kser_panel_stack(self, top_panel_name: str, bottom_panel_name: str):
        """
        Calcul la raideur du panneau en pile.

        Args:
            top_panel_name (str): nom du panneau supérieur.
            bottom_panel_name (str): nom du panneau inférieur.
        """
        k2 = self.data_Kser_panels.loc[top_panel_name]["Kser,p unit."]
        k1 = self.data_Kser_panels.loc[bottom_panel_name]["Kser,p unit."]
        k_equi = (k1 * k2) / (k1 + k2)
        return k_equi

    def _k_panel(self):
        """
        Calcul la raideur des panneaux.
        """
        df_kp = pd.DataFrame(
            columns=[
                "Étage",
                "Type",
                "Épaisseur",
                "Largeur",
                "Hauteur",
                "Nombre",
                "Kser,p unit.",
                "Kser,p,ref",
                "Coeff. équivalence",
            ]
        )
        df_kp.index.name = "Panneau"
        K_p_ref = {}
        for etage, sys in self.data.items():
            for sys_name, sys_wall in sys.items():
                h_sys = sys_wall["Hauteur du système de mur"]
                for wall_name, wall in sys_wall["Murs"].items():
                    if wall.get("Kser ext", None) is not None:
                        kser_ext = wall["Kser ext"]
                        couturage_ext_rive = wall["Couturage ext. rive"]
                        couturage_ext_inter = wall["Couturage ext. intermédiaire"]
                    if wall.get("Kser int", None) is not None:
                        kser_int = wall["Kser int"]
                        couturage_int_rive = wall["Couturage int. rive"]
                        couturage_int_inter = wall["Couturage int. intermédiaire"]

                    for panel_name, panel in wall["Panneaux"].items():
                        # On limite la largeur efficace de contreventement des panneaux
                        if not panel["Contreventement"]:
                            continue

                        if panel["Position"] == "Extérieur":
                            k_ser_r = kser_ext / couturage_ext_rive
                            k_ser_i = kser_ext / couturage_ext_inter
                        else:
                            k_ser_r = kser_int / couturage_int_rive
                            k_ser_i = kser_int / couturage_int_inter

                        h = panel["Hauteur"]
                        b = panel["Largeur"]
                        K_p = self._get_Kser_panel(b, h, k_ser_r, k_ser_i)

                        df_kp.loc[panel_name] = [
                            etage,
                            panel["Type"],
                            panel["Epaisseur"],
                            b,
                            h,
                            panel["Nombre"],
                            K_p,
                            None,
                            None,
                        ]
                        if K_p > K_p_ref.get(etage, 0 * si.N / si.mm):
                            K_p_ref[etage] = K_p

        # Trier par index (nommé 'Panneau') au lieu d'une colonne inexistante
        df_kp.sort_index(inplace=True)
        for panel_name, panel in df_kp.iterrows():
            etage = panel["Étage"]
            df_kp.loc[panel_name, "Kser,p,ref"] = K_p_ref.get(etage, 0 * si.N / si.mm)
            C_eq = panel["Kser,p unit."] / K_p_ref.get(etage, 0 * si.N / si.mm)
            df_kp.loc[panel_name, "Coeff. équivalence"] = C_eq

        return df_kp

    def get_Kser_wall(self) -> "pd.DataFrame":
        """Calcule et retourne la raideur de glissement ELS de chaque système de mur.

        Agrège les raideurs de panneaux (Kser,p) et de murs selon la méthode de la
        section §9.2.4 de l'EN 1995-1-1. Peuple également les attributs
        data_Kser_panels et data_Kser_walls.

        Returns:
            pd.DataFrame: Tableau data_Kser_sys_walls indexé par le nom du système de mur,
                avec les colonnes ["Étage", "Kser système de mur", "Coeff. équivalence", "Kser,p,ref"]
                (Kser en N/mm).
        """
        self.data_Kser_panels = self._k_panel()
        # Colonnes alignées avec les valeurs calculées: [Étage, Kser, Coeff, Kser,p,ref]
        self.data_Kser_sys_walls = pd.DataFrame(
            columns=["Étage", "Kser système de mur", "Coeff. équivalence", "Kser,p,ref"]
        )
        self.data_Kser_sys_walls.index.name = "Système de mur"
        self.data_Kser_walls = pd.DataFrame(columns=["Étage", "Kser mur", "Coeff. équivalence", "Kser,p,ref"])
        self.data_Kser_walls.index.name = "Mur"

        for etage, sys in self.data.items():
            for sys_name, sys_wall in sys.items():
                k_sys_wall = [etage, 0, 0, 0]

                for wall_name, wall in sys_wall["Murs"].items():
                    k_wall = [etage, 0, 0, 0]

                    # on gère d'abord les panneaux en pile
                    already_reported = []
                    for top_panel_name in self._list_on_top_of:
                        if top_panel_name not in wall["Panneaux"].keys():
                            continue

                        top_panel = wall["Panneaux"][top_panel_name]
                        bottom_panel_name = top_panel["On top of"]
                        bottom_panel = wall["Panneaux"][bottom_panel_name]

                        # on vérifie que les deux panneaux sont des panneaux de contreventement, si ce n'est pas le cas alors,
                        # les panneaux sont considérés comme des panneaux simples
                        if not top_panel["Contreventement"] or not bottom_panel["Contreventement"]:
                            continue
                        kser_equi = self._get_Kser_panel_stack(top_panel_name, bottom_panel_name)
                        Kp_ref = self.data_Kser_panels.loc[bottom_panel_name]["Kser,p,ref"]
                        coef_equi = kser_equi / Kp_ref
                        k_wall[1] = k_wall[1] + bottom_panel["Nombre"] * kser_equi
                        k_wall[2] = k_wall[2] + bottom_panel["Nombre"] * coef_equi
                        k_wall[3] = Kp_ref

                        already_reported.append(bottom_panel_name)
                        already_reported.append(top_panel_name)

                    # on gère ensuite les panneaux simples
                    for panel_name, panel in wall["Panneaux"].items():
                        if not panel["Contreventement"] or panel_name in already_reported:
                            continue
                        k_wall[1] = k_wall[1] + panel["Nombre"] * self.data_Kser_panels.loc[panel_name]["Kser,p unit."]
                        k_wall[2] = (
                            k_wall[2] + panel["Nombre"] * self.data_Kser_panels.loc[panel_name]["Coeff. équivalence"]
                        )
                        k_wall[3] = self.data_Kser_panels.loc[panel_name]["Kser,p,ref"]
                        k_sys_wall[3] = self.data_Kser_panels.loc[panel_name]["Kser,p,ref"]
                    self.data_Kser_walls.loc[wall_name] = k_wall
                    k_sys_wall[1] = k_sys_wall[1] + k_wall[1]
                    k_sys_wall[2] = k_sys_wall[2] + k_wall[2]
                self.data_Kser_sys_walls.loc[sys_name] = k_sys_wall

        # Trier par index (nommés 'Système de mur' et 'Mur')
        self.data_Kser_sys_walls.sort_index(inplace=True)
        self.data_Kser_walls.sort_index(inplace=True)
        return self.data_Kser_sys_walls

    def _FivRd(
        self,
        b_i: si.mm,
        h: si.mm,
        s: si.mm,
        F_f_Rk: si.kN,
        coef_double_panel: float = 1,
        K_mod: float = 1,
        gamma_M: float = 1.3,
    ):
        """
        Calcul la valeur de calcul (design) de la capacité résistante d'un panneau au contreventement.

        Args:
            b_i (si.mm): Largeur du panneau en mm.
            h (si.mm): Hauteur du panneau en mm.
            s (si.mm): Entraxe des organes en mm.
            F_f_Rk (si.kN): Capacité résistante caractéristique de cisaillement de l'organe en kN.
            K_mod (float, optional): Coefficient de modification Kmod. Defaults to 1.
            gamma_M (float, optional): Coefficient de matériaux. Par défaut 1.3 pour les assemblages.
        """
        if b_i >= h / 2:
            c_i = 1
        else:
            c_i = b_i / (h / 2)
        F_f_Rd = F_f_Rk * K_mod / gamma_M
        F_iv_Rd = (1.2 * coef_double_panel * F_f_Rd * b_i * c_i) / s
        return F_iv_Rd, c_i, coef_double_panel

    def _FtEd(self, wall_name: str):
        """
        Calcul l'effort de soulèvement d'un mur au contreventement.

        Args:
            wall_name (str): Nom du mur.
        """
        sys_wall = self.data_sys_walls[wall_name.split("_")[0]]
        wall = self.data_walls[wall_name]
        Fv_i_Ed_wall = self.data_walls_loads.loc[wall_name]["Fv,Ed,i ELU"]
        return Fv_i_Ed_wall * sys_wall["Hauteur du système de mur"] / wall["Longueur"]

    def taux_walls(self) -> "pd.DataFrame":
        """Calcule les taux de travail ELS et ELU en contreventement pour tous les éléments MOB.

        Peuple trois tableaux d'attribut accessibles après l'appel :

        - taux_FvRd_panels : taux de cisaillement de chaque panneau (F_iv,Ed / F_iv,Rd en %).
        - taux_FvRd_walls : taux de cisaillement de chaque mur (F_w,Ed / F_w,Rd en %).
        - taux_FvRd_sys_walls : taux de cisaillement et de déplacement ELS de chaque système de mur.

        Note:
            Appelle get_Kser_wall en interne. Les efforts peuvent être issus de
            calculate_loads_and_deformations ou fournis directement dans add_sys_wall.

        Returns:
            pd.DataFrame: Tableau taux_FvRd_walls indexé par le nom du mur, avec les colonnes
                ["Étage", "Fw,v,Ed", "Fw,v,Rd", "Taux cisaillement (%)"].
        """
        self.taux_FvRd_panels = pd.DataFrame(
            columns=["Étage", "Fiv,Ed", "ci", "Coeff. panneau double", "Fiv,Rd", "Taux cisaillement (%)"]
        )
        self.taux_FvRd_panels.index.name = "Panneau"
        self.taux_FvRd_walls = pd.DataFrame(columns=["Étage", "Fw,v,Ed", "Fw,v,Rd", "Taux cisaillement (%)"])
        self.taux_FvRd_walls.index.name = "Mur"
        self.taux_FvRd_sys_walls = pd.DataFrame(
            columns=[
                "Étage",
                "Fsw,v,Ed",
                "Fsw,v,Rd",
                "Taux cisaillement (%)",
                "Déplacement",
                "Critère ELS",
                "Taux déplacement (%)",
            ]
        )
        self.taux_FvRd_sys_walls.index.name = "Système de mur"

        self.get_Kser_wall()

        for etage, sys in self.data.items():
            for sys_name, sys_wall in sys.items():
                K_sw = self.data_Kser_sys_walls.loc[sys_name]["Kser système de mur"]
                Fv_Rd_sys_walls = 0

                # Si des charges manuelles sont données
                direct_loads = sys_wall.get("Fv_Ed_ELU", None)
                if direct_loads:
                    F_sw_Ed = direct_loads
                    F_sw_Ed_ELS = sys_wall.get("Fv_Ed_ELS", 0 * si.kN)
                    delta_v_sw_Ed = abs(F_sw_Ed_ELS / K_sw)
                else:
                    F_sw_Ed = self.data_sys_walls_loads.loc[sys_name]["Fv,Ed,i ELU"]
                    delta_v_sw_Ed = abs(self.data_sys_walls_loads.loc[sys_name]["delta,v,Ed,i ELS"])

                for wall_name, wall in sys_wall["Murs"].items():
                    K_w = self.data_Kser_walls.loc[wall_name]["Kser mur"]
                    Fv_Rd_walls = 0

                    if direct_loads:
                        F_w_Ed = direct_loads * K_w / K_sw
                        if not hasattr(self, "data_walls_loads"):
                            self.data_walls_loads = pd.DataFrame(
                                index=self.data_walls.keys(),
                                columns=["Étage", "Fv,Ed,i ELU", "Ft,Ed,i ELU (soulèvement)"],
                            )
                        self.data_walls_loads.loc[wall_name, "Fv,Ed,i ELU"] = F_w_Ed
                        self.data_walls_loads.loc[wall_name, "Ft,Ed,i ELU (soulèvement)"] = self._FtEd(wall_name)
                    else:
                        F_w_Ed = self.data_walls_loads.loc[wall_name]["Fv,Ed,i ELU"]

                    for panel_name, panel in wall["Panneaux"].items():
                        if not panel["Contreventement"]:
                            continue

                        # Si des charges manuelles sont données
                        if direct_loads:
                            F_f_Ed = (
                                F_w_Ed
                                * self.data_Kser_panels.loc[panel_name]["Kser,p unit."]
                                / self.data_Kser_walls.loc[wall_name]["Kser mur"]
                            )
                        else:
                            F_f_Ed = self.data_panels_loads.loc[panel_name]["Fv,Ed,i ELU"]

                        # Connecteur avec distance de couturage et capacité résistante caractéristique de cisaillement de l'organe
                        con = self._connecteurs[panel["Position"]]
                        if panel["Position"] == "Extérieur":
                            couturage_rive = wall["Couturage ext. rive"]
                            F_f_Rk = wall["Ff,Rk ext"]
                            K_mod = sqrt(con.beam_1._get_k_mod("Instantanee") * con.beam_2._get_k_mod("Instantanee"))
                        else:
                            couturage_rive = wall["Couturage int. rive"]
                            F_f_Rk = wall["Ff,Rk int"]
                            K_mod = sqrt(con.beam_1._get_k_mod("Instantanee") * con.beam_2._get_k_mod("Instantanee"))

                        # Calcul de la valeur de calcul (design) de la capacité résistante d'un panneau au contreventement
                        F_f_Rd, c_i, coef_double_panel = self._FivRd(
                            panel["Largeur"],
                            panel["Hauteur"],
                            couturage_rive,
                            F_f_Rk,
                            self.COEFF_DOUBLE_PANEL[panel["Panneau double"]],
                            K_mod,
                            gamma_M=1.3,
                        )

                        # Taux de travail panneau
                        self.taux_FvRd_panels.loc[panel_name] = [
                            etage,
                            F_f_Ed,
                            c_i,
                            coef_double_panel,
                            F_f_Rd,
                            round(abs(F_f_Ed * 100 / F_f_Rd), 2),
                        ]
                        Fv_Rd_walls = Fv_Rd_walls + F_f_Rd * panel["Nombre"]

                    # Taux de travail mur
                    self.taux_FvRd_walls.loc[wall_name] = [
                        etage,
                        F_w_Ed,
                        Fv_Rd_walls,
                        round(abs(F_w_Ed * 100 / Fv_Rd_walls), 2),
                    ]
                    Fv_Rd_sys_walls = Fv_Rd_sys_walls + Fv_Rd_walls

                # critère ELS
                crit_ELS = round(sys_wall["Hauteur étage"] / self.CRITERE_ELS, 2)
                # Taux de travail système de mur
                self.taux_FvRd_sys_walls.loc[sys_name] = [
                    etage,
                    F_sw_Ed,
                    Fv_Rd_sys_walls,
                    round(abs(F_sw_Ed * 100 / Fv_Rd_sys_walls), 2),
                    delta_v_sw_Ed,
                    crit_ELS,
                    round(abs(delta_v_sw_Ed * 100 / crit_ELS), 2),
                ]

        return self.taux_FvRd_walls

    def taux_ancrage(
        self,
        montant: object,
        lisse_basse: object,
        lisse_impl: object,
        Ft_Rk_wood: si.kN = 0,
        Ft_Rd_concrete: si.kN = 0,
        Fv_Rk_anc_lb: si.kN = 0,
        e_anc_lb: si.mm = 0,
        Fv_Rk_anc_li: si.kN = 0,
        e_anc_li: si.mm = 0,
        etage: str = None,
    ) -> "pd.DataFrame":
        """Calcule les taux de travail des ancrages pour tous les murs vérifiés.

        Effectue trois vérifications pour chaque mur :

        1. Équerre de soulèvement (traction) : F_t,Ed / min(F_t,Rd,bois, F_t,Rd,béton).
        2. Ancrage lisse basse sur lisse d'implantation : F_v,Ed / F_v,anc,lb,Rd.
        3. Ancrage lisse d'implantation au sol : F_v,Ed / F_v,anc,li,Rd.

        Note:
            taux_walls doit être appelé avant taux_ancrage pour peupler data_walls_loads.

        Args:
            montant (Barre): Montant de l'ossature bois portant l'équerre de traction,
                issu de Barre ou dérivé (EC5_Element_droit). Fournit le Kmod.
            lisse_basse (Barre): Lisse basse, issue de Barre ou dérivé. Fournit le Kmod.
            lisse_impl (Barre): Lisse d'implantation, issue de Barre ou dérivé. Fournit le Kmod.
            Ft_Rk_wood (float): Capacité résistante caractéristique de traction de l'équerre
                dans le montant en kN (nef déjà pris en compte). Defaults to 0.
            Ft_Rd_concrete (float): Capacité résistante de calcul en traction de l'équerre
                dans le béton en kN. Defaults to 0.
            Fv_Rk_anc_lb (float): Capacité résistante caractéristique de cisaillement de l'ancrage
                unitaire de la lisse basse en kN (nef déjà pris en compte). Defaults to 0.
            e_anc_lb (float): Entraxe de l'ancrage unitaire de la lisse basse en mm. Defaults to 0.
            Fv_Rk_anc_li (float): Capacité résistante caractéristique de cisaillement de l'ancrage
                unitaire de la lisse d'implantation en kN (nef déjà pris en compte). Defaults to 0.
            e_anc_li (float): Entraxe de l'ancrage unitaire de la lisse d'implantation en mm. Defaults to 0.
            etage (str, optional): Nom de l'étage à vérifier. Si None, tous les étages
                sont vérifiés avec ces mêmes valeurs. Defaults to None.

        Returns:
            pd.DataFrame: Tableau taux_ancr_walls indexé par le nom du mur, avec les colonnes
                ["Étage", "Longueur mur", "Ft,Ed", "Ft,Rd", "Taux ancr. soulèvement (%)",
                "Fv,anc,Ed", "Fv,anc,lb,Rd", "Nbr ancr,lb", "entraxe anc,lb",
                "Taux ancr. lisse basse (%)", "Fv,anc,li,Rd", "Nbr ancr,li",
                "entraxe anc,li", "Taux ancr. lisse impl. (%)"].
        """
        Ft_Rk_wood = Ft_Rk_wood * si.kN
        Ft_Rd_concrete = Ft_Rd_concrete * si.kN
        Fv_Rk_anc_lb = Fv_Rk_anc_lb * si.kN
        Fv_Rk_anc_li = Fv_Rk_anc_li * si.kN
        e_anc_lb = e_anc_lb * si.mm
        e_anc_li = e_anc_li * si.mm
        self.taux_ancr_walls = pd.DataFrame(
            columns=[
                "Étage",
                "Longueur mur",
                "Ft,Ed",
                "Ft,Rd",
                "Taux ancr. soulèvement (%)",
                "Fv,anc,Ed",
                "Fv,anc,lb,Rd",
                "Nbr ancr,lb",
                "entraxe anc,lb",
                "Taux ancr. lisse basse (%)",
                "Fv,anc,li,Rd",
                "Nbr ancr,li",
                "entraxe anc,li",
                "Taux ancr. lisse impl. (%)",
            ]
        )

        for wall_name, wall in self.data_walls_loads.iterrows():
            if etage and wall["Étage"] != etage:
                continue
            Ft_Ed = abs(wall["Ft,Ed,i ELU (soulèvement)"])
            Fv_Ed = abs(wall["Fv,Ed,i ELU"])
            K_mod_montant = montant._get_k_mod("Instantanee")
            K_mod_lb = lisse_basse._get_k_mod("Instantanee")
            K_mod_li = lisse_impl._get_k_mod("Instantanee")
            long_wall = self.data_walls[wall_name]["Longueur"]
            nb_anc_lb = floor(long_wall / e_anc_lb)
            nb_anc_li = floor(long_wall / e_anc_li)
            if nb_anc_lb == 0:
                raise ValueError(
                    f"Mur '{wall_name}' (L={long_wall}) : e_anc_lb={e_anc_lb} est supérieur "
                    f"à la longueur du mur — aucun ancrage lisse basse possible (nb_anc_lb=0)."
                )
            if nb_anc_li == 0:
                raise ValueError(
                    f"Mur '{wall_name}' (L={long_wall}) : e_anc_li={e_anc_li} est supérieur "
                    f"à la longueur du mur — aucun ancrage lisse d'implantation possible (nb_anc_li=0)."
                )
            Ft_Rd_wood = Ft_Rk_wood * K_mod_montant / Assemblage.GAMMA_M_ASS
            Fv_Rd_anc_lb = (Fv_Rk_anc_lb * K_mod_lb / Assemblage.GAMMA_M_ASS) * nb_anc_lb
            Fv_Rd_anc_li = (Fv_Rk_anc_li * K_mod_li / Assemblage.GAMMA_M_ASS) * nb_anc_li

            taux_soul = round(Ft_Ed * 100 / min(Ft_Rd_wood, Ft_Rd_concrete), 2)
            taux_lb = round(Fv_Ed * 100 / Fv_Rd_anc_lb, 2)
            taux_li = round(Fv_Ed * 100 / Fv_Rd_anc_li, 2)
            self.taux_ancr_walls.loc[wall_name] = [
                wall["Étage"],
                long_wall,
                Ft_Ed,
                min(Ft_Rd_wood, Ft_Rd_concrete),
                taux_soul,
                Fv_Ed,
                Fv_Rd_anc_lb,
                nb_anc_lb,
                e_anc_lb,
                taux_lb,
                Fv_Rd_anc_li,
                nb_anc_li,
                e_anc_li,
                taux_li,
            ]
        return self.taux_ancr_walls
