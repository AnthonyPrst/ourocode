# coding in UTF-8
# by Anthony PARISOT
import difflib
import math as mt
from math import sqrt
import pandas as pd

import forallpeople as si

si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.core.projet import Projet


class Plat(Projet):
    GAMMA_M = {
        "gamma_M0": 1,
        "gamma_M1": 1,
        "gamma_M2": 1.25,
        "gamma_M3": 1.1,
        "gamma_M3_ser": 1.25,
        "gamma_M4": 1,
        "gamma_M5": 1,
        "gamma_M6_ser": 1,
        "gamma_M7": 1.1,
    }
    E = 210000 * si.MPa
    CLASSE_STEEL = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_acier.csv").index)
    _PROFIL_DATA = Projet._data_from_csv(Projet, "profil_acier.csv", index_col=None)
    TYPES_PROFIL = tuple(_PROFIL_DATA["type"].unique())
    DESIGNATIONS = tuple(_PROFIL_DATA["designation"])

    def __init__(
        self,
        t: si.mm = 0,
        h: si.mm = 0,
        b: si.mm = 0,
        classe_acier: str = CLASSE_STEEL,
        classe_transv: int = ("1", "2", "3", "4"),
        designation: str = None,
        **kwargs,
    ):
        """Initialise un élément acier selon l'EN 1993-1-1.

        Définit la géométrie et les caractéristiques mécaniques (fy, fu) en fonction
        de la classe d'acier et de l'épaisseur. Classe de base dont héritent Profil, Traction,
        Compression, Cisaillement et Flexion.

        Sans `designation`, la section est rectangulaire t × b (Iy = t·b³/12). Avec une
        `designation` de profilé (ex : "IPE 300"), la géométrie est chargée depuis
        `data/profil_acier.csv` et t, h, b sont ignorés (t = tw). Les caractéristiques de
        section (A, Avz, Avy, Iy, Iz, Wely, Welz, Wply, Wplz, iy, iz) sont ainsi disponibles
        automatiquement pour les classes de vérification.

        Args:
            t (int, optional): Épaisseur du plat en mm. Defaults to 0.
            h (int, optional): Hauteur du plat en mm. Defaults to 0.
            b (int, optional): Largeur du plat en mm. C'est cette dimension qui, avec t,
                détermine le moment quadratique pour les vérifications. Defaults to 0.
            classe_acier (str, optional): Classe d'acier selon EN 1993-1-1 Tableau 3.1
                (ex. "S235", "S275", "S355"). Defaults to "S235".
            classe_transv (int, optional): Classe transversale de la section (1, 2 ou 3)
                selon EN 1993-1-1 §5.5. La classe 4 n'est pas implémentée. Defaults to 1.
            designation (str, optional): Désignation d'un profilé du catalogue (voir Profil).
                Defaults to None (section rectangulaire).
            **kwargs: Arguments transmis à la classe parent Projet.

        Raises:
            ValueError: Si classe_acier n'est pas dans CLASSE_STEEL.
            ValueError: Si classe_transv n'est pas 1, 2 ou 3.
            ValueError: Si designation n'est pas dans DESIGNATIONS.
        """
        super().__init__(**kwargs)
        self._check_dimension(t, "t", strict=False)
        self._check_dimension(h, "h", strict=False)
        self._check_dimension(b, "b", strict=False)
        classe_acier = self._check_choice(classe_acier, self.CLASSE_STEEL, "classe_acier")
        classe_transv = self._check_choice(
            int(classe_transv) if not isinstance(classe_transv, (list, tuple)) else int(classe_transv[0]),
            (1, 2, 3),
            "classe_transv",
        )
        self.classe_acier = classe_acier
        self.classe_transv = int(classe_transv)
        self.designation = designation
        if designation:
            self._load_profil(designation)
        else:
            self.t = t * si.mm
            self.h = h * si.mm
            self.b = b * si.mm
            self._geometrie_rectangulaire()
        self.__fy_fu()

    def _load_profil(self, designation: str):
        """Charge la géométrie d'un profilé laminé depuis `profil_acier.csv` (valeurs en mm)."""
        if designation not in self.DESIGNATIONS:
            proches = difflib.get_close_matches(designation, self.DESIGNATIONS, n=5)
            raise ValueError(
                f"Désignation '{designation}' invalide. Désignations proches : {proches}. "
                f"Types disponibles : {self.TYPES_PROFIL}"
            )
        df = self._PROFIL_DATA
        p = df[df["designation"] == designation].iloc[0]

        def val(name, unit):
            return float(p[name]) * unit if pd.notna(p[name]) else None

        self.type_profil = p["type"]
        self.h = val("h", si.mm)
        self.b = val("b", si.mm)
        self.tw = val("tw", si.mm)
        self.tf = val("tf", si.mm)
        self.r = val("r", si.mm)
        self.t = self.tw
        self.A = val("A", si.mm**2)
        self.Avz = val("Avz", si.mm**2)
        self.Avy = self.A - (self.h - 2 * self.tf) * self.tw
        self.Iy = val("Iy", si.mm**4)
        self.Iz = val("Iz", si.mm**4)
        self.Wely = val("Wely", si.mm**3)
        self.Welz = val("Welz", si.mm**3)
        self.Wply = val("Wply", si.mm**3)
        self.Wplz = val("Wplz", si.mm**3)
        self.iy = val("iy", si.mm)
        self.iz = val("iz", si.mm)
        self.masse = val("masse", si.kg / si.m)

    def _geometrie_rectangulaire(self):
        """Caractéristiques d'une section rectangulaire t × b (b selon z, flexion forte autour de y)."""
        self.type_profil = "Plat"
        if not (self.t.value and self.b.value):
            self.A = self.Avz = self.Avy = self.Iy = self.Iz = None
            self.Wely = self.Welz = self.Wply = self.Wplz = self.iy = self.iz = None
            return
        t, b = self.t, self.b
        self.A = t * b
        self.Avz = self.Avy = self.A
        self.Iy = t * b**3 / 12
        self.Iz = b * t**3 / 12
        self.Wely = t * b**2 / 6
        self.Welz = b * t**2 / 6
        self.Wply = t * b**2 / 4
        self.Wplz = b * t**2 / 4
        self.iy = b / sqrt(12)
        self.iz = t / sqrt(12)

    @property
    def _t_ref(self):
        """Épaisseur de référence pour fy/fu : max(tw, tf) pour un profilé, t sinon."""
        return max(self.tw, self.tf) if self.designation else self.t

    @property
    def __classe_acier(self):
        """Retourne le dataframe de la classe d'acier définit"""
        df = self._data_from_csv("caracteristique_meca_acier.csv")
        df = df.loc[self.classe_acier]
        return df

    def __fy_fu(self):
        """Défini fy (résistance élastique) et fu (résistance plastique) en MPa fonction de la classe d'acier choisi"""
        t = self._t_ref
        if t <= 40 * si.mm:
            self.fy = self.__classe_acier.loc["t<= 40  fy"] * si.MPa
            self.fu = self.__classe_acier.loc["t<= 40  fu"] * si.MPa
        elif t <= 80 * si.mm:
            self.fy = self.__classe_acier.loc["40<t<= 80  fy"] * si.MPa
            self.fu = self.__classe_acier.loc["40<t<= 80  fu"] * si.MPa
        else:
            raise ValueError(f"Épaisseur {t} > 80 mm : non couverte par le tableau 3.1 de l'EN 1993-1-1.")

    @property
    def _inertie(self):
        """Retourne les moments quadratiques [I_y, I_z] de la section en mm⁴ (avec unité si.mm⁴).

        Section rectangulaire : I_y = t × b³ / 12 et I_z = b × t³ / 12 ; profilé : valeurs du catalogue.
        """
        return [self.Iy, self.Iz]
