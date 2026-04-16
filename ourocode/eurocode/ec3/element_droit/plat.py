# coding in UTF-8
# by Anthony PARISOT
import math as mt
from math import sqrt
import pandas as pd

import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core.renderer import handcalc

from ourocode.eurocode.core.projet import Projet

class Plat(Projet):
    GAMMA_M = {"gamma_M0": 1, "gamma_M1": 1, "gamma_M2": 1.25, "gamma_M3": 1.1, "gamma_M3_ser": 1.25,
             "gamma_M4": 1, "gamma_M5": 1, "gamma_M6_ser": 1, "gamma_M7": 1.1}
    E = 210000 * si.MPa
    CLASSE_STEEL = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_acier.csv").index)
             
    def __init__(self, t: si.mm=0, h: si.mm=0, b: si.mm=0, classe_acier: str=CLASSE_STEEL, classe_transv: int=("1","2","3","4"), **kwargs):
        """Initialise un élément acier de section rectangulaire selon l'EN 1993-1-1.

        Définit la géométrie (t, b, h) et les caractéristiques mécaniques (fy, fu) en fonction
        de la classe d'acier et de l'épaisseur. Classe de base dont héritent Traction, Compression,
        Cisaillement et Flexion.

        Args:
            t (int, optional): Épaisseur du plat en mm. Defaults to 0.
            h (int, optional): Hauteur du plat en mm. Defaults to 0.
            b (int, optional): Largeur du plat en mm. C'est cette dimension qui, avec t,
                détermine le moment quadratique pour les vérifications. Defaults to 0.
            classe_acier (str, optional): Classe d'acier selon EN 1993-1-1 Tableau 3.1
                (ex. "S235", "S275", "S355"). Defaults to "S235".
            classe_transv (int, optional): Classe transversale de la section (1, 2 ou 3)
                selon EN 1993-1-1 §5.5. La classe 4 n'est pas implémentée. Defaults to 1.
            **kwargs: Arguments transmis à la classe parent Projet.

        Raises:
            ValueError: Si classe_acier n'est pas dans CLASSE_STEEL.
            ValueError: Si classe_transv n'est pas 1, 2 ou 3.
        """
        super().__init__(**kwargs)
        if classe_acier not in self.CLASSE_STEEL:
            raise ValueError(
                f"Classe d'acier '{classe_acier}' invalide. Valeurs acceptées : {self.CLASSE_STEEL}"
            )
        if int(classe_transv) not in (1, 2, 3):
            raise ValueError(
                f"Classe transversale '{classe_transv}' invalide. Valeurs acceptées : 1, 2, 3 "
                "(la classe 4 n'est pas développée)"
            )
        self.t = t * si.mm
        self.h = h * si.mm
        self.b = b * si.mm
        self.classe_acier = classe_acier
        self.classe_transv = int(classe_transv)
        self.__fy_fu()
    

    @property
    def __classe_acier(self):
        """Retourne le dataframe de la classe d'acier définit 
        """
        df = self._data_from_csv("caracteristique_meca_acier.csv")
        df = df.loc[self.classe_acier]
        return df


    def __fy_fu(self):
        """Défini fy (résistance élastique) et fu (résistance plastique) en MPa fonction de la classe d'acier choisi
        """
        if self.t <= 40*si.mm:
            self.fy = self.__classe_acier.loc["t<= 40  fy"] * si.MPa
            self.fu = self.__classe_acier.loc["t<= 40  fu"] * si.MPa
        elif self.t > 40*si.mm and self.t<= 80*si.mm:
            self.fy = self.__classe_acier.loc["40<t<= 80  fy"] * si.MPa
            self.fu = self.__classe_acier.loc["40<t<= 80  fu"] * si.MPa

    
    @property
    def _inertie(self):
        """Retourne les moments quadratiques [I_y, I_z] d'une section rectangulaire en mm⁴.

        Calcule I_y = t × b³ / 12 et I_z = b × t³ / 12 si b et t sont renseignés.
        Retourne les valeurs pré-calculées sinon.

        Returns:
            list: [I_y, I_z] en mm⁴ (avec unité si.mm⁴).
        """
        if self.t and self.b:
            self.Iy = (self.t * self.b**3)/12
            self.Iz = (self.b * self.t**3)/12
            return [self.Iy, self.Iz]

        elif self.Iy and self.Iz:
            return [self.Iy, self.Iz]


