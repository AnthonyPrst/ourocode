# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.assemblage.pointe import Pointe

# 8.4 Assemblage par Agrafe

class Agrafe(Pointe):
    TYPE_ASSEMBLAGE = ("Bois/Bois", ("CP", "Panneau dur", "PP/OSB"), "Bois/Métal")
    QUALITE_ACIER = ('8.8', '9.8', '10.9', '12.9')
    def __init__(self, d:si.mm, b_agrafe:si.mm, l:si.mm, qualite: str=QUALITE_ACIER, n: int=1, angle_sup_30: bool=["True", "False"], alpha1: float=0, alpha2: float=0, **kwargs):
        """Initialise une agrafe selon EN 1995-1-1 §8.4.

        Hérite de Pointe. Vérifie automatiquement les dimensions minimales (dos ≥ 6d,
        pénétration t2 ≥ 14d).

        Args:
            d (float): Diamètre de l'agrafe en mm. Pour section rectangulaire : racine carrée
                du produit des deux dimensions (EN 1995-1-1 §8.4(2)).
            b_agrafe (float): Dimension du dos de l'agrafe en mm (doit être ≥ 6d).
            l (float): Longueur sous tête en mm.
            qualite (str): Qualité de l'acier selon EN ISO 898-2. Defaults to "8.8".
            n (int): Nombre d'agrafes dans une file. Defaults to 1.
            angle_sup_30 (bool): True si l'angle entre la tête de l'agrafe et le fil du bois
                est supérieur à 30° (entraîne un facteur multiplicateur de 2 sur F_v,Rk,ass),
                False sinon (facteur 0.7). Defaults to True.
            alpha1 (float, optional): Angle entre l'effort et le fil du bois 1 en °. Defaults to 0.
            alpha2 (float, optional): Angle entre l'effort et le fil du bois 2 en °. Defaults to 0.
            **kwargs: Arguments transmis à la classe parent Pointe / Assemblage.
        """
        super().__init__(d=d, dh=0, l=l, qualite=qualite, n=n, alpha1=alpha1, alpha2=alpha2, type_organe="Agrafe", percage=False, **kwargs)
        self.b_agrafe = b_agrafe * si.mm
        self.angle_sup_30 = angle_sup_30
        self._dimension_min()

    def _dimension_min(self):
        """ Vérifie si les dimensions minimales sont respectées """
        if self.b_agrafe < self.d*6:
            raise ValueError(f"Erreur, la dimension du dos de l'agrafe est inférieur à 6 fois le diamètre de l'agrafe, le dos de l'agrafe minimal est de {self.d*6}")
        if self.t2 < self.d*14:
            raise ValueError(f"Erreur, la longueur de pénétration t2 est inférieur à 14 fois le diamètre de l'agrafe. La longueur d'agrafe minimal est de {self.t1 + self.d*14}")
        
    @property
    def MyRk(self) -> tuple:
        """Calcule le moment d'écoulement plastique d'une agrafe M_y,Rk selon EN 1995-1-1 §8.4(3).

        Formule : M_y,Rk = 150 × d³ en N.mm.

        Returns:
            tuple: (latex_string, M_y_Rk) où M_y_Rk est le moment plastique en N.mm (avec unité si.N×si.mm).
        """
        d = self.d.value * 10**3
        @handcalc(override="short", precision=2, left="\\[", right="\\]")
        def val():
            M_y_Rk = 150 * d**3 # N.mm
            return M_y_Rk * si.N*si.mm
        return val()

    def nef(self) -> int:
        """Retourne le nombre efficace d'agrafes n_ef dans une file.

        Pour les agrafes, n_ef = n (pas de réduction par effet de groupe, EN 1995-1-1 §8.4).

        Returns:
            int: Nombre efficace d'agrafes dans une file (= n).
        """
        self._nef = self.n
        return self._nef

    @property
    def pince(self) -> dict:
        """Retourne les pinces minimales pour une agrafe selon EN 1995-1-1 §8.4.

        Les pinces dépendent de l'angle alpha entre l'effort et le fil du bois,
        et de l'option angle_sup_30.

        Returns:
            dict: Dictionnaire par barre ("barre 1", "barre 2") contenant les pinces
                a1, a2, a3t, a3c, a4t, a4c en mm (avec unité si.mm).
        """
        dict_pince = {}
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if self._type_beam[i] not in self.TYPE_BOIS_ASSEMBLAGE:
                continue
            alpha = self.alpha[i]
            if not self.angle_sup_30:
                a1 = round((10 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
            else:
                a1 = round((15 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
            a2 = round(15 * self.d, 1)
            a3t = round((15 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
            a3c = round(15 * self.d, 1)
            a4t = round((15 + 5 * mt.sin(mt.radians(alpha))) * self.d, 1)
            a4c = round(10 * self.d, 1)
            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince
