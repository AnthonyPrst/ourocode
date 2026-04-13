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
        """
        Créer une classe Agrafe hérité de la classe Assemblage du module EC5_Assemblage.py.
        
        Args:
            d (float): diamètre de l'agrafe en mm, si l'agrafe est de section rectangulaire alors c'est la racine carrée du produit des 2 dimensions selon EN1995 §8.4(2).
            b_agrafe (float): dimension du dos de l'agrafe en mm.
            l (int): longueur sous la tête en mm
            qualite (str): qualité de l'acier
            n (int): nombre d'agrafe dans une file 
            angle_sup_30 (str): Si l'angle entre la tête de l'agrafe et le fil du bois est supérieur à 30° alors True sinon False. Defaults to True.
            alpha1 (float, optional): angle entre l'effort de l'agrafe et le fil du bois 1 en °. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'agrafe et le fil du bois 2 en °. Defaults to 0.
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
        """ Défini le moment d'écoulement plastique d'une pointe en N.mm avec"""
        d = self.d.value * 10**3
        @handcalc(override="short", precision=2, left="\\[", right="\\]")
        def val():
            M_y_Rk = 150 * d**3 # N.mm
            return M_y_Rk * si.N*si.mm
        return val()

    def nef(self) -> int:
        """Retourne le nombre efficace d'organe dans une file"""
        self._nef = self.n
        return self._nef

    @property
    def pince(self) -> dict:
        """
        Défini les différentes pinces minimales pour une pointe en mm.

        Args:
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace de la pointe ou du tire fond si d<=6mm en mm
        """
        dict_pince = {}
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if not self._type_beam[i] in self.TYPE_BOIS_ASSEMBLAGE:
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
