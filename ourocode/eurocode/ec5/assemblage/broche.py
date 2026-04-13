# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.assemblage.boulon import Boulon

# 8.6 Assemblage par broche

class Broche(Boulon):
    def __init__(self, d:float, qualite: str=Boulon.QUALITE_ACIER, n: int=1, alpha1: float=0, alpha2: float=0, t1: int=0, t2: int=0, **kwargs):
        """
        Créer une classe Broche hérité de la classe Assemblage du module EC5_Assemblage.py.

        Args:
            d (int): diamètre efficace de la broche ( entre 6 et 30 mm) en mm
            qualite (str): qualite de la broche
            n (int): nombre de broche dans une file
            alpha1 (float, optional): angle entre l'effort de l'organe et le fil du bois en ° pour la barre 1. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'organe et le fil du bois en ° pour la barre 2. Defaults to 0.
            t1 (int, optional): longueur de contacte avec la tige  pour la pièce 1 en mm. 
                ATTENTION : Cet argument n'est pas obligatoire par défaut, il est calculer par le type de tige utilisée.
                Il n'est nécessaire de le remplir que si vous avez un t1 spécifique, par exemple avec une chapelle réduisant ainsi la portance local à une longueur inférieur à celle de l'épaisseur de la pièce 1.
            t2 (int, optional): longueur de contacte avec la tige  pour la pièce 2 en mm.
                ATTENTION : Même chose que pour t1 mais pour la pièce 2.
        """
        super().__init__(d, qualite, n, alpha1, alpha2, **kwargs)
        self.type_organe = "Broche"
        self.FaxRk = 0

    @property
    def Fax_Rk(self):
        return self.FaxRk
        
    @property
    def pince(self):
        """Défini les différentes pinces minimales pour une broche en mm avec :
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm """
        dict_pince = {}
        for i, alpha in enumerate(self.alpha):
            a1 = round((3 + 2 * (mt.cos(mt.radians(alpha)))) * self.d, 1)
            a2 = round(3 * self.d, 1)
            a3t = round(max(7 * self.d, 80*si.mm), 1)

            if alpha <= 150 and alpha < 210:

                a3c = round(max(3.5 * self.d, 40*si.mm), 1)

            else:

                a3c = round((a3t * mt.sin(mt.radians(alpha))), 1)

            a4t = round(max((2 + 2 * mt.sin(mt.radians(alpha))) * self.d, 3 * self.d), 1)
            a4c = round(3 * self.d, 1)
            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince

