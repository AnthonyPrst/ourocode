# coding in UTF-8
# by Anthony PARISOT
import math as mt
from math import sqrt

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec3.element_droit.plat import Plat


class Cisaillement(Plat):
    def __init__(self, Av: si.mm**2, *args, **kwargs):
        """Defini une classe permettant le calcul d'un élément métallique en cisaillement selon l'EN 1993-1-1 §6.2.6.
        Cette classe est hérité de la classe Plat du module E3_Element_droit.py.

        Args:
            Av (float): Aire de cisaillemment en mm²
        """
        super().__init__(*args, **kwargs)
        self.Av = Av * si.mm**2

    @property
    def Vpl_Rd(self):
        """Calcul la résistance du cisaillement plastique en N (équa 6.18)
        """
        A_v = self.Av
        f_y = self.fy
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            Vpl_Rd = (A_v * f_y/sqrt(3)) / gamma_M0 #(équa 6.18)
            return Vpl_Rd
        return val() 

