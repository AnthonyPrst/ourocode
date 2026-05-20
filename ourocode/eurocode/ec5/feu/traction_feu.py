# coding in UTF-8
# by Anthony PARISOT
import sys
import os
import warnings
import math as mt
from math import sqrt, pi
import pandas as pd

import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.ec5.feu.feu import Feu
from ourocode.eurocode.ec5.element_droit.traction import Traction



class Traction_feu(Feu, Traction):
    def __init__(self, **kwargs):
        """Classe permettant le calcul de la Traction d'un élément bois selon l'EN 1995.
        Cette classe est hérité de la classe Feu, provenant du module EC5_Feu.py.
        """
        super().__init__(**kwargs)
    
    def f_t_0_d(self) -> tuple:
        """Retourne la résistance f,t,0,d au feu de l'élément en MPa
        """
        return super()._f_type_d("ft0k")

    def taux_t_0_d(self) -> tuple:
        """Retourne le taux de travail en traction axial.

        Returns:
            float: taux de travail en %
        """
        self.taux_t_0_rd = {}
        K_h_y = self.K_h["y"]
        K_h_z = self.K_h["z"]
        sigma_t_0_d_fi = self.sigma_t_0_rd
        f_t_0_d_fi = self.f_type_rd

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            K_h = min(K_h_y, K_h_z)
            taux_6_1 = sigma_t_0_d_fi / (K_h * f_t_0_d_fi)  # equ6.1
            return taux_6_1

        value = val()
        self.taux_t_0_rd["equ6.1"] = value[1]
        synthese = [
            ["Traction bois", None, self.taux_t_0_rd['equ6.1']],
        ]
        self._add_synthese_taux_travail(synthese)
        return value
