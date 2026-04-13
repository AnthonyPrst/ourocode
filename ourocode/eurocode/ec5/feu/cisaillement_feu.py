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
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.feu.feu import Feu
from ourocode.eurocode.ec5.element_droit.cisaillement import Cisaillement


class Cisaillement_feu(Cisaillement, Feu):
    def __init__(self, **kwargs):
        """Classe qui permet de calculer le cisaillement d'une poutre comme décrit à l'EN 1995 §6.1.7 et §6.5.
        Cette classe est hérité de la classe Feu, provenant du module EC5_Feu.py.
        """
        super().__init__(**kwargs)

    def Kv(self, hef: si.mm, x: si.mm, i_lo: si.mm, ent=("Dessous", "Dessus")):
        """Retourne le facteur d'entaille Kv pour une entaille au niveau d'un appuis

        Args:
            hef (int): Hauteur efficace à froid de la poutre (hauteur - hauteur de l'entaille) en mm
            x (int):Distance entre le centre de réaction à l'appuis et le coin de l'entaille en mm
            i_lo (float): longueur horizontal de l'entaille en mm
            ent (tuple, optional): Entaille sur le dessus ou dessous de la poutre.

        Returns:
            float: facteur Kv
        """
        h_ef = hef*si.mm - self._def["Bas"][1] - self._def["Haut"][1]  # Réduction de la section
        return super().Kv(h_ef.value*10**3, x, i_lo, ent)
    
    def f_v_d(self):
        """Retourne la résistance f,v,d au feu de l'élément en MPa
        """
        return super()._f_type_d("fvk")

    def taux_tau_d(self):
        """Retourne le taux de travail en cisaillement en %"""
        self.taux_tau_rd = {}
        tau_d_fi = self.tau_rd
        f_v_d_fi = self.f_type_rd
        K_v = self.K_v

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_6_13 = tau_d_fi / f_v_d_fi
            taux_6_60 = tau_d_fi / (K_v * f_v_d_fi)
            return taux_6_13, taux_6_60

        value = val()
        self.taux_tau_rd["equ6.13"] = value[1][0]
        self.taux_tau_rd["equ6.60"] = value[1][1]
        return value
