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
from ourocode.eurocode.ec5.element_droit.compression import Compression



class Compression_feu(Feu, Compression):
    def __init__(self, lo_y: si.mm, lo_z: si.mm, type_appuis: str = Compression.COEF_LF, **kwargs):
        """Classe permettant le calcul de la Compression d'un élément bois selon l'EN 1995.
        Cette classe est hérité de la classe Feu, provenant du module EC5_Feu.py.

        Args:
            lo : Longueur de flambement suivant l'axe de rotation (y ou z) en mm si pas de flambement alors 0
            type_appuis : Coefficient multiplicateur de la longueur pour obtenir la longeur efficace de flambement en
                        fonction des du type d'appuis :
                                                        Encastré 1 côté : 2
                                                        Rotule - Rotule : 1
                                                        Encastré - Rotule : 0.7
                                                        Encastré - Encastré : 0.5
                                                        Encastré - Rouleau : 1
        """
        super().__init__(lo_y=lo_y, lo_z=lo_z, type_appuis=type_appuis, **kwargs)

    @property
    def lamb_rel_Axe(self) -> tuple:
        """Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation"""
        lamb_y = self.lamb[1]["y"]
        lamb_z = self.lamb[1]["z"]
        f_c0k = float(self.caract_meca.loc["fc0k"]) * si.MPa
        E_0_05 = int(self.caract_meca.loc["E005"]) * si.MPa
        K_fi = self.K_fi

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            lamb_rel_y_fi = (lamb_y / pi) * sqrt(f_c0k * K_fi / (E_0_05 * K_fi))
            lamb_rel_z_fi = (lamb_z / pi) * sqrt(f_c0k * K_fi / (E_0_05 * K_fi))
            return {"y": lamb_rel_y_fi, "z": lamb_rel_z_fi}

        return val()

    @property
    def k_Axe(self) -> tuple:
        """Retourne le facteur Ky ou Kz (fonction de l'axe de flambement)"""
        beta_C = self.beta_C
        lamb_rel_y_fi = self.lamb_rel_Axe[1]["y"]
        lamb_rel_z_fi = self.lamb_rel_Axe[1]["z"]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            k_y_fi = 0.5 * (1 + beta_C * (lamb_rel_y_fi - 0.3) + lamb_rel_y_fi**2)
            k_z_fi = 0.5 * (1 + beta_C * (lamb_rel_z_fi - 0.3) + lamb_rel_z_fi**2)
            return {"y": k_y_fi, "z": k_z_fi}

        return val()

    @property
    def kc_Axe(self) -> tuple:
        """Retourne le coefficient multiplicateur KcAxe  (axe = y ou z suivant axe de rotation en flambement) de fc,0,d"""
        k_y_fi = self.k_Axe[1]["y"]
        k_z_fi = self.k_Axe[1]["z"]
        lamb_rel_y_fi = self.lamb_rel_Axe[1]["y"]
        lamb_rel_z_fi = self.lamb_rel_Axe[1]["z"]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            k_c_y = 1 / (k_y_fi + sqrt(k_y_fi**2 - lamb_rel_y_fi**2))
            k_c_z = 1 / (k_z_fi + sqrt(k_z_fi**2 - lamb_rel_z_fi**2))
            return {"y": min(k_c_y, 1), "z": min(k_c_z, 1)}

        return val()
    
    def f_c_0_d(self) -> tuple:
        """Retourne la résistance f,c,0,d au feu de l'élément en MPa
        """
        return super()._f_type_d("fc0k")

    def taux_c_0_d(self, flexion: object = None) -> tuple:
        """Retourne les taux de travaux de la compression axial.
        Si l'élement est un poteau (donc avec un travail principalement en compression) et de la flexion combinée (EN 1995-1-1 §6.3.2),
        il est possible d'ajouter l'objet flexion et de vérifier cette combinaison.

        Args:
            flexion (object, optional): L'objet Flexion avec ces taux de travaux préalablement calculés. Default to None.

        Returns:
            list: retourne la liste des taux de travaux en %
        """
        self.taux_c_0_rd = {}
        sigma_c_0_d_fi = self.sigma_c_0_rd
        f_c_0_d_fi = self.f_type_rd
        lamb_rel_y = self.lamb_rel_Axe[1]["y"]
        lamb_rel_z = self.lamb_rel_Axe[1]["z"]
        K_c_y_fi = self.kc_Axe[1]["y"]
        K_c_z_fi = self.kc_Axe[1]["z"]

        from ourocode.eurocode.ec5.feu.flexion_feu import Flexion_feu as _Flexion_feu
        if flexion and isinstance(flexion, _Flexion_feu):
            taux_6_11 = flexion.taux_m_rd["equ6.11"]
            taux_6_12 = flexion.taux_m_rd["equ6.12"]
        else:
            taux_6_11 = 0
            taux_6_12 = 0

        if lamb_rel_y < 0.3 and lamb_rel_z < 0.3:

            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_6_2 = sigma_c_0_d_fi / f_c_0_d_fi  # equ6.2
                taux_6_19 = (sigma_c_0_d_fi / (f_c_0_d_fi * K_c_y_fi)) ** 2 + taux_6_11  # equ6.19
                taux_6_20 = (sigma_c_0_d_fi / (f_c_0_d_fi * K_c_z_fi)) ** 2 + taux_6_12  # equ6.20
                return taux_6_2, taux_6_19, taux_6_20

            value = val()
            self.taux_c_0_rd["equ6.19"] = value[1][1]
            self.taux_c_0_rd["equ6.20"] = value[1][2]
        else:

            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_6_2 = sigma_c_0_d_fi / f_c_0_d_fi  # equ6.2
                taux_6_23 = (sigma_c_0_d_fi / (f_c_0_d_fi * K_c_y_fi) + taux_6_11)  # equ6.23
                taux_6_24 = (sigma_c_0_d_fi / (f_c_0_d_fi * K_c_z_fi) + taux_6_12)  # equ6.24
                return taux_6_2, taux_6_23, taux_6_24

            value = val()
            self.taux_c_0_rd["equ6.23"] = value[1][1]
            self.taux_c_0_rd["equ6.24"] = value[1][2]

        self.taux_c_0_rd["equ6.2"] = value[1][0]
        return value
