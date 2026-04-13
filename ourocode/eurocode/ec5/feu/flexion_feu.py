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
from ourocode.eurocode.ec5.element_droit.flexion import Flexion



class Flexion_feu(Feu, Flexion):
    def __init__(
        self,
        lo_rel_y:si.mm, lo_rel_z:si.mm, 
        coeflef_y: float=0.9, 
        coeflef_z: float=0.9,
        pos: str = Flexion.LOAD_POS,
        **kwargs,
    ):
        """Classe permettant le calcul de la flexion d'une poutre bois selon l'EN 1995 §6.1.6, §6.2.3, §6.2.4 et §6.3.3.
        Cette classe est hérité de la classe Feu, provenant du module EC5_Feu.py.

        Args:
            lo_rel_y/z (int): longueur de déversemment autour de l'axe défini en mm
            coeflef_y/z (float): appuis simple :
                                            Moment constant : 1
                                            Charge répartie constante : 0.9
                                            Charge concentrée au milieu de la portée : 0.8
                                porte à faux :
                                            Charge répartie constante : 0.5
                                            Charge concentrée agissant à l'extrémité libre : 0.8.
            pos (str): positionnement de la charge sur la hauteur de poutre
        """
        super().__init__(lo_rel_y=lo_rel_y, lo_rel_z=lo_rel_z, coeflef_y=coeflef_y, coeflef_z=coeflef_z, pos=pos, **kwargs)

    @property
    def sigma_m_crit(self):
        """ Retourne sigma m,crit pour la prise en compte du déversement d'une poutre """
        self.l_ef_y = self.lo_rel_y * self.coeflef['y']
        self.l_ef_z = self.lo_rel_z * self.coeflef['z']
        if self.pos == "Charge sur fibre comprimée":
            self.l_ef_y = self.l_ef_y + 2 * self.h_calcul
            self.l_ef_z = self.l_ef_z + 2 * self.h_calcul
        elif self.pos == "Charge sur fibre tendue":
            self.l_ef_y = self.l_ef_y - 0.5 * self.h_calcul
            self.l_ef_z = self.l_ef_z - 0.5 * self.h_calcul
        
        self.l_ef = {"y": self.l_ef_y, "z": self.l_ef_z}
        b_calcul = self.b_calcul
        h_calcul = self.h_calcul
        l_ef_y = self.l_ef['y']
        l_ef_z = self.l_ef['z']
        E_0_05 = int(self.caract_meca.loc['E005']) * si.MPa
        K_fi = self.K_fi
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_m_crit_y_fi = (0.78 * b_calcul ** 2 * E_0_05 * K_fi) / (h_calcul * l_ef_y)
            sigma_m_crit_z_fi = (0.78 * h_calcul ** 2 * E_0_05 * K_fi) / (b_calcul * l_ef_z)
            return {"y": sigma_m_crit_y_fi, "z": sigma_m_crit_z_fi}
        return val()

    @property
    def lamb_rel_m(self):
        """ Retourne l'élancement relatif de la section avec pour argument """
        f_m0k = float(self.caract_meca.loc['fm0k']) *si.MPa
        sigma_m_crit_y_fi = self.sigma_m_crit[1]['y']
        sigma_m_crit_z_fi = self.sigma_m_crit[1]['z']
        K_fi = self.K_fi

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            lamb_rel_m_y_fi = sqrt(f_m0k * K_fi / sigma_m_crit_y_fi)
            lamb_rel_m_z_fi = sqrt(f_m0k * K_fi / sigma_m_crit_z_fi)
            return {"y": lamb_rel_m_y_fi, "z": lamb_rel_m_z_fi}
        return val()

    @property
    def K_crit(self):
        """ Retourne K,crit le coef. de minoration de la résistance en flexion au déversement"""
        lamb_rel_m_y = self.lamb_rel_m[1]['y']
        lamb_rel_m_z = self.lamb_rel_m[1]['z']
        result = [None, {"y": None, "z": None}]

        for axe in ["y", "z"]:
            lamb_rel_m_fi = lamb_rel_m_y if axe == "y" else lamb_rel_m_z
            if lamb_rel_m_fi <= 0.75:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_crit_fi = 1
                    axe
                    return K_crit_fi
            elif 0.75 < lamb_rel_m_fi <= 1.4:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_crit_fi = 1.56 - 0.75 * lamb_rel_m_fi
                    axe
                    return K_crit_fi
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_crit_fi = 1 / (lamb_rel_m_fi ** 2)
                    axe
                    return K_crit_fi
            kcrit_axe = val()
            if result[0]:
                result[0] = result[0] + kcrit_axe[0]
            else:
                result[0] = kcrit_axe[0]
            result[1][axe] = kcrit_axe[1]
        result = (result[0], result[1])
        return result
    
    def f_m_d(self):
        """Retourne la résistance f,m,d au feu de l'élément en MPa
        """
        return super()._f_type_d("fm0k")

    def taux_m_d(self, compression: object = None, traction: object = None):
        """Retourne les différents taux de travaux en flexion.
        Si l'élement est une poutre (donc avec un travail principalement en flexion) et de la compression (EN 1995-1-1 §6.3.3) ou de la traction (EN 1995-1-1 §6.2.3) combinée,
        il est possible d'ajouter l'objet Compression et Traction et de vérifier ces combinaisons.

        Args:
            compression (object, optional): L'objet Compression avec ces taux de travaux préalablement calculés. Defaults to None.
            traction (object, optional): L'objet Traction avec ces taux de travaux préalablement calculés. Defaults to None.

        Returns:
            list: retourne la liste des taux de travaux en %"""
        self.taux_m_rd = {}

        sigma_my_d_fi = self.sigma_m_rd["y"]
        sigma_mz_d_fi = self.sigma_m_rd["z"]
        f_m_d_fi = self.f_type_rd
        K_h_y = self.K_h["y"]
        K_h_z = self.K_h["z"]
        K_m = self.K_m
        K_crit_y_fi = self.K_crit[1]["y"]
        K_crit_z_fi = self.K_crit[1]["z"]

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def base():
            taux_6_11 = sigma_my_d_fi / (f_m_d_fi * K_h_y) + K_m * sigma_mz_d_fi / (f_m_d_fi * K_h_z)  # equ6.11
            taux_6_12 = K_m * sigma_my_d_fi / (f_m_d_fi * K_h_y) + sigma_mz_d_fi / (f_m_d_fi * K_h_z)  # equ6.12
            taux_6_33y = sigma_my_d_fi / (f_m_d_fi * K_h_y * K_crit_y_fi)  # equ6.33
            taux_6_33z = sigma_mz_d_fi / (f_m_d_fi * K_h_z * K_crit_z_fi)  # equ6.33
            return taux_6_11, taux_6_12, taux_6_33y, taux_6_33z

        base_val = base()
        latex = base_val[0]
        self.taux_m_rd["equ6.11"] = base_val[1][0]
        self.taux_m_rd["equ6.12"] = base_val[1][1]
        self.taux_m_rd["equ6.33y"] = base_val[1][2]
        self.taux_m_rd["equ6.33z"] = base_val[1][3]

        from ourocode.eurocode.ec5.feu.compression_feu import Compression_feu as _Compression_feu
        from ourocode.eurocode.ec5.feu.traction_feu import Traction_feu as _Traction_feu
        if compression and isinstance(compression, _Compression_feu):
            sigma_c_0_d = compression.sigma_c_0_rd
            f_c_0_d = compression.f_type_rd
            K_c_y = compression.kc_Axe[1]["y"]
            K_c_z = compression.kc_Axe[1]["z"]
            taux_6_2 = compression.taux_c_0_rd["equ6.2"]

            # taux_6_23 = compression.taux_c_0_rd['equ6.23']
            # taux_6_24 = compression.taux_c_0_rd['equ6.24']
            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def comp(taux_6_11, taux_6_12, taux_6_2, taux_6_33y, taux_6_33z):
                taux_6_19 = taux_6_2**2 + taux_6_11  # equ6.19
                taux_6_20 = taux_6_2**2 + taux_6_12  # equ6.20
                taux_6_23 = sigma_c_0_d / (f_c_0_d * K_c_y)  # equ6.23
                taux_6_24 = sigma_c_0_d / (f_c_0_d * K_c_z)  # equ6.24
                taux_6_35zyz = (taux_6_33y**2 + (sigma_mz_d_fi / (f_m_d_fi * K_h_z)) + taux_6_24)  # equ6.35
                taux_6_35yzz = (taux_6_33y + (sigma_mz_d_fi / (f_m_d_fi * K_h_z)) ** 2 + taux_6_24)  # equ6.35 interprétation
                taux_6_35yzy = (taux_6_33z**2 + (sigma_my_d_fi / (f_m_d_fi * K_h_y)) + taux_6_23)  # equ6.35
                taux_6_35zyy = (taux_6_33z + (sigma_my_d_fi / (f_m_d_fi * K_h_y)) ** 2 + taux_6_23)  # equ6.35 interprétation
                return (taux_6_19, taux_6_20, taux_6_35zyz, taux_6_35yzz, taux_6_35yzy, taux_6_35zyy)

            compression_val = comp(
                self.taux_m_rd["equ6.11"],
                self.taux_m_rd["equ6.12"],
                taux_6_2,
                self.taux_m_rd["equ6.33y"],
                self.taux_m_rd["equ6.33z"],
            )
            latex = latex + compression_val[0]
            self.taux_m_rd["equ6.19"] = compression_val[1][0]
            self.taux_m_rd["equ6.20"] = compression_val[1][1]
            self.taux_m_rd["equ6.35zyz"] = compression_val[1][2]  # 1er item axe de flexion pas au carré, 2eme item axe de flexion au carré, 3eme axe de compression
            self.taux_m_rd["equ6.35yzz"] = compression_val[1][3]
            self.taux_m_rd["equ6.35yzy"] = compression_val[1][4]
            self.taux_m_rd["equ6.35zyy"] = compression_val[1][5]

        if traction and isinstance(traction, _Traction_feu):
            taux_6_1 = traction.taux_t_0_rd["equ6.1"]

            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def tract(taux_6_11, taux_6_12):
                taux_6_17 = taux_6_11 + taux_6_1  # equ6.17
                taux_6_18 = taux_6_12 + taux_6_1  # equ6.18
                return taux_6_17, taux_6_18

            traction_val = tract(self.taux_m_rd["equ6.11"], self.taux_m_rd["equ6.12"])
            latex = latex + traction_val[0]
            self.taux_m_rd["equ6.11"] = traction_val[1][0]
            self.taux_m_rd["equ6.12"] = traction_val[1][1]

        return (latex, self.taux_m_rd)
