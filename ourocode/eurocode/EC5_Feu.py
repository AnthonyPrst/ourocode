#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import math as mt
from math import sqrt, pi, cos, sin, radians
import numpy as np

import forallpeople as si

si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.A0_Projet import Projet
from ourocode.eurocode.EC5_Element_droit import (
    Barre,
    Flexion,
    Traction,
    Compression,
    Cisaillement,
)

def interpolation_lineaire(x, xa, xb, ya, yb):
    """Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb"""
    y = ya + (x - xa) * ((yb - ya) / (xb - xa))
    return y
# ================================ GLOBAL ==================================


class Feu(Barre):
    PROTECTION = (
        "Aucune protection",
        "1 plaque de platre type A joints comblés",
        "1 plaque de platre type A joints vides",
        "1 plaque de platre type F joints comblés",
        "1 plaque de platre type F joints vides",
        "1 plaque de platre type H joints comblés",
        "1 plaque de platre type H joints vides",
        "Fibre de roche",
        "Panneautage bois",
        "Contreplaqué",
        "Panneaux de fibres ou de particules",
    )
    D0 = 7*si.mm
    KMOD_FI = 1

    def __init__(
        self,
        t_expo: si.min=30,
        protection_haut: str = PROTECTION,
        protection_bas: str = PROTECTION,
        protection_gauche: str = PROTECTION,
        protection_droite: str = PROTECTION,
        double_couches: bool = False,
        hp: si.mm=0,
        rho_k_protect: float=0,
        tf: si.min=None,
        **kwargs,
    ):
        """Classe qui définit les caractéristiques d'un élément droit au feu.
        Cette classe est hérité de la classe Barre du module EC5_Element_droit.py.

        Args:
            t_expo (int, optional): Durée d'exposition au feu en minutes. Defaults to 30.
            
            Attention : 
                Pour la détermination des joints vides ou comblés en cas de protection rapportée:
                    un joint est considéré comme comblé si le vide est <= à 2mm.

            protection_haut (str): Type de protection au feu sur le haut de l'élément.
            protection_bas (str): Type de protection au feu sur le bas de l'élément.
            protection_gauche (str): Type de protection au feu sur le gauche de l'élément.
            protection_droite (str): Type de protection au feu sur le droite de l'élément.
            hp (int)): épaisseur totale (si double couches) des panneaux de protection en mm. Defaults to 0.
            rho_k_protect (float): A définir uniquement si il y a un ou des panneaux de protection bois ou en fibre de roche. 
                Masse volumique des panneaux de protection en kg/m3. Defaults to 0.
            tf (int, optional): Durée avant rupture du matériaux de protection au feu en minutes. Defaults to None.
                Attention uniquement pour les plaques de platre de type F ou les panneaux de fibres de roche.
                On ne peut par conséquent pas définir des plaques de platre de type F en même temps que des panneaux de fibres de roche.
        """
        super().__init__(**kwargs)
        self.t_expo = t_expo * si.min
        self.protection = {
            "haut": protection_haut,
            "bas": protection_bas,
            "gauche": protection_gauche,
            "droite": protection_droite,
        }
        self.double_couches = double_couches
        self.hp = hp * si.mm
        self.rho_k_protect = rho_k_protect
        self.t_f = tf
        self._get_section_reduction()

    def __convert_latex_ftyped(self, latex: str, type_caract: str):
        end_index_rk = type_caract.find("k")
        type_caract = type_caract[1:end_index_rk]
        latex = latex.replace("f_{type_{Rd}", "f_{" + type_caract + "_{" + "d_fi}")
        latex = latex.replace("f_{type_{k}", "f_{" + type_caract + "_{" + "k_fi}")
        return latex

    @property
    def K_fi(self):
        data_csv_kfi = self._data_from_csv("kfi.csv")
        return float(data_csv_kfi.loc[self.type_bois]["kfi"])

    def _get_bar_beta0_and_betan(self):
        beta_0 = None
        beta_n = None
        rho_k = int(self.caract_meca.loc["rhok"])
        if self.classe[0:1] == "C" and rho_k >= 290:
            beta_0 = 0.65
            beta_n = 0.7
        elif self.classe[0:2] == "GL" and rho_k >= 290:
            beta_0 = 0.65
            beta_n = 0.8
        elif self.classe[0:1] == "D":
            if rho_k >= 290 and rho_k < 450:
                beta_0 = interpolation_lineaire(rho_k, 290, 450, 0.65, 0.5)
                beta_n = interpolation_lineaire(rho_k, 290, 450, 0.7, 0.55)
            elif rho_k >= 450:
                beta_0 = 0.5
                beta_n = 0.55
        elif self.type_bois == "LVL" and rho_k >= 480:
            beta_0 = 0.65
            beta_n = 0.7
        elif self.b_calcul.value*10**3 >= 20 and rho_k >= 450:
            if self.type_bois == "CP":
                beta_0 = 1
            else:
                beta_0 = 0.9
        elif self.b_calcul.value*10**3 < 20 or rho_k < 450:
            kp = mt.sqrt(450 / rho_k)
            kh = mt.sqrt(20 / self.b_calcul.value * 10**3)
            if self.type_bois == "CP":
                beta_0 = 1
            else:
                beta_0 = 0.9
            beta_0 = beta_0 * kp * kh
        return beta_0, beta_n
    
    def _get_wood_protect_beta0_and_betan(self, type_wood_protect: str):
        beta_0 = None
        beta_n = None
        if self.hp.value*10**3 >= 20 and self.rho_k_protect >= 450:
            if self.type_bois == "Contreplaqué":
                beta_0 = 1
            else:
                beta_0 = 0.9
        elif self.hp.value*10**3 < 20 or self.rho_k_protect < 450:
            kp = mt.sqrt(450 / self.rho_k_protect)
            kh = mt.sqrt(20 / self.hp.value * 10**3)
            if self.type_bois == "Contreplaqué":
                beta_0 = 1
            else:
                beta_0 = 0.9
            beta_0 = beta_0 * kp * kh
        return beta_0, beta_n
        
        
    def d_char_0(self, t:int, beta_0:float):
        """Retourne la valeur de la profondeur de carbonisation en mm pour une carbonisation uni-dimensionnelle
        Args:
            t (int): le temps approprié d'exposition au feu en minutes.
        """
        @handcalc(
            override="short",
            precision=2,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            d_char_0 = beta_0 * t
            return d_char_0
        return val()
    
    def _d_char_n(self, t:int, beta_n:float):
        """Retourne la valeur de la profondeur de carbonisation fictive qui tient compte de l'effet des arrondis en coins
        Args:
            t (int): le temps approprié d'exposition au feu en minutes.
        """
        return beta_n * t

    def _k0(self, t:int):
        if t < 20:
            return t/20
        else:
            return 1
        
    def _d_ef(self, d_char_n: si.mm, t:si.min):
        """Retourne la valeur de la profondeur de carbonisation effective en mm
        Args:
            d_char_n (si.mm): profondeur de carbonisation fictive
            t (int): le temps approprié d'exposition au feu en minutes.
        """
        d_0 = self.D0
        k_0 = self._k0(t)
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            d_char_n = d_char_n + k_0 * d_0
            return d_char_n
        return val()
    

    def _get_tch_and_tf(self, orientation: str, beta_0: float):
        """Retorune le temps de démarrage de la carbonisation en minute (tch) et le temps de rupture de la protection (tf) en minute

        Args:
            orientation (str, optional): Orientation sur la barre. Defaults to "Bas".
        """
        h_p = self.hp.value*10**3
        t_f = self.t_f

        if self.protection[orientation] in ("Panneautage bois","Contreplaqué","Panneaux de fibres ou de particules"):
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_ch =  self.hp.value*10**3 / beta_0 # 3.10 (min)
                t_f = t_ch # 3.14 (min)
                return t_ch, t_f

        # platre à joint vide A,H,F
        elif "joints comblés" in self.protection[orientation]:
            if self.double_couches:
                if "Type F" in self.protection[orientation]:
                    h_p = (h_p / 2) * 1.8
                else:
                    h_p = (h_p / 2) * 1.5
            
            if "Type A" in self.protection[orientation] or "Type H" in self.protection[orientation]:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 14 # 3.11 (min)
                    t_f = t_ch # 3.15 (min)
                    return t_ch, t_f
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 14 # 3.11 (min)
                    t_f # Temps de rupture de la plaque de platre type F
                    return t_ch, t_f

        # platre à joint comblé A,H,F
        elif "joints vides" in self.protection[orientation]:
            if "Type A" in self.protection[orientation] or "Type H" in self.protection[orientation]:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 23 # 3.12 (min)
                    t_f = t_ch # 3.15 (min)
                    return t_ch, t_f
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 23 # 3.12 (min)
                    t_f # Temps de rupture de la plaque de platre type F
                    return t_ch, t_f
        # Fibre de roche
        elif self.protection[orientation] == "Fibre de roche":
            rho_k_protect = self.rho_k_protect
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_ch = 0.07 * (h_p-20) * sqrt(rho_k_protect) # 3.13 (min)
                t_f # Temps de rupture de la plaque en fibre de roche
                return t_ch, t_f
        # Aucune protection
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_ch = 0 # 3.12 (min)
                t_f = t_ch # 3.15 (min)
                return t_ch, t_f
        return val()
    
    def _get_k2(self, orientation:str):
        k2 = 1
        hp = self.hp.value*10**3
        if "Type F" in self.protection[orientation]:
            if self.double_couches:
                hp = self.hp.value*10**3 / 2
            k2 = 1 - 0.018 * hp # 3.7
        elif self.protection[orientation] == "Fibre de roche":
            if 45 > hp >= 20:
                k2 = interpolation_lineaire(hp, 20,  45, 1, 0.6)
            elif hp >= 45:
                k2 = 0.6
        return k2

    def _get_ta(self, t_ch: int, t_f: int, beta_n: float, k_2:float):
        """Retourne le temps de rupture de la protection en minutes (ta)
        Args:
            t_ch (int): Temps de démarrage de la carbonisation en minutes
            t_f (int): Temps de rupture de la protection en minutes
            beta_n (float): Coefficient de carbonisation fictif
            k_2 (float): Coefficient d'accelération de la carbonisation'"""
        if t_f:
            if t_ch < t_f :
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_a = (25 - (t_f - t_ch) * k_2 * beta_n) / (2 * beta_n) + t_f # 3.9
                    return t_a
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_a = min(2*t_f, 12.5/beta_n + t_f) # 3.8
                    return t_a
            return val()
        else:
            return ("", 0)
    
    def _get_section_reduction(self):
        self._def = {}
        for key, orientation in self.protection.items():
            d_char_n = 0
            beta_0, beta_n = self._get_bar_beta0_and_betan(orientation)
            beta_0_protect, beta_n_protect = self._get_wood_protect_beta0_and_betan(orientation)
            latex_tch_tf, res_tch_tf = self._get_tch_and_tf(orientation, beta_0_protect)
            tch = res_tch_tf[0]
            tf = res_tch_tf[1]
            k2 = self._get_k2(orientation)
            latex_ta, ta = self._get_ta(tch, tf, beta_n, k2)
            if tch < tf:
                d_char_n = self._d_char_n(tf-tch, beta_n*k2)
            k3 = 2
            d_char_n += self._d_char_n(ta-tf, beta_n*k3)
            d_char_n += self._d_char_n(self.t_expo-ta, beta_n*k3)
            d_ef = self._d_ef(d_char_n, self.t_expo)
            self._def[orientation] = (latex_tch_tf + latex_ta + d_ef[0], d_ef[1])
        print(self._def[orientation])
    

    def _f_type_d(self, typeCarac=Barre.CARACTERISTIQUE[0:6]):
        """Méthode donnant la résistance de calcul de l'élément fonction de la vérification

        Args:
            typeCarac (str, optional): Type de résistance caractéristique (flexion = "fm0k", compression = "fc0k" etc.). Defaults to "fm0k".
            loadtype (str, optional): Durée de chargement (Permanente, Court terme etc.). Defaults to "Permanente".
            typecombi (str, optional): Type de combinaison étudiée ("Fondamentales" ou " Accidentelles"). Defaults to "Fondamentales".

        Returns:
            float: Résistance de calcul en N/mm2 du type de vérification étudié.
        """
        gamma_M_fi = self._get_gamma_M("Accidentelles")
        K_mod_fi = self.KMOD_FI  # en situation à chaud
        K_fi = self.K_fi
        f_type_k = float(self.caract_meca.loc[typeCarac]) * si.MPa

        if typeCarac == "fm0k" and self.k_sys > 1:
            k_sys = self.k_sys

            @handcalc(
                override="short",
                precision=2,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                f_type_Rd = k_sys * f_type_k * K_fi * K_mod_fi / gamma_M_fi
                return f_type_Rd

        else:

            @handcalc(
                override="short",
                precision=2,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                f_type_Rd = f_type_k * K_fi * K_mod_fi / gamma_M_fi
                return f_type_Rd

        value = val()
        latex = self.__convert_latex_ftyped(value[0], typeCarac)
        self.f_type_rd = value[1]
        return (latex, value[1])

    def _K_h(self):
        """Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k"""
        kh = {}
        dim = {"y": self.h_calcul.value * 10**3, "z": self.b_calcul.value * 10**3}

        for cle, valeur in dim.items():
            if self.type_bois == "Massif":
                if valeur < 150:
                    kh[cle] = min((150 / valeur) ** 0.2, 1.3)
                else:
                    kh[cle] = 1
            elif self.type_bois == "BLC":
                if valeur < 600:
                    kh[cle] = min((600 / valeur) ** 0.1, 1.1)
                else:
                    kh[cle] = 1
            else:
                print("LVL non pris en compte dans cette fonction")
                kh[cle] = 1
        return kh

    def Emean_fin(self, psy_2: float):
        """Renvoie le E,mean,fin en fonction du Kdef et du psy2"""
        self.psy_2 = psy_2
        E0_mean = int(self.caract_meca.loc["E0mean"]) * si.MPa
        K_fi = self.K_fi
        K_def = self.K_def

        @handcalc(
            override="short",
            precision=2,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            E_mean_fin_fi = E0_mean * K_fi / (1 + psy_2 * K_def)
            return E_mean_fin_fi

        value = val()
        self.E_mean_fin = value[1]
        return value


# ================================ FLEXION ==================================


class Flexion_feu(Feu, Flexion):
    COEF_LEF = {"Appuis simple": [1, 0.9, 0.8], "Porte à faux": [0.5, 0.8]}
    LOAD_POS = {
        "Charge sur fibre comprimée": 0,
        "Charge sur fibre neutre": 1,
        "Charge sur fibre tendu": 2,
    }

    def __init__(
        self,
        lo: si.mm,
        coeflef: float = COEF_LEF["Appuis simple"][1],
        pos: str = LOAD_POS,
        **kwargs,
    ):
        """Classe permettant le calcul de la flexion d'une poutre bois selon l'EN 1995 §6.1.6, §6.2.3, §6.2.4 et §6.3.3.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.

        Args:
            lo (int): longueur de déversemment en mm
            coeflef (float): appuis simple :
                                                    Moment constant : 1
                                                    Charge répartie constante : 0.9
                                                    Charge concentrée au milieu de la portée : 0.8
                                    porte à faux :
                                                    Charge répartie constante : 0.5
                                                    Charge concentrée agissant à l'extrémité libre : 0.8.

            pos (str): positionnement de la charge sur la hauteur de poutre
        """
        Feu.__init__(**kwargs)
        Flexion.__init__(lo, coeflef, pos, **kwargs)

    @property
    def sigma_m_crit(self):
        """Retourne sigma m,crit pour la prise en compte du déversement d'une poutre"""
        self.l_ef = self.lo * self.coeflef
        if self.pos == 0:
            self.l_ef = self.l_ef + 2 * self.h_calcul
        elif self.pos == 2:
            self.l_ef = self.l_ef - 0.5 * self.h_calcul

        b_calcul = self.b_calcul
        h_calcul = self.h_calcul
        l_ef = self.l_ef
        E_0_05 = int(self.caract_meca.loc["E005"]) * si.MPa
        K_fi = self.K_fi

        @handcalc(
            override="short",
            precision=2,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            sigma_m_crit_fi = (0.78 * b_calcul**2 * E_0_05 * K_fi) / (h_calcul * l_ef)
            return sigma_m_crit_fi

        return val()

    @property
    def lamb_rel_m(self):
        """Retourne l'élancement relatif de la section avec pour argument"""
        f_m0k = float(self.caract_meca.loc["fm0k"]) * si.MPa
        K_fi = self.K_fi
        sigma_m_crit = self.sigma_m_crit[1]

        @handcalc(
            override="short",
            precision=2,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            lamb_rel_m_fi = sqrt(f_m0k * K_fi / sigma_m_crit)
            return lamb_rel_m_fi

        return val()

    @property
    def K_crit(self):
        """Retourne K,crit le coef. de minoration de la résistance en flexion au déversement"""
        lamb_rel_m_fi = self.lamb_rel_m[1]

        if lamb_rel_m_fi <= 0.75:

            @handcalc(
                override="short",
                precision=2,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                K_crit_fi = 1
                return K_crit_fi

        elif 0.75 < lamb_rel_m_fi <= 1.4:

            @handcalc(
                override="short",
                precision=2,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                K_crit_fi = 1.56 - 0.75 * lamb_rel_m_fi
                return K_crit_fi

        else:

            @handcalc(
                override="short",
                precision=2,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                K_crit_fi = 1 / (lamb_rel_m_fi**2)
                return K_crit_fi

        return val()

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
        K_crit_fi = self.K_crit[1]

        @handcalc(
            override="short",
            precision=3,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def base():
            taux_6_11 = sigma_my_d_fi / (f_m_d_fi * K_h_y) + K_m * sigma_mz_d_fi / (
                f_m_d_fi * K_h_z
            )  # equ6.11
            taux_6_12 = K_m * sigma_my_d_fi / (f_m_d_fi * K_h_y) + sigma_mz_d_fi / (
                f_m_d_fi * K_h_z
            )  # equ6.12
            taux_6_33y = sigma_my_d_fi / (f_m_d_fi * K_h_y * K_crit_fi)  # equ6.33
            taux_6_33z = sigma_mz_d_fi / (f_m_d_fi * K_h_z * K_crit_fi)  # equ6.33
            return taux_6_11, taux_6_12, taux_6_33y, taux_6_33z

        base_val = base()
        latex = base_val[0]
        self.taux_m_rd["equ6.11"] = base_val[1][0]
        self.taux_m_rd["equ6.12"] = base_val[1][1]
        self.taux_m_rd["equ6.33y"] = base_val[1][2]
        self.taux_m_rd["equ6.33z"] = base_val[1][3]

        if compression and isinstance(compression, Compression):
            sigma_c_0_d = compression.sigma_c_0_rd
            f_c_0_d = compression.f_type_rd
            K_c_y = compression.kc_Axe[1]["y"]
            K_c_z = compression.kc_Axe[1]["z"]
            taux_6_2 = compression.taux_c_0_rd["equ6.2"]

            # taux_6_23 = compression.taux_c_0_rd['equ6.23']
            # taux_6_24 = compression.taux_c_0_rd['equ6.24']
            @handcalc(
                override="short",
                precision=3,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def comp(taux_6_11, taux_6_12, taux_6_2, taux_6_33y, taux_6_33z):
                taux_6_19 = taux_6_2**2 + taux_6_11  # equ6.19
                taux_6_20 = taux_6_2**2 + taux_6_12  # equ6.20
                taux_6_23 = sigma_c_0_d / (f_c_0_d * K_c_y)  # equ6.23
                taux_6_24 = sigma_c_0_d / (f_c_0_d * K_c_z)  # equ6.24
                taux_6_35zyz = (
                    taux_6_33y**2 + (sigma_mz_d_fi / (f_m_d_fi * K_h_z)) + taux_6_24
                )  # equ6.35
                taux_6_35yzz = (
                    taux_6_33y + (sigma_mz_d_fi / (f_m_d_fi * K_h_z)) ** 2 + taux_6_24
                )  # equ6.35 interprétation
                taux_6_35yzy = (
                    taux_6_33z**2 + (sigma_my_d_fi / (f_m_d_fi * K_h_y)) + taux_6_23
                )  # equ6.35
                taux_6_35zyy = (
                    taux_6_33z + (sigma_my_d_fi / (f_m_d_fi * K_h_y)) ** 2 + taux_6_23
                )  # equ6.35 interprétation
                return (
                    taux_6_19,
                    taux_6_20,
                    taux_6_35zyz,
                    taux_6_35yzz,
                    taux_6_35yzy,
                    taux_6_35zyy,
                )

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
            self.taux_m_rd["equ6.35zyz"] = compression_val[1][
                2
            ]  # 1er item axe de flexion pas au carré, 2eme item axe de flexion au carré, 3eme axe de compression
            self.taux_m_rd["equ6.35yzz"] = compression_val[1][3]
            self.taux_m_rd["equ6.35yzy"] = compression_val[1][4]
            self.taux_m_rd["equ6.35zyy"] = compression_val[1][5]

        if traction and isinstance(traction, Traction):
            taux_6_1 = traction.taux_t_0_rd["equ6.1"]

            @handcalc(
                override="short",
                precision=3,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def tract(taux_6_11, taux_6_12):
                taux_6_17 = taux_6_11 + taux_6_1  # equ6.17
                taux_6_18 = taux_6_12 + taux_6_1  # equ6.18
                return taux_6_17, taux_6_18

            traction_val = tract(self.taux_m_rd["equ6.11"], self.taux_m_rd["equ6.12"])
            latex = latex + traction_val[0]
            self.taux_m_rd["equ6.11"] = traction_val[1][0]
            self.taux_m_rd["equ6.12"] = traction_val[1][1]

        return (latex, self.taux_m_rd)


# ================================ Traction ==================================


class Traction(Feu, Traction):
    def __init__(self, **kwargs):
        """Classe permettant le calcul de la Traction d'un élément bois selon l'EN 1995.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.
        """
        Feu.__init__(**kwargs)
        Traction.__init__(**kwargs)

    def taux_t_0_d(self):
        """Retourne le taux de travail en traction axial.

        Returns:
            float: taux de travail en %
        """
        self.taux_t_0_rd = {}
        K_h_y = self.K_h["y"]
        K_h_z = self.K_h["z"]
        sigma_t_0_d_fi = self.sigma_t_0_rd
        f_t_0_d_fi = self.f_type_rd

        @handcalc(
            override="short",
            precision=3,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            K_h = min(K_h_y, K_h_z)
            taux_6_1 = sigma_t_0_d_fi / (K_h * f_t_0_d_fi)  # equ6.1
            return taux_6_1

        value = val()

        self.taux_t_0_rd["equ6.1"] = value[1]
        return value


# ================================ Compression ==================================


class Compression_feu(Feu, Compression):
    COEF_LF = {
        "Encastré 1 côté": 2,
        "Rotule - Rotule": 1,
        "Encastré - Rotule": 0.7,
        "Encastré - Encastré": 0.5,
        "Encastré - Rouleau": 1,
    }

    def __init__(self, lo_y: si.mm, lo_z: si.mm, type_appuis: str = COEF_LF, **kwargs):
        """Classe permettant le calcul de la Compression d'un élément bois selon l'EN 1995.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.

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
        Feu.__init__(**kwargs)
        Compression.__init__(lo_y, lo_z, type_appuis, **kwargs)

    @property
    def lamb_rel_Axe(self):
        """Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation"""
        lamb_y = self.lamb[1]["y"]
        lamb_z = self.lamb[1]["z"]
        f_c0k = float(self.caract_meca.loc["fc0k"]) * si.MPa
        E_0_05 = int(self.caract_meca.loc["E005"]) * si.MPa
        K_fi = self.K_fi

        @handcalc(
            override="short",
            precision=3,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            lamb_rel_y_fi = (lamb_y / pi) * sqrt(f_c0k * K_fi / (E_0_05 * K_fi))
            lamb_rel_z_fi = (lamb_z / pi) * sqrt(f_c0k * K_fi / (E_0_05 * K_fi))
            return {"y": lamb_rel_y_fi, "z": lamb_rel_z_fi}

        return val()

    @property
    def k_Axe(self):
        """Retourne le facteur Ky ou Kz (fonction de l'axe de flambement)"""
        beta_C = self.beta_C
        lamb_rel_y_fi = self.lamb_rel_Axe[1]["y"]
        lamb_rel_z_fi = self.lamb_rel_Axe[1]["z"]

        @handcalc(
            override="short",
            precision=2,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            k_y_fi = 0.5 * (1 + beta_C * (lamb_rel_y_fi - 0.3) + lamb_rel_y_fi**2)
            k_z_fi = 0.5 * (1 + beta_C * (lamb_rel_z_fi - 0.3) + lamb_rel_z_fi**2)
            return {"y": k_y_fi, "z": k_z_fi}

        return val()

    @property
    def kc_Axe(self):
        """Retourne le coefficient multiplicateur KcAxe  (axe = y ou z suivant axe de rotation en flambement) de fc,0,d"""
        k_y_fi = self.k_Axe[1]["y"]
        k_z_fi = self.k_Axe[1]["z"]
        lamb_rel_y_fi = self.lamb_rel_Axe[1]["y"]
        lamb_rel_z_fi = self.lamb_rel_Axe[1]["z"]

        @handcalc(
            override="short",
            precision=2,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            k_c_y = 1 / (k_y_fi + sqrt(k_y_fi**2 - lamb_rel_y_fi**2))
            k_c_z = 1 / (k_z_fi + sqrt(k_z_fi**2 - lamb_rel_z_fi**2))
            return {"y": min(k_c_y, 1), "z": min(k_c_z, 1)}

        return val()

    def taux_c_0_d(self, flexion: object = None):
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

        if flexion and isinstance(flexion, Flexion):
            taux_6_11 = flexion.taux_m_rd["equ6.11"]
            taux_6_12 = flexion.taux_m_rd["equ6.12"]
        else:
            taux_6_11 = 0
            taux_6_12 = 0

        if lamb_rel_y < 0.3 and lamb_rel_z < 0.3:

            @handcalc(
                override="short",
                precision=3,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                taux_6_2 = sigma_c_0_d_fi / f_c_0_d_fi  # equ6.2
                taux_6_19 = (
                    sigma_c_0_d_fi / (f_c_0_d_fi * K_c_y_fi)
                ) ** 2 + taux_6_11  # equ6.19
                taux_6_20 = (
                    sigma_c_0_d_fi / (f_c_0_d_fi * K_c_z_fi)
                ) ** 2 + taux_6_12  # equ6.20
                return taux_6_2, taux_6_19, taux_6_20

            value = val()
            self.taux_c_0_rd["equ6.19"] = value[1][1]
            self.taux_c_0_rd["equ6.20"] = value[1][2]
        else:

            @handcalc(
                override="short",
                precision=3,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                taux_6_2 = sigma_c_0_d_fi / f_c_0_d_fi  # equ6.2
                taux_6_23 = (
                    sigma_c_0_d_fi / (f_c_0_d_fi * K_c_y_fi) + taux_6_11
                )  # equ6.23
                taux_6_24 = (
                    sigma_c_0_d_fi / (f_c_0_d_fi * K_c_z_fi) + taux_6_12
                )  # equ6.24
                return taux_6_2, taux_6_23, taux_6_24

            value = val()
            self.taux_c_0_rd["equ6.23"] = value[1][1]
            self.taux_c_0_rd["equ6.24"] = value[1][2]

        self.taux_c_0_rd["equ6.2"] = value[1][0]
        return value


# ================================ CISAILLEMENT ==================================


class Cisaillement_feu(Feu, Cisaillement):
    def __init__(self, **kwargs):
        """Classe qui permet de calculer le cisaillement d'une poutre comme décrit à l'EN 1995 §6.1.7 et §6.5.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.
        """
        super(Feu, self).__init__(**kwargs)
        super(Cisaillement, self).__init__(**kwargs)

    def Kv(self, hef: int, x: int, i_lo: int, ent=("Dessous", "Dessus")):
        """Retourne le facteur d'entaille Kv pour une entaille au niveau d'un appuis

        Args:
            hef (int): Hauteur efficace de la poutre (hauteur - hauteur de l'entaille) en mm
            x (int):Distance entre le centre de réaction à l'appuis et le coin de l'entaille en mm
            i_lo (float): longueur horizontal de l'entaille en mm
            ent (tuple, optional): Entaille sur le dessus ou dessous de la poutre.

        Returns:
            float: facteur Kv
        """
        K_n = self.DICT_KN[self.type_bois]
        x = x * si.mm
        h_ef = hef * si.mm
        i = i_lo * si.mm / h_ef
        h_calcul = self.h_calcul

        if ent == "Dessus":
            self.K_v = 1
            return self.K_v
        else:

            @handcalc(
                override="long",
                precision=3,
                jupyter_display=self.JUPYTER_DISPLAY,
                left="\\[",
                right="\\]",
            )
            def val():
                alpha = h_ef / h_calcul
                K_v = min(
                    1,
                    (K_n * (1 + (1.1 * i**1.5) / sqrt(h_calcul)))
                    / (
                        sqrt(h_calcul)
                        * (
                            sqrt(alpha * (1 - alpha))
                            + 0.8 * x / h_calcul * sqrt(1 / alpha - alpha**2)
                        )
                    ),
                )
                return K_v

            self.h_ef = h_ef
            value = val()
            self.K_v = value[1]
            return value

    def taux_tau_d(self):
        """Retourne le taux de travail en cisaillement en %"""
        self.taux_tau_rd = {}
        tau_d_fi = self.tau_rd
        f_v_d_fi = self.f_type_rd
        K_v = self.K_v

        @handcalc(
            override="short",
            precision=3,
            jupyter_display=self.JUPYTER_DISPLAY,
            left="\\[",
            right="\\]",
        )
        def val():
            taux_6_13 = tau_d_fi / f_v_d_fi
            taux_6_60 = tau_d_fi / (K_v * f_v_d_fi)
            return taux_6_13, taux_6_60

        value = val()
        self.taux_tau_rd["equ6.13"] = value[1][0]
        self.taux_tau_rd["equ6.60"] = value[1][1]
        return value


if __name__=='__main__':
    beam = Barre(60,200,"Rectangulaire", classe="C24", cs=1)
#     beam3 = Barre(60,100,"Rectangulaire", classe="C24", cs=1)
#     beam_ass = Poutre_assemblee_meca(beam_2=beam2, l=5000, disposition="Latérale", recouvrement=[0,120], Kser=[None,None,700], entraxe=[None,None,250], psy_2=0, beam_3=beam3)
#     pole_ass = Poteau_assemble_meca._from_parent_class(beam_ass, lo_y=5000, lo_z=5000, type_appuis="Rotule - Rotule")
#     # print(pole_ass.__dict__)
#     # print(pole_ass.lamb)
#     print(pole_ass.kc_Axe)
