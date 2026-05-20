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
from ourocode.eurocode.mixins.math_utils import MathUtilsMixin

from ourocode.eurocode.ec5.element_droit.barre import Barre

# ================================ GLOBAL ==================================


class Feu(Barre):
    EXPOSITION = (
        "Pas d'exposition",
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
    D0 = 7 * si.mm
    KMOD_FI = 1
    ORIENTATION = ("Haut", "Bas", "Gauche", "Droite")

    def __init__(
        self,
        t_expo: int = 30,
        haut: str = EXPOSITION,
        bas: str = EXPOSITION,
        gauche: str = EXPOSITION,
        droite: str = EXPOSITION,
        double_couches: bool = ("False", "True"),
        hp: si.mm = 0,
        rho_k_protect: float = 0,
        tf: int = None,
        **kwargs,
    ):
        """Classe qui définit les caractéristiques d'un élément droit au feu.
        Cette classe est hérité de la classe Barre du module EC5_Element_droit.py.

        Args:
            t_expo (int, optional): Durée d'exposition au feu en minutes. Defaults to 30.
            haut (str): Type d'expostion ou de protection au feu sur le haut de l'élément.
            bas (str): Type d'expostion ou de protection au feu sur le bas de l'élément.
            gauche (str): Type d'expostion ou de protection au feu sur le gauche de l'élément.
            droite (str): Type d'expostion ou de protection au feu sur le droite de l'élément.
            hp (int): Épaisseur totale (si double couches) des panneaux de protection en mm. Defaults to 0.
            rho_k_protect (float): Masse volumique des panneaux de protection bois ou fibre de roche en kg/m3. Defaults to 0.
            tf (int, optional): Durée avant rupture du matériau de protection au feu en minutes. Defaults to None.
                Uniquement pour les plaques de plâtre de type F ou les panneaux de fibres de roche.

        Note:
            Pour la détermination des joints vides ou comblés en cas de protection rapportée :
            un joint est considéré comme comblé si le vide est ≤ 2 mm.
        """
        super().__init__(**kwargs)
        self.t_expo = t_expo
        self.haut = haut
        self.bas = bas
        self.gauche = gauche
        self.droite = droite
        self.protection = {
            "Haut": haut,
            "Bas": bas,
            "Gauche": gauche,
            "Droite": droite,
        }
        self.double_couches = double_couches
        self.hp = hp * si.mm
        self.rho_k_protect = rho_k_protect
        self.tf = tf
        self._get_section_reduction()

    def __convert_latex_ftyped(self, latex: str, type_caract: str):
        end_index_rk = type_caract.find("k")
        type_caract = type_caract[1:end_index_rk]
        latex = latex.replace("f_{type_{Rd}", "f_{" + type_caract + "_{" + "d_fi}")
        latex = latex.replace("f_{type_{k}", "f_{" + type_caract + "_{" + "k_fi}")
        return latex

    @property
    def K_fi(self) -> float:
        data_csv_kfi = self._data_from_csv("kfi.csv")
        return float(data_csv_kfi.loc[self.type_bois]["kfi"])

    def _get_bar_beta0_and_betan(self):
        beta_0 = None
        beta_n = None
        rho_k = int(self.caract_meca.loc["rhok"])
        if self.classe[0:1] == "C" and rho_k >= 290:
            beta_0 = 0.65
            beta_n = 0.8
        elif self.classe[0:2] == "GL" and rho_k >= 290:
            beta_0 = 0.65
            beta_n = 0.7
        elif self.classe[0:1] == "D":
            if rho_k >= 290 and rho_k < 450:
                beta_0 = MathUtilsMixin.interpolation_lineaire(rho_k, 290, 450, 0.65, 0.5)
                beta_n = MathUtilsMixin.interpolation_lineaire(rho_k, 290, 450, 0.7, 0.55)
            elif rho_k >= 450:
                beta_0 = 0.5
                beta_n = 0.55
        elif self.type_bois == "LVL" and rho_k >= 480:
            beta_0 = 0.65
            beta_n = 0.7
        elif self.b_calcul.value * 10**3 >= 20 and rho_k >= 450:
            if self.type_bois == "CP":
                beta_0 = 1
            else:
                beta_0 = 0.9
        elif self.b_calcul.value * 10**3 < 20 or rho_k < 450:
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
        if type_wood_protect in (
            "Panneautage bois",
            "Contreplaqué",
            "Panneaux de fibres ou de particules",
        ):
            if self.hp.value * 10**3 >= 20 and self.rho_k_protect >= 450:
                if type_wood_protect == "Contreplaqué":
                    beta_0 = 1
                else:
                    beta_0 = 0.9
            elif self.hp.value * 10**3 < 20 or self.rho_k_protect < 450:
                kp = mt.sqrt(450 / self.rho_k_protect)
                kh = mt.sqrt(20 / self.hp.value * 10**3)
                if type_wood_protect == "Contreplaqué":
                    beta_0 = 1
                else:
                    beta_0 = 0.9
                beta_0 = beta_0 * kp * kh
        return beta_0, beta_n

    def d_char_0(self, t: int, beta_0: float) -> tuple:
        """Retourne la valeur de la profondeur de carbonisation en mm pour une carbonisation uni-dimensionnelle
        Args:
            t (int): le temps approprié d'exposition au feu en minutes.
        """

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            d_char_0 = beta_0 * t
            return d_char_0

        return val()

    def _d_char_n(self, t: int, beta_n: float):
        """Retourne la valeur de la profondeur de carbonisation fictive qui tient compte de l'effet des arrondis en coins
        Args:
            t (int): le temps approprié d'exposition au feu en minutes.
        """
        return beta_n * t * si.mm

    def _k0(self, t: int):
        if t < 20:
            return t / 20
        else:
            return 1

    def _d_ef(self, d_char_n: si.mm, t: int):
        """Retourne la valeur de la profondeur de carbonisation effective en mm
        Args:
            d_char_n (si.mm): profondeur de carbonisation fictive
            t (int): le temps approprié d'exposition au feu en minutes.
        """
        d_0 = self.D0
        k_0 = self._k0(t)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            d_ef = d_char_n + k_0 * d_0
            return d_ef
        return val()

    def _get_tch_and_tf(self, orientation: str, beta_0: float):
        """Retorune le temps de démarrage de la carbonisation en minute (tch) et le temps de rupture de la protection (tf) en minute

        Args:
            orientation (str, optional): Orientation sur la barre. Defaults to "Bas".
        """
        h_p = self.hp.value * 10**3
        t_f = self.tf

        if self.protection[orientation] in (
            "Panneautage bois",
            "Contreplaqué",
            "Panneaux de fibres ou de particules",
        ):

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_ch = self.hp.value * 10**3 / beta_0  # en min/equa3.10
                t_f = t_ch  # en min/equa3.14
                return t_ch, t_f

        # platre à joint vide A,H,F
        elif "joints comblés" in self.protection[orientation]:
            if self.double_couches:
                if "type F" in self.protection[orientation]:
                    h_p = (h_p / 2) * 1.8
                else:
                    h_p = (h_p / 2) * 1.5

            if (
                "type A" in self.protection[orientation]
                or "type H" in self.protection[orientation]
            ):

                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 14  # en min/equa3.11
                    t_f = t_ch  # en min/equa3.15
                    return t_ch, t_f

            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 14  # en min/equa3.11
                    t_f # Temps de rupture de la plaque de platre type F
                    return t_ch, t_f

        # platre à joint comblé A,H,F
        elif "joints vides" in self.protection[orientation]:
            if (
                "type A" in self.protection[orientation]
                or "type H" in self.protection[orientation]
            ):
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 23  # en min/equa3.12
                    t_f = t_ch  # en min/equa3.15
                    return t_ch, t_f

            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_ch = 2.8 * h_p - 23  # en min/equa3.12
                    t_f  # Temps de rupture de la plaque de platre type F
                    return t_ch, t_f

        # Fibre de roche
        elif self.protection[orientation] == "Fibre de roche":
            rho_k_protect = self.rho_k_protect

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_ch = 0.07 * (h_p - 20) * sqrt(rho_k_protect)  # en min/equa3.13
                t_f  # Temps de rupture de la plaque en fibre de roche
                return t_ch, t_f

        # Aucune protection
        else:

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_ch = 0  # en min/equa3.12
                t_f = t_ch  # en min/equa3.15
                return t_ch, t_f

        return val()

    def _get_k2(self, orientation: str):
        k2 = 1
        hp = self.hp.value * 10**3
        if "type F" in self.protection[orientation]:
            if self.double_couches:
                hp = self.hp.value * 10**3 / 2
            k2 = 1 - 0.018 * hp  # 3.7
        elif self.protection[orientation] == "Fibre de roche":
            if 45 > hp >= 20:
                k2 = MathUtilsMixin.interpolation_lineaire(hp, 20, 45, 1, 0.6)
            elif hp >= 45:
                k2 = 0.6
        return k2

    def _get_ta(self, t_ch: int, t_f: int, beta_n: float, k_2: float):
        """Retourne le temps de rupture de la protection en minutes (ta)
        Args:
            t_ch (int): Temps de démarrage de la carbonisation en minutes
            t_f (int): Temps de rupture de la protection en minutes
            beta_n (float): Coefficient de carbonisation fictif
            k_2 (float): Coefficient d'accelération de la carbonisation'"""
        if t_f:
            if t_ch < t_f:

                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_a = (25 - (t_f - t_ch) * k_2 * beta_n) / (2 * beta_n) + t_f  # en min/equa3.9
                    return t_a

            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    t_a = min(2 * t_f, 12.5 / beta_n + t_f)  # en min/equa3.8
                    return t_a

            return val()
        else:
            return ("", 0)

    def _get_section_reduction(self):
        self._base_b_calcul = self.b_calcul
        self._base_h_calcul = self.h_calcul
        self._def = {}
        # Index des lignes
        index = ["0", "tch", "tf", "ta", "t_expo", "d_ef"]
        # Colonnes multi-index : (orientation, type de valeur)
        cols = pd.MultiIndex.from_product(
            [("Haut", "Bas", "Gauche", "Droite"), ("Temps(min)", "d_char_n(mm)")],
            names=["Orientation", "Type"],
        )
        # DataFrame vide
        self.d_ef = pd.DataFrame(index=index, columns=cols)

        for orientation, protection in self.protection.items():
            if protection != "Pas d'exposition":
                d_char_n = 0
                d_char_n_values = {"0": 0, "tch": 0, "tf": 0, "ta": 0, "t_expo": 0}

                beta_0, beta_n = self._get_bar_beta0_and_betan()
                beta_0_protect, beta_n_protect = self._get_wood_protect_beta0_and_betan(
                    protection
                )
                latex_tch_tf, res_tch_tf = self._get_tch_and_tf(
                    orientation, beta_0_protect
                )
                tch = res_tch_tf[0]
                tf = res_tch_tf[1]
                k2 = self._get_k2(orientation)
                latex_ta, ta = self._get_ta(tch, tf, beta_n, k2)

                # Insertion des temps
                self.d_ef.loc["0", (orientation, "Temps(min)")] = 0
                self.d_ef.loc["tch", (orientation, "Temps(min)")] = tch
                self.d_ef.loc["tf", (orientation, "Temps(min)")] = tf
                self.d_ef.loc["ta", (orientation, "Temps(min)")] = ta
                self.d_ef.loc["t_expo", (orientation, "Temps(min)")] = self.t_expo

                # Calcul de d_char_n sur chaque phase
                if self.t_expo > tf:
                    if tch < tf:
                        d1 = self._d_char_n(tf - tch, beta_n * k2)
                        d_char_n = d_char_n + d1
                        d_char_n_values["tf"] = d_char_n

                    k3 = 2
                    if ta < self.t_expo:
                        d2 = self._d_char_n(ta - tf, beta_n * k3)
                        d_char_n = d_char_n + d2
                        d_char_n_values["ta"] = d_char_n

                        d3 = self._d_char_n(self.t_expo - ta, beta_n)
                        d_char_n = d_char_n + d3
                        d_char_n_values["t_expo"] = d_char_n
                    else:
                        d2 = self._d_char_n(ta - tf, beta_n * k3)
                        d_charn_sup_texpo = d_char_n + d2
                        d_char_n_values["ta"] = d_charn_sup_texpo

                        d3 = self._d_char_n(self.t_expo - tf, beta_n * k3)
                        d_char_n = d_char_n + d3
                    d_char_n_values["t_expo"] = d_char_n

                    # Stockage dans le DataFrame des d_char_n
                    for key, dcharn in d_char_n_values.items():
                        self.d_ef.loc[key, (orientation, "d_char_n(mm)")] = dcharn

                    # Calcul de d_ef final
                    d_ef = self._d_ef(d_char_n, self.t_expo)
                    self.d_ef.loc["d_ef", (orientation, "Temps(min)")] = self.t_expo
                    self.d_ef.loc["d_ef", (orientation, "d_char_n(mm)")] = d_ef[1]
                    self._def[orientation] = (latex_tch_tf + latex_ta + d_ef[0], d_ef[1])
            else:
                self._def[orientation] = ("\\[\text{Pas d'exposition}\\]", 0 * si.mm)
            
        self.b_calcul = self.b_calcul - self._def["Gauche"][1] - self._def["Droite"][1] # Réduction de la section
        self.h_calcul = self.h_calcul - self._def["Haut"][1] - self._def["Bas"][1] # Réduction de la section

    def get_def(self, orientation: str = ORIENTATION) -> tuple:
        """Retourne la profondeur de carbonisation effective en mm suivant l'orientation donnée"""
        return self._def[orientation]

    def _f_type_d(self, typeCarac=Barre.CARACTERISTIQUE[0:6]):
        """Méthode donnant la résistance de calcul de l'élément fonction de la vérification

        Args:
            typeCarac (str, optional): type de résistance caractéristique (flexion = "fm0k", compression = "fc0k" etc.). Defaults to "fm0k".
            loadtype (str, optional): Durée de chargement (Permanente, Court terme etc.). Defaults to "Permanente".
            typecombi (str, optional): type de combinaison étudiée ("Fondamentales" ou " Accidentelles"). Defaults to "Fondamentales".

        Returns:
            float: Résistance de calcul en N/mm2 du type de vérification étudié.
        """
        gamma_M_fi = self._get_gamma_M("Accidentelles")
        K_mod_fi = self.KMOD_FI  # en situation à chaud
        K_fi = self.K_fi
        f_type_k = float(self.caract_meca.loc[typeCarac]) * si.MPa

        if typeCarac == "fm0k" and self.k_sys > 1:
            k_sys = self.k_sys

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_type_Rd = k_sys * f_type_k * K_fi * K_mod_fi / gamma_M_fi
                return f_type_Rd

        else:

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
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
        dim = {"y": self._base_h_calcul.value * 10**3, "z": self._base_b_calcul.value * 10**3}

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
                warnings.warn("LVL non pris en compte dans cette fonction")
                kh[cle] = 1
        return kh

    def Emean_fin(self, psy_2: float) -> tuple:
        """Renvoie le E,mean,fin en fonction du Kdef et du psy2"""
        self.psy_2 = psy_2
        E0_mean = int(self.caract_meca.loc["E0mean"]) * si.MPa
        K_fi = self.K_fi
        K_def = self.K_def

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            E_mean_fin_fi = E0_mean * K_fi / (1 + psy_2 * K_def)
            return E_mean_fin_fi

        value = val()
        self.E_mean_fin = value[1]
        return value


# ================================ FLEXION ==================================
