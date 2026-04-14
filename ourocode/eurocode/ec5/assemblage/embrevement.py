# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.assemblage.assemblage import Assemblage
from ourocode.eurocode.ec5.element_droit.barre import Barre
from ourocode.eurocode.ec5.element_droit.compression import Compression_perpendiculaire, Compression_inclinees

class Embrevement(Assemblage):
    TYPE_EMBREVEMENT = ("Bissectrice", "Equerre à la pièce 2", "Equerre à la pièce 1")
    def __init__(
        self, 
        alpha: float=45, 
        type_embrevement: str=TYPE_EMBREVEMENT, 
        prof_embrevement: si.mm=0, 
        l_talon: si.mm=500,
        l1: si.mm=10000,
        **kwargs):
        """Créer une classe Embrevement hérité de la classe Assemblage du module EC5_Assemblage.py.
        
        Cette classe permet de vérifier un embrèvement simple avec comme i = 1 la pièce avec l'embrèvement en bout et comme i = 2 la pièce entaillée.
        Attention cette classe ne vérifie pas la section réduite de la pièce 2 !
        
        Args:
            alpha (float): angle en degrés entre le fil de la pièce 1 et le fil de la pièce 2 (Ne peux pas être supérieur à 90°). Defaults to 45.
            type_embrevement (str): type d'embrèvement "Bissectrice", "Equerre à la pièce 2", "Equerre à la pièce 1". Defaults to "Bissectrice".
            prof_embrevement (si.mm): profondeur de l'embrèvement en mm parallèle au fil de la pièce 2. Si ce dernier est égale à 0,
                alors la profondeur théorique est calculé automatiquement suivant les règles ci-dessous:
                    - 25°<= alpha < 50° : tv <= h2/4 
                    - 50°<= alpha < 60° : tv <= h2*(2/3 - alpha/120)
                    - alpha >= 60° : tv <= h2/6
                    avec h2 = hauteur de la pièce 2
            l_talon (si.mm): longueur du talon en mm sur la pièce 2. Defaults to 500.
            l1 (si.mm): distance de la face l'appuis à l'embrèvement sur la pièce 2 (l et l) (si pas de l1 ne rien mettre). Defaults to 10000.
        """
        super().__init__(**kwargs)
        self.alpha = alpha
        self.type_embrevement = type_embrevement
        self.prof_embrevement = prof_embrevement * si.mm
        self.l_talon = l_talon * si.mm
        self.l1 = l1 * si.mm
        self.dict_taux_embrevement = {}

        if self.prof_embrevement == 0*si.mm:
            h2 = self.beam_2.h_calcul
            def prof_embrevement():
                if 25 <= self.alpha <= 50:
                    prof_embrevement = h2/4
                elif 50 < self.alpha <= 60:
                    prof_embrevement = h2*(2/3 - self.alpha/120)
                elif self.alpha > 60:
                    prof_embrevement = h2/6
                return round(prof_embrevement, 0)
            self.prof_embrevement = prof_embrevement()
    
    def _compression_about(self, N_c1_d: si.kN, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Vérifie la compression d'about de l'embrevement"""
        if self.type_embrevement == "Bissectrice":
            alpha_about = self.alpha / 2
            tv = self.prof_embrevement / cos(radians(alpha_about))
        elif self.type_embrevement == "Equerre à la pièce 1":
            alpha_about = 0
            tv = self.prof_embrevement / cos(radians(self.alpha))
        elif self.type_embrevement == "Equerre à la pièce 2":
            alpha_about = self.alpha
            tv = self.prof_embrevement
        
        Fcad = N_c1_d * cos(radians(alpha_about))

        c_alpha = Compression_inclinees._from_parent_class(
            self.beam_1, 
            b_appuis=self.beam_1.b_calcul.value*10**3, 
            l_appuis=tv.value*10**3, 
            l1d=10000, 
            l1g=10000, 
            ad=10000, 
            ag=10000, 
            type_appuis_90="Appuis discret", 
            alpha=alpha_about)
        c_alpha.K_c90 = 1 # On bride le Kc90 à 1 de manière sécuritaire

        sigma_c_alpha_d = c_alpha.sigma_c_alpha_d(Fcad)
        taux_c_alpha_d = c_alpha.taux_c_alpha_d(loadtype, typecombi)
        latex = (sigma_c_alpha_d[0] + taux_c_alpha_d[0])
        self.dict_taux_embrevement["taux compression about"] = taux_c_alpha_d[1]
        return (latex, taux_c_alpha_d[1])

    def _compression_transversale(self, N_c1_d: si.kN, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Vérifie la compression transversale de l'embrevement"""
        l_appuis = self.beam_1.h_calcul / sin(radians(self.alpha))
        c_90 = Compression_perpendiculaire._from_parent_class(
            self.beam_2, 
            b_appuis=self.beam_1.b_calcul.value*10**3, 
            l_appuis=l_appuis.value*10**3, 
            l1d=self.l1.value*10**3, 
            l1g=10000, 
            ad=10000, 
            ag=10000, 
            type_appuis_90="Appuis discret")
        
        Fc90d = N_c1_d * sin(radians(self.alpha))
        f_c_90_d = c_90.f_c_90_d(loadtype, typecombi)
        sigma_c_90_d = c_90.sigma_c_90_d(Fc90d)
        taux_c_90_d = c_90.taux_c_90_d()
        latex = (f_c_90_d[0] + sigma_c_90_d[0] + taux_c_90_d[0])
        self.dict_taux_embrevement["taux compression transversale"] = taux_c_90_d[1]
        return (latex, taux_c_90_d[1])

    def _cisaillement_talon(self, N_c1_d: si.kN, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Vérifie le cisaillement du talon de la pièce 2, le talon est limité pour prendre en compte l'effet ciseau à bois"""
        N_c1_d = N_c1_d * si.kN
        l_talon = self.l_talon
        prof_embr = self.prof_embrevement
        b_embr = self.beam_1.b_calcul
        k_cr = 0.67 # On bride le Kcr à 0.67 de manière sécuritaire
        f_v_d = self.beam_2._f_type_d("fvk", loadtype, typecombi)[1]
        alpha = self.alpha

        if self.beam_1.b_calcul < self.beam_2.b_calcul:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                l_v_ef = min(l_talon, prof_embr * 8)
                L_v = (b_embr + 2 * prof_embr) * k_cr # périmètre transversal cisaillé
                N_c2_d = N_c1_d * cos(radians(alpha))
                taux_talon = N_c2_d / (f_v_d * L_v * l_v_ef)
                return taux_talon
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                l_v_ef = min(l_talon, prof_embr * 8)
                L_v = b_embr * k_cr # périmètre transversal cisaillé
                N_c2_d = N_c1_d * cos(radians(alpha))
                taux_talon = N_c2_d / (f_v_d * L_v * l_v_ef)
                return taux_talon
        value = val()
        self.dict_taux_embrevement["taux cisaillement talon"] = value[1]
        return value

    def taux_embrevement(self, N_c1_d: si.kN, loadtype: str=Barre.LOAD_TIME, typecombi: str=Barre.TYPE_ACTION) -> tuple:
        """Calcule les taux de travail de l'embrèvement simple selon EN 1995-1-1 §8.8.

        Effectue trois vérifications :

        1. Compression d'about de la pièce 1 (angle alpha_about selon type_embrevement).
        2. Compression transversale au fil de la pièce 2 (90° au fil).
        3. Cisaillement du talon de la pièce 2 (longueur efficace ≤ min(l_talon, 8×t_v)).

        Args:
            N_c1_d (float): Effort normal de calcul dans la pièce 1 (force diagonale) en kN.
            loadtype (str): Classe de durée de chargement selon EN 1995-1-1 §2.3.1.2.
                Defaults to "Permanente".
            typecombi (str): Type de combinaison (fondamentale ou accidentelle). Defaults to fondamentale.

        Returns:
            tuple: (latex_string, dict_taux_embrevement) où dict_taux_embrevement contient :
                - "taux compression about" : taux de compression inclinée de la pièce 1.
                - "taux compression transversale" : taux de compression perpendiculaire de la pièce 2.
                - "taux cisaillement talon" : taux de cisaillement du talon de la pièce 2.
        """
        about =self._compression_about(N_c1_d, loadtype, typecombi)
        transversale = self._compression_transversale(N_c1_d, loadtype, typecombi)
        talon = self._cisaillement_talon(N_c1_d, loadtype, typecombi)
        synthese = [
            ["Embrèvement avec compression d'about de la pièce 1", about[1], None],
            ["Embrèvement avec compression transversale de la pièce 2", transversale[1], None],
            ["Embrèvement avec cisaillement du talon de la pièce 2", talon[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return (about[0] + transversale[0] + talon[0], self.dict_taux_embrevement)



