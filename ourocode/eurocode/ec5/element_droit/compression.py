# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import matplotlib.pyplot as plt

import math as mt
from math import sqrt, pi, cos, sin, radians
import numpy as np

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.element_droit.barre import Barre


class Compression(Barre):
    """Classe de vérification des éléments bois en compression axiale selon l'EC5 §6.2 et §6.3.2.

    Effectue les calculs de résistance, élancement et flambement selon
    l'Eurocode 5 - Partie 1-1. Hérite de Barre pour les caractéristiques
    géométriques et mécaniques.

    Vérifie :
    - La résistance en compression axiale (§6.2.2, Eq. 6.2)
    - Le flambement avec coefficient kc (§6.3.2, Eq. 6.23-6.24)
    - Les combinaisons flexo-compression (si objet Flexion fourni)
    """

    COEF_LF = {"Encastré 1 côté" : 2,
                "Rotule - Rotule" : 1,
                "Encastré - Rotule" : 0.7,
                "Encastré - Encastré" : 0.5,
                "Encastré - Rouleau" : 1}

    def __init__(self, lo_y: si.mm, lo_z: si.mm, type_appuis: str=COEF_LF, *args, **kwargs):
        """Initialise un objet de vérification en compression axiale.

        Définit les longueurs de flambement et le coefficient de longueur efficace
        selon les conditions d'appui pour le calcul du flambement.

        Args:
            lo_y (si.mm): Longueur de flambement suivant l'axe y en mm.
                Mettre 0 si pas de risque de flambement selon cet axe.
            lo_z (si.mm): Longueur de flambement suivant l'axe z en mm.
                Mettre 0 si pas de risque de flambement selon cet axe.
            type_appuis (str): Type de conditions d'appui pour le coefficient β.
                Détermine la longueur efficace lf = β · lo.
                Valeurs possibles (voir COEF_LF):
                - "Encastré 1 côté" : β = 2.0 (console)
                - "Rotule - Rotule" : β = 1.0 (articulé-articulé)
                - "Encastré - Rotule" : β = 0.7
                - "Encastré - Encastré" : β = 0.5
                - "Encastré - Rouleau" : β = 1.0 (encastré-glissière)
                Defaults to "Rotule - Rotule".
            *args: Arguments positionnels transmis à Barre.
            **kwargs: Arguments nommés transmis à Barre (b, h, classe, etc.).

        Note:
            La longueur efficace de flambement lf est calculée par :
            lf = lo × coef_lef
        """
        super().__init__(*args, **kwargs)
        self.lo_comp = {"y":lo_y * si.mm, "z":lo_z * si.mm}
        self.lo_y = self.lo_comp['y']
        self.lo_z = self.lo_comp['z']
        self.type_appuis = type_appuis
        self.coef_lef = self.COEF_LF[type_appuis]
        self._Anet = self.aire

    @property
    def lamb(self) -> tuple:
        """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lo_y = self.lo_comp['y'].value * 10**3
        lo_z = self.lo_comp['z'].value * 10**3
        coef_lef = self.coef_lef
        I_y = self.inertie[0].value * 10**12
        I_z = self.inertie[1].value * 10**12
        A = self._Anet.value * 10**6

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            lamb_y = (lo_y * coef_lef) / sqrt(I_y / A)
            lamb_z = (lo_z * coef_lef) / sqrt(I_z / A)
            return {'y': lamb_y, 'z': lamb_z}
        return val()
    
    
    @property
    def lamb_rel_Axe(self) -> tuple:
        """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb_y = self.lamb[1]['y']
        lamb_z = self.lamb[1]['z']
        f_c0k = float(self.caract_meca.loc['fc0k']) * si.MPa
        E_0_05 = int(self.caract_meca.loc['E005']) * si.MPa

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            lamb_rel_y = (lamb_y / pi) * sqrt(f_c0k / E_0_05)
            lamb_rel_z = (lamb_z / pi) * sqrt(f_c0k / E_0_05)
            return {'y': lamb_rel_y, 'z': lamb_rel_z}
        return val()
    

    @property
    def beta_C(self) -> float:
        if self.type_bois == 'Massif':
            betaC = 0.2
        else:
            betaC = 0.1
        return betaC
    
    
    @property
    def k_Axe(self):
        """ Retourne le facteur Ky ou Kz (fonction de l'axe de flambement) """
        beta_C = self.beta_C
        lamb_rel_y = self.lamb_rel_Axe[1]['y']
        lamb_rel_z = self.lamb_rel_Axe[1]['z']

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            k_y = 0.5 * (1 + beta_C * (lamb_rel_y - 0.3) + lamb_rel_y ** 2)
            k_z = 0.5 * (1 + beta_C * (lamb_rel_z - 0.3) + lamb_rel_z ** 2)
            return {'y': k_y, 'z': k_z}
        return val()
    

    @property
    def kc_Axe(self) -> tuple:
        """ Retourne le coefficient multiplicateur KcAxe  (axe = y ou z suivant axe de rotation en flambement) de fc,0,d """
        k_y = self.k_Axe[1]['y']
        k_z = self.k_Axe[1]['z']
        lamb_rel_y = self.lamb_rel_Axe[1]['y']
        lamb_rel_z = self.lamb_rel_Axe[1]['z']

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            k_c_y = 1 / (k_y + sqrt(k_y** 2 - lamb_rel_y ** 2))
            k_c_z = 1 / (k_z + sqrt(k_z** 2 - lamb_rel_z ** 2))
            return {'y': min(k_c_y, 1), 'z': min(k_c_z, 1)}
        return val()


    def f_c_0_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul en compression axiale f_c,0,d selon l'EC5 §6.2.2.

        Détermine la résistance à partir de la résistance caractéristique fc,0,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul fc,0,d en MPa avec unité (si.MPa).

        Note:
            Cette valeur est réduite par le coefficient kc en cas de flambement.
        """
        return super()._f_type_d("fc0k", loadtype, typecombi)
    
    
    def sigma_c_0_d(self, Fc0d: si.kN, Anet: si.mm**2=None) -> tuple:
        """Calcule la contrainte de compression axiale sigma_c,0,d selon l'EC5 §6.2.2.

        Détermine la contrainte normale due à l'effort de compression axial.
        Prend en compte une section nette réduite si spécifiée.

        Args:
            Fc0d (si.kN): Effort de compression axial en kN.
            Anet (si.mm**2, optional): Aire nette de la section en mm² si réduction
                (entailles, perçages). Doit être ≤ aire brute. Defaults to None.

        Returns:
            tuple: (latex_string, valeur) où valeur est sigma_c,0,d en MPa avec unité.

        Raises:
            ValueError: Si Anet > aire brute de la section.

        Note:
            La valeur est stockée dans l'attribut sigma_c_0_rd.
        """
        self.Fc_0_d = Fc0d * si.kN
        Fc_0_d = self.Fc_0_d
        if Anet and Anet * si.mm**2 <= self.aire:
            self._Anet = Anet * si.mm**2
        A = self._Anet

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_c_0_d = Fc_0_d / A
            return sigma_c_0_d
        value = val()
        self.sigma_c_0_rd = value[1]
        return value
    
    
    def taux_c_0_d(self, flexion: object=None) -> tuple:
        """Calcule les taux de travail en compression axiale selon l'EC5 §6.2.2 et §6.3.2.

        Vérifie les critères de résistance :
        - Equ. 6.2 : Compression simple (sigma_c,0,d / f_c,0,d)
        - Equ. 6.23-6.24 : Flambement (sigma_c,0,d / (kc · f_c,0,d))
        - Equ. 6.19-6.20 : Flexo-compression (si objet Flexion fourni)

        Args:
            flexion (Flexion, optional): Objet Flexion déjà calculé
                pour les combinaisons flexo-compression. Defaults to None.

        Returns:
            tuple: (latex_string, taux_dict) où taux_dict contient :
                - "equ6.2" : Compression simple sans flambement
                - "equ6.23", "equ6.24" : Compression avec flambement selon y et z
                - "equ6.19", "equ6.20" : Flexo-compression (si flexion fournie)
                Valeurs en pourcentage (0.85 = 85%).

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_c_0_rd = {}
        sigma_c_0_d = self.sigma_c_0_rd
        f_c_0_d = self.f_type_rd
        lamb_rel_y  = self.lamb_rel_Axe[1]['y']
        lamb_rel_z  = self.lamb_rel_Axe[1]['z']
        K_c_y = self.kc_Axe[1]['y']
        K_c_z = self.kc_Axe[1]['z']

        from ourocode.eurocode.ec5.element_droit.flexion import Flexion as _Flexion
        if flexion and isinstance(flexion, _Flexion):
            taux_6_11 = flexion.taux_m_rd["equ6.11"]
            taux_6_12 = flexion.taux_m_rd["equ6.12"]
        else:
            taux_6_11 = 0
            taux_6_12 = 0
        
        if lamb_rel_y <= 0.3 and lamb_rel_z <= 0.3:
            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_6_2 = sigma_c_0_d / f_c_0_d # equ6.2
                taux_6_19 = (sigma_c_0_d / (f_c_0_d * K_c_y))**2 + taux_6_11 # equ6.19
                taux_6_20 = (sigma_c_0_d / (f_c_0_d * K_c_z))**2 + taux_6_12 # equ6.20
                return taux_6_2, taux_6_19, taux_6_20
            value = val()
            self.taux_c_0_rd['equ6.19'] = value[1][1]
            self.taux_c_0_rd['equ6.20'] = value[1][2]
        else:      
            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_6_2 = sigma_c_0_d / f_c_0_d # equ6.2
                taux_6_23 = sigma_c_0_d / (f_c_0_d * K_c_y) + taux_6_11 # equ6.23
                taux_6_24 = sigma_c_0_d / (f_c_0_d * K_c_z) + taux_6_12 # equ6.24
                return taux_6_2, taux_6_23, taux_6_24
            value = val()
            self.taux_c_0_rd['equ6.23'] = value[1][1]
            self.taux_c_0_rd['equ6.24'] = value[1][2]

        self.taux_c_0_rd['equ6.2'] = value[1][0]
        max_taux = max([v for v in self.taux_c_0_rd.values()])
        synthese = [
            ["Compression bois", max_taux, None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value

class Compression_perpendiculaire(Barre):
    """Classe de vérification des éléments bois en compression perpendiculaire selon l'EC5 §6.1.5.

    Effectue les calculs de résistance et de contrainte en compression perpendiculaire
    au fil du bois (appuis de poutres, abouts de pieux, etc.) selon l'Eurocode 5.
    Hérite de Barre pour les caractéristiques géométriques et mécaniques.

    La vérification utilise l'équation 6.3 avec le coefficient K_c,90 :
    σ_c,90,d / (K_c,90 · f_c,90,d) ≤ 1
    """

    TYPE_APPUIS = ("Appuis discret", "Appuis continu")
    def __init__(
        self, 
        b_appuis:si.mm, 
        l_appuis: si.mm, 
        l1d: si.mm=10000, 
        l1g: si.mm=10000, 
        ad: si.mm=0, 
        ag: si.mm=0, 
        type_appuis_90: str=TYPE_APPUIS, 
        *args, 
        **kwargs
        ):
        """Initialise un objet de vérification en compression perpendiculaire.

        Args:
            b_appuis (si.mm): Largeur d'appuis en mm.
            l_appuis (si.mm): Longueur de l'appuis en mm.
            l1d (si.mm, optional): Distance entre les charges en mm (l et l) (si pas de l1d ne rien mettre). Defaults to 10000.
            l1g (si.mm, optional): Distance entre les charges en mm (l et l) (si pas de l1g ne rien mettre). Defaults to 10000.
            ad (si.mm, optional): Distance depuis le bord jusqu'à l'appuis à droite (l) en mm (si pas de ad et au bord ne rien mettre). Defaults to 0.
            ag (si.mm, optional): Distance depuis le bord jusqu'à l'appuis à gauche (l) en mm (si pas de ad et au bord ne rien mettre). Defaults to 0.
            type_appuis_90 (str, optional): Type d'appuis (Appui continu, Appui discret). Defaults to TYPE_APPUIS.
            *args: Arguments transmis à Barre (classe, etc.).
            **kwargs: Arguments transmis à Barre (b, h, etc.).

        Note:
            La longueur efficace de compression perpendiculaire l_ef est calculée
            en fonction de la configuration des appuis et des distances entre charges.
        """
        super().__init__(*args, **kwargs)
        self.b_appuis = b_appuis * si.mm
        self.l_appuis = l_appuis * si.mm
        self.l1d = l1d * si.mm
        self.l1g = l1g * si.mm
        self.ad = ad * si.mm
        self.ag = ag * si.mm
        self.type_appuis_90 = type_appuis_90

    @property
    def K_c90(self):
        """Retourne le facteur K_c,90 qui tient compte de la configuration de chargement, du fendage et de la déformation
            en compression avec pour argument :
            h : Hauteur de l'élement subissant la compression en mm
            lO : Longeur de l'appuis en compression en mm
            l1 : Distance la plus petite entre deux appuis en mm (l et l)
        """
        try:
            return self._setter_K_c90
        except AttributeError:
            if self.l1d == 0 and self.l1g > 0:
                l1 = self.l1g
                self.ag = self.l1g
            elif self.l1g == 0 and self.l1d > 0:
                l1 = self.l1d
                self.ad = self.l1d
            else:
                l1 = min(self.l1d, self.l1g)

            if self.type_appuis_90 == "Appuis discret":
                if self.type_bois == 'Massif':
                    if self.h_calcul.value * 10**3 <= 300:
                        if l1 >= 2 * self.h_calcul:
                            kc_90 = 1.5
                        else:
                            kc_90 = 1
                    else:
                        kc_90 = 1.5
                elif self.type_bois == 'BLC':
                    if self.h_calcul.value * 10**3 <= 300 and self.l_appuis.value * 10**3 <= 400:
                        if l1 >= 2 * self.h_calcul:
                            kc_90 = 1.75
                        else:
                            kc_90 = 1
                    elif self.h_calcul.value * 10**3 <= 300 and self.l_appuis.value * 10**3 > 400:
                        kc_90 = 1
                    else:
                        kc_90 = 1.75
                else:
                    if self.h_calcul.value * 10**3 > 300:
                        kc_90 = 1.75
                    else:
                        kc_90 = 1
            else:
                if self.type_bois == 'Massif':
                    if self.h_calcul.value * 10**3 <= 300:
                        if l1 >= 2 * self.h_calcul:
                            kc_90 = 1.25
                        else:
                            kc_90 = 1
                    else:
                        kc_90 = 1.5
                elif self.type_bois == 'BLC':
                    if self.h_calcul.value * 10**3 <= 300:
                        if l1 >= 2 * self.h_calcul:
                            kc_90 = 1.5
                        else:
                            kc_90 = 1
                    else:
                        kc_90 = 1.75
                else:
                    if self.h_calcul.value * 10**3 > 300:
                        kc_90 = 1.75
                    else:
                        kc_90 = 1
            return kc_90

    @K_c90.setter
    def K_c90(self, value):
        self._setter_K_c90 = value


    def f_c_90_d(self, loadtype: str=Barre.LOAD_TIME, typecombi: str=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul en compression perpendiculaire f_c,90,d selon l'EC5 §6.1.5.

        Détermine la résistance à partir de la résistance caractéristique fc,90,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul fc,90,d en MPa avec unité (si.MPa).
        """
        return super()._f_type_d("fc90k", loadtype, typecombi)
    
    
    def sigma_c_90_d(self, Fc90d: si.kN) -> tuple:
        """Calcule la contrainte de compression perpendiculaire sigma_c,90,d selon l'EC5 §6.1.5.

        Détermine la contrainte en tenant compte de l'aire effective d'appui,
        qui inclut une diffusion des efforts sur 30 mm ou jusqu'aux bords/entraxe.

        Args:
            Fc90d (si.kN): Effort de compression perpendiculaire en kN.

        Returns:
            tuple: (latex_string, valeur) où valeur est sigma_c,90,d en MPa avec unité.

        Note:
            L'aire effective a_ef = (l_appuis + min(30mm, distance_bord, l_appuis, 0.5×entraxe)) × b_appuis
            La valeur est stockée dans l'attribut sigma_c_90_rd.
        """
        
        self.Fc90d = Fc90d * si.kN
        Fc_90_d = self.Fc90d
        l_appuis = self.l_appuis
        b_appuis = self.b_appuis
        l_1d = self.l1d
        l_1g = self.l1g
        ad = self.ad
        ag = self.ag
        mm = si.mm

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            a_ef = (l_appuis + min(30*mm, ad, l_appuis, 0.5*l_1d) + min(30*mm, ag, l_appuis, 0.5*l_1g)) * b_appuis
            sigma_c_90_d = Fc_90_d / a_ef
            return sigma_c_90_d
        value = val()
        self.sigma_c_90_rd = value[1]
        return value


    def taux_c_90_d(self) -> tuple:
        """Calcule le taux de travail en compression perpendiculaire selon l'EC5 §6.1.5 (Eq. 6.3).

        Vérifie le critère : σ_c,90,d / (K_c,90 · f_c,90,d) ≤ 1

        Le coefficient K_c,90 (déterminé par la propriété K_c90) tient compte :
        - De la configuration des appuis (discrets ou continus)
        - Du type de bois (massif, BLC, etc.)
        - De la hauteur de l'élément et des distances entre appuis

        Returns:
            tuple: (latex_string, valeur) où valeur est le taux en pourcentage
                (0.75 = 75%). Stocké dans taux_c_90_rd['equ6.3'].

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_c_90_rd = {}
        sigma_c_90_d = self.sigma_c_90_rd
        K_c90 = self.K_c90
        f_c_90_d = self.f_type_rd

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_6_3 = sigma_c_90_d / (K_c90 * f_c_90_d) # equ6.3
            return taux_6_3
        
        value = val()
        self.taux_c_90_rd['equ6.3'] = value[1]
        synthese = [
            ["Compression perpendiculaire bois", self.taux_c_90_rd['equ6.3'], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value
    
    
    def show_c90(self) -> None:
        """Affiche l'image des caractéristiques d'une compression perpendiculaire
        """
        self._show_element("C90_def.png")



class Compression_inclinees(Compression_perpendiculaire):
    """Classe de vérification des éléments bois en compression inclinée selon l'EC5 §6.2.2.

    Effectue les calculs de résistance et de contrainte en compression inclinée
    par rapport au fil du bois selon l'Eurocode 5 - Partie 1-1, article 6.2.2.
    Hérite de Compression_perpendiculaire pour la gestion des appuis et des coefficients.

    La vérification utilise l'équation 6.16 avec la formule de Hankinson :
    σ_c,α,d ≤ f_c,0,d / [(f_c,0,d/(K_c,90·f_c,90,d))·sin²(α) + cos²(α)]
    """

    def __init__(self, alpha: float=45, **kwargs):
        """Initialise un objet de vérification en compression inclinée.

        Args:
            alpha (float, optional): Angle d'inclinaison de la compression par rapport
                au fil du bois en degrés. 0° = compression axiale, 90° = compression
                perpendiculaire. Defaults to 45°.
            **kwargs: Arguments transmis à Compression_perpendiculaire
                (b_appuis, l_appuis, etc.).
        """
        super().__init__(**kwargs)
        self.alpha = alpha
    
    def sigma_c_alpha_d(self, Fcad: si.kN) -> tuple:
        """Calcule la contrainte de compression inclinée sigma_c,alpha,d.

        Détermine la contrainte normale due à l'effort de compression inclinée,
        en utilisant l'aire brute d'appui (sans diffusion).

        Args:
            Fcad (si.kN): Effort de compression inclinée en kN.

        Returns:
            tuple: (latex_string, valeur) où valeur est sigma_c,alpha,d en MPa avec unité.

        Note:
            La valeur est stockée dans l'attribut sigma_c_alpha_rd.
        """
        b_appuis = self.b_appuis
        l_appuis = self.l_appuis
        self.Fc_alpha_d = Fcad * si.kN
        Fc_alpha_d  = self.Fc_alpha_d
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            A = b_appuis * l_appuis
            sigma_c_alpha_d = Fc_alpha_d / A
            return sigma_c_alpha_d
        
        value = val()
        self.sigma_c_alpha_rd = value[1]
        return value
    

    def taux_c_alpha_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION) -> tuple:
        """Calcule le taux de travail en compression inclinée selon l'EC5 §6.2.2 (Eq. 6.16).

        Vérifie le critère de Hankinson : σ_c,α,d ≤ f_c,α,d
        où f_c,α,d = f_c,0,d / [(f_c,0,d/(K_c,90·f_c,90,d))·sin²(α) + cos²(α)]

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            tuple: (latex_string, valeur) où valeur est le taux en pourcentage
                (0.75 = 75%). Stocké dans taux_c_alpha_rd['equ6.16'].

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_c_alpha_rd = {}
        f_c_0_d = self._f_type_d("fc0k", loadtype, typecombi)[1]
        f_c_90_d = self.f_c_90_d(loadtype, typecombi)[1]
        alpha = self.alpha
        K_c90 = self.K_c90
        sigma_c_alpha_d = self.sigma_c_alpha_rd

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_6_16 = sigma_c_alpha_d / (f_c_0_d / ((f_c_0_d /(K_c90 * f_c_90_d)) * sin(radians(alpha))**2 + cos(radians(alpha))**2)) # equ6.16
            return taux_6_16
         
        value = val()
        self.taux_c_alpha_rd['equ6.16'] = value[1]
        synthese = [
            ["Compression inclinée bois", self.taux_c_alpha_rd['equ6.16'], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value
    


