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
from ourocode.eurocode.ec5.element_droit.compression import Compression
from ourocode.eurocode.ec5.element_droit.traction import Traction


class Flexion(Barre):
    """Classe de vérification à la flexion selon l'EN 1995-1-1 §6.1.6, §6.2.3, §6.2.4 et §6.3.3.

    Cette classe effectue les vérifications de résistance à la flexion et
    de stabilité au déversement (flambement latéral) pour des poutres en bois.

    Elle hérite de la classe Barre pour récupérer les caractéristiques
    géométriques et mécaniques, et utilise le pattern _from_parent_class
    pour l'enchaînement des vérifications.
    """

    COEF_LEF = {"Appuis simple" : [1, 0.9, 0.8], 
                "Porte à faux": [0.5, 0.8]}
    LOAD_POS = (
        "Charge sur fibre comprimée",
        "Charge sur fibre neutre",
        "Charge sur fibre tendue"
        )

    def __init__(self, 
        lo_rel_y:si.mm, 
        lo_rel_z:si.mm, 
        coeflef_y: float=0.9, 
        coeflef_z: float=0.9, 
        pos: str=LOAD_POS, 
        *args, **kwargs):
        """Initialise une vérification en flexion avec paramètres de déversement.

        Args:
            lo_rel_y (si.mm): Longueur de déversement effective autour de l'axe Y
                (entre appuis latéraux), en millimètres.
            lo_rel_z (si.mm): Longueur de déversement effective autour de l'axe Z,
                en millimètres.
            coeflef_y (float, optional): Coefficient de longueur efficace selon l'EC5.
                - Appuis simple : 1.0 (moment constant), 0.9 (charge répartie),
                  0.8 (charge concentrée centrale)
                - Porte-à-faux : 0.5 (charge répartie), 0.8 (charge concentrée bout)
                Defaults to 0.9.
            coeflef_z (float, optional): Idem pour l'axe Z. Defaults to 0.9.
            pos (str, optional): Position de la charge verticale sur la hauteur.
                "Charge sur fibre comprimée": charge au-dessus de l'axe neutre
                    (aggrave le déversement, +2h sur l_ef)
                "Charge sur fibre neutre": charge au centre de gravité
                "Charge sur fibre tendue": charge en dessous de l'axe neutre
                    (favorise la stabilité, -0.5h sur l_ef)
            *args: Arguments transmis à la classe parent Barre.
            **kwargs: Arguments nommés transmis à Barre (b, h, classe, etc.).

        Note:
            La longueur efficace de déversement l_ef est calculée par :
            l_ef = lo_rel × coeflef (+ correction selon pos)
        """
        super().__init__(*args, **kwargs)
        self.lo_rel_y = lo_rel_y* si.mm
        self.lo_rel_z = lo_rel_z* si.mm
        self.lo_rel = {"y": lo_rel_y, "z": lo_rel_z}
        self.coeflef_y = coeflef_y
        self.coeflef_z = coeflef_z
        self.coeflef = {"y": coeflef_y, "z": coeflef_z}
        self.pos = pos

    @property
    def K_h(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        return self._K_h()

    @property
    def K_m(self):
        """Coefficient de distribution des contraintes K_m selon l'EC5 §6.1.6.

        Ce coefficient réduit la contrainte de flexion calculée pour les
        sections rectangulaires en bois massif, BLC ou LVL afin de tenir
        compte de la redistribution plastique des contraintes.

        Returns:
            float: Valeur de K_m.
                - 0.7 pour les sections rectangulaires en bois massif, BLC, LVL
                - 1.0 pour les sections circulaires ou les panneaux dérivés
        """
        if self.type_bois == "Massif" or self.type_bois == "BLC" or self.type_bois == "LVL":
            if self.section == "Rectangulaire":
                km = 0.7
            else:
                km = 1
        else:
            km = 1
        return km

    @property
    def sigma_m_crit(self):
        """Contrainte critique de déversement sigma_m,crit selon l'EC5 §6.3.3.

        Calculée par la formule de l'EC5 : sigma_m,crit = (0.78 × b² × E_0,05) / (h × l_ef)

        Cette contrainte caractérise la stabilité latérale de la poutre.
        Elle est corrigée en fonction de la position de la charge (pos).

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {'y': ..., 'z': ...}
                avec les contraintes critiques pour chaque direction.
        """
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
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_m_crit_y = (0.78 * b_calcul ** 2 * E_0_05) / (h_calcul * l_ef_y)
            sigma_m_crit_z = (0.78 * h_calcul ** 2 * E_0_05) / (b_calcul * l_ef_z)
            return {"y": sigma_m_crit_y, "z": sigma_m_crit_z}
        return val()

    @property
    def lamb_rel_m(self):
        """Élancement relatif en flexion lambda_rel,m selon l'EC5 §6.3.3.

        Rapport entre la résistance caractéristique et la contrainte critique :
        lambda_rel,m = sqrt(f_m,k / sigma_m,crit)

        Cet élancement caractérise le risque de déversement :
        - lambda_rel,m <= 0.75 : pas de risque de déversement (K_crit = 1)
        - 0.75 < lambda_rel,m <= 1.4 : zone de transition
        - lambda_rel,m > 1.4 : risque élevé de déversement

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {'y': ..., 'z': ...}
                avec les élancements relatifs pour chaque direction.
        """
        f_m0k = float(self.caract_meca.loc['fm0k']) *si.MPa
        sigma_m_crit_y = self.sigma_m_crit[1]['y']
        sigma_m_crit_z = self.sigma_m_crit[1]['z']

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            lamb_rel_m_y = sqrt(f_m0k / sigma_m_crit_y)
            lamb_rel_m_z = sqrt(f_m0k / sigma_m_crit_z)
            return {"y": lamb_rel_m_y, "z": lamb_rel_m_z}
        return val()

    @property
    def K_crit(self):
        """Coefficient de déversement K_crit selon l'EC5 §6.3.3.

        Ce coefficient minore la résistance à la flexion pour tenir compte
du risque de déversement latéral. Il dépend de l'élancement relatif :

        - lambda_rel,m <= 0.75 : K_crit = 1 (pas de déversement)
        - 0.75 < lambda_rel,m <= 1.4 : K_crit = 1.56 - 0.75 × lambda_rel,m
        - lambda_rel,m > 1.4 : K_crit = 1 / lambda_rel,m²

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {'y': ..., 'z': ...}
                avec les coefficients de déversement pour chaque direction.

        Note:
            La vérification finale utilise : sigma_m,d <= K_crit × f_m,d
        """
        lamb_rel_m_y = self.lamb_rel_m[1]['y']
        lamb_rel_m_z = self.lamb_rel_m[1]['z']
        result = [None, {"y": None, "z": None}]

        for axe in ["y", "z"]:
            lamb_rel_m = lamb_rel_m_y if axe == "y" else lamb_rel_m_z
            if lamb_rel_m <= 0.75:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_crit = 1
                    axe
                    return K_crit
            elif 0.75 < lamb_rel_m <= 1.4:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_crit = 1.56 - 0.75 * lamb_rel_m
                    axe
                    return K_crit
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_crit = 1 / (lamb_rel_m ** 2)
                    axe
                    return K_crit
            kcrit_axe = val()
            if result[0]:
                result[0] = result[0] + kcrit_axe[0]
            else:
                result[0] = kcrit_axe[0]
            result[1][axe] = kcrit_axe[1]
        result = (result[0], result[1])
        return result
    
    def f_m_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul en flexion f_m,d selon l'EC5 §6.1.6.

        La résistance est déterminée à partir de la résistance caractéristique fm,0,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement (permanent, long terme, etc.).
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul fm,d en MPa avec unité (si.MPa).
        """
        return self._f_type_d("fm0k", loadtype, typecombi)
    
    
    def sigma_m_d(self, My: si.kN*si.m, Mz: si.kN*si.m):
        """Calcule les contraintes de flexion sigma_m,d selon l'EC5 §6.1.6.

        Détermine les contraintes normales dues aux moments fléchissants My et Mz
        en utilisant la formule de Navier : σ = M·y/I

        Args:
            My (si.kN*m): Moment fléchissant autour de l'axe y (moment vertical)
                en kN·m. Mettre 0 si pas de flexion selon cet axe.
            Mz (si.kN*m): Moment fléchissant autour de l'axe z (moment horizontal)
                en kN·m. Mettre 0 si pas de flexion selon cet axe.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dictionnaire :
                {"y": sigma_my_d, "z": sigma_mz_d} en MPa avec unités.

        Note:
            Les valeurs sont stockées dans l'attribut sigma_m_rd.
            Pour une section rectangulaire : sigma = M·h/(2·I) = 6·M/(b·h²)
        """
        self.Md = {'y': My * si.kN*si.m, 'z': Mz * si.kN*si.m}
        self.sigma_m_rd = {'y': 0 * si.MPa, 'z': 0 * si.MPa}
        
        Iy = self.inertie[0]
        h_calcul = self.h_calcul
        Iz = self.inertie[1]
        b_calcul = self.b_calcul
        
        M_y = self.Md['y']
        M_z = self.Md['z']
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_my_d = M_y * h_calcul / (Iy * 2)
            sigma_mz_d = M_z * b_calcul / (Iz * 2)
            return {"y": sigma_my_d, "z": sigma_mz_d}
        value = val()
        self.sigma_m_rd["y"] = value[1]["y"]
        self.sigma_m_rd["z"] = value[1]["z"]
        return value
    

    def taux_m_d(self, compression: object=None, traction: object=None):
        """Calcule les taux de travail en flexion selon l'EC5 §6.1.6, §6.2.3 et §6.3.3.

        Vérifie les critères de résistance en flexion pure, flexion déviée,
        flexo-compression et flexo-traction selon les équations :
        - 6.11 et 6.12 : Flexion déviée avec K_m (facteur de distribution)
        - 6.33 : Flexion avec déversement (K_crit)
        - 6.17-6.20 : Combinaisons flexion + traction/compression
        - 6.35 : Flexo-compression avec risque de déversement

        Args:
            compression (Compression, optional): Objet Compression déjà calculé
                pour les combinaisons flexo-compression. Defaults to None.
            traction (Traction, optional): Objet Traction déjà calculé
                pour les combinaisons flexo-traction. Defaults to None.

        Returns:
            tuple: (latex_string, taux_dict) où taux_dict contient :
                - "equ6.11", "equ6.12" : Flexion déviée
                - "equ6.33y", "equ6.33z" : Flexion avec déversement
                - "equ6.17", "equ6.18" : Flexion + traction (si traction fournie)
                - "equ6.19", "equ6.20" : Flexion + compression (si compression fournie)
                - "equ6.23-6.35" : Combinaisons avancées (si compression fournie)
                Valeurs en pourcentage (0.85 = 85%).

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_m_rd = {}

        sigma_my_d = self.sigma_m_rd['y']
        sigma_mz_d = self.sigma_m_rd['z']
        f_m_d = self.f_type_rd
        K_h_y = self.K_h['y']
        K_h_z = self.K_h['z']
        K_m = self.K_m
        K_crit_y = self.K_crit[1]["y"]
        K_crit_z = self.K_crit[1]["z"]

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def base():
            taux_6_11 = sigma_my_d / (f_m_d * K_h_y) + K_m * sigma_mz_d / (f_m_d * K_h_z) # equ6.11
            taux_6_12 = K_m * sigma_my_d / (f_m_d * K_h_y) + sigma_mz_d / (f_m_d * K_h_z) # equ6.12
            taux_6_33y = sigma_my_d / (f_m_d * K_h_y * K_crit_y) # equ6.33
            taux_6_33z = sigma_mz_d / (f_m_d * K_h_z * K_crit_z) # equ6.33
            return taux_6_11, taux_6_12, taux_6_33y, taux_6_33z
        
        base_val = base()
        latex = base_val[0]
        self.taux_m_rd['equ6.11'] = base_val[1][0]
        self.taux_m_rd['equ6.12'] = base_val[1][1]
        self.taux_m_rd['equ6.33y'] = base_val[1][2]
        self.taux_m_rd['equ6.33z'] = base_val[1][3]
        

        if compression and isinstance(compression, Compression):
            sigma_c_0_d = compression.sigma_c_0_rd
            f_c_0_d = compression.f_type_rd
            K_c_y = compression.kc_Axe[1]['y']
            K_c_z = compression.kc_Axe[1]['z']
            taux_6_2 = compression.taux_c_0_rd['equ6.2']
            # taux_6_23 = compression.taux_c_0_rd['equ6.23']
            # taux_6_24 = compression.taux_c_0_rd['equ6.24']
            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def comp(taux_6_11, taux_6_12, taux_6_2, taux_6_33y, taux_6_33z):
                taux_6_19 = taux_6_2**2 + taux_6_11 # equ6.19
                taux_6_20 = taux_6_2**2 + taux_6_12 # equ6.20
                taux_6_23 = sigma_c_0_d / (f_c_0_d * K_c_y) # equ6.23
                taux_6_24 = sigma_c_0_d / (f_c_0_d * K_c_z) #equ6.24
                taux_6_35zyz = taux_6_33y** 2 + (sigma_mz_d / (f_m_d * K_h_z)) + taux_6_24 # equ6.35
                taux_6_35yzz = taux_6_33y  + (sigma_mz_d / (f_m_d * K_h_z)) ** 2 + taux_6_24 # equ6.35 interprétation
                taux_6_35yzy = taux_6_33z** 2 + (sigma_my_d / (f_m_d * K_h_y)) + taux_6_23 # equ6.35
                taux_6_35zyy = taux_6_33z + (sigma_my_d / (f_m_d * K_h_y)) ** 2 + taux_6_23 # equ6.35 interprétation
                return taux_6_19, taux_6_20, taux_6_35zyz, taux_6_35yzz, taux_6_35yzy, taux_6_35zyy
            
            compression_val = comp(self.taux_m_rd['equ6.11'], self.taux_m_rd['equ6.12'], taux_6_2, self.taux_m_rd['equ6.33y'], self.taux_m_rd['equ6.33z'])
            latex = latex + compression_val[0]
            self.taux_m_rd['equ6.19'] = compression_val[1][0]
            self.taux_m_rd['equ6.20'] = compression_val[1][1]
            self.taux_m_rd['equ6.35zyz'] = compression_val[1][2] # 1er item axe de flexion pas au carré, 2eme item axe de flexion au carré, 3eme axe de compression
            self.taux_m_rd['equ6.35yzz'] = compression_val[1][3]
            self.taux_m_rd['equ6.35yzy'] = compression_val[1][4]
            self.taux_m_rd['equ6.35zyy'] = compression_val[1][5]

        if traction and isinstance(traction, Traction):
            taux_6_1 = traction.taux_t_0_rd['equ6.1']
            @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def tract(taux_6_11, taux_6_12):
                taux_6_17 = taux_6_11 + taux_6_1 # equ6.17
                taux_6_18 = taux_6_12 + taux_6_1 # equ6.18
                return taux_6_17, taux_6_18
            
            traction_val = tract(self.taux_m_rd['equ6.11'], self.taux_m_rd['equ6.12'])
            latex = latex + traction_val[0]
            self.taux_m_rd['equ6.11'] = traction_val[1][0]
            self.taux_m_rd['equ6.12'] = traction_val[1][1]

        max_taux = max([v for v in self.taux_m_rd.values()])
        synthese = [
            ["Flexion bois", max_taux, None],
        ]
        self._add_synthese_taux_travail(synthese)
        return (latex, self.taux_m_rd)

