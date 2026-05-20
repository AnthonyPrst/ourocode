# coding in UTF-8
import os
import sys
from math import *
import pandas as pd

import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.ec3.assemblage.tige import Tige
from ourocode.eurocode.ec3.element_droit.plat import Plat
from ourocode.eurocode.ec5.element_droit.barre import Barre


class Platine_assise_compression_beton(Tige):
    GAMMA_C = {"Fondamentale": 1.5, "Accidentelle": 1.2}
    CLASSE_CONCRETE = list(Tige._data_from_csv(Tige, "caracteristique_meca_beton.csv").index)[2:]

    def __init__(self, gorge: si.mm, classe_beton: str=CLASSE_CONCRETE, *args, **kwargs):
        """Crée une platine d'assise comprimé avec support béton selon le §6.2.5 de l'EN1993-1-8.
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py
        Cette classe ne prend pas en compte pour le moment les dimensions connues de la fondation ce qui place l'utilisateur du côté de la sécurité.

        Attention:
        Si l’appui est constitué de mortier de scellement, il devra respecter les dispositions suivantes :  
            - De manière générale, l'épaisseur du scellement ne doit pas excéder 0,2 fois la plus petite largeur de la plaque 
            métallique d'assise, et sa résistance caractéristique ne doit pas être inférieure à 0,2 fois la résistance caractéristique 
            fck du béton de fondation. 
            - Dans les cas où l'épaisseur du scellement est supérieure à 50 mm, il convient que la résistance caractéristique du 
            scellement soit au moins égale à celle du béton de fondation.

        On détermine tout d'abord la la dimension effective additionnelle (fonction "c") depuis la face d'une semelle ou d'une ame métallique en mm 
        selon l'EN1993-1-8 §6.2.5(4).
        Puis on utilise la fonction "taux_compression" avec l'aire efficace de compression sur le béton Aef voir EN1993-1-8 §6.2.5(6) 
        et le guide de dimensionnement des assemblages acier édité par le CODIFAB.
        
        Args:
            gorge (int): gorge de la soudure qui relie le T du tronçon
            classe_beton (str): classe du béton
        """
        super().__init__(*args, **kwargs)
        self.gorge = gorge * si.mm
        self.classe_beton = classe_beton
        caract_meca = self._data_from_csv("caracteristique_meca_beton.csv").loc[self.classe_beton]
        self.fck = float(caract_meca.loc["fck"]) * si.MPa

    @property
    def f_jd(self):
        """Retourne la valeur de calcul de la résistance à la compression localisée f_jd en MPa selon EN 1993-1-8 §6.2.5.

        Formule : f_jd = β_j × (α_cc × f_ck / γ_c) × α_bf.
        Coefficients conservateurs par défaut : β_j = 1, α_cc = 1, α_bf = 1.

        Returns:
            tuple: (latex_string, f_jd) en MPa (avec unité si.MPa).
        """
        alpha_cc = 1
        alpha_bf = 1
        beta_i = 1
        f_ck = self.fck
        gamma_c = self.GAMMA_C["Fondamentale"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_jd = beta_i * (alpha_cc * f_ck / gamma_c) * alpha_bf
            return f_jd
        return val()

    def c(self):
        """Retourne la largeur effective complémentaire C depuis la face d'une semelle ou d'une âme selon EN 1993-1-8 §6.2.5(4).

        Formule : C = t × √(f_y / (3 × f_jd × γ_M0)).

        Returns:
            tuple: (latex_string, C) en mm (avec unité si.mm).
        """
        t = self.t
        f_y = self.fy
        f_jd = self.f_jd[1]
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            C = t * sqrt(f_y / (3 * f_jd * gamma_M0))
            return C
        return val()

    def taux_compression(self, N_c_Ed:si.kN, Aef: si.mm**2):
        """Retourne le taux de travail en compression de la plaque d'assise sur le béton selon EN 1993-1-8 §6.2.5(6).

        Formule : taux = N_c,Ed / (A_ef × f_jd).

        Args:
            N_c_Ed (float): Effort de compression de calcul sur la plaque d'assise en kN.
            Aef (float): Aire efficace de compression sur le béton en mm²
                (tenant compte de la largeur effective complémentaire C obtenue via la méthode c()).

        Returns:
            tuple: (latex_string, taux_compression) sans unité.
        """
        N_c_Ed = N_c_Ed * si.kN
        A_ef = Aef * si.mm**2
        f_jd = self.f_jd[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_compression = N_c_Ed / (A_ef * f_jd)
            return taux_compression
        result = val()
        self._add_synthese_taux_travail([["Compression sur plaque d'assise et support béton", result[1], None]])
        return result


class Platine_assise_compression_bois(Tige):
    CLASSE_WOOD = list(Tige._data_from_csv(Tige, "caracteristique_meca_bois.csv").index)[2:]

    def __init__(self, wood_beam: object, gorge: si.mm, alpha_wood: float=0, *args, **kwargs):
        """Crée une platine d'assise comprimé avec support bois selon le §6.2.5 de l'EN1993-1-8.
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py
        Pour ce placer en sécurité cette fonction considère un Kc90 = 1 pour le bois.

        On détermine tout d'abord la la dimension effective additionnelle (fonction "c") depuis la face d'une semelle ou d'une ame métallique en mm 
        selon l'EN1993-1-8 §6.2.5(4).
        Puis on utilise la fonction "taux_compression" avec l'aire efficace de compression sur le bois Aef voir EN1993-1-8 §6.2.5(6) 
        et le guide de dimensionnement des assemblages acier édité par le CODIFAB.
        
        Args:
            wood_beam (object): objet correspondant à la classe Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
                             ou bien objet Element ou dérivé de cet objet provenant du module EC3_Element_droit.py
            gorge (int): gorge de la soudure qui relie le T du tronçon
            alpha_wood (str): angle du fil du bois vis à vis de l'effort de compression sur la plaque d'assise.
        """
        super().__init__(*args, **kwargs)
        self.wood_beam = wood_beam
        self.gorge = gorge * si.mm
        self.alpha_wood = alpha_wood

    def _f_jd(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Valeur de calcul de la résistance à la compression localisée de la plaque d'assise sur le bois
        """
        f_c_0_d = self.wood_beam._f_type_d("fc0k", loadtype, typecombi)[1]
        f_c_90_d = self.wood_beam._f_type_d("fc90k", loadtype, typecombi)[1]
        K_c90 = 1 # On se place du coté de la sécurité
        alpha_wood = self.alpha_wood
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_jd =  (f_c_0_d / ((f_c_0_d /(K_c90 * f_c_90_d)) * sin(radians(alpha_wood))**2 + cos(radians(alpha_wood))**2))
            return f_jd
        result = val()
        self.f_jd = result
        return result

    def c(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Retourne la largeur effective complémentaire C pour une platine sur support bois selon EN 1993-1-8 §6.2.5(4).

        Formule : C = t × √(f_y / (3 × f_jd × γ_M0)) + a_ef
        où a_ef = 0.8 × (a × √2) est la largeur efficace de la soudure.

        Args:
            loadtype: Classe de durée de charge du bois (via Barre.LOAD_TIME). Defaults to Barre.LOAD_TIME.
            typecombi: Type de combinaison d'actions (via Barre.TYPE_ACTION). Defaults to Barre.TYPE_ACTION.

        Returns:
            tuple: (latex_string, C) en mm (avec unité si.mm).
        """
        t = self.t
        f_y = self.fy
        f_jd = self._f_jd(loadtype, typecombi)[1]
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        gorge = self.gorge
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            a_ef = 0.8 * (gorge * sqrt(2)) # largeur efficace d'une soudure
            C = t * sqrt(f_y / (3 * f_jd * gamma_M0)) + a_ef
            return C
        return val()

    def taux_compression(self, N_c_Ed:si.kN, Aef: si.mm**2):
        """Retourne le taux de travail en compression de la plaque d'assise sur le bois selon EN 1993-1-8 §6.2.5(6).

        Formule : taux = N_c,Ed / (A_ef × f_jd).
        Appeler d'abord self._f_jd() puis self.c() pour obtenir A_ef.

        Args:
            N_c_Ed (float): Effort de compression de calcul sur la plaque d'assise en kN.
            Aef (float): Aire efficace de compression sur le bois en mm²
                (tenant compte de la largeur effective complémentaire C obtenue via la méthode c()).

        Returns:
            tuple: (latex_string, taux_compression) sans unité.
        """
        N_c_Ed = N_c_Ed * si.kN
        A_ef = Aef * si.mm**2
        latex_fjd, f_jd = self.f_jd
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_compression = N_c_Ed / (A_ef * f_jd)
            return taux_compression
        result = val()
        self._add_synthese_taux_travail([["Compression sur plaque d'assise et support bois", result[1], None]])
        return result
        
    
class Platine_assise_traction(Tige):
    def __init__(self, gorge: si.mm, effet_levier: bool=("False", "True"), *args, **kwargs):
        """Crée une platine d'about fléchie en traction selon le §6.2.6.5 de l'EN1993-1-8 et du CNC2M §4.3.3
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py
        Attention cette classe ne permet le calcul que des platines avec 2 boulons par rangée de tronçon en T équvalent et non 4 comme décrit dans l'annexe A.

        Pour utiliser cette classe on doit d'abord ajouter les rangées isolées de chaque tronçon équivalent avec les fonctions: 
            - ajouter_rangee_simple_unique
            - ajouter_rangee_interieure
            - ajouter_rangee_exterieure_non_raidie
            - ajouter_rangee_centrale_non_raidie.

        Ensuite on ajoute les groupes de rangées pour la vérification des modes de rupture globale avec les fonctions:
            - ajouter_groupe_rangees_interieures_centrales
            - ajouter_groupe_rangees_centrales
            - ajouter_groupe_rangees_interieures

        Puis on vérifie les taux de travail avec la fonction:
            - taux_traction        
        Args:
            gorge (int): gorge de la soudure qui relie le T du tronçon
            effet_levier (bool, optional): définit si un effet levier est appliqué sur la platine. 
                Il est louable de considéré qu'une platine reposant sur un support bois n'a pas d'effet levier, 
                en effet la raideur en compression perpendiculaire du bois est faible. Defaults to "False".
        """
        super().__init__(*args, **kwargs)
        self.gorge = gorge * si.mm
        self.effet_levier = effet_levier
        self._bolts_rows = []
        
    def _m(self, l:si.mm):
        """Dimension m selon la figure 6.2 de l'EN1993-1-8 §6.2.4.1.
        Args:
            l (int): longueur en mm entre la face du T/raidisseur et l'axe de la tige d'assemblage

        Returns:
            int: dimension m en mm
        """
        m = l - 0.8 * self.gorge * sqrt(2)
        return m


    def show_l_eff(self):
        """Affiche le schéma des rangées issue du CNC2M §4.3.3 (4) tab15 et 16 pour définir les longeurs efficaces"""
        file1 = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "CNC2M_4.3.3_tab15.png"
        )
        file2 = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "CNC2M_4.3.3_tab16.png"
        )
        image1 = Image.open(file1)
        image2 = Image.open(file2)
        image1.show()
        image2.show()

    def show_fig6_11(self):
        """Affiche la figure 6.11 de l'EN1993-1-8 pour définir le facteur alpha"""
        file = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "EN1993-1-8_fig6_11.png"
        )
        image = Image.open(file)
        image.show()

    def ajouter_rangee_simple_unique(self, ml:si.mm, e:si.mm, bp:si.mm):
        """Ajoute une rangée isolée simple, unique (pas plus de rangée) et non raidie sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 1 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            ml (si.mm): distance en mm entre la face du T/raidisseur et l'axe de la tige d'assemblage
            e (si.mm): distance de bord en mm
            bp (si.mm): longueur de la platine en mm
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        bp = bp * si.mm

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2 * pi * m
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(4*m + 1.25*e, bp)
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Simple non raidie", "m": m, "e_min": e, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_rangee_interieure(self, ml:si.mm, m2l:si.mm, e:si.mm, alpha: float):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 2 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
            m2l (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
            e (si.mm): distance de bord en mm
            alpha (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda 1 et 2 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        m2l = m2l * si.mm
        m2 = self._m(m2l)
        e = e * si.mm
        lamb1 = m / (m + e)
        lamb2 = m2 / (m + e)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2 * pi * m
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(4*m + 1.25*e, alpha*m)
            return l_eff_nc
        l_eff_nc = leff_nc()

        row = {"Nom": name, "type": "Interieur", "m": m, "e_min": e, "lamb1": lamb1, "lamb2": lamb2, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row
    
    def ajouter_rangee_exterieure_non_raidie(self, w: si.mm, mxl:si.mm, e:si.mm, ex:si.mm, bp:si.mm):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 3 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            w (si.mm): distance d'entraxe des tiges en mm
            mxl (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
            e (si.mm): distance de bord parallèle à la semelle en mm
            ex (si.mm): distance de bord peprpendiculaire à la semelle en mm
            bp (si.mm): longueur de la platine en mm
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        mxl = mxl * si.mm 
        mx = self._m(mxl)
        e = e * si.mm
        ex = ex * si.mm

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = min(2 * pi * mx, pi * mx + w, pi * mx + 2*e)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(4*mx + 1.25*e, 2*mx + 0.625*ex + w/2, 2*mx + 0.625*ex + e, bp/2)
            return l_eff_nc
        l_eff_nc = leff_nc()

        row = {"Nom": name, "type": "Extérieur non raidie", "m": mx, "e_min": min(e, ex), "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_rangee_centrale_non_raidie(self, ml:si.mm, e:si.mm):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 4 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
            e (si.mm): distance de bord en mm
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2 * pi * m
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = 4*m + 1.25*e
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Centrale non raidie", "m": m, "e_min": e, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_rangee_exterieure_raidie(self, mxl:si.mm, ex:si.mm, m2xl:si.mm, e2x:si.mm, alpha1: float, alpha2: float):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 5 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            mxl (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
            ex (si.mm): distance de bord peprpendiculaire à la semelle en mm
            m2xl (si.mm): distance en mm entre la face du raidisseur et l'axe de la tige d'assemblage
            e2x (si.mm): distance de bord peprpendiculaire au raidisseur en mm
            alpha1 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda11 et lambda21 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
            alpha2 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda12 et lambda22 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        mxl = mxl * si.mm 
        mx = self._m(mxl)
        m2xl = m2xl * si.mm
        m2x = self._m(m2xl)
        ex = ex * si.mm
        e2x = e2x * si.mm

        lamb11 = mx / (mx + ex)
        lamb21 = m2x / (mx + ex)
        lamb12 = m2x / (m2x + e2x)
        lamb22 = mx / (m2x + e2x)

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp1 = min(2 * pi * mx, pi * mx + 2*e2x)
            l_eff_cp2 = min(2 * pi * m2x, pi * m2x + 2*ex)
            l_eff_cp = min(l_eff_cp1, l_eff_cp2)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc1 = min(alpha1 * mx, alpha1 * mx - 2*mx - 0.625*ex + e2x)
            l_eff_nc2 = min(alpha2 * m2x, alpha2 * m2x - 2*mx - 0.625*e2x + ex)
            l_eff_nc = min(l_eff_nc1, l_eff_nc2)
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Extérieur raidie", "m": min(mx,m2x), "e_min": min(ex, e2x), "lamb11": lamb11, "lamb21": lamb21, "lamb12": lamb12, "lamb22": lamb22, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    ########### MODE DE RUPTURE GLOBALE ###################

    def ajouter_groupe_rangees_interieures_centrales(self, ml:si.mm, e:si.mm, m2l:si.mm, sum_p:si.mm, alpha: float):
        """Ajoute une rangée isolée intérieur sur la platine 
            pour le calcul du l,eff des modes circulaires et non circulaires en traction.
            Schéma 1 dans le tableau 16 du CNC2M (voir la fonction show_l_eff).

            Args:
                ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
                e (si.mm): distance de bord perpendiculaire à l'ame en mm
                m2l (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
                sum_p (si.mm): somme des entraxes p entres les tiges d'assemblages intérieur/central
                alpha (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                    Pour déterminer alpha on peut récupérer les lambda 1 et 2 en exécutant cette fonction avec un lambda aléatoire, puis
                    en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon groupé " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        m2l = m2l * si.mm
        m2 = self._m(m2l)
        sum_p = sum_p * si.mm

        lamb1 = m / (m + e)
        lamb2 = m2 / (m + e)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2*(pi * m + sum_p)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(alpha*m + sum_p, 4*m + 1.25*e + sum_p)
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Intérieure et centrale(s) en groupe", "m": m, "e_min": e, "lamb1": lamb1, "lamb2": lamb2, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_groupe_rangees_centrales(self, ml:si.mm, e:si.mm, sum_p:si.mm):
        """Ajoute une rangée isolée intérieur sur la platine 
            pour le calcul du l,eff des modes circulaires et non circulaires en traction.
            Schéma 2 dans le tableau 16 du CNC2M (voir la fonction show_l_eff).

            Args:
                ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
                e (si.mm): distance de bord perpendiculaire à l'ame en mm
                sum_p (si.mm): somme des entraxes p entres les tiges d'assemblages intérieur/central
        """
        name = "Tronçon groupé " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        sum_p = sum_p * si.mm

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2*(pi * m + sum_p)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = 4*m + 1.25*e + sum_p
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Centrales en groupe", "m": m, "e_min": e, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_groupe_rangees_interieures(self, ml:si.mm, e:si.mm, m21l:si.mm, m22l:si.mm, sum_p:si.mm, alpha1: float, alpha2: float):
        """Ajoute une rangée isolée intérieur sur la platine pour le calcul du l,eff 
        des modes circulaires et non circulaires en traction.
        Schéma 3 dans le tableau 16 du CNC2M (voir la fonction show_l_eff).

        Attention: Ce cas ne se rencontre que pour les assemblages de pièces tendues ou pour des encastrements avec raidisseurs entre rangées

        Args:
            ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
            e (si.mm): distance de bord perpendiculaire à l'ame en mm
            m21l (si.mm): distance en mm entre la face de la semelle haute et l'axe de la tige d'assemblage
            m22l (si.mm): distance en mm entre la face de la semelle basse et l'axe de la tige d'assemblage
            sum_p (si.mm): somme des entraxes p entres les tiges d'assemblages intérieur/central
            alpha1 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda11 et lambda21 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
            alpha2 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda12 et lambda22 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon groupé " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        m21l = m21l * si.mm
        m21 = self._m(m21l)
        m22l = m22l * si.mm
        m22 = self._m(m22l)
        sum_p = sum_p * si.mm

        lamb11 = m / (m + e)
        lamb21 = m21 / (m + e)
        lamb12 = m / (m + e)
        lamb22 = m22 / (m22 + e)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2*(pi * m + sum_p)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = (alpha1 + alpha2)*m -4*m + 1.25*e + sum_p
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Intérieures en groupe", "m": m, "e_min": e, "lamb11": lamb11, "lamb21": lamb21, "lamb12": lamb12, "lamb22": lamb22, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def taux_traction(self, N_T_Ed:si.kN):
        """Retourne le taux de travail global en traction de la platine d'about selon EN 1993-1-8 §6.2.6.5 et CNC2M §4.3.3.

        Parcourt toutes les rangées/groupes de tronçons ajoutés préalablement pour calculer F_T,Rd
        selon les trois modes de rupture (modes 1, 2, 3) et le résultat minimisant.
        La résistance globale est min(somme des F_T,Rd isolés, F_T,Rd groupé).

        Args:
            N_T_Ed (float): Effort de traction de calcul sur la platine en kN.

        Returns:
            tuple: (latex_string, taux_traction) sans unité.
        """
        N_T_Ed = N_T_Ed * si.kN
        F_T_Rd_isole = 0 * si.kN
        F_T_Rd_groupe = 0 * si.kN
        n_bl = 2
        t_f = self.t
        f_y = self.fy
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        latex_FtRd, F_t_Rd = self.FtRd
        for row in self._bolts_rows:
            name = row["Nom"]
            m = row["m"]
            e_min = row["e_min"]
            l_eff_cp = row["l_eff_cp"][1]
            l_eff_nc = row["l_eff_nc"][1]

            if self.effet_levier:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    name
                    effet_levier = True
                    l_eff_1 = min(l_eff_cp, l_eff_nc)
                    M_pl1_Rd = (l_eff_1 * t_f**2 * f_y) / (4 * gamma_M0)
                    F_T1_Rd = (4 * M_pl1_Rd) / m
                    
                    l_eff_2 = l_eff_nc
                    n = min(e_min, 1.25*si.m)
                    M_pl2_Rd = (l_eff_2 * t_f**2 * f_y) / (4 * gamma_M0)
                    F_T2_Rd = (2 * M_pl2_Rd + n * n_bl * F_t_Rd) / (m + n)

                    F_T3_Rd = n_bl * F_t_Rd # n,bl est le nombre de boulon dans la rangée en traction

                    F_T_Rd = min(F_T1_Rd, F_T2_Rd, F_T3_Rd)
                    return F_T_Rd
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    name
                    effet_levier = False
                    l_eff_1 = min(l_eff_cp, l_eff_nc)
                    M_pl1_Rd = (l_eff_1 * t_f**2 * f_y) / (4 * gamma_M0)
                    F_T12_Rd = (2 * M_pl1_Rd) / m

                    F_T3_Rd = n_bl * F_t_Rd # n,bl est le nombre de boulon dans la rangée en traction

                    F_T_Rd = min(F_T12_Rd, F_T3_Rd)
                    return F_T_Rd
            row["F,T,Rd"] = val()
            if "groupé" in name:
                if F_T_Rd_groupe == 0*si.kN:
                    F_T_Rd_groupe = row["F,T,Rd"][1]
                else:
                    F_T_Rd_groupe = min(F_T_Rd_groupe, row["F,T,Rd"][1])
            else:
                F_T_Rd_isole = F_T_Rd_isole + row["F,T,Rd"][1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def taux():
            F_T_Rd = min(F_T_Rd_isole, F_T_Rd_groupe)
            taux_traction = N_T_Ed / F_T_Rd
            return taux_traction
        latex = "".join((row["F,T,Rd"][0] for row in self._bolts_rows))
        taux = taux()
        self._add_synthese_taux_travail([["Traction d'une plaque d'assise avec attache en T", taux[1], None]])
        return (latex + taux[0], taux[1])


if __name__ == "__main__":
    # Création d'un pied de poteau HE200B
    pied_poteau = Platine_about(
        h=200,     # hauteur du profilé en mm
        b=200,     # largeur des ailes en mm
        tw=9,      # épaisseur de l'âme en mm
        tf=15,     # épaisseur des ailes en mm
        r=18,      # rayon de raccordement âme-aile en mm
        type_profil="HE"
    )

    # Ajout de boulons M20
    positions_boulons = [
        (-80, 50),   # boulon en haut à gauche
        (80, 50),    # boulon en haut à droite
        (-80, -50),  # boulon en bas à gauche
        (80, -50)    # boulon en bas à droite
    ]
    pied_poteau.ajouter_boulons(positions_boulons, 20)

    # Affichage
    pied_poteau.afficher()
