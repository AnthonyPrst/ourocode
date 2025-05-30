#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os
from re import M
from _pytest.stash import T
import matplotlib.pyplot as plt
from math import sqrt
from PySide6.QtWidgets import QFileDialog
import numpy as np
import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc
from ourocode.eurocode.A0_Projet import Batiment

class Sismique(Batiment):
    CAT_IMPORTANCE = tuple(Batiment._data_from_csv(Batiment, os.path.join("sismique", "categorie_importance.csv")).index)
    CLASSE_SOL = tuple(Batiment._data_from_csv(Batiment, os.path.join("sismique", "classe_sol.csv")).index)
    TYPE_SPECTRE = ("Type 2", "Type 1")
    AGR = {"Zone 1": 0.4, "Zone 2": 0.7, "Zone 3": 1.1, "Zone 4": 1.6, "Zone 5": 3}
    CLASSE_DUCTILITE = {
        "Consoles, poutres, arcs avec 2 ou 3 assemblages brochés": {"q":1.5, "Classe de ductilité": "DCL"},
        "Panneaux de Murs avec diaphragmes collés, assemblés par clous et boulons": {"q":2, "Classe de ductilité": "DCM"},
        "Treillis assemblés par broches ou boulons": {"q":2, "Classe de ductilité": "DCM"},
        "Structures mixtes Ossature bois (contreventement) + remplissage non porteur": {"q":2, "Classe de ductilité": "DCM"},
        "Portiques hyperstatiques assemblés par broches ou boulons": {"q":2.5, "Classe de ductilité": "DCM"},
        "Panneaux de murs avec diaphragmes cloués, assemblés par clous et boulons": {"q":3, "Classe de ductilité": "DCH"},
    }

    def __init__(
        self, 
        Kbx: float,
        Kby: float,
        type_constructif_x: str=CLASSE_DUCTILITE,
        type_constructif_y: str=CLASSE_DUCTILITE,
        regulier_plan: bool=("True", "False"), 
        regulier_elevation: bool=("True", "False"), 
        cat_importance: str=CAT_IMPORTANCE, 
        classe_sol: str=CLASSE_SOL, 
        type_spectre: str=TYPE_SPECTRE,
        **kwargs
        ):
        """
        Créer une classe qui permet de calculer l'action sismique pour les bâtiments bois, selon la méthode des forces latérales.
        Attention, tout les bâtiments ne ce prêtent pas à ce type d'étude (voir EN 1998).
        Cette classe est hérité de la classe Batiment du module A0_Projet.py.

        Args:
            Kbx (float): Raideur du bâtiment selon la direction x en kN/m.
            Kby (float): Raideur du bâtiment selon la direction y en kN/m.
            type_constructif_x (str): Type de système constructif pour les éléments dissipatifs dans la longeur x du bâtiment.
            type_constructif_y (str): Type de système constructif pour les éléments dissipatifs dans la largeur y du bâtiment.
            regulier_plan (bool): Détermine si le bâtiment est régulier en plan.
            regulier_elevation (bool): Détermine si le bâtiment est régulier en élévation.
            cat_importance (str): Catégorie d'importance du bâtiment:
                - I : Bâtiments dans lesquels il n'y a aucune activité humaine nécessitant un séjour de longue durée.
                - II : 
                    - Habitations individuelles.
                    - Établissements recevant du public (ERP) de catégories 4 et 5.
                    - Habitations collectives de hauteur inférieure à 28 m.
                    - Bureaux ou établissements commerciaux non ERP, h ≤ 28 m, max. 300 pers. 
                    - Bâtiments industriels pouvant accueillir au plus 300 personnes.
                    - Parcs de stationnement ouverts au public.
                - III : 
                    - ERP de catégories 1, 2 et 3.
                    - Habitations collectives et bureaux, h > 28 m.
                    - Bâtiments pouvant accueillir plus de 300 personnes. 
                    - Établissements sanitaires et sociaux.
                    - Centres de production collective d'énergie.
                    - Établissements scolaires.
                - IV : 
                    - Bâtiments indispensables à la sécurité civile, la défense nationale et le maintien de l'ordre public.
                    - Bâtiments assurant le maintien des communications, la production et le stockage d'eau potable, la distribution publique de l'énergie.
                    - Bâtiments assurant le contrôle de la sécurité aérienne.
                    - Établissements de santé nécessaires à la gestion de crise.
                    - Centres météorologiques.
            classe_sol (str): Classe du sol selon EN 1998-1 §3.1.2.
            type_spectre (str): Type de spectre. Le spectre en France métropolitaine est de type 2 (magnitude < 6) pour le reste type 1 (magnitude >= 6).
        """
        super().__init__(**kwargs)
        self.cat_importance = cat_importance
        self.classe_sol = classe_sol
        self.type_spectre = type_spectre
        self.type_constructif_x = type_constructif_x
        self.type_constructif_y = type_constructif_y
        self.type_constructif = {"x": type_constructif_x, "y": type_constructif_y}
        self.coeff_comportement = {"x": self.CLASSE_DUCTILITE[type_constructif_x], "y": self.CLASSE_DUCTILITE[type_constructif_y]}
        self.Kbx = Kbx * si.kN / si.m
        self.Kby = Kby * si.kN / si.m
        self.K_b = {"Raideur du bâtiment dans la direction x": self.Kbx, "Raideur du bâtiment dans la direction y": self.Kby}
        self.regulier_plan = regulier_plan
        self.regulier_elevation = regulier_elevation
        self._has_to_be_analyzed()
        self._is_ok_for_lateral_force_method()

    def _has_to_be_analyzed(self):
        """Retourne True si le bâtiment doit être analysé selon la catégorie d'importance"""
        if self.region_sismique == "Zone 1":
            raise ValueError("Le bâtiment n'est pas à analyser au niveau sismique")
        elif self.cat_importance == "I":
            raise ValueError("Le bâtiment n'est pas à analyser au niveau sismique")
        elif self.cat_importance == "II" and self.region_sismique == "Zone 2":
            raise ValueError("Le bâtiment n'est pas à analyser au niveau sismique")
        else: 
            return True
    
    def _is_ok_for_lateral_force_method(self):
        T1 = self.T1[1]
        TC = float(self.type_spectre_table["TC"])
        if self.cat_importance == "IV":
            raise ValueError("Le bâtiment ne peut pas être analysé avec la méthode des forces latérales car le bâtiment est de catégorie IV")
        elif not self.regulier_elevation:
            raise ValueError("Le bâtiment ne peut pas être analysé avec la méthode des forces latérales car le bâtiment n'est pas régulier en élévation")
        elif T1["x"] > 4 * TC or T1["x"] > 2  or T1["y"] > 4 * TC or T1["y"] > 2:
            raise ValueError("Le bâtiment ne peut pas être analysé avec la méthode des forces latérales car le bâtiment est soumis à des modes de vibrations\
                de rang plus élevé que le mode fondamentale dans chaque direction principale (EN 1998-1 §4.3.3.2.1)")
        else:
            return True 

    @property
    def region_sismique(self):
        """Retourne la région sismique du bâtiment"""
        file = "carte_action_region.csv"
        df = self._data_from_csv(file, index_col=1)
        return df.loc[str(self.code_INSEE)]["Alea_sismique"]
    
    @property
    def cat_importance_table(self):
        """Retourne le dataframe de la catégorie d'importance choisi"""
        file = os.path.join("sismique", "categorie_importance.csv")
        data_cat_imp = self._data_from_csv(file)
        return data_cat_imp.loc[self.cat_importance]

    @property
    def classe_sol_table(self):
        """Retourne le dataframe de la classe de sol choisi"""
        file = os.path.join("sismique", "classe_sol.csv")
        data_classe_sol = self._data_from_csv(file)
        return data_classe_sol.loc[self.classe_sol]

    @property
    def type_spectre_table(self):
        """Retourne le dataframe du spectre choisi"""
        if self.type_spectre == "Type 2":
            file = os.path.join("sismique", "spectre_eleastique_h_type2.csv")
        else:
            file = os.path.join("sismique", "spectre_eleastique_h_type1.csv")
        return self._data_from_csv(file).loc[self.classe_sol]

    @property
    def type_constructif_table(self):
        """Retourne les classes de ductilité bois"""
        return self.CLASSE_DUCTILITE

    @property
    def a_gr(self):
        """Retourne l'accélération de base pour un sol de classe A"""
        return self.AGR[self.region_sismique] *si.m / si.s**2

    @property
    def a_g(self):
        """Retourne l'accélération de calcul pour un sol de classe A"""
        a_gr = self.a_gr
        gamma_1 = self.gamma_1
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            a_g = a_gr * gamma_1
            return a_g
        return val()

    @property
    def gamma_1(self):
        """Retourne le gamma_1 fonction de la catégorie d'importance"""
        return self.cat_importance_table["gamma_1"]
    

    def _spectre_elastique_calcul(self, T1: float, q: float, dir: str) -> tuple[str, dict[str, float]]:
        """
        Retourne le spectre elastique de calcul selon EN 1998-1 §3.2.2.5

        Args:
            T1 (float): période du bâtiment suivant la direction considéré x ou y.
            q (float): coefficient de comportement pour la direction considéré.
        """
        if dir == "x":
            dir = "x (sens de la longeur)"
        else:
            dir = "y (sens de la largeur)"
        beta = 0.2
        a_g = self.a_g[1]
        S = float(self.type_spectre_table["S"])
        TB = float(self.type_spectre_table["TB"])
        TC = float(self.type_spectre_table["TC"])
        TD = float(self.type_spectre_table["TD"])
        if 0 <= T1 <= TB:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                direction = dir
                S_d_t1 = a_g * S * (2/3 + T1/TB *(2.5/q - 2/3))
                return S_d_t1
            
        elif TB <= T1 <= TC:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                direction = dir
                S_d_t1 = a_g * S * 2.5/q
                return S_d_t1

        elif TC <= T1 <= TD:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                direction = dir
                res1 = a_g * S * 2.5/q * (TC/T1)
                res2 = a_g * beta
                S_d_t1 = max(res1, res2)
                return S_d_t1

        elif TD <= T1:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                direction = dir
                res1 = a_g * S * 2.5/q * ((TC*TD)/T1**2)
                res2 = a_g * beta
                S_d_t1 = max(res1, res2)
                return S_d_t1
        return val()

    @property
    def T1(self):
        """
        Retourne les periodes de calcul selon EN 1998-1 §3.2.2.5
        """
        K_b_x = self.Kbx.value
        K_b_y = self.Kby.value
        m_total = 100000
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            T_1x = 2 * sqrt(m_total/K_b_x) # en s
            T_1y = 2 * sqrt(m_total/K_b_y) # en s
            return {"x": T_1x, "y": T_1y}
        return val()

    @property
    def Sd_t(self):
        """
        Retourne le spectre elastique de calcul selon EN 1998-1 §3.2.2.5
        """
        Sd_t = {}
        latex = ""
        latex_periode, T1 = self.T1
        latex += latex_periode
        for dir, T in T1.items():
            latex_spectre, Sd_t[dir] = self._spectre_elastique_calcul(T, self.coeff_comportement[dir]["q"], dir)
            latex += latex_spectre
        return (latex, Sd_t)

    @property
    def Fb(self):
        """
        Retourne l'effort tranchant à la base de la structure selon EN 1998-1 §4.3.3.2.2
        """
        m = 100000
        Sd_t = self.Sd_t[1]
        S_d_T1_x = Sd_t["x"].value
        S_d_T1_y = Sd_t["y"].value
        T1 = self.T1[1]
        TC = float(self.type_spectre_table["TC"])
        lamb_x = 0.85 if T1["x"] <= 2 * TC else 1
        lamb_y = 0.85 if T1["y"] <= 2 * TC else 1
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_bx = S_d_T1_x * m * lamb_x * si.N
            F_by = S_d_T1_y * m * lamb_y * si.N
            return {"x": F_bx, "y": F_by}
        return val()
            
    def coeff_torsion_accidentelle(self, x: si.m, Le: si.m):
        """
        Retourne le coefficient de torsion accidentelle qui est à déterminer de cette manière si les raideurs latérales et de la masses sont symétriques.
        Selon EN 1998-1 §4.3.3.2.4
        
        Args:
            x (si.m): est la distance en plan de l'élément considéré au centre de masse du bâtiment en plan,
                mesurée perpendiculairement à la direction de l'action sismique considérée.
            Le (si.m): est la distance entre les deux éléments de contreventement extrêmes, mesurée perpendiculairement à la direction de
                l'action sismique considérée.
        """
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            delta = 1 + 0.6 * (x/Le)
            return delta
        return val()
            
            
 