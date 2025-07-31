#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os
from PySide6.QtWidgets import QFileDialog
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import pandas as pd
import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc
from ourocode.eurocode.A0_Projet import Batiment

class Sismique(Batiment):
    ACTION = (
        "Permanente G",
        "Exploitation Q",
        "Neige normale Sn"
    )
    OCCUPATION = {"Étages à occupations corrélées": 0.8, "Étages à occupations indépendantes": 0.5, "Toiture": 1, "Autres": 1}
    CAT_IMPORTANCE = tuple(Batiment._data_from_csv(Batiment, os.path.join("sismique", "categorie_importance.csv")).index)
    CAT_IMPORTANCE_NS = tuple(Batiment._data_from_csv(Batiment, os.path.join("sismique", "categorie_importance_ns.csv")).index)
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
    CLASSE_DUCTILITE_NS = {
        "Garde-corps ou ornements": {"q":1, "Classe de ductilité": "DCL"},
        "Signalisations et panneaux d'affichage": {"q":1, "Classe de ductilité": "DCL"},
        "Cheminées, mâts et réservoirs sur poteaux consoles non cvt sur plus de la moitié de leur h totale": {"q":1, "Classe de ductilité": "DCL"},
        "Cheminées, mâts et réservoirs sur poteaux consoles non cvt sur moins de la moitié de leur h totale": {"q":2, "Classe de ductilité": "DCM"},
        "Murs de façade et intermédiaires": {"q":2, "Classe de ductilité": "DCM"},
        "Cloisons et façades": {"q":2, "Classe de ductilité": "DCM"},
        "Éléments de fixations des meubles lourds et des biblio. supportés par les planchers": {"q":2, "Classe de ductilité": "DCM"},
        "Éléments de fixations des faux-plafonds et autres dispositifs légers de fixation": {"q":2, "Classe de ductilité": "DCM"},
    }
    TYPE_DOMMAGES = ("Fragiles", "Ductiles", "Désolidarisé/Pas de risque")
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
        self.gravity_loads = {}
        self._loads = {"default": {"Zi": 0*si.m, "load": 0*si.kN}}
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
    
    def add_gravity_load(self,name: str, load: float, surface: si.m**2, etage: str=Batiment.ETAGE, z_i: float=0, action: str=ACTION, categorie_Q: str=Batiment.CAT_TYPE, occupations: str=OCCUPATION, comment: str=""):
        """Ajoute une charge gravitaire au bâtiment, cela permet de considérer la masse par niveau sur le bâtiment.
        Attention ne pas oublier les charges G de mur, de menuiserie, d'élément technique et autre.

        Args:
            name (str): nom de la charge.
            load (float): charge gravitaire en kN/m².
            surface (float): surface d'application de la charge.
            etage (str): étage auquel est appliquée la charge.
            z_i (float): est la hauteur de l'étage i en mètres depuis les fondations ou le sommet d'un soubassement rigide.
            action (str): type d'action de la charge.
            categorie_Q (str): catégorie d'exploitation (valable uniquement pour les charges d'exploitation Q, mettre "Aucune" pour les autres charges).
            occupations (str): type d'occupation du bâtiment (valable uniquement pour les charges d'exploitation Q, mettre "Autres" pour les autres charges).
            comment (str): commentaire sur la charge.
        """
        def coef_psy(cat=None):
            """Retourne les caractéristiques psy sous forme de dictionnaire"""
            dict_psy = {"Vent": {}, "Température": {}}
            if cat:
                dict_psy[cat] = {}
            if self.alt.value > 1000:
                dict_psy["Neige > 1000m"] = {}
            else:
                dict_psy["Neige <= 1000m"] = {}

            data_csv_psy = self._data_from_csv("coeff_psy.csv")
            psy_columns = data_csv_psy.columns.to_list()

            for psy_i in psy_columns:
                for key in dict_psy.keys():
                    dict_psy[key][psy_i] = data_csv_psy.loc[key].loc[psy_i]
            return dict_psy

        def _key_action_psy(action_variable):
            if action_variable == "Exploitation Q":
                index = categorie_Q
            elif action_variable == "Neige normale Sn":
                if self.alt.value > 1000:
                    index = "Neige > 1000m"
                else:
                    index = "Neige <= 1000m"
            return index
            
        load = abs(load) * si.kN / si.m**2
        surface = surface * si.m**2
        value = {
            "Zi": z_i * si.m,
            "Charge gravitaire": load,
            "Surface": surface,
            "Action": action,
            "Catégorie": categorie_Q,
            "Occupations": occupations,
            "Commentaire": comment
        }
        coeff_occupation = 1
        psy2 = 1
        if self._loads.get("default"):
            self._loads.pop("default")
            
        if action != "Permanente G":
            if action == "Exploitation Q":
                psy2 = coef_psy(categorie_Q)[_key_action_psy(action)]["psy2"]
                coeff_occupation = self.OCCUPATION[occupations]
            else:
                psy2 = coef_psy()[_key_action_psy(action)]["psy2"]
        load_f = load * surface * coeff_occupation * psy2
        if self._loads.get(etage):
            self._loads[etage]["load"] = self._loads[etage]["load"] + load_f
        else:
            self._loads[etage] = {"load": load_f, "Zi": z_i * si.m}
        if self.gravity_loads.get(etage):
            self.gravity_loads[etage][name] = value
        else:
            self.gravity_loads[etage] = {name: value}
        return value

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
    
    def show_spectre_elastique_calcul(
        self, 
        direction: str=("x", "y"), 
        screenshot: bool = ("False", "True"),
        filepath: str=None
        ):
        """
        Affiche le spectre de calcul pour l'analyse élastique.

        Args:
            direction (str): direction du spectre ("x" ou "y")
            screenshot (bool): si True, enregistre le graphique
            filepath (str): chemin d'enregistrement du graphique, si ce dernier est vide, 
                alors une boite de dialogue s'ouvre pour choisir le chemin.
        """
        q = self.coeff_comportement[direction]["q"]
        array = np.array([])
        a_g = self.a_g[1]
        S = float(self.type_spectre_table["S"])
        TB = float(self.type_spectre_table["TB"])
        TC = float(self.type_spectre_table["TC"])
        TD = float(self.type_spectre_table["TD"])
        beta = 0.2
        for T1 in np.arange(0, 4, 0.01):
            if 0 <= T1 <= TB:
                S_d_t1 = a_g * S * (2/3 + T1/TB * (2.5/q - 2/3))
            elif TB <= T1 <= TC:
                S_d_t1 = a_g * S * 2.5/q
            elif TC <= T1 <= TD:
                res1 = a_g * S * 2.5/q * (TC/T1)
                res2 = a_g * beta
                S_d_t1 = max(res1, res2)
            elif TD <= T1:
                res1 = a_g * S * 2.5/q * ((TC*TD)/T1**2)
                res2 = a_g * beta
                S_d_t1 = max(res1, res2)
            array = np.append(array, S_d_t1.value)
        color = "blue"
        if direction == "y":
            color = "red"
        plt.figure(figsize=(10, 5))
        plt.plot(np.arange(0, 4, 0.01), array, color=color)
        plt.title(f"Spectre élastique de calcul / classe de sol {self.classe_sol} / direction {direction} / q={q}")
        plt.xlabel("Période T (s)")
        plt.ylabel("Accélération Sd,T1 (m/s^2)")
        plt.fill_between(np.arange(0, 4, 0.01), array, color=color, alpha=0.2)
        plt.grid()
        if screenshot:
            if not filepath:
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
            plt.savefig(filepath)
            return filepath
        else:
            plt.show()

    @property
    def T1(self):
        """
        Retourne les periodes de calcul selon EN 1998-1 §4.3.3.2.2
        """
        K_b_x = self.Kbx.value
        K_b_y = self.Kby.value
        m_total = sum(load["load"].value/10 for load in self._loads.values())
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            T_1_x = 2 * sqrt(m_total/K_b_x) # en s
            T_1_y = 2 * sqrt(m_total/K_b_y) # en s
            return {"x": T_1_x, "y": T_1_y}
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
        m_total = sum(load["load"].value/9.81 for load in self._loads.values()) * si.kg
        Sd_t = self.Sd_t[1]
        S_d_T1_x = Sd_t["x"]
        S_d_T1_y = Sd_t["y"]
        T1 = self.T1[1]
        TC = float(self.type_spectre_table["TC"])
        lamb_x = 0.85 if T1["x"] <= 2 * TC and len(self._loads.keys()) > 2 else 1
        lamb_y = 0.85 if T1["y"] <= 2 * TC and len(self._loads.keys()) > 2 else 1
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_b_x = S_d_T1_x * m_total * lamb_x # form4.5
            F_b_y = S_d_T1_y * m_total * lamb_y # form4.5
            return {"x": F_b_x, "y": F_b_y}
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
            delta = 1 + 1.2 * (x/Le) # form4.12
            return delta
        return val()

    def Fi(self, etage: str=Batiment.ETAGE):
        """
        Retourne l'effort horizontal équivalent à l'étage i selon EN 1998-1 §4.3.3.2.3.
        Attention cette formule ne fonctionne que si les déplacements horizontaux croissent linéairement suivant la hauteur.
        Les efforts n'intègrent pas le coefficient de torsion accidentelle !
        
        Args:
            etage (str): est le nom de l'étage considéré.
        """
        sum_zj_mj = 0
        for level, load in self._loads.items():
            z_i = load["Zi"]
            sum_zj_mj = sum_zj_mj + z_i * load["load"].value/9.81*si.kg
        m_i = self._loads.get(etage)["load"].value/9.81*si.kg
        z_i = self._loads.get(etage)["Zi"]
        F_b_x = self.Fb[1]["x"]
        F_b_y = self.Fb[1]["y"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            niveau = etage
            F_i_x = F_b_x * (z_i*m_i)/sum_zj_mj  # form4.13
            F_i_y = F_b_y * (z_i*m_i)/sum_zj_mj # form4.13
            return {"x": F_i_x, "y": F_i_y}
        return val()
    
    @property
    def Fi_table(self):
        """
        Retourne les efforts horizontaux équivalents à chaque étage selon EN 1998-1 §4.3.3.2.3.
        Attention cette méthodes ne fonctionne que si les déplacements horizontaux croissent linéairement suivant la hauteur.
        Les efforts n'intègrent pas le coefficient de torsion accidentelle !
        """
        dict_Fi = pd.DataFrame(columns=["Fi,x", "Fi,y"])
        F_b_x = self.Fb[1]["x"]
        F_b_y = self.Fb[1]["y"]
        sum_zj_mj = 0
        for level, load in self._loads.items():
            z_j = load["Zi"]
            sum_zj_mj = sum_zj_mj + z_j * load["load"].value/9.81*si.kg
        for level, load in self._loads.items():
            if level != "default":
                m_i = load["load"].value/9.81*si.kg
                z_i = load["Zi"]
                F_i_x = F_b_x * (z_i*m_i)/sum_zj_mj  # form4.13
                F_i_y = F_b_y * (z_i*m_i)/sum_zj_mj # form4.13
                dict_Fi.loc[level] = [F_i_x, F_i_y]
        return dict_Fi

    def ds(self, de: float, direction: str=("x", "y")):
        """
        Retourne le déplacement de calcul dû à l'action sismique de calcul avec prise en compte du coefficient de comportement
        conformément à EN 1998-1 §4.3.4.

        Args:
            de (float): déplacement du même point déterminé par une analyse linéaire basée sur le spectre de réponse de calcul en mm, 
                conformément à EN 1998-1 §3.2.2.5
            direction (str): direction du déplacement selon l'axe x ou y du bâtiment.
        """
        d_e = de * si.mm
        q_d = self.coeff_comportement[direction]["q"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            d_s = d_e * q_d # form4.23
            return d_s
        return val()

    def coeff_second_ordre(self, dr: float, V_tot: si.kN, etage: str=Batiment.ETAGE):
        """
        Retourne le coefficient de second ordre selon EN 1998-1 §4.3.5.2.2
        
        Args:
            dr (float): déplacement relatif de calcul entre étages en mm, pris comme la diférence de 
                déplacement latéral entre le bas et le haut du niveau considéré. Calculé conformément à EN 1998-1 §4.3.4
            V_tot (float): effort tranchant sismique total au niveau de l'étage considéré en kN.
                Ne pas oublier de multiplier par le coefficient de torsion accidentelle !
            etage (str): nom de l'étage considéré.
        """
        P_tot = 0
        V_tot = V_tot * si.kN
        d_r = dr * si.mm
        z_i = self._loads.get(etage)["Zi"]
        h_lvl = self.h_bat
        for level, load in self._loads.items():
            if level != "default":
                z_j = load["Zi"]
                if z_i <= z_j:
                    P_tot = P_tot + load["load"]
                if h_lvl < abs(z_j - z_i) and level != etage:
                    h_lvl = abs(z_j - z_i)
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            theta_r = (P_tot * d_r) / (V_tot * h_lvl) # form4.28
            if theta_r <= 0.1:
                coeff_P_delta = 1
            elif 0.1 < theta_r <= 0.2:
                coeff_P_delta = 1 /( 1 - theta_r)
            return coeff_P_delta
        return val()
    
    def taux_limitations_dommages(self, dr: float, etage: str=Batiment.ETAGE, type_dommages: str=TYPE_DOMMAGES):
        """
        Retourne le taux de limitation des dommages selon EN 1998-1 §4.4.3.2.
        
        Args:
            dr (float): déplacement relatif de calcul entre étages en mm, pris comme la diférence de 
                déplacement latéral entre le bas et le haut du niveau considéré. Calculé conformément à EN 1998-1 §4.3.4
            etage (str): nom de l'étage considéré.
            type_dommages (str): type de dommage considéré:
                Fragiles: pour les bâtiments ayant des éléments non structuraux composés de matériaux fragiles fixés à la structure
                Ductiles: pour les bâtiments ayant des éléments non structuraux ductiles
                Désolidarisé/Pas de risque: pour les bâtiments ayant des éléments non structuraux fixés de manière à ne pas interférer 
                    avec les déformations de la structure ou n'ayant pas d'éléments non structuraux
        """
        nu = 0.4
        d_r = dr * si.mm
        h_lvl = self.h_bat
        z_i = self._loads.get(etage)["Zi"]
        for level, load in self._loads.items():
            if level != "default":
                z_j = load["Zi"]
                if h_lvl < abs(z_j - z_i) and level != etage:
                    h_lvl = abs(z_j - z_i)
        if type_dommages == "Fragiles":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_4_31 = (d_r * nu) / (h_lvl / 200) # equ4.31
                return taux_4_31
        elif type_dommages == "Ductiles":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_4_32 = (d_r * nu) / (h_lvl / 133) # equ4.32
                return taux_4_32
        elif type_dommages == "Désolidarisé/Pas de risque":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_4_33 = (d_r * nu) / (h_lvl / 100) # equ4.33
                return taux_4_33
        return val()

    def Fa(self, ma: float, Ta_x: float, Ta_y: float, z: float, type_element_ns: str=CLASSE_DUCTILITE_NS, cat_importance_ns: str=CAT_IMPORTANCE_NS):
        """
        Retourne l'effort sismique horizontal Fa à appliquer au centre de gravité des éléments non structuraux selon EN 1998-1 §4.3.5.2.

        Args:
            ma (float): masse de l'élément en kg.
            Ta_(x/y) (float): période fondamentale de vibration de l'élément non structural en seconde suivant la direction x ou y.
            z (float): hauteur de l'élément non structural au-dessus du niveau d'application de l'action sismique en m.
        """
        z = z * si.m
        h_bat = self.h_bat
        T_a_x = Ta_x * si.s
        T_a_y = Ta_y * si.s
        m_a = ma * si.kg
        a_g = self.a_g[1]
        g = 9.81 * si.m/si.s**2
        S = float(self.type_spectre_table["S"]) * si.m*si.s**-2
        T_1_x = self.T1[1]["x"] * si.s
        T_1_y = self.T1[1]["y"] * si.s

        file = os.path.join("sismique", "categorie_importance_ns.csv")
        gamma_a = self._data_from_csv(file).loc[cat_importance_ns]["gamma_a"]
        q_a = self.CLASSE_DUCTILITE_NS[type_element_ns]["q"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            alpha = a_g / g
            S_a_x = alpha * S *(3*(1 + z/h_bat)/(1 + (1-T_a_x/T_1_x)**2)-0.5)
            S_a_x = max(S_a_x, alpha * S)
            F_a_x = (S_a_x * m_a * gamma_a) / q_a # form4.24 dans le sens x
            S_a_y = alpha * S *(3*(1 + z/h_bat)/(1 + (1-T_a_y/T_1_y)**2)-0.5)
            S_a_y = max(S_a_y, alpha * S)
            F_a_y = (S_a_y * m_a * gamma_a) / q_a # form4.24 dans le sens y
            return {"x": F_a_x, "y": F_a_y}
        return val()

    def F_sismique_final_capacite(
        self,
        etage: str=Batiment.ETAGE, 
        gamma_d_x: str=("Rupture fragile", "Rupture ductile"), 
        gamma_d_y: str=("Rupture fragile", "Rupture ductile"), 
        Omega_x: float=1, 
        Omega_y: float=1, 
        eta_torsion_x: float=1, 
        eta_torsion_y: float=1, 
        P_delta_x: float=1, 
        P_delta_y: float=1
        ):
        """
        Retourne l'effort sismique final pour un dimensionnement bois des éléments en capacités.

        Args:
            etage (str): nom de l'étage considéré.
            gamma_d (str): type de rupture considéré:
                Fragile: tels que l'effort tranchant dans les diaphragmes en béton, cisaillant dans les embrèvements, 
                    assemblages collés, connecteurs par plaques embouties, instabilité de flambement déversement.
                Ductile: tels que l'effort tranchant dans les diaphragmes en panneaux dérivés du bois cloué, 
                    assemblages par pointes ou tiges avec raideurs faibles.
            omega (float): facteur de sur résistance des éléments dissipatifs (Rd/Ed).
            eta_torsion (float): facteur de prise en compte de la torsion accidentelle.
            P_delta (float): facteur de prise en compte des effets du second ordre.
        """
        F_i_x = self.Fi_table.loc[etage]["Fi,x"]
        F_i_y = self.Fi_table.loc[etage]["Fi,y"]
        if gamma_d_x == "Rupture fragile":
            gamma_d_x = 1.3
        elif gamma_d_x == "Rupture ductile":
            gamma_d_x = 1.1
        if gamma_d_y == "Rupture fragile":
            gamma_d_y = 1.3
        elif gamma_d_y == "Rupture ductile":
            gamma_d_y = 1.1
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_s_fin_x = F_i_x * gamma_d_x * Omega_x * eta_torsion_x * P_delta_x
            F_s_fin_y = F_i_y * gamma_d_y * Omega_y * eta_torsion_y * P_delta_y
            return {"x": F_s_fin_x, "y": F_s_fin_y}
        return val()

    def F_sismique_final_dissipatif(
        self,
        etage: str=Batiment.ETAGE, 
        eta_torsion_x: float=1, 
        eta_torsion_y: float=1, 
        P_delta_x: float=1, 
        P_delta_y: float=1
        ):
        """
        Retourne l'effort sismique final pour un dimensionnement bois des éléments dissipatifs.

        Args:
            etage (str): nom de l'étage considéré.
            eta_torsion (float): facteur de prise en compte de la torsion accidentelle.
            P_delta (float): facteur de prise en compte des effets du second ordre.
        """
        F_i_x = self.Fi_table.loc[etage]["Fi,x"]
        F_i_y = self.Fi_table.loc[etage]["Fi,y"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_s_fin_x = F_i_x * eta_torsion_x * P_delta_x
            F_s_fin_y = F_i_y * eta_torsion_y * P_delta_y
            return {"x": F_s_fin_x, "y": F_s_fin_y}
        return val()
        
