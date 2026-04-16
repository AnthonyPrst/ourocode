#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import pandas as pd
import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core.renderer import handcalc
from ourocode.eurocode.core.batiment import Batiment

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
        """Initialise le calcul sismique d'un bâtiment bois selon la méthode des forces latérales (EN 1998-1 §4.3.3.2).

        Attention : tous les bâtiments ne se prêtent pas à cette méthode (voir EN 1998-1 §4.3.3.2.1).
        La classe vérifie automatiquement la nécessité d'une analyse sismique et la validité de la
        méthode des forces latérales lors de l'instanciation.

        Args:
            Kbx (float): Raideur latérale globale du bâtiment selon la direction x (sens de la longueur) en kN/m.
            Kby (float): Raideur latérale globale du bâtiment selon la direction y (sens de la largeur) en kN/m.
            type_constructif_x (str): Système constructif dissipatif dans la direction x.
                Doit être une clé de CLASSE_DUCTILITE. Détermine le coefficient de comportement q et la
                classe de ductilité (DCL, DCM, DCH) selon l'EN 1998.
            type_constructif_y (str): Système constructif dissipatif dans la direction y.
                Même valeurs que type_constructif_x.
            regulier_plan (bool): True si le bâtiment est régulier en plan selon EN 1998-1 §4.2.3.
            regulier_elevation (bool): True si le bâtiment est régulier en élévation selon EN 1998-1 §4.2.3.
            cat_importance (str): Catégorie d'importance du bâtiment selon EN 1998-1 §4.2.5 :
                - "I"  : Bâtiments sans activité humaine prolongée.
                - "II" : Habitations individuelles, ERP cat. 4/5, logements collectifs < 28 m,
                          bureaux/commerces non ERP ≤ 28 m (≤ 300 pers.), parcs de stationnement.
                - "III": ERP cat. 1-3, logements collectifs/bureaux > 28 m, bâtiments > 300 pers.,
                          établissements sanitaires, scolaires, centres énergétiques.
                - "IV" : Bâtiments indispensables à la sécurité civile, défense, communications,
                          eau potable, énergie, sécurité aérienne, santé de crise, météorologie.
            classe_sol (str): Classe de sol selon EN 1998-1 §3.1.2.
            type_spectre (str): Type de spectre de réponse.
                "Type 2" : France métropolitaine (magnitude < 6).
                "Type 1" : Zones à forte sismicité (magnitude ≥ 6).
            **kwargs: Arguments transmis à la classe parent Batiment.

        Raises:
            ValueError: Si le bâtiment n'est pas soumis à une analyse sismique (zone 1, cat. I,
                cat. II en zone 2) ou si la méthode des forces latérales n'est pas applicable
                (cat. IV, irrégularité en élévation, T1 > min(4TC, 2s)).
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
    
    def save_gravity_load_data(self, path: str=None) -> None:
        """Sauvegarde les données des charges gravitaires dans un fichier JSON.

        Args:
            path (str, optional): Chemin absolu du fichier de destination (.json).
                Si None, une boîte de dialogue s'ouvre pour sélectionner le fichier.
                Defaults to None.
        """
        super().save_data(self.gravity_loads, type_data="JSON", path=path)
    
    def load_gravity_load_data(self, path: str=None) -> dict:
        """Charge les données des charges gravitaires depuis un fichier JSON et reconstruit l'état interne.

        Relit le fichier produit par save_gravity_load_data et appelle _set_loads pour
        reconstituer le dictionnaire _loads utilisé dans les calculs sismiques.

        Args:
            path (str, optional): Chemin absolu du fichier source (.json).
                Si None, une boîte de dialogue s'ouvre pour sélectionner le fichier.
                Defaults to None.

        Returns:
            dict: Dictionnaire gravity_loads mis à jour, structuré par étage puis par nom de charge.
        """
        data = super().load_data(type_data="JSON", path=path)
        self.gravity_loads = data
        for etage, loads in data.items():
            for load_name, load in loads.items():
                self._set_loads(
                    load["Charge gravitaire"], 
                    load["Surface"], 
                    etage, 
                    load["Zi"], 
                    load["Action"], 
                    load["Catégorie"], 
                    load["Occupations"]
                    )
        return self.gravity_loads

    def _set_loads(self, load: float, surface: float, etage: str, z_i: float, action: str, categorie_Q: str, occupations: str):
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
                    if key != "Aucune":
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
            self._loads[etage] = {"load": load_f, "Zi": z_i}
        

    def add_gravity_load(
        self,
        name: str, 
        load: si.kN / si.m**2, 
        surface: si.m**2, 
        etage: str=Batiment.ETAGE, 
        z_i: float=0, 
        action: str=ACTION, 
        categorie_Q: str=Batiment.CAT_TYPE, 
        occupations: str=OCCUPATION, 
        comment: str=""
        ) -> dict:
        """Ajoute une charge gravitaire à un étage pour le calcul de la masse sismique.

        La charge est convertie en masse sismique selon la combinaison quasi-permanente
        (G + psi_2 × Q) conformément à EN 1998-1 §3.2.4. Ne pas oublier les charges permanentes
        de murs, menuiseries, éléments techniques et autres éléments non structuraux.

        Args:
            name (str): Identifiant unique de la charge dans l'étage considéré.
            load (float): Intensité de la charge gravitaire en kN/m².
            surface (float): Surface d'application de la charge en m².
            etage (str): Nom de l'étage auquel est appliquée la charge.
                Doit correspondre à une clé de ETAGE définie dans Batiment.
            z_i (float): Hauteur de l'étage i en mètres depuis les fondations
                ou le sommet d'un soubassement rigide. Defaults to 0.
            action (str): Type d'action de la charge selon ACTION.
                "Permanente G", "Exploitation Q" ou "Neige normale Sn".
            categorie_Q (str): Catégorie d'exploitation selon l'EC1-1-1.
                Utilisée uniquement pour les charges Q (psi_2). Defaults to "Aucune".
            occupations (str): Type d'occupation pour la réduction de la charge Q
                (coefficient phi selon EN 1998-1 §4.2.4).
                "Étages à occupations corrélées", "Étages à occupations indépendantes",
                "Toiture", ou "Autres". Defaults to "Autres".
            comment (str): Commentaire descriptif libre sur la charge. Defaults to "".

        Returns:
            dict: Dictionnaire décrivant la charge ajoutée avec les clés :
                "Zi", "Charge gravitaire", "Surface", "Action", "Catégorie",
                "Occupations", "Commentaire".
        """
            
        load = abs(load) * si.kN / si.m**2
        surface = surface * si.m**2
        z_i = z_i * si.m
        value = {
            "Zi": z_i,
            "Charge gravitaire": load,
            "Surface": surface,
            "Action": action,
            "Catégorie": categorie_Q,
            "Occupations": occupations,
            "Commentaire": comment
        }

        self._set_loads(
            load, 
            surface, 
            etage, 
            z_i, 
            action, 
            categorie_Q, 
            occupations
            )

        if self.gravity_loads.get(etage):
            self.gravity_loads[etage][name] = value
        else:
            self.gravity_loads[etage] = {name: value}
        return value

    @property
    def region_sismique(self) -> str:
        """Retourne la zone d'aléa sismique du bâtiment selon l'AN français (carte par code INSEE).

        Returns:
            str: Zone sismique parmi "Zone 1" à "Zone 5".
        """
        file = "carte_action_region.csv"
        df = self._data_from_csv(file, index_col=1)
        return df.loc[str(self.code_INSEE)]["Alea_sismique"]
    
    @property
    def cat_importance_table(self) -> 'pd.Series':
        """Retourne les données normatives de la catégorie d'importance sélectionnée.

        Returns:
            pd.Series: Ligne du fichier categorie_importance.csv correspondant à cat_importance,
                contenant notamment le coefficient d'importance gamma_I.
        """
        file = os.path.join("sismique", "categorie_importance.csv")
        data_cat_imp = self._data_from_csv(file)
        return data_cat_imp.loc[self.cat_importance]

    @property
    def classe_sol_table(self) -> 'pd.Series':
        """Retourne les paramètres de sol de la classe sélectionnée selon EN 1998-1 §3.1.2.

        Returns:
            pd.Series: Ligne du fichier classe_sol.csv pour la classe_sol choisie,
                contenant les paramètres S, TB, TC, TD du spectre de réponse.
        """
        file = os.path.join("sismique", "classe_sol.csv")
        data_classe_sol = self._data_from_csv(file)
        return data_classe_sol.loc[self.classe_sol]

    @property
    def type_spectre_table(self) -> 'pd.Series':
        """Retourne les paramètres du spectre élastique horizontal selon EN 1998-1 §3.2.2.

        Sélectionne le fichier de spectre type 1 ou type 2 selon type_spectre,
        puis retourne la ligne correspondant à la classe de sol.

        Returns:
            pd.Series: Paramètres S, TB, TC, TD pour le couple (type_spectre, classe_sol).
        """
        if self.type_spectre == "Type 2":
            file = os.path.join("sismique", "spectre_eleastique_h_type2.csv")
        else:
            file = os.path.join("sismique", "spectre_eleastique_h_type1.csv")
        return self._data_from_csv(file).loc[self.classe_sol]

    @property
    def type_constructif_table(self) -> dict:
        """Retourne le dictionnaire complet des classes de ductilité bois disponibles.

        Returns:
            dict: CLASSE_DUCTILITE : clés = description du système constructif,
                valeurs = dict {"q": coefficient de comportement, "Classe de ductilité": "DCL"/"DCM"/"DCH"}.
        """
        return self.CLASSE_DUCTILITE

    @property
    def a_gr(self) -> float:
        """Retourne l'accélération de référence a_gR pour un sol de classe A selon l'AN français.

        Returns:
            float: Accélération de référence en m/s² avec unité (si.m/si.s²),
                issue du dictionnaire AGR pour la zone sismique du projet.
        """
        return self.AGR[self.region_sismique] *si.m / si.s**2

    @property
    def a_g(self) -> tuple:
        """Calcule l'accélération de calcul a_g = a_gR × gamma_I selon EN 1998-1 §3.2.1.

        Returns:
            tuple: (latex_string, valeur) où valeur est l'accélération de calcul en m/s².
        """
        a_gr = self.a_gr
        gamma_I = self.gamma_I
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            a_g = a_gr * gamma_I
            return a_g
        return val()

    @property
    def gamma_I(self) -> float:
        """Retourne le coefficient d'importance gamma_I de la catégorie d'importance selon EN 1998-1 §4.2.5.

        Returns:
            float: Coefficient gamma_I (ex. : 0.8, 1.0, 1.2, 1.4 pour cat. I à IV).
        """
        return self.cat_importance_table["gamma_I"]
    

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
        ) -> None:
        """Affiche le spectre élastique de calcul S_d(T) selon EN 1998-1 §3.2.2.5.

        Trace la courbe S_d(T) sur [0, 4 s] pour la direction choisie avec matplotlib.

        Args:
            direction (str): Direction du spectre. "x" (bleu) ou "y" (rouge). Defaults to "x".
            screenshot (bool): Si True, enregistre le graphique au lieu de l'afficher.
                Defaults to False.
            filepath (str, optional): Chemin du fichier image de sortie (.png).
                Si None et screenshot=True, une boîte de dialogue s'ouvre.
                Defaults to None.
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
                from PySide6.QtWidgets import QFileDialog
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
            plt.savefig(filepath)
            return filepath
        else:
            plt.show()

    @property
    def T1(self) -> tuple:
        """Calcule les périodes propres de vibration T1 selon EN 1998-1 §4.3.3.2.2.

        Formule : T_1 = 2 * sqrt(m_total / K_b) pour chaque direction.
        La masse totale est déduite des charges gravitaires ajoutées via add_gravity_load.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": T_1_x, "y": T_1_y}
                avec les périodes en secondes (sans unité forallpeople).
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
    def Sd_t(self) -> tuple:
        """Retourne les ordonnées spectrales de calcul S_d(T1) selon EN 1998-1 §3.2.2.5.

        Calcule S_d(T1) pour les périodes T1x et T1y dans chaque direction en appelant
        _spectre_elastique_calcul, avec le coefficient de comportement q propre à chaque direction.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": Sd_x, "y": Sd_y}
                avec les accélérations spectrales en m/s².
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
    def Fb(self) -> tuple:
        """Calcule l'effort tranchant sismique à la base F_b selon EN 1998-1 §4.3.3.2.2 (éq. 4.5).

        F_b = S_d(T1) × m_total × lambda, avec lambda = 0.85 si T1 ≤ 2TC et plus de 2 niveaux, 1 sinon.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": F_b_x, "y": F_b_y}
                avec les efforts de base en kN (avec unité si.kN).
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
            
    def coeff_torsion_accidentelle(self, x: si.m, Le: si.m) -> tuple:
        """Calcule le coefficient de torsion accidentelle delta selon EN 1998-1 §4.3.3.2.4 (éq. 4.12).

        Applicable lorsque les raideurs latérales et les masses sont symétriques.
        Formule : delta = 1 + 1.2 × (x / Le).

        Args:
            x (float): Distance en plan de l'élément au centre de masse du bâtiment,
                mesurée perpendiculairement à la direction de l'action sismique, en m.
            Le (float): Distance entre les deux éléments de contreventement extrêmes,
                mesurée perpendiculairement à la direction de l'action sismique, en m.

        Returns:
            tuple: (latex_string, delta) où delta est le coefficient de torsion accidentelle (sans unité).
        """
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            delta = 1 + 1.2 * (x/Le) # form4.12
            return delta
        return val()

    def Fi(self, etage: str=Batiment.ETAGE) -> tuple:
        """Calcule l'effort horizontal équivalent F_i à l'étage i selon EN 1998-1 §4.3.3.2.3 (éq. 4.11).

        Formule : F_i = F_b × (z_i × m_i) / sum(z_j × m_j).
        Valable uniquement si les déplacements horizontaux croissent linéairement avec la hauteur.
        Le coefficient de torsion accidentelle n'est pas inclus (à appliquer séparément).

        Args:
            etage (str): Nom de l'étage considéré, tel que défini lors de add_gravity_load.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": F_i_x, "y": F_i_y}
                avec les efforts horizontaux en kN (avec unité si.kN).
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
    def Fi_table(self) -> 'pd.DataFrame':
        """Retourne le tableau des efforts horizontaux équivalents à chaque étage selon EN 1998-1 §4.3.3.2.3.

        Valable uniquement si les déplacements horizontaux croissent linéairement avec la hauteur.
        Le coefficient de torsion accidentelle n'est pas inclus.

        Returns:
            pd.DataFrame: DataFrame indexé par le nom d'étage avec les colonnes ["Fi,x", "Fi,y"]
                contenant les efforts en kN (avec unité si.kN).
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

    def ds(self, de: float, direction: str=("x", "y")) -> tuple:
        """Calcule le déplacement de calcul d_s selon EN 1998-1 §4.3.4 (éq. 4.23).

        Formule : d_s = q_d × d_e, avec q_d = q (coefficient de comportement de la direction).

        Args:
            de (float): Déplacement élastique du même point issu d'une analyse linéaire
                sur le spectre de calcul (EN 1998-1 §3.2.2.5), en mm.
            direction (str): Direction du déplacement. "x" ou "y". Defaults to "x".

        Returns:
            tuple: (latex_string, d_s) où d_s est le déplacement de calcul en mm (avec unité si.mm).
        """
        d_e = de * si.mm
        q_d = self.coeff_comportement[direction]["q"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            d_s = d_e * q_d # form4.23
            return d_s
        return val()

    def coeff_second_ordre(self, dr: float, V_tot: si.kN, etage: str=Batiment.ETAGE) -> tuple:
        """Calcule le coefficient d'effet P-delta theta_r selon EN 1998-1 §4.3.5.2.2 (éq. 4.28).

        Formule : theta_r = (P_tot × d_r) / (V_tot × h_lvl).
        Si theta_r ≤ 0.1 : coeff_P_delta = 1 (effets du 2nd ordre négligeables).
        Si 0.1 < theta_r ≤ 0.2 : coeff_P_delta = 1 / (1 - theta_r).
        Si theta_r > 0.2 : non couvert (restructuration requise).

        Args:
            dr (float): Déplacement relatif de calcul entre niveaux en mm (différence des
                déplacements latéraux entre le bas et le haut du niveau). Calculé via ds().
            V_tot (float): Effort tranchant sismique total au niveau de l'étage en kN.
                Inclure le coefficient de torsion accidentelle le cas échéant.
            etage (str): Nom de l'étage considéré, tel que défini lors de add_gravity_load.

        Returns:
            tuple: (latex_string, coeff_P_delta) où coeff_P_delta est le facteur d'amplification
                du 2nd ordre (sans unité).
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
    
    def taux_limitations_dommages(self, dr: float, etage: str=Batiment.ETAGE, type_dommages: str=TYPE_DOMMAGES) -> tuple:
        """Calcule le taux de limitation des dommages selon EN 1998-1 §4.4.3.2.

        Vérifie que le déplacement inter-étages d_r × nu reste dans la limite normative
        selon le type d'éléments non structuraux (éq. 4.31, 4.32 ou 4.33).
        Le coefficient de réduction nu = 0.4 est fixe (sismicité modérée AN français).

        Args:
            dr (float): Déplacement relatif de calcul entre niveaux en mm. Calculé via ds().
            etage (str): Nom de l'étage considéré, tel que défini lors de add_gravity_load.
            type_dommages (str): Type d'éléments non structuraux :
                "Fragiles" : éléments en matériaux fragiles fixés à la structure (limite h/200, éq. 4.31).
                "Ductiles" : éléments non structuraux ductiles (limite h/133, éq. 4.32).
                "Désolidarisé/Pas de risque" : éléments sans interaction avec la structure (limite h/100, éq. 4.33).

        Returns:
            tuple: (latex_string, taux) où taux est le taux de travail en limitation de dommages (sans unité,
                valeur ≤ 1 si la vérification est satisfaite).
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

    def Fa(
        self, 
        ma: float, 
        Ta_x: float, 
        Ta_y: float, 
        z: float, 
        type_element_ns: str=CLASSE_DUCTILITE_NS, 
        cat_importance_ns: str=CAT_IMPORTANCE_NS
        ) -> tuple:
        """Calcule l'effort sismique horizontal F_a sur un élément non structural selon EN 1998-1 §4.3.5.2 (éq. 4.24).

        Formule : F_a = (S_a × m_a × gamma_a) / q_a, avec S_a basé sur le spectre de plancher.

        Args:
            ma (float): Masse de l'élément non structural en kg ou kg/m² en fonction du besoin.
            Ta_x (float): Période fondamentale de vibration de l'élément suivant la direction x en s.
            Ta_y (float): Période fondamentale de vibration de l'élément suivant la direction y en s.
            z (float): Hauteur de l'élément non structural au-dessus du niveau d'application
                de l'action sismique (base du bâtiment ou sommet du soubassement) en m.
            type_element_ns (str): Type d'élément non structural selon CLASSE_DUCTILITE_NS.
                Détermine le coefficient de comportement q_a.
            cat_importance_ns (str): Catégorie d'importance de l'élément non structural selon
                CAT_IMPORTANCE_NS. Détermine le coefficient d'importance gamma_a.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": F_a_x, "y": F_a_y}
                avec les efforts sismiques en kN (avec unité si.kN).
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
        ) -> tuple:
        """Calcule l'effort sismique de dimensionnement des éléments non dissipatifs par capacité.

        Formule : F_s = F_i × gamma_d × Omega × eta_torsion × P_delta,
        avec gamma_d = 1.3 (rupture fragile) ou 1.1 (rupture ductile).

        Args:
            etage (str): Nom de l'étage considéré, tel que défini lors de add_gravity_load.
            gamma_d_x (str): Type de rupture dans la direction x :
                "Rupture fragile" (gamma_d = 1.3) : diaphragmes béton, embrèvements, assemblages collés,
                    connecteurs par plaques embouties, instabilité de flambement/déversement.
                "Rupture ductile" (gamma_d = 1.1) : diaphragmes en panneaux bois cloués,
                    assemblages par pointes ou tiges à faible raideur.
            gamma_d_y (str): Type de rupture dans la direction y. Mêmes valeurs que gamma_d_x.
            Omega_x (float): Facteur de sur-résistance des éléments dissipatifs dans la direction x
                (Rd/Ed). Defaults to 1.
            Omega_y (float): Facteur de sur-résistance dans la direction y. Defaults to 1.
            eta_torsion_x (float): Coefficient de torsion accidentelle dans la direction x,
                calculé via coeff_torsion_accidentelle(). Defaults to 1.
            eta_torsion_y (float): Coefficient de torsion accidentelle dans la direction y. Defaults to 1.
            P_delta_x (float): Coefficient d'effet du second ordre dans la direction x,
                calculé via coeff_second_ordre(). Defaults to 1.
            P_delta_y (float): Coefficient d'effet du second ordre dans la direction y. Defaults to 1.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": F_s_x, "y": F_s_y}
                avec les efforts sismiques de dimensionnement en kN (avec unité si.kN).
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
        ) -> tuple:
        """Calcule l'effort sismique de dimensionnement des éléments dissipatifs.

        Formule : F_s = F_i × eta_torsion × P_delta (sans amplification par capacité).
        Les éléments dissipatifs sont dimensionnés directement avec l'effort de l'analyse
        spectrale, sans majoration par gamma_d ni Omega.

        Args:
            etage (str): Nom de l'étage considéré, tel que défini lors de add_gravity_load.
            eta_torsion_x (float): Coefficient de torsion accidentelle dans la direction x,
                calculé via coeff_torsion_accidentelle(). Defaults to 1.
            eta_torsion_y (float): Coefficient de torsion accidentelle dans la direction y. Defaults to 1.
            P_delta_x (float): Coefficient d'effet du second ordre dans la direction x,
                calculé via coeff_second_ordre(). Defaults to 1.
            P_delta_y (float): Coefficient d'effet du second ordre dans la direction y. Defaults to 1.

        Returns:
            tuple: (latex_string, valeurs) où valeurs est un dict {"x": F_s_x, "y": F_s_y}
                avec les efforts sismiques de dimensionnement en kN (avec unité si.kN).
        """
        F_i_x = self.Fi_table.loc[etage]["Fi,x"]
        F_i_y = self.Fi_table.loc[etage]["Fi,y"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_s_fin_x = F_i_x * eta_torsion_x * P_delta_x
            F_s_fin_y = F_i_y * eta_torsion_y * P_delta_y
            return {"x": F_s_fin_x, "y": F_s_fin_y}
        return val()
        
