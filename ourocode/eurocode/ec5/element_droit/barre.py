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

from ourocode.eurocode.core.projet import Projet


class Barre(Projet):
    """Classe définissant les caractéristiques d'un élément droit en bois.

    Cette classe décrit la géométrie, la classe de résistance et les conditions
    d'exploitation d'une barre (poutre, colonne) selon l'EN 1995.

    Elle calcule automatiquement les dimensions de section en fonction de
    l'humidité de pose (retrait/gonflement) et donne accès aux caractéristiques
    mécaniques normatives du matériau.
    """

    LIST_SECTION = ["Rectangulaire","Circulaire"]
    LIST_TYPE_B = ["Massif", "BLC", "LVL", "OSB 2", "OSB 3/4", "CP"]
    CLASSE_WOOD = list(Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").index)[2:]
    CLASSE_PANEL = list(Projet._data_from_csv(Projet, "caracteristique_meca_panel.csv").index)[2:]
    CLASSE = CLASSE_WOOD + CLASSE_PANEL
    CS = ["1","2","3"]
    CARACTERISTIQUE = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").columns)
    LOAD_TIME = tuple(Projet._data_from_csv(Projet, "kmod.csv").columns)[1:]
    TYPE_ACTION = ("Fondamentales" ,"Accidentelles")
    TYPE_BAT = ("Bâtiments courants", "Bâtiments agricoles et similaires")
    TYPE_ELE = tuple(Projet._data_from_csv(Projet, "limite_fleche.csv").index.unique())
    B90 = 0.25

    def __init__(self, b:si.mm, h:si.mm, section: str=LIST_SECTION, Hi: int=12, Hf: int=12, classe: str=CLASSE, cs: int=CS, effet_systeme: bool=("False", "True"), **kwargs):
        """Initialise un élément droit en bois avec ses caractéristiques.

        Args:
            b (si.mm): Largeur de pose de la pièce en millimètres (dimension brute).
            h (si.mm): Hauteur de pose de la pièce en millimètres (dimension brute).
            section (str, optional): Type de section transversale.
                "Rectangulaire" ou "Circulaire". Defaults to "Rectangulaire".
            Hi (int, optional): Humidité initiale de pose en %.
                Humidité au moment de la fabrication/pose. Defaults to 12.
            Hf (int, optional): Humidité finale d'équilibre en % selon l'AN (Hf = 12).
                Humidité en service. Defaults to 12.
            classe (str, optional): Classe de résistance du bois selon l'EC5.
                Ex: "C24", "GL28h", "LVL". Defaults to "C24".
            cs (int, optional): Classe de service selon l'EC5 §2.3.1.3.
                1 = intérieur chauffé, 2 = couvert non chauffé, 3 = extérieur.
                Defaults to 1.
            effet_systeme (bool, optional): Active l'effet de système (k_sys = 1.1)
                pour les éléments permettant une redistribution des charges
                (solives avec répartition continue). Defaults to False.
            **kwargs: Arguments supplémentaires transmis à Projet.

        Note:
            Les dimensions de calcul (b_calcul, h_calcul) sont automatiquement
            ajustées pour tenir compte du retrait si Hi > Hf (AN B90 = 0.25%).

        Raises:
            ValueError: Si la section, la classe ou la classe de service est invalide.
        """
        super().__init__(**kwargs)
        if section not in self.LIST_SECTION:
            raise ValueError(
                f"Section '{section}' invalide. Valeurs acceptées : {self.LIST_SECTION}"
            )
        if classe not in self.CLASSE:
            raise ValueError(
                f"Classe '{classe}' invalide. Valeurs acceptées : {self.CLASSE}"
            )
        if int(cs) not in (1, 2, 3):
            raise ValueError(
                f"Classe de service '{cs}' invalide. Valeurs acceptées : 1, 2, 3"
            )
        self.b = b * si.mm
        self.h = h * si.mm
        self.section = section
        self.Hi = Hi
        self.Hf = Hf
        self.classe = classe
        self.cs = int(cs)
        self.effet_systeme = effet_systeme
        self._sectionCalcul()


    def _sectionCalcul(self):
        """Calcule les dimensions de section de calcul avec correction humidité.

        Ajuste les dimensions brutes (b, h) en fonction de la variation
        d'humidité entre la pose et l'équilibre en service selon la formule
        de l'Annexe Nationale française (coefficient B90 = 0.25 %).

        Formule : dimension_calcul = dimension_pose × (1 - B90/100 × (Hi - Hf))

        Cette correction est appliquée pour le calcul des aires et inerties,
        mais les résistances caractéristiques sont basées sur les dimensions
        de pose selon les règles de l'EC5.

        Attributs modifiés:
            b_calcul (si.mm): Largeur corrigée pour le calcul des inerties.
            h_calcul (si.mm): Hauteur corrigée pour le calcul des inerties.
        """
        self.b_calcul = self.b * (1 - self.B90 / 100 * (self.Hi - self.Hf))
        self.h_calcul = self.h * (1 - self.B90 / 100 * (self.Hi - self.Hf))

    
    @property
    def aire(self):
        if self.section == self.LIST_SECTION[0]:
            return self.b_calcul * self.h_calcul
        else:
            return mt.pi * (self.b_calcul/2)**2
        

    @property
    def inertie(self):
        """ Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
            b ou d : Largeur ou diamètre de la poutre en mm
            h : Hauteur de la poutre en mm """
        if self.section == "Rectangulaire":
            self.I_y = (self.b_calcul * self.h_calcul**3)/12
            self.I_z = (self.h_calcul * self.b_calcul**3)/12
            return [self.I_y, self.I_z]

        elif hasattr(self, "Iy") and hasattr(self, "Iz"):
            return [self.I_y * si.mm**4, self.I_z * si.mm**4]

        else:
            self.I_y = (mt.pi * self.b_calcul ** 4) / 64
            self.I_z = self.I_y
            return [self.I_y, self.I_z]
        
    
    @property
    def caract_meca(self):
        """ Retourne les caractéristiques méca du bois sous forme de dataframe pandas """
        if self.classe in self.CLASSE_WOOD:
            data_csv_meca = self._data_from_csv("caracteristique_meca_bois.csv")
            return data_csv_meca.loc[self.classe]
        elif self.classe in self.CLASSE_PANEL:
            data_csv_meca = self._data_from_csv("caracteristique_meca_panel.csv")
            return data_csv_meca.loc[self.classe]
    

    @property
    def gamma_M_table(self):
        """Retourne le tableau des gamma M pour le type de bois sélectionné
        """
        data_csv_gammaM = self._data_from_csv("gammaM.csv")
        return data_csv_gammaM.loc[self.type_bois]
    
    
    def _get_gamma_M(self, typecombi=TYPE_ACTION):
        self.gamma_M = self.gamma_M_table.loc[typecombi]
        return self.gamma_M
    

    @property
    def K_def(self):
        data_csv_kdef = self._data_from_csv("kdef.csv")
        kdef = float(data_csv_kdef.loc[self.type_bois][str(self.cs)])
        if self.Hi > 20:
            kdef += 1
        return kdef
    
    
    @property
    def K_mod_table(self):
        """ Retourne le tableau des Kmod du bois
        """
        data_csv_kmod = self._data_from_csv("kmod.csv")
        data_kmod = data_csv_kmod.loc[self.type_bois]
        return data_kmod.loc[data_kmod["CS"]==self.cs]
    

    def _get_k_mod(self, loadtype=LOAD_TIME):
        self.K_mod = self.K_mod_table[loadtype].iloc[0]
        return self.K_mod
    
    
    @property
    def k_sys(self):
        """Détermine le Ksys d'un élément si celui-ci permet une redistribution des charges continues.
        """
        if self.effet_systeme:
            return 1.1
        else:
            return 1

    
    @property
    def type_bois(self):
        if self.classe[0:1] == "C" or self.classe[0:1] == "D":
            type_b = self.LIST_TYPE_B[0]
        elif self.classe[0:2] == "GL":
            type_b = self.LIST_TYPE_B[1]
        elif self.classe[0:3] == "LVL":
            type_b = self.LIST_TYPE_B[2]
        elif self.classe[0:5] == "OSB/2":
            type_b = self.LIST_TYPE_B[3]
        elif self.classe[0:5] == "OSB/3" or self.classe[0:5] == "OSB/4":
            type_b = self.LIST_TYPE_B[4]
        else:
            type_b = self.LIST_TYPE_B[5]
        
        return type_b


    def __convert_latex_ftyped(self, latex: str, type_caract: str):
        end_index_rk = type_caract.find("k")
        type_caract = type_caract[1:end_index_rk]
        latex = latex.replace("f_{type_{Rd}", "f_{"+type_caract+"_{"+"d}")
        latex = latex.replace("f_{type_{k}", "f_{"+type_caract+"_{"+"k}")
        return latex
    
    
    def _f_type_d(self, typeCarac=CARACTERISTIQUE[0:6], loadtype=LOAD_TIME, typecombi=TYPE_ACTION):
        """Calcule la résistance de calcul f_d selon l'EC5 §2.4.1.

        Applique la formule : f_d = k_mod × k_sys × f_k / gamma_M
        où k_mod dépend de la classe de service et de la durée de chargement,
        et k_sys est le coefficient d'effet de système (1.1 si activé, 1.0 sinon).

        Args:
            typeCarac (str, optional): Type de résistance caractéristique.
                Valeurs possibles : "fm0k" (flexion), "fc0k" (compression),
                "ft0k" (traction), "fvk" (cisaillement), etc.
                Defaults to "fm0k".
            loadtype (str, optional): Durée de chargement selon l'EC5 Tableau 3.1.
                Valeurs : "Permanente", "Long terme", "Moyen terme", "Court terme",
                "Instantanée".
            typecombi (str, optional): Type de combinaison.
                "Fondamentales" ou "Accidentelles". Defaults to "Fondamentales".

        Returns:
            tuple: (latex_string, valeur) où valeur est la résistance de calcul
                f_d en N/mm² (MPa) prête à être comparée aux contraintes.

        Note:
            Cette méthode utilise @handcalc pour générer la justification LaTeX.
        """
        gamma_M = self._get_gamma_M(typecombi)
        K_mod = self._get_k_mod(loadtype)
        f_type_k = float(self.caract_meca.loc[typeCarac]) * si.MPa

        if typeCarac == "fm0k" and self.k_sys > 1:
            k_sys = self.k_sys
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_type_Rd = k_sys * f_type_k * K_mod / gamma_M
                return f_type_Rd
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_type_Rd = f_type_k * K_mod / gamma_M
                return f_type_Rd
        value = val()
        latex = self.__convert_latex_ftyped(value[0], typeCarac)
        self.f_type_rd = value[1]
        return (latex, value[1])
    

    def _K_h(self):
        """Calcule le coefficient de taille K_h selon l'EC5 §3.2 et §3.3.

        Ce coefficient majore la résistance caractéristique à la flexion et
        à la traction pour les petites sections, selon le type de bois :

        - Bois massif : K_h = min((150/h)^0.2, 1.3) où h < 150 mm
        - BLC : K_h = min((600/h)^0.1, 1.1) où h < 600 mm
        - LVL : K_h = 1.0 (non applicable)

        Le coefficient s'applique séparément aux deux dimensions (hauteur
        et largeur) pour les sections rectangulaires.

        Returns:
            dict: Dictionnaire {'y': K_h_y, 'z': K_h_z} avec les coefficients
                pour chaque direction de flexion.

        Note:
            Ce coefficient ne s'applique qu'à fm,k et ft,k selon l'EC5.
        """
        kh = {}
        dim = {'y': self.h_calcul.value *10**3, 'z': self.b_calcul.value *10**3}

        for cle, valeur in dim.items():
            if self.type_bois == "Massif":
                if valeur < 150:
                    kh[cle] = min((150 / valeur) ** 0.2, 1.3)
                else :
                    kh[cle] = 1
            elif self.type_bois == "BLC":
                if valeur < 600:
                    kh[cle] = min((600 / valeur) ** 0.1, 1.1)
                else :
                    kh[cle] = 1
            else:
                warnings.warn("LVL non pris en compte dans cette fonction")
                kh[cle] = 1
        return kh       
    
    
    def Emean_fin(self, psy_2: float):
        """Calcule le module de Young final E_mean,fin selon l'EC5 §2.3.2.2.

        Le module final tient compte du fluage par la formule :
        E_mean,fin = E_0,mean / (1 + psi_2 × k_def)

        où psi_2 est le coefficient de combinaison quasi-permanente et
        k_def dépend de la classe de service et du type de bois.

        Args:
            psy_2 (float): Coefficient psi_2 de la combinaison quasi-permanente.
                0 pour le court terme, 1 pour le long terme, ou valeur calculée.

        Returns:
            tuple: (latex_string, valeur) où valeur est E_mean,fin en MPa.

        Note:
            Ce module final est utilisé pour les calculs de flèche en ELS
            selon l'EC5 §2.2.3 et §7.2.
        """
        self.psy_2 = psy_2
        E0_mean = int(self.caract_meca.loc["E0mean"]) * si.MPa
        K_def = self.K_def

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            E_mean_fin = E0_mean / (1 + psy_2 * K_def)
            return E_mean_fin
            
        value = val()
        self.E_mean_fin = value[1]
        return value   
    

    def fleche(self, long:si.mm, Ed_WinstQ:si.mm=0, Ed_Wnetfin:si.mm=0, Ed_Wfin:si.mm=0, Ed_W2:si.mm=0, limit_W2:int=500, type_ele=TYPE_ELE, type_bat=TYPE_BAT):
        """Vérifie les taux de travail des flèches selon l'EC5 §7.2.

        Compare les flèches calculées aux limites normatives pour différents
        critères ELS. Génère automatiquement un tableau de synthèse avec
        _add_synthese_taux_travail.

        Args:
            long (si.mm): Portée entre appuis à vérifier en millimètres.
            Ed_WinstQ (si.mm, optional): Flèche instantanée sous charge variable Q seule.
                Defaults to 0.
            Ed_Wnetfin (si.mm, optional): Flèche nette finale (sous combinaison quasi-permanente).
                Defaults to 0.
            Ed_Wfin (si.mm, optional): Flèche finale totale (y compris fluage).
                Defaults to 0.
            Ed_W2 (si.mm, optional): Flèche w2 tenant compte du phasage de pose
                pour éléments fragiles (cloisons, platrerie). Defaults to 0.
            limit_W2 (int, optional): Limite de flèche w2 pour éléments fragiles.
                Valeur courante : 500 (L/500). Defaults to 500.
            type_ele (str, optional): Type d'élément selon limite_fleche.csv.
            type_bat (str, optional): Type de bâtiment.

        Returns:
            tuple: (latex_string, valeurs) où valeurs contient les taux de travail.

        Note:
            Les limites de flèche sont définies dans le fichier limite_fleche.csv
            selon les recommandations de l'Annexe Nationale française.
        """
        data_csv_fleche = self._data_from_csv("limite_fleche.csv")
        self.data_fleche= data_csv_fleche.loc[type_ele]
        self.data_fleche = self.data_fleche.loc[self.data_fleche["Type bâtiment"]==type_bat]
        self.taux_ELS = {}

        long = long * si.mm
        Ed_W_inst_Q = Ed_WinstQ * si.mm
        Ed_W_net_fin = Ed_Wnetfin * si.mm
        Ed_W_fin = Ed_Wfin * si.mm
        Ed_W2 = Ed_W2 * si.mm
        limit_W2 = int(limit_W2)

        limit_W_inst_Q = self.data_fleche['Winst(Q)'].iloc[0]
        limit_W_net_fin = int(self.data_fleche['Wnet,fin'].iloc[0])
        limit_W_fin = int(self.data_fleche['Wfin'].iloc[0])
        limit_U_fin_max = self.data_fleche['Ufin,max'].iloc[0]
        
        if np.isnan(limit_U_fin_max):
            limit_U_fin_max = long / limit_W_fin
        else:
            limit_U_fin_max = int(limit_U_fin_max)

        if not np.isnan(limit_W_inst_Q):
            limit_W_inst_Q = int(limit_W_inst_Q)
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                Rd_W_inst_Q = long / limit_W_inst_Q
                Rd_W_net_fin = long / limit_W_net_fin
                Rd_W_fin = min(long / limit_W_fin, limit_U_fin_max)
                Rd_W2 = long / limit_W2

                taux_W_inst_Q = Ed_W_inst_Q / Rd_W_inst_Q
                taux_W_net_fin = Ed_W_net_fin / Rd_W_net_fin
                taux_W_fin = Ed_W_fin / Rd_W_fin
                taux_W2 = Ed_W2 / Rd_W2
                return taux_W_inst_Q, taux_W_net_fin, taux_W_fin, taux_W2
            
            value = val()
            self.taux_ELS["Winst(Q)"] = value[1][0]
            self.taux_ELS["Wnet,fin"] = value[1][1]
            self.taux_ELS["Wfin"] = value[1][2]
            self.taux_ELS["W2"] = value[1][3]
            synthese = [
                ["Flèche W,inst(Q)", self.taux_ELS["Winst(Q)"], None],
                ["Flèche W,net,fin", self.taux_ELS["Wnet,fin"], None],
                ["Flèche W,fin", self.taux_ELS["Wfin"], None],
                ["Flèche W2", self.taux_ELS["W2"], None],
            ]
            self._add_synthese_taux_travail(synthese)

        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                Rd_W_net_fin = long / limit_W_net_fin
                Rd_W_fin = min(long / limit_W_fin, limit_U_fin_max)
                Rd_W2 = long / limit_W2

                taux_W_net_fin = Ed_W_net_fin / Rd_W_net_fin
                taux_W_fin = Ed_W_fin / Rd_W_fin
                taux_W2 = Ed_W2 / Rd_W2
                return taux_W_net_fin, taux_W_fin, taux_W2
            
            value = val()
            self.taux_ELS["Wnet,fin"] = value[1][0]
            self.taux_ELS["Wfin"] = value[1][1]
            self.taux_ELS["W2"] = value[1][2]
            synthese = [
                ["Flèche W,net,fin", self.taux_ELS["Wnet,fin"], None],
                ["Flèche W,fin", self.taux_ELS["Wfin"], None],
                ["Flèche W2", self.taux_ELS["W2"], None],
            ]
            self._add_synthese_taux_travail(synthese)
        return value 


