#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
from copy import deepcopy
import warnings
import matplotlib.pyplot as plt

import math as mt
from math import sqrt, pi, cos, sin, radians
import numpy as np

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.A0_Projet import Projet


# ================================ GLOBAL ==================================

class Barre(Projet):
    """Classe définissant les caractéristiques d'un élément droit en bois.

    Cette classe décrit la géométrie, la classe de résistance et les conditions
    d'exploitation d'une barre (poutre, colonne, liteau) selon l'EN 1995.

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


# ================================ FLEXION ==================================

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

class Traction(Barre):
    """Classe de vérification des éléments bois en traction axiale selon l'EC5 §6.1.2.

    Effectue les calculs de résistance et de contrainte en traction axiale selon
    l'Eurocode 5 - Partie 1-1. Hérite de Barre pour les caractéristiques
    géométriques et mécaniques.

    La vérification principale est le taux de travail en traction (équation 6.1):
    σ_t,0,d / (f_t,0,d · k_h) ≤ 1
    """

    def __init__(self, *args, **kwargs):
        """Initialise un objet de vérification en traction axiale.

        Hérite de toutes les caractéristiques de Barre (section, classe de bois,
        etc.). Aucun paramètre supplémentaire requis à l'initialisation.

        Args:
            *args: Arguments positionnels transmis à Barre.
            **kwargs: Arguments nommés transmis à Barre (b, h, classe, etc.).
        """
        super().__init__(*args, **kwargs)


    @property
    def K_h(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        return self._K_h()


    def f_t_0_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul en traction axiale f_t,0,d selon l'EC5 §6.1.2.

        Détermine la résistance à partir de la résistance caractéristique ft,0,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul ft,0,d en MPa avec unité (si.MPa).
        """
        return super()._f_type_d("ft0k", loadtype, typecombi)


    def sigma_t_0_d(self, Ft0d: si.kN, Anet: si.mm**2=None):
        """Calcule la contrainte de traction axiale sigma_t,0,d selon l'EC5 §6.1.2.

        Détermine la contrainte normale due à l'effort de traction axial.
        Prend en compte une section nette réduite (perçages, entailles) si spécifiée.

        Args:
            Ft0d (si.kN): Effort de traction axial en kN.
            Anet (si.mm**2, optional): Aire nette de la section en mm² si réduction
                (perçages, entailles). Doit être ≤ aire brute. Defaults to None.

        Returns:
            tuple: (latex_string, valeur) où valeur est sigma_t,0,d en MPa avec unité.

        Raises:
            ValueError: Si Anet > aire brute de la section.

        Note:
            La valeur est stockée dans l'attribut sigma_t_0_rd.
            Pour les assemblages boulonnés, utiliser Anet pour tenir compte des trous.
        """
        self.Ft_0_d = Ft0d * si.kN
        Ft_0_d = self.Ft_0_d
        if Anet and Anet * si.mm**2<= self.aire:
            A = Anet * si.mm**2
        else:
            A = self.aire

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_t_0_d = Ft_0_d / A
            return sigma_t_0_d
        value = val()
        self.sigma_t_0_rd = value[1]
        return value


    def taux_t_0_d(self):
        """Calcule le taux de travail en traction axiale selon l'EC5 §6.1.2 (Eq. 6.1).

        Vérifie le critère : σ_t,0,d / (k_h · f_t,0,d) ≤ 1

        Le coefficient k_h (effet de hauteur) est pris comme le minimum des
        valeurs selon y et z pour être conservateur.

        Returns:
            tuple: (latex_string, valeur) où valeur est le taux en pourcentage
                (0.75 = 75%). Stocké dans taux_t_0_rd['equ6.1'].

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_t_0_rd = {}
        K_h_y = self.K_h['y']
        K_h_z = self.K_h['z']
        sigma_t_0_d = self.sigma_t_0_rd
        f_t_0_d = self.f_type_rd

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            K_h = min(K_h_y, K_h_z)
            taux_6_1 = sigma_t_0_d / (K_h * f_t_0_d) # equ6.1
            return taux_6_1
        value = val()

        self.taux_t_0_rd['equ6.1'] = value[1]
        synthese = [
            ["Traction bois", self.taux_t_0_rd['equ6.1'], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value


# ================================ Compression ==================================

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
    def lamb(self):
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
    def lamb_rel_Axe(self):
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
    def beta_C(self):
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
    def kc_Axe(self):
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
    
    
    def sigma_c_0_d(self, Fc0d: si.kN, Anet: si.mm**2=None):
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
    
    
    def taux_c_0_d(self, flexion: object=None):
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

        if flexion and isinstance(flexion, Flexion):
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
    
    
    def sigma_c_90_d(self, Fc90d: si.kN):
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


    def taux_c_90_d(self):
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
    
    
    def show_c90(self):
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
    
    def sigma_c_alpha_d(self, Fcad: si.kN):
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
    

    def taux_c_alpha_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
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
    


# ================================ CISAILLEMENT ==================================

class Cisaillement(Barre):
    """Classe de vérification des éléments bois au cisaillement selon l'EC5 §6.1.7 et §6.5.

    Effectue les calculs de contrainte et taux de travail au cisaillement
    longitudinal pour des poutres en bois selon l'Eurocode 5.
    Hérite de Barre pour les caractéristiques géométriques et mécaniques.

    Vérifie :
    - La résistance au cisaillement (§6.1.7, Eq. 6.13)
    - Le cisaillement avec entaille (§6.5, Eq. 6.60) avec facteur K_v
    """

    DICT_KN = {"Massif": 5,"BLC": 6.5, "LVL": 4.5}
    def __init__(self, **kwargs):
        """Initialise un objet de vérification au cisaillement.

        Hérite de toutes les caractéristiques de Barre. Initialise K_v à 1
        (pas d'entaille) et h_ef à la hauteur totale.

        Args:
            **kwargs: Arguments transmis à Barre (classe, etc.).
        """
        super().__init__(**kwargs)
        self.K_v = 1
        self.h_ef = self.h_calcul

    @property
    def K_cr(self):
        """Facteur de réduction de largeur K_cr selon l'EC5 §6.1.7.

        Tient compte des fissures de séchage dans le bois massif et BLC.
        Réduction de 33% (K_cr = 0.67) pour les sections fragilisées.

        Returns:
            float: Valeur de K_cr.
                - 1.0 : Pas de réduction (bois sans fissuration significative)
                - 0.67 : Réduction pour bois massif h > 150mm (CS1/CS2) ou BLC (CS2/CS3)

        Note:
            CS1/CS2 : K_cr = 0.67 si h > 150mm et bois massif
            CS2 : K_cr = 0.67 pour BLC
            CS3 : K_cr = 0.67 pour tout bois
        """
        if self.cs == 1:
            if self.h_calcul.value * 10**3 > 150 and self.type_bois == "Massif":
                return 0.67
            else:
                return 1
        elif self.cs == 2:
            if self.h_calcul.value * 10**3 > 150 and self.type_bois == "Massif":
                return 0.67
            elif self.type_bois == "BLC":
                return 0.67
            else:
                return 1
        else:
            return 0.67


    def Kv(self, hef:si.mm, x:si.mm, i_lo:si.mm, ent=("Dessous", "Dessus")):
        """Calcule le facteur de réduction d'entaille K_v selon l'EC5 §6.5.

        Ce coefficient réduit la résistance au cisaillement en présence d'une
        entaille au niveau d'un appui (entaille en dessous ou au dessus de la poutre).

        Formule EC5 : K_v = min(1, [K_n(1 + 1.1·i^1.5/√h)] / [√h·(√(α(1-α)) + 0.8·x/h·√(1/α - α²))])
        où α = h_ef/h et i = i_lo/h_ef

        Args:
            hef (si.mm): Hauteur efficace de la poutre (hauteur - profondeur entaille) en mm.
            x (si.mm): Distance entre le centre de réaction à l'appui et le coin de l'entaille en mm.
            i_lo (si.mm): Longueur horizontale de l'entaille en mm.
            ent (str, optional): Position de l'entaille : "Dessous" ou "Dessus".
                Defaults to "Dessous".

        Returns:
            tuple: (latex_string, valeur) où valeur est le facteur K_v (≤ 1).
                Si ent="Dessus", retourne K_v = 1 (pas de réduction).

        Note:
            La valeur est stockée dans l'attribut K_v.
            K_n dépend du type de bois (voir DICT_KN).
        """
        K_n = self.DICT_KN[self.type_bois]
        x = x * si.mm
        h_ef = hef * si.mm
        i = i_lo * si.mm / h_ef
        h_calcul = self.h_calcul

        self.h_ef = h_ef
        if ent == "Dessus":
            self.K_v = 1
            return self.K_v
        else:
            @handcalc(override="long", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                alpha = h_ef / h_calcul
                K_v = min(1,(K_n * (1 + (1.1 * i ** 1.5) / sqrt(h_calcul))) / (sqrt(h_calcul) * (sqrt(alpha * (1 - alpha)) + 0.8 * x / h_calcul * sqrt(1 / alpha - alpha ** 2))))
                return K_v
            value = val()
            self.K_v = value[1]
            return value

    def f_v_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul au cisaillement f_v,d selon l'EC5 §6.1.7.

        Détermine la résistance à partir de la résistance caractéristique f_v,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul f_v,d en MPa avec unité (si.MPa).
        """
        return super()._f_type_d("fvk", loadtype, typecombi)

    def tau_d(self, Vd:si.kN):
        """Calcule la contrainte de cisaillement tau_d selon l'EC5 §6.1.7.

        Détermine la contrainte de cisaillement longitudinal pour une poutre
        rectangulaire : τ = 1.5·V/(b_ef·h_ef)

        Args:
            Vd (si.kN): Effort tranchant (cisaillement) sur la poutre en kN.

        Returns:
            tuple: (latex_string, valeur) où valeur est tau_d en MPa avec unité.

        Note:
            La valeur est stockée dans l'attribut tau_rd.
            Prend en compte K_cr (fissuration) et K_v (entaille si défini).
        """
        self.V_d = Vd * si.kN
        V_d = self.V_d
        K_cr = self.K_cr
        b_calcul = self.b_calcul
        h_ef = self.h_ef

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            b_ef = K_cr * b_calcul
            tau_d = (1.5 * V_d) / (b_ef * h_ef)
            return tau_d

        value = val()
        self.tau_rd = value[1]
        return value

    def taux_tau_d(self):
        """Calcule les taux de travail au cisaillement selon l'EC5 §6.1.7 et §6.5.

        Vérifie les critères :
        - Equ. 6.13 : Cisaillement simple (tau_d / f_v,d)
        - Equ. 6.60 : Cisaillement avec entaille (tau_d / (K_v · f_v,d))

        Returns:
            tuple: (latex_string, valeurs) où valeurs contient les taux pour
                equ6.13 (sans entaille) et equ6.60 (avec entaille) en pourcentage.
                Stockés dans taux_tau_rd['equ6.13'] et ['equ6.60'].

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_tau_rd = {}
        tau_d = self.tau_rd
        f_v_d = self.f_type_rd
        K_v  = self.K_v

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_6_13 = tau_d / f_v_d
            taux_6_60 = tau_d / (K_v * f_v_d)
            return taux_6_13, taux_6_60

        value = val()
        self.taux_tau_rd['equ6.13'] = value[1][0]
        self.taux_tau_rd['equ6.60'] = value[1][1]
        synthese = [
            ["Cisaillement bois", max(value[1][0], value[1][1]), None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value
    
    
    def show_Kv(self):
        """Affiche l'image des caractéristiques d'une entaille au cisaillement
        """
        self._show_element("Kv_def.png")


# ================================ Barre assemblées mécaniquement Annexe B ==================================

# class Poutre_assemblee_meca(Projet):
#     def __init__(self, beam_2:object, l: si.mm, disposition: str=["Latérale", "Dessus / Dessous"], recouvrement: list=[0,0], Ki: list=[0,None,0], entraxe: list=[1, None, 1], psy_2: int|float=0, **kwargs):
#         """Classe définitssant une poutre composée d'au maximum 3 éléments connectés entre eux par liaisons mécanique 
#         suivant la théorie de HEIMESHOFF Annexe B de l'EN 1995
#         Args:
#             beam_2 (object): objet contenant l'élément centrale de la poutre i=2, cette objet doit être issu de la classe élément ou de son héritage
#             l (int | float): longueur de la poutre en mm
#             disposition (str): Disposition des éléments supplémentaires à l'élément 2, "Latérale" ou "Dessus / Dessous"
#             recouvrement (int): Si disposition latérale alors donner un recouvrement en mm. 
#                                 Cela correspond à la distance entre le centre géométrique de la pièce 2 et celui de la pièce i.
#                                 ATTENTION z local est vers le bas quand vous donnez le recouvrement.
                                
#             Ki (list, optional): Rigidité de liaison par connecteur entre les éléments entre i=1/2 et i=2/3, en N/mm. Soit Kser soit Ku en fonction du type de vérification
#                                     S'il n'y a que 2 éléments connectés, laisser l'index correspondant vide (ex: [0, 2000]). Defaults to [0,None,0].
                                    
#             entraxe (list, optional): Entraxe entre connecteur en mm suivant i=1 ou i=3. Defaults to [1, None, 1].
#             type_action (int | float, optional): Psy 2 qui permet de prendre en compte le fluage dans le temps,
#                                         si calcul en instantanée alors 0, 
#                                         si intermédiaire = voir dans data/coeff_psy.csv 
#                                         et enfin temps infini = 1. 
#                                         Defaults to 0.
                                        
#             **kwargs (object): beam_1 et ou beam_3
#         """
#         super().__init__(**kwargs)
#         self.beam = [None , beam_2, None]
#         self.l = l * si.mm
#         self.disposition = disposition
#         self.recouvrement = [recouvrement[0]*si.mm, recouvrement[1]*si.mm]
#         self.entraxe = []
#         self.Ki = []
#         for i in range(3):
#             if entraxe[i]:
#                 self.entraxe.append(entraxe[i] * si.mm)
#             else:
#                 self.entraxe.append(None)
#             if Ki[i]:
#                 self.Ki.append(Ki[i] * si.N / si.mm)
#             else:
#                 self.Ki.append(None)
            
#         for key, value in kwargs.items():
#             match key[0:4]:
#                 case "beam" if key[-1] == "2":
#                     print("L'attribut ne peut pas être nommé beam_2, il est déjà pris par l'élément centrale ! beam_1 ou beam_3 possible !")
#                 case "beam" if key[-1] != "2" and 1<= int(key[-1]) <=3 :
#                     self.beam[int(key[-1])-1] = value
#         for index, beam in enumerate(self.beam):
#             if beam is not None:
#                 beam.Emean_fin(psy_2)

   
#     @property
#     def K_def(self):
#         kdef = (None, 0)
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 K_def_i = beam.K_def
#                 K_def_2 = self.beam[1].K_def
#                 if K_def_i != K_def_2:
#                     @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                     def val():
#                         K_def = 2 * sqrt(K_def_i * K_def_2)
#                         return K_def
#                 else:
#                     @handcalc(override="params", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                     def val():
#                         K_def= K_def_2
#                         return K_def
#                 value = val()
#                 if value[1] > kdef[1]:
#                     kdef = value
#         return kdef
    
        
#     @property
#     def Ki_fin(self):
#         """Renvoie le Ki en fonction des Kdef des pièces assemblées et du psy2"""
#         ki_fin = {}
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 K_i = self.Ki[index]
#                 psy_2 = beam.psy_2
#                 K_def = self.K_def[1]
                
#                 @handcalc(override="long", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     K_i_fin = K_i / (1 + psy_2 * K_def)
#                     return K_i_fin
                
#                 if index == 0:
#                     ki_fin["Ki fin 1-2"] = val()
#                 else:
#                     ki_fin["Ki fin 2-3"] = val()
#         return ki_fin
    
    
#     @property
#     def gamma_i(self):
#         """Renvoie le gamma"""
#         gamma = [0, 1, 0]
#         gamma = {"gamma 1": None, "gamma 2": 1, "gamma 3": None}
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 E_mean_fin = beam.E_mean_fin.value * 10 ** -6
#                 A = beam.aire.value * 10 ** 6
#                 entraxe = self.entraxe[index].value*10**3
#                 lo = self.l.value * 10 ** 3
#                 if index == 0:
#                    K_i_fin = self.Ki_fin["Ki fin 1-2"][1].value * 10**-3
#                 else:
#                    K_i_fin = self.Ki_fin["Ki fin 2-3"][1].value * 10**-3
                
#                 @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     gamma_i = (1 + pi ** 2 * E_mean_fin * A * entraxe / (K_i_fin * lo ** 2)) ** (-1)
#                     return gamma_i
#                 gamma["gamma "+str(index+1)] = val()
#         return gamma

    
#     @property
#     def distance_ai(self):
#         """Renvoie la distance à l'axe neutre de la pièce 2"""
#         denominateur = 0
#         if self.disposition == "Latérale":
#             d1 = self.recouvrement[0]
#             d3 = self.recouvrement[1]
            
#         else:
#             h1 = self.beam[0].h_calcul
#             h2 = self.beam[1].h_calcul
#             h3 = self.beam[2].h_calcul
#             d1 = -(h1 + h2)/2
#             d3 = (h3 + h2)/2
            
#         for index, beam in enumerate(self.beam):
#             if beam is not None:
#                 if index == 1:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)]
#                 else:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)][1]
                    
#                 denominateur = denominateur + (gamma_i * beam.E_mean_fin * beam.aire)
                
#         if self.beam[0] == None or self.beam[2] == None:
#             di = [d1, None, d3]
#             for index, beam in enumerate(self.beam):
#                 if beam is not None and index != 1:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)][1]
#                     numerateur = gamma_i * beam.E_mean_fin * beam.aire * di[index]
#         else:
#             numerateur = (self.gamma_i["gamma 1"][1] * self.beam[0].E_mean_fin * self.beam[0].aire * d1 
#                         - self.gamma_i["gamma 3"][1] * self.beam[2].E_mean_fin  * self.beam[2].aire * d3)
#         a2 =  numerateur / denominateur
#         ai = [None, -a2, None]
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 di = [d1, None, d3]
#                 ai[index] = di[index] - a2
#         return ai
       
#     @property
#     def EI_eff (self):
#         """Renvoie la rigidité de la poutre connectée"""
#         EIeff_latex = ""
#         EIeff_value = 0
#         for index, beam in enumerate(self.beam):
#             if beam is not None:
#                 if index == 1:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)]
#                 else:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)][1]
                    
#                 E_mean_fin = beam.E_mean_fin
#                 inertie = beam.inertie[0]
#                 aire = self.beam[index].aire
#                 distance_ai = self.distance_ai[index]
                
#                 @handcalc(override="long", precision=1, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():                    
#                     EI_eff_i = E_mean_fin * inertie + gamma_i * E_mean_fin * aire * distance_ai**2
#                     return EI_eff_i
#                 inter = val()
#                 EIeff_latex += f'EI eff {index+1}: '
#                 EIeff_latex += inter[0]
#                 EIeff_value = EIeff_value + inter[1]
                
#         @handcalc(override="params", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def sum_EI(): 
#             EI_eff_global = EIeff_value # Somme des EI efficace
#             return EI_eff_global
#         sumEI = sum_EI()
#         return (EIeff_latex + sumEI[0], sumEI[1])


#     def tau_2_max(self, Vz:float):
#         """Calcul de la contrainte de cisaillement maximale dans l'élément 2 selon Annexe B.4 de l'EN 1995

#         Args:
#             vz (float): Effort de cisaillement Vz en kN
#         """
        
#         beam_2 = deepcopy(self.beam[1])
#         K_cr = Cisaillement._from_parent_class(beam_2).K_cr
        
#         V_z = Vz * si.kN
#         h = self.beam[1].h_calcul / 2 + self.distance_ai[1]
#         gamma_3 = self.gamma_i["gamma 3"][1]
#         E_mean_fin_3 = self.beam[2].E_mean_fin
#         E_mean_fin_2 = self.beam[1].E_mean_fin
#         aire_3 = self.beam[2].aire
#         a_3 = self.distance_ai[2]
#         b_2 = self.beam[1].b_calcul
#         EI_eff = self.EI_eff[1]
        
#         @handcalc(override="long", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             tau_2 = V_z * (gamma_3 * E_mean_fin_3 * aire_3 * a_3 + 0.5 * E_mean_fin_2 * b_2 * h**2) / (b_2 * K_cr * EI_eff)
#             return tau_2
#         return val()
    

#     def sigma_i(self, My:float, beam: int=[1, 2, 3]):
#         """Contrainte de compression parallèle au fil selon Annexe B.3 de l'EN 1995

#         Args:
#             My (float): Moment selon l'axe y en kN.m
#             beam (int, optional): Beam 1 = Barre 1
#                                   Beam 2 = Barre 2 
#                                   Beam 3 = Barre 3.
#         """
#         Mf_z = My * si.kN*si.m
#         index = beam
#         if self.beam[index-1] is not None:
#             if index == 2:
#                 gamma_i = self.gamma_i["gamma "+str(index)]
#             else:
#                 gamma_i = self.gamma_i["gamma "+str(index)][1]
#             E_mean_fin = self.beam[index-1].E_mean_fin
#             distance_ai = self.distance_ai[index-1]
#             EI_eff = self.EI_eff[1]
            
#             @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 sigma_i = Mf_z * gamma_i * E_mean_fin * distance_ai / EI_eff
                
#                 return sigma_i
#             return val()
            
    
#     def sigma_mi(self, My:float, beam: int=[1, 2, 3]):
#         """Contrainte de flexion selon l'axe y selon Annexe B.3 de l'EN 1995

#         Args:
#             My (float): Moment selon l'axe y en kN.m
#             beam (int, optional): Beam 1 = Barre 1
#                                   Beam 2 = Barre 2 
#                                   Beam 3 = Barre 3.
#         """
#         M_y = My * si.kN*si.m
#         index = beam
#         if self.beam[index-1] is not None:
#             E_mean_fin = self.beam[index-1].E_mean_fin
#             h_i = self.beam[index-1].h_calcul
#             EI_eff = self.EI_eff[1]
            
#             @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 sigma_m_i = 0.5 * M_y * E_mean_fin * h_i / EI_eff
                
#                 return sigma_m_i
#             return val()
            
            
#     def F_i(self, Vz:float, connecteur: int=[1,2]):
#         """Permet de calculer la charge à reprendre par organe d'assemblage selon Annexe B.5 de l'EN 1995

#         Args:
#             Vz (float): Effort de cisaillement Vz en kN
#             connecteur (int, optional): Position des connecteurs:
#                         1 : connecteur entre les planches 1 et 2
#                         2 : connecteur entre les planches 2 et 3. 
#                         Defaults to [1,2].
#         """
#         V_z = Vz * si.kN
#         if connecteur == 1 :
#             connecteur = 0
#         gamma_i = self.gamma_i["gamma "+str(connecteur + 1)][1]
#         E_mean_fin = self.beam[connecteur].E_mean_fin
#         aire = self.beam[connecteur].aire
#         distance_ai = self.distance_ai[connecteur]
#         entraxe = self.entraxe[connecteur]
#         EI_eff = self.EI_eff[1]
        
#         @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             F_i = gamma_i * E_mean_fin * aire * distance_ai * entraxe * V_z / EI_eff
#             return F_i
#         return val()
    

# class Poteau_assemble_meca(Poutre_assemblee_meca):
#     def __init__(self, lo_y: si.mm, lo_z: si.mm, type_appuis: str=Compression.COEF_LF, **kwargs):
#         super(Poutre_assemblee_meca, self).__init__(**kwargs)
#         self.lo_comp = {"y":lo_y * si.mm, "z":lo_z * si.mm}
#         self.pole = []
#         if not isinstance(self.l, si.Physical):
#             self.l = self.l * si.mm
#         for i, beam in enumerate(self.beam):
#             if beam is not None:
#                 compression = Compression._from_parent_class(beam, lo_y=lo_y, lo_z=lo_z, type_appuis=type_appuis)
#                 self.pole.append(compression)
#             else:
#                 self.pole.append(None)
    
#     @property
#     def aire(self):
#         """Détermine la surface total du poteau assemblé mécaniquement
#         """
#         aire = 0
#         for beam in self.beam:
#             if beam is not None:
#                 aire = aire + beam.aire
#         return aire
    
#     @property
#     def Ief(self):
#         """Détermine l'inertie efficace à partir d'un module d'élasticité moyen Emean et 
#         de la rigidité équivalente en flexion EIef d'une poutre assemblée mécaniquement.
#         """
#         E_mean_beams = []
#         for pole in self.pole:
#             if pole is not None:
#                 E_mean_beams.append(pole.E_mean_fin)
#         E_mean = np.mean(E_mean_beams)
#         EI_ef = self.EI_eff[1]

#         @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def inertie_ef():
#             I_ef = EI_ef / E_mean # equC.4
#             return I_ef
#         return inertie_ef()
    
#     @property
#     def lamb_ef(self):
#         """ Retourne l'élancement d'un poteau assemblé mécaniquement en compression avec risque de flambement suivant l'axe de rotation z donc une direction de flèche suivant y"""
#         lo_z = self.lo_comp['z'].value * 10**3
#         I_z_ef = self.Ief[1]
#         A_tot = self.aire

#         @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             lamb_z_ef = lo_z * sqrt(A_tot / I_z_ef) # equC.3
#             return lamb_z_ef
#         return val()
    

#     @property
#     def lamb(self):
#         """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
#         dict_lamb = {"lamb,ef,z": self.lamb_ef}
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 lo_y = self.lo_comp['y'].value * 10**3
#                 coef_lef = pole.coef_lef
#                 I_y = pole.inertie[0].value * 10**12
#                 A = pole.aire.value * 10**6

#                 @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     lamb_y = (lo_y * coef_lef) / sqrt(I_y / A)
#                     return lamb_y
#                 dict_lamb["Pole"+str(i+1) + " lamb,y"] = val()
#         return dict_lamb
    

#     @property
#     def lamb_rel_Axe(self):
#         """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
#         lamb_ef_z = self.lamb['lamb,ef,z'][1]
#         E_0_05_beams, f_c0k_beams = [], []
#         for pole in self.pole:
#             if pole is not None:
#                 f_c0k_beams.append(float(pole.caract_meca.loc['fc0k'])* si.MPa)
#                 E_0_05_beams.append(float(pole.caract_meca.loc['E005'])* si.MPa)
#         f_c0k_mean = np.mean(f_c0k_beams)
#         E_0_05_mean = np.mean(E_0_05_beams)

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def lamb_rel_z():
#             lamb_rel_ef_z = (lamb_ef_z / pi) * sqrt(f_c0k_mean / E_0_05_mean)
#             return lamb_rel_ef_z
#         dict_lamb_rel = {"lamb,rel,ef,z": lamb_rel_z()}

#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 f_c0k = float(pole.caract_meca.loc['fc0k']) * si.MPa
#                 E_0_05 = int(pole.caract_meca.loc['E005']) * si.MPa
#                 lamb_y = self.lamb["Pole"+str(i+1) + " lamb,y"][1]

#                 @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     lamb_rel_y = (lamb_y / pi) * sqrt(f_c0k / E_0_05)
#                     return lamb_rel_y
#                 dict_lamb_rel["Pole"+str(i+1) + " lamb,rel,y"] = val()
#         return dict_lamb_rel
    

#     @property
#     def k_Axe(self):
#         """ Retourne le facteur Ky ou Kz (fonction de l'axe de flambement) """
#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def ky(beta_C:float, lamb_rel_y:float):
#             k_y = 0.5 * (1 + beta_C * (lamb_rel_y - 0.3) + lamb_rel_y ** 2)
#             return k_y
        
#         dict_k_axe = {}
#         beta_C_poles = 0
#         lamb_rel_ef_z = self.lamb_rel_Axe["lamb,rel,ef,z"][1]
        
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 beta_C = pole.beta_C
#                 if beta_C > beta_C_poles:
#                     beta_C_poles = beta_C
#                 lamb_rel_y = self.lamb_rel_Axe["Pole"+str(i+1) + " lamb,rel,y"][1]
#                 dict_k_axe["Pole"+str(i+1) + " ky"] =  ky(beta_C, lamb_rel_y)

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def kz():
#             k_ef_y = 0.5 * (1 + beta_C_poles * (lamb_rel_ef_z - 0.3) + lamb_rel_ef_z ** 2)
#             return k_ef_y
#         dict_k_axe["k,ef,z"] =  kz()
#         return dict_k_axe
    

#     @property
#     def kc_Axe(self):
#         """ Retourne le coefficient multiplicateur KcAxe  (axe = y ou z suivant axe de rotation en flambement) de fc,0,d """
#         lamb_rel_ef_z = self.lamb_rel_Axe["lamb,rel,ef,z"][1]
#         k_ef_z = self.k_Axe["k,ef,z"][1]

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def kcz():
#             k_c_ef_z = 1 / (k_ef_z + sqrt(k_ef_z** 2 - lamb_rel_ef_z ** 2))
#             return k_c_ef_z

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def kcy(k_y, lamb_rel_y):
#             k_c_y = 1 / (k_y + sqrt(k_y** 2 - lamb_rel_y ** 2))
#             return k_c_y
        
#         dict_kc_axe = {"kc,ef,z": kcz()}
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 k_y = self.k_Axe["Pole"+str(i+1) + " ky"][1]
#                 lamb_rel_y = self.lamb_rel_Axe["Pole"+str(i+1) + " lamb,rel,y"][1]
#                 dict_kc_axe["Pole"+str(i+1) + " kc,y"] =  kcy(k_y, lamb_rel_y)
#         return dict_kc_axe
    
    
#     def sigma_c_0_d(self, Fc0d: float):
#         """ Retourne la contrainte de compression axial en MPa avec:
#             Fc0d : la charge en kN de compression """
#         self.Fc_0_d = Fc0d * si.kN
#         Fc_0_d = self.Fc_0_d
#         A_tot = self.aire

#         @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             sigma_c_0_d = Fc_0_d / A_tot # equC.2
#             return sigma_c_0_d
#         value = val()
#         self.sigma_c_0_rd = value[1]
#         return value
    

#     def Vd_organe(self, Fc0d: float):
#         """ Retourne l'effort de cisaillement à prendr en compte sur un organe d'assemblage suivant l'annexe C §C.2.2 avec:
#             Fc0d : la charge en kN de compression """
#         self.Fc_0_d = Fc0d * si.kN
#         Fc_0_d = self.Fc_0_d
#         K_c_eff_z = self.kc_Axe["kc,ef,z"][1]
#         lamb_ef_z = self.lamb_ef[1]

#         if lamb_ef_z < 30:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def vd():
#                 V_d = Fc_0_d / (120 * K_c_eff_z)
#                 return V_d
#         elif 30 <= lamb_ef_z < 60:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def vd():
#                 V_d = Fc_0_d * lamb_ef_z / (3600 * K_c_eff_z)
#                 return V_d
#         else:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def vd():
#                 V_d = Fc_0_d / (60 * K_c_eff_z)
#                 return V_d
#         return vd()


#     def f_c_0_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
#         """Retourne le dictionnaire des résistances f,c,0,d de l'élément assemblé mécaniquement en MPa

#         Args:
#             loadtype (str): chargement de plus courte durée sur l'élément.
#             typecombi (str): type de combinaison, fondamentale ou accidentelle.

#         Returns:
#             float: f,c,0,d en MPa
#         """
#         self.dict_fc0d = {}
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 self.dict_fc0d["Pole"+str(i+1) + " fc,0,d"] = pole._f_type_d("fc0k", loadtype, typecombi)
#         self.dict_fc0d["fc,0,ef,d"] = np.mean([value[1] for value in self.dict_fc0d.values()])
#         return self.dict_fc0d
    

#     def taux_c_0_d(self, flexion: object=None):
#         """Retourne les taux de travaux de la compression axial.
#         Si l'élement est un poteau (donc avec un travail principalement en compression) et de la flexion combinée (EN 1995-1-1 §6.3.2), 
#         il est possible d'ajouter l'objet flexion et de vérifier cette combinaison.

#         Args:
#             flexion (object, optional): L'objet Flexion avec ces taux de travaux préalablement calculés. Default to None.

#         Returns:
#             list: retourne la liste des taux de travaux en %
#         """
#         self.taux_c_0_rd = {}
#         sigma_c_0_d = self.sigma_c_0_rd

#         f_c_0_ef_d = self.dict_fc0d["fc,0,ef,d"]
#         f_c0_d_tot = sum([value[1] for key, value in self.dict_fc0d.items() if key != "fc,0,ef,d"])
#         f_c0d_tot_Kcy = 0
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 kcy = self.kc_Axe["Pole"+str(i+1) + " kc,y"][1]
#                 fcod_y = self.dict_fc0d["Pole"+str(i+1) + " fc,0,d"][1]
#                 f_c0d_tot_Kcy = f_c0d_tot_Kcy + kcy * fcod_y
        
#         lamb_rel_y  = max([value[1] for key, value in self.lamb_rel_Axe.items() if key != "lamb,rel,ef,z"])
#         lamb_rel_z  = self.lamb_rel_Axe["lamb,rel,ef,z"][1]

#         K_c_eff_z = self.kc_Axe["kc,ef,z"][1]

#         if flexion and isinstance(flexion, Flexion):
#             taux_6_11 = flexion.taux_m_rd["equ6.11"]
#             taux_6_12 = flexion.taux_m_rd["equ6.12"]
#         else:
#             taux_6_11 = 0
#             taux_6_12 = 0
        
#         if lamb_rel_y < 0.3 and lamb_rel_z < 0.3:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 taux_6_2 = sigma_c_0_d / f_c0_d_tot # equ6.2
#                 taux_6_19 = (sigma_c_0_d / (f_c0_d_tot))**2 + taux_6_11 # equ6.19
#                 taux_6_20 = sigma_c_0_d / (f_c0_d_tot)**2 + taux_6_12 # equ6.20
#                 return taux_6_2, taux_6_19, taux_6_20
#             value = val()
#             self.taux_c_0_rd['equ6.19'] = value[1][3]
#             self.taux_c_0_rd['equ6.20'] = value[1][4]
#         else:      
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 taux_6_2 = sigma_c_0_d / f_c0_d_tot # equ6.2
#                 taux_6_23 = (sigma_c_0_d / (f_c0d_tot_Kcy)) + taux_6_11 # equ6.23
#                 taux_6_24 = sigma_c_0_d / (f_c_0_ef_d * K_c_eff_z) + taux_6_12 # equC.1
#                 return taux_6_2, taux_6_23, taux_6_24
#             value = val()
#             self.taux_c_0_rd['equ6.23'] = value[1][1]
#             self.taux_c_0_rd['equ6.24'] = value[1][2]

#         self.taux_c_0_rd['equ6.2'] = value[1][0]
#         return value



# if __name__=='__main__':
#     beam2 = Barre(60,200,"Rectangulaire", classe="C24", cs=1)
#     beam3 = Barre(60,100,"Rectangulaire", classe="C24", cs=1)
#     beam_ass = Poutre_assemblee_meca(beam_2=beam2, l=5000, disposition="Latérale", recouvrement=[0,120], Kser=[None,None,700], entraxe=[None,None,250], psy_2=0, beam_3=beam3)
#     pole_ass = Poteau_assemble_meca._from_parent_class(beam_ass, lo_y=5000, lo_z=5000, type_appuis="Rotule - Rotule")
#     # print(pole_ass.__dict__)
#     # print(pole_ass.lamb)
#     print(pole_ass.kc_Axe)