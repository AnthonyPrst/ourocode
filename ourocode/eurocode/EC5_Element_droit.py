#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT



# Code permettant le calcul des principaux coef. EC5 (EN 1995) ainsi que les calculs simple de vérification au
# élément droit.
from copy import deepcopy
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
    LIST_SECTION = ["Rectangulaire","Circulaire"]
    LIST_TYPE_B = ["Massif", "BLC", "LVL", "OSB 2", "OSB 3/4", "CP"]
    CLASSE_WOOD = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").index)[2:]
    CS = ["1","2","3"]
    CARACTERISTIQUE = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").columns)
    LOAD_TIME = tuple(Projet._data_from_csv(Projet, "kmod.csv").columns)[1:]
    TYPE_ACTION = ("Fondamentales" ,"Accidentelles")
    TYPE_BAT = ("Bâtiments courants", "Bâtiments agricoles et similaires")
    TYPE_ELE = tuple(Projet._data_from_csv(Projet, "limite_fleche.csv").index.unique())
    B90 = 0.25

    def __init__(self, b:si.mm, h:si.mm, section: str=LIST_SECTION, Hi: int=12, Hf: int=12, classe: str=CLASSE_WOOD, cs: int=CS, effet_systeme: bool=("False", "True"), **kwargs):
        """Classe qui définit les caractéristiques d'un élément droit. 
        Cette classe est hérité de la classe Projet du module A0_Projet.py.

        Args:
            b (int): largeur de pose de la pièce en mm
            h (int): hauteur de pose de la pièce en mm
            section (str, optional): Type de section. Defaults to "Rectangulaire".
            Hi (int, optional): Humidité initiale de pose en %. Defaults to 12.
            Hf (int, optional): Humidité finale de pose en %. Defaults to 12.
            classe (str, optional): Classe mécanique du bois. Defaults to 'C24'.
            cs (int, optional): Classe de service de l'élément. Defaults to 1.
            effet_systeme: Détermine si l'effet système s'applique.
        """
        super().__init__(**kwargs)
        self.b = b * si.mm
        self.h = h * si.mm
        self.section = section
        self.Hi = Hi
        self.Hf = Hf
        self.classe = classe
        self.cs = cs
        self.effet_systeme = effet_systeme
        self._sectionCalcul()


    def _sectionCalcul(self):
        """ Retourne la section de calcul en fonction de l'humidité de pose et celle d'utilisation avec pour argument:
                Hi : Humidité de pose en %
                Hf : Humidité finale en % selon AN Hf = 12%
                B90 : Coefficient de correction de section selon AN B90 = 0.25 %
                cote : Largeur ou hauteur de la section initiale en mm """
        self.b_calcul = self.b * (1 - self.B90 / 100 * (self.Hi - self.Hf))
        self.h_calcul = self.h * (1 - self.B90 / 100 * (self.Hi - self.Hf))

    
    @property
    def aire(self):
        if self.section == __class__.LIST_SECTION[0]:
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
        data_csv_meca = self._data_from_csv("caracteristique_meca_bois.csv")
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
            type_b = __class__.LIST_TYPE_B[0]
        elif self.classe[0:2] == "GL":
            type_b = __class__.LIST_TYPE_B[1]
        elif self.classe[0:3] == "LVL":
            type_b = __class__.LIST_TYPE_B[2]
        elif self.classe[0:5] == "OSB/2":
            type_b = __class__.LIST_TYPE_B[3]
        elif self.classe[0:5] == "OSB/3" or self.classe[0:5] == "OSB/4":
            type_b = __class__.LIST_TYPE_B[4]
        else:
            type_b = __class__.LIST_TYPE_B[5]
        
        return type_b


    def __convert_latex_ftyped(self, latex: str, type_caract: str):
        end_index_rk = type_caract.find("k")
        type_caract = type_caract[1:end_index_rk]
        latex = latex.replace("f_{type_{Rd}", "f_{"+type_caract+"_{"+"d}")
        latex = latex.replace("f_{type_{k}", "f_{"+type_caract+"_{"+"k}")
        return latex
    
    
    def _f_type_d(self,typeCarac=CARACTERISTIQUE[0:6], loadtype=LOAD_TIME, typecombi=TYPE_ACTION):
        """Méthode donnant la résistance de calcul de l'élément fonction de la vérification

        Args:
            typeCarac (str, optional): Type de résistance caractéristique (flexion = "fm0k", compression = "fc0k" etc.). Defaults to "fm0k".
            loadtype (str, optional): Durée de chargement (Permanente, Court terme etc.). Defaults to "Permanente".
            typecombi (str, optional): Type de combinaison étudiée ("Fondamentales" ou " Accidentelles"). Defaults to "Fondamentales".

        Returns:
            float: Résistance de calcul en N/mm2 du type de vérification étudié.
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
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
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
                print("LVL non pris en compte dans cette fonction")
                kh[cle] = 1
        return kh       
    
    
    def Emean_fin(self, psy_2: float):
        """Renvoie le E,mean,fin en fonction du Kdef et du psy2"""
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
    

    def fleche(self, long:si.mm, Ed_WinstQ:si.mm=0, Ed_Wnetfin:si.mm=0, Ed_Wfin:si.mm=0, type_ele=TYPE_ELE, type_bat=TYPE_BAT):
        """Retourne le taux de travail de la flèche en % avec pour argument:

        Args:
            long (int): La longueur entre appuis à vérifier en mm
            Ed_WinstQ (float, optional): La flèche instanténée sous charge variable Q en mm. Defaults to 0.
            Ed_Wnetfin (float, optional): La flèche net finale en mm. Defaults to 0.
            Ed_Wfin (float, optional): La flèche finale en mm. Defaults to 0.
            type_ele (_type_, optional): Le type d'élément à vérifier. Defaults to TYPE_ELE.
            type_bat (_type_, optional): Le type de bâtiment sur lequel on vérifie notre élémennt. Defaults to TYPE_BAT.

        Returns:
            dict: Retourne le dictionnaire des taux de travail.
        """
        data_csv_fleche = self._data_from_csv("limite_fleche.csv")
        self.data_fleche= data_csv_fleche.loc[type_ele]
        self.data_fleche = self.data_fleche.loc[self.data_fleche["Type bâtiment"]==type_bat]
        self.taux_ELS = {}

        long = long * si.mm
        Ed_W_inst_Q = Ed_WinstQ * si.mm
        Ed_W_net_fin = Ed_Wnetfin * si.mm
        Ed_W_fin = Ed_Wfin * si.mm

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

                taux_W_inst_Q = Ed_W_inst_Q / Rd_W_inst_Q * 100 #%
                taux_W_net_fin = Ed_W_net_fin / Rd_W_net_fin * 100 #%
                taux_W_fin = Ed_W_fin / Rd_W_fin * 100 #%
                return taux_W_inst_Q, taux_W_net_fin, taux_W_fin
            
            value = val()
            self.taux_ELS["Winst(Q)"] = value[1][0]
            self.taux_ELS["Wnet,fin"] = value[1][1]
            self.taux_ELS["Wfin"] = value[1][2]

        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                Rd_W_net_fin = long / limit_W_net_fin
                Rd_W_fin = min(long / limit_W_fin, limit_U_fin_max)

                taux_W_net_fin = Ed_W_net_fin / Rd_W_net_fin * 100 #%
                taux_W_fin = Ed_W_fin / Rd_W_fin * 100 #%
                return taux_W_net_fin, taux_W_fin
            
            value = val()
            self.taux_ELS["Wnet,fin"] = value[1][0]
            self.taux_ELS["Wfin"] = value[1][1]
        return value 


# ================================ FLEXION ==================================

class Flexion(Barre):
    COEF_LEF = {"Appuis simple" : [1, 0.9, 0.8], 
                "Porte à faux": [0.5, 0.8]}
    LOAD_POS = (
        "Charge sur fibre comprimée",
        "Charge sur fibre neutre",
        "Charge sur fibre tendue"
        )

    def __init__(self, lo:si.mm, coeflef: float=COEF_LEF['Appuis simple'][1], pos: str=LOAD_POS, *args, **kwargs):
        """Classe permettant le calcul de la flexion d'une poutre bois selon l'EN 1995 §6.1.6, §6.2.3, §6.2.4 et §6.3.3.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.

        Args:
            lo (int): longueur de déversemment en mm
            coeflef (float): appuis simple :
                                            Moment constant : 1
                                            Charge répartie constante : 0.9
                                            Charge concentrée au milieu de la portée : 0.8
                            porte à faux :
                                            Charge répartie constante : 0.5
                                            Charge concentrée agissant à l'extrémité libre : 0.8.

            pos (str): positionnement de la charge sur la hauteur de poutre
        """
        super().__init__(*args, **kwargs)
        self.lo = lo * si.mm
        self.coeflef = coeflef
        self.pos = pos


    @property
    def K_h(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        return self._K_h()


    @property
    def K_m(self):
        """ Retourne le coef. Km qui reduit les contrainte d'une poutre scié en flexion """
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
        """ Retourne sigma m,crit pour la prise en compte du déversement d'une poutre """
        self.l_ef = self.lo * self.coeflef
        if self.pos == "Charge sur fibre comprimée":
            self.l_ef = self.l_ef + 2 * self.h_calcul
        elif self.pos == "Charge sur fibre tendue":
            self.l_ef = self.l_ef - 0.5 * self.h_calcul

        b_calcul = self.b_calcul
        h_calcul = self.h_calcul
        l_ef = self.l_ef
        E_0_05 = int(self.caract_meca.loc['E005']) * si.MPa
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_m_crit = (0.78 * b_calcul ** 2 * E_0_05) / (h_calcul * l_ef)
            return sigma_m_crit
        return val()
    

    @property
    def lamb_rel_m(self):
        """ Retourne l'élancement relatif de la section avec pour argument """
        f_m0k = float(self.caract_meca.loc['fm0k']) *si.MPa
        sigma_m_crit = self.sigma_m_crit[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            lamb_rel_m = sqrt(f_m0k / sigma_m_crit)
            return lamb_rel_m
        return val()
    

    @property
    def K_crit(self):
        """ Retourne K,crit le coef. de minoration de la résistance en flexion au déversement"""
        lamb_rel_m = self.lamb_rel_m[1]

        if lamb_rel_m <= 0.75:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_crit = 1
                return K_crit
        elif 0.75 < lamb_rel_m <= 1.4:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_crit = 1.56 - 0.75 * lamb_rel_m
                return K_crit
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_crit = 1 / (lamb_rel_m ** 2)
                return K_crit
        return val()
    
    
    def f_m_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Retourne la résistance f,m,d de l'élément en MPa

        Args:
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.

        Returns:
            float: f,m,d en MPa
        """
        return self._f_type_d("fm0k", loadtype, typecombi)
    
    
    def sigma_m_d(self, M: float, axe=['y', 'z']):
        """ Retourne la contrainte sigma,m,d suivant sont axes de flexion avec :
            M : Moment max dans la barre en kN.m
            axe : Axe d'inertie quadratique à considérer
        """
        if not hasattr(self, "Md"):
            self.Md = {'y': 0* si.kN*si.m, 'z': 0 * si.kN*si.m}
            self.sigma_m_rd = {'y': 0 * si.MPa, 'z': 0 * si.MPa}
        
        self.Md[axe] = M * si.kN*si.m

        if axe == 'y':
            I = self.inertie[0]
            h_calcul = self.h_calcul
        elif axe == 'z':
            I = self.inertie[1]
            h_calcul = self.b_calcul
        #prend en compte une section circulaire si aucun des deux axes
        else: 
            I = self.inertie
        
        M_d = self.Md[axe]
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_m_d = M_d * h_calcul / (I * 2)
            return sigma_m_d
        value = val()
        self.sigma_m_rd[axe] = value[1]
        return value
    

    def taux_m_d(self, compression: object=None, traction: object=None):
        """Retourne les différents taux de travaux en flexion.
        Si l'élement est une poutre (donc avec un travail principalement en flexion) et de la compression (EN 1995-1-1 §6.3.3) ou de la traction (EN 1995-1-1 §6.2.3) combinée, 
        il est possible d'ajouter l'objet Compression et Traction et de vérifier ces combinaisons.

        Args:
            compression (object, optional): L'objet Compression avec ces taux de travaux préalablement calculés. Defaults to None.
            traction (object, optional): L'objet Traction avec ces taux de travaux préalablement calculés. Defaults to None.

        Returns:
            list: retourne la liste des taux de travaux en %"""
        self.taux_m_rd = {}

        sigma_my_d = self.sigma_m_rd['y']
        sigma_mz_d = self.sigma_m_rd['z']
        f_m_d = self.f_type_rd
        K_h_y = self.K_h['y']
        K_h_z = self.K_h['z']
        K_m = self.K_m
        K_crit = self.K_crit[1]

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def base():
            taux_6_11 = sigma_my_d / (f_m_d * K_h_y) + K_m * sigma_mz_d / (f_m_d * K_h_z) # equ6.11
            taux_6_12 = K_m * sigma_my_d / (f_m_d * K_h_y) + sigma_mz_d / (f_m_d * K_h_z) # equ6.12
            taux_6_33y = sigma_my_d / (f_m_d * K_h_y * K_crit) # equ6.33
            taux_6_33z = sigma_mz_d / (f_m_d * K_h_z * K_crit) # equ6.33
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

        return (latex, self.taux_m_rd)



# ================================ Traction ==================================

class Traction(Barre):
    def __init__(self, *args, **kwargs):
        """Classe permettant le calcul de la Traction d'un élément bois selon l'EN 1995.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.
        """
        super().__init__(*args, **kwargs)


    @property
    def K_h(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        return self._K_h()
    
    
    def f_t_0_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Retourne la résistance f,t,0,d de l'élément en MPa

        Args:
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.

        Returns:
            float: f,t,0,d en MPa
        """
        return super()._f_type_d("ft0k", loadtype, typecombi)
    
    
    def sigma_t_0_d(self, Ft0d: si.kN):
        """ Retourne la contrainte de traxion axial en MPa avec:
            Ft0d : la charge en kN de compression """
        self.Ft_0_d = Ft0d * si.kN
        Ft_0_d = self.Ft_0_d
        A = self.aire

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_t_0_d = Ft_0_d / A
            return sigma_t_0_d
        value = val()
        self.sigma_t_0_rd = value[1]
        return value


    def taux_t_0_d(self):
        """Retourne le taux de travail en traction axial.

        Returns:
            float: taux de travail en %
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
        return value


# ================================ Compression ==================================
 
class Compression(Barre):
    COEF_LF = {"Encastré 1 côté" : 2,
                "Rotule - Rotule" : 1,
                "Encastré - Rotule" : 0.7,
                "Encastré - Encastré" : 0.5,
                "Encastré - Rouleau" : 1}
    def __init__(self, lo_y: si.mm, lo_z: si.mm, type_appuis: str=COEF_LF, *args, **kwargs):
        """ Classe permettant le calcul de la Compression d'un élément bois selon l'EN 1995.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.
        
        Args:
            lo : Longueur de flambement suivant l'axe de rotation (y ou z) en mm si pas de flambement alors 0
            type_appuis : Coefficient multiplicateur de la longueur pour obtenir la longeur efficace de flambement en
                        fonction des du type d'appuis :
                                                        Encastré 1 côté : 2
                                                        Rotule - Rotule : 1
                                                        Encastré - Rotule : 0.7
                                                        Encastré - Encastré : 0.5
                                                        Encastré - Rouleau : 1
        """

        super().__init__(*args, **kwargs)
        self.lo_comp = {"y":lo_y * si.mm, "z":lo_z * si.mm}
        self.lo_y = self.lo_comp['y']
        self.lo_z = self.lo_comp['z']
        self.type_appuis = type_appuis
        self.coef_lef = __class__.COEF_LF[type_appuis]

    @property
    def lamb(self):
        """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lo_y = self.lo_comp['y'].value * 10**3
        lo_z = self.lo_comp['z'].value * 10**3
        coef_lef = self.coef_lef
        I_y = self.inertie[0].value * 10**12
        I_z = self.inertie[1].value * 10**12
        A = self.aire.value * 10**6

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
        """Retourne la résistance f,c,0,d de l'élément en MPa

        Args:
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.

        Returns:
            float: f,c,0,d en MPa
        """
        return super()._f_type_d("fc0k", loadtype, typecombi)
    
    
    def sigma_c_0_d(self, Fc0d: si.kN):
        """ Retourne la contrainte de compression axial en MPa avec:
            Fc0d : la charge en kN de compression """
        self.Fc_0_d = Fc0d * si.kN
        Fc_0_d = self.Fc_0_d
        A = self.aire

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_c_0_d = Fc_0_d / A
            return sigma_c_0_d
        value = val()
        self.sigma_c_0_rd = value[1]
        return value
    
    
    def taux_c_0_d(self, flexion: object=None):
        """Retourne les taux de travaux de la compression axial.
        Si l'élement est un poteau (donc avec un travail principalement en compression) et de la flexion combinée (EN 1995-1-1 §6.3.2), 
        il est possible d'ajouter l'objet flexion et de vérifier cette combinaison.

        Args:
            flexion (object, optional): L'objet Flexion avec ces taux de travaux préalablement calculés. Default to None.

        Returns:
            list: retourne la liste des taux de travaux en %
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
        
        if lamb_rel_y < 0.3 and lamb_rel_z < 0.3:
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
        return value

    

    # ================================ COMPRESSION PERPENDICULAIRE ==================================

class Compression_perpendiculaire(Barre):
    TYPE_APPUIS = ("Appuis discret", "Appuis continu")
    def __init__(self, b_appuis:si.mm, l_appuis: si.mm, l1d: si.mm=10000, l1g: si.mm=10000, ad: si.mm=0, ag: si.mm=0, type_appuis_90: str=TYPE_APPUIS, *args, **kwargs):
        """Classe intégrant les formules de compression perpendiculaire selon l'EN 1995 §6.1.5.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.

        Args:
            b_appuis(int): largeur d'appuis en mm.
            l_appuis(int): longeur de l'appuis en mm.
            l1d(int) : Distance entre les charges en mm (l et l) (si pas de l1d ne rien mettre).
            l1g(int) : Distance entre les charges en mm (l et l) (si pas de l1g ne rien mettre).
            ad(int) : Distance depuis le bord jusqu'à l'appuis à droite (l) en mm (si pas de ad et au bord ne rien mettre).
            ag(int) : Distance depuis le bord jusqu'à l'appuis à gauche (l) en mm (si pas de ad et au bord ne rien mettre).
            type_appuis_90(str) : Type d'appuis (Appui continu, Appui discret)
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
        """ Retourne le facteur Kc,90 qui tient compte de la configuration de chargement, du fendage et de la déformation
            en compression avec pour argument :
            h : Hauteur de l'élement subissant la compression en mm
            lO : Longeur de l'appuis en compression en mm
            l1 : Distance la plus petite entre deux appuis en mm (l et l)
            """
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


    def f_c_90_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Retourne la résistance f,c,90,d de l'élément en MPa

        Args:
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.

        Returns:
            float: f,c,90,d en MPa
        """
        return super()._f_type_d("fc90k", loadtype, typecombi)
    
    
    def sigma_c_90_d(self, Fc90d: si.kN):
        """ Retourne la contrainte normal de compression à 90 degrés en MPa avec pour argument:

            Fc90d : Charge en compression perpendiculaire en kN
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
        """ Retourne le taux de travail de la compression perpendiculaire """
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
        return value
    
    
    def show_c90(self):
        """Affiche l'image des caractéristiques d'une compression perpendiculaire
        """
        self._show_element("C90_def.png")



class Compression_inclinees(Compression, Compression_perpendiculaire):
    def __init__(self, alpha: float=45, **kwargs):
        """Classe qui permet de calculer la compression inclinées par rapport au fil comme décrit à l'EN 1995 §6.2.2.
        Cette classe est hérité de la classe Compression et Compression_perpendiculaire, les deux provenants du module EC5_Element_droit.py.

        Args:
            alpha (float, optional): angle d'inclinaison en degrés de la compression. Defaults to 0.
        """
        # super(Compression, self).__init__(**kwargs)
        # super(Compression_perpendiculaire, self).__init__(**kwargs)
        super().__init__(**kwargs)
        self.alpha = alpha
    
    def sigma_c_alpha_d(self, Fcad: si.kN):
        """ Retourne la contrainte de compression inclinée en MPa avec:
            Fcad : la charge en kN de compression inclinée """
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
        """ Retourne le taux de travail de la compression inclinées par rapport au fil """
        self.taux_c_alpha_rd = {}
        f_c_0_d = self.f_c_0_d(loadtype, typecombi)[1]
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
        return value
    


# ================================ CISAILLEMENT ==================================

class Cisaillement(Barre):
    DICT_KN = {"Massif": 5,"BLC": 6.5, "LVL": 4.5}
    def __init__(self, **kwargs):
        """Classe qui permet de calculer le cisaillement d'une poutre comme décrit à l'EN 1995 §6.1.7 et §6.5.
        Cette classe est hérité de la classe Barre, provenant du module EC5_Element_droit.py.
        """
        super().__init__(**kwargs)
        self.K_v = 1
        self.h_ef = self.h_calcul

    @property
    def K_cr(self):
        """ Retourne le facteur de réduction de largeur Kcr avec pour argument:
            cs: Classe de service de la poutre
                                            CS 1 : 1
                                            CS 2 : 2
                                            CS 3 : 3
            h : Hauteur en mm
            type_bois : Type de bois
                                Massif : 0
                                BLC : 1
                                Autre : 2 """
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
        """Retourne le facteur d'entaille Kv pour une entaille au niveau d'un appuis

        Args:
            hef (int): Hauteur efficace de la poutre (hauteur - hauteur de l'entaille) en mm
            x (int):Distance entre le centre de réaction à l'appuis et le coin de l'entaille en mm
            i_lo (float): longueur horizontal de l'entaille en mm
            ent (tuple, optional): Entaille sur le dessus ou dessous de la poutre.

        Returns:
            float: facteur Kv
        """
        K_n = self.DICT_KN[self.type_bois]
        x = x * si.mm
        h_ef = hef * si.mm
        i = i_lo * si.mm / h_ef
        h_calcul = self.h_calcul
        print(h_ef, h_calcul)

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
        """Retourne la résistance f,v,d de l'élément en MPa

        Args:
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.

        Returns:
            float: f,v,d en MPa
        """
        return super()._f_type_d("fvk", loadtype, typecombi)
       
    
    def tau_d(self, Vd:si.kN):
        """ Retourne la contrainte tau en  MPa pour le cisaillement longitudinale d'une poutre rectangulaire
              Vd : Effort de cisaillement sur la poutre en kN"""
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
        """ Retourne le taux de travail en cisaillement en % """
        self.taux_tau_rd = {}
        tau_d = self.tau_rd
        f_v_d = self.f_type_rd
        K_v  = self.K_v

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_6_13 = tau_d / f_v_d
            taux_6_60 = tau_d/ (K_v * f_v_d)
            return taux_6_13, taux_6_60
         
        value = val()
        self.taux_tau_rd['equ6.13'] = value[1][0]
        self.taux_tau_rd['equ6.60'] = value[1][1]
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