#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT



# Code permettant le calcul des principaux coef. EC5 (EN 1995) ainsi que les calculs simple de vérification au
# élément droit.
import os
import sys

import math as mt
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Project


# ================================ GLOBAL ==================================


class Beam(Project):
    LIST_SECTION = ["Rectangulaire","Circulaire"]
    LIST_TYPE_B = ["Massif","BLC", "LVL", "OSB 2", "OSB 3/4", "CP"]

    def __init__(self, b:int, h:int, section: str=LIST_SECTION[0], Hi: int=12, Hf: int=12, classe: str='C24', cs: int=1, **kwargs):
        """Classe qui définis les caractéristiques d'un élément droit. 
        Cette classe est hérité de la classe Project du module A0_Projet.py.

        Args:
            b (int): largeur de pose de la pièce en mm
            h (int): hauteur de pose de la pièce en mm
            section (str, optional): Type de section. Defaults to "Rectangulaire".
            Hi (int, optional): Humidité initiale de pose en %. Defaults to 12.
            Hf (int, optional): Humidité finale de pose en %. Defaults to 12.
            classe (str, optional): Classe mécanique du bois. Defaults to 'C24'.
            cs (int, optional): Classe de service de l'élément. Defaults to 1.
        """
        super().__init__(**kwargs)
        self.b = b
        self.h = h
        self.section = section
        self.Hi = Hi
        self.Hf = Hf
        self.classe = classe
        self.cs = cs

        for key, value in kwargs.items():
            setattr(self, key, value)
            
        self.sectionCalcul()
    

    def __data_from_csv(self, data_file:str):
        """ Retourne un dataframe d'un fichier CSV """
        repertory = os.getcwd() + "/data/" + data_file
        data_csv = pd.read_csv(repertory, sep=';', index_col=0)
        return data_csv
    
    
    def sectionCalcul(self, B90=0.25):
        """ Retourne la section de calcul en fonction de l'humidité de pose et celle d'utilisation avec pour argument:
                Hi : Humidité de pose en %
                Hf : Humidité finale en % selon AN Hf = 12%
                B90 : Coefficient de correction de section selon AN B90 = 0.25 %
                cote : Largeur ou hauteur de la section initiale en mm """
        self.b_calcul = self.b * (1 - B90 / 100 * (self.Hi - self.Hf))
        self.h_calcul = self.h * (1 - B90 / 100 * (self.Hi - self.Hf))
        
        
    @property
    def inertie(self):
        """ Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
            b ou d : Largeur ou diamètre de la poutre en mm
            h : Hauteur de la poutre en mm """
        if self.section == "Rectangulaire":
            iy = (self.b_calcul * self.h_calcul**3)/12
            iz = (self.h_calcul * self.b_calcul**3)/12
            return [iy, iz]

        elif self.Iy and self.Iz:
            return [self.Iy, self.Iz]

        else:
            i = (mt.pi * self.b_calcul ** 4) / 64
            return i
        
    
    @property
    def caract_meca(self):
        """ Retourne les caractéristiques méca du bois sous forme de dataframe pandas """
        data_csv_meca = self.__data_from_csv("caracteristique_meca_bois.csv")
        return data_csv_meca.loc[self.classe]
    
    
    @property
    def kmod(self):
        """ Retourne le kmod du bois avec pour argument:
            cs : classe de service (1, 2 ou 3) 
            duree_charge : temps de chargement ("instantanée" etc) """
        data_csv_kmod = self.__data_from_csv("kmod.csv")
        data_kmod = data_csv_kmod.loc[self.type_bois]
        return data_kmod.loc[data_kmod["CS"]==self.cs]


    @property
    def gamma_M(self):
        data_csv_gammaM = self.__data_from_csv("gammaM.csv")
        return data_csv_gammaM.loc[self.type_bois]
    
    
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
    

    @property
    def Kdef(self):
        data_csv_kdef = self.__data_from_csv("kdef.csv")
        kdef = float(data_csv_kdef.loc[self.type_bois][str(self.cs)])
        if self.Hi > 20:
            kdef += 1
        return kdef

    
    def f_type_d(self,typeCarac="fm0k", loadtype = "Permanente", typecombi="Fondamentales"):
        """Méthode donnant la résistance de calcul de l'élément fonction de la vérification

        Args:
            typeCarac (str, optional): Type de résistance caractéristique (flexion = "fm0k", compression = "fc0k" etc.). Defaults to "fm0k".
            loadtype (str, optional): Durée de chargement (Permanente, Court terme etc.). Defaults to "Permanente".
            typecombi (str, optional): Type de combinaison étudiée ("Fondamentales" ou " Accidentelles"). Defaults to "Fondamentales".

        Returns:
            float: Résistance de calcul en N/mm2 du type de vérification étudié.
        """
        self.Kmod = self.kmod[loadtype].iloc[0]
        self.Gm = self.gamma_M.loc[typecombi]
        self.f_type_rd = (float(self.caract_meca.loc[typeCarac]) * self.Kmod) / self.Gm
        return self.f_type_rd
    

    def fleche(self, long, Ed_WinstQ, Ed_Wnetfin, Ed_Wfin, type_ele='Element structuraux', type_bat='Batiments courants'):
        """ Retourne le taux de travail de la flèche en % avec pour argument:
            """
        data_csv_fleche = self.__data_from_csv("val_limite_fleche.csv")
        self.data_fleche= data_csv_fleche.loc[type_ele]
        self.data_fleche = self.data_fleche.loc[self.data_fleche["Type batiment"]==type_bat]

        self.Ed_WinstQ = Ed_WinstQ
        self.Ed_Wnetfin = Ed_Wnetfin 
        self.Ed_Wfin = Ed_Wfin
  
        try:
            self.Rd_WinstQ = long / int(self.data_fleche['Winst(Q)'].iloc[0])
        except ValueError:
            self.Ed_WinstQ = 0
            self.Rd_WinstQ = 1

        self.Rd_Wnetfin = long / int(self.data_fleche['Wnet,fin'].iloc[0])
        self.Rd_Wfin = long / int(self.data_fleche['Wfin'].iloc[0])

        self.taux_ELS = {}
        self.taux_ELS["Winst(Q)"] = self.Ed_WinstQ / self.Rd_WinstQ
        self.taux_ELS["Wnet,fin"] = self.Ed_Wnetfin / self.Rd_Wnetfin
        self.taux_ELS["Wfin"] = self.Ed_Wfin / self.Rd_Wfin

        return self.taux_ELS


    def module_young_mean_fin (self, psy2: int|float) -> float:
        """renvoie le E,mean,fin en fonction du Kdef et du psy2"""
        self.psy2 = psy2
        module = int(self.caract_meca.loc["E0mean"])
        self.E_mean_fin = module / (1 + self.psy2 * self.Kdef)
        return self.E_mean_fin


# ================================ FLEXION ==================================

class Flexion(Beam):
    """ Classe intégrant les formules de flexion à l'EC5 """
    COEF_LEF = {"Appuis simple" : [1, 0.9, 0.8], 
                "Porte à faux": [0.5, 0.8]}
    def __init__(self, lo, coeflef = COEF_LEF['Appuis simple'][1], pos = 1, *args, **kwargs):
        """ 
        b ou d : Largeur ou diamètre de la poutre en mm
        lo : Longeur de déversement en mm de la poutre
        h : Hauteur en mm
        Hi : Humidité de pose en %

        coeflef :
                    appuis simple :
                                    Moment constant : 1
                                    Charge répartie constante : 0.9
                                    Charge concentrée au milieu de la portée : 0.8
                    Porte à faux :
                                    Charge répartie constante : 0.5
                                    Charge concentrée agissant à l'extrémité libre : 0.8
        pos : Positionnement de la charge sur la hauteur de poutre
                                    
                                    si charge sur fibre comprimée pos = 0
                                    si charge sur fibre neutre pos = 1
                                    si charge sur fibre tendu pos = 2 
        
        type_bois : Défini le type de bois utilisé (Massif, BLC, LVL)
        section : Forme de la section (Rectangulaire, Circulaire etc.)
        classe : Classe mécanique du bois (C24 etc.) 
        """
        super().__init__(*args, **kwargs)
        self.lo = lo
        self.coeflef = coeflef
        self.pos = pos
  

    @property
    def Kh(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        kh = {}
        dim = {'y': self.h_calcul, 'z': self.b_calcul}

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

    @property
    def Km(self):
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
    def sig_m_crit(self):
        """ Retourne sigma m,crit pour la prise en compte du déversement d'une poutre """

        self.lef = self.lo * self.coeflef
        if self.pos == 0:
            self.lef = self.lef + 2 * self.h_calcul
        elif self.pos == 2:
            self.lef = self.lef - 0.5 * self.h_calcul
        return (0.78 * self.b_calcul ** 2 * int(self.caract_meca.loc['E005'])) / (self.h_calcul * self.lef)

    @property
    def lamb_rel_m(self):
        """ Retourne l'élancement relatif de la section avec pour argument """
        return mt.sqrt(float(self.caract_meca.loc['fm0k']) / self.sig_m_crit)

    @property
    def Kcrit(self):
        """ Retourne Kcrit le coef. de minoration de la résistance en flexion au déversement"""
        if self.lamb_rel_m <= 0.75:
            k_crit = 1
        elif 0.75 < self.lamb_rel_m <= 1.4:
            k_crit = 1.56 - 0.75 * self.lamb_rel_m
        else:
            k_crit = 1 / (self.lamb_rel_m ** 2)
        return k_crit
    
    
    def f_m_d(self, loadtype="Permanente", typecombi="Fondamentales"):
        return self.f_type_d("fm0k", loadtype, typecombi)
    
    
    def sigma_m_d(self, M, axe='y'):
        """ Retourne la contrainte sigma,m,d suivant sont axes de flexion avec :
            M : Momment max dans la barre en kN.m
            axe : Axe d'inertie quadratique à considérer"""
        self.Md = {'y': 0, 'z': 0}
        self.sigma_m_rd = {'y': 0, 'z': 0}
        
        if axe == 'y':
            inertie = self.inertie[0]
        elif axe == 'z':
            inertie = self.inertie[1]
        #prend en compte une section circulaire si aucun des deux axes
        else: 
            inertie = self.inertie
            
        self.Md[axe] = M 
        self.sigma_m_rd[axe] = (self.Md[axe] * 10**6) * (self.h_calcul/2) / inertie
        return self.sigma_m_rd
    

    def taux_m_d(self, taux_c_0_dy=0, taux_c_0_dz=0, taux_t_0_d=0):
        """ retourne les différents taux de travaux en flexion """
        self.taux_m_rd = {}
        self.taux_m_rd['equ6.11'] = (self.sigma_m_rd['y'] / (self.f_type_rd * self.Kh['y'])) + (self.Km * self.sigma_m_rd['z'] / (self.f_type_rd * self.Kh['z']))
        self.taux_m_rd['equ6.12'] = (self.Km * self.sigma_m_rd['y'] / (self.f_type_rd * self.Kh['y'])) + (self.sigma_m_rd['z'] / (self.f_type_rd * self.Kh['z']))
        self.taux_m_rd['equ6.33y'] = self.sigma_m_rd['y'] / (self.f_type_rd * self.Kh['y'] * self.Kcrit)
        self.taux_m_rd['equ6.33z'] = self.sigma_m_rd['z'] / (self.f_type_rd * self.Kh['z'] * self.Kcrit)
  
        if taux_c_0_dy or taux_c_0_dz: 
            self.taux_m_rd['equ6.35zyz'] = self.taux_m_rd['equ6.33y'] ** 2 + (self.sigma_m_rd['z'] / (self.f_type_rd * self.Kh['z'])) + taux_c_0_dz # 1er item axe de flexion pas au carré, 2eme item axe de flexion au carré, 3eme axe de compression
            self.taux_m_rd['equ6.35yzz'] = self.taux_m_rd['equ6.33y'] + (self.sigma_m_rd['z'] / (self.f_type_rd * self.Kh['z'])) ** 2 + taux_c_0_dz
            self.taux_m_rd['equ6.35yzy'] = self.taux_m_rd['equ6.33z'] ** 2 + (self.sigma_m_rd['y'] / (self.f_type_rd * self.Kh['y'])) + taux_c_0_dy
            self.taux_m_rd['equ6.35zyy'] = self.taux_m_rd['equ6.33z'] + (self.sigma_m_rd['y'] / (self.f_type_rd * self.Kh['y'])) ** 2 + taux_c_0_dy
        if taux_t_0_d:
            self.taux_m_rd['equ6.17'] = self.taux_m_rd['equ6.11'] + taux_t_0_d
            self.taux_m_rd['equ6.18'] = self.taux_m_rd['equ6.12'] + taux_t_0_d
        return self.taux_m_rd



# ================================ Traction ==================================

class Traction(Beam):
    """ Classe intégrant les formules de compression et d'instabilité au flambement à l'EC5 """

    def __init__(self, *args, **kwargs):
        """ 
        b ou d : Largeur ou diamètre de la poutre en mm
        h : Hauteur en mm
        section : Forme de la section (Rectangulaire, Circulaire etc.)
        classe : Classe mécanique du bois (C24 etc.)"""

        super().__init__(*args, **kwargs)


    @property
    def Kh(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        kh = {}
        dim = {'y': self.h_calcul, 'z': self.b_calcul}

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
        return kh
    
    
    def f_type_d(self, loadtype="Permanente", typecombi="Fondamentales"):
        return super().f_type_d("ft0k", loadtype, typecombi)
    
    
    def sigma_t_0_d(self, Ft0d):
        """ Retourne la contrainte de traxion axial en MPa avec:
            Ft0d : la charge en N de compression """
        self.Ft0d = Ft0d
        self.sigma_t_0_rd = self.Ft0d / (self.b_calcul * self.h_calcul)
        return self.sigma_t_0_rd


    def taux_t_0_d(self, taux_m_dy=0, taux_m_dz=0):
        """ Retourne le taux de travail en traction axial """
        kh = min(self.Kh['z'], self.Kh['y'])
        self.taux_t_0_rd = {}
        self.taux_t_0_rd['equ6.1'] = self.sigma_t_0_rd / (kh * self.f_type_rd)
        if taux_m_dy or taux_m_dz:
            self.taux_t_0_rd['equ6.17'] = self.taux_t_0_rd['equ6.1'] + taux_m_dy
            self.taux_t_0_rd['equ6.18'] = self.taux_t_0_rd['equ6.1'] + taux_m_dz
        
        return self.taux_t_0_rd
    

# ================================ Compression ==================================
 
class Compression(Beam):
    """ Classe intégrant les formules de compression et d'instabilité au flambement à l'EC5 """

    def __init__(self, lo = {'y':0, 'z':0}, coeflf=1, *args, **kwargs):
        """ 
        b ou d : Largeur ou diamètre de la poutre en mm
        h : Hauteur en mm
        lo : Longueur de flambement suivant l'axe de rotation (y ou z) en mm

        coeflf : Coefficient multiplicateur de la longueur pour obtenir la longeur efficace de flambement en
                    fonction des du type d'appuis :
                                                    Encastré 1 côté : 2
                                                    Rotule - Rotule : 1
                                                    Encastré - Rotule : 0.7
                                                    Encastré - Encastré : 0.5
                                                    Encastré - Rouleau : 1
        
        type_bois : Défini le type de bois utilisé (Massif, BLC, LVL)
        section : Forme de la section (Rectangulaire, Circulaire etc.)
        classe : Classe mécanique du bois (C24 etc.)"""

        super().__init__(*args, **kwargs)
        self.lo = lo
        self.coeflf = coeflf

    @property
    def lamb(self):
        """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb = {'y':0, 'z':0}
        lamb['y'] = (self.lo['y'] * self.coeflf) / mt.sqrt(self.inertie[1] / (self.b_calcul * self.h_calcul))
        lamb['z'] = (self.lo['z'] * self.coeflf) / mt.sqrt(self.inertie[0] / (self.b_calcul * self.h_calcul))
        return lamb
    
    @property
    def lamb_rel_Axe(self):
        """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb_rel_Axe = {'y':0, 'z':0}
        for cle, value in lamb_rel_Axe.items():
            lamb_rel_Axe[cle] = (self.lamb[cle] / mt.pi) * mt.sqrt(float(self.caract_meca.loc['fc0k']) / int(self.caract_meca.loc['E005']))
        return lamb_rel_Axe
    
    @property
    def betaC(self):
        if self.type_bois == 'Massif':
            betaC = 0.2
        else:
            betaC = 0.1
        return betaC
    
    @property
    def kAxe(self):
        """ Retourne le facteur Ky ou Kz (fonction de l'axe de flambement) """
        kA = {'y':0, 'z':0}
        for cle, value in kA.items():
            kA[cle] = 0.5 * (1 + self.betaC * (self.lamb_rel_Axe[cle] - 0.3) + self.lamb_rel_Axe[cle] ** 2)
        return kA

    @property
    def kc_Axe(self):
        """ Retourne le coefficient multiplicateur KcAxe  (axe = y ou z suivant axe de rotation en flambement) de fc,0,d """
        kc_Axe = {'y':0, 'z':0}
        for cle, value in kc_Axe.items():
            kc_Axe[cle] = 1 / (self.kAxe[cle] + mt.sqrt(self.kAxe[cle] ** 2 - self.lamb_rel_Axe[cle] ** 2))
        return kc_Axe


    def f_type_d(self, loadtype="Permanente", typecombi="Fondamentales"):
        return super().f_type_d("fc0k", loadtype, typecombi)
    
    
    def sigma_c_0_d(self, Fc0d):
        """ Retourne la contrainte de compression axial en MPa avec:
            Fc0d : la charge en N de compression """
        aire = self.b_calcul * self.h_calcul
        self.Fc0d = Fc0d
        self.sigma_c_0_rd = self.Fc0d / aire
        return self.sigma_c_0_rd
    
    
    def taux_c_0_d(self, taux_m_dy=0, taux_m_dz=0):
        """ Retourne le taux de travail de la compression axial """
        self.taux_c_0_rd = {}
        self.taux_c_0_rd['equ6.2'] = self.sigma_c_0_rd / self.f_type_rd

        if self.lamb_rel_Axe['y'] < 0.3 and self.lamb_rel_Axe['z'] < 0.3:
            self.taux_c_0_rd['equ6.19'] = (self.sigma_c_0_rd / (self.f_type_rd * self.kc_Axe['y']))**2 + taux_m_dy
            self.taux_c_0_rd['equ6.20'] = self.sigma_c_0_rd / (self.f_type_rd * self.kc_Axe['z'])**2 + taux_m_dz

        self.taux_c_0_rd['equ6.23'] = (self.sigma_c_0_rd / (self.f_type_rd * self.kc_Axe['y'])) + taux_m_dy
        self.taux_c_0_rd['equ6.24'] = self.sigma_c_0_rd / (self.f_type_rd * self.kc_Axe['z']) + taux_m_dz
        return self.taux_c_0_rd

    

    # ================================ COMPRESSION PERPENDICULAIRE ==================================

class Compression_perpendiculaire(Beam):
    """ Classe intégrant les formules de compression perpendiculaire à l'EC5 """

    def __init__(self, typea= 1, *args, **kwargs):
        """ 
        b ou d : Largeur ou diamètre de la poutre en mm
        h : Hauteur en mm
        Hi : Humidité de pose en %

        typea : Type d'appui
                                Appui continu : 0
                                Appui discret : 1
                                
        section : Forme de la section (Rectangulaire, Circulaire etc.)
        type_bois : Défini le type de bois utilisé (Massif, BLC, LVL)
        classe : Classe mécanique du bois (C24 etc.)"""

        super().__init__(*args, **kwargs)
        self.typea = typea
    

    def kc90(self, l0, l1):
        """ Retourne le facteur Kc,90 qui tient compte de la configuration de chargement, du fendage et de la déformation
            en compression avec pour argument :
            h : Hauteur de l'élement subissant la compression en mm
            lO : Longeur de l'appuis en compression en mm
            l1 : Distance entre deux appuis en mm (l et l)
            """
        if self.typea:
            if self.type_bois == 'Massif':
                if self.h_calcul <= 300:
                    if l1 >= 2 * self.h_calcul:
                        self.kc_90 = 1.5
                    else:
                        self.kc_90 = 1
                else:
                    self.kc_90 = 1.5
            elif self.type_bois == 'BLC':
                if self.h_calcul <= 300 and l0 <= 400:
                    if l1 >= 2 * self.h_calcul:
                        self.kc_90 = 1.75
                    else:
                        self.kc_90 = 1
                elif self.h_calcul <= 300 and l0 > 400:
                    self.kc_90 = 1
                else:
                    self.kc_90 = 1.75
            else:
                if self.h_calcul > 300:
                    self.kc_90 = 1.75
                else:
                    self.kc_90 = 1
        else:
            if self.type_bois == 'Massif':
                if self.h_calcul <= 300:
                    if l1 >= 2 * self.h_calcul:
                        self.kc_90 = 1.25
                    else:
                        self.kc_90 = 1
                else:
                    self.kc_90 = 1.5
            elif self.type_bois == 'BLC':
                if self.h_calcul <= 300:
                    if l1 >= 2 * self.h_calcul:
                        self.kc_90 = 1.5
                    else:
                        self.kc_90 = 1
                else:
                    self.kc_90 = 1.75
            else:
                if self.h_calcul > 300:
                    self.kc_90 = 1.75
                else:
                    self.kc_90 = 1

        return self.kc_90


    def f_type_d(self, loadtype="Permanente", typecombi="Fondamentales"):
        return super().f_type_d("fc90k", loadtype, typecombi)
    
    
    def sigma_c_90_d(self, Fc90d, l0, l1d=10000, l1g=10000, ad=0, ag=0):
        """ Retourne la contrainte normal de compression à 90 degrés en MPa avec pour argument:

            Fc90d : Charge en compression perpendiculaire en N
            l1d : Distance entre les charges en mm (l et l) (si pas de l1d ne rien mettre)
            l1g : Distance entre les charges en mm (l et l) (si pas de l1g ne rien mettre)
            ad : Distance depuis le bord jusqu'à l'appuis à droite (l) en mm (si pas de ad ne rien mettre)
            ag : Distance depuis le bord jusqu'à l'appuis à gauche (l) en mm (si pas de ad ne rien mettre)"""

        if l1d == 0 and l1g > 0:
            l1 = l1g
        elif l1g == 0 and l1d > 0:
            l1 = l1d
        else:
            l1 = min(l1d, l1g)
            
        self.kc90(l0, l1)

        self.aef = (l0 + min(30.0, ad, l0, 0.5*l1d) + min(30.0, ag, l0, 0.5*l1g)) * self.b_calcul
        # self.aef = (l0 + min(30.0, ad, 0.5 * l1d) + min(30.0, ag, 0.5 * l1g)) * self.b_calcul
        self.Fc90d = Fc90d
        self.sigma_c_90_rd = self.Fc90d / self.aef

        return self.sigma_c_90_rd


    def taux_c_90_d(self):
        """ Retourne le taux de travail de la compression perpendiculaire """
        self.taux_c_90_rd = {}
        self.taux_c_90_rd['equ6.3'] = (self.sigma_c_90_rd / (self.kc_90 * self.f_type_rd))

        return self.taux_c_90_rd


# ================================ CISAILLEMENT ==================================

class Cisaillement(Beam):
    """ Classe intégrant les formules de cisaillement à l'EC5 """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Kv = 1

    def kv(self, hef, x, i, ent="dessous"):
        """ Retourne le facteur d'entaille kv avec pour argument :
            h : Hauteur de la poutre en mm
            hef : Hauteur efficace de la poutre (hauteur - hauteur de l'entaille) en mm
            x : Distance entre le centre de réaction à l'appuis et le coin de l'entaille en mm
            i : pente de l'entaille (ex : i45° = 1) en % non multiplié par 100
            ent: Entaille sur le dessus ou dessous de la poutre
            type_bois : Type de bois
                                Massif
                                BLC
                                LVL """
        dictkn = {"Massif": 5,"BLC": 6.5, "LVL": 4.5}
        kn = dictkn[self.type_bois]
        a = hef / self.h_calcul
        if ent == "dessus" or ent == "Dessus":
            self.Kv = 1
        else:
            self.Kv = min(1, (kn * (1 + ((1.1 * i ** 1.5) / mt.sqrt(self.h_calcul)))) / (
                    mt.sqrt(self.h_calcul) * (mt.sqrt(a * (1 - a)) + 0.8 * (x / self.h_calcul) * mt.sqrt((1 / a) - a ** 2))))
        self.h_calcul = hef
        
        return self.Kv

    @property
    def Kcr(self):
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
            if self.h_calcul > 150 and self.type_bois == "Massif":
                k_cr = 0.67
            else:
                k_cr = 1
        elif self.cs == 2:
            if self.h_calcul > 150 and self.type_bois == "Massif":
                k_cr = 0.67
            elif self.type_bois == "BLC":
                k_cr = 0.67
            else:
                k_cr = 1
        else:
            k_cr = 0.67

        return k_cr

    
    def f_type_d(self, loadtype="Permanente", typecombi="Fondamentales"):
        return super().f_type_d("fvk", loadtype, typecombi)
       
    
    def tau_d(self, Vd):
        """ Retourne la contrainte tau en  MPa pour le cisaillement longitudinale d'une poutre rectangle
              Vd : Effort de cisaillement sur la poutre en N"""
        self.Vd = Vd
        bef = self.Kcr * self.b_calcul
        self.tau_rd = (1.5 * self.Vd) / (bef * self.h_calcul)
        return self.tau_rd


    def taux_tau_d(self):
        """ Retourne le taux de travail en cisaillement en % avec pour argument:
            Vd : Effort de cisaillement sur la poutre en N
            b : Largeur en mm
            kcr : coef de réduction de largeur avec prise en compte de la fissuration
            kv : coef prenant en compte la réduction de résistance liés à une entaille."""
        self.taux_tau_rd = {}
        self.taux_tau_rd['equ6.13'] = self.tau_rd / self.f_type_rd
        self.taux_tau_rd['equ6.60'] = self.tau_rd/ (self.Kv * self.f_type_rd)
        
        return self.taux_tau_rd


# ================================ Poutre assemblées mécaniquement Annexe B ==================================

class Poutre_assemblee_meca(Project):
    def __init__(self, beam_2:object, l: int|float, Kser: list=[0,None,0], entraxe: list=[1, None, 1], psy2: int|float=0, **kwargs):
        """Classe définissant une poutre composée d'au maximum 3 éléments connectés entre eux par liaisons mécanique 
        suivant la théorie de HEIMESHOFF Annexe B de l'EN 1995

        Args:
            beam_2 (object): objet contenant l'élément centrale de la poutre i=2, cette objet doit être issu de la classe élément ou de son héritage
            l (int | float): longueur de la poutre en mm

            Kser (list, optional): Rigidité de liaison par connecteur entre les éléments entre i=1/2 et i=2/3, en N/mm. 
                                    S'il n'y a que 2 éléments connectés, laisser le l'index correspondant vide (ex: [0, 2000]). Defaults to [0,None,0].

            entraxe (list, optional): Entraxe entre connecteur en mm suivant i=1 ou i=3. Defaults to [1, None, 1].

            type_action (int | float, optional): Psy 2 qui permet de prendre en compte le fluage dans le temps,
                                         si calcul en instantanée alors 0, 
                                         si intermédiaire = voir dans data/coeff_psy.csv 
                                         et enfin temps infini = 1. 
                                         Defaults to 0.
            **kwargs (object): beam_1 et ou beam_3
        """
        super().__init__(**kwargs)
        self.beam = [None , beam_2, None]
        self.l = l

        for key, value in kwargs.items():
            match key[0:4]:
                case "beam" if key[-1] == "2":
                    print("L'attribut ne peut pas être nommé beam_2, il est déjà pris par l'élément centrale ! beam_1 ou beam_3 possible !")
                case "beam" if key[-1] != "2" and 1<= int(key[-1]) <=3 :
                    self.beam[int(key[-1])-1] = value

        for index, beam in enumerate(self.beam):
            if beam is not None:
                beam.module_young_mean_fin(psy2)
                setattr(beam, "aire", beam.b_calcul * beam.h_calcul)
     
        self.Kser = Kser
        self.entraxe = entraxe

    
    @property
    def Kser_fin(self):
        """renvoie le Kser en fonction des Kdef des pièces assemblées et du psy2"""
        kser_fin = [None] *3
        for index, beam in enumerate(self.beam):
            if beam is not None and index != 1:
                kser_fin[index] = self.Kser[index] / (1 + self.beam[1].psy2 * (2 * (beam.Kdef * self.beam[1].Kdef)**0.5))
        return kser_fin


    @property
    def gamma_i(self):
        """renvoie le gamma"""
        gamma = [0, 1, 0]
        for index, beam in enumerate(self.beam):
            if beam is not None and index != 1:
                gamma[index] = (1 + mt.pi**(2) * beam.E_mean_fin * beam.aire * self.entraxe[index] / (self.Kser_fin[index] * self.l**(2)))**(-1)
        return gamma


    @property
    def distance_ai(self):
        """renvoie la distance à l'axe neutre de la pièce 2"""
        denominateur = 0
        for index, beam in enumerate(self.beam):
            if beam is not None:
                denominateur += (self.gamma_i[index] * beam.E_mean_fin * beam.aire)

        denominateur = 2 * denominateur 
        if self.beam[0] == None or self.beam[2] == None:
            for index, beam in enumerate(self.beam):
                if beam is not None and index != 1:
                    numerateur = self.gamma_i[index] * beam.E_mean_fin * beam.aire * (beam.h_calcul + self.beam[1].h_calcul)
        else:
            numerateur = (self.gamma_i[0] * self.beam[0].E_mean_fin * self.beam[0].aire * (self.beam[0].h_calcul + self.beam[1].h_calcul) 
                        - self.gamma_i[2] * self.beam[2].E_mean_fin  * self.beam[2].aire * (self.beam[1].h_calcul + self.beam[2].h_calcul))
        a2 =  numerateur / denominateur

        ai = [None, a2, None]
        for index, beam in enumerate(self.beam):
            if beam is not None and index != 1:
                ai[index] = (beam.h_calcul + self.beam[1].h_calcul)/2 - (-a2)
                if a2 >= 0 and index == 0 :
                    ai[index] -= a2
                elif a2 < 0 and index == 2:
                    ai[index] -= a2
                else:
                    ai[index] += a2
        return ai
        


    @property
    def EI_eff (self):
        """renvoie la rigidité de la poutre connectée"""
        ei_eff = 0
        for index, beam in enumerate(self.beam):
            if beam is not None:
                ei_eff += (beam.E_mean_fin * beam.inertie[0] + self.gamma_i[index] * beam.E_mean_fin * self.beam[index].aire * self.distance_ai[index]**2)
        return ei_eff


    # def verif_cisaillement(self, qz, l, EIeff):
    #     vz = qz * l / 2
    #     h = self.beam[1].h_calcul / 2 + self.distance_ai[0]
    #     kcr = 0.67
    #     tau = vz * (self.g3 * self.E3 * self.A3 * self.a3 + 0.5 * self.E2 * self.Bois2.b * h**2) / ( self.Bois2.b * kcr *EIeff)
    #     return tau  #fonctionne


    # def verif_flexion(self, qz, Ei, ai, EIeff, gammai, l, hi):
    #     my = qz * l**2 / 8
    #     sigma_i = my * gammai * Ei * ai / EIeff
    #     sigma_mi = 0.5 * Ei * hi *my / EIeff
    #     return sigma_i, sigma_mi #fonctionne


    # def verif_els(self, qz_qp, qz_c, EIeff_inst, EIeff_fin, l):
    #     u_inst_qp = 5 * qz_qp * l**4 / (384 * EIeff_inst)
    #     u_fin_qp = 5 * qz_qp * l**4 / (384 * EIeff_fin)
    #     u_creep = u_fin_qp - u_inst_qp
    #     u_inst_c = u_inst_qp * qz_c / qz_qp
    #     u_fin = u_inst_c + u_creep
    #     return u_fin



class Calcul_EC5_ELU(object):
    
    def __init__(self,  b: int, h: int, long: int, l0_flamb = {'y':0, 'z':0}, coeflf: float = 1, l0_dev: int = 0, coeflef=0.9, 
                   pos: int =0, section: str ="Rectangulaire", Hi: int = 12, Hf: int = 12, classe: str = 'C24', cs: int = 1,
                load_type: str = "Moyen terme", effort : list = [] , list_appuis: list = []):
        """Fait un  calcul standart à l'EC5 ELU STR pour vérifier les taux de travaux de la traction, compression, flexion, cisaillement, C90 et retourne ainsi les objets crée dans une liste

        Args:
            b (int): largeur de pose de la pièce en mm
            h (int): hauteur de pose de la pièce en mm
            long (int): longeur de la pièce en mm
            l0_flamb (dict, optional): longeur de flambement en mm suivant les deux axes y et z. Defaults to {'y':0, 'z':0}.
            coeflf (int, optional): Coefficient multiplicateur de la longueur pour obtenir la longeur efficace de flambement en
                    fonction des du type d'appuis :
                                                    Encastré 1 côté : 2
                                                    Rotule - Rotule : 1
                                                    Encastré - Rotule : 0.7
                                                    Encastré - Encastré : 0.5
                                                    Encastré - Rouleau : 1 . Defaults to 1.
                    
            l0_dev (int, optional): longeur de déversement en mm. Defaults to 0.
            coeflef (int, optional): Coefficient multiplicateur de la longeur pour obtenir la longueur efficace de déversement en fonction:
   
                       Appuis simple :
                                    Moment constant : 1
                                    Charge répartie constante : 0.9
                                    Charge concentrée au milieu de la portée : 0.8
                    Porte à faux :
                                    Charge répartie constante : 0.5
                                    Charge concentrée agissant à l'extrémité libre : 0.8 . Defaults to 0.9.. Defaults to 0.
         
            pos (int, optional): positionnement de la charge sur la hauteur de poutre
                        
                                    si charge sur fibre comprimée pos = 0
                                    si charge sur fibre neutre pos = 1
                                    si charge sur fibre tendu pos = 2  Defaults to 0.
         
            section (str, optional): Type de section. Defaults to "Rectangulaire".
            Hi (int, optional): Humidité initiale de pose en %. Defaults to 12.
            Hf (int, optional): Humidité finale de pose en %. Defaults to 12.
            classe (str, optional): Classe mécanique du bois. Defaults to 'C24'.
            cs (int, optional): Classe de service de l'élément. Defaults to 1.
            load_type (str, optional): Type de chargment de plus courte durée:
   
                       "Permanente",
                    "Moyen terme",
                    "Court terme",
                    "Instantanee", Defaults to "Moyen terme".
     
            effort (list, optional): Liste des efforts [{'t': 0, 'c': 0}, Vz, My, [Fc,90,d]] en N. Defaults to [0,0,0,[0,0]].
            list_appuis (list, optional): Liste des appuis [N°, type d'appuis, position mm, longueur d'appuis mm]. Defaults to [""].

        Returns:
            list: Liste des objets crées pour la compression, traction, flexion, cisaillment et C90
        """
        self.b = b
        self.h= h
        self.long = long
        self.section = section 
        self.classe = classe
        self.Hi = Hi
        self.Hf = Hf
        self.cs = cs
        self.load_type = load_type
        self.effort = effort
        self.list_appuis = list_appuis
  
        self.taux_c_0_dy, self.taux_c_0_dz, self.taux_t_0_d = 0, 0, 0
        self.Beam_C0, self.Beam_T0 = 0, 0
        
        if self.effort[0]['c'] > 0:
                self.calcul_compression(l0_flamb, coeflf)

        if self.effort[0]['t'] > 0:
            self.calcul_traction()

        if self.effort[1] > 0:
            self.calcul_cisaillement()

        if self.effort[2] > 0:
            self.calcul_flexion(l0_dev, coeflef, pos)

        if self.effort[3][0] > 0:
            self.calcul_compression90()
   
  
    def calcul_compression(self, l0_flamb, coeflf):
        self.Beam_C0 = Compression(self.b, self.h, l0_flamb, coeflf, self.section, self.Hi, self.Hf, self.classe, self.cs)
        self.Beam_C0.sigma_c_0_d(self.effort[0]['c'] * 10**3)
        self.Beam_C0.f_type_d("fc0k", self.load_type, "Fondamentales")
        self.Beam_C0.taux_c_0_d()
        self.taux_c_0_dy = self.Beam_C0.taux_c_0_rd['equ6.23']
        self.taux_c_0_dz = self.Beam_C0.taux_c_0_rd['equ6.24']
  
  
    def calcul_traction(self):
        self.Beam_T0 = Traction(self.b, self.h, self.section, self.Hi, self.Hf, self.classe, self.cs)
        self.Beam_T0.sigma_t_0_d(self.effort[0]['t'] * 10**3)
        self.Beam_T0.f_type_d("ft0k", self.load_type, "Fondamentales")
        self.Beam_T0.taux_t_0_d()
        self.taux_t_0_d = self.Beam_T0.taux_t_0_rd['equ6.1']
  
  
    def calcul_flexion(self, l0_dev, coeflef, pos):
        self.Beam_flexion = Flexion(self.b, self.h, l0_dev, coeflef, pos, self.section, self.Hi, self.Hf, self.classe, self.cs)
        self.Beam_flexion.sigma_m_d(self.effort[2])
        self.Beam_flexion.f_type_d("fm0k", self.load_type, "Fondamentales")
        self.Beam_flexion.taux_m_d(self.taux_c_0_dy, self.taux_c_0_dz, self.taux_t_0_d)


    def calcul_cisaillement(self):
        self.Beam_cisail = Cisaillement(self.b, self.h, self.section, self.Hi, self.Hf, self.classe, self.cs)
        self.Beam_cisail.tau_d(self.effort[1] * 10**3)
        self.Beam_cisail.f_type_d("fvk", self.load_type, "Fondamentales")
        self.Beam_cisail.taux_tau_d()
  
  
    def calcul_compression90(self):
        i, appuis_zero = 0, 0
        self.Beam_C90 = [0] * len(self.list_appuis)
        
        for appuis in self.list_appuis:
            Fc90d = self.effort[3][i] * 10**3
            
            posl1g = 0
            posl1d = 10**6
            ag = ad = 0

            for pos in self.list_appuis:
                if pos[2] == 0:
                    appuis_zero = 1
                if pos[2] < appuis[2] or pos[2] > appuis[2]:
                    if posl1g < pos[2] < appuis[2]:
                        posl1g = pos[2]
                    if posl1d > pos[2] > appuis[2]:
                        posl1d = pos[2]

            if posl1d != 10**6:
                l1d = posl1d - appuis[2]
            else:
                l1d = 10**6
            if posl1g or (appuis[2] != 0  and appuis_zero == 1):
                l1g = appuis[2] - posl1g
            else:
                l1g = 10**6

            if appuis[2] == 0 or appuis[2] == self.long:
                pass
            else:
                if l1g == 10**6 and appuis_zero == 0:
                    ag = appuis[2]-(appuis[3]/2)
                if l1d == 10**6:
                    ad = self.long-appuis[2]-(appuis[3]/2)

            self.Beam_C90[i] = Compression_perpendiculaire(self.b, self.h, 1, self.section, self.Hi, self.Hf, self.classe, self.cs)
            self.Beam_C90[i].sigma_c_90_d(Fc90d, appuis[3], l1d, l1g, ad, ag)
            self.Beam_C90[i].f_type_d("fc90k", self.load_type, "Fondamentales")
            self.Beam_C90[i].taux_c_90_d()
            # print('carac c90:', l1d, l1g, ad, ag, self.Beam_C90[i].aef)
            i += 1
  
    def get_object(self):
        return [self.Beam_flexion, self.Beam_cisail, self.Beam_C90, self.Beam_C0, self.Beam_T0]
            

if __name__=='__main__':
    beam2 = Beam(200, 100, Hi=20, classe="C18", cs=1)
    beam1 = Beam(200, 100, Hi=20, classe="C18", cs=1)
    beam3 = Beam(200, 100, Hi=20, classe="C18", cs=1)

    poutre_compo = Poutre_assemblee_meca(beam_2=beam2, beam_1=beam1, beam_3=beam3, l=5000, Kser=[30400, None, 30400], entraxe=[357, None, 357], psy2=0)
    print(poutre_compo.Kser_fin)
    print(poutre_compo.EI_eff)
    