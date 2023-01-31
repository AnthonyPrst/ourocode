# coding in UTF-8 

############# Le but de ce fichier est de regrouper toute les fonctions d'assemblage par organe métalique dans l'EN-1993 #############

import os
import sys

import math as mt
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from EC3_Element_droit import Element

#======================================================= BOULON =========================================================

class Boulon(Element):
    QUALITE_ACIER = tuple(Element._data_from_csv(Element, "qualite_acier.csv").index)

    def __init__(self, d:float, d0:float, qualite: float=QUALITE_ACIER, verif_filetage: bool=("False", "True"), filetage_EN1090: bool=("True", "False"), *args, **kwargs):
        """Configure un objet boulon permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Element du fichier EC3_Element_droit.py.

        Args:
            d (int): le diamètre de la tige en mm
            d0 (int): le diamètre de percage en mm
            qualite (float): classe d'acier de la tige (ex: 4.8)
            verif_filetage (bool, optional): défini si le filetage du boulon doit être vérifier, si c'est le cas alors True. Defaults to False.
            filetage_EN1090 (bool, optional): défini si le filetage est conforme à l'EN 1090, soit matricé. Si filetage usiné alors False. Defaults to True.
        """
        super().__init__(*args, **kwargs)
        self.d = d
        self.d0 = d0
        self.qualite = qualite
        self.verif_filetage = verif_filetage

        if filetage_EN1090:
            self.filetage_EN1090 = 1
        else:
            self.filetage_EN1090 = 0.85

        self.fyb = self.__qualite_acier.loc["fyb"]
        self.fub = self.__qualite_acier.loc["fub"]
        
        if type(self.__section_boulon) is pd.core.series.Series:
            self.As = self.__section_boulon.loc["As"]
            self.An = self.__section_boulon.loc["An"]
        else:
            self.As = 0
            self.An = mt.pi * (self.d/2)**2
        
            


    @property
    def __qualite_acier(self):
        df = self._data_from_csv("qualite_acier.csv")
        df = df.loc[self.qualite]
        return df
    

    @property
    def __section_boulon(self):
        try:
            df = self._Element__data_from_csv("section_boulon.csv")
            df = df.loc[self.d]
            return df
        except KeyError:
            print("Le diamètre ne peut pas être celui d'un boulon, vérifier As = 0, An = aire de d")
        


    def pince_metal_boulon(self, trous_oblongs: bool=("False", "True"), corrosion: bool=("False", "True")):
        """ Retourne les pinces du métal minimum dans un assemblage constitué de boulon (EN 1993-1-8 §3.5) avec:
            self.t : epaisseur de la pièce extérieur la plus mince de l'assemblage
            trous_oblongs : si les trous oblongs True sinon False
            corrosion : assemblage exposé à la corrosion = True sinon False 
            en_10025_5 : structure réalisées en acier conformes à l'EN 10025-5, True si vrai sinon False

            NE PRENDS PAS EN COMPTE P1,0 ou P1,i ou P2 diminué quand les boulons sont en quinconce (voir §3.5 -5)
        """
        
        en_10025_5 = False
        
        if self._Element__classe_acier.loc["norme"] == "EN 10025-5": 
            en_10025_5 = True
        
        
        pince = {}
        
        e1 = {"e1_min": round(1.2 * self.d0, 1)}
        e2 = {"e2_min": e1["e1_min"]}
        e3 = round(1.5 * self.d0, 1)
        e4 = e3
        p1 = {"p1_min": round(2.2 * self.d0, 1)}
        p2 = {"p2_min": round(2.4 * self.d0, 1)}
        

        if not en_10025_5:
            p1["p1_max"] = round(min(14 * self.t, 200))
            p2["p2_max"] = p1["p1_max"]

            if corrosion:
                e1["e1_max_corro"] = round(4 * self.t + 40)
                e2["e2_max_corro"] = e1["e1_max_corro"]

        else:
            e1["e1_max"] = round(max(8 * self.t, 125))
            e2["e2_max"] = e1["e1_max"]
            p1["p1_max"] = round(min(14 * self.t, 175))
            p2["p2_max"] = p1["p1_max"]
                
            
        if trous_oblongs :
            pince["e3"] = e3
            pince["e4"] = e4  
        else:
            pince["e1"] = e1
            pince["e2"] = e2
            
        pince["p1"] = p1      
        pince["p2"] = p2      

        return pince
    

    @property
    def FvRd(self) -> float:
        """Retourne la résistance en cisaillement de la partie fileté et lisse du boulon par plan en N
        """
        quali_1 = [4.6,5.6,8.8]
        if self.qualite in quali_1:
            fvrd_filete = 0.6 * self.fub * self.As / __class__.GAMMA_M["gamma_M2"]
        else:
            fvrd_filete = 0.5 * self.fub * self.As / __class__.GAMMA_M["gamma_M2"]
        
        fvrd_lisse = 0.6 * self.fub * self.An / __class__.GAMMA_M["gamma_M2"]
        return {"filetage": fvrd_filete * self.filetage_EN1090, "lisse": fvrd_lisse * self.filetage_EN1090}
    
    
    @property
    def FtRd(self) -> float:
        """Retourne la résistance en traction du boulon en N
        """
        return  0.9 * self.fub * self.As / __class__.GAMMA_M["gamma_M2"] * self.filetage_EN1090

    
    #Combinaison des efforts
    def taux_FvEd_FtEd(self, FvEd: int, FtEd: int) -> dict:
        """Retourne les taux de travail en cisaillement, en traction et combiné du boulon

        Args:
            FvEd (int): effort à reprendre en cisaillment en N
            FtEd (int): effort de traction à reprendre en N 

        Returns:
            dict: dictionnaire des taux de travail slon tab. 3.4 de l'EN 1993-1-8
        """
        self.FvEd = FvEd
        self.FtEd = FtEd
        self.taux_bl = {}

        self.taux_bl["taux_t_d"] = self.FtEd / self.FtRd

        self.taux_bl["taux_v_d_lisse"] = self.FvEd / self.FvRd["lisse"]

        if self.verif_filetage:
            self.taux_bl["taux_v_d_filetage"] = self.FvEd / self.FvRd["filetage"]
            if FvEd and FtEd:
                self.taux_bl["taux_v_t_d"] = self.FvEd / min(self.FvRd["lisse"], self.FvRd["filetage"]) + self.FtEd / (1.4 * self.FtRd)
        else:
            self.taux_bl["taux_v_t_d"] = self.FvEd / self.FvRd["lisse"] + self.FtEd / (1.4 * self.FtRd)
            
        return self.taux_bl
    
    
    def BpRd(self, d_ecrou: int, d_head_bl: int, *args) -> float:
        """Retourne la résistance au poinçonnement de la plaque en N

        Args:
            d_ecrou (int): diamètre extérieur de l'écrou en mm
            d_head_bl (int): diamètre de la tête de boulon en mm
            tp, i (int, args): épaisseur des plaques dans l'assemblage, sinon récupère le t de la classe Element

        Returns:
            float: résistance de calcul en N
        """
        dm = (d_ecrou + d_head_bl) / 2
        tp = min(self.t, *args)
        return (0.6 * mt.pi * dm * tp * self.fu) / __class__.GAMMA_M["gamma_M2"]

    
    def FbRd(self, e1: float ,e2: float , p1: float, p2: float) -> float:
        """Retourne la pression diamétrale en N. 
           ATTENTION: ne prends pas en compte les réductions de résistance lié au critère de l'assemblage (voir §3.6.1-10 et tab 3.4)
                        - jeux non normalisés -> 0.8 * Fb,Rd
                        - Trous oblongs avec axe perpendiculaire à la direction de la charge -> 0.6 * Fb,Rd
                        - Assemblage à une seule rangée de boulon en simple cisaillement -> rondelles + limitation de Fb,Rd <= 1.5*fu*d*t/gamma_M2

        Args:
            e1 (float): pince e1 en mm
            e2 (float): pince e2 en mm
            p1 (float): pince p1 en mm
            p2 (float): pince p2 en mm

        Returns:
            float: résistance de calcul en N 
        """
        self._alpha_b = min(e1 / (3 * self.d0), (p1 / (3 * self.d0)) - 0.25, self.fub / self.fu, 1)
        self._k1 = min(2.8 * e2 / self.d0 - 1.7, 1.4 * p2 / self.d0 - 1.7, 2.5)
        # print(alpha_b , k1)
        return self._k1 * self._alpha_b * self.fu * self.d * self.t / __class__.GAMMA_M["gamma_M2"]




class Soudure(Element):
    def __init__(self, gorge: int, l: int, retour_soudure: bool=("False", "True"), alpha: float=90, *args, **kwargs):
        """Configure un objet soudure permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Element du fichier EC3_Element_droit.py.  

        Args:
            gorge (int): dimension de la gorge "a" en mm
            l (int): longueur de soudure "brute" sans cratère en mm
            retour_soudure (bool, optional): _description_. Defaults to False.
            alpha (int | float, optional): angle en degré de la de la pièce 2 sur la pièce 1. Defaults to 90.
        """
     
        super().__init__(*args, **kwargs)
        self.gorge = gorge
        self.l = l
        self.retour_soudure = retour_soudure
        self.alpha = alpha

        self.verif_soudure()


    @property
    def beta_w(self):
        return float(self._Element__classe_acier.loc["betaW"])


    @property
    def lef(self):
        if not self.retour_soudure:
            return self.l - 2 * self.gorge
        else:
            return self.l


    def verif_soudure(self):
        if 60 <= self.alpha <= 120:
            if self.gorge >= 3:
                if self.l > max(30, 6*self.gorge):
                    return True
                else:
                    print(f"La longueur du cordon de soudure est trop petite, elle doit être supérieur à {min(30, 6*self.gorge)} mm")
            else:
                print("La gorge doit être au minimum de 3mm")
        else:
            print("L'angle entre les deux pièces à souder doit être compris entre 60° et 120°")
            return False
        
    
    def beta_Lw1(self, Lj: int) -> float:
        """Calcul le facteur beta_Lw,1 qui dimminue la résistance pour des cordons de soudure des assemblages par recouvrement (à plat)

        Args:
            Lj (int): Longueur de recouvrement des plats en mm
        """
        return min((1.2 - 0.2 * Lj) / (150 * self.gorge), 1)

    
    def cordon_frontal(self, Ned: float):
        """Calcul un cordon de soudure frontale et retourne le taux de travail.

        Args:
            Ned (int | float): Effort de traction en N.
        """
        return (self.beta_w * self.GAMMA_M["gamma_M2"] * (Ned * mt.sqrt(2)) / self.fu) / (self.gorge * self.lef)


    def cordon_laterale(self, Ned: float):
        """Calcul un cordon de soudure latérale et retourne le taux de travail.

        Args:
            Ned (int | float): Effort de cisaillement du cordon en N.
        """
        return (self.beta_w * self.GAMMA_M["gamma_M2"] * (Ned * mt.sqrt(3)) / self.fu) / (self.gorge * self.lef)


    def cordon_oblique(self, alpha_cordon: float, Ned: float):
        """Calcul un cordon de soudure oblique et retourne le taux de travail.

        Args:
            Ned (int | float): Effort de traction en N.
        """
        self.alpha_cordon = alpha_cordon
        return (self.beta_w * self.GAMMA_M["gamma_M2"] * (Ned * mt.sqrt(3 - mt.sin(mt.radians(self.alpha_cordon))**2)) / self.fu) / (self.gorge * self.lef)


    def cordon_pieces_obliques(self, Ned: float):
        """Calcul un cordon de soudure sur des pièces à positionnement obliques et retourne le taux de travail.

        Args:
            Ned (int | float): Effort de traction en N.
        """
        if self.alpha < 90:
            return (self.beta_w * self.GAMMA_M["gamma_M2"] * (Ned * mt.sqrt(2 - mt.sin(mt.radians(self.alpha)))) / self.fu) / (self.gorge * self.lef)
        elif self.alpha > 90:
            return (self.beta_w * self.GAMMA_M["gamma_M2"] * (Ned * mt.sqrt(2 + mt.sin(mt.radians(self.alpha)))) / self.fu) / (self.gorge * self.lef)


    def critere_generale(self, FvEd: float, FaxEd: float) -> float:
        """Calcul le critère générale de Von Mises d'une soudure et retourne le taux.

        Args:
            FvEd (int | float): Effort de cisaillement sur la en N
            FaxEd (int | float): Effort de traction sur la soudure en N

        Returns:
            (float) : Taux de travail de la soudure
        """
        self.tau_para = FvEd / (self.gorge * self.lef)
        self.sigma_perpend = (FaxEd * mt.cos(mt.radians(self.alpha/2)))/ (self.gorge * self.lef)
        self.tau_perpend = (FaxEd * mt.cos(mt.radians(self.alpha/2)))/ (self.gorge * self.lef)
        self.sigma_e = mt.sqrt(self.sigma_perpend**2 + 3 * (self.tau_perpend**2 + self.tau_para**2))

        self.sigma_Rd = self.fu / (self.beta_w * self.GAMMA_M["gamma_M2"])
        self.sigma_perpend_Rd = (0.9 * self.fu) / self.GAMMA_M["gamma_M2"]
        return max(self.sigma_e / self.sigma_Rd, self.sigma_perpend / self.sigma_perpend_Rd)
    
    

    def soudure_discontinue(self, b: int, b1: int, t1: int, corrosion: bool=("False", "True")):
        """_summary_

        Args:
            b (int): voir EC3 1-8
            b1 (int): hauteur en mm de la pièce 2 soudé sur la pièce 1
            t (int): défini dans la classe Element
            t1 (int): épaisseur en mm de la piece 2 soudé sur pièce 1
            corrosion (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if corrosion:
            print("Il n'est pas possible d'avoir une soudure discontinue en ambiance corrosive")
            return False

        lwe = max(0.75 * b, 0.75 * b1)
        l1 = min(16 * self.t, 16 * t1, 200)
        l2 = min(12 * self.t, 12 * t1, 0.25 * b, 200)
        return {"Lwe": lwe, "L1": l1, "L2": l2}




if __name__ == "__main__":
    soudure = Soudure(gorge=4, l=140, retour_soudure=True, alpha=90, classe_acier="S235")
    
    print(soudure.critere_generale(0, 100135))