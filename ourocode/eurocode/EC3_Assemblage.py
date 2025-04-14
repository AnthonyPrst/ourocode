# coding in UTF-8 

############# Le but de ce fichier est de regrouper toute les fonctions d'assemblage par organe métalique dans l'EN-1993 #############
from math import *
import pandas as pd

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.EC3_Element_droit import Element

#======================================================= Tige =========================================================

class Tige(Element):
    QUALITE_ACIER = tuple(str(key) for key in Element._data_from_csv(Element, "qualite_acier.csv").index)

    def __init__(self, d:si.mm, d0:si.mm, qualite: str=QUALITE_ACIER, verif_filetage: bool=("False", "True"), filetage_EN1090: bool=("True", "False"), *args, **kwargs):
        """Configure un objet Tige permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Element du fichier EC3_Element_droit.py.

        Args:
            d (int): le diamètre de la tige en mm
            d0 (int): le diamètre de percage en mm
            qualite (float): classe d'acier de la tige (ex: 4.8)
            verif_filetage (bool, optional): définit si le filetage du boulon doit être vérifier, si c'est le cas alors True. Defaults to False.
            filetage_EN1090 (bool, optional): définit si le filetage est conforme à l'EN 1090, soit matricé. Si filetage usiné alors False. Defaults to True.
        """
        super().__init__(*args, **kwargs)
        self.d = d * si.mm
        self.d0 = d0 * si.mm
        self.qualite = qualite
        self.verif_filetage = verif_filetage

        if filetage_EN1090:
            self.filetage_EN1090 = 1
        else:
            self.filetage_EN1090 = 0.85

        self.fyb = self.__qualite_acier.loc["fyb"] * si.MPa
        self.fub = self.__qualite_acier.loc["fub"] * si.MPa
        
        if type(self.__section_boulon) is pd.core.series.Series:
            self.As = self.__section_boulon.loc["As"] * si.mm**2
            self.An = self.__section_boulon.loc["An"] * si.mm**2
        else:
            self.As = 0 * si.mm**2
            self.An = pi * (self.d.value*10**3/2)**2 * si.mm**2
            

    @property
    def __qualite_acier(self):
        df = self._data_from_csv("qualite_acier.csv")
        df = df.loc[self.qualite]
        return df
    

    @property
    def __section_boulon(self):
        try:
            df = self._data_from_csv("section_boulon.csv")
            df = df.loc[self.d.value*10**3]
            return df
        except KeyError:
            print("Le diamètre ne peut pas être celui d'un boulon, vérifier As = 0, An = aire de d")
        


    def pince_metal_boulon(self, trous_oblongs: bool=("False", "True"), corrosion: bool=("False", "True")):
        """ Retourne les pinces du métal minimum dans un assemblage constitué de tige (EN 1993-1-8 §3.5) avec:
            trous_oblongs : si les trous oblongs True sinon False
            corrosion : assemblage exposé à la corrosion = True sinon False 
            en_10025_5 : structure réalisées en acier conformes à l'EN 10025-5, True si vrai sinon False

            NE PRENDS PAS EN COMPTE P1,0 ou P1,i ou P2 diminué quand les boulons sont en quinconce (voir §3.5 -5)
        """
        
        en_10025_5 = False
        
        if self._Element__classe_acier.loc["norme"] == "EN 10025-5": 
            en_10025_5 = True
        
        
        pince = {}
        
        e1 = {"e1_min": round(1.2 * self.d0.value*10**3, 1)}
        e2 = {"e2_min": e1["e1_min"]}
        e3 = round(1.5 * self.d0.value*10**3, 1)
        e4 = e3
        p1 = {"p1_min": round(2.2 * self.d0.value*10**3 , 1)}
        p2 = {"p2_min": round(2.4 * self.d0.value*10**3, 1)}
        

        if not en_10025_5:
            p1["p1_max"] = round(min(14 * self.t.value*10**3, 200))
            p2["p2_max"] = p1["p1_max"]

            if corrosion:
                e1["e1_max_corro"] = round(4 * self.t.value*10**3 + 40)
                e2["e2_max_corro"] = e1["e1_max_corro"]

        else:
            e1["e1_max"] = round(max(8 * self.t.value*10**3, 125))
            e2["e2_max"] = e1["e1_max"]
            p1["p1_max"] = round(min(14 * self.t.value*10**3, 175))
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
        """Retourne la résistance en cisaillement de la partie fileté et lisse d'une tige par plan en N
        """
        quali_1 = [4.6,5.6,8.8]
        if self.qualite in quali_1:
            fvrd_filete = 0.6 * self.fub * self.As / self.GAMMA_M["gamma_M2"]
        else:
            fvrd_filete = 0.5 * self.fub * self.As / self.GAMMA_M["gamma_M2"]
        
        fvrd_lisse = 0.6 * self.fub * self.An / self.GAMMA_M["gamma_M2"]
        return {"filetage": fvrd_filete * self.filetage_EN1090, "lisse": fvrd_lisse * self.filetage_EN1090}
    
    
    @property
    def FtRd(self) -> float:
        """Retourne la résistance en traction de la tige en N
        """
        f_ub = self.fub
        A_s = self.As
        gamma_M2 =  self.GAMMA_M["gamma_M2"]
        coef_filetage_EN1090 = self.filetage_EN1090
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            F_t_Rd = 0.9 * f_ub * A_s /  gamma_M2 * coef_filetage_EN1090
            return F_t_Rd
        return val()

    
    #Combinaison des efforts
    def taux_FvEd_FtEd(self, FvEd: float=0, FtEd: float=0):
        """Retourne les taux de travail en cisaillement, en traction et combiné d'une tige

        Args:
            FvEd (int): effort à reprendre en cisaillment en kN
            FtEd (int): effort de traction à reprendre en kN 

        Returns:
            dict: dictionnaire des taux de travail slon tab. 3.4 de l'EN 1993-1-8
        """
        self.FvEd = FvEd * si.kN
        self.FtEd = FtEd * si.kN
        self.taux_bl = {}

        self.taux_bl["taux_t_d"] = self.FtEd / self.FtRd[1]

        self.taux_bl["taux_v_d_lisse"] = self.FvEd / self.FvRd["lisse"]

        if self.verif_filetage:
            self.taux_bl["taux_v_d_filetage"] = self.FvEd / self.FvRd["filetage"]
            if FvEd and FtEd:
                self.taux_bl["taux_v_t_d"] = self.FvEd / min(self.FvRd["lisse"], self.FvRd["filetage"]) + self.FtEd / (1.4 * self.FtRd[1])
        else:
            self.taux_bl["taux_v_t_d"] = self.FvEd / self.FvRd["lisse"] + self.FtEd / (1.4 * self.FtRd[1])
            
        return self.taux_bl
    
    
    def Bp_Rd(self, d_ecrou: int, d_head_bl: int, *args):
        """Retourne la résistance au poinçonnement de la plaque en N

        Args:
            d_ecrou (int): diamètre extérieur de l'écrou en mm
            d_head_bl (int): diamètre de la tête de boulon en mm
            tp, i (int, args): épaisseur des plaques dans l'assemblage, sinon récupère le t de la classe Element

        Returns:
            float: résistance de calcul en N
        """
        d_ecrou = d_ecrou * si.mm
        d_head_bl = d_head_bl * si.mm
        dm = (d_ecrou + d_head_bl) / 2
        tp = min(self.t.value*10**3, *args) * si.mm
        return (0.6 * pi * dm * tp * self.fu) / __class__.GAMMA_M["gamma_M2"]

    
    def Fb_Rd(self, e1: float ,e2: float , p1: float, p2: float):
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
        e_1 = e1 * si.mm
        e_2 = e2 * si.mm
        p_1 = p1 * si.mm
        p_2 = p2 * si.mm
        f_u = self.fu
        f_u_b = self.fub
        d_0 = self.d0
        d = self.d
        t = self.t
        gamma_M_2 = __class__.GAMMA_M["gamma_M2"]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            alpha_b = min(e_1 / (3 * d_0), (p_1 / (3 * d_0)) - 0.25, f_u_b / f_u, 1)
            k_1 = min(2.8 * e_2 / d_0 - 1.7, 1.4 * p_2 / d_0 - 1.7, 2.5)

            F_b_Rd = k_1 * alpha_b * f_u_b * d * t / gamma_M_2
            return F_b_Rd
        return val()


    def Veff_Rd(self, Lnt: float, Lvt: float, effort: str=("Centré", "Excentré")):
        """Retourne la résistance en cisaillement bloc de l'assemblage selon §3.10.2.

        Args:
            Lnt (float, optional): longueur nette soumise à la traction déduction faite des trous de boulons en mm.
            Lvt (float, optional): longueur nette soumise au cisaillement déduction faite des trous de boulons en mm.
            effort (str, optional): définit si l'effort est centré ou excentré par rapport au cisaillement de bloc. Defaults to ("Centré", "Excentré").

        Returns:
            float: résistance de calcul en N 
        """
        A_nt = Lnt * si.mm * self.t
        A_vt = Lvt * si.mm * self.t
        f_y = self.fy
        f_u = self.fu
        k_ex = 1
        gamma_M_2 = __class__.GAMMA_M["gamma_M2"]

        if effort == "Excentré":
            k_ex = 0.5

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            V_eff_nt_Rd = (A_nt * f_u) / gamma_M_2
            V_eff_nv_Rd = 1/sqrt(3) * (A_vt * f_y) / gamma_M_2
            V_eff_Rd = k_ex * V_eff_nt_Rd + V_eff_nv_Rd
            return V_eff_Rd
        return val()
    
    
    def taux_Veff_d(self, N_Ed: float=0, N_Veff_Rd: float=0, V_Ed: float=0, V_Veff_Rd: float=0):
        """Retourne le taux de travail du cisaillement bloc avec un effort combiné ou non selon l'EN 1993 §3.10.2 et le CNC2M-N0175-REC §2.1(10).    

        Args:
            N_Ed (float, optional): effort normal en kN. Defaults to 0.
            V_Ed (float, optional): effort de cisaillement en kN. Defaults to 0.
            V_Veff_Rd (float, optional): résistance de l'effort de traction en rupture bloc en kN, si pas d'effort ne rien mettre. Defaults to 0.
            V_Veff_Rd (float, optional): résistance de l'effort de cisaillement en kN ,si pas d'effort ne rien mettre. Defaults to 0.

        """
        N_Ed = N_Ed * si.kN
        N_Veff_Rd = N_Veff_Rd * si.kN
        V_Ed = V_Ed * si.kN
        V_Veff_Rd = V_Veff_Rd * si.kN

        # selon CNC2M §2.1(10)
        if N_Ed.value and V_Ed.value:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                taux = N_Ed / N_Veff_Rd + V_Ed / V_Veff_Rd
                return taux
            return val()
        elif N_Ed.value:
            taux = N_Ed / N_Veff_Rd
        else:
            taux = V_Ed / V_Veff_Rd
        return taux



class Soudure(Element):
    def __init__(self, t2: si.mm, gorge: si.mm, l: si.mm, retour_soudure: bool=("False", "True"), alpha: float=90, *args, **kwargs):
        """Configure un objet soudure permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Element du fichier EC3_Element_droit.py.  

        Args:
            t2 (float): épaisseur de la pièce 2 sur laquelle on soude en mm
            gorge (int): dimension de la gorge "a" en mm
            l (int): longueur de soudure "brute" sans cratère en mm
            retour_soudure (bool, optional): détermine si un retour de la soudure est fait si oui alors True. Defaults to False.
            alpha (int | float, optional): angle en degré de la de la pièce 2 sur la pièce 1. Defaults to 90.
        """
     
        super().__init__(*args, **kwargs)
        self.t2 = t2 * si.mm
        self.gorge = gorge * si.mm
        self.l = l * si.mm
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
        """Vérifie que la soudure répond aux critères de l'EC3

        Returns:
            bool: si la soudure est correctement définie, alors True sinon False
        """
        if 60 <= self.alpha <= 120:
            # selon CNC2M-N0175-REC §3.3
            tmin = min(self.t, self.t2)
            tmax = max(self.t, self.t2)
            amin = max(3*si.mm, (sqrt(tmax)-0.5)*si.mm)
            amax = 0.7 * tmin

            if amin <= self.gorge <= amax:
                if self.l > max(30*si.mm, 6*self.gorge):
                    return True
                else:
                    print(f"La longueur du cordon de soudure est trop petite, elle doit être supérieur à {min(30*si.mm, 6*self.gorge)}")
                    return False
            else:
                print(f"La gorge doit être au minimum de {amin} et au maximum de {amax}")
                return False
        else:
            print("L'angle entre les deux pièces à souder doit être compris entre 60° et 120°")
            return False
        
    
    def beta_Lw1(self, Lj: int) -> float:
        """Calcul le facteur beta_Lw,1 qui dimminue la résistance pour des cordons de soudure des assemblages par recouvrement (à plat)

        Args:
            Lj (int): Longueur de recouvrement des plats en mm
        """
        Lj = Lj * si.mm
        return min((1.2 - 0.2 * Lj) / (150 * self.gorge), 1)

    
    def cordon_frontal(self, N_Ed: float):
        """Calcul un cordon de soudure frontale et retourne le taux de travail.

        Args:
            N_Ed (float): Effort de traction en kN.
        """
        N_Ed = N_Ed * si.kN
        return (self.beta_w * self.GAMMA_M["gamma_M2"] * (N_Ed * sqrt(2)) / self.fu) / (self.gorge * self.lef)


    def cordon_laterale(self, V_Ed: float):
        """Calcul un cordon de soudure latérale et retourne le taux de travail.

        Args:
            V_Ed (float): Effort de cisaillement du cordon en kN.
        """
        V_Ed = V_Ed * si.kN
        return (self.beta_w * self.GAMMA_M["gamma_M2"] * (V_Ed * sqrt(3)) / self.fu) / (self.gorge * self.lef)


    def cordon_oblique(self, alpha_cordon: float, N_Ed: float):
        """Calcul un cordon de soudure oblique et retourne le taux de travail.

        Args:
            N_Ed (float): Effort de traction en kN.
        """
        self.alpha_cordon = alpha_cordon
        N_Ed = N_Ed * si.kN
        return (self.beta_w * self.GAMMA_M["gamma_M2"] * (N_Ed * sqrt(3 - sin(radians(self.alpha_cordon))**2)) / self.fu) / (self.gorge * self.lef)


    def cordon_pieces_obliques(self, N_Ed: float):
        """Calcul un cordon de soudure sur des pièces à positionnement obliques et retourne le taux de travail.

        Args:
            N_Ed (float): Effort de traction en kN.
        """
        N_Ed = N_Ed * si.kN
        if self.alpha < 90:
            return (self.beta_w * self.GAMMA_M["gamma_M2"] * (N_Ed * sqrt(2 - sin(radians(self.alpha)))) / self.fu) / (self.gorge * self.lef)
        elif self.alpha > 90:
            return (self.beta_w * self.GAMMA_M["gamma_M2"] * (N_Ed * sqrt(2 + sin(radians(self.alpha)))) / self.fu) / (self.gorge * self.lef)


    def critere_generale(self, FvEd: float, FaxEd: float):
        """Calcul le critère générale de Von Mises d'une soudure et retourne le taux.

        Args:
            FvEd (int | float): Effort de cisaillement sur la en kN
            FaxEd (int | float): Effort de traction sur la soudure en kN

        Returns:
            (float) : Taux de travail de la soudure
        """
        Fv_Ed = FvEd * si.kN
        Fax_Ed = FaxEd * si.kN
        self.tau_para = Fv_Ed / (self.gorge * self.lef)
        self.sigma_perpend = (Fax_Ed * cos(radians(self.alpha/2)))/ (self.gorge * self.lef)
        self.tau_perpend = (Fax_Ed * cos(radians(self.alpha/2)))/ (self.gorge * self.lef)
        self.sigma_e = sqrt(self.sigma_perpend**2 + 3 * (self.tau_perpend**2 + self.tau_para**2))

        self.sigma_Rd = self.fu / (self.beta_w * self.GAMMA_M["gamma_M2"])
        self.sigma_perpend_Rd = (0.9 * self.fu) / self.GAMMA_M["gamma_M2"]
        return max(self.sigma_e / self.sigma_Rd, self.sigma_perpend / self.sigma_perpend_Rd)
    
    

    def soudure_discontinue(self, b: int, b1: int, t1: int, corrosion: bool=("False", "True")):
        """_summary_

        Args:
            b (int): voir EC3 1-8
            b1 (int): hauteur en mm de la pièce 2 soudé sur la pièce 1
            t (int): définit dans la classe Element
            t2 (int): épaisseur en mm de la piece 2 soudé sur pièce 1
            corrosion (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if corrosion:
            print("Il n'est pas possible d'avoir une soudure discontinue en ambiance corrosive")
            return False
        b = b * si.mm
        lwe = max(0.75 * b, 0.75 * b1*si.mm)
        l1 = min(16 * self.t, 16 * self.t2, 200)
        l2 = min(12 * self.t, 12 * self.t2, 0.25 * b, 200)
        return {"Lwe": lwe, "L1": l1, "L2": l2}




if __name__ == "__main__":
    ele = Element(classe_acier="S235", t=6, h=200, classe_transv=1)
    soudure = Soudure._from_parent_class(ele, t2=6, gorge=4, l=140, retour_soudure=True, alpha=90)
    bl = Tige._from_parent_class(ele, d=12,d0=14,qualite="A2-50",verif_filetage=True, filetage_EN1090=True)
    print(soudure.critere_generale(0, 100135))
    N_rd = bl.Veff_Rd(300, 300, "Centré")
    V_Rd = bl.Veff_Rd(300, 300, "Centré")
    taux = bl.taux_Veff_d(12, N_rd[1], 25, V_Rd[1])
    print(taux)