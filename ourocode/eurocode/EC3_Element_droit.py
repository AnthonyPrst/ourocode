#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import math as mt
from math import sqrt
import pandas as pd

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.A0_Projet import Projet

class Element(Projet):

    GAMMA_M = {"gamma_M0": 1, "gamma_M1": 1, "gamma_M2": 1.25, "gamma_M3": 1.1, "gamma_M3_ser": 1.25,
             "gamma_M4": 1, "gamma_M5": 1, "gamma_M6_ser": 1, "gamma_M7": 1.1}
    E = 210000 * si.MPa
    CLASSE_STEEL = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_acier.csv").index)
             
    def __init__(self, t: si.mm=0, h: si.mm=0, classe_acier: str=CLASSE_STEEL, classe_transv: int=("1","2","3","4"), **kwargs):
        """Configure un objet Element pour vérifier un élément acier suivant l'EN 1993-1-1. 

        Args:
            t (int, optional): épaisseur de la plaque en mm. Defaults to 0.
            h (int, optional): hauteur de la plaque en mm. Defaults to 0.
            classe_acier (str, optional): classe d'acier. Defaults to "S235".
            classe_transv (int, optional): classe transversale de la section en fonction de sa capacité de plastification. Defaults to 1.
        """
        super().__init__(**kwargs)
        self.t = t * si.mm
        self.h = h * si.mm
        self.classe_acier = classe_acier
        self.classe_transv = classe_transv
        self.__fy_fu()
    

    @property
    def __classe_acier(self):
        """Retourne le dataframe de la classe d'acier définit 
        """
        df = self._data_from_csv("caracteristique_meca_acier.csv")
        df = df.loc[self.classe_acier]
        return df


    def __fy_fu(self):
        """Défini fy (résistance élastique) et fu (résistance plastique) en MPa fonction de la classe d'acier choisi
        """
        print(self.t)
        if self.t <= 40*si.mm:
            self.fy = self.__classe_acier.loc["t<= 40  fy"] * si.MPa
            self.fu = self.__classe_acier.loc["t<= 40  fu"] * si.MPa
        elif self.t > 40*si.mm and self.t<= 80*si.mm:
            self.fy = self.__classe_acier.loc["40<t<= 80  fy"] * si.MPa
            self.fu = self.__classe_acier.loc["40<t<= 80  fu"] * si.MPa

    
    @property
    def _inertie(self):
        """ Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
            b ou d : Largeur ou diamètre de la poutre en mm
            h : Hauteur de la poutre en mm """
        if self.t and self.h:
            self.Iy = (self.t * self.h**3)/12
            self.Iz = (self.h * self.t**3)/12
            return [self.Iy, self.Iz]

        elif self.Iy and self.Iz:
            return [self.Iy, self.Iz]



class Traction(Element):
    def __init__(self, A: si.mm*2, Anet: si.mm*2=0, ass_cat_C: bool=("False", "True"), *args, **kwargs):
        """Défini une classe traction permettant le calcul d'un élément métallique à l'EN 1993-1-1 §6.2.3.
        Cette classe est hérité de la classe Element du module EC3_Element_droit.py.

        Args:
            A (float): Aire brute de la section en mm².
            Anet (float, optional): Aire nette au droit droit des trous de fixation selon §6.2.2.2 en mm². Defaults to 0.
            ass_cat_C (bool, optional): Si assemblage de catégorie C alors True sinon False, voir EN 1993-1-8 §3.4.1.(1). Defaults to ("False", "True").
        """
        super().__init__(*args, **kwargs)
        self.A = A * si.mm*2
        self.Anet = Anet * si.mm*2
        self.ass_cat_C = ass_cat_C 

    @property
    def _Npl_Rd(self):
        """Calcul la résistance plastique en traction de la section transversale brute en N (équa 6.6)
        """
        return (self.A * self.fy)/__class__.GAMMA_M["gamma_M0"] #(équa 6.6)

    @property
    def _Nu_Rd(self):
        """Calcul la résistance ultime en traction de la section transversale nette en N (équa 6.7)
        """
        return 0.9*(self.Anet * self.fu)/__class__.GAMMA_M["gamma_M2"] #(équa 6.7)

    @property
    def _Nnet_Rd(self):
        """Calcul la résistance ultime en traction de la section transversale nette en N
        lorsque l'assemblage est de catégorie C (équa 6.8)
        """
        return (self.Anet * self.fy)/__class__.GAMMA_M["gamma_M0"] #(équa 6.8)
    
    @property
    def Nt_Rd (self):
        if self.ass_cat_C:
            return self._Nnet_Rd
        if self.Anet:
           return min(self._Npl_Rd, self._Nu_Rd)
        else:
            return self._Npl_Rd  
        


class Compression(Element):
    """Classe intégrant les formules de compression et d'instabilité au flambement à l'EC3

    Args:
        Element (class): hérite des propriétés de la classe Element à l'EC3
    """

    FACTEUR_ALPHA = {"a0": 0.13, "a": 0.21, "b": 0.34, "c": 0.49, "d": 0.76}
    def __init__(self, A: si.mm**2, lo_y: si.mm=0, lo_z: si.mm=0, courbe_flamb: dict="{'y':'c', 'z':'c'}", coeflf: float=1, *args, **kwargs):
        """
        Args:
            A (float | int): Aire brute si classe 1,2 ou 3 et Aeff si classe 4 en mm²
            lo (int, optional): Longueur de flambement suivant l'axe de rotation (y ou z) en mm. Defaults to {'y':0, 'z':0}.
            coeflf (float, optional): Coefficient multiplicateur de la longueur pour obtenir la longeur efficace de flambement en
                fonction des du type d'appuis :
                            Encastré 1 côté : 2
                            Rotule - Rotule : 1
                            Encastré - Rotule : 0.7
                            Encastré - Encastré : 0.5
                            Encastré - Rouleau : 1. Defaults to 1.
        """
        super().__init__(*args, **kwargs)
        self.A = A * si.mm**2
        self.lo ={'y': lo_y*si.mm, 'z': lo_z*si.mm}
        self.courbe_flamb = courbe_flamb
        self.coeflf = coeflf


    
    @property
    def Nc_Rd(self):
        """Calcul la résistance en compression de la section transversale en N (équa 6.10 et 6.11)
        """
        return (self.A * self.fy)/__class__.GAMMA_M["gamma_M0"]

    @property
    def lamb(self):
        """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb = {'y':0, 'z':0}
        lamb['y'] = (self.lo['y'].value * 10**3 * self.coeflf) / mt.sqrt(self._inertie[0] / (self.A))
        lamb['z'] = (self.lo['z'].value * 10**3 * self.coeflf) / mt.sqrt(self._inertie[1] / (self.A))
        return lamb

    @property
    def _lamb_rel_Axe(self):
        """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb_rel_Axe = {'y':0, 'z':0}
        for cle, value in lamb_rel_Axe.items():
            lamb_rel_Axe[cle] = (self.lamb[cle] / mt.pi) * mt.sqrt(self.fy / __class__.E)
        return lamb_rel_Axe

    @property
    def _alpha(self):
        """Détermine le facteur d'imperfection fonction des courbes de flambement
        """
        a = {}
        for key, value in self.courbe_flamb.items():
            a[key] = __class__.FACTEUR_ALPHA[value]
        return a

    @property
    def _phi(self):
        """Détermine le facteur phi (fonction de l'axe de flambement)
        """
        phi = {'y':0, 'z':0}
        for key in phi.keys():
            phi[key] = 0.5 * (1 + self._alpha[key] * (self._lamb_rel_Axe[key] - 0.2) + self._lamb_rel_Axe[key]**2)
        return phi

    @property
    def _chi(self):
        """Détermine le facteur ki (fonction de l'axe de flambement)
        """
        ki = {'y':0, 'z':0}
        for key in ki.keys():
            ki[key] = 1 / (self._phi[key] + mt.sqrt(self._phi[key]**2 - self._lamb_rel_Axe[key]**2))
        return ki
    
    @property
    def Nb_Rd(self):
        """Renvoie la capacité résitante en compression avec flambement en N (fonction de l'axe de flambement)
        """
        NbRd = {'y':0, 'z':0}
        for key in NbRd.keys():
            NbRd[key] = (self._chi[key] * self.A * self.fy)/__class__.GAMMA_M["gamma_M0"]
        return NbRd



class Cisaillement(Element):
    def __init__(self, Av: si.mm**2, *args, **kwargs):
        """Defini une classe permettant le calcul d'un élément métalique en cisaillement selon l'EN 1993-1-1 §6.2.6.
        Cette classe est hérité de la classe Element du module E3_Element_droit.py.

        Args:
            Av (float): Aire de cisaillemment en mm²
        """
        super().__init__(*args, **kwargs)
        self.Av = Av * si.mm**2

    @property
    def Vpl_Rd(self):
        """Calcul la résistance du cisaillment plastique en N (équa 6.18)
        """
        A_v = self.Av
        f_y = self.fy
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            Vpl_Rd = (A_v * f_y/sqrt(3)) / gamma_M0 #(équa 6.18)
            return Vpl_Rd
        return val() 


class Flexion(Element):
    def __init__(self, W: si.mm**3, *args, **kwargs):
        """Defini une classe permettant le calcul en flexion d'un élément acier suivant l'EN 1993-1-1 §6.2.5.
        Cette classe est hérité de la classe Element du module EC3_Element_droit.py.

        Args:
            Wpl (float, optional): Module de flexion plastique (pour les sections transversales de classe 1 et 2) de la section en mm3.
            Wel_min (float, optional): Module de flexion élastique (pour les sections transversales de classe 3) de la section en mm3.
            Weff_min (float, optional): Module de flexion efficace (pour les sections transversales de classe 4) de la section en mm3.
        """
        super().__init__(*args, **kwargs)
        self.W = W * si.mm**3

    @property
    def Mc_Rd(self):
        """Calcul la résistance du moment fléchissant de la section transversale en N.mm (équa 6.13, 6.14, et 6.15)
        """
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        f_y = self.fy
        match self.classe_transv:
            case 1|2:
                W_pl = self.W
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    M_pl_Rd = W_pl * f_y/ gamma_M0 #équa 6.13
                    return M_pl_Rd
            case 3:
                W_el_min = self.W
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    M_el_Rd = W_el_min * f_y/ gamma_M0 #équa 6.14
                    return M_el_Rd
            case 4:
                W_eff_min = self.W
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    M_c_Rd = W_eff_min * f_y/ gamma_M0 #équa 6.15
                    return M_c_Rd
        return val()
            
        
    def Mc_V_Rd(self, Av: float, V_Ed: float):
        """Calcul la résistance du moment fléchissant avec prise en compte de l'incidence du cisaillement sur ce dernier EN 1993-1-1 §6.2.8
            (équa 6.29)

        Args:
            Av (float): Aire de cisaillemment en mm².
            V_Ed (float): Effort de cisaillement en kN.
        """
        A_v = Av * si.mm**2
        V_Ed = V_Ed * si.kN
        cis = Cisaillement._from_parent_class(self, Av=Av)
        Vpl_Rd = cis.Vpl_Rd[1]
        Mc_Rd = self.Mc_Rd[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            if V_Ed/Vpl_Rd > 0.5:
                rho = ((2*V_Ed) / Vpl_Rd-1)**2
                My_V_Rd = rho * Mc_Rd
            elif V_Ed/Vpl_Rd <= 0.5:
                My_V_Rd = Mc_Rd
            return My_V_Rd 
        return val()


if __name__ == "__main__":
    aire = 10*112
    calcul = Flexion(t=10, h=112, classe_acier="S235", classe_transv=1,  W=2000000)
    print(calcul.Mc_V_Rd(200, 12000))
    val =10* si.mm
    print(val.value)