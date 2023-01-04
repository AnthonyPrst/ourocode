#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import os
import sys

import math as mt
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Project

class Element(Project):

    GAMMA_M = {"gamma_M0": 1, "gamma_M1": 1, "gamma_M2": 1.25, "gamma_M3": 1.1, "gamma_M3_ser": 1.25,
             "gamma_M4": 1, "gamma_M5": 1, "gamma_M6_ser": 1, "gamma_M7": 1.1}
    E = 210000
             
    def __init__(self, t: int=0, h: int=0, classe_acier: str="S235", classe_transv: int=1, **kwargs):
        """Configure un objet Element pour vérifier un élément acier suivant l'EN 1993-1-1. 

        Args:
            t (int, optional): épaisseur de la plaque en mm. Defaults to 0.
            h (int, optional): hauteur de la plaque en mm. Defaults to 0.
            classe_acier (str, optional): _description_. Defaults to "S235".
            classe_transv (int, optional): _description_. Defaults to 1.
        """
        super().__init__(**kwargs)
        self.t = t
        self.h = h
        self.classe_acier = classe_acier
        self.classe_transv = classe_transv
        self.__fy_fu()
        for key, value in kwargs.items():
            setattr(self, key, value)
    

    def __data_from_csv(self, data_file: str):
            """ Retourne un dataframe d'un fichier CSV """
            repertory = os.path.join(os.getcwd(), "data", data_file)
            data_csv = pd.read_csv(repertory, sep=';', index_col=0)
            return data_csv

    @property
    def __classe_acier(self):
        """Retourne le dataframe de la classe d'acier défini 
        """
        df = self.__data_from_csv("caracteristique_meca_acier.csv")
        df = df.loc[self.classe_acier]
        return df


    def __fy_fu(self):
        """Défini fy (résistance élastique) et fu (résistance plastique) en MPa fonction de la classe d'acier choisi
        """
        if self.t <= 40:
            self.fy = self.__classe_acier.loc["t<= 40  fy"]
            self.fu = self.__classe_acier.loc["t<= 40  fu"]
        elif self.t > 40 and self.t <= 80:
            self.fy = self.__classe_acier.loc["40<t<= 80  fy"]
            self.fu = self.__classe_acier.loc["40<t<= 80  fu"]

    
    @property
    def _inertie(self):
        """ Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
            b ou d : Largeur ou diamètre de la poutre en mm
            h : Hauteur de la poutre en mm """
        if self.t and self.h:
            iy = (self.t * self.h**3)/12
            iz = (self.h * self.t**3)/12
            return [iy, iz]

        elif self.Iy and self.Iz:
            return [self.Iy, self.Iz]



class Traction(Element):
    def __init__(self, A: float|int, Anet: float|int=0, ass_cat_C: bool=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.A = A
        self.Anet = Anet
        self.ass_cat_C = ass_cat_C

    @property
    def _Npl_Rd(self):
        """Calcul la résistance plastique en traction de la section transversale brute en N (équa 6.6)
        """
        return (self.A * self.fy)/__class__.GAMMA_M["gamma_M0"]

    @property
    def _Nu_Rd(self):
        """Calcul la résistance ultime en traction de la section transversale nette en N (équa 6.7)
        """
        return 0.9*(self.Anet * self.fu)/__class__.GAMMA_M["gamma_M2"]

    @property
    def _Nnet_Rd(self):
        """Calcul la résistance ultime en traction de la section transversale nette en N
        lorsque l'assemblage est de catégorie C (équa 6.8)
        """
        return (self.Anet * self.fy)/__class__.GAMMA_M["gamma_M0"]
    
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
    def __init__(self, A: float|int, lo={'y':0, 'z':0}, courbe_flamb={'y':'c', 'z':'c'}, coeflf=1, *args, **kwargs):
        """
        Args:
            A (float | int): Aire brute si classe 1,2 ou 3 et Aeff si classe 4
            lo (dict, optional): Longueur de flambement suivant l'axe de rotation (y ou z) en mm. Defaults to {'y':0, 'z':0}.
            coeflf (int, optional): Coefficient multiplicateur de la longueur pour obtenir la longeur efficace de flambement en
                                    fonction des du type d'appuis :
                                                Encastré 1 côté : 2
                                                Rotule - Rotule : 1
                                                Encastré - Rotule : 0.7
                                                Encastré - Encastré : 0.5
                                                Encastré - Rouleau : 1. Defaults to 1.
        """
        super().__init__(*args, **kwargs)
        self.A = A
        self.lo = lo
        self.courbe_flamb = courbe_flamb
        self.coeflf = coeflf


    
    @property
    def Nc_Rd(self):
        """Calcul la résistance en compression de la section transversale en N (équa 6.10 et 6.11)
        """
        return (self.A * self.fy)/__class__.GAMMA_M["gamma_M0"]

    @property
    def _lamb(self):
        """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb = {'y':0, 'z':0}
        lamb['y'] = (self.lo['y'] * self.coeflf) / mt.sqrt(self._inertie[0] / (self.A))
        lamb['z'] = (self.lo['z'] * self.coeflf) / mt.sqrt(self._inertie[1] / (self.A))
        return lamb

    @property
    def _lamb_rel_Axe(self):
        """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb_rel_Axe = {'y':0, 'z':0}
        for cle, value in lamb_rel_Axe.items():
            lamb_rel_Axe[cle] = (self._lamb[cle] / mt.pi) * mt.sqrt(self.fy / __class__.E)
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
    def __init__(self, Av: int|float, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Av = Av

    @property
    def Vpl_Rd(self):
        """Calcul la résistance du cisaillment plastique en N (équa 6.18)
        """
        return (self.Av * (self.fy/mt.sqrt(3)))/__class__.GAMMA_M["gamma_M0"]       



class Flexion(Element):
    def __init__(self, Wpl: float|int=0, Wel_min: float|int=0, Weff_min: float|int=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Wpl = Wpl
        self.Wel_min = Wel_min
        self.Weff_min = Weff_min

    @property
    def Mc_Rd(self):
        """Calcul la résistance du moment fléchissant de la section transversale en N (équa 6.13, 6.14, et 6.15)
        """
        match self.classe_transv:
            case (1|2):
                return (self.Wpl* self.fy)/__class__.GAMMA_M["gamma_M0"]
            case 3:
                return (self.Wel_min * self.fy)/__class__.GAMMA_M["gamma_M0"]
            case 4:
                return (self.Weff_min * self.fy)/__class__.GAMMA_M["gamma_M0"]
            
        
  
    

if __name__ == "__main__":
    aire = 10*112
    calcul = Compression(t=10, h=112, classe_acier="S235", classe_transv=1, A=aire, lo={'y':228, 'z':228}, coeflf=2)
    print(calcul.fy)