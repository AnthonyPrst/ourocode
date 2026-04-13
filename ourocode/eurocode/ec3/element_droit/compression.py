# coding in UTF-8
# by Anthony PARISOT
import math as mt
from math import sqrt

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec3.element_droit.plat import Plat


class Compression(Plat):
    FACTEUR_ALPHA = {"a0": 0.13, "a": 0.21, "b": 0.34, "c": 0.49, "d": 0.76}
    COEF_LF = {"Encastré 1 côté" : 2,
                "Rotule - Rotule" : 1,
                "Encastré - Rotule" : 0.7,
                "Encastré - Encastré" : 0.5,
                "Encastré - Rouleau" : 1}
    def __init__(self, A: si.mm**2, Iy:si.mm**4, Iz:si.mm**4, lo_y: si.mm=0, lo_z: si.mm=0, courbe_flamb: dict="{'y':'c', 'z':'c'}", type_appuis: str=COEF_LF, *args, **kwargs):
        """
        Classe intégrant les formules de compression et d'instabilité au flambement à l'EC3.
        Cette classe est hérité de la classe Plat du module EC3_Element_droit.py.

        Args:
            A (float): Aire brute si classe 1,2 ou 3 et Aeff si classe 4 en mm²
            Iy (float): Moment quadratique suivant l'axe de rotation y en mm4
            Iz (float): Moment quadratique suivant l'axe de rotation z en mm4
            lo (int): Longueur de flambement suivant l'axe de rotation (y ou z) en mm. Defaults to {'y':0, 'z':0}.
            type_appuis (str): Permet de déterminé la forme du flambement en fonction des types d'appui:
                Encastré 1 côté : 2
                Rotule - Rotule : 1
                Encastré - Rotule : 0.7
                Encastré - Encastré : 0.5
                Encastré - Rouleau : 1
        """
        super().__init__(*args, **kwargs)
        self.A = A * si.mm**2
        self.Iy = Iy * si.mm**4
        self.Iz = Iz * si.mm**4
        self.lo ={'y': lo_y*si.mm, 'z': lo_z*si.mm}
        self.courbe_flamb = courbe_flamb
        self.type_appuis = type_appuis
        self.coef_lef = self.COEF_LF[type_appuis]


    
    @property
    def Nc_Rd(self):
        """Calcul la résistance en compression de la section transversale en N (équa 6.10 et 6.11)
        """
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        A = self.A
        fy = self.fy
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            N_c_Rd = A * fy/gamma_M0
            return N_c_Rd
        return val()

    @property
    def lamb(self):
        """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb = {'y':0, 'z':0}
        lamb['y'] = (self.lo['y'].value * 10**3 * self.coef_lef) / mt.sqrt(self.Iy / (self.A))
        lamb['z'] = (self.lo['z'].value * 10**3 * self.coef_lef) / mt.sqrt(self.Iz / (self.A))
        return lamb

    @property
    def _lamb_rel_Axe(self):
        """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
        lamb_rel_Axe = {'y':0, 'z':0}
        for cle, value in lamb_rel_Axe.items():
            lamb_rel_Axe[cle] = (self.lamb[cle] / mt.pi) * mt.sqrt(self.fy / self.E)
        return lamb_rel_Axe

    @property
    def _alpha(self):
        """Détermine le facteur d'imperfection fonction des courbes de flambement
        """
        a = {}
        for key, value in self.courbe_flamb.items():
            a[key] = self.FACTEUR_ALPHA[value]
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
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        f_y = self.fy
        A = self.A
        chi_y = self._chi['y']
        chi_z = self._chi['z']
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            N_b_y_Rd = chi_y * A * f_y / gamma_M0
            N_b_z_Rd = chi_z * A * f_y / gamma_M0
            return {"y": N_b_y_Rd, "z": N_b_z_Rd}
        return val()

