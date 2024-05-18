# coding in UTF-8 
# by Anthony PARISOT

############# Le but de ce fichier est de regrouper toute les fonctions du BLC dans l'EN-1995 #############
import os
import sys

import math as mt

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from EC5_Element_droit import Flexion


class Poutre_simple_decroissance(Flexion):
    """ Défini une classe poutre à simple décroissance hérité à partir de la classe Beam du fichier EC5_Element_droit.py.
        Avec pour argument :
            alpha : angle de la fibre coupé par rapport à la fibre neutre """

    def __init__(self, alpha, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.alpha = alpha


    def km_alpha(self, type_contrainte=("Traction", "Compression"), loadtype=Flexion.LOAD_TIME, typecombi=Flexion.TYPE_ACTION):
        fvd = self.f_type_d("fvk", loadtype, typecombi)[1]
        ft90d = self.f_type_d("ft90k", loadtype, typecombi)[1]
        fc90d = self.f_type_d("fc90k", loadtype, typecombi)[1]
        fmd = self.f_m_d(loadtype, typecombi)[1]

        if type_contrainte == "Traction":
            self.Km_alpha = 1 / mt.sqrt(1 + ((fmd / (0.75 * fvd)) * mt.tan(mt.radians(self.alpha)))**2 + ((fmd / ft90d) * mt.tan(mt.radians(self.alpha))**2)**2)
        else:
            self.Km_alpha = 1 / mt.sqrt(1 + ((fmd / (1.5 * fvd)) * mt.tan(mt.radians(self.alpha)))**2 + ((fmd / fc90d) * mt.tan(mt.radians(self.alpha))**2)**2)
        return self.Km_alpha
    
    # N'est plus à jour à modifier !
    def taux_m_alpha_d(self, taux_c_0_dy: float=0, taux_c_0_dz: float=0, taux_t_0_d: float=0):
        """ Retourne le taux de travail de la compression perpendiculaire """
        self.f_type_rd = self.Km_alpha * self.f_type_rd
        self.taux_m_alpha_rd = {}
        self.taux_m_alpha_rd['equ6.38y'] = (self.sigma_m_rd["y"] / (self.f_type_rd))
        self.taux_m_alpha_rd['equ6.38z'] = (self.sigma_m_rd["z"] / (self.f_type_rd))
        return self.taux_m_d(taux_c_0_dy, taux_c_0_dz, taux_t_0_d), self.taux_m_alpha_rd
    
if __name__ == "__main__":
    Beam = Poutre_simple_decroissance(b=100, h=200, lo=4000, coeflef=0.9, pos=0, alpha=2.86, section="Rectangulaire", Hi=12, classe="GL24h", cs=2)
    print(Beam.km_alpha())
    