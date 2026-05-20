# coding in UTF-8
# by Anthony PARISOT
import math as mt
from math import sqrt

import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.ec3.element_droit.plat import Plat


class Traction(Plat):
    def __init__(self, A: si.mm**2, Anet: si.mm**2=0, ass_cat_C: bool=("False", "True"), *args, **kwargs):
        """Défini une classe traction permettant le calcul d'un élément métallique à l'EN 1993-1-1 §6.2.3.
        Cette classe est hérité de la classe Plat du module EC3_Element_droit.py.

        Args:
            A (float): Aire brute de la section en mm².
            Anet (float, optional): Aire nette au droit droit des trous de fixation selon §6.2.2.2 en mm².
            ass_cat_C (bool, optional): Si assemblage de catégorie C alors True sinon False, voir EN 1993-1-8 §3.4.1.(1). Defaults to False.
        """
        super().__init__(*args, **kwargs)
        self.A = A * si.mm**2
        self.Anet = Anet * si.mm**2
        self.ass_cat_C = ass_cat_C 

    @property
    def _Npl_Rd(self):
        """Calcul la résistance plastique en traction de la section transversale brute en N (équa 6.6)
        """
        A = self.A
        fy = self.fy
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            N_pl_Rd = (A * fy)/gamma_M0 #(équa 6.6)
            return N_pl_Rd
        return val()

    @property
    def _Nu_Rd(self):
        """Calcul la résistance ultime en traction de la section transversale nette en N (équa 6.7)
        """
        Anet = self.Anet
        fu = self.fu
        gamma_M2 = self.GAMMA_M["gamma_M2"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            N_u_Rd = 0.9*(Anet * fu)/gamma_M2 #(équa 6.7)
            return N_u_Rd
        return val()

    @property
    def _Nnet_Rd(self):
        """Calcul la résistance ultime en traction de la section transversale nette en N
        lorsque l'assemblage est de catégorie C (équa 6.8)
        """
        Anet = self.Anet
        fy = self.fy
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            N_net_Rd = (Anet * fy)/gamma_M0 #(équa 6.8)
            return N_net_Rd
        return val()
    
    @property
    def Nt_Rd(self):
        """Retourne la résistance de la section transversale en traction N_t,Rd en N selon EN 1993-1-1 §6.2.3.

        Sélectionne la formule appropriée selon la présence de trous et la catégorie d'assemblage :
        - Assemblage cat. C (glissement à l'ELU) : N_net,Rd = A_net × f_y / γ_M0 (eq. 6.8).
        - Avec trous : N_t,Rd = min(N_pl,Rd, N_u,Rd) (eq. 6.6 et 6.7).
        - Sans trous : N_t,Rd = N_pl,Rd (eq. 6.6).

        Returns:
            tuple: (latex_string, N_t_Rd) où N_t_Rd est la résistance en traction en N (avec unité si.N).
        """
        if self.ass_cat_C:
            return self._Nnet_Rd
        if self.Anet:
            latex_Npl_Rd, N_pl_Rd = self._Npl_Rd
            latex_Nu_Rd, N_u_Rd = self._Nu_Rd
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                N_t_Rd = min(N_pl_Rd, N_u_Rd)
                return N_t_Rd
            value = val()
            return (latex_Npl_Rd + latex_Nu_Rd + value[0], value[1])
        else:
            return self._Npl_Rd  
        

