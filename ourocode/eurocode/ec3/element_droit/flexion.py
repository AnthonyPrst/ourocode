# coding in UTF-8
# by Anthony PARISOT
import math as mt
from math import sqrt

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec3.element_droit.plat import Plat
from ourocode.eurocode.ec3.element_droit.cisaillement import Cisaillement


class Flexion(Plat):
    def __init__(self, W: si.mm**3, *args, **kwargs):
        """Defini une classe permettant le calcul en flexion d'un élément acier suivant l'EN 1993-1-1 §6.2.5.
        Cette classe est hérité de la classe Plat du module EC3_Element_droit.py.

        Args:
            Wpl (float, optional): Module de flexion plastique (pour les sections transversales de classe 1 et 2) de la section en mm3.
            Wel_min (float, optional): Module de flexion élastique (pour les sections transversales de classe 3) de la section en mm3.
            Weff_min (float, optional): Module de flexion efficace (pour les sections transversales de classe 4) de la section en mm3.
        """
        super().__init__(*args, **kwargs)
        self.W = W * si.mm**3

    @property
    def Mc_Rd(self):
        """Retourne la résistance au moment fléchissant M_c,Rd en N.mm selon EN 1993-1-1 §6.2.5.

        Sélectionne la formule selon la classe transversale :
        - Classe 1 ou 2 : M_pl,Rd = W_pl × f_y / γ_M0 (eq. 6.13).
        - Classe 3 : M_el,Rd = W_el,min × f_y / γ_M0 (eq. 6.14).
        - Classe 4 : M_c,Rd = W_eff,min × f_y / γ_M0 (eq. 6.15).

        Returns:
            tuple: (latex_string, M_c_Rd) où M_c_Rd est la résistance en flexion en N.mm (avec unité si.N×si.mm).
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
            
        
    def Mc_V_Rd(self, Av: si.mm**2, V_Ed: si.kN):
        """Retourne la résistance réduite au moment fléchissant avec cisaillement M_V,Rd selon EN 1993-1-1 §6.2.8.

        Si V_Ed > 0.5 × V_pl,Rd, le moment résistant est réduit par le facteur ρ (eq. 6.29).
        Sinon, M_V,Rd = M_c,Rd (pas d'interaction).

        Args:
            Av (float): Aire de cisaillement de la section en mm².
            V_Ed (float): Effort de cisaillement de calcul en kN.

        Returns:
            tuple: (latex_string, M_V_Rd) où M_V_Rd est le moment réduit en N.mm (avec unité si.N×si.mm).
        """
        V_Ed = V_Ed * si.kN
        cis = Cisaillement._from_parent_class(self, Av=Av)
        Vpl_Rd = cis.Vpl_Rd[1]
        Mc_Rd = self.Mc_Rd[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            if V_Ed/Vpl_Rd > 0.5:
                rho = (2*(V_Ed/Vpl_Rd)-1)**2
                My_V_Rd = (1-rho) * Mc_Rd
            elif V_Ed/Vpl_Rd <= 0.5:
                My_V_Rd = Mc_Rd
            return My_V_Rd 
        return val()


