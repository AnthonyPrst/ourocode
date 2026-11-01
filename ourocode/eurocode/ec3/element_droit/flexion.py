# coding in UTF-8
# by Anthony PARISOT
import logging
import math as mt
from math import sqrt

import forallpeople as si

si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.ec3.element_droit.plat import Plat
from ourocode.eurocode.ec3.element_droit.cisaillement import Cisaillement

logger = logging.getLogger(__name__)


class Flexion(Plat):
    def __init__(self, W: si.mm**3 = None, axe: str = ("y", "z"), *args, **kwargs):
        """Defini une classe permettant le calcul en flexion d'un élément acier suivant l'EN 1993-1-1 §6.2.5.
        Cette classe est hérité de la classe Plat du module EC3_Element_droit.py.

        Args:
            W (float, optional): Module de flexion de la section en mm³.
                Utiliser Wpl pour les sections de classe 1 ou 2, Wel_min pour la classe 3, Weff_min pour la classe 4.
                Si None, il est déduit de la section selon la classe transversale et l'axe
                (Wpl pour les classes 1 et 2, Wel pour la classe 3). Defaults to None.
            axe (str, optional): Axe de flexion ("y" : axe fort, "z" : axe faible). Defaults to "y".
        """
        super().__init__(*args, **kwargs)
        self.axe = self._check_choice(axe, ("y", "z"), "axe")
        if W is not None:
            self.W = W * si.mm**3
            return
        W_el = getattr(self, f"Wel{self.axe}")
        W_pl = getattr(self, f"Wpl{self.axe}")
        if self.classe_transv in (1, 2) and W_pl is None and W_el is not None:
            logger.warning("Wpl%s non disponible pour %s : Wel%s retenu.", self.axe, self.designation, self.axe)
        self.W = W_pl if self.classe_transv in (1, 2) and W_pl is not None else W_el
        if self.W is None:
            raise ValueError("W non défini : renseigner W ou la géométrie (t, b ou designation).")

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
            case 1 | 2:
                W_pl = self.W

                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    M_pl_Rd = W_pl * f_y / gamma_M0  # équa 6.13
                    return M_pl_Rd
            case 3:
                W_el_min = self.W

                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    M_el_Rd = W_el_min * f_y / gamma_M0  # équa 6.14
                    return M_el_Rd
            case 4:
                W_eff_min = self.W

                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    M_c_Rd = W_eff_min * f_y / gamma_M0  # équa 6.15
                    return M_c_Rd

        return val()

    def Mc_V_Rd(self, Av: si.mm**2 = None, V_Ed: si.kN = 0):
        """Retourne la résistance réduite au moment fléchissant avec cisaillement M_V,Rd selon EN 1993-1-1 §6.2.8.

        Si V_Ed > 0.5 × V_pl,Rd, le moment résistant est réduit par le facteur ρ (eq. 6.29).
        Sinon, M_V,Rd = M_c,Rd (pas d'interaction).

        Args:
            Av (float, optional): Aire de cisaillement de la section en mm². Si None, Avz
                (flexion autour de y) ou Avy (flexion autour de z) de la section. Defaults to None.
            V_Ed (float): Effort de cisaillement de calcul en kN. Defaults to 0.

        Returns:
            tuple: (latex_string, M_V_Rd) où M_V_Rd est le moment réduit en N.mm (avec unité si.N×si.mm).
        """
        V_Ed = V_Ed * si.kN
        cis = Cisaillement._from_parent_class(self, Av=Av, axe_effort="z" if self.axe == "y" else "y")
        Vpl_Rd = cis.Vpl_Rd[1]
        Mc_Rd = self.Mc_Rd[1]

        if V_Ed / Vpl_Rd > 0.5:

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                rho = (2 * (V_Ed / Vpl_Rd) - 1) ** 2
                My_V_Rd = (1 - rho) * Mc_Rd
                return My_V_Rd
        else:

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                My_V_Rd = Mc_Rd
                return My_V_Rd

        return val()
