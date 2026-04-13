# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import matplotlib.pyplot as plt

import math as mt
from math import sqrt, pi, cos, sin, radians
import numpy as np

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.element_droit.barre import Barre


class Traction(Barre):
    """Classe de vérification des éléments bois en traction axiale selon l'EC5 §6.1.2.

    Effectue les calculs de résistance et de contrainte en traction axiale selon
    l'Eurocode 5 - Partie 1-1. Hérite de Barre pour les caractéristiques
    géométriques et mécaniques.

    La vérification principale est le taux de travail en traction (équation 6.1):
    σ_t,0,d / (f_t,0,d · k_h) ≤ 1
    """

    def __init__(self, *args, **kwargs):
        """Initialise un objet de vérification en traction axiale.

        Hérite de toutes les caractéristiques de Barre (section, classe de bois,
        etc.). Aucun paramètre supplémentaire requis à l'initialisation.

        Args:
            *args: Arguments positionnels transmis à Barre.
            **kwargs: Arguments nommés transmis à Barre (b, h, classe, etc.).
        """
        super().__init__(*args, **kwargs)


    @property
    def K_h(self):
        """ Retourne le coef. Kh qui peut augmenter la resistance caractéristique fm,k et ft,k """
        return self._K_h()


    def f_t_0_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul en traction axiale f_t,0,d selon l'EC5 §6.1.2.

        Détermine la résistance à partir de la résistance caractéristique ft,0,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul ft,0,d en MPa avec unité (si.MPa).
        """
        return super()._f_type_d("ft0k", loadtype, typecombi)


    def sigma_t_0_d(self, Ft0d: si.kN, Anet: si.mm**2=None) -> tuple:
        """Calcule la contrainte de traction axiale sigma_t,0,d selon l'EC5 §6.1.2.

        Détermine la contrainte normale due à l'effort de traction axial.
        Prend en compte une section nette réduite (perçages, entailles) si spécifiée.

        Args:
            Ft0d (si.kN): Effort de traction axial en kN.
            Anet (si.mm**2, optional): Aire nette de la section en mm² si réduction
                (perçages, entailles). Doit être ≤ aire brute. Defaults to None.

        Returns:
            tuple: (latex_string, valeur) où valeur est sigma_t,0,d en MPa avec unité.

        Raises:
            ValueError: Si Anet > aire brute de la section.

        Note:
            La valeur est stockée dans l'attribut sigma_t_0_rd.
            Pour les assemblages boulonnés, utiliser Anet pour tenir compte des trous.
        """
        self.Ft_0_d = Ft0d * si.kN
        Ft_0_d = self.Ft_0_d
        if Anet and Anet * si.mm**2<= self.aire:
            A = Anet * si.mm**2
        else:
            A = self.aire

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            sigma_t_0_d = Ft_0_d / A
            return sigma_t_0_d
        value = val()
        self.sigma_t_0_rd = value[1]
        return value


    def taux_t_0_d(self) -> tuple:
        """Calcule le taux de travail en traction axiale selon l'EC5 §6.1.2 (Eq. 6.1).

        Vérifie le critère : σ_t,0,d / (k_h · f_t,0,d) ≤ 1

        Le coefficient k_h (effet de hauteur) est pris comme le minimum des
        valeurs selon y et z pour être conservateur.

        Returns:
            tuple: (latex_string, valeur) où valeur est le taux en pourcentage
                (0.75 = 75%). Stocké dans taux_t_0_rd['equ6.1'].

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_t_0_rd = {}
        K_h_y = self.K_h['y']
        K_h_z = self.K_h['z']
        sigma_t_0_d = self.sigma_t_0_rd
        f_t_0_d = self.f_type_rd

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            K_h = min(K_h_y, K_h_z)
            taux_6_1 = sigma_t_0_d / (K_h * f_t_0_d) # equ6.1
            return taux_6_1
        value = val()

        self.taux_t_0_rd['equ6.1'] = value[1]
        synthese = [
            ["Traction bois", self.taux_t_0_rd['equ6.1'], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value


