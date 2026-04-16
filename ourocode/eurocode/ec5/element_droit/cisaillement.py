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
from ourocode.eurocode.core.renderer import handcalc

from ourocode.eurocode.ec5.element_droit.barre import Barre


class Cisaillement(Barre):
    """Classe de vérification des éléments bois au cisaillement selon l'EC5 §6.1.7 et §6.5.

    Effectue les calculs de contrainte et taux de travail au cisaillement
    longitudinal pour des poutres en bois selon l'Eurocode 5.
    Hérite de Barre pour les caractéristiques géométriques et mécaniques.

    Vérifie :
    - La résistance au cisaillement (§6.1.7, Eq. 6.13)
    - Le cisaillement avec entaille (§6.5, Eq. 6.60) avec facteur K_v
    """

    DICT_KN = {"Massif": 5,"BLC": 6.5, "LVL": 4.5}
    def __init__(self, **kwargs):
        """Initialise un objet de vérification au cisaillement.

        Hérite de toutes les caractéristiques de Barre. Initialise K_v à 1
        (pas d'entaille) et h_ef à la hauteur totale.

        Args:
            **kwargs: Arguments transmis à Barre (classe, etc.).
        """
        super().__init__(**kwargs)
        self.K_v = 1
        self.h_ef = self.h_calcul

    @property
    def K_cr(self):
        """Facteur de réduction de largeur K_cr selon l'EC5 §6.1.7.

        Tient compte des fissures de séchage dans le bois massif et BLC.
        Réduction de 33% (K_cr = 0.67) pour les sections fragilisées.

        Returns:
            float: Valeur de K_cr.
                - 1.0 : Pas de réduction (bois sans fissuration significative)
                - 0.67 : Réduction pour bois massif h > 150mm (CS1/CS2) ou BLC (CS2/CS3)

        Note:
            CS1/CS2 : K_cr = 0.67 si h > 150mm et bois massif
            CS2 : K_cr = 0.67 pour BLC
            CS3 : K_cr = 0.67 pour tout bois
        """
        if self.cs == 1:
            if self.h_calcul.value * 10**3 > 150 and self.type_bois == "Massif":
                return 0.67
            else:
                return 1
        elif self.cs == 2:
            if self.h_calcul.value * 10**3 > 150 and self.type_bois == "Massif":
                return 0.67
            elif self.type_bois == "BLC":
                return 0.67
            else:
                return 1
        else:
            return 0.67


    def Kv(self, hef:si.mm, x:si.mm, i_lo:si.mm, ent=("Dessous", "Dessus")) -> tuple:
        """Calcule le facteur de réduction d'entaille K_v selon l'EC5 §6.5.

        Ce coefficient réduit la résistance au cisaillement en présence d'une
        entaille au niveau d'un appui (entaille en dessous ou au dessus de la poutre).

        Formule EC5 : K_v = min(1, [K_n(1 + 1.1·i^1.5/√h)] / [√h·(√(α(1-α)) + 0.8·x/h·√(1/α - α²))])
        où α = h_ef/h et i = i_lo/h_ef

        Args:
            hef (si.mm): Hauteur efficace de la poutre (hauteur - profondeur entaille) en mm.
            x (si.mm): Distance entre le centre de réaction à l'appui et le coin de l'entaille en mm.
            i_lo (si.mm): Longueur horizontale de l'entaille en mm.
            ent (str, optional): Position de l'entaille : "Dessous" ou "Dessus".
                Defaults to "Dessous".

        Returns:
            tuple: (latex_string, valeur) où valeur est le facteur K_v (≤ 1).
                Si ent="Dessus", retourne K_v = 1 (pas de réduction).

        Note:
            La valeur est stockée dans l'attribut K_v.
            K_n dépend du type de bois (voir DICT_KN).
        """
        K_n = self.DICT_KN[self.type_bois]
        x = x * si.mm
        h_ef = hef * si.mm
        i = i_lo * si.mm / h_ef
        h_calcul = self.h_calcul

        self.h_ef = h_ef
        if ent == "Dessus":
            self.K_v = 1
            return self.K_v
        else:
            @handcalc(override="long", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                alpha = h_ef / h_calcul
                K_v = min(1,(K_n * (1 + (1.1 * i ** 1.5) / sqrt(h_calcul))) / (sqrt(h_calcul) * (sqrt(alpha * (1 - alpha)) + 0.8 * x / h_calcul * sqrt(1 / alpha - alpha ** 2))))
                return K_v
            value = val()
            self.K_v = value[1]
            return value

    def f_v_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcule la résistance de calcul au cisaillement f_v,d selon l'EC5 §6.1.7.

        Détermine la résistance à partir de la résistance caractéristique f_v,k
        et des coefficients de modification (kmod, γM).

        Args:
            loadtype (str): Classe de durée de chargement.
                Voir Barre.LOAD_TIME pour les valeurs possibles.
            typecombi (str): Type de combinaison d'actions.
                "fondamentale" ou "accidentelle". Defaults to "fondamentale".

        Returns:
            float: Résistance de calcul f_v,d en MPa avec unité (si.MPa).
        """
        return super()._f_type_d("fvk", loadtype, typecombi)

    def tau_d(self, Vd:si.kN) -> tuple:
        """Calcule la contrainte de cisaillement tau_d selon l'EC5 §6.1.7.

        Détermine la contrainte de cisaillement longitudinal pour une poutre
        rectangulaire : τ = 1.5·V/(b_ef·h_ef)

        Args:
            Vd (si.kN): Effort tranchant (cisaillement) sur la poutre en kN.

        Returns:
            tuple: (latex_string, valeur) où valeur est tau_d en MPa avec unité.

        Note:
            La valeur est stockée dans l'attribut tau_rd.
            Prend en compte K_cr (fissuration) et K_v (entaille si défini).
        """
        self.V_d = Vd * si.kN
        V_d = self.V_d
        K_cr = self.K_cr
        b_calcul = self.b_calcul
        h_ef = self.h_ef

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            b_ef = K_cr * b_calcul
            tau_d = (1.5 * V_d) / (b_ef * h_ef)
            return tau_d

        value = val()
        self.tau_rd = value[1]
        return value

    def taux_tau_d(self) -> tuple:
        """Calcule les taux de travail au cisaillement selon l'EC5 §6.1.7 et §6.5.

        Vérifie les critères :
        - Equ. 6.13 : Cisaillement simple (tau_d / f_v,d)
        - Equ. 6.60 : Cisaillement avec entaille (tau_d / (K_v · f_v,d))

        Returns:
            tuple: (latex_string, valeurs) où valeurs contient les taux pour
                equ6.13 (sans entaille) et equ6.60 (avec entaille) en pourcentage.
                Stockés dans taux_tau_rd['equ6.13'] et ['equ6.60'].

        Note:
            Cette méthode met à jour automatiquement la synthèse des taux
            de travail via _add_synthese_taux_travail.
        """
        self.taux_tau_rd = {}
        tau_d = self.tau_rd
        f_v_d = self.f_type_rd
        K_v  = self.K_v

        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_6_13 = tau_d / f_v_d
            taux_6_60 = tau_d / (K_v * f_v_d)
            return taux_6_13, taux_6_60

        value = val()
        self.taux_tau_rd['equ6.13'] = value[1][0]
        self.taux_tau_rd['equ6.60'] = value[1][1]
        synthese = [
            ["Cisaillement bois", max(value[1][0], value[1][1]), None],
        ]
        self._add_synthese_taux_travail(synthese)
        return value
    
    
    def show_Kv(self) -> None:
        """Affiche l'image des caractéristiques d'une entaille au cisaillement
        """
        self._show_element("Kv_def.png")


# ================================ Barre assemblées mécaniquement Annexe B ==================================

# class Poutre_assemblee_meca(Projet):
#     def __init__(self, beam_2:object, l: si.mm, disposition: str=["Latérale", "Dessus / Dessous"], recouvrement: list=[0,0], Ki: list=[0,None,0], entraxe: list=[1, None, 1], psy_2: int|float=0, **kwargs):
#         """Classe définitssant une poutre composée d'au maximum 3 éléments connectés entre eux par liaisons mécanique 
#         suivant la théorie de HEIMESHOFF Annexe B de l'EN 1995
#         Args:
#             beam_2 (object): objet contenant l'élément centrale de la poutre i=2, cette objet doit être issu de la classe élément ou de son héritage
#             l (int | float): longueur de la poutre en mm
#             disposition (str): Disposition des éléments supplémentaires à l'élément 2, "Latérale" ou "Dessus / Dessous"
#             recouvrement (int): Si disposition latérale alors donner un recouvrement en mm. 
#                                 Cela correspond à la distance entre le centre géométrique de la pièce 2 et celui de la pièce i.
#                                 ATTENTION z local est vers le bas quand vous donnez le recouvrement.
                                
#             Ki (list, optional): Rigidité de liaison par connecteur entre les éléments entre i=1/2 et i=2/3, en N/mm. Soit Kser soit Ku en fonction du type de vérification
#                                     S'il n'y a que 2 éléments connectés, laisser l'index correspondant vide (ex: [0, 2000]). Defaults to [0,None,0].
                                    
#             entraxe (list, optional): Entraxe entre connecteur en mm suivant i=1 ou i=3. Defaults to [1, None, 1].
#             type_action (int | float, optional): Psy 2 qui permet de prendre en compte le fluage dans le temps,
#                                         si calcul en instantanée alors 0, 
#                                         si intermédiaire = voir dans data/coeff_psy.csv 
#                                         et enfin temps infini = 1. 
#                                         Defaults to 0.
                                        
#             **kwargs (object): beam_1 et ou beam_3
#         """
#         super().__init__(**kwargs)
#         self.beam = [None , beam_2, None]
#         self.l = l * si.mm
#         self.disposition = disposition
#         self.recouvrement = [recouvrement[0]*si.mm, recouvrement[1]*si.mm]
#         self.entraxe = []
#         self.Ki = []
#         for i in range(3):
#             if entraxe[i]:
#                 self.entraxe.append(entraxe[i] * si.mm)
#             else:
#                 self.entraxe.append(None)
#             if Ki[i]:
#                 self.Ki.append(Ki[i] * si.N / si.mm)
#             else:
#                 self.Ki.append(None)
            
#         for key, value in kwargs.items():
#             match key[0:4]:
#                 case "beam" if key[-1] == "2":
#                     print("L'attribut ne peut pas être nommé beam_2, il est déjà pris par l'élément centrale ! beam_1 ou beam_3 possible !")
#                 case "beam" if key[-1] != "2" and 1<= int(key[-1]) <=3 :
#                     self.beam[int(key[-1])-1] = value
#         for index, beam in enumerate(self.beam):
#             if beam is not None:
#                 beam.Emean_fin(psy_2)

   
#     @property
#     def K_def(self):
#         kdef = (None, 0)
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 K_def_i = beam.K_def
#                 K_def_2 = self.beam[1].K_def
#                 if K_def_i != K_def_2:
#                     @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                     def val():
#                         K_def = 2 * sqrt(K_def_i * K_def_2)
#                         return K_def
#                 else:
#                     @handcalc(override="params", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                     def val():
#                         K_def= K_def_2
#                         return K_def
#                 value = val()
#                 if value[1] > kdef[1]:
#                     kdef = value
#         return kdef
    
        
#     @property
#     def Ki_fin(self):
#         """Renvoie le Ki en fonction des Kdef des pièces assemblées et du psy2"""
#         ki_fin = {}
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 K_i = self.Ki[index]
#                 psy_2 = beam.psy_2
#                 K_def = self.K_def[1]
                
#                 @handcalc(override="long", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     K_i_fin = K_i / (1 + psy_2 * K_def)
#                     return K_i_fin
                
#                 if index == 0:
#                     ki_fin["Ki fin 1-2"] = val()
#                 else:
#                     ki_fin["Ki fin 2-3"] = val()
#         return ki_fin
    
    
#     @property
#     def gamma_i(self):
#         """Renvoie le gamma"""
#         gamma = [0, 1, 0]
#         gamma = {"gamma 1": None, "gamma 2": 1, "gamma 3": None}
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 E_mean_fin = beam.E_mean_fin.value * 10 ** -6
#                 A = beam.aire.value * 10 ** 6
#                 entraxe = self.entraxe[index].value*10**3
#                 lo = self.l.value * 10 ** 3
#                 if index == 0:
#                    K_i_fin = self.Ki_fin["Ki fin 1-2"][1].value * 10**-3
#                 else:
#                    K_i_fin = self.Ki_fin["Ki fin 2-3"][1].value * 10**-3
                
#                 @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     gamma_i = (1 + pi ** 2 * E_mean_fin * A * entraxe / (K_i_fin * lo ** 2)) ** (-1)
#                     return gamma_i
#                 gamma["gamma "+str(index+1)] = val()
#         return gamma

    
#     @property
#     def distance_ai(self):
#         """Renvoie la distance à l'axe neutre de la pièce 2"""
#         denominateur = 0
#         if self.disposition == "Latérale":
#             d1 = self.recouvrement[0]
#             d3 = self.recouvrement[1]
            
#         else:
#             h1 = self.beam[0].h_calcul
#             h2 = self.beam[1].h_calcul
#             h3 = self.beam[2].h_calcul
#             d1 = -(h1 + h2)/2
#             d3 = (h3 + h2)/2
            
#         for index, beam in enumerate(self.beam):
#             if beam is not None:
#                 if index == 1:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)]
#                 else:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)][1]
                    
#                 denominateur = denominateur + (gamma_i * beam.E_mean_fin * beam.aire)
                
#         if self.beam[0] == None or self.beam[2] == None:
#             di = [d1, None, d3]
#             for index, beam in enumerate(self.beam):
#                 if beam is not None and index != 1:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)][1]
#                     numerateur = gamma_i * beam.E_mean_fin * beam.aire * di[index]
#         else:
#             numerateur = (self.gamma_i["gamma 1"][1] * self.beam[0].E_mean_fin * self.beam[0].aire * d1 
#                         - self.gamma_i["gamma 3"][1] * self.beam[2].E_mean_fin  * self.beam[2].aire * d3)
#         a2 =  numerateur / denominateur
#         ai = [None, -a2, None]
#         for index, beam in enumerate(self.beam):
#             if beam is not None and index != 1:
#                 di = [d1, None, d3]
#                 ai[index] = di[index] - a2
#         return ai
       
#     @property
#     def EI_eff (self):
#         """Renvoie la rigidité de la poutre connectée"""
#         EIeff_latex = ""
#         EIeff_value = 0
#         for index, beam in enumerate(self.beam):
#             if beam is not None:
#                 if index == 1:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)]
#                 else:
#                     gamma_i = self.gamma_i["gamma "+str(index+1)][1]
                    
#                 E_mean_fin = beam.E_mean_fin
#                 inertie = beam.inertie[0]
#                 aire = self.beam[index].aire
#                 distance_ai = self.distance_ai[index]
                
#                 @handcalc(override="long", precision=1, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():                    
#                     EI_eff_i = E_mean_fin * inertie + gamma_i * E_mean_fin * aire * distance_ai**2
#                     return EI_eff_i
#                 inter = val()
#                 EIeff_latex += f'EI eff {index+1}: '
#                 EIeff_latex += inter[0]
#                 EIeff_value = EIeff_value + inter[1]
                
#         @handcalc(override="params", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def sum_EI(): 
#             EI_eff_global = EIeff_value # Somme des EI efficace
#             return EI_eff_global
#         sumEI = sum_EI()
#         return (EIeff_latex + sumEI[0], sumEI[1])


#     def tau_2_max(self, Vz:float):
#         """Calcul de la contrainte de cisaillement maximale dans l'élément 2 selon Annexe B.4 de l'EN 1995

#         Args:
#             vz (float): Effort de cisaillement Vz en kN
#         """
        
#         beam_2 = deepcopy(self.beam[1])
#         K_cr = Cisaillement._from_parent_class(beam_2).K_cr
        
#         V_z = Vz * si.kN
#         h = self.beam[1].h_calcul / 2 + self.distance_ai[1]
#         gamma_3 = self.gamma_i["gamma 3"][1]
#         E_mean_fin_3 = self.beam[2].E_mean_fin
#         E_mean_fin_2 = self.beam[1].E_mean_fin
#         aire_3 = self.beam[2].aire
#         a_3 = self.distance_ai[2]
#         b_2 = self.beam[1].b_calcul
#         EI_eff = self.EI_eff[1]
        
#         @handcalc(override="long", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             tau_2 = V_z * (gamma_3 * E_mean_fin_3 * aire_3 * a_3 + 0.5 * E_mean_fin_2 * b_2 * h**2) / (b_2 * K_cr * EI_eff)
#             return tau_2
#         return val()
    

#     def sigma_i(self, My:float, beam: int=[1, 2, 3]):
#         """Contrainte de compression parallèle au fil selon Annexe B.3 de l'EN 1995

#         Args:
#             My (float): Moment selon l'axe y en kN.m
#             beam (int, optional): Beam 1 = Barre 1
#                                   Beam 2 = Barre 2 
#                                   Beam 3 = Barre 3.
#         """
#         Mf_z = My * si.kN*si.m
#         index = beam
#         if self.beam[index-1] is not None:
#             if index == 2:
#                 gamma_i = self.gamma_i["gamma "+str(index)]
#             else:
#                 gamma_i = self.gamma_i["gamma "+str(index)][1]
#             E_mean_fin = self.beam[index-1].E_mean_fin
#             distance_ai = self.distance_ai[index-1]
#             EI_eff = self.EI_eff[1]
            
#             @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 sigma_i = Mf_z * gamma_i * E_mean_fin * distance_ai / EI_eff
                
#                 return sigma_i
#             return val()
            
    
#     def sigma_mi(self, My:float, beam: int=[1, 2, 3]):
#         """Contrainte de flexion selon l'axe y selon Annexe B.3 de l'EN 1995

#         Args:
#             My (float): Moment selon l'axe y en kN.m
#             beam (int, optional): Beam 1 = Barre 1
#                                   Beam 2 = Barre 2 
#                                   Beam 3 = Barre 3.
#         """
#         M_y = My * si.kN*si.m
#         index = beam
#         if self.beam[index-1] is not None:
#             E_mean_fin = self.beam[index-1].E_mean_fin
#             h_i = self.beam[index-1].h_calcul
#             EI_eff = self.EI_eff[1]
            
#             @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 sigma_m_i = 0.5 * M_y * E_mean_fin * h_i / EI_eff
                
#                 return sigma_m_i
#             return val()
            
            
#     def F_i(self, Vz:float, connecteur: int=[1,2]):
#         """Permet de calculer la charge à reprendre par organe d'assemblage selon Annexe B.5 de l'EN 1995

#         Args:
#             Vz (float): Effort de cisaillement Vz en kN
#             connecteur (int, optional): Position des connecteurs:
#                         1 : connecteur entre les planches 1 et 2
#                         2 : connecteur entre les planches 2 et 3. 
#                         Defaults to [1,2].
#         """
#         V_z = Vz * si.kN
#         if connecteur == 1 :
#             connecteur = 0
#         gamma_i = self.gamma_i["gamma "+str(connecteur + 1)][1]
#         E_mean_fin = self.beam[connecteur].E_mean_fin
#         aire = self.beam[connecteur].aire
#         distance_ai = self.distance_ai[connecteur]
#         entraxe = self.entraxe[connecteur]
#         EI_eff = self.EI_eff[1]
        
#         @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             F_i = gamma_i * E_mean_fin * aire * distance_ai * entraxe * V_z / EI_eff
#             return F_i
#         return val()
    

# class Poteau_assemble_meca(Poutre_assemblee_meca):
#     def __init__(self, lo_y: si.mm, lo_z: si.mm, type_appuis: str=Compression.COEF_LF, **kwargs):
#         super(Poutre_assemblee_meca, self).__init__(**kwargs)
#         self.lo_comp = {"y":lo_y * si.mm, "z":lo_z * si.mm}
#         self.pole = []
#         if not isinstance(self.l, si.Physical):
#             self.l = self.l * si.mm
#         for i, beam in enumerate(self.beam):
#             if beam is not None:
#                 compression = Compression._from_parent_class(beam, lo_y=lo_y, lo_z=lo_z, type_appuis=type_appuis)
#                 self.pole.append(compression)
#             else:
#                 self.pole.append(None)
    
#     @property
#     def aire(self):
#         """Détermine la surface total du poteau assemblé mécaniquement
#         """
#         aire = 0
#         for beam in self.beam:
#             if beam is not None:
#                 aire = aire + beam.aire
#         return aire
    
#     @property
#     def Ief(self):
#         """Détermine l'inertie efficace à partir d'un module d'élasticité moyen Emean et 
#         de la rigidité équivalente en flexion EIef d'une poutre assemblée mécaniquement.
#         """
#         E_mean_beams = []
#         for pole in self.pole:
#             if pole is not None:
#                 E_mean_beams.append(pole.E_mean_fin)
#         E_mean = np.mean(E_mean_beams)
#         EI_ef = self.EI_eff[1]

#         @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def inertie_ef():
#             I_ef = EI_ef / E_mean # equC.4
#             return I_ef
#         return inertie_ef()
    
#     @property
#     def lamb_ef(self):
#         """ Retourne l'élancement d'un poteau assemblé mécaniquement en compression avec risque de flambement suivant l'axe de rotation z donc une direction de flèche suivant y"""
#         lo_z = self.lo_comp['z'].value * 10**3
#         I_z_ef = self.Ief[1]
#         A_tot = self.aire

#         @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             lamb_z_ef = lo_z * sqrt(A_tot / I_z_ef) # equC.3
#             return lamb_z_ef
#         return val()
    

#     @property
#     def lamb(self):
#         """ Retourne l'élancement d'un poteau en compression avec risque de flambement suivant son axe de rotation """
#         dict_lamb = {"lamb,ef,z": self.lamb_ef}
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 lo_y = self.lo_comp['y'].value * 10**3
#                 coef_lef = pole.coef_lef
#                 I_y = pole.inertie[0].value * 10**12
#                 A = pole.aire.value * 10**6

#                 @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     lamb_y = (lo_y * coef_lef) / sqrt(I_y / A)
#                     return lamb_y
#                 dict_lamb["Pole"+str(i+1) + " lamb,y"] = val()
#         return dict_lamb
    

#     @property
#     def lamb_rel_Axe(self):
#         """ Retourne l'élancement relatif d'un poteau en compression avec risque de flambement suivant son axe de rotation """
#         lamb_ef_z = self.lamb['lamb,ef,z'][1]
#         E_0_05_beams, f_c0k_beams = [], []
#         for pole in self.pole:
#             if pole is not None:
#                 f_c0k_beams.append(float(pole.caract_meca.loc['fc0k'])* si.MPa)
#                 E_0_05_beams.append(float(pole.caract_meca.loc['E005'])* si.MPa)
#         f_c0k_mean = np.mean(f_c0k_beams)
#         E_0_05_mean = np.mean(E_0_05_beams)

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def lamb_rel_z():
#             lamb_rel_ef_z = (lamb_ef_z / pi) * sqrt(f_c0k_mean / E_0_05_mean)
#             return lamb_rel_ef_z
#         dict_lamb_rel = {"lamb,rel,ef,z": lamb_rel_z()}

#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 f_c0k = float(pole.caract_meca.loc['fc0k']) * si.MPa
#                 E_0_05 = int(pole.caract_meca.loc['E005']) * si.MPa
#                 lamb_y = self.lamb["Pole"+str(i+1) + " lamb,y"][1]

#                 @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#                 def val():
#                     lamb_rel_y = (lamb_y / pi) * sqrt(f_c0k / E_0_05)
#                     return lamb_rel_y
#                 dict_lamb_rel["Pole"+str(i+1) + " lamb,rel,y"] = val()
#         return dict_lamb_rel
    

#     @property
#     def k_Axe(self):
#         """ Retourne le facteur Ky ou Kz (fonction de l'axe de flambement) """
#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def ky(beta_C:float, lamb_rel_y:float):
#             k_y = 0.5 * (1 + beta_C * (lamb_rel_y - 0.3) + lamb_rel_y ** 2)
#             return k_y
        
#         dict_k_axe = {}
#         beta_C_poles = 0
#         lamb_rel_ef_z = self.lamb_rel_Axe["lamb,rel,ef,z"][1]
        
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 beta_C = pole.beta_C
#                 if beta_C > beta_C_poles:
#                     beta_C_poles = beta_C
#                 lamb_rel_y = self.lamb_rel_Axe["Pole"+str(i+1) + " lamb,rel,y"][1]
#                 dict_k_axe["Pole"+str(i+1) + " ky"] =  ky(beta_C, lamb_rel_y)

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def kz():
#             k_ef_y = 0.5 * (1 + beta_C_poles * (lamb_rel_ef_z - 0.3) + lamb_rel_ef_z ** 2)
#             return k_ef_y
#         dict_k_axe["k,ef,z"] =  kz()
#         return dict_k_axe
    

#     @property
#     def kc_Axe(self):
#         """ Retourne le coefficient multiplicateur KcAxe  (axe = y ou z suivant axe de rotation en flambement) de fc,0,d """
#         lamb_rel_ef_z = self.lamb_rel_Axe["lamb,rel,ef,z"][1]
#         k_ef_z = self.k_Axe["k,ef,z"][1]

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def kcz():
#             k_c_ef_z = 1 / (k_ef_z + sqrt(k_ef_z** 2 - lamb_rel_ef_z ** 2))
#             return k_c_ef_z

#         @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def kcy(k_y, lamb_rel_y):
#             k_c_y = 1 / (k_y + sqrt(k_y** 2 - lamb_rel_y ** 2))
#             return k_c_y
        
#         dict_kc_axe = {"kc,ef,z": kcz()}
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 k_y = self.k_Axe["Pole"+str(i+1) + " ky"][1]
#                 lamb_rel_y = self.lamb_rel_Axe["Pole"+str(i+1) + " lamb,rel,y"][1]
#                 dict_kc_axe["Pole"+str(i+1) + " kc,y"] =  kcy(k_y, lamb_rel_y)
#         return dict_kc_axe
    
    
#     def sigma_c_0_d(self, Fc0d: float):
#         """ Retourne la contrainte de compression axial en MPa avec:
#             Fc0d : la charge en kN de compression """
#         self.Fc_0_d = Fc0d * si.kN
#         Fc_0_d = self.Fc_0_d
#         A_tot = self.aire

#         @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#         def val():
#             sigma_c_0_d = Fc_0_d / A_tot # equC.2
#             return sigma_c_0_d
#         value = val()
#         self.sigma_c_0_rd = value[1]
#         return value
    

#     def Vd_organe(self, Fc0d: float):
#         """ Retourne l'effort de cisaillement à prendr en compte sur un organe d'assemblage suivant l'annexe C §C.2.2 avec:
#             Fc0d : la charge en kN de compression """
#         self.Fc_0_d = Fc0d * si.kN
#         Fc_0_d = self.Fc_0_d
#         K_c_eff_z = self.kc_Axe["kc,ef,z"][1]
#         lamb_ef_z = self.lamb_ef[1]

#         if lamb_ef_z < 30:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def vd():
#                 V_d = Fc_0_d / (120 * K_c_eff_z)
#                 return V_d
#         elif 30 <= lamb_ef_z < 60:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def vd():
#                 V_d = Fc_0_d * lamb_ef_z / (3600 * K_c_eff_z)
#                 return V_d
#         else:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def vd():
#                 V_d = Fc_0_d / (60 * K_c_eff_z)
#                 return V_d
#         return vd()


#     def f_c_0_d(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
#         """Retourne le dictionnaire des résistances f,c,0,d de l'élément assemblé mécaniquement en MPa

#         Args:
#             loadtype (str): chargement de plus courte durée sur l'élément.
#             typecombi (str): type de combinaison, fondamentale ou accidentelle.

#         Returns:
#             float: f,c,0,d en MPa
#         """
#         self.dict_fc0d = {}
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 self.dict_fc0d["Pole"+str(i+1) + " fc,0,d"] = pole._f_type_d("fc0k", loadtype, typecombi)
#         self.dict_fc0d["fc,0,ef,d"] = np.mean([value[1] for value in self.dict_fc0d.values()])
#         return self.dict_fc0d
    

#     def taux_c_0_d(self, flexion: object=None):
#         """Retourne les taux de travaux de la compression axial.
#         Si l'élement est un poteau (donc avec un travail principalement en compression) et de la flexion combinée (EN 1995-1-1 §6.3.2), 
#         il est possible d'ajouter l'objet flexion et de vérifier cette combinaison.

#         Args:
#             flexion (object, optional): L'objet Flexion avec ces taux de travaux préalablement calculés. Default to None.

#         Returns:
#             list: retourne la liste des taux de travaux en %
#         """
#         self.taux_c_0_rd = {}
#         sigma_c_0_d = self.sigma_c_0_rd

#         f_c_0_ef_d = self.dict_fc0d["fc,0,ef,d"]
#         f_c0_d_tot = sum([value[1] for key, value in self.dict_fc0d.items() if key != "fc,0,ef,d"])
#         f_c0d_tot_Kcy = 0
#         for i, pole in enumerate(self.pole):
#             if pole is not None:
#                 kcy = self.kc_Axe["Pole"+str(i+1) + " kc,y"][1]
#                 fcod_y = self.dict_fc0d["Pole"+str(i+1) + " fc,0,d"][1]
#                 f_c0d_tot_Kcy = f_c0d_tot_Kcy + kcy * fcod_y
        
#         lamb_rel_y  = max([value[1] for key, value in self.lamb_rel_Axe.items() if key != "lamb,rel,ef,z"])
#         lamb_rel_z  = self.lamb_rel_Axe["lamb,rel,ef,z"][1]

#         K_c_eff_z = self.kc_Axe["kc,ef,z"][1]

#         if flexion and isinstance(flexion, Flexion):
#             taux_6_11 = flexion.taux_m_rd["equ6.11"]
#             taux_6_12 = flexion.taux_m_rd["equ6.12"]
#         else:
#             taux_6_11 = 0
#             taux_6_12 = 0
        
#         if lamb_rel_y < 0.3 and lamb_rel_z < 0.3:
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 taux_6_2 = sigma_c_0_d / f_c0_d_tot # equ6.2
#                 taux_6_19 = (sigma_c_0_d / (f_c0_d_tot))**2 + taux_6_11 # equ6.19
#                 taux_6_20 = sigma_c_0_d / (f_c0_d_tot)**2 + taux_6_12 # equ6.20
#                 return taux_6_2, taux_6_19, taux_6_20
#             value = val()
#             self.taux_c_0_rd['equ6.19'] = value[1][3]
#             self.taux_c_0_rd['equ6.20'] = value[1][4]
#         else:      
#             @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
#             def val():
#                 taux_6_2 = sigma_c_0_d / f_c0_d_tot # equ6.2
#                 taux_6_23 = (sigma_c_0_d / (f_c0d_tot_Kcy)) + taux_6_11 # equ6.23
#                 taux_6_24 = sigma_c_0_d / (f_c_0_ef_d * K_c_eff_z) + taux_6_12 # equC.1
#                 return taux_6_2, taux_6_23, taux_6_24
#             value = val()
#             self.taux_c_0_rd['equ6.23'] = value[1][1]
#             self.taux_c_0_rd['equ6.24'] = value[1][2]

#         self.taux_c_0_rd['equ6.2'] = value[1][0]
#         return value



# if __name__=='__main__':
#     beam2 = Barre(60,200,"Rectangulaire", classe="C24", cs=1)
#     beam3 = Barre(60,100,"Rectangulaire", classe="C24", cs=1)
#     beam_ass = Poutre_assemblee_meca(beam_2=beam2, l=5000, disposition="Latérale", recouvrement=[0,120], Kser=[None,None,700], entraxe=[None,None,250], psy_2=0, beam_3=beam3)
#     pole_ass = Poteau_assemble_meca._from_parent_class(beam_ass, lo_y=5000, lo_z=5000, type_appuis="Rotule - Rotule")
#     # print(pole_ass.__dict__)
#     # print(pole_ass.lamb)
#     print(pole_ass.kc_Axe)
