# coding in UTF-8 
import os
import sys
from PIL import Image
from math import *
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from PySide6.QtWidgets import QApplication, QInputDialog
from PySide6.QtCore import Qt

# Calculs structurels
import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

# Import local
from ourocode.eurocode.EC3_Element_droit import Plat
from ourocode.eurocode.EC5_Element_droit import Barre

#======================================================= Tige =========================================================

class Tige(Plat):
    QUALITE_ACIER = tuple(str(key) for key in Plat._data_from_csv(Plat, "qualite_acier.csv").index)

    def __init__(
        self, 
        d:si.mm, 
        d0:si.mm, 
        trou_oblong: bool=("False", "True"),
        qualite: str=QUALITE_ACIER, 
        tete_fraisee: bool=("False", "True"),
        prof_tete_fraisee: si.mm=0,
        verif_filetage: bool=("False", "True"), 
        filetage_EN1090: bool=("True", "False"),
        **kwargs
        ):
        """Configure un objet Tige permettant les vérification suivant l'EN 1993-1-8. 
        Cette classe est hérité de la classe Plat du fichier EC3_Element_droit.py.

        Args:
            d (int): le diamètre de la tige en mm
            d0 (int): le diamètre de percage en mm
            trou_oblong (bool): définit si le trou est oblong, si c'est le cas alors True. Defaults to False.
            qualite (str): classe d'acier de la tige (ex: "4.8")
            tete_fraisee (bool): définit si la tige à tete fraisee, si c'est le cas alors True. Defaults to False.
            prof_tete_fraisee (int): profondeur de la tête fraisée dans la plaque en mm. Si pas de tête fraisée alors laisser 0. Defaults to 0.
            verif_filetage (bool): définit si le filetage du boulon doit être vérifier, si c'est le cas alors True. Defaults to False.
            filetage_EN1090 (bool): définit si le filetage est conforme à l'EN 1090, soit matricé. Si filetage usiné alors False. Defaults to True.
        """
        super().__init__(**kwargs)
        self.d = d * si.mm
        self.d0 = d0 * si.mm
        self.trou_oblong = trou_oblong
        self.qualite = qualite
        self.tete_fraisee = tete_fraisee
        self.prof_tete_fraisee = prof_tete_fraisee * si.mm
        self.verif_filetage = verif_filetage
        self._verif_trou_surdimensionne()

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
            raise KeyError("Le diamètre ne correspond pas à celui d'un boulon standard, vérifier As = 0, An = aire de d")
        
        
    def _verif_trou_surdimensionne(self):
        self.percage_surdimensionne = False
        if self.d.value*10**3 in [12, 14]:
            if self.d0.value*10**3 >= self.d.value*10**3+3:
                self.percage_surdimensionne = True
                raise ValueError((f"Le percage est surdimensionné, pour revenir à un perçage normal, "
                               f"réduire le diamètre de perçage à {self.d.value*10**3 + 1} mm "
                               "selon CNC2M §2.1(1) et NF EN 1090-2.\nIl est toutefois possible de réaliser un trou d+2mm si Fb,Rd <= Fv,Rd"))
        elif self.d.value*10**3 < 27:
            if self.d0.value*10**3 >= self.d.value*10**3+4:
                self.percage_surdimensionne = True
                raise ValueError((f"Le percage est surdimensionné, pour revenir à un perçage normal, "
                               f"réduire le diamètre de perçage à {self.d.value*10**3 + 2} mm "
                               "selon CNC2M §2.1(1) et NF EN 1090-2"))
        elif self.d.value*10**3 >= 27:
            if self.d0.value*10**3 >= self.d.value*10**3+8:
                self.percage_surdimensionne = True
                raise ValueError((f"Le percage est surdimensionné, pour revenir à un perçage normal, "
                               f"réduire le diamètre de perçage à {self.d.value*10**3 + 3} mm "
                               "selon CNC2M §2.1(1) et NF EN 1090-2"))  
        return self.percage_surdimensionne
            

    def pince_metal_boulon(self, trous_oblongs: bool=("False", "True"), corrosion: bool=("False", "True")):
        """ Retourne les pinces du métal minimum dans un assemblage constitué de tige (EN 1993-1-8 §3.5) avec:
            trous_oblongs : si les trous oblongs True sinon False
            corrosion : assemblage exposé à la corrosion = True sinon False 
            en_10025_5 : structure réalisées en acier conformes à l'EN 10025-5, True si vrai sinon False

            NE PRENDS PAS EN COMPTE P1,0 ou P1,i ou P2 diminué quand les boulons sont en quinconce (voir §3.5 -5)
        """
        en_10025_5 = False
        if self._Plat__classe_acier.loc["norme"] == "EN 10025-5": 
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
        k_filetage_EN1090 = self.filetage_EN1090
        A_s = self.As
        A_n = self.An
        f_ub = self.fub
        gamma_M2 = self.GAMMA_M["gamma_M2"]

        # On vérifie la clause §3.6.1(5) si les boulons M12 et M14 sont à d+2mm
        if self.d.value*10**3 in [12, 14] and self.qualite in ["4.8", "5.8", "6.8", "8.8", "10.9"] and self.d0.value*10**3 == self.d.value*10**3+2:
            k_d0 = 0.85
        else:
            k_d0 = 1
            
        if self.qualite in ["4.6", "5.6", "8.8"]:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def filetage():
                alpha_v = 0.6
                F_v_Rd_filete = k_d0 * alpha_v * k_filetage_EN1090 * (f_ub * A_s) / gamma_M2
                return F_v_Rd_filete
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def filetage():
                alpha_v = 0.5
                F_v_Rd_filete = k_d0 * alpha_v * k_filetage_EN1090 * (f_ub * A_s) / gamma_M2
                return F_v_Rd_filete

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def lisse():
            alpha_v = 0.6
            F_v_Rd_lisse = k_d0 * alpha_v * (f_ub * A_n) / gamma_M2
            return F_v_Rd_lisse
        return {"filetage": filetage(), "lisse": lisse()}
    
    
    @property
    def FtRd(self) -> float:
        """Retourne la résistance en traction de la tige en N
        """
        f_ub = self.fub
        A_s = self.As
        gamma_M2 =  self.GAMMA_M["gamma_M2"]
        k_filetage_EN1090 = self.filetage_EN1090
        if self.tete_fraisee:
            k_2 = 0.63
        else:
            k_2 = 0.9
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_t_Rd = k_filetage_EN1090 * (k_2 * f_ub * A_s) /  gamma_M2
            return F_t_Rd
        return val()

    
    #Combinaison des efforts
    def taux_FvEd_FtEd(self, FvEd: si.kN=0, FtEd: si.kN=0):
        """Retourne les taux de travail en cisaillement, en traction et combiné d'une tige

        Args:
            FvEd (int): effort à reprendre en cisaillment en kN
            FtEd (int): effort de traction à reprendre en kN 

        Returns:
            dict: dictionnaire des taux de travail slon tab. 3.4 de l'EN 1993-1-8
        """
        self.FvEd = FvEd * si.kN
        F_v_Ed = self.FvEd
        self.FtEd = FtEd * si.kN
        F_t_Ed = self.FtEd
        self.taux_bl = {}
        F_v_Rd_lisse = self.FvRd["lisse"][1]
        F_v_Rd_filetage = self.FvRd["filetage"][1]
        F_t_Rd = self.FtRd[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_t_d = F_t_Ed / F_t_Rd
            taux_v_d_lisse = F_v_Ed / F_v_Rd_lisse
            return taux_t_d, taux_v_d_lisse

        value1 = val()
        self.taux_bl["taux_t_d"] = value1[1][0]
        self.taux_bl["taux_v_d_lisse"] = value1[1][1]
        synthese = [
            ["Traction de la tige", self.taux_bl["taux_t_d"], None],
            ["Cisaillement sur la partie lisse de la tige", self.taux_bl["taux_v_d_lisse"], None]
        ]
        self._add_synthese_taux_travail(synthese)

        if self.verif_filetage:
            if FvEd and FtEd:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val2():
                    taux_v_d_filetage = F_v_Ed / F_v_Rd_filetage
                    taux_v_t_d = F_v_Ed / min(F_v_Rd_lisse, F_v_Rd_filetage) + F_t_Ed / (1.4 * F_t_Rd)
                    return taux_v_d_filetage, taux_v_t_d
                value2 = val2()
                self.taux_bl["taux_v_d_filetage"] = value2[1][0]
                self.taux_bl["taux_v_t_d"] = value2[1][1]
                synthese2 = [
                    ["Cisaillement au droit du filetage de la tige", self.taux_bl["taux_v_d_filetage"], None],
                    ["Cisaillement et traction combiné sur la tige", self.taux_bl["taux_v_t_d"], None]
                ]
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val2():
                    taux_v_d_filetage = F_v_Ed / F_v_Rd_filetage
                    return taux_v_d_filetage
                value2 = val2()
                self.taux_bl["taux_v_d_filetage"] = value2[1]
                synthese2 = [
                    ["Cisaillement au droit du filetage de la tige", self.taux_bl["taux_v_d_filetage"], None],
                ]

        elif not self.verif_filetage and FvEd and FtEd:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                taux_v_t_d = F_v_Ed / F_v_Rd_lisse + F_t_Ed / (1.4 * F_t_Rd)
                return taux_v_t_d
            value2 = val2()
            self.taux_bl["taux_v_t_d"] = value2[1]
            synthese2 = [
                ["Cisaillement et traction combiné sur la tige", self.taux_bl["taux_v_t_d"], None],
            ]

        self._add_synthese_taux_travail(synthese2)
        latex = value1[0] + value2[0]
        return (latex, self.taux_bl)
    
    
    def Bp_Rd(self, d_ecrou: si.mm, d_head_bl: si.mm):
        """Retourne la résistance au poinçonnement de la plaque en N

        Args:
            d_ecrou (int): diamètre extérieur de l'écrou en mm
            d_head_bl (int): diamètre de la tête de boulon en mm

        Returns:
            float: résistance de calcul en N
        """
        d_ecrou = d_ecrou * si.mm
        d_head_bl = d_head_bl * si.mm
        d_m = (d_ecrou + d_head_bl) / 2
        t_p = self.t
        gamma_M2 = self.GAMMA_M["gamma_M2"]
        f_u = self.fu
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            B_p_Rd = (0.6 * pi * d_m * t_p * f_u) / gamma_M2
            return B_p_Rd
        return val()

    def _Kb(self):
        """Détermine le facteur kb selon CNC2M §2.1 tab2"""
        kb = 1
        if self.trou_oblong:
            kb = 0.6
        elif self.percage_surdimensionne:
            kb = 0.8
        return kb
    
    def Fb_Rd(self, e1: si.mm ,e2: si.mm , p1: si.mm, p2: si.mm):
        """Retourne la pression diamétrale en N. 
           ATTENTION: ne prends pas en compte les réductions de résistance lié au critère de l'assemblage (voir §3.6.1(10))
                - Assemblage à une seule rangée de boulon en simple cisaillement -> rondelles + limitation de Fb,Rd <= 1.5*fu*d*t/gamma_M2

        Args:
            e1 (float): pince e1 en mm parallèlle à l'effort si alpha = 0°
            e2 (float): pince e2 en mm perpendiculaire à l'effort si alpha = 0°
            p1 (float): pince p1 en mm parallèlle à l'effort si alpha = 0°
            p2 (float): pince p2 en mm perpendiculaire à l'effort si alpha = 0°
            alpha (float): angle de l'effort en degrés

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
        gamma_M_2 = self.GAMMA_M["gamma_M2"]
        k_b = self._Kb()
        if self.tete_fraisee:
            e_tete_fraisee = self.prof_tete_fraisee
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_eff = t - e_tete_fraisee / 2
                alpha_bx = min(e_1 / (3 * d_0), (p_1 / (3 * d_0)) - 0.25, f_u_b / f_u, 1)
                k_1x = min(2.8 * e_2 / d_0 - 1.7, 1.4 * p_2 / d_0 - 1.7, 2.5)
                F_bx_Rd = k_b * k_1x * alpha_bx * f_u_b * d * t_eff / gamma_M_2

                alpha_by = min(e_2 / (3 * d_0), (p_2 / (3 * d_0)) - 0.25, f_u_b / f_u, 1)
                k_1y = min(2.8 * e_1 / d_0 - 1.7, 1.4 * p_1 / d_0 - 1.7, 2.5)
                F_by_Rd = k_b * k_1y * alpha_by * f_u_b * d * t_eff / gamma_M_2
                return {"F_bx_Rd": F_bx_Rd, "F_by_Rd": F_by_Rd}
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                alpha_bx = min(e_1 / (3 * d_0), (p_1 / (3 * d_0)) - 0.25, f_u_b / f_u, 1)
                k_1x = min(2.8 * e_2 / d_0 - 1.7, 1.4 * p_2 / d_0 - 1.7, 2.5)
                F_bx_Rd = k_b * k_1x * alpha_bx * f_u_b * d * t / gamma_M_2

                alpha_by = min(e_2 / (3 * d_0), (p_2 / (3 * d_0)) - 0.25, f_u_b / f_u, 1)
                k_1y = min(2.8 * e_1 / d_0 - 1.7, 1.4 * p_1 / d_0 - 1.7, 2.5)
                F_by_Rd = k_b * k_1y * alpha_by * f_u_b * d * t / gamma_M_2
                return {"F_bx_Rd": F_bx_Rd, "F_by_Rd": F_by_Rd}
        return val()
    

    def taux_Fb_d(self, FvEd: si.kN=0, alpha: float=0, FbxRd: si.kN=0, FbyRd: si.kN=0):
        """Retourne le taux de travail de la pression diamétrale avec un effort combiné ou non selon l'EN 1993 §3.6.1 tab3.4 et le CNC2M-N0175-REC §2.1(5).    

        Args:
            FvEd (float, optional): effort de cisaillement sur la tige en kN.
            alpha (float, optional): angle de l'effort en degrés suivant la direction parallèle à la pince e1.
            FbxRd (float, optional): résistance de calcul de la pression diamétrale suivant la direction parallèle à la pince e1 en kN.
            FbyRd (float, optional): résistance de calcul de la pression diamétrale suivant la direction parallèle à la pince e2 en kN.

        Returns:
            float: taux de travail
        """
        F_v_Ed = FvEd * si.kN
        F_bx_Rd = FbxRd * si.kN
        F_by_Rd = FbyRd * si.kN
        
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_vx_Ed = F_v_Ed * cos(radians(alpha))
            F_vy_Ed = F_v_Ed * sin(radians(alpha))
            taux_Fbx_d = F_vx_Ed / F_bx_Rd
            taux_Fby_d = F_vy_Ed / F_by_Rd
            taux_Fb_d = (F_vx_Ed / F_bx_Rd)**2 + (F_vy_Ed / F_by_Rd)**2
            return taux_Fbx_d, taux_Fby_d, taux_Fb_d

        result = val()
        synthese = [
            ["Pression diamétrale dans la direction x", result[1][0], None],
            ["Pression diamétrale dans la direction y", result[1][1], None],
            ["Pression diamétrale combinée", result[1][2], None]
        ]
        self._add_synthese_taux_travail(synthese)
        return result

    def Veff_Rd(self, Lnt: si.mm, Lvt: si.mm, effort: str=("Centré", "Excentré")):
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
        gamma_M_2 = self.GAMMA_M["gamma_M2"]

        if effort == "Excentré":
            k_ex = 0.5

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            V_eff_nt_Rd = (A_nt * f_u) / gamma_M_2
            V_eff_nv_Rd = 1/sqrt(3) * (A_vt * f_y) / gamma_M_2
            V_eff_Rd = k_ex * V_eff_nt_Rd + V_eff_nv_Rd
            return V_eff_Rd
        return val()
    
    
    def taux_Veff_d(self, N_Ed: si.kN=0, N_Veff_Rd: si.kN=0, V_Ed: si.kN=0, V_Veff_Rd: si.kN=0):
        """Retourne le taux de travail du cisaillement bloc avec un effort combiné ou non selon l'EN 1993 §3.10.2 et le CNC2M-N0175-REC §2.1(10).    

        Args:
            N_Ed (float, optional): effort normal en kN.
            N_Veff_Rd (float, optional): résistance de l'effort de traction en rupture bloc en kN, si pas d'effort ne rien mettre. Defaults to 0.
            V_Ed (float, optional): effort de cisaillement en kN.
            V_Veff_Rd (float, optional): résistance de l'effort de cisaillement en kN ,si pas d'effort ne rien mettre. Defaults to 0.

        """
        N_Ed = N_Ed * si.kN
        N_Veff_Rd = N_Veff_Rd * si.kN
        V_Ed = V_Ed * si.kN
        V_Veff_Rd = V_Veff_Rd * si.kN

        # selon CNC2M §2.1(10)
        if N_Ed.value and V_Ed.value:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux = N_Ed / N_Veff_Rd + V_Ed / V_Veff_Rd
                return taux
        elif N_Ed.value:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_N_eff = N_Ed / N_Veff_Rd
                return taux_N_eff
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_V_eff = V_Ed / V_Veff_Rd
                return taux_V_eff
        result = val()
        synthese = [
            ["Cisaillement bloc du plat", result[1], None]
        ]
        self._add_synthese_taux_travail(synthese)
        return result

    def axe_articulation(self, FvEd:si.kN=0, jeu:si.mm=0, t_flasque:si.mm=0, t_ame:si.mm=0):
        """Retourne le taux de travail en double cisaillement et en flexion d'un axe d'articulation selon EN 1993-1-8 §3.13

        Args:
            FvEd (si.kN): effort de cisaillement sur l'axe en kN provenant du plat en âme.
            jeu (si.mm): jeu entre plat en mm.
            t_flasque (si.mm): épaisseur de la flasque en mm.
            t_ame (si.mm): épaisseur de l'ame en mm.

        Returns:
            tuple: (taux_v_d, taux_m_d, taux_m_v_d, taux_f_b_flasque_d, taux_f_b_ame_d)
        """
        self.taux_axe_articulation = {}
        F_v_Ed = FvEd * si.kN
        jeu = jeu * si.mm
        t_flasque = t_flasque * si.mm
        t_ame = t_ame * si.mm
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        gamma_M2 = self.GAMMA_M["gamma_M2"]
        f_y = self.fy
        f_y_b = self.fyb
        f_u_b = self.fub
        d = self.d
        if self.verif_filetage:
            A = self.As
        else:
            A = self.An

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_v_Rd = 0.6 * A * f_u_b / gamma_M2
            taux_v_d = F_v_Ed / F_v_Rd

            M_Ed = F_v_Ed * (t_ame + 4 * jeu + 2 * t_flasque)
            M_Rd = 1.5 * pi * d**3 * f_y_b / (32 * gamma_M0)
            taux_m_d = M_Ed / M_Rd
            taux_m_v_d = (M_Ed / M_Rd)**2 + (F_v_Ed / F_v_Rd)**2

            F_b_flasque_Rd = 1.5 * t_flasque * d * f_y / gamma_M0
            taux_f_b_flasque_d = 0.5 * F_v_Ed / F_b_flasque_Rd

            F_b_ame_Rd = 1.5 * t_ame * d * f_y / gamma_M0
            taux_f_b_ame_d = F_v_Ed / F_b_ame_Rd

            return taux_v_d, taux_m_d, taux_m_v_d, taux_f_b_flasque_d, taux_f_b_ame_d
        result = val()
        self.taux_axe_articulation["taux_v_d"] = result[1][0]
        self.taux_axe_articulation["taux_m_d"] = result[1][1]
        self.taux_axe_articulation["taux_m_v_d"] = result[1][2]
        self.taux_axe_articulation["taux_f_b_flasque_d"] = result[1][3]
        self.taux_axe_articulation["taux_f_b_ame_d"] = result[1][4]

        synthese = [
            ["Cisaillement de l'axe d'articulation", self.taux_axe_articulation["taux_v_d"], None],
            ["Flexion de l'axe d'articulation", self.taux_axe_articulation["taux_m_d"], None],
            ["Cisaillement et flexion combinée de l'axe d'articulation", self.taux_axe_articulation["taux_m_v_d"], None],
            ["Pression diamétrale des flasques de l'axe d'articulation", self.taux_axe_articulation["taux_f_b_flasque_d"], None],
            ["Pression diamétrale de l'ame de l'axe d'articulation", self.taux_axe_articulation["taux_f_b_ame_d"], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result


class _Tige_scellement(Tige):
    GAMMA_C = {"Fondamentale": 1.5, "Accidentelle": 1.2}
    CLASSE_CONCRETE = list(Tige._data_from_csv(Tige, "caracteristique_meca_beton.csv").index)[2:]
    def __init__(self, n:int=2, n_traction:int=2, classe_beton: str=CLASSE_CONCRETE, **kwargs):
        """Initialise une tige de prescellement avec plaque circulaire d'après l'EN1993-1-8 et le guide du CNC2M.
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py.

        Args:
            n (int, optional): nombre de tiges d'ancrage. Defaults to 2.
            n_traction (int, optional): nombre de tiges d'ancrage en traction dans le cas
             d'une reprise de moment en pieds par exemple sinon égale à n. Defaults to 2.
            classe_beton (str): classe du béton.
        """
        kwargs["filetage_EN1090"] = False
        super().__init__(**kwargs)
        self._test_qualite_acier()
        self.n = n
        self.n_traction = n_traction
        self.classe_beton = classe_beton
        caract_meca = self._data_from_csv("caracteristique_meca_beton.csv").loc[self.classe_beton]
        self.fck = float(caract_meca.loc["fck"]) * si.MPa

    def _test_qualite_acier(self):
        acier_valide = ["4.6", "4.8", "5.6", "5.8", "6.8", "8.8", "9.8", "10.9"]
        if self.qualite not in acier_valide:
            raise ValueError("Veuillez entrer une qualité d'acier pour l'ancrage valide {}".format(acier_valide))

    @property
    def FvRd(self) -> float:
        """Retourne la résistance en cisaillement de la partie fileté et lisse d'une tige par plan en N dans le cas 
        d'une reprise du cisaillement par les tiges d'ancrage selon l'EN1993-1-8 §6.2.2(7)
        """
        k_filetage_EN1090 = self.filetage_EN1090
        A_s = self.As
        f_yb = self.fyb.value * 10**-6
        f_ub = self.fub
        gamma_M2 = self.GAMMA_M["gamma_M2"]
        if 235*si.MPa <= self.fyb <= 640*si.MPa:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                alpha_bc = 0.44 - 0.0003 * f_yb
                F_v_Rd_filete = alpha_bc * k_filetage_EN1090 * f_ub * A_s / gamma_M2
                return F_v_Rd_filete
            return val()
        else:
            raise ValueError("La limite elastique fyb doit être comprise entre 235 MPa et 640 MPa pour être considéré en cisaillement")

    def taux_ancrage(self, FtEd: si.kN=0, FvEd: si.kN=0):
        """Calcul le taux de travail d'une tige d'ancrage d'après l'EN1993-1-8 et le guide du CNC2M.

        Args:
            FtEd (float): effort de traction sur la platine de scellement en kN.
            FvEd (float): effort de cisaillement sur la platine de scellement en kN. 
                Il est en général déconseillé de faire passer de grand effort de cisaillement par les tiges, privilégier une bèche.

        Returns:
            tuple: le tuple des taux de travail
        """
        F_t_Ed = FtEd * si.kN
        F_v_Ed = FvEd * si.kN
        latex_F_tc_Rd, F_tc_Rd = self.FtcRd
        latex_F_tb_Rd, F_tb_Rd = self.FtRd
        n = self.n
        n_traction = self.n_traction
        if not F_v_Ed or F_v_Ed.value == 0:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_t_ancr_Rd = min(F_tc_Rd, F_tb_Rd)
                taux_ancrage = F_t_Ed / (F_t_ancr_Rd * n_traction) #CNC2M §6(8)
                return taux_ancrage
            result = val()
            synthese = [
                ["Traction de la tige d'ancrage", result[1], None],
            ]
            self._add_synthese_taux_travail(synthese)
            return (latex_F_tb_Rd + latex_F_tc_Rd + result[0], result[1])

        elif not self.percage_surdimensionne:
            latex_F_v_Rd, F_v_Rd = self.FvRd
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_t_ancr_Rd = min(F_tc_Rd, F_tb_Rd)
                taux_t = F_t_Ed / (F_t_ancr_Rd * n_traction) #CNC2M §6(8)
                taux_t_v = F_v_Ed / (F_v_Rd * n) + F_t_Ed / (1.4 * F_t_ancr_Rd * n_traction)
                return taux_t, taux_t_v
            result = val()
            synthese = [
                ["Traction de la tige d'ancrage", result[1][0], None],
                ["Cisaillement et traction combiné de la tige d'ancrage", result[1][1], None],
            ]
            self._add_synthese_taux_travail(synthese)
            return (latex_F_v_Rd + latex_F_tb_Rd + latex_F_tc_Rd + result[0], result[1])

        else:
            t_p = self.t
            d = self.d
            latex_F_v_Rd, F_v_Rd = self.FvRd
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                e = t_p + d / 2
                F_t_Ed_eq = F_v_Ed * (e / d) * (5 * pi)/6 #CNC2M §6(9)
                F_t_ancr_Rd = min(F_tc_Rd, F_tb_Rd)
                taux_t = (F_t_Ed + F_t_Ed_eq) / (F_t_ancr_Rd * n_traction) #CNC2M §6(9)
                taux_t_v = F_v_Ed / (F_v_Rd * n) + (F_t_Ed + F_t_Ed_eq) / (1.4 * F_t_ancr_Rd * n_traction)
                return taux_t, taux_t_v
            result = val()
            synthese = [
                ["Traction de la tige d'ancrage", result[1][0], None],
                ["Cisaillement et traction combiné de la tige d'ancrage", result[1][1], None],
            ]
            self._add_synthese_taux_travail(synthese)
            return (latex_F_v_Rd + latex_F_tb_Rd + latex_F_tc_Rd + result[0], result[1])


class Tige_scellement_plaque_ancrage(_Tige_scellement):
    def __init__(self, tr: si.mm, dr: si.mm, l_ancr: si.mm, d1: si.mm, e_tiges: si.mm, n:int=2, n_traction:int=2, classe_beton: str=_Tige_scellement.CLASSE_CONCRETE, **kwargs):
        """Initialise une tige de scellement avec plaque circulaire d'après l'EN1993-1-8 et le guide du CNC2M.
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py.

        Args:
            tr (si.mm): épaisseur de la plaque d'ancrage en mm.
            dr (si.mm): diamètre de la plaque d'ancrage en mm.
            l_ancr (si.mm): longueur de l'ancrage en mm (distance entre le nu béton et la face de la plaque d'ancrage).
            d1 (si.mm): distance perpendiculaire du bord béton à l'axe de la tige en mm.
            e_tiges (si.mm): entraxe entre les tiges d'ancrage en mm.
            n (int, optional): nombre de tiges d'ancrage. Defaults to 2.
            n_traction (int, optional): nombre de tiges d'ancrage en traction dans le cas
                d'une reprise de moment en pieds par exemple sinon égale à n. Defaults to 2.
            classe_beton (str): classe du béton.
        """
        super().__init__(n=n, n_traction=n_traction, classe_beton=classe_beton,**kwargs)
        self.tr = tr * si.mm
        self.dr = dr * si.mm
        self.l_ancr = l_ancr * si.mm
        self.d1 = d1 * si.mm
        self.e_tiges = e_tiges * si.mm
        self.type_ancrage = "Plaque d'ancrage circulaire"

        if self.tr.value < 0.3*self.dr.value/2:
            raise ValueError("L'épaisseur tr doit être supérieur ou égal à {} mm".format(ceil(0.3*self.dr.value*10**3/2)))

    @property
    def FtcRd(self):
        d = self.d
        d_r = self.dr
        l_ancr = self.l_ancr
        d_1 = self.d1
        e_tiges = self.e_tiges
        gamma_c = self.GAMMA_C["Fondamentale"]
        f_ck = self.fck
        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            r_r = d_r / 2
            v = min(l_ancr, d_1, e_tiges)
            A_ancr = 2.55 * pi * (r_r**2 - d**2/4) * (1-r_r/v)
            F_tc_Rd = A_ancr * f_ck / gamma_c #CNC2M §6(6) tab19
            return F_tc_Rd
        return val()

class Tige_scellement_crochet(_Tige_scellement):
    def __init__(self, l1: si.mm, l2: si.mm, rayon_crochet: si.mm, n:int=2, n_traction:int=2, classe_beton: str=_Tige_scellement.CLASSE_CONCRETE, **kwargs):
        """Initialise une tige de scellement avec crochet d'après l'EN1993-1-8 et le guide du CNC2M.
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py.

        Args:
            l1 (si.mm): longueur d'ancrage l1 selon tab19 du CNC2M en mm (distance entre le nu béton et le début du cintrage).
            l2 (si.mm): longueur d'ancrage l2 selon tab19 du CNC2M en mm (petite distance de l'extrémité du crochet à la fin du cintrage).
            rayon_crochet (si.mm): rayon du crochet en mm.
            n (int, optional): nombre de tiges d'ancrage. Defaults to 2.
            n_traction (int, optional): nombre de tiges d'ancrage en traction dans le cas
                d'une reprise de moment en pieds par exemple sinon égale à n. Defaults to 2.
            classe_beton (str): classe du béton.
        """
        super().__init__(n=n, n_traction=n_traction, classe_beton=classe_beton,**kwargs)
        self.l1 = l1 * si.mm
        self.l2 = l2 * si.mm
        self.rayon_crochet = rayon_crochet * si.mm
        self.type_ancrage = "Crochet"

        if self.rayon_crochet.value*10**3 < 3*self.d.value*10**3:
            raise ValueError("Le rayon du crochet doit être supérieur ou égal à {} mm".format(ceil(3*self.d.value*10**3)))
        if self.l2.value*10**3 < 1.5*self.d.value*10**3 or self.l2.value*10**3 > 2*self.d.value*10**3:
            raise ValueError("La longueur l2 doit être comprise entre {} mm <= l2 <= {} mm".format(ceil(1.5*self.d.value*10**3), floor(2*self.d.value*10**3)))

    @property
    def FtcRd(self):
        d = self.d
        l_1 = self.l1
        l_2 = self.l2
        rayon = self.rayon_crochet
        gamma_c = self.GAMMA_C["Fondamentale"]
        f_c_k = self.fck
        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            A_ancr = pi * d * (l_1 + 6.4 * rayon + 3.5 * l_2)
            f_b_d = 0.36 * sqrt(f_c_k) / gamma_c * si.MPa
            F_tc_Rd = A_ancr * f_b_d #CNC2M §6(6) tab19
            return F_tc_Rd
        return val()

class Soudure(Plat):
    def __init__(self, t2: si.mm, gorge: si.mm, l: si.mm, retour_soudure: bool=("False", "True"), alpha: float=90, **kwargs):
        """Configure un objet soudure permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Plat du fichier EC3_Element_droit.py.  

        Args:
            t2 (float): épaisseur de la pièce 2 sur laquelle on soude en mm
            gorge (int): dimension de la gorge "a" en mm
            l (int): longueur de soudure "brute" sans cratère en mm
            retour_soudure (bool, optional): détermine si un retour de la soudure est fait si oui alors True. Defaults to False.
            alpha (int | float, optional): angle en degré de la de la pièce 2 sur la pièce 1. Defaults to 90.
        """
     
        super().__init__(**kwargs)
        self.t2 = t2 * si.mm
        self.gorge = gorge * si.mm
        self.l = l * si.mm
        self.retour_soudure = retour_soudure
        self.alpha = alpha

        self.verif_soudure()


    @property
    def beta_w(self):
        return float(self._Plat__classe_acier.loc["betaW"])


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
                    raise ValueError(f"La longueur du cordon de soudure est trop petite, elle doit être supérieur à {min(30*si.mm, 6*self.gorge)}")
            else:
                raise ValueError(f"La gorge doit être au minimum de {amin} et au maximum de {amax}")
        else:
            raise ValueError("L'angle entre les deux pièces à souder doit être compris entre 60° et 120°")
        
    
    def beta_Lw1(self, Lj: si.mm) -> float:
        """Calcul le facteur beta_Lw,1 qui dimminue la résistance pour des cordons de soudure des assemblages par recouvrement (à plat)

        Args:
            Lj (int): Longueur de recouvrement des plats en mm
        """
        Lj = Lj * si.mm
        return min((1.2 - 0.2 * Lj) / (150 * self.gorge), 1)

    
    def cordon_frontal(self, N_Ed: si.kN=0):
        """Calcul un cordon de soudure frontale et retourne le taux de travail.

        Args:
            N_Ed (float): Effort de traction en kN.
        """
        N_Ed = N_Ed * si.kN
        gamma_M_2 = self.GAMMA_M["gamma_M2"]
        beta_w = self.beta_w
        a = self.gorge
        l_ef = self.lef
        fu = self.fu
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            cordon_frontal = (beta_w * gamma_M_2 * (N_Ed * sqrt(2)) / fu) / (a * l_ef)
            return cordon_frontal
        result = val()
        synthese = [
            ["Cordon de soudure frontal", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result


    def cordon_lateral(self, V_Ed: si.kN=0):
        """Calcul un cordon de soudure latérale et retourne le taux de travail.

        Args:
            V_Ed (float): Effort de cisaillement du cordon en kN.
        """
        V_Ed = V_Ed * si.kN
        gamma_M_2 = self.GAMMA_M["gamma_M2"]
        beta_w = self.beta_w
        a = self.gorge
        l_ef = self.lef
        fu = self.fu
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            cordon_laterale = (beta_w * gamma_M_2 * (V_Ed * sqrt(3)) / fu) / (a * l_ef)
            return cordon_laterale
        result = val()
        synthese = [
            ["Cordon de soudure latéral", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result


    def cordon_oblique(self, alpha_cordon: float, N_Ed: si.kN=0):
        """Calcul un cordon de soudure oblique et retourne le taux de travail.

        Args:
            N_Ed (float): Effort de traction en kN.
        """
        self.alpha_cordon = alpha_cordon
        N_Ed = N_Ed * si.kN
        gamma_M_2 = self.GAMMA_M["gamma_M2"]
        beta_w = self.beta_w
        a = self.gorge
        l_ef = self.lef
        fu = self.fu
        alpha_cordon = self.alpha_cordon
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            cordon_oblique = (beta_w * gamma_M_2 * (N_Ed * sqrt(3 - sin(radians(alpha_cordon))**2)) / fu) / (a * l_ef)
            return cordon_oblique
        result = val()
        synthese = [
            ["Cordon de soudure oblique", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result


    def cordon_pieces_obliques(self, N_Ed: si.kN=0):
        """Calcul un cordon de soudure sur des pièces à positionnement obliques et retourne le taux de travail.

        Args:
            N_Ed (float): Effort de traction en kN.
        """
        N_Ed = N_Ed * si.kN
        gamma_M_2 = self.GAMMA_M["gamma_M2"]
        beta_w = self.beta_w
        a = self.gorge
        l_ef = self.lef
        fu = self.fu
        alpha = self.alpha
        if alpha < 90:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                cordon_pieces_obliques = (beta_w * gamma_M_2 * (N_Ed * sqrt(2 - sin(radians(alpha)))) / fu) / (a * l_ef)
                return cordon_pieces_obliques
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                cordon_pieces_obliques = (beta_w * gamma_M_2 * (N_Ed * sqrt(2 + sin(radians(alpha)))) / fu) / (a * l_ef)
                return cordon_pieces_obliques
        result = val()
        synthese = [
            ["Cordon de soudure d'une pièce oblique", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result


    def critere_generale(self, FvEd: si.kN=0, FaxEd: si.kN=0):
        """Calcul le critère générale de Von Mises d'une soudure et retourne le taux.

        Args:
            FvEd (int | float): Effort de cisaillement sur la en kN
            FaxEd (int | float): Effort de traction sur la soudure en kN

        Returns:
            (float) : Taux de travail de la soudure
        """
        F_v_Ed = FvEd * si.kN
        F_ax_Ed = FaxEd * si.kN
        f_u = self.fu
        gamma_M2 = self.GAMMA_M["gamma_M2"]
        beta_w = self.beta_w
        a = self.gorge
        l_ef = self.lef
        alpha = self.alpha
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            tau_para = F_v_Ed / (a * l_ef)
            sigma_perpend = (F_ax_Ed * cos(radians(alpha/2)))/ (a * l_ef)
            tau_perpend = (F_ax_Ed * cos(radians(alpha/2)))/ (a * l_ef)
            sigma_e = sqrt(sigma_perpend**2 + 3 * (tau_perpend**2 + tau_para**2))
            sigma_Rd = f_u / (beta_w * gamma_M2)
            sigma_perpend_Rd = (0.9 * f_u) / gamma_M2
            taux_von_mises = max(sigma_e / sigma_Rd, sigma_perpend / sigma_perpend_Rd)
            return taux_von_mises
        result = val()
        synthese = [
            ["Cordon de soudure vérification du critère de Von Mises", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result
    

    def soudure_discontinue(self, b: si.mm, b1: si.mm, corrosion: bool=("False", "True")):
        """Détermine les dimensions des cordons de soudure discontinue

        Args:
            b (int): voir EC3 1-8
            b1 (int): hauteur en mm de la pièce 2 soudé sur la pièce 1
            corrosion (bool, optional): _description_. Defaults to False.

        Returns:
            dict: dimensions des cordons de soudure discontinue
        """
        if corrosion:
            raise ValueError("Il n'est pas possible d'avoir une soudure discontinue en ambiance corrosive")
        b = b * si.mm
        lwe = max(0.75 * b, 0.75 * b1*si.mm)
        l1 = min(16 * self.t, 16 * self.t2, 200)
        l2 = min(12 * self.t, 12 * self.t2, 0.25 * b, 200)
        return {"Lwe": lwe, "L1": l1, "L2": l2}

# class Platine_about(Plat):
#     """
#     Classe pour la visualisation d'un pied de poteau en acier avec profilé H ou I.
#     Permet de représenter graphiquement le profilé et l'emplacement des boulons.
#     """
    
#     def __init__(self, h_profil: float, b_profil: float, tw: float, tf: float, r: float, 
#                  type_profil: str = "HEB", echelle: float = 1.0,
#                  decalage_y: float = 0.0,
#                  **kwargs):
#         """
#         Initialise la visualisation d'un pied de poteau.
        
#         Args:
#             h_profil (float): Hauteur totale du profilé en mm
#             b_profil (float): Largeur des ailes du profilée en mm
#             tw (float): Épaisseur de l'âme du profilée en mm
#             tf (float): Épaisseur des ailes du profilée en mm
#             r (float): Rayon de raccordement âme-aile en mm
#             type_profil (str): Type de profilé ('HE', 'IPE', 'INP', etc.)
#             echelle (float): Échelle de visualisation (utile pour les très grands/petits profils)
#             decalage_y (float, optional): Décalage vertical du profilé par rapport au centre de la platine. Par défaut 0.0.
#         """
#         super().__init__(**kwargs)
#         self.h_profil = h_profil
#         self.b_profil = b_profil
#         self.tw = tw
#         self.tf = tf
#         self.r = r
#         self.type_profil = type_profil
#         self.echelle = echelle
        
#         # Décalage du profilé par rapport au centre de la platine
#         self.decalage_y = decalage_y
        
#         self.boulons = []
#         self.fig, self.ax = plt.subplots(figsize=(12, 10))
#         self.ax.set_aspect('equal')
#         self.ax.set_title(f"Platine d'about - Profile {type_profil}")
#         self.ax.set_xlabel('Largeur (mm)')
#         self.ax.set_ylabel('Hauteur (mm)')
#         self.ax.grid(True, linestyle='--', alpha=0.7)
    
#     def add_tiges(self, positions: list[tuple[float, float]], diametre: float):
#         """
#         Ajoute des tiges à la visualisation.
        
#         Args:
#             positions (list): Liste de tuples [(x, y)] des positions des tiges (en mm) vis à vis du centre de gravité de la platine
#             diametre (float): Diamètre des boulons en mm
#         """
#         self.boulons.extend([(x, y, diametre) for x, y in positions])
    
#     def _dessiner_profil(self):
#         """Dessine le profilé du poteau"""
#         # Position du profilé avec décalage
#         x_center = 0
#         y_center = self.decalage_y
        
#         # Dessin du profilé en H
#         if self.type_profil.upper() in ['HEB', 'HEA', 'IPE', 'INP']:
#             # Aile supérieure
#             aile_sup = Rectangle((x_center - self.b_profil/2, y_center + (self.h_profil/2 - self.tf)), 
#                                 self.b_profil, self.tf, 
#                                 linewidth=1, edgecolor='black', facecolor='lightgray')
            
#             # Aile inférieure
#             aile_inf = Rectangle((x_center - self.b_profil/2, y_center - self.h_profil/2), 
#                                self.b_profil, self.tf, 
#                                linewidth=1, edgecolor='black', facecolor='lightgray')
            
#             # Âme
#             ame = Rectangle((x_center - self.tw/2, y_center - self.h_profil/2 + self.tf), 
#                            self.tw, self.h_profil - 2*self.tf, 
#                            linewidth=1, edgecolor='black', facecolor='lightgray')
            
#             # Ajout des éléments au graphique
#             for patch in [aile_sup, aile_inf, ame]:
#                 self.ax.add_patch(patch)
        
#         # Ajout des cotes
#         self._ajouter_cotes()
    
#     def _dessiner_boulons(self):
#         """Dessine les boulons sur le dessin"""
#         # Ajout des boulons avec le décalage du profilé
#         for x, y, diametre in self.boulons:
#             boulon = Circle((x, y), diametre/2, 
#                            facecolor='yellow', edgecolor='black', linewidth=1, zorder=5)
#             self.ax.add_patch(boulon)
#             # Ajout d'un point au centre pour les petits boulons
#             self.ax.plot(x, y, 'k+', markersize=5, zorder=6)
    
#     def _ajouter_cotes(self):
#         """Ajoute les cotes principales au dessin"""

#           # Cote décalage Y si nécessaire
#         if abs(self.decalage_y) > 0.1:  # Seulement si le décalage est significatif
#             y_platine = -self.h.value*10**3/2  # Bord inférieur de la platine
#             y_debut = y_platine - 30
#             y_fin = -self.h_profil/2 + self.decalage_y
            
#             self.ax.annotate('', 
#                             xy=(10, 0), 
#                             xytext=(10, self.decalage_y),
#                             arrowprops=dict(arrowstyle='<->', color='red'))
#             self.ax.text(15, self.decalage_y, 
#                         f'décalage Y={abs(self.decalage_y):.0f} mm', 
#                         va='center', ha='center', rotation=90, color='red')

#         # Cote hauteur totale
#         self.ax.annotate('', 
#                         xy=(-self.b_profil/2 - 5, -self.h_profil/2 + self.decalage_y), 
#                         xytext=(-self.b_profil/2 - 5, self.h_profil/2 + self.decalage_y),
#                         arrowprops=dict(arrowstyle='<->'))
#         self.ax.text(-self.b_profil/2 - 15, self.decalage_y, f'h prof={self.h_profil} mm', 
#                     va='center', ha='center', rotation=90)
        
#         # Cote largeur aile
#         self.ax.annotate('', 
#                         xy=(-self.b_profil/2, self.h_profil/2 + self.decalage_y + 5), 
#                         xytext=(self.b_profil/2, self.h_profil/2 + self.decalage_y + 5),
#                         arrowprops=dict(arrowstyle='<->'))
#         self.ax.text(0, self.h_profil/2 + self.decalage_y + 15, f'b prof={self.b_profil} mm', 
#                     ha='center')
        
#         # Cote épaisseur aile
#         self.ax.annotate('', 
#                         xy=(self.b_profil/2 + 5, self.h_profil/2 + self.decalage_y), 
#                         xytext=(self.b_profil/2 + 5, self.h_profil/2 - self.tf + self.decalage_y),
#                         arrowprops=dict(arrowstyle='<->'))
#         self.ax.text(self.b_profil/2 + 15, self.h_profil/2 - self.tf/2 + self.decalage_y, 
#                     f'tf={self.tf} mm', 
#                     va='center')
        
#         # Cote épaisseur âme
#         self.ax.annotate('', 
#                         xy=(-self.tf/2, self.decalage_y), 
#                         xytext=(self.tf/2, self.decalage_y),
#                         arrowprops=dict(arrowstyle='<->'))
#         self.ax.text(-20, self.decalage_y + 10, 
#                     f'tw={self.tw} mm', 
#                     va='center', rotation=0)
    
#     def _dessiner_platine(self):
#         """Dessine la platine sous le profilé en vue de face"""
#         # Position verticale de la platine (sous le profilé)
#         y_platine = -self.h/2
        
#         # Création du rectangle de la platine (vue de face)
#         platine = Rectangle(
#             xy=(-self.b/2, y_platine),
#             width=self.b,
#             height=self.h,
#             linewidth=1,
#             edgecolor='black',
#             facecolor='#E0E0E0',  # Gris clair
#             zorder=1,  # Pour que la platine soit en arrière-plan
#             alpha=0.7  # Légère transparence pour voir le profilé derrière
#         )
#         self.ax.add_patch(platine)
#         # Ajout d'un point au centre de la platine
#         self.ax.plot(0, 0, 'k+', markersize=6, zorder=6, color='red')
        
#         # Ajout des cotes de la platine
#         self._ajouter_cotes_platine()
        
#         return platine
    
#     def _ajouter_cotes_platine(self):
#         """Ajoute les cotes spécifiques à la platine"""
#         y_platine = -self.h.value*10**3/2
        
#         # Cote hauteur platine
#         self.ax.annotate(
#             '',
#             xy=(-self.b.value*10**3/2 - 10, y_platine),
#             xytext=(-self.b.value*10**3/2 - 10, y_platine + self.h.value*10**3),
#             arrowprops=dict(arrowstyle='<->', color='blue')
#         )
#         self.ax.text(
#             -self.b.value*10**3/2 - 15, 
#             y_platine + self.h.value*10**3/2, 
#             f'hp={self.h}',
#             va='center',
#             ha='center',
#             rotation=90,
#             color='blue'
#         )
        
#         # Cote largeur platine
#         self.ax.annotate(
#             '',
#             xy=(-self.b.value*10**3/2, y_platine - 10),
#             xytext=(self.b.value*10**3/2, y_platine - 10),
#             arrowprops=dict(arrowstyle='<->', color='blue')
#         )
#         self.ax.text(
#             0, 
#             y_platine - 15, 
#             f'bp={self.b}',
#             ha='center',
#             color='blue'
#         )
    
#     def calculer_distances_boulons(self):
#         """
#         Calcule les distances minimales entre les boulons et les éléments du profilé/platine.
        
#         Returns:
#             list: Liste de dictionnaires contenant les distances pour chaque boulon
#         """
#         if not hasattr(self, 'boulons') or not self.boulons:
#             return []
            
#         # Trier les boulons par coordonnée y décroissante (du haut vers le bas)
#         boulons_tries = sorted(self.boulons, key=lambda x: x[1], reverse=True)
        
#         # Grouper les boulons par rangée (même y)
#         rangees = {}
#         for x, y, diametre in boulons_tries:
#             y_rounded = round(y, 1)  # Arrondi pour gérer les erreurs de précision
#             rangees.setdefault(y_rounded, []).append((x, y, diametre))
        
#         # Calculer les entraxes pour chaque rangée
#         entraxes = {}
#         for y, boulons_rangee in rangees.items():
#             # Trier les boulons de la rangée par x croissant
#             boulons_rangee.sort()
            
#             # Calculer les entraxes pour cette rangée
#             for i in range(len(boulons_rangee) - 1):
#                 x1, y1, d1 = boulons_rangee[i]
#                 x2, y2, d2 = boulons_rangee[i+1]
#                 distance_x = abs(x2 - x1)
#                 distance_y = abs(y2 - y1)
                
#                 # Stocker l'entraxe pour les deux boulons concernés
#                 for x, y, d in [(x1, y1, d1), (x2, y2, d2)]:
#                     if (x, y) not in entraxes:
#                         entraxes[(x, y)] = []
#                     entraxes[(x, y)].append((distance_x, distance_y))
        
#         # Calculer les distances pour chaque boulon
#         distances = []
#         for i, (x, y, diametre) in enumerate(self.boulons, 1):
#             # Distances aux bords de la platine
#             dist_gauche = abs(x - (-self.b.value*10**3/2))
#             dist_droite = abs(self.b.value*10**3/2 - x)
#             dist_bas = abs(y - (-self.h.value*10**3/2))
#             dist_haut = abs(self.h.value*10**3/2 - y)
            
#             # Distances par rapport au profilé (si le boulon est à proximité)
#             dist_semelle = None
#             dist_ame = None
#             y_relatif = y - self.decalage_y  # Position relative par rapport au centre du profilé
            
#             # Vérification si le boulon est dans la zone du profilé
#             if (-self.b_profil/2 <= x <= self.b_profil/2 and 
#                 -self.h_profil/2 <= y_relatif <= self.h_profil/2):
                
#                 # Distance à l'aile supérieure
#                 if y_relatif > 0:
#                     dist_semelle = abs(y_relatif - (self.h_profil/2 - self.tf))
#                 # Distance à l'aile inférieure
#                 else:
#                     dist_semelle = abs(-self.h_profil/2 + self.tf - y_relatif)
                
#                 # Distance à l'âme
#                 dist_ame = min(abs(x - self.tw/2), abs(x + self.tw/2))
            
#             # Récupérer les entraxes pour ce boulon
#             entraxes_boulon = entraxes.get((x, y), [])
            
#             distances.append({
#                 'boulon': i,
#                 'x': x,
#                 'y': y,
#                 'diametre': diametre,
#                 'distance_bord_gauche': dist_gauche,
#                 'distance_bord_droit': dist_droite,
#                 'distance_bord_bas': dist_bas,
#                 'distance_bord_haut': dist_haut,
#                 'distance_semelle': dist_semelle if dist_semelle is not None else "N/A",
#                 'distance_ame': dist_ame if dist_ame is not None else "N/A",
#                 'entraxes': entraxes_boulon if entraxes_boulon else "N/A"
#             })
        
#         return distances

#     def afficher_distances_boulons(self):
#         """Affiche les distances minimales des boulons dans la console"""
#         if not hasattr(self, 'boulons') or not self.boulons:
#             print("Aucun boulon défini.")
#             return
            
#         distances = self.calculer_distances_boulons()
#         print("\n" + "="*70)
#         print("DISTANCES MINIMALES DES BOULONS".center(70))
#         print("="*70)
        
#         for dist in distances:
#             print(f"\nBoulon {dist['boulon']} (x={dist['x']:.1f} mm, y={dist['y']:.1f} mm, Ø={dist['diametre']} mm):")
#             print("-"*60)
#             print(f"  Distance bord gauche : {dist['distance_bord_gauche']:7.1f} mm")
#             print(f"  Distance bord droit  : {dist['distance_bord_droit']:7.1f} mm")
#             print(f"  Distance bord bas    : {dist['distance_bord_bas']:7.1f} mm")
#             print(f"  Distance bord haut   : {dist['distance_bord_haut']:7.1f} mm")
            
#             if dist['distance_semelle'] != "N/A":
#                 print(f"  Distance semelle     : {dist['distance_semelle']:7.1f} mm")
#             if dist['distance_ame'] != "N/A":
#                 print(f"  Distance âme         : {dist['distance_ame']:7.1f} mm")
                
#             if dist['entraxes'] != "N/A":
#                 print(f"  Entraxes             : {', '.join(f'{e} mm' for e in dist['entraxes'])}")
        
#         print("\n" + "="*70 + "\n")
    
#     def show(self, montrer_distances: bool = True):
#         """Affiche la visualisation finale
        
#         Args:
#             montrer_distances (bool): Si True, affiche les distances minimales dans la console
#         """
#         # Ajustement des limites pour inclure tous les éléments
#         margin = max(self.b.value*10**3, self.h.value*10**3) * 0.4
#         self.ax.set_xlim(-self.b.value*10**3/2 - margin/2, self.b.value*10**3/2 + margin/2)
#         self.ax.set_ylim(-self.h.value*10**3/2 - margin/2, self.h.value*10**3/2 + margin/2)
        
#         # Dessin des éléments (d'abord la platine, puis le profilé, puis les boulons)
#         self._dessiner_platine()
#         self._dessiner_profil()
#         self._dessiner_boulons()
        
#         # Affichage des distances si demandé
#         if montrer_distances and hasattr(self, 'boulons') and self.boulons:
#             self.afficher_distances_boulons()
        
#         # Ajustement du layout et affichage
#         plt.tight_layout()
#         plt.show()
    
#     def sauvegarder(self, nom_fichier: str, format_fichier: str = 'png', dpi: int = 300):
#         """
#         Sauvegarde la visualisation dans un fichier.
        
#         Args:
#             nom_fichier (str): Nom du fichier de sortie (sans extension)
#             format_fichier (str): Format de l'image ('png', 'jpg', 'pdf', 'svg', etc.)
#             dpi (int): Résolution en points par pouce
#         """
#         # S'assurer que la figure est à jour
#         self.show()
        
#         # Construction du nom de fichier complet
#         nom_complet = f"{nom_fichier}.{format_fichier}"
        
#         # Sauvegarde
#         self.fig.savefig(nom_complet, dpi=dpi, bbox_inches='tight')
#         print(f"Visualisation sauvegardée sous : {nom_complet}")


class Platine_assise_compression_beton(Tige):
    GAMMA_C = {"Fondamentale": 1.5, "Accidentelle": 1.2}
    CLASSE_CONCRETE = list(Tige._data_from_csv(Tige, "caracteristique_meca_beton.csv").index)[2:]

    def __init__(self, gorge: si.mm, classe_beton: str=CLASSE_CONCRETE, *args, **kwargs):
        """Crée une platine d'assise comprimé avec support béton selon le §6.2.5 de l'EN1993-1-8.
        Cette classe ne prend pas en compte pour le moment les dimensions connues de la fondation ce qui place l'utilisateur du côté de la sécurité.

        Attention:
        Si l’appui est constitué de mortier de scellement, il devra respecter les dispositions suivantes :  
            - De manière générale, l'épaisseur du scellement ne doit pas excéder 0,2 fois la plus petite largeur de la plaque 
            métallique d'assise, et sa résistance caractéristique ne doit pas être inférieure à 0,2 fois la résistance caractéristique 
            fck du béton de fondation. 
            - Dans les cas où l'épaisseur du scellement est supérieure à 50 mm, il convient que la résistance caractéristique du 
            scellement soit au moins égale à celle du béton de fondation.

        On détermine tout d'abord la la dimension effective additionnelle (fonction "c") depuis la face d'une semelle ou d'une ame métallique en mm 
        selon l'EN1993-1-8 §6.2.5(4).
        Puis on utilise la fonction "taux_compression" avec l'aire efficace de compression sur le béton Aef voir EN1993-1-8 §6.2.5(6) 
        et le guide de dimensionnement des assemblages acier édité par le CODIFAB.
        
        Args:
            gorge (int): gorge de la soudure qui relie le T du tronçon
            classe_beton (str): classe du béton
        """
        super().__init__(*args, **kwargs)
        self.gorge = gorge * si.mm
        self.classe_beton = classe_beton
        caract_meca = self._data_from_csv("caracteristique_meca_beton.csv").loc[self.classe_beton]
        self.fck = float(caract_meca.loc["fck"]) * si.MPa

    @property
    def f_jd(self):
        """Valeur de calcul de la résistance à la compression localisée de la plaque d'assise sur le béton
        """
        alpha_cc = 1
        alpha_bf = 1
        beta_i = 1
        f_ck = self.fck
        gamma_c = self.GAMMA_C["Fondamentale"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_jd = beta_i * (alpha_cc * f_ck / gamma_c) * alpha_bf
            return f_jd
        return val()

    def c(self):
        """Détermine la largeur effective complémentaire depuis la face d'une semelle ou d'une ame métallique en mm.
        """
        t = self.t
        f_y = self.fy
        f_jd = self.f_jd[1]
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            C = t * sqrt(f_y / (3 * f_jd * gamma_M0))
            return C
        return val()

    def taux_compression(self, N_c_Ed:si.kN, Aef: si.mm**2):
        """Détermine le taux de travail en compression de la plaque d'assise sur le béton.

        Args:
            N_c_Ed (si.kN): Effort de compression sur la plaque d'assise
            Aef (si.mm): Aire efficace de compression sur le béton tenant compte de la largeur effective complémentaire C
        """
        N_c_Ed = N_c_Ed * si.kN
        A_ef = Aef * si.mm**2
        f_jd = self.f_jd[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_compression = N_c_Ed / (A_ef * f_jd)
            return taux_compression
        result = val()
        self._add_synthese_taux_travail([["Compression sur plaque d'assise et support béton", result[1], None]])
        return result


class Platine_assise_compression_bois(Tige):
    CLASSE_WOOD = list(Tige._data_from_csv(Tige, "caracteristique_meca_bois.csv").index)[2:]

    def __init__(self, wood_beam: object, gorge: si.mm, alpha_wood: float=0, *args, **kwargs):
        """Crée une platine d'assise comprimé avec support bois selon le §6.2.5 de l'EN1993-1-8.
        Pour ce placer en sécurité cette fonction considère un Kc90 = 1 pour le bois.

        On détermine tout d'abord la la dimension effective additionnelle (fonction "c") depuis la face d'une semelle ou d'une ame métallique en mm 
        selon l'EN1993-1-8 §6.2.5(4).
        Puis on utilise la fonction "taux_compression" avec l'aire efficace de compression sur le bois Aef voir EN1993-1-8 §6.2.5(6) 
        et le guide de dimensionnement des assemblages acier édité par le CODIFAB.
        
        Args:
            wood_beam (object): objet correspondant à la classe Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
                             ou bien objet Element ou dérivé de cet objet provenant du module EC3_Element_droit.py
            gorge (int): gorge de la soudure qui relie le T du tronçon
            alpha_wood (str): angle du fil du bois vis à vis de l'effort de compression sur la plaque d'assise.
        """
        super().__init__(*args, **kwargs)
        self.wood_beam = wood_beam
        self.gorge = gorge * si.mm
        self.alpha_wood = alpha_wood

    def _f_jd(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Valeur de calcul de la résistance à la compression localisée de la plaque d'assise sur le bois
        """
        f_c_0_d = self.wood_beam._f_type_d("fc0k", loadtype, typecombi)[1]
        f_c_90_d = self.wood_beam._f_type_d("fc90k", loadtype, typecombi)[1]
        K_c90 = 1 # On se place du coté de la sécurité
        alpha_wood = self.alpha_wood
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_jd =  (f_c_0_d / ((f_c_0_d /(K_c90 * f_c_90_d)) * sin(radians(alpha_wood))**2 + cos(radians(alpha_wood))**2))
            return f_jd
        result = val()
        self.f_jd = result
        return result

    def c(self, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Détermine la largeur effective complémentaire depuis la face d'une semelle ou d'une ame métallique en mm.
        """
        t = self.t
        f_y = self.fy
        f_jd = self._f_jd(loadtype, typecombi)[1]
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        gorge = self.gorge
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            a_ef = 0.8 * (gorge * sqrt(2)) # largeur efficace d'une soudure
            C = t * sqrt(f_y / (3 * f_jd * gamma_M0)) + a_ef
            return C
        return val()

    def taux_compression(self, N_c_Ed:si.kN, Aef: si.mm**2):
        """Détermine le taux de travail en compression de la plaque d'assise sur le bois.

        Args:
            N_c_Ed (si.kN): Effort de compression sur la plaque d'assise
            Aef (si.mm): Aire efficace de compression sur le bois tenant compte de la largeur effective complémentaire C
        """
        N_c_Ed = N_c_Ed * si.kN
        A_ef = Aef * si.mm**2
        latex_fjd, f_jd = self.f_jd
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_compression = N_c_Ed / (A_ef * f_jd)
            return taux_compression
        result = val()
        self._add_synthese_taux_travail([["Compression sur plaque d'assise et support bois", result[1], None]])
        return result
        
    
class Platine_assise_traction(Tige):
    def __init__(self, gorge: si.mm, effet_levier: bool=("False", "True"), *args, **kwargs):
        """Crée une platine d'about fléchie en traction selon le §6.2.6.5 de l'EN1993-1-8 et du CNC2M §4.3.3
        Attention cette classe ne permet le calcul que des platines avec 2 boulons par rangée de tronçon en T équvalent et non 4 comme décrit dans l'annexe A.

        Pour utiliser cette classe on doit d'abord ajouter les rangées isolées de chaque tronçon équivalent avec les fonctions: 
            - ajouter_rangee_simple_unique
            - ajouter_rangee_interieure
            - ajouter_rangee_exterieure_non_raidie
            - ajouter_rangee_centrale_non_raidie.

        Ensuite on ajoute les groupes de rangées pour la vérification des modes de rupture globale avec les fonctions:
            - ajouter_groupe_rangees_interieures_centrales
            - ajouter_groupe_rangees_centrales
            - ajouter_groupe_rangees_interieures

        Puis on vérifie les taux de travail avec la fonction:
            - taux_traction        
        Args:
            gorge (int): gorge de la soudure qui relie le T du tronçon
            effet_levier (bool, optional): définit si un effet levier est appliqué sur la platine. 
                Il est louable de considéré qu'une platine reposant sur un support bois n'a pas d'effet levier, 
                en effet la raideur en compression perpendiculaire du bois est faible. Defaults to "False".
        """
        super().__init__(*args, **kwargs)
        self.gorge = gorge * si.mm
        self.effet_levier = effet_levier
        self._bolts_rows = []
        
    def _m(self, l:si.mm):
        """Dimension m selon la figure 6.2 de l'EN1993-1-8 §6.2.4.1.
        Args:
            l (int): longueur en mm entre la face du T/raidisseur et l'axe de la tige d'assemblage

        Returns:
            int: dimension m en mm
        """
        m = l - 0.8 * self.gorge * sqrt(2)
        return m


    def show_l_eff(self):
        """Affiche le schéma des rangées issue du CNC2M §4.3.3 (4) tab15 et 16 pour définir les longeurs efficaces"""
        file1 = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "CNC2M_4.3.3_tab15.png"
        )
        file2 = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "CNC2M_4.3.3_tab16.png"
        )
        image1 = Image.open(file1)
        image2 = Image.open(file2)
        image1.show()
        image2.show()

    def show_fig6_11(self):
        """Affiche la figure 6.11 de l'EN1993-1-8 pour définir le facteur alpha"""
        file = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "EN1993-1-8_fig6_11.png"
        )
        image = Image.open(file)
        image.show()

    def ajouter_rangee_simple_unique(self, ml:si.mm, e:si.mm, bp:si.mm):
        """Ajoute une rangée isolée simple, unique (pas plus de rangée) et non raidie sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 1 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            ml (si.mm): distance en mm entre la face du T/raidisseur et l'axe de la tige d'assemblage
            e (si.mm): distance de bord en mm
            bp (si.mm): longueur de la platine en mm
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        bp = bp * si.mm

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2 * pi * m
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(4*m + 1.25*e, bp)
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Simple non raidie", "m": m, "e_min": e, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_rangee_interieure(self, ml:si.mm, m2l:si.mm, e:si.mm, alpha: float):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 2 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
            m2l (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
            e (si.mm): distance de bord en mm
            alpha (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda 1 et 2 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        m2l = m2l * si.mm
        m2 = self._m(m2l)
        e = e * si.mm
        lamb1 = m / (m + e)
        lamb2 = m2 / (m + e)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2 * pi * m
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(4*m + 1.25*e, alpha*m)
            return l_eff_nc
        l_eff_nc = leff_nc()

        row = {"Nom": name, "type": "Interieur", "m": m, "e_min": e, "lamb1": lamb1, "lamb2": lamb2, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row
    
    def ajouter_rangee_exterieure_non_raidie(self, w: si.mm, mxl:si.mm, e:si.mm, ex:si.mm, bp:si.mm):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 3 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            w (si.mm): distance d'entraxe des tiges en mm
            mxl (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
            e (si.mm): distance de bord parallèle à la semelle en mm
            ex (si.mm): distance de bord peprpendiculaire à la semelle en mm
            bp (si.mm): longueur de la platine en mm
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        mxl = mxl * si.mm 
        mx = self._m(mxl)
        e = e * si.mm
        ex = ex * si.mm

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = min(2 * pi * mx, pi * mx + w, pi * mx + 2*e)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(4*mx + 1.25*e, 2*mx + 0.625*ex + w/2, 2*mx + 0.625*ex + e, bp/2)
            return l_eff_nc
        l_eff_nc = leff_nc()

        row = {"Nom": name, "type": "Extérieur non raidie", "m": mx, "e_min": min(e, ex), "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_rangee_centrale_non_raidie(self, ml:si.mm, e:si.mm):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 4 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
            e (si.mm): distance de bord en mm
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2 * pi * m
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = 4*m + 1.25*e
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Centrale non raidie", "m": m, "e_min": e, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_rangee_exterieure_raidie(self, mxl:si.mm, ex:si.mm, m2xl:si.mm, e2x:si.mm, alpha1: float, alpha2: float):
        """Ajoute une rangée isolée intérieur sur la platine 
        pour le calcul du l,eff des modes circulaires et non circulaires en traction.
        Schéma 5 dans le tableau 15 du CNC2M (voir la fonction show_l_eff).

        Args:
            mxl (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
            ex (si.mm): distance de bord peprpendiculaire à la semelle en mm
            m2xl (si.mm): distance en mm entre la face du raidisseur et l'axe de la tige d'assemblage
            e2x (si.mm): distance de bord peprpendiculaire au raidisseur en mm
            alpha1 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda11 et lambda21 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
            alpha2 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda12 et lambda22 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon " + str(len(self._bolts_rows)+1)
        mxl = mxl * si.mm 
        mx = self._m(mxl)
        m2xl = m2xl * si.mm
        m2x = self._m(m2xl)
        ex = ex * si.mm
        e2x = e2x * si.mm

        lamb11 = mx / (mx + ex)
        lamb21 = m2x / (mx + ex)
        lamb12 = m2x / (m2x + e2x)
        lamb22 = mx / (m2x + e2x)

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp1 = min(2 * pi * mx, pi * mx + 2*e2x)
            l_eff_cp2 = min(2 * pi * m2x, pi * m2x + 2*ex)
            l_eff_cp = min(l_eff_cp1, l_eff_cp2)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc1 = min(alpha1 * mx, alpha1 * mx - 2*mx - 0.625*ex + e2x)
            l_eff_nc2 = min(alpha2 * m2x, alpha2 * m2x - 2*mx - 0.625*e2x + ex)
            l_eff_nc = min(l_eff_nc1, l_eff_nc2)
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Extérieur raidie", "m": min(mx,m2x), "e_min": min(ex, e2x), "lamb11": lamb11, "lamb21": lamb21, "lamb12": lamb12, "lamb22": lamb22, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    ########### MODE DE RUPTURE GLOBALE ###################

    def ajouter_groupe_rangees_interieures_centrales(self, ml:si.mm, e:si.mm, m2l:si.mm, sum_p:si.mm, alpha: float):
        """Ajoute une rangée isolée intérieur sur la platine 
            pour le calcul du l,eff des modes circulaires et non circulaires en traction.
            Schéma 1 dans le tableau 16 du CNC2M (voir la fonction show_l_eff).

            Args:
                ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
                e (si.mm): distance de bord perpendiculaire à l'ame en mm
                m2l (si.mm): distance en mm entre la face de la semelle et l'axe de la tige d'assemblage
                sum_p (si.mm): somme des entraxes p entres les tiges d'assemblages intérieur/central
                alpha (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                    Pour déterminer alpha on peut récupérer les lambda 1 et 2 en exécutant cette fonction avec un lambda aléatoire, puis
                    en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon groupé " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        m2l = m2l * si.mm
        m2 = self._m(m2l)
        sum_p = sum_p * si.mm

        lamb1 = m / (m + e)
        lamb2 = m2 / (m + e)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2*(pi * m + sum_p)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = min(alpha*m + sum_p, 4*m + 1.25*e + sum_p)
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Intérieure et centrale(s) en groupe", "m": m, "e_min": e, "lamb1": lamb1, "lamb2": lamb2, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_groupe_rangees_centrales(self, ml:si.mm, e:si.mm, sum_p:si.mm):
        """Ajoute une rangée isolée intérieur sur la platine 
            pour le calcul du l,eff des modes circulaires et non circulaires en traction.
            Schéma 2 dans le tableau 16 du CNC2M (voir la fonction show_l_eff).

            Args:
                ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
                e (si.mm): distance de bord perpendiculaire à l'ame en mm
                sum_p (si.mm): somme des entraxes p entres les tiges d'assemblages intérieur/central
        """
        name = "Tronçon groupé " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        sum_p = sum_p * si.mm

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2*(pi * m + sum_p)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = 4*m + 1.25*e + sum_p
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Centrales en groupe", "m": m, "e_min": e, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def ajouter_groupe_rangees_interieures(self, ml:si.mm, e:si.mm, m21l:si.mm, m22l:si.mm, sum_p:si.mm, alpha1: float, alpha2: float):
        """Ajoute une rangée isolée intérieur sur la platine pour le calcul du l,eff 
        des modes circulaires et non circulaires en traction.
        Schéma 3 dans le tableau 16 du CNC2M (voir la fonction show_l_eff).

        Attention: Ce cas ne se rencontre que pour les assemblages de pièces tendues ou pour des encastrements avec raidisseurs entre rangées

        Args:
            ml (si.mm): distance en mm entre la face de l'ame et l'axe de la tige d'assemblage
            e (si.mm): distance de bord perpendiculaire à l'ame en mm
            m21l (si.mm): distance en mm entre la face de la semelle haute et l'axe de la tige d'assemblage
            m22l (si.mm): distance en mm entre la face de la semelle basse et l'axe de la tige d'assemblage
            sum_p (si.mm): somme des entraxes p entres les tiges d'assemblages intérieur/central
            alpha1 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda11 et lambda21 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
            alpha2 (float): facteur définit dans la figure 6.11 de l'EN1993-1-8. 
                Pour déterminer alpha on peut récupérer les lambda12 et lambda22 en exécutant cette fonction avec un lambda aléatoire, puis
                en observant le resultat de sortie pour ensuite redéfinir alpha (voir la fonction show_fig6_11)
        """
        name = "Tronçon groupé " + str(len(self._bolts_rows)+1)
        ml = ml * si.mm 
        m = self._m(ml)
        e = e * si.mm
        m21l = m21l * si.mm
        m21 = self._m(m21l)
        m22l = m22l * si.mm
        m22 = self._m(m22l)
        sum_p = sum_p * si.mm

        lamb11 = m / (m + e)
        lamb21 = m21 / (m + e)
        lamb12 = m / (m + e)
        lamb22 = m22 / (m22 + e)

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_cp():
            l_eff_cp = 2*(pi * m + sum_p)
            return l_eff_cp
        l_eff_cp = leff_cp()

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def leff_nc():
            l_eff_nc = (alpha1 + alpha2)*m -4*m + 1.25*e + sum_p
            return l_eff_nc
        l_eff_nc = leff_nc()
        
        row = {"Nom": name, "type": "Intérieures en groupe", "m": m, "e_min": e, "lamb11": lamb11, "lamb21": lamb21, "lamb12": lamb12, "lamb22": lamb22, "l_eff_cp": l_eff_cp, "l_eff_nc": l_eff_nc}
        self._bolts_rows.append(row)
        return row

    def taux_traction(self, N_T_Ed:si.kN):
        N_T_Ed = N_T_Ed * si.kN
        F_T_Rd_isole = 0 * si.kN
        F_T_Rd_groupe = 0 * si.kN
        n_bl = 2
        t_f = self.t
        f_y = self.fy
        gamma_M0 = self.GAMMA_M["gamma_M0"]
        latex_FtRd, F_t_Rd = self.FtRd
        for row in self._bolts_rows:
            name = row["Nom"]
            m = row["m"]
            e_min = row["e_min"]
            l_eff_cp = row["l_eff_cp"][1]
            l_eff_nc = row["l_eff_nc"][1]

            if self.effet_levier:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    name
                    effet_levier = True
                    l_eff_1 = min(l_eff_cp, l_eff_nc)
                    M_pl1_Rd = (l_eff_1 * t_f**2 * f_y) / (4 * gamma_M0)
                    F_T1_Rd = (4 * M_pl1_Rd) / m
                    
                    l_eff_2 = l_eff_nc
                    n = min(e_min, 1.25*si.m)
                    M_pl2_Rd = (l_eff_2 * t_f**2 * f_y) / (4 * gamma_M0)
                    F_T2_Rd = (2 * M_pl2_Rd + n * n_bl * F_t_Rd) / (m + n)

                    F_T3_Rd = n_bl * F_t_Rd # n,bl est le nombre de boulon dans la rangée en traction

                    F_T_Rd = min(F_T1_Rd, F_T2_Rd, F_T3_Rd)
                    return F_T_Rd
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    name
                    effet_levier = False
                    l_eff_1 = min(l_eff_cp, l_eff_nc)
                    M_pl1_Rd = (l_eff_1 * t_f**2 * f_y) / (4 * gamma_M0)
                    F_T12_Rd = (2 * M_pl1_Rd) / m

                    F_T3_Rd = n_bl * F_t_Rd # n,bl est le nombre de boulon dans la rangée en traction

                    F_T_Rd = min(F_T12_Rd, F_T3_Rd)
                    return F_T_Rd
            row["F,T,Rd"] = val()
            if "groupé" in name:
                if F_T_Rd_groupe == 0*si.kN:
                    F_T_Rd_groupe = row["F,T,Rd"][1]
                else:
                    F_T_Rd_groupe = min(F_T_Rd_groupe, row["F,T,Rd"][1])
            else:
                F_T_Rd_isole = F_T_Rd_isole + row["F,T,Rd"][1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def taux():
            F_T_Rd = min(F_T_Rd_isole, F_T_Rd_groupe)
            taux_traction = N_T_Ed / F_T_Rd
            return taux_traction
        latex = "".join((row["F,T,Rd"][0] for row in self._bolts_rows))
        taux = taux()
        self._add_synthese_taux_travail([["Traction d'une plaque d'assise avec attache en T", taux[1], None]])
        return (latex + taux[0], taux[1])


if __name__ == "__main__":
    # Création d'un pied de poteau HE200B
    pied_poteau = Platine_about(
        h=200,     # hauteur du profilé en mm
        b=200,     # largeur des ailes en mm
        tw=9,      # épaisseur de l'âme en mm
        tf=15,     # épaisseur des ailes en mm
        r=18,      # rayon de raccordement âme-aile en mm
        type_profil="HE"
    )

    # Ajout de boulons M20
    positions_boulons = [
        (-80, 50),   # boulon en haut à gauche
        (80, 50),    # boulon en haut à droite
        (-80, -50),  # boulon en bas à gauche
        (80, -50)    # boulon en bas à droite
    ]
    pied_poteau.ajouter_boulons(positions_boulons, 20)

    # Affichage
    pied_poteau.afficher()