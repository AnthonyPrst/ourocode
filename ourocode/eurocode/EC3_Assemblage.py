# coding in UTF-8 

from math import *
import pandas as pd

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.EC3_Element_droit import Plat

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
        *args, 
        **kwargs
        ):
        """Configure un objet Tige permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Plat du fichier EC3_Element_droit.py.

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
        super().__init__(*args, **kwargs)
        self.d = d * si.mm
        self.d0 = d0 * si.mm
        self.trou_oblong = trou_oblong
        self.qualite = qualite
        self.tete_fraisee = tete_fraisee
        self.prof_tete_fraisee = prof_tete_fraisee * si.mm
        self.verif_filetage = verif_filetage
        self.__verif_trou_surdimensionne()

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
        
        
    def __verif_trou_surdimensionne(self):
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
        else:
            self.percage_surdimensionne = False
        return self.percage_surdimensionne

    @property
    def _Kb(self):
        """Détermine le facteur kb selon CNC2M §2.1 tab2"""
        kb = 1
        if self.trou_oblong:
            kb = 0.6
        elif self.percage_surdimensionne:
            kb = 0.8
        return kb
            


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
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val2():
                    taux_v_d_filetage = F_v_Ed / F_v_Rd_filetage
                    return taux_v_d_filetage
                value2 = val2()
                self.taux_bl["taux_v_d_filetage"] = value2[1]

        elif not self.verif_filetage and FvEd and FtEd:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                taux_v_t_d = F_v_Ed / F_v_Rd_lisse + F_t_Ed / (1.4 * F_t_Rd)
                return taux_v_t_d
            value2 = val2()
            self.taux_bl["taux_v_t_d"] = value2[1]

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
        k_b = self._Kb
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
        return val()

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
        return val()

    def axe_articulation(self, FvEd:si.kN=0, jeu:si.mm=0, t_flasque:si.mm=0, t_ame:si.mm=0):
        """Retourne le taux de travail en double cisaillement et en flexion d'un axe d'articulation selon EN 1993-1-8 §3.13

        Args:
            FvEd (si.kN): effort de cisaillement sur l'axe en kN.
            jeu (si.mm): jeu de l'axe en mm.
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
        return result


class Prescellement_plaque_ancrage(Tige):
    GAMMA_C = {"Fondamentale": 1.5, "Accidentelle": 1.2}
    CLASSE_CONCRETE = list(Tige._data_from_csv(Tige, "caracteristique_meca_beton.csv").index)[2:]
    def __init__(self, tr: si.mm, dr: si.mm, l_ancr: si.mm, d1: si.mm, e_tiges: si.mm, classe_beton: str=CLASSE_CONCRETE, *args, **kwargs):
        """Initialise une tige de prescellement avec plaque circulaire d'après l'EN1993-1-8 et le guide du CNC2M.
        Cette classe est hérité de la classe Tige du module EC3_Assemblage.py.

        Args:
            tr (si.mm): épaisseur de la plaque d'ancrage en mm.
            dr (si.mm): diamètre de la plaque d'ancrage en mm.
            l_ancr (si.mm): longueur de l'ancrage en mm (distance entre le nu béton et la face de la plaque d'ancrage).
            d1 (si.mm): distance perpendiculaire du bord béton à l'axe de la tige en mm.
            e_tiges (si.mm): entraxe des tiges d'ancrage en mm.
            classe_beton (str): classe du béton.
        """
        kwargs["filetage_EN1090"] = False
        super().__init__(*args, **kwargs)
        self.tr = tr * si.mm
        self.dr = dr * si.mm
        self.l_ancr = l_ancr * si.mm
        self.d1 = d1 * si.mm
        self.e_tiges = e_tiges * si.mm
        self.classe_beton = classe_beton
        self.type_ancrage = "Plaque d'ancrage circulaire"

        caract_meca = self._data_from_csv("caracteristique_meca_beton.csv").loc[self.classe_beton]
        self.fck = float(caract_meca.loc["fck"]) * si.MPa

        if self.tr.value < 0.3*self.dr.value/2:
            raise ValueError("L'épaisseur tr doit être supérieur ou égal à {} mm".format(ceil(0.3*self.dr.value*10**3/2)))

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
            raise ValueError("La limite elastique Fyb doit être comprise entre 235 MPa et 640 MPa pour être considéré en cisaillement")

    @property
    def FtcRd(self):
        d = self.d
        d_r = self.dr
        l_ancr = self.l_ancr
        d_1 = self.d1
        e_tiges = self.e_tiges
        gamma_c = self.GAMMA_C["Fondamentale"]
        f_c_k = self.fck
        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            r_r = d_r / 2
            v = min(l_ancr, d_1, e_tiges)
            A_ancr = 2.55 * pi * (r_r**2 - d**2/4) * (1-r_r/v)
            F_tc_Rd = A_ancr * f_c_k / gamma_c #CNC2M §6(6) tab19
            return F_tc_Rd
        return val()

    def taux_ancrage(self, FtEd: si.kN=0, FvEd: si.kN=0):
        """Calcul le taux de travail d'une tige d'ancrage avec plaque circulaire d'après l'EN1993-1-8 et le guide du CNC2M.

        Args:
            FtEd (float): effort de traction sur la tige en kN.
            FvEd (float): effort de cisaillement sur la tige en kN. 
                Il est en général déconseillé de faire passer de grand effort de cisaillement par les tiges, privilégier une bèche.

        Returns:
            tuple: le tuple des taux de travail
        """
        F_t_Ed = FtEd * si.kN
        F_v_Ed = FvEd * si.kN
        latex_F_tc_Rd, F_tc_Rd = self.FtcRd
        latex_F_tb_Rd, F_tb_Rd = self.FtRd
        if not F_v_Ed or F_v_Ed.value == 0:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_ancrage = F_t_Ed / min(F_tc_Rd, F_tb_Rd)
                return taux_ancrage
            result = val()
            return (latex_F_tb_Rd + latex_F_tc_Rd + result[0], result[1])
        elif not self.percage_surdimensionne:
            latex_F_v_Rd, F_v_Rd = self.FvRd
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_t_ancr_Rd = min(F_tc_Rd, F_tb_Rd)
                taux_t = F_t_Ed / F_t_ancr_Rd #CNC2M §6(8)
                taux_t_v = F_v_Ed / F_v_Rd + F_t_Ed / (1.4 * F_t_ancr_Rd)
                return taux_t, taux_t_v
            result = val()
            return (latex_F_v_Rd + latex_F_tb_Rd + latex_F_tc_Rd + result[0], result[1])
        else:
            t_r = self.tr
            d = self.d
            latex_F_v_Rd, F_v_Rd = self.FvRd
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                e = t_r + d / 2
                F_t_Ed_eq = F_v_Ed * (e / d) * (5 * pi)/6 #CNC2M §6(9)
                F_t_ancr_Rd = min(F_tc_Rd, F_tb_Rd)
                taux_t = (F_t_Ed + F_t_Ed_eq) / F_t_ancr_Rd #CNC2M §6(9)
                taux_t_v = F_v_Ed / F_v_Rd + (F_t_Ed + F_t_Ed_eq) / (1.4 * F_t_ancr_Rd)
                return taux_t, taux_t_v
            result = val()
            return (latex_F_v_Rd + latex_F_tb_Rd + latex_F_tc_Rd + result[0], result[1])

        


class Soudure(Plat):
    def __init__(self, t2: si.mm, gorge: si.mm, l: si.mm, retour_soudure: bool=("False", "True"), alpha: float=90, *args, **kwargs):
        """Configure un objet soudure permettant les vérification suivant l'EN 1993-1-8. Cette classe est hérité de la classe Plat du fichier EC3_Element_droit.py.  

        Args:
            t2 (float): épaisseur de la pièce 2 sur laquelle on soude en mm
            gorge (int): dimension de la gorge "a" en mm
            l (int): longueur de soudure "brute" sans cratère en mm
            retour_soudure (bool, optional): détermine si un retour de la soudure est fait si oui alors True. Defaults to False.
            alpha (int | float, optional): angle en degré de la de la pièce 2 sur la pièce 1. Defaults to 90.
        """
     
        super().__init__(*args, **kwargs)
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
        return val()


    def cordon_laterale(self, V_Ed: si.kN=0):
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
        return val()


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
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            cordon_oblique = (beta_w * gamma_M_2 * (N_Ed * sqrt(3 - sin(radians(self.alpha_cordon))**2)) / fu) / (a * l_ef)
            return cordon_oblique
        return val()


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
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            if self.alpha < 90:
                cordon_pieces_obliques = (beta_w * gamma_M_2 * (N_Ed * sqrt(2 - sin(radians(self.alpha)))) / fu) / (a * l_ef)
            elif self.alpha > 90:
                cordon_pieces_obliques = (beta_w * gamma_M_2 * (N_Ed * sqrt(2 + sin(radians(self.alpha)))) / fu) / (a * l_ef)
            return cordon_pieces_obliques
        return val()


    def critere_generale(self, FvEd: si.kN=0, FaxEd: si.kN=0):
        """Calcul le critère générale de Von Mises d'une soudure et retourne le taux.

        Args:
            FvEd (int | float): Effort de cisaillement sur la en kN
            FaxEd (int | float): Effort de traction sur la soudure en kN

        Returns:
            (float) : Taux de travail de la soudure
        """
        Fv_Ed = FvEd * si.kN
        Fax_Ed = FaxEd * si.kN
        self.tau_para = Fv_Ed / (self.gorge * self.lef)
        self.sigma_perpend = (Fax_Ed * cos(radians(self.alpha/2)))/ (self.gorge * self.lef)
        self.tau_perpend = (Fax_Ed * cos(radians(self.alpha/2)))/ (self.gorge * self.lef)
        self.sigma_e = sqrt(self.sigma_perpend**2 + 3 * (self.tau_perpend**2 + self.tau_para**2))

        self.sigma_Rd = self.fu / (self.beta_w * self.GAMMA_M["gamma_M2"])
        self.sigma_perpend_Rd = (0.9 * self.fu) / self.GAMMA_M["gamma_M2"]
        return max(self.sigma_e / self.sigma_Rd, self.sigma_perpend / self.sigma_perpend_Rd)
    
    

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
            print("Il n'est pas possible d'avoir une soudure discontinue en ambiance corrosive")
            return False
        b = b * si.mm
        lwe = max(0.75 * b, 0.75 * b1*si.mm)
        l1 = min(16 * self.t, 16 * self.t2, 200)
        l2 = min(12 * self.t, 12 * self.t2, 0.25 * b, 200)
        return {"Lwe": lwe, "L1": l1, "L2": l2}

    
class Platine_about(Plat):
    def __init__(self, r_or_a: int, l: int, attache_T: str=("Soudée", "Profil formé à froid/chaud"), effet_levier: bool=("False", "True"), *args, **kwargs):
        """Crée une platine flechie

        Args:
            r_or_a (int): rayon du profilé ou gorge de la soudure qui relie le T du tronçon
            l (int): longueur en mm entre la face du T et l'axe de la tige d'assemblage
            attache_T (str, optional): attache de la platine. Defaults to "Soudée".
            effet_levier (bool, optional): défini si un effet levier est appliqué sur la platine. Defaults to "False".
        """
        super().__init__(*args, **kwargs)
        self.attache_T = attache_T
        self.r_or_a = r_or_a
        self.l = l
        self.effet_levier = effet_levier
        
    @property
    def _m(self):
        if self.attache_T == "Soudée":
            m = self.l - 0.8 * self.r_or_a * sqrt(2)
        else:
            m = self.l - 0.8 * self.r_or_a
        return m


    @property
    def l_eff(self):
        """Longueur effective de la platine pour une platine d'about"""
        dict_l_eff = {
            "Rangée de boulons prise séparément": {"Mécanismes circulaires": {}, "Mécanismes non circulaires": {}},
            "Rangée de boulons prise ensemble": {"Mécanismes circulaires": {}, "Mécanismes non circulaires": {}}
            }


if __name__ == "__main__":
    ele = Plat(classe_acier="S235", t=6, h=200, classe_transv=1)
    soudure = Soudure._from_parent_class(ele, t2=6, gorge=4, l=140, retour_soudure=True, alpha=90)
    bl = Tige._from_parent_class(ele, d=12,d0=14,qualite="A2-50",verif_filetage=True, filetage_EN1090=True)
    print(soudure.critere_generale(0, 100135))
    N_rd = bl.Veff_Rd(300, 300, "Centré")
    V_Rd = bl.Veff_Rd(300, 300, "Centré")
    taux = bl.taux_Veff_d(12, N_rd[1], 25, V_Rd[1])
    print(taux)