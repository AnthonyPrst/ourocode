# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.core.projet import Projet
from ourocode.eurocode.ec5.element_droit.barre import Barre
from ourocode.eurocode.ec5.element_droit.cisaillement import Cisaillement
from ourocode.eurocode.ec5.element_droit.compression import Compression_perpendiculaire, Compression_inclinees

def interpolationLineaire(x, xa, xb, ya, yb) -> float:
    """Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
    y = ya + (x - xa) * ((yb - ya)/(xb - xa))
    return y
    
# ====================================================== GENERAL =========================================================
# 8.1 Généralité
class Assemblage(Projet):
    GAMMA_M_ASS = 1.3
    DICO_COEF_LIMITE = {"Pointe circulaire lisse": 0.15, "Agrafe": 0.15,"Pointe carrée lisse": 0.25,
                     "Boulon": 0.25, "Autres pointes": 0.5, "Tirefond": 1}
    TYPE_BOIS_ASSEMBLAGE = ("Bois","PP/OSB", "CP", "Panneau dur")
    TYPE_ASSEMBLAGE = ("Bois/Bois", "Bois/Métal")

    def __init__(self, beam_1:object, beam_2:object, nfile: int=1, nCis: int=("1","2"), **kwargs):
        """Créer un objet Assemblage qui permet de calculer un assemblage bois/bois ou bois/métal à l'EN 1995.
        Cette classe est hérité de la classe Projet du module A0_Project.py

        Args:
            beam_1 (object): objet correspondant à i=1, Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
                             ou bien objet Element ou dérivé de cet objet provenant du module EC3_Element_droit.py
                             
            beam_2 (object): objet correspondant à i=2, Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
                             ou bien objet Element ou dérivé de cet objet provenant du module EC3_Element_droit.py
                             
            nfile (int, optional): le nombre de file dans l'assemblage en considérant i=1. Defaults to 1.
            nCis (int, optional): Nombre de plan cisaillé entre 1 et 2. Defaults to ["1","2"].
        """
        
        super().__init__(**kwargs)
        self.beam_1 = beam_1
        self.beam_2 = beam_2
        self.nfile = nfile
        self.nCis = nCis
        self.__type_assemblage()
    

    def __type_assemblage(self):
        self._type_beam = []
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            try:
                if beam.type_bois:
                    if beam.type_bois in ["Massif","BLC", "LVL"]:
                        self._type_beam.append(self.TYPE_BOIS_ASSEMBLAGE[0])
                    elif beam.type_bois in ["OSB 2", "OSB 3/4"]:
                        self._type_beam.append(self.TYPE_BOIS_ASSEMBLAGE[1])
                    elif beam.type_bois == "CP":
                        self._type_beam.append(self.TYPE_BOIS_ASSEMBLAGE[2])
                    elif beam.type_bois == "Panneau dur":
                        self._type_beam.append(self.TYPE_BOIS_ASSEMBLAGE[3])

                    beam.rho_k = int(beam.caract_meca.loc["rhok"])
            except AttributeError:
                self._type_beam.append("Métal")

        # Détermine le type d'assemblage Bois/Bois ou Bois/métal  
        if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
            if self._type_beam[1] in self.TYPE_BOIS_ASSEMBLAGE:
                self.type_assemblage = self.TYPE_ASSEMBLAGE[0]
            else:
                self.type_assemblage = self.TYPE_ASSEMBLAGE[1]

        else:
            self.type_assemblage = self.TYPE_ASSEMBLAGE[1]
        

    @property
    def rho_mean_ass(self) -> float:
        if self.type_assemblage == self.TYPE_ASSEMBLAGE[0]:
            rho_m1 = int(self.beam_1.caract_meca.loc["rhomean"])
            rho_m2 = int(self.beam_2.caract_meca.loc["rhomean"])
            return mt.sqrt(rho_m1 * rho_m2)
        else:
            if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
                return int(self.beam_1.caract_meca.loc["rhomean"])
            else:
                return int(self.beam_2.caract_meca.loc["rhomean"])
           
            
    def _min_nef(self, list_nef: list):
        """Retourne le nef minimum entre les 2

        Args:
            list_nef (list): list des nef en format handcalcs

        Returns:
            _type_: _description_
        """
        nef_1 = list_nef[0]
        nef_2 = list_nef[1]
        
        min_nef = min(self.nfile * nef_1[1], self.n * nef_2[1])
        if min_nef == (self.nfile * nef_1[1]):
            return list_nef[0]
        else:
            n_file = deepcopy(self.nfile)
            self.nfile = self.n
            self.n = n_file
            warnings.warn("Le sens de traitement de l'assemblage a été changé car le nef mini ce trouve sur la pièce 2 et non sur la pièce 1.\nAttention aux efforts de calcul à prendre en compte.")
            return list_nef[1]
        
        
    # 7.1 Glissement des assemblages
    def Kser(self) -> tuple:
        """Retourne le module de glissement de l'organe pour les états limites de services."""
        rho_mean = self.rho_mean_ass
        if self.type_organe == "Anneau" or self.type_organe == "Crampon C10/C11":
            dc = self.dc.value*10**3
        else: 
            d = self.d.value*10**3
            
        if self.type_organe == "Boulon" or self.type_organe == "Broche" or self.type_organe == "Tirefond":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_ser = rho_mean**1.5 * d / 23 # N/mm
                return K_ser * si.N / si.mm
        elif self.type_organe == "Pointe circulaire lisse" or self.type_organe == "Pointe carrée lisse" or self.type_organe == "Autres pointes":
            if self.percage:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_ser = rho_mean**1.5 * d / 23 # N/mm
                    return K_ser * si.N / si.mm
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    K_ser = rho_mean**1.5 * d**0.8 / 30 # N/mm
                    return K_ser * si.N / si.mm
        elif self.type_organe == "Agrafe":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_ser = 2 * rho_mean**1.5 * self.d**0.8 / 80 # N/mm
                return K_ser * si.N / si.mm
        elif self.type_organe == "Anneau" or self.type_organe == "Crampon C10/C11":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_ser = rho_mean * dc / 2 # N/mm
                return K_ser * si.N / si.mm
        elif self.type_organe == "Crampon C1/C9":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                K_ser = 1.5 * rho_mean * d / 4 # N/mm
                return K_ser * si.N / si.mm
        return val()

    def Kser_ass(self) -> tuple:
        """Retourne le module de glissement de l'assemblage pour les états limites de services."""
        n_file = self.nfile
        K_ser = self.Kser()[1]
        n = self.n
        n_Cis = self.nCis
        k_type = 1
        if self.type_assemblage == self.TYPE_ASSEMBLAGE[1]:
            k_type = 2
        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():    
            K_ser_ass = K_ser * n_file * n * n_Cis * k_type
            return K_ser_ass
        return val()

    def Ku(self) -> tuple:
        """Retourne le module de glissement de l'organe pour les états limites ultimes."""
        K_ser = self.Kser()[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            K_u = K_ser * 2 / 3
            return K_u
        return val()

    def Ku_ass(self) -> tuple:
        """Retourne le module de glissement de l'assemblage pour les états limites ultimes."""
        n_file = self.nfile
        K_u = self.Ku()[1]
        n = self.n
        n_Cis = self.nCis
        k_type = 1
        if self.type_assemblage == self.TYPE_ASSEMBLAGE[1]:
            k_type = 2
        @handcalc(override="short", precision=3, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():    
            K_u_ass = K_u * n_file * n * n_Cis * k_type
            return K_u_ass
        return val()

    # 8.1.4 Effort d'assemblage inclinés par rapport au fil

    def _w(self, wpl: int=0):
        """Calcul le facteur de modification pour le calcul de la valeur caractéristique au fendage F90,Rk avec:
            wpl : largeur de la plaque métallique emboutie parallèlement au fil en mm
            type_organe : type d'organe utilisé, pour les plaques métalliques embouties : "plaque", pour les autres : "autres" """
        if self.type_organe == "plaques métaliques embouties":
            w = max((wpl / 100) ** 0.35, 1)
        else:
            w = 1
        return w


    def F90Rk(self, b:si.mm, h:si.mm, he:si.mm, w: int=1) -> tuple:
        """Calcul la valeur caractérisque de la capacité au fendage en N

        Args:
            b (int) : l'épaisseur de l'élément en mm
            h (int) : la hauteur de l'élément en mm
            he (int) : la distance de rive chargée vis à vis du centre de l'organe le plus plus éloigné ou du bord de la plaque
            w (int) : facteur de modification """
        h_e = he
        @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_90_Rk = 14 * b * w * sqrt(h_e / (1 - (h_e / h))) # N
            return F_90_Rk * si.N
        value = val()
        self.F90_Rk = value[1]
        return value

    # 8.2 Capacité résistante latérale pour les organes métaliques de type tige
    # 8.2.2 Assemblage bois/bois bois/panneau


    def _FvRk_BoisBois(self, effet_corde: bool):
        """Calcul la capacité résistante en cisaillement de la tige en N par plan de cisaillement avec
            t1 : valeur minimale entre epaisseur de l'élément bois latéral et la profondeur de pénétration en mm
            t2 : epaisseur de l'élément bois central en mm """
            
        f_h1k = self.fh1k[1].value * 10**-6
        f_h2k = self.fh2k[1].value * 10**-6
        t_1 = self.t1.value * 10**3
        t_2 = self.t2.value * 10**3
        diam = self.d.value * 10**3
        M_y_Rk = self.MyRk[1].value * 10**3

        coef_limit_Johansen = self.DICO_COEF_LIMITE.get(self.type_organe, 0)
        if effet_corde:
            F_ax_Rk = self.FaxRk.value
        else:
            F_ax_Rk = 0

        if self.nCis == 1:
            
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                effet_corde = F_ax_Rk/4 # N
                beta = f_h2k / f_h1k
                a = f_h1k * t_1 * diam # N
                b = f_h2k * t_2 * diam # N
                c_johansen = a/(1 + beta) * (sqrt(beta + 2 * beta**2 * (1 + t_2 / t_1 + (t_2 / t_1)**2) + beta**3 * (t_2 / t_1)**2) - beta * (1 + t_2 / t_1)) # N
                c = c_johansen + min(effet_corde, c_johansen*coef_limit_Johansen) # N
                d_johansen = 1.05 * a / (2 + beta) * (sqrt(2 * beta * (1 + beta) + (4 * beta * (2 + beta) * M_y_Rk)/(f_h1k * t_1**2 * diam)) - beta) # N
                d = d_johansen + min(effet_corde, d_johansen*coef_limit_Johansen) # N
                e_johansen = 1.05 * (f_h1k * t_2 * diam) / (1 + 2 * beta) * (sqrt(2 * beta**2 * (1 + beta) + (4 * beta * (1 + 2 * beta) * M_y_Rk)/(f_h1k * t_2**2 * diam)) - beta) # N
                e = e_johansen + min(effet_corde, e_johansen*coef_limit_Johansen) # N
                f_johansen = 1.15 * sqrt((2 * beta)/(1 + beta)) * sqrt(2 * M_y_Rk * f_h1k * diam) # N
                f = f_johansen + min(effet_corde, f_johansen*coef_limit_Johansen) # N
                return a, b, c, d, e, f
            
            calcul = val()
            a = calcul[1][0] * si.N
            b = calcul[1][1] * si.N
            c = calcul[1][2] * si.N
            d = calcul[1][3] * si.N
            e = calcul[1][4] * si.N
            f = calcul[1][5] * si.N
            dicoRupture = {a: "A", b: "B", c: "C", d: "D", e: "E", f: "F"}
            FvRk = min(a, b, c, d, e, f)
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(a, b, c, d, e, f)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()

        else:
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                effet_corde = F_ax_Rk/4 # N
                beta = f_h2k / f_h1k
                g = f_h1k * t_1 * diam # N
                h = 0.5 * f_h2k * t_2 * diam # N
                j_johansen = 1.05 * g / (2 + beta) * (sqrt(2 * beta * (1 + beta) + (4 * beta * (2 + beta) * M_y_Rk)/ (f_h1k * t_1**2 * diam)) - beta) # N
                j = j_johansen + min(effet_corde, j_johansen*coef_limit_Johansen) # N
                k_johansen = 1.15 * sqrt((2 * beta)/(1 + beta)) * sqrt(2 * M_y_Rk * f_h1k * diam) # N
                k = k_johansen + min(effet_corde, k_johansen*coef_limit_Johansen) # N
                return g, h, j, k
            
            calcul = val()
            g = calcul[1][0] * si.N
            h = calcul[1][1] * si.N
            j = calcul[1][2] * si.N
            k = calcul[1][3] * si.N
            FvRk = min(g, h, j, k)
            dicoRupture = {g: "G", h: "H", j: "J", k: "K"}
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(g, h, j, k)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()
            
        return (calcul[0] + calcul2[0], calcul2[1])

    def _type_plaque(self):
        if self.t <= 0.5 * self.d:
            self.type_plaque = "mince"

        elif self.d <= self.t:
            self.type_plaque = "epaisse"
            
        else:
            self.type_plaque = "intermédiaire"
        return self.type_plaque

    # 8.2.3 Assemblage bois métal
    def _FvRk_BoisMetal(self, effet_corde: bool, test_type_plaque:bool=True):
        """Calcul la capacité résistante en cisaillement de la tige en N par plan de cisaillement avec
            t1 : valeur minimale entre epaisseur de l'élément bois latéral et la profondeur de pénétration en mm
            t2 : epaisseur de l'élément bois central en mm
            """
        diam = self.d.value * 10**3
        M_y_Rk = self.MyRk[1].value * 10**3
        coef_limit_Johansen = self.DICO_COEF_LIMITE.get(self.type_organe, 0)
        if effet_corde:
            F_ax_Rk = self.FaxRk.value
        else:
            F_ax_Rk = 0

        if self._type_beam[0] == "Métal":
            f_h2k = self.fh2k[1].value * 10**-6
            t_2 = self.t2.value * 10**3
            self.t = self.beam_1.t
            self.pos_plaq = "externe"
        else:
            f_h1k = self.fh1k[1].value * 10**-6
            t_1 = self.t1.value * 10**3
            self.t = self.beam_2.t
            if self.nCis == 2:
                self.pos_plaq = "centrale"
            else:
                self.pos_plaq = "externe"

        # On détecte le type de plaque
        if test_type_plaque:
            self._type_plaque()

        # Si la plaque est intermédiaire, alors on fait une interpolation linéaire entre la valeur critique d'une plaque mince et d'une plaque epaisse
        if self.type_plaque == "intermédiaire" and self.pos_plaq != "centrale":
            t = self.t
            t_mince = 0.5 * self.d
            t_epaisse = self.d
            self.type_plaque = "mince"
            FvRk_mince = self._FvRk_BoisMetal(effet_corde,test_type_plaque=False)
            F_v_Rk_mince = FvRk_mince[1][0]
            self.type_plaque = "epaisse"
            FvRk_epaisse = self._FvRk_BoisMetal(effet_corde,test_type_plaque=False)
            F_v_Rk_epaisse = FvRk_epaisse[1][0]
            self.type_plaque = "intermédiaire"
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_v_Rk_inter = F_v_Rk_mince + (t - t_mince) * ((F_v_Rk_epaisse - F_v_Rk_mince)/(t_epaisse - t_mince))  # interpolation entre plaque mince et épaisse
                return F_v_Rk_inter
            FvRk_inter = val()
            return (FvRk_mince[0] + FvRk_epaisse[0] + FvRk_inter[0], (FvRk_inter[1], (FvRk_mince[1][1], FvRk_epaisse[1][1])))

        if self.type_plaque == "mince" and self.nCis == 1:
            if self._type_beam[0] == "Métal":
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    effet_corde = F_ax_Rk/4 # N
                    a = 0.4 * f_h2k * t_2 * diam # N
                    b_johansen = 1.15 * sqrt(2 * M_y_Rk * f_h2k * diam) # N
                    b = b_johansen + min(effet_corde, b_johansen*coef_limit_Johansen) # N
                    return a, b
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    effet_corde = F_ax_Rk/4 # N
                    a = 0.4 * f_h1k * t_1 * diam # N
                    b_johansen = 1.15 * sqrt(2 * M_y_Rk * f_h1k * diam) # N
                    b = b_johansen + min(effet_corde, b_johansen*coef_limit_Johansen) # N
                    return a, b
            
            calcul = val()
            a = calcul[1][0] * si.N
            b = calcul[1][1] * si.N
            dicoRupture = {a: "A", b: "B"}
            FvRk = min(a, b)
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(a, b)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()

        elif self.type_plaque == "epaisse" and self.nCis == 1:
            if self._type_beam[0] == "Métal":
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    effet_corde = F_ax_Rk/4 # N
                    c = f_h2k * t_2 * diam # N
                    d_johansen = c * (sqrt(2 + (4 * M_y_Rk) / (f_h2k * diam * t_2 ** 2)) - 1) # N
                    d = d_johansen + min(effet_corde, d_johansen*coef_limit_Johansen) # N
                    e_johansen = 2.3 * sqrt(M_y_Rk * f_h2k * diam) # N
                    e = e_johansen + min(effet_corde, e_johansen*coef_limit_Johansen) # N
                    return c, d, e
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    effet_corde = F_ax_Rk/4 # N
                    c = f_h1k * t_1 * diam # N
                    d_johansen = c * (sqrt(2 + (4 * M_y_Rk) / (f_h1k * diam * t_1 ** 2)) - 1) # N
                    d = d_johansen + min(effet_corde, d_johansen*coef_limit_Johansen) # N
                    e_johansen = 2.3 * sqrt(M_y_Rk * f_h1k * diam) # N
                    e = e_johansen + min(effet_corde, e_johansen*coef_limit_Johansen) # N
                    return c, d, e
            
            calcul = val()
            c = calcul[1][0] * si.N
            d = calcul[1][1] * si.N
            e = calcul[1][2] * si.N
            dicoRupture = {c: "C", d: "D", e: "E"}
            FvRk = min(c, d, e)
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(c, d, e)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()

        elif self.nCis == 2 and self.pos_plaq == "centrale":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                effet_corde = F_ax_Rk/4 # N
                f = f_h1k * t_1 * diam # N
                g_johansen = f * (sqrt(2 + (4 * M_y_Rk) / (f_h1k * diam * t_1 ** 2)) - 1) # N
                g = g_johansen + min(effet_corde, g_johansen*coef_limit_Johansen) # N
                h_johansen = 2.3 * sqrt(M_y_Rk * f_h1k * diam) # N
                h = h_johansen + min(effet_corde, h_johansen*coef_limit_Johansen) # N
                return f, g, h

            calcul = val()
            f = calcul[1][0] * si.N
            g = calcul[1][1] * si.N
            h = calcul[1][2] * si.N
            dicoRupture = {f: "F", g: "G", h: "H"}
            FvRk = min(f, g, h)
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(f, g, h)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()

        elif self.type_plaque == "mince" and self.nCis == 2 and self.pos_plaq != "centrale":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                effet_corde = F_ax_Rk/4 # N
                j = 0.5 * f_h2k * t_2 * diam # N
                k_johansen = 1.15 * sqrt(2 * M_y_Rk * f_h2k * diam) # N
                k = k_johansen + min(effet_corde, k_johansen*coef_limit_Johansen) # N
                return j, k

            calcul = val()
            j = calcul[1][0] * si.N
            k = calcul[1][1] * si.N
            dicoRupture = {j: "J", k: "K"}
            FvRk = min(j, k)
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(j, k)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()

        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                effet_corde = F_ax_Rk/4 # N
                l = 0.5 * f_h2k * t_2 * diam # N
                m_johansen = 2.3 * sqrt(M_y_Rk * f_h2k * diam) # N
                m = m_johansen + min(effet_corde, m_johansen*coef_limit_Johansen) # N
                return l, m

            calcul = val()
            l = calcul[1][0] * si.N
            m = calcul[1][1] * si.N
            dicoRupture = {l: "L", m: "M"}
            FvRk = min(l, m)
            mode_rupture = dicoRupture[FvRk]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_v_Rk = min(l, m)
                mode_rupture
                return F_v_Rk, mode_rupture
            calcul2 = val2()

        return (calcul[0] + calcul2[0], calcul2[1])

    
        # 8.1.2 Assemblage par organe multiple

    def FvRk(self, effet_corde: bool=("True", "False")) -> tuple:
        """Calcul la valeur de calcul caractéristique de résistance au cisaillement de l'assemblage en N

        Args:
            effet_corde (bool): prise en compte de l'effet de corde, si oui alors True.
            Attention: pour que l'effet de corde soit prit en compte, il faut que le FaxRk du type de tige est été calculé préalablement.

        Returns:
            float: effort de reprise caractéristique de l'assemblage en N
        """
        #     Fvrktot : capacité résistante en cisaillement caractéristique avec la partie de Johansen + l'effet de corde en N
        if self.type_assemblage == self.TYPE_ASSEMBLAGE[0]:
            latex, self.Fv_Rk = self._FvRk_BoisBois(effet_corde)
        else:
            latex, self.Fv_Rk = self._FvRk_BoisMetal(effet_corde)

        F_v_Rk = self.Fv_Rk[0]
        n_file = self.nfile
        n_ef = self._nef
        n_cisaillement = self.nCis
        if self.type_organe == "Agrafe":
            if not self.angle_sup_30:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    F_v_Rk_ass = F_v_Rk * 0.7 * n_file * n_ef * n_cisaillement
                    return F_v_Rk_ass
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    F_v_Rk_ass = F_v_Rk * 2 * n_file * n_ef * n_cisaillement
                    return F_v_Rk_ass
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_v_Rk_ass = F_v_Rk * n_file * n_ef * n_cisaillement
                return F_v_Rk_ass
        FvRkass_latex, self.Fv_Rk_ass = val()
        latex = latex + FvRkass_latex
        return (latex, self.Fv_Rk_ass)


    def F_Rd(self, F_Rk: si.kN, loadtype=Barre.LOAD_TIME) -> tuple:
        """Calcul la valeur de calcul (design) de résistance de l'assemblage en N avec :

        Args:
            F_rk (float): capacité résistante caractéristique de l'organe en kN
            loadtype (str, optional): Durée de chargement (Permanente, Court terme etc.). Defaults to "Permanente"."""
        F_Rk = F_Rk * si.kN
        gamma_M = self.GAMMA_M_ASS
        if "Métal" in self._type_beam:
            if self._type_beam[0] == "Métal":
                k_mod = self.beam_2._get_k_mod(loadtype)
            else:
                k_mod = self.beam_1._get_k_mod(loadtype)
        else:
            k_mod = sqrt(self.beam_1._get_k_mod(loadtype) * self.beam_2._get_k_mod(loadtype))
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_Rd = (F_Rk * k_mod) /gamma_M # Valeur de calcul (design)
            return F_Rd
        return val()
    
    
    # Annexe A : Cisaillement de bloc
    def FbsRk(self, dp:si.mm, a1:si.mm, a2:si.mm, a3t:si.mm, Kcr: float=0.67, num_beam: int=("1", "2")) -> tuple:
        """Calcul la valeur caractéristique en cisaillement de bloc en N pour l'élément 1 ou 2 de l'assemblage.
        Attention pour que cette fonction puisse s'éxecuter il faut avoir préalablement éxécuté la fontion FvRk.

        Args:
            dp (float): diamètre de perçage en mm
            a1 (float): pince longitudinale en mm
            a2 (float):  pince perpendiculaire en mm
            a3t (float): pince en bord chargée suivant le file en mm
            kcr (float, optional): coeff de réduction largeur en cisaillement. 
                Si aucune valeur n'est rentrée alors il sera automatiquement calculé. Defaults to 0.67.
            num_beam (int, optional): numéro de l'élément à vérifier. Defaults to 1.
        """
        diam = self.d
        diam_percage = dp * si.mm
        
        a_1 = a1 * si.mm
        a_2 = a2 * si.mm
        a_3_t = a3t * si.mm
        M_y_Rk = self.MyRk[1]
        n = self.n
        n_file = self.nfile

        mode_rupture = self.Fv_Rk[1]

        # on calcul les valeur caractéristique de base de notre éléméent
        if num_beam == 1:
            parent_beam = self.beam_1
            f_hk = self.fh1k[1]
            t_1 = self.t1
        else:
            parent_beam = self.beam_2
            f_hk = self.fh2k[1]
            t_1 = self.t2
        f_t0_k = float(parent_beam.caract_meca.loc["ft0k"]) * si.MPa
        f_v_k = float(parent_beam.caract_meca.loc["fvk"]) * si.MPa

        # On calcul le kcr si aucune valeur n'est donnée en argument
        if not Kcr:
            cis_beam = Cisaillement._from_parent_class(parent_beam)
            K_cr = cis_beam.K_cr
        else:
            K_cr = Kcr

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def l_net():
            l_v_i = a_1 - diam_percage # Distance entre perçage dans le sens fil
            L_net_v = 2 * (l_v_i * (n - 1) + a_3_t - (diam_percage / 2)) # Longueur résiduelle de la surface de rupture en cisaillement
            
            l_t_i = a_2 - diam_percage # Distance entre perçage dans le sens perpendiculaire au fil
            L_net_t = l_t_i * (n_file - 1) # Largeur résiduelle de la section perpendiculaire au fil

            mode_rupture
            return L_net_v, L_net_t
        
        L_net_value = l_net()
        latex = L_net_value[0]
        L_net_v, L_net_t = L_net_value[1]

        if not isinstance(mode_rupture, (list, tuple)):
            mode_rupture = [mode_rupture]

        def _Anet_v(mode: str):
            """ détermine la surface Anet_v en fonction du mode de rupture """
            if mode in ("C", "F", "J", "L", "K", "M"):
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def a_net_v():
                    A_net_v = L_net_v * (t_1 * K_cr)
                    return A_net_v
                a_net_v_result = a_net_v()
                return a_net_v_result
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def a_net_v(t_ef):
                    A_net_v = L_net_v / 2 * (L_net_t + 2 * (t_ef * K_cr))
                    return A_net_v
                
                match mode:
                    case "A":
                        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                        def tef():
                            t_ef = 0.4*t_1 # Épaisseur efficace
                            return t_ef
                    case "B":
                        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                        def tef():
                            t_ef = 1.4 * sqrt(M_y_Rk / (f_hk * diam)) # Épaisseur efficace
                            return t_ef
                    case "D"|"G":
                        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                        def tef():
                            t_ef = t_1 * (sqrt(2 + (4 * M_y_Rk ) / (f_hk * diam * t_1**2)-1)) # Épaisseur efficace
                            return t_ef
                    case "E"|"H":
                        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                        def tef():
                            t_ef = 2 * sqrt(M_y_Rk / (f_hk * diam)) # Épaisseur efficace
                            return t_ef
                tef_result = tef()
                a_net_v_result = a_net_v(tef_result[1])
                a_net_v_result = (tef_result[0] + a_net_v_result[0], a_net_v_result[1])
                return a_net_v_result

        list_A_net_v = []
        for mode in mode_rupture:
            a_net_v = _Anet_v(mode)
            list_A_net_v.append(a_net_v)
        
        a_net_v_result = min(list_A_net_v, key=lambda x: x[1])

        latex = latex + a_net_v_result[0]    
        A_net_v = a_net_v_result[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def f_bs_Rk():
            A_net_t = L_net_t * t_1
            F_bs_Rk = max(0.7 * A_net_v * f_v_k, 1.5 * A_net_t * f_t0_k)
            return F_bs_Rk
        
        result = f_bs_Rk()
        self.F_bs_Rk = result[1]
        return (latex + result[0], self.F_bs_Rk)

    def taux_F_bs_Ed(self, Fv_Ed: si.kN=0, FbsRk_1: si.kN=0, FbsRk_2: si.kN=0, loadtype=Barre.LOAD_TIME) -> tuple:
        """Détermine le taux de la rupture de bloc de l'assemblage

        Args:
            Fv_Ed (si.kN, optional): effort de cisaillement total sur l'assemblage en kN. Defaults to 0.
            FbsRk_1 (si.kN, optional): Résistance caractéristique en cisaillement de bloc de l'élément 1 en kN. Defaults to 0.
            FbsRk_2 (si.kN, optional): Résistance caractéristique en cisaillement de bloc de l'élément 2 en kN. Defaults to 0.
            loadtype (str, optional): type de durée de chargement.
        """
        Fv_Ed = abs(Fv_Ed) * si.kN
        Fbs_Rd_1 = self.F_Rd(abs(FbsRk_1), loadtype)[1]
        Fbs_Rd_2 = self.F_Rd(abs(FbsRk_2), loadtype)[1]
        alpha_1 = self.alpha[0]
        alpha_2 = self.alpha[1]
        n_cis = self.nCis
        if not Fbs_Rd_1.value and not Fbs_Rd_2.value:
            raise ValueError("Aucune des résistances caractéristiques en cisaillement de bloc n'est fournie.")
        else:
            if Fbs_Rd_1.value and Fbs_Rd_2.value:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    taux_rupture_bloc_1 = ((Fv_Ed / n_cis) * cos(radians(alpha_1))) / Fbs_Rd_1
                    taux_rupture_bloc_2 = Fv_Ed * cos(radians(alpha_2)) / Fbs_Rd_2
                    return taux_rupture_bloc_1, taux_rupture_bloc_2
                result = val()
                max_taux = max(result[1][0], result[1][1])
            elif Fbs_Rd_1.value:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    taux_rupture_bloc_1 = ((Fv_Ed / n_cis) * cos(radians(alpha_1))) / Fbs_Rd_1
                    return taux_rupture_bloc_1
                result = val()
                max_taux = result[1]
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    taux_rupture_bloc_2 = Fv_Ed * cos(radians(alpha_2)) / Fbs_Rd_2
                    return taux_rupture_bloc_2
                result = val()
                max_taux = result[1]
            synthese = [
                [f"Rupture de bloc de l'assemblage {self.type_assemblage} avec {self.type_organe}", max_taux, None],
            ]
            self._add_synthese_taux_travail(synthese)
            return result

    
    def taux_F_v_Ed(self, Fv_Ed: si.kN=0, Fax_Ed: si.kN=0, loadtype=Barre.LOAD_TIME) -> tuple:
        """Détermine le taux de cisaillement et éventuellement de la traction combinée de l'assemblage

        Args:
            Fv_Ed (float): effort de cisaillement total sur l'assemblage en kN
            Fax_Ed (float, optional): effort axialement à reprendre en kN. Defaults to 0.
            loadtype (str, optional): type de durée de chargement.
        """
        Fv_Rd_ass = self.F_Rd(self.Fv_Rk_ass.value*10**-3, loadtype)[1]
        Fv_Ed = abs(Fv_Ed) * si.kN
        Fax_Ed = abs(Fax_Ed) * si.kN
        if Fax_Ed:
            if self.type_organe in ("Pointe circulaire lisse", "Pointe carrée lisse", "Agrafe", "Boulon"):
                self.Fax_Rk_ass = self.FaxRk * self.nfile * self.n
                Fax_Rd_ass = self.F_Rd(self.Fax_Rk_ass.value*10**-3, loadtype)[1]
                @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    taux_combi = Fax_Ed / Fax_Rd_ass + Fv_Ed / Fv_Rd_ass
                    return taux_combi
            elif self.type_organe in ("Tirefond", "Autres pointes"):
                Fax_Rd_ass = self.F_Rd(self.Fax_Rk_ass.value*10**-3, loadtype)[1]
                @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    taux_combi = (Fax_Ed / Fax_Rd_ass)**2 + (Fv_Ed / Fv_Rd_ass)**2
                    return taux_combi
            else:
                raise f"L'organe {self.type_organe} ne peut pas être pris en compte pour une vérification de l'élément chargées à la fois axialement et latéralement"
        else:
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                taux_cisaillement = Fv_Ed / Fv_Rd_ass
                return taux_cisaillement
        result = val()
        synthese = [
            [f"Cisaillement assemblage {self.type_assemblage} avec {self.type_organe}", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result

    def taux_F_90_Ed(self, Fv_Ed: si.kN=0, loadtype=Barre.LOAD_TIME) -> tuple:
        """Détermine le taux de rupture par fendage de l'assemblage

        Args:
            Fv_Ed (float): effort de cisaillement en kN.
            loadtype (str, optional): type de durée de chargement.
        """
        F_90_Rd = self.F_Rd(self.F90_Rk.value*10**-3, loadtype)[1]
        Fv_Ed = abs(Fv_Ed) * si.kN
        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            taux_cisaillement = Fv_Ed / F_90_Rd
            return taux_cisaillement
        result = val()
        synthese = [
            [f"Fendage de la pièce bois", result[1], None],
        ]
        self._add_synthese_taux_travail(synthese)
        return result
    

