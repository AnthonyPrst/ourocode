# coding in UTF-8
# by Anthony PARISOT

############# Le but de ce fichier est de regrouper toute les fonctions d'assemblage par organe métalique dans l'EN-1995 #############
from copy import deepcopy
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.A0_Projet import Projet
from ourocode.eurocode.EC5_Element_droit import Barre, Cisaillement, Compression_perpendiculaire, Compression_inclinees

def interpolationLineaire(x, xa, xb, ya, yb):
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

    def __init__(self,beam_1:object, beam_2:object, nfile: int=1, nCis: int=("1","2"), **kwargs):
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
    def rho_mean_ass(self):
        if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[0]:
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
            print("Le sens de traitement de l'assemblage a été changé car le nef mini ce trouve sur la pièce 2 et non sur la pièce 1.\nAttention aux efforts de calcul à prendre en compte.")
            return list_nef[1]
        
        
    # 7.1 Glissement des assemblages
    def Kser(self):
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

    def Kser_ass(self):
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

    def Ku(self):
        """Retourne le module de glissement de l'organe pour les états limites ultimes."""
        K_ser = self.Kser()[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            K_u = K_ser * 2 / 3
            return K_u
        return val()

    def Ku_ass(self):
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


    def F90Rk(self, b:int, h:int, he:int, w: int=1):
        """Calcul la valeur caractérisque de la capacité au fendage en N avec :
            b : l'épaisseur de l'élément en mm
            h : la hauteur de l'élément en mm
            he : la distance de rive chargée vis à vis du centre de l'organe le plus plus éloigné ou du bord de la plaque
            w : facteur de modification """
        h_e = he
        @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_90_Rk = 14 * b * w * sqrt(h_e / (1 - (h_e / h))) # N
            return F_90_Rk * si.N
        return val()

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
            if self.nCis == 2:
                self.pos_plaq = "externe"
        else:
            f_h1k = self.fh1k[1].value * 10**-6
            t_1 = self.t1.value * 10**3
            self.t = self.beam_2.t
            if self.nCis == 2:
                self.pos_plaq = "centrale"

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

    def FvRk(self, effet_corde: bool=("True", "False")):
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


    def F_Rd(self, F_Rk: si.kN, loadtype=Barre.LOAD_TIME):
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
    def FbsRk(self, dp:si.mm, a1:si.mm, a2:si.mm, a3t:si.mm, Kcr: float=0.67, num_beam: int=("1", "2")):
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
    
    def taux_cisaillement(self, Fv_Ed: si.kN=0, Fax_Ed: si.kN=0, loadtype=Barre.LOAD_TIME):
        """Détermine le taux de cisaillement ou du chargement combiné de l'assemblage

        Args:
            Fv_Ed (float): effort de cisaillement à reprendre en kN
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
        return val()
    

class Embrevement(Assemblage):
    TYPE_EMBREVEMENT = ("Bissectrice", "Equerre à la pièce 2", "Equerre à la pièce 1")
    def __init__(
        self, 
        alpha: float=45, 
        type_embrevement: str=TYPE_EMBREVEMENT, 
        prof_embrevement: si.mm=0, 
        l_talon: si.mm=500,
        l1: si.mm=10000,
        **kwargs):
        """Créer une classe Embrevement hérité de la classe Assemblage du module EC5_Assemblage.py.
        
        Cette classe permet de vérifier un embrèvement simple avec comme i = 1 la pièce avec l'embrèvement en bout et comme i = 2 la pièce entaillée.
        Attention cette classe ne vérifie pas la section réduite de la pièce 2 !
        
        Args:
            alpha (float): angle en degrés entre le fil de la pièce 1 et le fil de la pièce 2 (Ne peux pas être supérieur à 90°). Defaults to 45.
            type_embrevement (str): type d'embrèvement "Bissectrice", "Equerre à la pièce 2", "Equerre à la pièce 1". Defaults to "Bissectrice".
            prof_embrevement (si.mm): profondeur de l'embrèvement en mm parallèle au fil de la pièce 2. Si ce dernier est égale à 0,
                alors la profondeur théorique est calculé automatiquement suivant les règles ci-dessous:
                    - 25°<= alpha < 50° : tv <= h2/4 
                    - 50°<= alpha < 60° : tv <= h2*(2/3 - alpha/120)
                    - alpha >= 60° : tv <= h2/6
                    avec h2 = hauteur de la pièce 2
            l_talon (si.mm): longueur du talon en mm sur la pièce 2. Defaults to 500.
            l1 (si.mm): distance de la face l'appuis à l'embrèvement sur la pièce 2 (l et l) (si pas de l1 ne rien mettre). Defaults to 10000.
        """
        super().__init__(**kwargs)
        self.alpha = alpha
        self.type_embrevement = type_embrevement
        self.prof_embrevement = prof_embrevement * si.mm
        self.l_talon = l_talon * si.mm
        self.l1 = l1 * si.mm
        self.dict_taux_embrevement = {}

        if self.prof_embrevement == 0*si.mm:
            h2 = self.beam_2.h_calcul
            def prof_embrevement():
                if 25 <= self.alpha <= 50:
                    prof_embrevement = h2/4
                elif 50 < self.alpha <= 60:
                    prof_embrevement = h2*(2/3 - self.alpha/120)
                elif self.alpha > 60:
                    prof_embrevement = h2/6
                return round(prof_embrevement, 0)
            self.prof_embrevement = prof_embrevement()
    
    def _compression_about(self, N_c1_d: si.kN, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Vérifie la compression d'about de l'embrevement"""
        if self.type_embrevement == "Bissectrice":
            alpha_about = self.alpha / 2
            tv = self.prof_embrevement / cos(radians(alpha_about))
        elif self.type_embrevement == "Equerre à la pièce 1":
            alpha_about = 0
            tv = self.prof_embrevement / cos(radians(self.alpha))
        elif self.type_embrevement == "Equerre à la pièce 2":
            alpha_about = self.alpha
            tv = self.prof_embrevement
        
        Fcad = N_c1_d * cos(radians(alpha_about))

        c_alpha = Compression_inclinees._from_parent_class(
            self.beam_1, 
            b_appuis=self.beam_1.b_calcul.value*10**3, 
            l_appuis=tv.value*10**3, 
            l1d=10000, 
            l1g=10000, 
            ad=10000, 
            ag=10000, 
            type_appuis_90="Appuis discret", 
            alpha=alpha_about)
        c_alpha.K_c90 = 1 # On bride le Kc90 à 1 de manière sécuritaire

        sigma_c_alpha_d = c_alpha.sigma_c_alpha_d(Fcad)
        taux_c_alpha_d = c_alpha.taux_c_alpha_d(loadtype, typecombi)
        latex = (sigma_c_alpha_d[0] + taux_c_alpha_d[0])
        self.dict_taux_embrevement["taux compression about"] = taux_c_alpha_d[1]
        return (latex, taux_c_alpha_d[1])

    def _compression_transversale(self, N_c1_d: si.kN, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Vérifie la compression transversale de l'embrevement"""
        l_appuis = self.beam_1.h_calcul / sin(radians(self.alpha))
        c_90 = Compression_perpendiculaire._from_parent_class(
            self.beam_2, 
            b_appuis=self.beam_1.b_calcul.value*10**3, 
            l_appuis=l_appuis.value*10**3, 
            l1d=self.l1.value*10**3, 
            l1g=10000, 
            ad=10000, 
            ag=10000, 
            type_appuis_90="Appuis discret")
        
        Fc90d = N_c1_d * sin(radians(self.alpha))
        f_c_90_d = c_90.f_c_90_d(loadtype, typecombi)
        sigma_c_90_d = c_90.sigma_c_90_d(Fc90d)
        taux_c_90_d = c_90.taux_c_90_d()
        latex = (f_c_90_d[0] + sigma_c_90_d[0] + taux_c_90_d[0])
        self.dict_taux_embrevement["taux compression transversale"] = taux_c_90_d[1]
        return (latex, taux_c_90_d[1])

    def _cisaillement_talon(self, N_c1_d: si.kN, loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Vérifie le cisaillement du talon de la pièce 2, le talon est limité pour prendre en compte l'effet ciseau à bois"""
        N_c1_d = N_c1_d * si.kN
        l_talon = self.l_talon
        prof_embr = self.prof_embrevement
        b_embr = self.beam_1.b_calcul
        k_cr = 0.67 # On bride le Kcr à 0.67 de manière sécuritaire
        f_v_d = self.beam_2._f_type_d("fvk", loadtype, typecombi)[1]
        alpha = self.alpha

        if self.beam_1.b_calcul < self.beam_2.b_calcul:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                l_v_ef = min(l_talon, prof_embr * 8)
                L_v = (b_embr + 2 * prof_embr) * k_cr # périmètre transversal cisaillé
                N_c2_d = N_c1_d * cos(radians(alpha))
                taux_talon = N_c2_d / (f_v_d * L_v * l_v_ef)
                return taux_talon
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                l_v_ef = min(l_talon, prof_embr * 8)
                L_v = b_embr * k_cr # périmètre transversal cisaillé
                N_c2_d = N_c1_d * cos(radians(alpha))
                taux_talon = N_c2_d / (f_v_d * L_v * l_v_ef)
                return taux_talon
        value = val()
        self.dict_taux_embrevement["taux cisaillement talon"] = value[1]
        return value

    def taux_embrevement(self, N_c1_d: si.kN, loadtype: str=Barre.LOAD_TIME, typecombi: str=Barre.TYPE_ACTION):
        """Retourne les taux de travail de l'embrèvement avec:
            - la compression d'about de la pièce 1
            - la compression transversale de la pièce 2
            - le cisaillement du talon de la pièce 2

        Args:
            N_c1_d (si.kN): effort normal de la pièce 1 en kN.
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.
        """
        about =self._compression_about(N_c1_d, loadtype, typecombi)
        transversale = self._compression_transversale(N_c1_d, loadtype, typecombi)
        talon = self._cisaillement_talon(N_c1_d, loadtype, typecombi)
        return (about[0] + transversale[0] + talon[0], self.dict_taux_embrevement)



# ======================================================= POINTE =========================================================
# 8.3 Assemblage par pointes

class Pointe(Assemblage):
    QUALITE_ACIER = ('6.8', '8.8', '9.8', '10.9', '12.9')
    TYPE_ORGANE = ("Pointe circulaire lisse", "Pointe carrée lisse", "Autres pointes")
    def __init__(self, d:si.mm, dh:si.mm, l:si.mm, qualite: str=QUALITE_ACIER, n: int=1, alpha1: float=0, alpha2: float=0, type_organe: str=TYPE_ORGANE, percage: bool=("False", "True"), *args, **kwargs):
        """ 
        Créer une classe Pointe hérité de la classe Assemblage du module EC5_Assemblage.py.
        
        Args:
            d (float): diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            dh (float): diamètre de la tête en mm
            l (int): longueur sous la tête en mm
            qualite (str): qualité de l'acier
            n (int): nombre d'organe dans une file 
            alpha1 (float, optional): angle entre l'effort de l'organe et le fil du bois 1 en °
            alpha2 (float, optional): angle entre l'effort de l'organe et le fil du bois 2 en °
            type_organe (str): type de l'organe "Pointe circulaire lisse", "Pointe carrée lisse", "Autres pointes" pour les pointes torsadées, annelées, crantées
            percage (bool): Si il y a un prépercage de la pointe alors True sinon False. Defaults to False".
        """
        super().__init__(*args, **kwargs)
        self.d = d * si.mm
        self.dh = dh * si.mm
        self.l = l * si.mm #longueur sous tête
        self.qualite = qualite
        self.n = n
        self._nef = n
        self.fu = self.__qualite_acier.loc["fub"] *si.MPa
        self.alpha = [alpha1, alpha2]
        self.percage = percage
        self.type_organe = type_organe
        self.__t1_t2()
        self._caracteristique_min()
        self._fhik()
    
    def _caracteristique_min(self):
        """ Vérifie si les caractéristiques minimales sont respectées """
        if self.fu < 600*si.MPa:
            raise ValueError("La résistance du fil en traction est inférieur à 600 MPa, vérifier vos données !")
        if not self.percage and self.d > 6*si.mm:
            raise ValueError(f"Erreur, le diamètre de la pointe est supérieur à 6mm, le prépercage est obligatoire")
        if self.type_assemblage == self.TYPE_ASSEMBLAGE[0]: #Si assemblage bois bois
            if self.type_organe in ("Pointe circulaire lisse", "Pointe carrée lisse") and self.t2 < 8*self.d:
                raise ValueError(f"Erreur, la longueur de pénétration t2 est inférieur à 8 fois le diamètre de la pointe. La longueur de pénétration minimal est de {self.d*8}")
            elif self.type_organe in ("Autres pointes", "Tirefond") and self.t2 < 6*self.d:
                raise ValueError(f"Erreur, la longueur de pénétration t2 est inférieur à 6 fois le diamètre de la pointe. La longueur de pénétration minimal est de {self.d*6}")


    @property
    def __qualite_acier(self):
        df = self._data_from_csv("qualite_acier.csv")
        df = df.loc[self.qualite]
        return df

    @property
    def _type_circulaire(self):
        if self.type_organe == "Pointe carrée lisse":
            return False
        else:
            return True


    def __t1_t2(self):
        """Retourne t1 et t2 en mm suivant l'EN 1995 §8.3.1.1
        """
        l_pointe = self.d #Longueur de la pointe à déduire
        if self.type_organe == "Tirefond":
            l_pointe = self.d_vis
            
        if self.nCis == 1:
            if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t1 = self.beam_1.b_calcul
            else:
                self.t1 = self.beam_1.t

            if self._type_beam[1] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t2 = self.l - self.t1 - l_pointe
            else:
                raise ValueError("Il n'est pas considéré possible d'avoir un assemblage simple cisaillement avec flasque métalique en position 2")

        else:
            if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
                if "Bois/Métal" in self.type_assemblage:
                    b_beam_2 = self.beam_2.t + 2
                else:
                    b_beam_2 = self.beam_2.b_calcul
                    self.t1 = min(self.beam_1.b_calcul, self.l - self.beam_1.b_calcul - b_beam_2 - l_pointe)
            else:
                raise ValueError("Il n'est pas considéré possible d'avoir un assemblage double cisaillement avec flasque métalique et un organe de type pointe!")
            if self._type_beam[1] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t2 = self.beam_2.b_calcul
            else:
                self.t2 = self.beam_2.t


    @property
    def MyRk(self):
        """ Défini le moment d'écoulement plastique d'une pointe en N.mm avec:
            d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            fu : la résistance caractéristique en traction du fil d'acier en N/mm2 """
        f_u = self.fu.value * 10**-6
        d = self.d.value * 10**3
        if self.fu >= 600:
            if self._type_circulaire == True:
                @handcalc(override="short", precision=2, left="\\[", right="\\]")
                def val():
                    M_y_Rk = 0.3 * f_u* d**2.6 # N.mm
                    return M_y_Rk * si.N*si.mm
            else:
                @handcalc(override="short", precision=2, left="\\[", right="\\]")
                def val():
                    M_y_Rk = 0.45 * f_u* d**2.6 # N.mm
                    return M_y_Rk * si.N*si.mm
            return val()
        else:
            raise ValueError("La résistance du fil en traction est inférieur à 600 MPa, vérifier vos données !")


    def _fhk_bois(self, beam:object):
        """ Calcul la portance locale des pointes inférieur à 8mm dans le bois et le LVL en MPa
            """
        d = self.d.value * 10**3
        rho_k = beam.rho_k
        if self.percage:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                fhk = 0.082 * (1 - 0.01 * d) * rho_k # MPa
                return fhk * si.MPa
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                fhk = 0.082 * rho_k * d**(-0.3) # MPa
                return fhk * si.MPa
        return val()


    def _fhk_panneau(self, beam:object):
        """Calcul la portance locale des pointes dans les panneaux en MPa

        Args:
            t (int): épaisseur du panneau en mm
            self.dh (int): diamètre de la tête de la pointe

        Returns:
            float: portance locale en MPa
        """
        d = self.d.value * 10**3
        b_calcul = beam.b_calcul.value * 10**3
        rho_k = beam.rho_k
        if self.dh >= 2*self.d or self.type_organe == "Agrafe":
            if beam.type_bois == "CP":
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    fhk = 0.11 * rho_k * d**-0.3 # MPa
                    return fhk * si.MPa
            elif beam.type_bois == "Panneau dur":
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    fhk = 30 * d**-0.3 * b_calcul**0.6 # MPa
                    return fhk * si.MPa
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    fhk = 65 * d**-0.7 * b_calcul**0.1 # MPa
                    return fhk * si.MPa
            return val()
        else:
            raise ValueError(f"La tête de la pointe doit être au moins égale à {2*d} mm selon l'EN 1995-1-1 §8.3.1.3 (3).")


    def _fhik(self) -> tuple:
        """Calcul la portance locale d'une pointe en MPa dans les deux éléments de l'assemblage

        ATTENTION: Pas de prise en compte des panneaux durs au sens de la norme EN 622-2

        Returns:
            tuple: (fh1k, fh2k) en MPa
        """
        dict_beam = {"1": {}, "2": {}}
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if self._type_beam[i] == "Bois":
                dict_beam[str(i+1)]["fhk"] = self._fhk_bois(beam)
                
            elif self._type_beam[i] == "CP" or self._type_beam[i] == "PP/OSB":
                dict_beam[str(i+1)]["fhk"] = self._fhk_panneau(beam)
                
            else:
                dict_beam[str(i+1)]["fhk"] = 0
        
            # if self._type_beam[i] != "Métal":
            #     if i:
            #         self.t2 = self.beam_2.b_calcul
            #     else:
            #         self.t1 = self.beam_1.b_calcul
        self.fh1k = dict_beam["1"]["fhk"]
        self.fh2k = dict_beam["2"]["fhk"]
        return self.fh1k, self.fh2k
    
    
    def Fax_Rk(self):
        """Calcul la valeur caractéristique de résistance axial entre la résistance caractéristique de la tête et du corps de la tige dans le bois en N.
        """
        if self._type_beam[1] == "Métal":
            rho_k_ax  = self.beam_1.rho_k
        else:
            rho_k_ax  = self.beam_2.rho_k
        rho_k_head  = self.beam_1.rho_k
        d = self.d
        d_h = self.dh
        t = self.t1
        t_pen = self.t2
        
        if self.type_organe == "Autres pointes":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_ax_k = 20 * 10**-6 * rho_k_ax**2 * si.MPa
                F_ax_a_Rk = f_ax_k * d * t_pen
                f_head_k = 70 * 10**-6 * rho_k_head**2 * si.MPa
                F_head_Rk = f_head_k * d_h**2
                F_ax_Rk = min(F_ax_a_Rk, F_head_Rk)
                return F_ax_Rk
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_ax_k = 20 * 10**-6 * rho_k_ax**2 * si.MPa
                F_ax_a_Rk = f_ax_k * d * t_pen
                f_head_k = 70 * 10**-6 * rho_k_head**2 * si.MPa
                F_head_Rk = f_ax_k * d * t + f_head_k * d_h**2
                F_ax_Rk = min(F_ax_a_Rk, F_head_Rk)
                return F_ax_Rk
        
        F_ax_Rk = val()
        if self.type_organe == "Autres pointes" and self.t2 < 8 * self.d:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_ax_Rk_min = F_ax_Rk * (t_pen / (2 * d) - 3)
                return F_ax_Rk_min
            minoration = val2()
            self.FaxRk = minoration[1]
            return (F_ax_Rk[0] + minoration[0], minoration[1])

        elif self.type_organe and self.t2 < 12 * self.d:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val2():
                F_ax_Rk_min = F_ax_Rk * (t_pen / (4 * d) - 2)
                return F_ax_Rk_min
            minoration = val2()
            self.FaxRk = minoration[1]
            return (F_ax_Rk[0] + minoration[0], minoration[1])
            
        else:
            self.FaxRk = F_ax_Rk[1]
            return F_ax_Rk


    def _kef(self, a1:int):
        """ coefficient donnée dans le tableau 8.1 fonction de a1 et du percage qui réduit le nombre efficace de pointe dans le sens 
            du fil avec :
                a1 : l'espacement entre tige dans le sens du fil du bois """
        if self.type_organe == "Tirefond":
            d = self.d_vis
        else:
            d = self.d
        listeTab = (d * 4, d * 7, d *
                    10, d * 14, ("x", 0.7, 0.85, 1))

        if a1 >= listeTab[3]:
            kef = 1
        j = 0
        for i in range(4):

            if listeTab[i] <= a1:

                if a1 == listeTab[i]:
                    kef = listeTab[4][i]

                    if kef == "x" and self.percage == False:
                        kef = 0

                    elif kef == "x" and self.percage == True:
                        kef = 0.5
                j += 1

            else:
                if listeTab[4][j-1] == "x" and self.percage == False:
                    kef = interpolationLineaire(
                        a1, listeTab[j-1], listeTab[j], 0, listeTab[4][j])

                elif listeTab[4][j-1] == "x" and self.percage == True:
                    kef = interpolationLineaire(
                        a1, listeTab[j-1], listeTab[j], 0.5, listeTab[4][j])

                else:
                    kef = interpolationLineaire(
                        a1, listeTab[j-1], listeTab[j], listeTab[4][j-1], listeTab[4][j])
        self.kef = kef
        return self.kef


    def nef(self, a1_beam1:int, a1_beam2:int):
        """Défini le nombre efficace d'organe dans une file avec :
            a1_beam1 : espacement entre les organes dans la barre 1 en mm
            a1_beam2 : espacement entre les organes dans la barre 2 en mm"""
        nef_list = []
        for i, _ in enumerate([self.beam_1, self.beam_2]):
            if i:
                a_1 = a1_beam2 * si.mm
                n = self.nfile
            else:
                a_1 = a1_beam1 * si.mm
                n = self.n
                
            k_ef = self._kef(a_1)
        
            if n == 1:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    n_ef = 1
                    return n_ef
                nef_list.append(val())
            else :
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    n_ef = n**k_ef
                    return n_ef
                nef_list.append(val())
        result = self._min_nef(nef_list)
        self._nef = result[1]
        return result

    @property
    def pince(self):
        """
        Défini les différentes pinces minimales pour une pointe en mm.

        Args:
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace de la pointe ou du tire fond si d<=6mm en mm
        """
        dict_pince = {}
        if self.type_organe == "Tirefond":
            self.d = self.d_vis
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if not self._type_beam[i] in self.TYPE_BOIS_ASSEMBLAGE:
                continue
            
            alpha = self.alpha[i]
            rho_k = beam.rho_k
            if self.percage == True:
                a1 = round((4 + mt.cos(mt.radians(alpha))) * self.d, 1)
                a2 = round((3 + mt.sin(mt.radians(alpha))) * self.d, 1)
                a3t = round((7 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
                a3c = round(7 * self.d, 1)
                a4c = round(3 * self.d, 1)

                if self.d < 5:
                    a4t = round((3 + 2 * mt.sin(mt.radians(alpha))) * self.d, 1)
                else:
                    a4t = round((3 + 4 * mt.sin(mt.radians(alpha))) * self.d, 1)
            else:
                if rho_k <= 420:
                    if self.d < 5:
                        print (mt.cos(mt.radians(alpha)))
                        a1 = round((5 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
                        a4t = round((5 + 2 * mt.sin(mt.radians(alpha))) * self.d, 1)
                    else:
                        a1 = round((5 + 7 * mt.cos(mt.radians(alpha))) * self.d, 1)
                        a4t = round((5 + 5 * mt.sin(mt.radians(alpha))) * self.d, 1)

                    a2 = round(5 * self.d, 1)
                    a3t = round((10 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
                    a3c = round(10 * self.d, 1)
                    a4c = a2

                elif (rho_k > 420 and rho_k <= 500) or self._type_beam[i] in self.TYPE_BOIS_ASSEMBLAGE[1:]:
                    a1 = round((7 + 8 * mt.cos(mt.radians(alpha))) * self.d, 1)
                    a2 = round(7 * self.d, 1)
                    a3t = round((15 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
                    a3c = round(15 * self.d, 1)
                    a4c = a2

                    if self.d < 5:
                        a4t = round((7 + 2 * mt.sin(mt.radians(alpha))) * self.d, 1)
                    else:
                        a4t = round((7 + 5 * mt.sin(mt.radians(alpha))) * self.d, 1)
                else:
                    raise ValueError("Il faut absolument prépercer au dessus de rho,k: 500 kg/m3")

            if self._type_beam[i] in self.TYPE_BOIS_ASSEMBLAGE[1:]: #Si assemblage bois panneau
                a1 = round(a1 * 0.85, 1)
                a2 = round(a2 * 0.85, 1)
                if beam.type_bois == "CP":
                    a3t = round((3 + 4 * mt.sin(mt.radians(alpha))) * self.d, 1)
                    a3c = round(3 * self.d, 1)
                    a4t = a3t
                    a4c = a3c

            elif self.type_assemblage == "Bois/Bois": #Si assemblage bois bois
                pass

            elif self.type_assemblage == "Bois/Métal": #Si assemblage bois métal
                a1 = round(a1 * 0.7, 1)
                a2 = round(a2 * 0.7, 1)

            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince


    def prepercage(self, beam: str=["1", "2"], sensible: bool=("False", "True")):
        """Retourne l'épaisseur en mm minimale pour éviter le pré-perçage des pointes

        Args:
            beam (str, optional): Défini la barre à calculer entre 1 et 2 selon l'EN1995. Defaults to "1".
            sensible (bool, optional): Défini si le bois utilisé est sensible à la fissuration (selon AN §8.3.1.2(7) douglas et pin maritime). Defaults to False.

        Returns:
            int: l'épaisseur minimale du bois en mm
        """
        d = self.d
        if beam == "1":
            rho_k  = self.beam_1.rho_k
        else:
            rho_k  = self.beam_2.rho_k
        if sensible:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_min = max(14 * d, (13 * d - 30) * (rho_k / 200))
                return t_min
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                t_min = max(7 * d, (13 * d - 30) * (rho_k / 400))
                return t_min
        return val() 
        

# ======================================================= Agrafe =========================================================
# 8.4 Assemblage par Agrafe

class Agrafe(Pointe):
    TYPE_ASSEMBLAGE = ("Bois/Bois", ("CP", "Panneau dur", "PP/OSB"), "Bois/Métal")
    QUALITE_ACIER = ('8.8', '9.8', '10.9', '12.9')
    def __init__(self, d:si.mm, b_agrafe:si.mm, l:si.mm, qualite: str=QUALITE_ACIER, n: int=1, angle_sup_30: bool=["True", "False"], alpha1: float=0, alpha2: float=0, **kwargs):
        """
        Créer une classe Agrafe hérité de la classe Assemblage du module EC5_Assemblage.py.
        
        Args:
            d (float): diamètre de l'agrafe en mm, si l'agrafe est de section rectangulaire alors c'est la racine carrée du produit des 2 dimensions selon EN1995 §8.4(2).
            b_agrafe (float): dimension du dos de l'agrafe en mm.
            l (int): longueur sous la tête en mm
            qualite (str): qualité de l'acier
            n (int): nombre d'agrafe dans une file 
            angle_sup_30 (str): Si l'angle entre la tête de l'agrafe et le fil du bois est supérieur à 30° alors True sinon False. Defaults to True.
            alpha1 (float, optional): angle entre l'effort de l'agrafe et le fil du bois 1 en °. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'agrafe et le fil du bois 2 en °. Defaults to 0.
        """
        super().__init__(d=d, dh=0, l=l, qualite=qualite, n=n, alpha1=alpha1, alpha2=alpha2, type_organe="Agrafe", percage=False, **kwargs)
        self.b_agrafe = b_agrafe * si.mm
        self.angle_sup_30 = angle_sup_30
        self._dimension_min()

    def _dimension_min(self):
        """ Vérifie si les dimensions minimales sont respectées """
        if self.b_agrafe < self.d*6:
            raise ValueError(f"Erreur, la dimension du dos de l'agrafe est inférieur à 6 fois le diamètre de l'agrafe, le dos de l'agrafe minimal est de {self.d*6}")
        if self.t2 < self.d*14:
            raise ValueError(f"Erreur, la longueur de pénétration t2 est inférieur à 14 fois le diamètre de l'agrafe. La longueur d'agrafe minimal est de {self.t1 + self.d*14}")
        
    @property
    def MyRk(self):
        """ Défini le moment d'écoulement plastique d'une pointe en N.mm avec"""
        d = self.d.value * 10**3
        @handcalc(override="short", precision=2, left="\\[", right="\\]")
        def val():
            M_y_Rk = 150 * d**3 # N.mm
            return M_y_Rk * si.N*si.mm
        return val()

    def nef(self):
        """Retourne le nombre efficace d'organe dans une file"""
        self._nef = self.n
        return self._nef

    @property
    def pince(self):
        """
        Défini les différentes pinces minimales pour une pointe en mm.

        Args:
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace de la pointe ou du tire fond si d<=6mm en mm
        """
        dict_pince = {}
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if not self._type_beam[i] in self.TYPE_BOIS_ASSEMBLAGE:
                continue
            alpha = self.alpha[i]
            if not self.angle_sup_30:
                a1 = round((10 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
            else:
                a1 = round((15 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
            a2 = round(15 * self.d, 1)
            a3t = round((15 + 5 * mt.cos(mt.radians(alpha))) * self.d, 1)
            a3c = round(15 * self.d, 1)
            a4t = round((15 + 5 * mt.sin(mt.radians(alpha))) * self.d, 1)
            a4c = round(10 * self.d, 1)
            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince

# ======================================================= BOULON =========================================================
# 8.5 Assemblage par boulon

class Boulon(Assemblage):
    QUALITE_ACIER = tuple(Assemblage._data_from_csv(Assemblage, "qualite_acier.csv").index)
    def __init__(self, d:si.mm, qualite: str=QUALITE_ACIER, n: int=1, alpha1: float=0, alpha2: float=0, t1: int=0, t2: int=0, **kwargs):
        """
        Créer une classe Boulon hérité de la classe Assemblage du module EC5_Assemblage.py.
        
        Args:
            d (int): diamètre efficace du boulon (ou du tire fond si >6mm) en mm
            qualite (str): qualité de l'acier
            n (int): nombre de boulons dans une file
            alpha1 (float, optional): angle entre l'effort de l'organe et le fil du bois en ° pour la barre 1. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'organe et le fil du bois en ° pour la barre 2. Defaults to 0.
            t1 (int, optional): longueur de contact avec la tige  pour la pièce 1 en mm. 
                ATTENTION : Cet argument n'est pas obligatoire par défaut, il est calculé par le type de tige utilisée.
                Il n'est nécessaire de le remplir uniquement si vous avez un t1 spécifique, par exemple avec une chapelle réduisant ainsi la portance local à une longueur inférieur à celle de l'épaisseur de la pièce 1.
            t2 (int, optional): longueur de contact avec la tige  pour la pièce 2 en mm. 
                ATTENTION : Même chose que pour t1 mais pour la pièce 2.
        """
        super().__init__(**kwargs)
        self.type_organe = "Boulon"
        if "type_organe" in kwargs.keys():
            self.type_organe = kwargs.pop("type_organe")
        self.d = d * si.mm
        self.qualite = qualite
        self.fuk = self.__qualite_acier.loc["fub"] *si.MPa
        self.n = n
        self._nef = n
        self.alpha = [alpha1, alpha2]
        self.t1 = t1
        self.t2 = t2

        self.__t1_t2()
        self._fhik()
    
    @property
    def __qualite_acier(self):
        df = self._data_from_csv("qualite_acier.csv")
        df = df.loc[self.qualite]
        return df

    # 8.5.1 Boulons chargés latéralement
    # 8.5.1.1 Généralité et assemblage bois/bois
    
    def __t1_t2(self):
        """Retourne t1 et t2 en mm suivant l'EN 1995 §8.3.1.1
        """
        for i, _ in enumerate([self.beam_1, self.beam_2]):
            if self._type_beam[i] != "Métal":
                if not i:
                    if not self.t1:
                        if self.type_organe == "Tirefond" and self.nCis == 2:
                            if "Bois/Métal" in self.type_assemblage:
                                b_beam_2 = self.beam_2.t + 2
                            else:
                                b_beam_2 = self.beam_2.b_calcul
                            l_pointe = self.d_vis
                            print(self.l , self.beam_1.b_calcul , b_beam_2 , l_pointe)
                            self.t1 = min(self.beam_1.b_calcul, self.l - self.beam_1.b_calcul - b_beam_2 - l_pointe)
                            continue
                        self.t1 = self.beam_1.b_calcul
                    else:
                        self.t1 = self.t1 * si.mm
                else:
                    l_pointe = 0
                    if self.type_organe == "Tirefond":
                        l_pointe = self.d_vis
                        l_vis = self.l  
                        if i:
                            if not self.t2:
                                self.t2 = min(self.beam_2.b_calcul,l_vis-l_pointe-self.t1)
                            else:
                                self.t2 = self.t2 * si.mm
                    else:
                        if i:
                            if not self.t2:
                                self.t2 = self.beam_2.b_calcul
                            else:
                                self.t2 = self.t2 * si.mm
                                
    
    def _fh0k(self, beam: object):
        """Calcul la portance locale d'un boulon bois/bois ou d'un tire fond si d>6mm avec :
            beam: poutre à calculer
        """
        rho_k = beam.rho_k
        d = self.d.value *10**3
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_h0k = 0.082 * (1 - 0.01 * d) * rho_k # MPa
            return f_h0k * si.MPa
        return val()
    
    
    def _K_90(self, beam: object):
        """Coef. modélisant la diminution de portance local quand un angle est donnée entre l'effort et le fil avec
            beam: poutre à calculer"""
            
        if beam.classe[0:1] == "C" or  beam.classe[0:2] == "GL":
                type_b = "C"
        elif beam.classe[0:3] == "LVL":
            type_b = "LVL"
        else: 
            type_b = "D"
                    
        if type_b == "C":
            ck90 = 1.35 + 0.015 * self.d.value *10**3
        elif type_b == "LVL":
            ck90 = 1.30 + 0.015 * self.d.value *10**3
        else:
            ck90 = 0.9 + 0.015 * self.d.value *10**3
        return ck90
    

    def _fhak(self, fh0k:float, k_90:float, alpha):
        """Calcul la portance locale d'un boulon bois/bois ou d'un tire fond si d>6mm par rapport à un effort donné à un angle
        du fil en MPa avec :
            fh0k : portance locale dans le sens du fil d'un boulon
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            k90 : coef. de réduction de la portance locale quand un effort à un angle par rapport au fil du bois"""
        f_h0k = fh0k[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_h_alpha_k = f_h0k / (k_90 * sin(radians(alpha)) ** 2 + cos(radians(alpha)) ** 2)
            return f_h_alpha_k
        return val()
    
    
    # 8.5.1.1 Généralité et assemblage bois/panneaux
    def _fhk(self, beam:object):
        """ Calcul la portance locale d'un boulon dans un assemblage bois/panneaux en MPa avec :
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm
            rho_k : masse volumique caractéristique du contreplaqué en kg/m3
            ep : epaisseur du panneau en mm """
        rho_k = beam.rho_k
        d = self.d.value *10**3
        ep = beam.b_calcul.value *10**3
        
        if beam.type_bois == "CP":
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_hk = 0.11 * (1 - 0.01 * d) * rho_k # MPa
                return f_hk
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_hk = 50 * (d**(-0.6)) * (ep**0.2) # MPa
                return f_hk
        return val()
    
    def _fhik(self):
        """Calcul la portance locale d'un boulon bois/bois ou d'un tire fond si d>6mm
        """
    
        dict_beam = {"1": {}, "2": {}}
            
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if self._type_beam[i] == "Bois":
                dict_beam[str(i+1)]["fh0k"] = self._fh0k(beam)
                dict_beam[str(i+1)]["K90"] = self._K_90(beam)
                dict_beam[str(i+1)]["fhik"] = self._fhak(dict_beam[str(i+1)]["fh0k"], dict_beam[str(i+1)]["K90"], self.alpha[i])
                dict_beam[str(i+1)]["fhik"] = (dict_beam[str(i+1)]["fh0k"][0] + dict_beam[str(i+1)]["fhik"][0], dict_beam[str(i+1)]["fhik"][1])
                
            elif self._type_beam[i] == "CP" or self._type_beam[i] == "PP/OSB":
                dict_beam[str(i+1)]["fhik"] = self._fhk(beam)
                
            else:
                dict_beam[str(i+1)]["fhik"] = 0
        
        self.fh1k = dict_beam["1"]["fhik"]
        self.fh2k = dict_beam["2"]["fhik"]
        return self.fh1k, self.fh2k


    @property
    def MyRk(self):
        """Défini le moment d'écoulement plastique d'un boulon en N.mm avec:
            fuk : la valeur caractéristique de résistance à la traction du boulon en N/mm2
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm"""
        f_uk = self.fuk.value * 10**-6
        d = self.d.value * 10**3
        @handcalc(override="short", precision=2, left="\\[", right="\\]")
        def val():
            M_y_Rk = 0.3 * f_uk * d ** 2.6 # N.mm
            return M_y_Rk * si.N*si.mm
        return val()


    @property
    def pince(self):
        """Défini les différentes pinces minimales pour un boulon en mm avec :
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm """
        dict_pince = {}
        if self.type_organe == "Tirefond":
            self.d = self.d_vis
        for i, alpha in enumerate(self.alpha):
            if not self._type_beam[i] in self.TYPE_BOIS_ASSEMBLAGE:
                continue
            a1 = round((4 + mt.cos(mt.radians(alpha))) * self.d, 1)
            a2 = round(4 * self.d, 1)
            a3t = round(max(7 * self.d, 80*si.mm),1)

            if alpha <= 30:
                a3c = round(4 * self.d, 1)
            else:
                a3c = round((1 + 6 * mt.sin(mt.radians(alpha))) * self.d, 1)

            a4t = round(max((2 + 2 * mt.sin(mt.radians(alpha))) * self.d, 3 * self.d), 1)
            a4c = round(3 * self.d, 1)
            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince
    
        
    def nef(self, a1_beam1:int, a1_beam2:int):
        """Défini le nombre efficace d'organe (boulon) dans une file avec :
            a1 : l'espacement entre boulon dans le sens du fil du bois
            n : nombre de boulons dans une file 
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en mm"""
        d = self.d
        if self.type_organe == "Tirefond":
            d = self.d_vis
            
        nef_list = []
        for i, alpha in enumerate(self.alpha):
            if i:
                a_1 = a1_beam2 * si.mm
                n = self.nfile
            else:
                a_1 = a1_beam1 * si.mm
                n = self.n
                
            if n == 1 or alpha == 90:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    n_ef = n
                    return n_ef
                nef_list.append(val())
            else:
                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    n_ef = min(n**0.9 * (a_1/(13 * d))**(1/4), n)
                    return n_ef
                nef_list.append(val())
            
                if alpha > 0 and alpha < 90: 
                    n_ef = nef_list[i][1]
                    @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                    def nef_a():
                        n_ef_a = n_ef + alpha * (n - n_ef)/90
                        return n_ef_a
                    value2 = nef_a()
                    nef_list[i] = (nef_list[i][0] + value2[0] , value2[1])
        result = self._min_nef(nef_list)
        self._nef = result[1]
        return result


    # 8.5.2 Boulons chargés axialement
    def Fax_Rk(self, d_int: float=0, d_ext: float=0, filetage_EN1090: bool=("True", "False")):
        """Calcul la résistance axial caractéristique d'un boulon chargé axialement à partir soit de la rondelle soit de la plaque métalique.

        Args:
            d_int (float, optional): diamètre intérieur de la rondelle en mm ou du trou de perçage dans la plaque métallique. Defaults to 0.
            d_ext (float, optional): diamètre extérieur de la rondelle en mm. Defaults to 0.
            filetage_EN1090 (bool, optional): définit si le filetage est conforme à l'EN 1090, soit matricé. Si filetage usiné alors False. Defaults to True.

        Returns:
            FaxRk: la résistance axial d'un boulon en N
        """
        from ourocode.eurocode.EC3_Assemblage import Tige
        d_int = d_int * si.mm
        d_ext = d_ext * si.mm

        if self._type_beam[0] == "Métal":
            FtRd = Tige(self.d.value*10**3, d_int, self.qualite, True, filetage_EN1090, t=self.beam_1.t.value*10**3, h=self.beam_1.h.value*10**3, classe_acier=self.beam_1.classe_acier, classe_transv=self.beam_1.classe_transv).FtRd
            fc_90_k = float(self.beam_2.caract_meca.loc["fc90k"]) * si.MPa
            if self.nCis == 1:
                d_ext = min(self.beam_1.t*12, 4*self.d, d_ext)
            else:
                d_ext = min(self.beam_1.t*12, 4*self.d)
        else:
            FtRd = Tige(self.d.value*10**3, d_int, self.qualite, True, filetage_EN1090, t=0, h=0, classe_acier="S235", classe_transv=3).FtRd
            fc_90_k = float(self.beam_1.caract_meca.loc["fc90k"]) * si.MPa

        Ft_Rd = FtRd[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            A_int = pi * (d_int / 2)**2 
            A_rondelle = pi * (d_ext / 2)**2 - A_int

            f_c90_k_rond = fc_90_k * 3
            F_c90_d = f_c90_k_rond * A_rondelle
            F_ax_Rk = min(F_c90_d, Ft_Rd)
            return F_ax_Rk
        FaxRk = val()
        self.FaxRk = FaxRk[1]
        return (FtRd[0] + FaxRk[0], FaxRk[1])
        

# ======================================================= BROCHE =========================================================
# 8.6 Assemblage par broche

class Broche(Boulon):
    def __init__(self, d:float, qualite: str=Boulon.QUALITE_ACIER, n: int=1, alpha1: float=0, alpha2: float=0, t1: int=0, t2: int=0, **kwargs):
        """
        Créer une classe Broche hérité de la classe Assemblage du module EC5_Assemblage.py.

        Args:
            d (int): diamètre efficace de la broche ( entre 6 et 30 mm) en mm
            qualite (str): qualite de la broche
            n (int): nombre de broche dans une file
            alpha (float): angle entre l'effort de l'organe et le fil du bois en °
            t1 (int, optional): longueur de contacte avec la tige  pour la pièce 1 en mm. 
                ATTENTION : Cet argument n'est pas obligatoire par défaut, il est calculer par le type de tige utilisée.
                Il n'est nécessaire de le remplir que si vous avez un t1 spécifique, par exemple avec une chapelle réduisant ainsi la portance local à une longueur inférieur à celle de l'épaisseur de la pièce 1.
            t2 (int, optional): longueur de contacte avec la tige  pour la pièce 2 en mm.
                ATTENTION : Même chose que pour t1 mais pour la pièce 2.
        """
        super().__init__(d, qualite, n, alpha1, alpha2, **kwargs)
        self.type_organe = "Broche"
        self.FaxRk = 0

    @property
    def Fax_Rk(self):
        return self.FaxRk
        
    @property
    def pince(self):
        """Défini les différentes pinces minimales pour une broche en mm avec :
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm """
        dict_pince = {}
        for i, alpha in enumerate(self.alpha):
            a1 = round((3 + 2 * (mt.cos(mt.radians(alpha)))) * self.d, 1)
            a2 = round(3 * self.d, 1)
            a3t = round(max(7 * self.d, 80*si.mm), 1)

            if alpha <= 150 and alpha < 210:

                a3c = round(max(3.5 * self.d, 40*si.mm), 1)

            else:

                a3c = round((a3t * mt.sin(mt.radians(alpha))), 1)

            a4t = round(max((2 + 2 * mt.sin(mt.radians(alpha))) * self.d, 3 * self.d), 1)
            a4c = round(3 * self.d, 1)
            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince


# ======================================================= TIREFOND =========================================================


class _Tirefond(object):
    def __init__(self, d:si.mm, d1:si.mm, ds:si.mm, dh:si.mm, l:si.mm, n: int, rho_a:float, fhead:float, ftensk:float, MyRk:float=0, alpha1: float=0, alpha2: float=0, percage: bool=("False", "True"), **kwargs):
        """Défini un object tirefond

        Args:
            d (si.mm): diamètre extérieur du filet en mm
            d1 (float): diamètre du noyaux en mm
            ds (float): diamètre de la tige en mm
            dh (float): diamètre de la tête en mm
            l (float): longueur total de la vis en mm
            n (int): le nombre de vis dans une file
            rho_a (float): masse volumique associée au tirefond en fax,k en kg/m3
            fhead (float): valeur caractéristique de traversée de la tête du tirefond à l'EN 14592 en Mpa
            ftensk (float): valeur caractéristique en traction du tirefond en N

            MyRk (float, optional): le moment d'écoulement plastique de la vis en N.mm. 
                Si cette attribut est remplie, alors on récupère la valeur founis sinon on le calcul à l'EC5. 
                Defaults to 0.

            alpha1 (float, optional): angle entre l'effort de l'organe et le fil du bois 1 en °. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'organe et le fil du bois 2 en °. Defaults to 0.
            percage (bool, optional): l'élément est-il préperçé ? Si oui alors True. Defaults to ("False", "True").
        """
        self.d_vis = d * si.mm
        self.d1 = d1 * si.mm
        self.ds = ds * si.mm
        self.dh = dh * si.mm
        # self.d_ef = min(d1*1.1, ds) A valider !!!!!
        self.d_ef = d1*1.1
        self.l = l * si.mm #longueur sous tête
        self.n = n
        self.percage = percage
        self.alpha = [alpha1, alpha2]
        self.type_organe = "Tirefond"
        self.rho_a = rho_a * si.kg/si.m**3
        self.fhead = fhead * si.MPa
        self.ftensk = ftensk * si.N
    

    # 8.7.2 Tirefond chargés axialement
    def pince_tirefond_axial(self, t: int):
        """
        Défini les pinces d'un tirefond en mm lorsqu'il est chargée axialement et l'epaisseur de bois supérieur à 12*d.
        Args:
            t(int) : epaisseur de bois en mm
        """
        if t >= 12 * self.d_vis:

            a1 = 7 * self.d_vis
            a2 = 5 * self.d_vis
            a1CG = 10 * self.d_vis
            a2CG = 4 * self.d_vis
        else:
            print("L'épaisseur de bois n'est pas suffisante, il faut un bois de {0} mm minimum !".format(
                12*self.d_vis))
        return {"a1": a1, "a2": a2, "a1CG": a1CG, "a2CG": a2CG}


    @property
    def nefTraction(self):
        """Renvoie le nombre efficace de tirefond quand ils sont solicités par une composante parallèle à la partie lisse."""
        n = self.n * self.nfile
        if n > 1:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                n_ef_traction = n**0.9
                return n_ef_traction
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                n_ef_traction = 1
                return n_ef_traction
        return val()


    def faxk(self, l_ef:int, beam:str=["1", "2"]):
        """
        Calcul la valeur caractéristique de la résistance à l'arrachement perpendiculaire au fil en N/mm2 si 6mm<=d<=12mm
        et 0.6<=d1/d<=0.75.

        Args:
            l_ef (int): longueur de pénétration de la partie filetée en mm
            beam (str): élément à calculer 1 ou 2
        """
        d = self.d_vis.value*10**3
        d1 = self.d1.value*10**3
        if beam == "1":
            rho_k  = self.beam_1.rho_k
        else:
            rho_k  = self.beam_2.rho_k
            
        if 6 <= d <= 12 and 0.6 <= d1 / d <= 0.75:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                f_ax_k = 0.52 * (d ** -0.5) * (l_ef ** -0.1) * (rho_k ** 0.8)
                return f_ax_k
            return val()

        else:
            print(
                "le diamètre ne répond pas aux spécifications demandées en 8.7.2(4) de l'EN 1995 partie assemblage")


    def _FaxaRk(self, faxk:float, l_ef:int, alpha:int):
        """
        Calcul la valeur caractéristique de la résistance à l'arrachement du tirefond à un angle alpha par rapport au fil en N.

        Args:
            faxk (float): Valeur caractéristique de résistance à l'arrachement perpendiculaire au fil en N/mm2
            l_ef (int): longueur de pénétration de la partie filetée en mm
            alpha (int): angle formé entre l'axe du tirefond et le fil du bois, doit être supérieur à 30°
        """
        d = self.d_vis.value*10**3
        d_1 = self.d1
        rho_a = self.rho_a.value
        if self._type_beam[1] == "Métal":
            rho_k  = self.beam_1.rho_k
        else:
            rho_k  = self.beam_2.rho_k
            
        if 6 <= d <= 12 and 0.6 <= (d_1 / d) <= 0.75:
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                kd = min((d / 8), 1)
                F_ax_a_Rk = (faxk * d * l_ef * kd) / (1.2 * cos(radians(alpha)) ** 2 + sin(radians(alpha)) ** 2)
                return F_ax_a_Rk * si.N
        else:
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_ax_a_Rk = ((faxk * d * l_ef) / (1.2 * (cos(radians(alpha))) ** 2 + (sin(radians(alpha))) ** 2)) * ((rho_k / rho_a) ** 0.8) #N
                return F_ax_a_Rk * si.N
        return val()

    def _FaxaRkHead(self):
        """
        Calcul la valeur caractéristique de résistance à la traversée de la tête du tirefond dans l'assemblage en N.
        """
        f_head = self.fhead.value*10**-6
        d_h = self.dh.value*10**3
        rho_a = self.rho_a.value
        if self._type_beam[0] != "Métal":
            rho_k = self.beam_1.rho_k
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_head_Rk = f_head * d_h**2 * ((rho_k/rho_a)**0.8) #N
                return F_head_Rk * si.N
        else:
            @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
            def val():
                F_head_Rk = 10**6 * si.kN #l'élément 1 est métallique donc infini pour ce calcul
                return F_head_Rk
        return val()
    
    
    def Fax_Rk(self, faxk:float, l_ef:int, alpha:int):
        """Calcul la valeur caractéristique de résistance du tirefond axialement entre la résistance caractéristique de la tête et du pas de vis dans le bois en N.
        Cette fonction détermine ensuite la valeur de résistance caractéristique de l'assemblage en prenant compte du nombre efficace de tirefond en traction.

        Args:
            faxk (float): Valeur caractéristique de résistance à l'arrachement perpendiculaire au fil en N/mm2
            l_ef (int): longueur de pénétration de la partie filetée en mm
            alpha (int): angle formé entre l'axe du tirefond et le fil du bois, doit être supérieur à 30°
        """
        F_ax_a_Rk_value = self._FaxaRk(faxk, l_ef, alpha)
        F_ax_a_Rk = F_ax_a_Rk_value[1]
        F_head_Rk_value = self._FaxaRkHead()
        F_head_Rk = F_head_Rk_value[1]
        
        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_ax_Rk = min(F_ax_a_Rk, F_head_Rk)
            return F_ax_Rk
        FaxRk = val()
        self.FaxRk = FaxRk[1]

        f_ax_Rk = FaxRk[1]
        n_ef_traction = self.nefTraction[1]

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val_ass():
            F_ax_Rk_ass = f_ax_Rk * n_ef_traction
            return F_ax_Rk_ass
        F_ax_Rk_ass = val_ass()
        self.Fax_Rk_ass = F_ax_Rk_ass[1]
        
        return (F_ax_a_Rk_value[0] + F_head_Rk_value[0] + FaxRk[0] + F_ax_Rk_ass[0], F_ax_Rk_ass[1])


    def FtRk(self):
        """ Calcul la résistance caractéristique en traction pur du fil des tirefonds dans l'assemblage en N."""
        n_ef_traction = self.nefTraction[1]
        f_tens_k = self.ftensk
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            f_tRk = n_ef_traction * f_tens_k
            return f_tRk
        return val()
    

class Tirefond_inf_7(_Tirefond, Pointe):
    def __init__(self, d:si.mm, d1:float, ds:float, dh:float, l:float, n:int, rho_a:float, fhead:float, ftensk:float, MyRk:float=0, alpha1: float=0, alpha2: float=0, percage: bool=("False", "True"), **kwargs):
        """
        Crée une classe Tirefond_inf_7 pour les tirefonds avec un diamètre efficace inférieur à 7mm.
        Cette classe hérite de la classe Assemblage du module EC5_Assemblage.py.

        Args:
            d (int): diamètre extérieur du filet en mm
            d1 (float): diamètre du noyaux en mm
            ds (float): diamètre de la tige en mm
            dh (float): diamètre de la tête en mm
            l (float): longueur total de la vis en mm
            n (int): le nombre de vis dans une file
            rho_a (float): masse volumique associée au tirefond en fax,k en kg/m3
            fhead (float): valeur caractéristique de traversée de la tête du tirefond à l'EN 14592 en Mpa
            ftensk (float): valeur caractéristique en traction du tirefond en N
            MyRk (float, optional): le moment d'écoulement plastique de la vis en N.mm. 
                Si cette attribut est remplie, alors on récupère la valeur founis sinon on le calcul à l'EC5. 
                Defaults to 0.
            alpha1 (float, optional): angle entre l'effort de l'organe et le fil du bois 1 en °. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'organe et le fil du bois 2 en °. Defaults to 0.
            percage (bool, optional): l'élément est-il préperçé ? Si oui alors True. Defaults to ("False", "True").
        """
        qualite = "6.8"
        self.d_vis = d * si.mm
        if "qualite" in kwargs.keys():
            qualite = kwargs.pop("qualite")
                
        if d1*1.1 <= 6:
            Pointe.__init__(self, d=d1*1.1, dh=dh, l=l, qualite=qualite, n=n, alpha1=alpha1, alpha2=alpha2, type_organe="Tirefond", percage=percage, **kwargs)
            _Tirefond.__init__(self, d, d1, ds, dh, l, n, rho_a, fhead, ftensk, MyRk, alpha1, alpha2, percage)
        else:
            raise ValueError("Erreur, le tirefond est considéré comme un boulon et non une pointe")
        
        if MyRk:
            self._MyRk_fourni = MyRk * si.N*si.mm
        
    @property
    def MyRk(self):
        if hasattr(self, "_MyRk_fourni"):
            M_y_Rk_fourni = self._MyRk_fourni
            @handcalc(override="short", precision=2, left="\\[", right="\\]")
            def val():
                M_y_Rk = M_y_Rk_fourni
                return M_y_Rk
            return val()
        else:
            return super().MyRk


class Tirefond_sup_6(_Tirefond, Boulon):
    def __init__(self, d:si.mm, d1:float, ds:float, dh:float, l:si.mm, n: int, rho_a:float, fhead:float, ftensk:float, MyRk:float=0, alpha1: float=0, alpha2: float=0, **kwargs):
        """
        Crée une classe Tirefond_sup_6 pour les tirefonds avec un diamètre efficace supérieur à 6mm
        Cette classe hérite de la classe Assemblage du module EC5_Assemblage.py.

        Args:
            d (int): diamètre extérieur du filet en mm
            d1 (float): diamètre du noyaux en mm
            ds (float): diamètre de la tige en mm
            dh (float): diamètre de la tête en mm
            l (float): longueur total de la vis en mm
            n (int): le nombre de vis dans une file
            rho_a (float): masse volumique associée au tirefond en fax,k en kg/m3
            fhead (float): valeur caractéristique de traversée de la tête du tirefond à l'EN 14592 en Mpa
            ftensk (float): valeur caractéristique en traction du tirefond en N
            MyRk (float, optional): le moment d'écoulement plastique de la vis en N.mm. 
                Si cette attribut est remplie, alors on récupère la valeur founis sinon on le calcul à l'EC5. 
                Defaults to 0.
            alpha1 (float, optional): angle entre l'effort de l'organe et le fil du bois 1 en °. Defaults to 0.
            alpha2 (float, optional): angle entre l'effort de l'organe et le fil du bois 2 en °. Defaults to 0.
        """
        qualite = "6.8"
        self.d_vis = d * si.mm
        self.l = l * si.mm
        if "qualite" in kwargs.keys():
            qualite = kwargs.pop("qualite")
                
        if d1*1.1 > 6:
            Boulon.__init__(self, d=d1*1.1, qualite=qualite, n=n, alpha1=alpha1, alpha2=alpha2, type_organe="Tirefond", **kwargs)
            _Tirefond.__init__(self, d, d1, ds, dh, l, n, rho_a, fhead, ftensk, MyRk, alpha1, alpha2)
            
        else:
            raise ValueError("Erreur, le tirefond est considéré comme une pointe et non un boulon")
        
        if MyRk:
            self._MyRk_fourni = MyRk * si.N*si.mm
        
    @property
    def MyRk(self):
        if hasattr(self, "_MyRk_fourni"):
            M_y_Rk_fourni = self._MyRk_fourni
            @handcalc(override="short", precision=2, left="\\[", right="\\]")
            def val():
                M_y_Rk = M_y_Rk_fourni
                return M_y_Rk
            return val()
        else:
            return super().MyRk
# ======================================================= ANNEAU =========================================================
# 8.9 Assemblage par anneaux


# class Annneau(object):
#     """Défini un objet anneau avec :"""

#     def __init__(self, dc: float, t1:float, t2:float, hc:float, typeA="bois"):
#         self.type_organe = "Anneau"
#         self.dc = dc
#         self.t1 = t1
#         self.t2 = t2
#         self.he = hc/2
#         self.typeA = typeA

#     def ki(self, nAss:int, a3t:float, rhok:int):
#         """ Donne les facteur ki (de 1 à 4) dans un dico avec:
#             nAss : nombre d'assemblage par plan de cisaillement
#             a3t = distance d'extrémité chargé (en traction) """
#         listk = [0.0]*4

#         if nAss > 1:
#             ka = 1
#         else:
#             ka = 1.25

#         listk[0] = min(1,
#                        self.t1 / (3 * self.he),
#                        self.t2 / (5 * self.he))

#         listk[1] = min(ka,
#                        a3t / (2 * self.dc))

#         listk[2] = min(1.75,
#                        rhok / 350)

#         if self.typeA == "bois":
#             k4 = 1
#         else:
#             k4 = 1.25
#         listk[3] = k4
#         dico = {}
#         for i in range(1, 5):
#             cle = "k" + str(i)
#             dico[cle] = listk[i-1]

#         return dico

#     def Fv0Rk(self, dicoKi:dict):
#         """ Retourne la résistance en cidaillement de l'anneau en N avec:
#             dicoKi = dictionnaire des facteurs ki (def ki)"""
#         fv0rk = min(dicoKi["k1"] * dicoKi["k2"] * dicoKi["k3"] * dicoKi["k4"] * (35 * self.dc**1.5),
#                     dicoKi["k1"] * dicoKi["k3"] * self.he * (31.5 * self.dc))
#         return fv0rk
    

#     def FvaRk(self):
#         pass



