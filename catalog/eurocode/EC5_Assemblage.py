# coding in UTF-8
# by Anthony PARISOT
# under Python  3.9.4

############# Le but de ce fichier est de regrouper toute les fonctions d'assemblage par organe métalique dans l'EN-1995 #############
import os
import sys

import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Projet
from EC5_Element_droit import Barre, Cisaillement

def interpolationLineaire(x, xa, xb, ya, yb):
    """Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
    y = ya + (x - xa) * ((yb - ya)/(xb - xa))
    return y
    
# ====================================================== GENERAL =========================================================
# 8.1 Généralité
class Assemblage(Projet):
    GAMMA_M_ASS = 1.3
    DICO_COEF_LIMITE = {"Pointe circulaire": 0.15, "Pointe carrée": 0.25,
                     "Boulon": 0.25, "Autres pointes": 0.5, "Tirefond": 1}
    TYPE_BOIS_ASSEMBLAGE = ("Bois","PP/OSB", "CP", "Panneau dur")
    TYPE_ASSEMBLAGE = ("Bois/Bois", "Bois/Métal")

    def __init__(self,beam_1:object, beam_2:object, nfile: int=1, nCis: int=("1","2"), **kwargs):
        """Créer un objet Assemblage qui permet de calculer un assemblage bois/bois ou bois/métal à l'EN 1995.
        Cette classe est dérivé de la classe Projet du module A0_Project.py

        Args:
            beam_1 (object): objet correspondant à i=1, Beam  ou dérivé de cet objet provenant du module EC5_Element_droit.py
                             ou bien objet Element ou dérivé de cet objet provenant du module EC3_Element_droit.py
                             
            beam_2 (object): objet correspondant à i=2, Beam ou dérivé de cet objet provenant du module EC5_Element_droit.py
                             ou bien objet Element ou dérivé de cet objet provenant du module EC3_Element_droit.py
                             
            nfile (int, optional): le nombre de file dans l'assemblage. Defaults to 1.
            nCis (int, optional): Nombre de plan cisaillé entre 1 et 2. Defaults to ["1","2"].
        """
        
        super().__init__(**kwargs)
        #type_organe : le type d'organe de l'assemblage "Pointe circulaire", "Pointe carrée", "Autres pointes", "Agrafe", "Tirefond", "Boulon"ou"Broche", "plaque", "Anneau", "Crampon C1/C9","Crampon C10/C11".
        self.beam_1 = beam_1
        self.beam_2 = beam_2
        self.nfile = nfile
        # self.FaxRk = FaxRk
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
            
        
    # 7.1 Glissement des assemblages
    def Kser(self, perc: bool=("False", "True")):
        """Calcul le kser du type d'un organe et de l'asssemblage en N/mm 

        Args:
            perc (bool, optional): True si pré perçage de fait (valable que pour les pointes). Defaults to False.

        Returns:
            list: retourne une liste du kser d'un organe par plan de cisaillement et celui de l'assemblage
        """
        if self.type_organe == "Boulon" or self.type_organe == "Broche" or self.type_organe == "Tirefond":
            kser = self.rho_mean_ass**1.5 * self.d / 23
        elif self.type_organe == "Pointe circulaire" or self.type_organe == "Pointe carrée" or self.type_organe == "Autres pointes":
            if perc:
                kser = self.rho_mean_ass**1.5 * self.d / 23
            else:
                kser = self.rho_mean_ass**1.5 * self.d**0.8 / 30
        elif self.type_organe == "Agrafe":
            kser = self.rho_mean_ass**1.5 * self.d**0.8 / 80
        elif self.type_organe == "Anneau" or self.type_organe == "Crampon C10/C11":
            kser = self.rho_mean_ass * self.dc / 2
        elif self.type_organe == "Crampon C1/C9":
            kser = 1.5 * self.rho_mean_ass * self.d / 4
        
        ktype = 1
        if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[1]:
            ktype = 2
            
        kser_ass = kser * self.nfile * self.n * self.nCis * ktype
        return kser, kser_ass



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
        @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            F_90_Rk = 14 * b * w * sqrt(h_e / (1 - (h_e / h))) # N
            return F_90_Rk * si.N
        return val()

    # 8.2 Capacité résistante latérale pour les organes métaliques de type tige
    # 8.2.2 Assemblage bois/bois bois/panneau


    def _FvRk_BoisBois_Johansen(self):
        """Calcul la capacité résistante en cisaillement de la tige en N par plan de cisaillement avec
            t1 : valeur minimale entre epaisseur de l'élément bois latéral et la profondeur de pénétration en mm
            t2 : epaisseur de l'élément bois central en mm """
            
        f_h1k = self.fh1k[1].value * 10**-6
        f_h2k = self.fh2k[1].value * 10**-6
        t_1 = self.t1.value * 10**3
        t_2 = self.t2.value * 10**3
        diam = self.d.value * 10**3
        M_y_Rk = self.MyRk[1].value * 10**-6

        if self.nCis == 1:
            
            @handcalc(override="long", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                beta = f_h2k / f_h1k
                a = f_h1k * t_1 * diam # N
                b = f_h2k * t_2 * diam # N
                c = a/(1 + beta) * (sqrt(beta + 2 * beta**2 * (1 + t_2 / t_1 + (t_2 / t_1)**2) + beta**3 * (t_2 / t_1)**2) - beta * (1 + t_2 / t_1)) # N
                d = 1.05 * a / (2 + beta) * (sqrt(2 * beta * (1 + beta) + (4 * beta * (2 + beta) * M_y_Rk)/(f_h1k * t_1**2 * diam)) - beta) # N
                e = 1.05 * (f_h1k * t_2 * diam) / (1 + 2 * beta) * (sqrt(2 * beta**2 * (1 + beta) + (4 * beta * (1 + 2 * beta) * M_y_Rk)/(f_h1k * t_2**2 * diam)) - beta) # N
                f = 1.15 * sqrt((2 * beta)/(1 + beta)) * sqrt(2 * M_y_Rk * f_h1k * diam) # N
                return a, b, c, d, e, f
            
            calcul = val()
            a = calcul[1][0] * si.N
            b = calcul[1][1] * si.N
            c = calcul[1][2] * si.N
            d = calcul[1][3] * si.N
            e = calcul[1][4] * si.N
            f = calcul[1][5] * si.N
            dicoRupture = {a: "A", b: "B", c: "C", d: "D", e: "E", f: "F"}
            F_v_Rk_johansen = min(a, b, c, d, e, f)
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(a, b, c, d, e, f)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()

        else:
            @handcalc(override="long", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                beta = f_h2k / f_h1k
                g = f_h1k * t_1 * diam # N
                h = 0.5 * f_h2k * t_2 * diam # N
                j = 1.05 * g / (2 + beta) * (sqrt(2 * beta * (1 + beta) + (4 * beta * (2 + beta) * M_y_Rk)/ (f_h1k * t_1**2 * diam)) - beta) # N
                k = 1.15 * sqrt((2 * beta)/(1 + beta)) * sqrt(2 * M_y_Rk * f_h1k * diam) # N
                return g, h, j, k
            
            calcul = val()
            g = calcul[1][0] * si.N
            h = calcul[1][1] * si.N
            j = calcul[1][2] * si.N
            k = calcul[1][3] * si.N
            F_v_Rk_johansen = min(g, h, j, k)
            dicoRupture = {g: "G", h: "H", j: "J", k: "K"}
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(g, h, j, k)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()
            
        return (calcul[0] + calcul2[0], calcul2[1])


    # 8.2.3 Assemblage bois métal
    def _FvRk_BoisMetal_Johansen(self):
        """Calcul la capacité résistante en cisaillement de la tige en N par plan de cisaillement avec
            t1 : valeur minimale entre epaisseur de l'élément bois latéral et la profondeur de pénétration en mm
            t2 : epaisseur de l'élément bois central en mm
            """
        diam = self.d.value * 10**3
        M_y_Rk = self.MyRk[1].value * 10**-6

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
            
        if self.t <= 0.5 * self.d:
            self.type_plaque = "mince"

        elif self.d <= self.t:
            self.type_plaque = "epaisse"
            
        else:
            self.type_plaque = "intermédiaire"
            print("ATTENTION interpolation linéaire à faire ! EC5-8.2.3.1")

        if self.type_plaque == "mince" and self.nCis == 1:
            @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                a = 0.4 * f_h1k * t_1 * diam # N
                b = 1.15 * sqrt(2 * M_y_Rk * f_h1k * diam) # N
                return a, b

            calcul = val()
            a = calcul[1][0] * si.N
            b = calcul[1][1] * si.N
            dicoRupture = {a: "A", b: "B"}
            F_v_Rk_johansen = min(a, b)
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(a, b)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()

        elif self.type_plaque == "epaisse" and self.nCis == 1:
            @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                c = f_h1k * t_1 * diam # N
                d = c * (sqrt(2 + (4 * M_y_Rk) / (f_h1k * diam * t_1 ** 2)) - 1) # N
                e = 2.3 * sqrt(M_y_Rk * f_h1k * diam) # N
                return c, d, e
            
            calcul = val()
            c = calcul[1][0] * si.N
            d = calcul[1][1] * si.N
            e = calcul[1][2] * si.N
            dicoRupture = {c: "C", d: "D", e: "E"}
            F_v_Rk_johansen = min(c, d, e)
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(c, d, e)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()

        elif self.nCis == 2 and self.pos_plaq == "centrale":
            @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                f = f_h1k * t_1 * diam # N
                g = f * (sqrt(2 + (4 * M_y_Rk) / (f_h1k * diam * t_1 ** 2)) - 1) # N
                h = 2.3 * sqrt(M_y_Rk * f_h1k * diam) # N
                return f, g, h

            calcul = val()
            f = calcul[1][0] * si.N
            g = calcul[1][1] * si.N
            h = calcul[1][2] * si.N
            dicoRupture = {f: "F", g: "G", h: "H"}
            F_v_Rk_johansen = min(f, g, h)
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(f, g, h)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()

        elif self.type_plaque == "mince" and self.nCis == 2 and self.pos_plaq != "centrale":
            @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                j = 0.5 * f_h2k * t_2 * diam # N
                k = 1.15 * sqrt(2 * M_y_Rk * f_h2k * diam) # N
                return j, k

            calcul = val()
            j = calcul[1][0] * si.N
            k = calcul[1][1] * si.N
            dicoRupture = {j: "J", k: "K"}
            F_v_Rk_johansen = min(j, k)
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(j, k)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()

        else:
            @handcalc(override="short", precision=0, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                l = 0.5 * f_h2k * t_2 * diam # N
                m = 2.3 * sqrt(M_y_Rk * f_h2k * diam) # N
                return l, m

            calcul = val()
            l = calcul[1][0] * si.N
            m = calcul[1][1] * si.N
            dicoRupture = {l: "L", m: "M"}
            F_v_Rk_johansen = min(l, m)
            mode_rupture = dicoRupture[F_v_Rk_johansen]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val2():
                F_v_Rk_johansen = min(l, m)
                mode_rupture
                return F_v_Rk_johansen, mode_rupture
            calcul2 = val2()

        return (calcul[0] + calcul2[0], calcul2[1])

        
    def _FvRk_total(self, ass_bois: bool= True):
        """Calcul la capacité résistante en cisaillement caractéristique avec la partie de Johansen + l'effet de corde si il existe dans le
        mode de rupture avec :
            ass_bois(bool): si assemblage bois/bois alors True sinon False
        """
        coeflimit = self.DICO_COEF_LIMITE.get(self.type_organe, 0)
        if ass_bois and self.FvRk_Johansen[1] in ["C", "D", "E", "F", "J", "K"] or \
            not ass_bois and self.FvRk_Johansen[1] in ["B", "D", "E", "G", "H", "K", "M"]:

            F_ax_Rk = self.Fax_Rk
            F_v_Rk_johansen = self.FvRk_Johansen[0]

            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                F_ax_Rk_reel = min(F_ax_Rk / 4, coeflimit * F_v_Rk_johansen)
                F_v_Rk = F_v_Rk_johansen + F_ax_Rk_reel
                return F_v_Rk
            return val()
        else:
            return self.FvRk_Johansen[0]
    
        # 8.1.2 Assemblage par organe multiple

    def FvRk(self, effet_corde: bool=("True", "False")):
        """Calcul la valeur de calcul caractéristique de résistance au cisaillement de l'assemblage en N

        Args:
            effet_corde (bool): prise en compte de l'effet de corde, si oui alors True.

        Returns:
            float: effort de reprise caractéristique de l'assemblage en N
        """
        #     Fvrktot : capacité résistante en cisaillement caractéristique avec la partie de Johansen + l'effet de corde en N
        if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[0]:
            latex, self.FvRk_Johansen = self._FvRk_BoisBois_Johansen()
            assemblage_bois = True
        else:
            assemblage_bois = False
            latex, self.FvRk_Johansen = self._FvRk_BoisMetal_Johansen()

        if effet_corde:
            FvRk_latex, self.Fv_Rk = self._FvRk_total(assemblage_bois)
            latex = latex + FvRk_latex

        else:
            self.Fv_Rk = self.FvRk_Johansen[0]

        F_v_Rk = self.Fv_Rk
        n_file = self.nfile
        n_ef = self._nef
        n_cisaillement = self.nCis
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            F_v_Rk_ass = F_v_Rk * n_file * n_ef * n_cisaillement
            return F_v_Rk_ass
        FvRkass_latex, self.Fv_Rk_ass = val()
        latex = latex + FvRkass_latex
        return (latex, self.Fv_Rk_ass)


    def F_Rd(self, F_rk: float):
        """Calcul la valeur de calcul (design) de résistance de l'assemblage en N avec :

        Args:
            F_rk (float): capacité résistante caractéristique de l'organe en N
            num_beam (int, optional): numéro de l'élément à vérifier. Defaults to 1."""
        gamma_M = self.GAMMA_M_ASS
        if "Métal" in self._type_beam:
            if self._type_beam[0] == "Métal":
                k_mod = self.beam_2.K_mod
            else:
                k_mod = self.beam_1.K_mod
        else:
            k_mod = sqrt(self.beam_1.K_mod * self.beam_2.K_mod)
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            F_rd = (F_rk * k_mod) /gamma_M # Valeur de calcul (design)
            return F_rd
        return val()
    
    
    # Annexe A : Cisaillement de bloc
    def FbsRk(self, dp:float, a1:float, a2:float, a3t:float, Kcr: float=0.67, num_beam: int=("1", "2")):
        """Calcul la valeur caractéristique en cisaillement de bloc en N pour l'élément 1 ou 2 de l'assemblage avec:

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

        mode_rupture = self.FvRk_Johansen[1]

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

        @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def l_net():
            l_v_i = a_1 - diam # Distance entre perçage dans le sens fil
            L_net_v = 2 * (l_v_i * (n - 1) + a_3_t - (diam_percage / 2)) # Longueur résiduelle de la surface de rupture en cisaillement
            
            l_t_i = a_2 - diam_percage # Distance entre perçage dans le sens perpendiculaire au fil
            L_net_t = l_t_i * (n_file - 1) # Largeur résiduelle de la section perpendiculaire au fil

            mode_rupture
            return L_net_v, L_net_t
        
        L_net_value = l_net()
        latex = L_net_value[0]
        L_net_v, L_net_t = L_net_value[1]

        
        if mode_rupture in ("C", "F", "J", "L", "K", "M"):
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def a_net_v():
                A_net_v = L_net_v * (t_1 * K_cr)
                return A_net_v
            a_net_v_result = a_net_v()
            latex = latex + tef_result[0] + a_net_v_result[0]
            A_net_v = a_net_v_result[1]
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def a_net_v(t_ef):
                A_net_v = L_net_v / 2 * (L_net_t + 2 * (t_ef * K_cr))
                return A_net_v
            
            match mode_rupture:
                case "A":
                    @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
                    def tef():
                        t_ef = 0.4*t_1 # Épaisseur efficace
                        return t_ef
                case "B":
                    @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
                    def tef():
                        t_ef = 1.4 * sqrt(M_y_Rk / (f_hk * diam)) # Épaisseur efficace
                        return t_ef
                case "D"|"G":
                    @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
                    def tef():
                        t_ef = t_1 * (sqrt(2 + (4 * M_y_Rk ) / (f_hk * diam * t_1**2)-1)) # Épaisseur efficace
                        return t_ef
                case "E"|"H":
                    @handcalc(override="long", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
                    def tef():
                        t_ef = 2 * sqrt(M_y_Rk / (f_hk * diam)) # Épaisseur efficace
                        return t_ef
            tef_result = tef()
            a_net_v_result = a_net_v(tef_result[1])
            latex = latex + tef_result[0] + a_net_v_result[0]
            A_net_v = a_net_v_result[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def f_bs_Rk():
            A_net_t = L_net_t * t_1
            F_bs_Rk = max(0.7 * A_net_v * f_v_k, 1.5 * A_net_t * f_t0_k)
            return F_bs_Rk
        
        result = f_bs_Rk()
        self.F_bs_Rk = result[1]
        return (latex + result[0], self.F_bs_Rk)
    

# ======================================================= POINTE =========================================================
# 8.3 Assemblage par pointes

class Pointe(Assemblage):
    """ Défini un objet pointe avec :
        d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
        n : nombre d'organe dans une file 
        fu : la résistance caractéristique en traction du fil d'acier en N/mm2
        alpha : angle entre l'effort de l'organe et le fil du bois en °
        pk : masse volumique caractéristique du bois en kg/m3
        type_pointe : "Carrée" = False "Circulaire" = True
        """

    TYPE_ASSEMBLAGE = ("Bois/Bois", ("CP", "Panneau dur", "PP/OSB"), "Bois/Métal")
    QUALITE_ACIER = ('6.8', '8.8', '9.8', '10.9', '12.9')
    TYPE_ORGANE = ("Pointe circulaire", "Pointe carrée", "Autres pointes")

    def __init__(self, d:float, d_tete:int, l:float, qualite: str=QUALITE_ACIER, n: int=1, alpha: float=0, type_organe: str=TYPE_ORGANE, percage: bool=("False", "True"), *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.d = d * si.mm
        self.d_tete = d_tete * si.mm
        self.l = l * si.mm #longueur sous tête
        self.qualite = qualite
        self.n = n
        self.fu = self.__qualite_acier.loc["fub"]
        self.alpha = mt.radians(alpha)
        self.percage = percage
        self.type_organe = type_organe

        self.__t1_t2()
        self._fhik()


    @property
    def __qualite_acier(self):
        df = self._data_from_csv("qualite_acier.csv")
        df = df.loc[self.qualite]
        return df

    @property
    def _type_circulaire(self):
        if self.type_organe == "Pointe carrée":
            return False
        else:
            return True


    def __t1_t2(self):
        """Retourne t1 et t2 en mm suivant l'EN 1995 §8.3.1.1
        """
        if self.nCis == 1:
            if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t1 = self.beam_1.b_calcul
            else:
                self.t1 = self.beam_1.t

            if self._type_beam[1] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t2 = self.l-self.t1

        else:
            if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t1 = min(self.beam_1.b_calcul, self.l-self.beam_1.b_calcul)
            else:
                print("Il n'est pas possible sauf erreur d'avoir un assemblage double cisaillement avec une pointe !")
            if self._type_beam[1] in self.TYPE_BOIS_ASSEMBLAGE:
                self.t2 = self.beam_2.b_calcul
            else:
                self.t2 = self.beam_2.t


    @property
    def MyRk(self):
        """ Défini le moment d'écoulement plastique d'une pointe en N.mm avec:
            d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            fu : la résistance caractéristique en traction du fil d'acier en N/mm2 """
        if self.fu >= 600:
            if self._type_circulaire == True:
                MyRk = 0.3 * self.fu * self.d**2.6
            else:
                MyRk = 0.45 * self.fu * self.d**2.6
            return MyRk
        else:
            print("La résistance du fil en traction est inférieur à 600 MPa, vérifier vos données !")


    def _fhk_bois(self, beam:object):
        """ Calcul la portance locale des pointes inférieur à 8mm dans le bois et le LVL en MPa
            """
        if self.percage:
            fhk = 0.082 * (1 - 0.01 * self.d) * beam.rho_k
        else:
            fhk = 0.082 * beam.rho_k * self.d**(-0.3)
        return fhk


    def _fhk_panneau(self, beam:object):
        """Calcul la portance locale des pointes dans les panneaux en MPa

        Args:
            t (int): épaisseur du panneau en mm
            self.d_tete (int): diamètre de la tête de la pointe

        Returns:
            float: portance locale en MPa
        """
        if self.d_tete >= 2*self.d:
            if beam.type_bois == "CP":
                fhk = 0.11 * beam.rho_k * self.d**-0.3
            elif beam.type_bois == "Panneau dur":
                fhk = 30 * self.d**-0.3 * beam.b_calcul**0.6
            else:
                fhk = 65 * self.d**-0.7 * beam.b_calcul**0.1
            return fhk
        else:
            print(f"La tête de la pointe doit être au moins égale à {2*self.d} mm")


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


    def _kef(self, a1:int):
        """ coefficient donnée dans le tableau 8.1 fonction de a1 et du percage qui réduit le nombre efficace de pointe dans le sens 
            du fil avec :
                a1 : l'espacement entre tige dans le sens du fil du bois """
        listeTab = (self.d * 4, self.d * 7, self.d *
                    10, self.d * 14, ("x", 0.7, 0.85, 1))

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


    def nef(self, a1:int):
        """Défini le nombre efficace d'organe dans une file avec :
            n : nombre de boulons dans une file 
            d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            percage : True or False
            kef : coefficient donnée dans le tableau 8.1 fonction de a1 et du percage """
        self._nef = self.n**self._kef(a1)
        return self._nef


    def _pince(self, beam:object):
        """Défini les différentes pinces minimales pour une pointe en mm"""

        rho_k = beam.rho_k

        if self.percage:
            a1 = round((4 + mt.cos(self.alpha)) * self.d, 1)
            a2 = round((3 + mt.sin(self.alpha)) * self.d, 1)
            a3t = round((7 + 5 * mt.cos(self.alpha)) * self.d, 1)
            a3c = round(7 * self.d, 1)
            a4c = round(3 * self.d, 1)

            if self.d < 5:
                a4t = round((3 + 2 * mt.sin(self.alpha)) * self.d, 1)
            else:
                a4t = round((3 + 4 * mt.sin(self.alpha)) * self.d, 1)
        else:
            if rho_k <= 420:
                if self.d < 5:
                    a1 = round((5 + 5 * mt.cos(self.alpha)) * self.d, 1)
                    a4t = round((5 + 2 * mt.sin(self.alpha)) * self.d, 1)
                else:
                    a1 = round((3 + 7 * mt.cos(self.alpha)) * self.d, 1)
                    a4t = round((5 + 5 * mt.sin(self.alpha)) * self.d, 1)

                a2 = round(5 * self.d, 1)
                a3t = round((10 + 5 * mt.cos(self.alpha)) * self.d, 1)
                a3c = round(10 * self.d, 1)
                a4c = a2

            elif rho_k > 420 and rho_k <= 500:
                a1 = round((7 + 8 * mt.cos(self.alpha)) * self.d, 1)
                a2 = round(7 * self.d, 1)
                a3t = round((15 + 5 * mt.cos(self.alpha)) * self.d, 1)
                a3c = round(15 * self.d, 1)
                a4c = a2

                if self.d < 5:
                    a4t = round((7 + 2 * mt.sin(self.alpha)) * self.d, 1)
                else:
                    a4t = round((7 + 5 * mt.sin(self.alpha)) * self.d, 1)
            else:
                print("Il faut absolument prépercer au dessus de rho,k: 500 kg/m3")
        
        if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[0]: #Si assemblage bois bois
            pass

        elif self.type_assemblage == __class__.TYPE_ASSEMBLAGE[1]: #Si assemblage bois métal
            a1 = round(a1 * 0.7, 1)
            a2 = round(a2 * 0.7, 1)

        else: #Si assemblage bois panneau
            a1 = round(a1 * 0.85, 1)
            a2 = round(a2 * 0.85, 1)
            if beam.type_bois == "CP":
                a3t = round((3 + 4 * mt.sin(self.alpha)) * self.d, 1)
                a3c = round(3 * self.d, 1)
                a4t = a3t
                a4c = a3c

        return {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}


    @property
    def pince_ass(self):
        pince = []
        if self._type_beam[0] in self.TYPE_BOIS_ASSEMBLAGE:
            pince.append(self._pince(self.beam_1))
        else:
            pince.append(None)

        if self._type_beam[1] in self.TYPE_BOIS_ASSEMBLAGE:
            pince.append(self._pince(self.beam_2))
        else:
            pince.append(None)

        return pince


    def prepercage(self, sensible: bool=("False", "True")):
        """Défini l'épaisseur en mm minimale pour éviter le pré-perçage des pointes

        Args:
            sensible (bool, optional): Défini si le bois utilisé est sensible à la fissuration (selon AN §8.3.1.2(7) douglas et pin maritime). Defaults to False.

        Returns:
            int: l'épaisseur du bois en mm
        """
        if sensible:
            t = max(14 * self.d,
                    (13 * self.d - 30) * (self.pk / 200))
        else:
            t = max(7 * self.d,
                    (13 * self.d - 30) * (self.pk / 400))
        print(f"l'épaisseur minimale pour éviter le pré-perçage est de {t} mm")
        return int(t) 
        
        
        



# ======================================================= BOULON =========================================================
# 8.5 Assemblage par boulon

class Boulon(Assemblage):
    """ Défini un objet boulon avec :
        d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm
        fuk : la valeur caractéristique de résistance à la traction du boulon en N/mm2
        n : nombre de boulons dans une file
        alpha : angle entre l'effort de l'organe et le fil du bois en ° par barre 

        t1 (int, optional): longueur de contact avec la tige  pour la pièce 1 en mm. 
                                ATTENTION : Cette argument n'est pas obligatoire par défaut, il est calculer par le type de tige utilisé.
                                Il n'est nécessaire de le remplir que si vous avez un t1 spécifique, par exemple avec une chapelle réduisant ainsi la portance local à une longueur inférieur à celle de l'épaisseur de la pièce 1.

        t2 (int, optional): longueur de contact avec la tige  pour la pièce 2 en mm. 
                            ATTENTION : Même chose que pour t1 mais pour la pièce 2.
    """
        
    QUALITE_ACIER = tuple(Assemblage._data_from_csv(Assemblage, "qualite_acier.csv").index)
    
    def __init__(self, d:float, qualite: str=QUALITE_ACIER, n: int=1, alpha1: float=0, alpha2: float=0, t1: int=0, t2: int=0, **kwargs):
        super().__init__(**kwargs)
        self.type_organe = "Boulon"
        self.d = d * si.mm
        self.qualite = qualite
        self.fuk = self.__qualite_acier.loc["fub"] *si.MPa
        self.n = n
        self._nef = n
        self.alpha = [alpha1, alpha2]
        self.t1 = t1
        self.t2 = t2

        self._fhik()
    
    @property
    def __qualite_acier(self):
        df = self._data_from_csv("qualite_acier.csv")
        df = df.loc[self.qualite]
        return df

        
    # 8.5.1 Boulons chargés latéralement
    # 8.5.1.1 Généralité et assemblage bois/bois
    
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
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
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
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
            def val():
                f_hk = 0.11 * (1 - 0.01 * d) * rho_k # MPa
                return f_hk
        else:
            @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
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

        
            if self._type_beam[i] != "Métal":
                if i:
                    if not self.t2:
                        self.t2 = self.beam_2.b_calcul
                    else:
                        self.t2 = self.t2 * si.mm
                else:
                    if not self.t1:
                        self.t1 = self.beam_1.b_calcul
                    else:
                        self.t1 = self.t1 * si.mm
        
        
        self.fh1k = dict_beam["1"]["fhik"]
        self.fh2k = dict_beam["2"]["fhik"]
        return self.fh1k, self.fh2k


    @property
    def MyRk(self):
        """Défini le moment d'écoulement plastique d'un boulon en N.mm avec:
            fuk : la valeur caractéristique de résistance à la traction du boulon en N/mm2
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm"""
        f_uk = self.fuk
        d = self.d.value * 10**3
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            M_y_Rk = 0.3 * f_uk * d ** 2.6
            return M_y_Rk
        return val()


    @property
    def pince(self):
        """Défini les différentes pinces minimales pour un boulon en mm avec :
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm """
        dict_pince = {}
        for i, alpha in enumerate(self.alpha):
            a1 = round((4 + mt.cos(mt.radians(alpha))) * self.d, 1)
            a2 = round(4 * self.d, 1)
            a3t = max(7 * self.d, 80*si.mm)

            if alpha <= 150 and alpha < 210:
                a3c = round(4 * self.d, 1)
            else:
                a3c = round((1 + 6 * mt.sin(mt.radians(alpha))) * self.d, 1)

            a4t = round(max((2 + 2 * mt.sin(mt.radians(alpha))) * self.d, 3 * self.d), 1)
            a4c = round(3 * self.d, 1)
            dict_pince["barre "+str(i+1)] = {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}
        return dict_pince


    def nef(self, a1:int):
        """Défini le nombre efficace d'organe (boulon) dans une file avec :
            a1 : l'espacement entre boulon dans le sens du fil du bois
            n : nombre de boulons dans une file 
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm"""
        n = self.n
        a1 = a1 * si.mm
        d = self.d
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            n_ef = min(n**0.9 * (a1/(13 * d))**(1/4), n)
            return n_ef
        value = val()
        self._nef = value[1]
        return value


    # 8.5.2 Boulons chargés axialement
    def FaxRk(self, d_int: float=0, d_ext: float=0, filetage_EN1090: bool=("True", "False"), loadtype=Barre.LOAD_TIME, typecombi=Barre.TYPE_ACTION):
        """Calcul la résistance axial caractéristique d'un boulon chargé axialement à partir soit de la rondelle soit de la plaque métalique.

        Args:
            d_int (float, optional): diamètre intérieur de la rondelle en mm ou du trou de perçage dans la plaque métallique. Defaults to 0.
            d_ext (float, optional): diamètre extérieur de la rondelle en mm. Defaults to 0.
            filetage_EN1090 (bool, optional): défini si le filetage est conforme à l'EN 1090, soit matricé. Si filetage usiné alors False. Defaults to True.
            loadtype (str): chargement de plus courte durée sur l'élément.
            typecombi (str): type de combinaison, fondamentale ou accidentelle.

        Returns:
            FaxRk: la résistance axial d'un boulon en N
        """
        from EC3_Assemblage import Tige
        d_int = d_int * si.mm
        d_ext = d_ext * si.mm

        if self._type_beam[0] == "Métal":
            FtRd = Tige(self.d.value*10**3, d_int, self.qualite, True, filetage_EN1090, t=self.beam_1.t.value*10**3, h=self.beam_1.h.value*10**3, classe_acier=self.beam_1.classe_acier, classe_transv=self.beam_1.classe_transv).FtRd
            d_ext = min(self.beam_1.t*12, 4*self.d)
            fc_90_d = self.beam_2._f_type_d("fc90k", loadtype, typecombi)[1]
        else:
            FtRd = Tige(self.d.value*10**3, d_int, self.qualite, True, filetage_EN1090, t=0, h=0, classe_acier="S235", classe_transv=3).FtRd
            fc_90_d = self.beam_1._f_type_d("fc90k", loadtype, typecombi)[1]

        Ft_Rd = FtRd[1]
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
        def val():
            A_int = pi * (d_int / 2)**2 
            A_rondelle = pi * (d_ext / 2)**2 - A_int

            f_c90_d_rond = fc_90_d * 3
            F_c90_d = f_c90_d_rond * A_rondelle
            F_ax_Rk = min(F_c90_d, Ft_Rd)
            return F_ax_Rk
        FaxRk = val()
        self.Fax_Rk = FaxRk[1]
        return (FtRd[0] + FaxRk[0], FaxRk[1])
        

# ======================================================= BROCHE =========================================================
# 8.6 Assemblage par broche

class Broche(Boulon):
    """ Créer une classe Broche.
        d : diamètre efficace de la broche ( entre 6 et 30 mm) en  mm
        qualite : qualite de la broche
        n : nombre de broche dans une file
        alpha : angle entre l'effort de l'organe et le fil du bois en °

        t1 (int, optional): longueur de contacte avec la tige  pour la pièce 1 en mm. 
                                ATTENTION : Cette argument n'est pas obligatoire par défaut, il est calculer par le type de tige utilisé.
                                Il n'est nécessaire de le remplir que si vous avez un t1 spécifique, par exemple avec une chapelle réduisant ainsi la portance local à une longueur inférieur à celle de l'épaisseur de la pièce 1.

        t2 (int, optional): longueur de contacte avec la tige  pour la pièce 2 en mm.
                            ATTENTION : Même chose que pour t1 mais pour la pièce 2.
    """

    def __init__(self, d:float, qualite: str=Boulon.QUALITE_ACIER, n: int=1, alpha1: float=0, alpha2: float=0, t1: int=0, t2: int=0, **kwargs):
        super().__init__(d, qualite, n, alpha1, alpha2, **kwargs)
        self.type_organe = "Broche"

    @property
    def FaxRk(self):
        self.Fax_Rk = 0
        
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


class Tirefond(object):
    """ Défini un object tirefond avec :
        d : diamètre extérieur du filet en mm
        d1 : diamètre du noyaux en mm
        dh : diamètre de la tête en mm
        pa : masse volumique associée au tirefond en fax,k en kg/m3
        fhead : valeur caractéristique de traversée de la tête du tirefond à l'EN14592 en Mpa
        ftensk : valeur caractéristique en traction du tirefond en MPa
        n : nombre de boulons dans une file
        alphaTirefond : angle formé entre l'axe du tirefond et le fil du bois, doit être supérieur à 30°"""

    def __init__(self, d:float, d1:float, dh:float,  pa:float, fhead:float, ftensk:float, n, alphaTirefond: float=90, **kwargs):
        self.d = d * si.mm
        self.d1 = d1 * si.mm
        self.d_ef = d1*1.1
        if self.d_ef <= 6:
            super(Pointe, self).__init__(d=self.d_ef, d_tete=dh, **kwargs)
        self.type_organe = "Tirefond"
        
        self.dh = dh * si.mm
        self.pa = pa * si.kg/si.m**3
        self.fhead = fhead * si.MPa
        self.ftensk = ftensk * si.MPa
        self.n = n
        self.alphaTirefond = mt.radians(alphaTirefond)


    # 8.7.2 Tirefond chargés axialement
    def pinceTirefondAxial(self, t: int):
        """ Défini les pinces d'un tirefond en mm lorsqu'il est chargée axialement et l'epaisseur de bois supérieur à 12*d avec:
        d : diamètre extérieur du filet en mm
        t : epaisseur de bois en mm """
        if t >= 12 * self.d:

            a1 = 7 * self.d
            a2 = 5 * self.d
            a1CG = 10 * self.d
            a2CG = 4 * self.d

        else:

            print("L'épaisseur de bois n'est pas suffisante, il faut un bois de {0} mm minimum !".format(
                12*self.d))
        return {"a1": a1, "a2": a2, "a1CG": a1CG, "a2CG": a2CG}


    def nefTraction(self):
        """ Renvoie le nombre efficace de tirefond quand ils sont solicité par une composante parallèle à la partie lisse avec :
            n = nombre de tirefond agissant conjointement dans l'assemblage """
        if self.n > 1:

            self._nefTraction = self.n**0.9

        else:

            self._nefTraction = 1

        return self._nefTraction


    def faxk(self, lef: int, pk: int=350):
        """Calcul la valeur caractéristique de la résistance à l'arrachement perpendiculaire au fil en N/mm2 si 6mm<=d<=12mm
        et 0.6<=d1/d<=0.75 avec :
            d : diamètre extérieur du filet en mm
            d1 : diamètre du noyaux en mm
            lef : longueur de pénétration de la partie filetée en mm
            pk : masse volumique caractéristique en kg/m3"""
        if 6 <= self.d <= 12 and 0.6 <= (self.d1 / self.d) <= 0.75:

            faxk = 0.52 * (self.d ** -0.5) * (lef ** -0.1) * (pk ** 0.8)
            return faxk

        else:
            print(
                "le diamètre ne répond pas au spécification demandée en 8.7.2(4) de l'EN 1995 partie assemblage")


    def FaxaRk(self, faxk:float, lef:int, pk: int=350):
        """Calcul la valeur caractéristique de la résistance à l'arrachement du tirefond à un angle alpha par rapport au fil en N avec :
                _nef : nombre efficace de tirefond en traction compression
                faxk : Valeur caractéristique de résistance à l'arrachement perpendiculaire au fil en N/mm2
                d : diamètre extérieur du filet en mm
                d1 : diamètre du noyaux en mm
                lef : longueur de pénétration de la partie filetée en mm
                pk : masse volumique caractéristique en kg/m3
                pa : masse volumique associée au tirefond en fax,k en kg/m3
                alpha : angle formé entre l'axe du tirefond et le fil du bois, doit être supérieur à 30°"""
        if 6 <= self.d <= 12 and 0.6 <= (self.d1 / self.d) <= 0.75:

            kd = min((self.d / 8), 1)
            faxark = (self._nefTraction() * faxk * self.d * lef * kd) / (1.2 * mt.cos(self.alphaTirefond) ** 2 +
                                                                           mt.sin(self.alphaTirefond) ** 2)

        else:

            faxark = ((self._nefTraction() * faxk * self.d * lef) / (1.2 * (mt.cos(self.alphaTirefond)) ** 2 +
                                                                       (mt.sin(self.alphaTirefond)) ** 2)) * ((pk / self.pa) ** 0.8)

        return faxark


    def FaxaRkHead(self, pk:int):
        """Calcul la valeur caractéristique de résistance à la traversée de l'assemblage en N avec :
            _nef : nombre efficace de tirefond en traction compression
            fhead : valeur caractéristique de traversée de la tête du tirefond à l'EN14592 en Mpa
            dh : diamètre de la tête en mm
            pk : masse volumique caractéristique en kg/m3
            pa : masse volumique associée au tirefond en fax,k en kg/m3 """
        FaxRkhead = self._nefTraction() * self.fhead * self.dh**2 * ((pk/self.pa)**0.8)
        return FaxRkhead


    def FtRk(self):
        """ Calcul la résistance caractéristique en traction du tirefond en N avec:
            _nef : nombre efficace de tirefond en traction compression
            ftensk : valeur caractéristique en traction du tirefond en MPa """
        ftRk = self._nefTraction() * self.ftensk
        return ftRk

# ======================================================= ANNEAU =========================================================
# 8.9 Assemblage par anneaux


class Annneau(object):
    """Défini un objet anneau avec :"""

    def __init__(self, dc: float, t1:float, t2:float, hc:float, typeA="bois"):
        self.type_organe = "Anneau"
        self.dc = dc
        self.t1 = t1
        self.t2 = t2
        self.he = hc/2
        self.typeA = typeA

    def ki(self, nAss:int, a3t:float, rhok:int):
        """ Donne les facteur ki (de 1 à 4) dans un dico avec:
            nAss : nombre d'assemblage par plan de cisaillement
            a3t = distance d'extrémité chargé (en traction) """
        listk = [0.0]*4

        if nAss > 1:
            ka = 1
        else:
            ka = 1.25

        listk[0] = min(1,
                       self.t1 / (3 * self.he),
                       self.t2 / (5 * self.he))

        listk[1] = min(ka,
                       a3t / (2 * self.dc))

        listk[2] = min(1.75,
                       rhok / 350)

        if self.typeA == "bois":
            k4 = 1
        else:
            k4 = 1.25
        listk[3] = k4
        dico = {}
        for i in range(1, 5):
            cle = "k" + str(i)
            dico[cle] = listk[i-1]

        return dico

    def Fv0Rk(self, dicoKi:dict):
        """ Retourne la résistance en cidaillement de l'anneau en N avec:
            dicoKi = dictionnaire des facteurs ki (def ki)"""
        fv0rk = min(dicoKi["k1"] * dicoKi["k2"] * dicoKi["k3"] * dicoKi["k4"] * (35 * self.dc**1.5),
                    dicoKi["k1"] * dicoKi["k3"] * self.he * (31.5 * self.dc))
        return fv0rk
    

    def FvaRk(self):
        pass


if __name__ == "__main__":
    from EC5_Element_droit import Barre
    from EC3_Element_droit import Element
    
    beam1 = Barre(60, 140, "Rectangulaire", classe="C24", cs=2)
    beam1._f_type_d("fm0k", "Instantanee", "Fondamentales")

    beam2 = Barre(160, 160, "Rectangulaire", classe="GL24h", cs=2)
    beam1._f_type_d("fm0k", "Instantanee", "Fondamentales")

    ass = Assemblage(beam1, beam2, nfile=1, nCis=2)
    bl = Boulon._from_parent_class(ass, d=12, qualite="4.6", n=2, alpha=0)
    bl.FaxRk(14,30,True,"Court terme", "Fondamentales")
    print(bl.FvRk(True))

    # pointe = Pointe._from_parent_class(ass, d=6, d_tete=12, l=120, qualite="10.9", n=1, alpha=0, type_organe="Autres pointes", percage="False")
    # print(pointe.pince_ass)

