# coding in UTF-8
# by Anthony PARISOT
# under Python  3.9.4

############# Le but de ce fichier est de regrouper toute les fonctions d'assemblage par organe métalique dans l'EN-1995 #############
import os
import sys

import math as mt

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Project

def interpolationLineaire(x, xa, xb, ya, yb):
    """Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
    y = ya + (x - xa) * ((yb - ya)/(xb - xa))
    return y
    
# ====================================================== GENERAL =========================================================
# 8.1 Généralité
class Assemblage(Project):
    GAMMA_M_ASS = 1.3
    DICO_COEF_LIMITE = {"Pointe circulaire": 0.15, "Pointe carrée": 0.25,
                     "Boulon": 0.25, "Autres pointes": 0.5, "Tirefond": 1}
    TYPE_ASSEMBLAGE = ["Bois/Bois","Bois/Métal"]

    def __init__(self,beam_1:object, beam_2:object, nfile: int=1, nCis: int=["1","2"], **kwargs):
        """Créer un objet Assemblage qui permet de calculer un assemblage bois/bois ou bois/métal à l'EN 1995.
        Cette classe est dérivé de la classe Project du module A0_Project.py

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
        # self.t1 = t1
        # self.t2 = t2
        # self.fh1k = fh1k
        # self.fh2k = fh2k
        self.nCis = nCis
        
        self.__type_assemblage()
    
    def __type_assemblage(self):
        self._type_beam = []
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            try:
                if beam.type_bois:
                    if beam.type_bois in ["Massif","BLC", "LVL"]:
                        self._type_beam.append("Bois")
                    elif beam.type_bois in ["OSB 2", "OSB 3/4"]:
                        self._type_beam.append("PP/OSB")
                    else:
                        self._type_beam.append("CP")
                    self.K_mod = beam.K_mod
            except AttributeError:
                self._type_beam.append("Métal")
                
        if self._type_beam[0] == "Bois":
            self.type_assemblage = self._type_beam[0] + "/" + self._type_beam[1]
        else:
            self.type_assemblage = self._type_beam[1] + "/" + self._type_beam[0]
        
    @property
    def rho_mean_ass(self):
        if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[0]:
            rho_m1 = int(self.beam_1.caract_meca.loc["rhomean"])
            rho_m2 = int(self.beam_2.caract_meca.loc["rhomean"])
            return mt.sqrt(rho_m1 * rho_m2)
        else:
            if self._type_beam[0] == "Bois":
                return int(self.beam_1.caract_meca.loc["rhomean"])
            else:
                return int(self.beam_2.caract_meca.loc["rhomean"])
            
        
    # 7.1 Glissement des assemblages
    def Kser(self, perc=False)->list:
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


    # 8.1.2 Assemblage par organe multiple

    def FvRd(self, effet_corde: bool=["True", "False"]) -> float:
        """Calcul la valeur de calcul (design) de résistance au cisaillement de l'assemblage en N

        Args:
            effet_corde (bool): prise en compte de l'effet de corde, si oui alors True.

        Returns:
            float: retourne Fv_Rd de l'assemblage en N
        """
        #     Fvrktot : capacité résistante en cisaillement caractéristique avec la partie de Johansen + l'effet de corde en N
        if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[0]:
            self.FvRk_Johansen = self._FvRk_BoisBois_Johansen()
            if effet_corde:
                self.FvRk = self._FvRk_BoisBois_Tot()
        else:
            self.FvRk_Johansen = self._FvRk_BoisMetal_Johansen()
            if effet_corde:
                self.FvRk = self._FvRk_BoisMetal_Tot()
        
        fvrd = (self.FvRk * self.nfile * self._nef * self.nCis * self.K_mod) / __class__.GAMMA_M_ASS
        return fvrd

    def F_Rd(self, F_rk):
        """Calcul la valeur de calcul (design) de résistance de l'assemblage en N avec :
            F_rk : capacité résistante caractéristique de l'organe en N"""
        f_rd = (F_rk * self.K_mod) / __class__.GAMMA_M_ASS
        return f_rd



    # 8.1.4 Effort d'assemblage inclinés par rapport au fil

    def _w(self, wpl=0):
        """Calcul le facteur de modification pour le calcul de la valeur caractéristique au fendage F90,Rk avec:
            wpl : largeur de la plaque métallique emboutie parallèlement au fil en mm
            type_organe : type d'organe utilisé, pour les plaques métalliques embouties : "plaque", pour les autres : "autres" """
        if self.type_organe == "plaques métaliques embouties":
            w = max((wpl / 100) ** 0.35, 1)
        else:
            w = 1
        return w


    def F90Rk(self, b, h, he, w=1):
        """Calcul la valeur caractérisque de la capacité au fendage en N avec :
            b : l'épaisseur de l'élément en mm
            h : la hauteur de l'élément en mm
            he : la distance de rive chargée vis à vis du centre de l'organe le plus plus éloigné ou du bord de la plaque
            w : facteur de modification """
        f90rk = 14 * b * w * mt.sqrt(he / (1 - (he / h)))
        return f90rk

    # 8.2 Capacité résistante latérale pour les organes métaliques de type tige
    # 8.2.2 Assemblage bois/bois bois/panneau


    def _FvRk_BoisBois_Johansen(self):
        """Calcul la capacité résistante en cisaillement de la tige en N par plan de cisaillement avec
            t1 : valeur minimale entre epaisseur de l'élément bois latéral et la profondeur de pénétration en mm
            t2 : epaisseur de l'élément bois central en mm """
            
        beta = self.fh2k/self.fh1k

        if self.nCis == 1:

            a = self.fh1k * self.t1 * self.d
            b = self.fh2k * self.t2 * self.d
            c = (a/(1 + beta)) * (mt.sqrt(beta + 2 * beta**2 * (1 + (self.t2 / self.t1) + (self.t2 / self.t1)**2) + beta**3 *
                                          (self.t2 / self.t1)**2) - beta * (1 + (self.t2 / self.t1)))
            d = 1.05 * (a / (2 + beta)) * (mt.sqrt(2 * beta * (1 + beta) + ((4 * beta * (2 + beta) * self.MyRk)/(self.fh1k *
                                                                                                                 self.t1**2 * self.d))) - beta)
            e = 1.05 * ((self.fh1k * self.t2 * self.d) / (1 + 2 * beta)) * (mt.sqrt(2 * beta**2 * (1 + beta) + ((4 * beta * (1 + 2 * beta) * self.MyRk)/(self.fh1k *
                                                                                                                                                            self.t2**2 * self.d))) - beta)
            f = 1.15 * mt.sqrt((2 * beta)/(1 + beta)) * \
                mt.sqrt(2 * self.MyRk * self.fh1k * self.d)

            modeVal = min(a, b, c, d, e, f)
            dicoRupture = {a: "A", b: "B", c: "C", d: "D", e: "E", f: "F"}
            self.mode_rupture = dicoRupture[modeVal]

        else:

            g = self.fh1k * self.t1 * self.d
            h = 0.5 * self.fh2k * self.t2 * self.d
            j = 1.05 * (g / (2 + beta)) * (mt.sqrt(2 * beta * (1 + beta) + ((4 * beta * (2 + beta) * self.MyRk)/(self.fh1k *
                                                                                                                 self.t1**2 * self.d))) - beta)
            k = 1.15 * mt.sqrt((2 * beta)/(1 + beta)) * \
                mt.sqrt(2 * self.MyRk * self.fh1k * self.d)

            modeVal = min(g, h, j, k)
            dicoRupture = {g: "G", h: "H", j: "J", k: "K"}
            self.mode_rupture = dicoRupture[modeVal]
        return modeVal, self.mode_rupture


    def _FvRk_BoisBois_Tot(self):
        """Calcul la capacité résistante en cisaillement caractéristique avec la partie de Johansen + l'effet de corde si il existe dans le
        mode de rupture avec :
            Fvrk : résistance en cisaillement (uniquement partie de Johansen !) de l'organe en N par plan"
            FaxRk : résistance à l'arrachement de l'organe en N
            mode : mode de rupture de l'organe "A","B","C","D","E","F","G","H","I","J","K","L" ou "M"
        """
        coeflimit = __class__.DICO_COEF_LIMITE.get(self.type_organe, 0)
        
        if self.mode_rupture in ["C", "D", "E", "F", "J", "K"]:

            self.FaxRk_reel = min((self.FaxRk / 4), (coeflimit * self.FvRk_Johansen[0]))
            fvrk = self.FvRk_Johansen[0] + self.FaxRk_reel

        else:
            fvrk = self.FvRk_Johansen[0]
        return fvrk


    # 8.2.3 Assemblage bois métal

    def _FvRk_BoisMetal_Johansen(self):
        """Calcul la capacité résistante en cisaillement de la tige en N par plan de cisaillement avec
            t1 : valeur minimale entre epaisseur de l'élément bois latéral et la profondeur de pénétration en mm
            t2 : epaisseur de l'élément bois central en mm
            """

        if self._type_beam[0] == "Métal":
            self.t = self.beam_1.t
            if self.nCis == 2:
                self.pos_plaq = "externe"
        else:
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

            a = 0.4 * self.fh1k * self.t1 * self.d
            b = 1.15 * mt.sqrt(2 * self.MyRk * self.fh1k * self.d)

            modeVal = min(a, b)
            dicoRupture = {a: "A", b: "B"}
            self.mode_rupture = dicoRupture[modeVal]

        elif self.type_plaque == "epaisse" and self.nCis == 1:

            c = self.fh1k * self.t1 * self.d
            d = c * (mt.sqrt(2 + (4 * self.MyRk) /
                     (self.fh1k * self.d * self.t1 ** 2)) - 1)
            e = 2.3 * mt.sqrt(self.MyRk * self.fh1k * self.d)

            modeVal = min(c, d, e)
            dicoRupture = {c: "C", d: "D", e: "E"}
            self.mode_rupture = dicoRupture[modeVal]

        elif self.nCis == 2 and self.pos_plaq == "centrale":

            f = self.fh1k * self.t1 * self.d
            g = f * (mt.sqrt(2 + (4 * self.MyRk) /
                     (self.fh1k * self.d * self.t1 ** 2)) - 1)
            h = 2.3 * mt.sqrt(self.MyRk * self.fh1k * self.d)

            modeVal = min(f, g, h)
            dicoRupture = {f: "F", g: "G", h: "H"}
            self.mode_rupture = dicoRupture[modeVal]

        elif self.type_plaque == "mince" and self.nCis == 2 and self.pos_plaq != "centrale":

            j = 0.5 * self.fh2k * self.t2 * self.d
            k = 1.15 * mt.sqrt(2 * self.MyRk * self.fh2k * self.d)

            modeVal = min(j, k)
            dicoRupture = {j: "J", k: "K"}
            self.mode_rupture = dicoRupture[modeVal]

        else:

            lmode = 0.5 * self.fh2k * self.t2 * self.d
            m = 2.3 * mt.sqrt(self.MyRk * self.fh2k * self.d)

            modeVal = min(lmode, m)
            dicoRupture = {lmode: "L", m: "M"}
            self.mode_rupture = dicoRupture[modeVal]

        return modeVal, self.mode_rupture


    def _FvRk_BoisMetal_Tot(self):
        """Calcul la capacité résistante en cisaillement caractéristique avec la partie de Johansen + l'effet de corde si il existe dans le
        mode de rupture avec"""
        coeflimit = __class__.DICO_COEF_LIMITE.get(self.type_organe, 0)

        if self.mode_rupture in ["B", "D", "E", "G", "H", "K", "M"]:

            self.FaxRk_reel = min((self.FaxRk / 4), (coeflimit * self.FvRk_Johansen[0]))
            fvrk = self.FvRk_Johansen[0] + self.FaxRk_reel
            
        else:
            fvrk = self.FvRk_Johansen[0]
        return fvrk
    
    
    # Annexe A : Cisaillement de bloc
    def FbsRk(self, dp, t1, a1, a2, a3t, mode, fhk, ft0k, fvk, kcr=0.67):
        """ Calcul la valeur caractéristique en cisaillement de bloc en N (la valeur max entre FbsRkT et FbsRkV est à prendre) avec:
            dp : diamètre de perçage en mm
            t1 : épaisseur du bois en mm
            a1 : pince longitudinale en mm
            a2 : pince perpendiculaire en mm
            a3t : pince en bord chargée suivant le file en mm
            kcr : coeff de réduction largeur en cisaillement 
            """
        lvi = a1 - self.d
        lnetV = 2 * (lvi * (self.n - 1) + a3t - (dp * 0.5))
        
        lti = a2 - dp
        lnetT = lti * (self.nfile - 1)
        
        if mode == "C" or mode == "F" or mode == "J" or mode == "L" or mode == "K" or mode == "M":
            anetV = lnetV * (t1 * kcr)
        else:
            dicotef = {"A": 0.4*t1,
                       "B": 1.4 * mt.sqrt(self.MyRk/(fhk*self.d)),
                       "D": t1 * (mt.sqrt(2 + (4*self.MyRk)/(fhk * self.d * t1**2)-1)),
                       "E": 2 * mt.sqrt(self.MyRk/(fhk*self.d)),
                       "G": t1 * (mt.sqrt(2 + (4*self.MyRk)/(fhk * self.d * t1**2)-1)),
                       "H": 2 * mt.sqrt(self.MyRk/(fhk*self.d))}
            anetV = (lnetV/2) * (lnetT + 2 * (dicotef[mode]*kcr))
            
        anetT = lnetT * t1
        
        fbsrkT = 1.5 * anetT * ft0k
        fbsrkV = 0.7 * anetV * fvk
        
        return fbsrkT, fbsrkV
    

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

    TYPE_ASSEMBLAGE = ["Bois/Bois", ["CP", "Panneau dur", "PP/OSB"], "Bois/Métal"]

    def __init__(self, d: float|int, fu: int, n: int, alpha: float, pk: int=350, type_assemblage: str=TYPE_ASSEMBLAGE[0], type_organe: str=["Pointe circulaire", "Pointe carrée", "Autres pointes"], percage: bool=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.d = d
        self.n = n
        self.fu = fu
        self.alpha = mt.radians(alpha)
        self.pk = pk
        self.type_assemblage = type_assemblage
        self.percage = percage

        if type_organe == "Pointe carrée":
            self.type_organe = __class__.DICO_COEF_LIMITE[type_organe]
        else:
            self.type_organe = __class__.DICO_COEF_LIMITE[type_organe]
            self.type_circulaire = True


    @property
    def MyRk(self):
        """ Défini le moment d'écoulement plastique d'une pointe en N.mm avec:
            d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            fu : la résistance caractéristique en traction du fil d'acier en N/mm2 """
        if self.fu >= 600:
            if self.type_circulaire == True:
                MyRk = 0.3 * self.fu * self.d**2.6
            else:
                MyRk = 0.45 * self.fu * self.d**2.6
            return MyRk
        else:
            print("La résistance du fil en traction est inférieur à 600 MPa, vérifier vos données !")


    def fhk_bois(self):
        """ Calcul la portance locale des pointes inférieur à 8mm dans le bois et le LVL en MPa
            """
        if self.percage:
            fhk = 0.082 * (1 - 0.01 * self.d) * self.pk
        else:
            fhk = 0.082 * self.pk * self.d**(-0.3)
        return fhk


    def fhk_panneau(self, t: int, diam_tete_pointe: int):
        """Calcul la portance locale des pointes dans les panneaux en MPa

        Args:
            t (int): épaisseur du panneau en mm
            diam_tete_pointe (int): diamètre de la tête de la pointe

        Returns:
            float: portance locale en MPa
        """
        if diam_tete_pointe >= 2*self.d:
            if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[1][0]:
                fhk = 0.11 * self.pk * self.d**-0.3
            elif self.type_assemblage == __class__.TYPE_ASSEMBLAGE[1][1]:
                fhk = 30 * self.d**-0.3 * t**0.6
            else:
                fhk = 65 * self.d**-0.7 * t**0.1
            return fhk
        else:
            print(f"La tête de la pointe doit être au moins égale à {2*self.d} mm")


    def _kef(self, a1):
        """ coefficient donnée dans le tableau 8.1 fonction de a1 et du percage qui réduit le nombre efficace de pointe dans le sens 
            du fil avec :
                a1 : l'espacement entre tige dans le sens du fil du bois """
        listeTab = [self.d * 4, self.d * 7, self.d *
                    10, self.d * 14, ["x", 0.7, 0.85, 1]]

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
        return kef


    def nef(self, a1):
        """Défini le nombre efficace d'organe dans une file avec :
            n : nombre de boulons dans une file 
            d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            percage : True or False
            kef : coefficient donnée dans le tableau 8.1 fonction de a1 et du percage """
        self._nef = self.n**self._kef(a1)
        return self._nef


    @property
    def pince(self):
        """Défini les différentes pinces minimales pour une pointe en mm"""
        
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
            if self.pk <= 420:
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

            elif self.pk > 420 and self.pk <= 500:
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

        elif self.type_assemblage == __class__.TYPE_ASSEMBLAGE[2]: #Si assemblage bois métal
            a1 = round(a1 * 0.7, 1)
            a2 = round(a2 * 0.7, 1)

        else: #Si assemblage bois panneau
            a1 = round(a1 * 0.85, 1)
            a2 = round(a2 * 0.85, 1)
            if self.type_assemblage == __class__.TYPE_ASSEMBLAGE[1][0]:
                a3t = round((3 + 4 * mt.sin(self.alpha)) * self.d, 1)
                a3c = round(3 * self.d, 1)
                a4t = a3t
                a4c = a3c

        return {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}


    def prepercage(self, sensible: bool=False):
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
        alpha : angle entre l'effort de l'organe et le fil du bois en ° """
        
    QUALITE_ACIER = ['4.6', '4.8', '5.6', '5.8', '6.8', '8.8', '9.8', '10.9', '12.9']
    
    def __init__(self, d, qualite: float=QUALITE_ACIER, n: int=1, alpha: int|float=0, **kwargs):
        super().__init__(**kwargs)
        self.type_organe = "Boulon"
        self.d = d
        self.qualite = qualite
        self.fuk = self.__qualite_acier.loc["fub"]
        self.n = n
        self._nef = n
        self.alpha = mt.radians(alpha)
    
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
        rho_k = int(beam.caract_meca.loc["rhok"])
        return 0.082 * (1 - 0.01 * self.d) * rho_k
    
    
    def _K_90(self, beam: object):
        """Coef. modélisant la diminution de portance local quand un angle est donnée entre l'effort et le fil avec
            beam: poutre à calculer"""
            
        if beam.classe[0:1] == "C":
                type_b = "C"
        elif beam.classe[0:3] == "LVL":
            type_b = "LVL"
        else: 
            type_b == "D"
                    
        if type_b == "C" or type_b == "c":
            ck90 = 1.35 + 0.015 * self.d
        elif type_b == "LVL" or type_b == "lvl":
            ck90 = 1.30 + 0.015 * self.d
        else:
            ck90 = 0.9 + 0.015 * self.d
        return ck90
    

    def _fhak(self, fh0k, k90):
        """Calcul la portance locale d'un boulon bois/bois ou d'un tire fond si d>6mm par rapport à un effort donné à un angle
        du fil en MPa avec :
            fh0k : portance locale dans le sens du fil d'un boulon
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            k90 : coef. de réduction de la portance locale quand un effort à un angle par rapport au fil du bois"""
        
        fhak = fh0k / (k90 * mt.sin(self.alpha) ** 2 + mt.cos(self.alpha) ** 2)
        return fhak
    
    
    # 8.5.1.1 Généralité et assemblage bois/panneaux
    def _fhk(self, beam:object):
        """ Calcul la portance locale d'un boulon dans un assemblage bois/panneaux en MPa avec :
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm
            pk : masse volumique caractéristique du contreplaqué en kg/m3
            ep : epaisseur du panneau en mm
            typeP : type de panneau utilisé :  "CP" ou "Autres" """
        if beam.type_bois == "CP":
            rho_k = int(beam.caract_meca.loc["rhok"])
            fhk = 0.11 * (1 - 0.01 * self.d) * rho_k
        else:
            ep = beam.b_calcul
            fhk = 50 * (self.d**(-0.6)) * (ep**0.2)
        return fhk
    
    
    
    def fhik(self):
        """Calcul la portance locale d'un boulon bois/bois ou d'un tire fond si d>6mm
        """
    
        dict_beam = {"1": {}, "2": {}}
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if self._type_beam[i] == "Bois":
                dict_beam[str(i+1)]["fh0k"] = self._fh0k(beam)
                dict_beam[str(i+1)]["K90"] = self._K_90(beam)
                dict_beam[str(i+1)]["fhik"] = self._fhak(dict_beam[str(i+1)]["fh0k"], dict_beam[str(i+1)]["K90"])
                
            elif self._type_beam[i] == "CP" or self._type_beam[i] == "PP/OSB":
                dict_beam[str(i+1)]["fhik"] = self._fhk(beam)
                
            else:
                dict_beam[str(i+1)]["fhik"] = 0
        
            if self._type_beam[i] != "Métal":
                if i:
                    self.t2 = self.beam_2.b_calcul
                else:
                    self.t1 = self.beam_1.b_calcul
        
        
        self.fh1k = dict_beam["1"]["fhik"]
        self.fh2k = dict_beam["2"]["fhik"]
        return self.fh1k, self.fh2k


    @property
    def MyRk(self):
        """Défini le moment d'écoulement plastique d'un boulon en N.mm avec:
            fuk : la valeur caractéristique de résistance à la traction du boulon en N/mm2
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm"""
        return 0.3 * self.fuk * self.d ** 2.6


    @property
    def pinceBoulon(self):
        """Défini les différentes pinces minimales pour un boulon en mm avec :
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm """
        a1 = round((4 + mt.cos(self.alpha)) * self.d, 1)
        a2 = round(4 * self.d, 1)
        a3t = max(7 * self.d, 80)

        if self.alpha <= 150 and self.alpha < 210:
            a3c = round(4 * self.d, 1)
        else:
            a3c = round((1 + 6 * mt.sin(self.alpha)) * self.d, 1)

        a4t = round(max((2 + 2 * mt.sin(self.alpha))
                    * self.d, 3 * self.d), 1)
        a4c = round(3 * self.d, 1)
        return {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}

    def nef(self, a1):
        """Défini le nombre efficace d'organe (boulon) dans une file avec :
            a1 : l'espacement entre boulon dans le sens du fil du bois
            n : nombre de boulons dans une file 
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm"""
        self._nef = min(self.n**0.9 * (a1/(13 * self.d))**(1/4), self.n)
        return self._nef

    # 8.5.2 Boulons chargés axialement
    @property
    def FaxRk(self):
        return 0

# ======================================================= BROCHE =========================================================
# 8.6 Assemblage par broche

class Broche(Boulon):
    """ Défini une objet broche avec : 
        d : diamètre efficace de la broche ( entre 6 et 30 mm) en  mm
        fuk : la valeur caractéristique de résistance à la traction de la broche en N/mm2
        n : nombre de broche dans une file
        alpha : angle entre l'effort de l'organe et le fil du bois en ° """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_organe = "Broche"


    @property
    def pinceBroche(self):
        """Défini les différentes pinces minimales pour une broche en mm avec :
            alpha : angle entre l'effort de l'organe et le fil du bois en °
            d : diamètre efficace du boulon (ou du tire fond si >6mm) en  mm """
        a1 = round((3 + 2 * (mt.cos(self.alpha))) * self.d, 1)
        a2 = round(3 * self.d, 1)
        a3t = round(max(7 * self.d, 80), 1)

        if self.alpha <= 150 and self.alpha < 210:

            a3c = round(max(3.5 * self.d, 40), 1)

        else:

            a3c = round((a3t * mt.sin(self.alpha)), 1)

        a4t = round(max((2 + 2 * mt.sin(self.alpha))
                    * self.d, 3 * self.d), 1)
        a4c = round(3 * self.d, 1)
        return {"a1": a1, "a2":a2, "a3t": a3t, "a3c": a3c, "a4t": a4t, "a4c": a4c}


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

    def __init__(self, d, d1, dh,  pa, fhead, ftensk, n, alphaTirefond=90):
        self.type_organe = "Tirefond"
        self.d = d
        self.d1 = d1
        self.dh = dh
        self.pa = pa
        self.fhead = fhead
        self.ftensk = ftensk
        self.n = n
        self.alphaTirefond = mt.radians(alphaTirefond)


    # 8.7.2 Tirefond chargés axialement
    def pinceTirefondAxial(self, t):
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


    def faxk(self, lef, pk=350):
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


    def FaxaRk(self, faxk, lef, pk=350):
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


    def FaxaRkHead(self, pk):
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

    def __init__(self, dc, t1, t2, hc, typeA="bois"):
        self.type_organe = "Anneau"
        self.dc = dc
        self.t1 = t1
        self.t2 = t2
        self.he = hc/2
        self.typeA = typeA

    def ki(self, nAss, a3t, rhok):
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

    def Fv0Rk(self, dicoKi):
        """ Retourne la résistance en cidaillement de l'anneau en N avec:
            dicoKi = dictionnaire des facteurs ki (def ki)"""
        fv0rk = min(dicoKi["k1"] * dicoKi["k2"] * dicoKi["k3"] * dicoKi["k4"] * (35 * self.dc**1.5),
                    dicoKi["k1"] * dicoKi["k3"] * self.he * (31.5 * self.dc))
        return fv0rk
    

    def FvaRk(self):
        pass


if __name__ == "__main__":
    from EC5_Element_droit import Beam
    from EC3_Element_droit import Element
    
    # beam1 = Beam(60, 200, "Rectangulaire")
    # beam1.f_type_d()
    beam1 = Element(6, 200)
    beam2 = Beam(60, 200, "Rectangulaire", classe="C24")
    beam2.f_type_d()
    ass= Assemblage(beam1, beam2, nfile=1, nCis=2)
    bl = Boulon._from_parent_class(ass, d=16, qualite=8.8, n=2)
    print(bl.fhik())
    print(bl.pinceBoulon)
    print(bl.nef(100))
    print(bl.FvRd(effet_corde=True))
    print(bl.Kser())