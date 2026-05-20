# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core._renderer import handcalc

from ourocode.eurocode.ec5.assemblage.assemblage import Assemblage, interpolationLineaire

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
            raise ValueError("Erreur, le diamètre de la pointe est supérieur à 6mm, le prépercage est obligatoire")
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
    def MyRk(self) -> tuple:
        """ Défini le moment d'écoulement plastique d'une pointe en N.mm avec:
            d : diamètre de la pointe en mm (pour les pointe carrée = coté de la pointe)
            fu : la résistance caractéristique en traction du fil d'acier en N/mm2 """
        f_u = self.fu.value * 10**-6
        d = self.d.value * 10**3
        if self.fu >= 600:
            if self._type_circulaire:
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
    
    
    def Fax_Rk(self) -> tuple:
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

                    if kef == "x" and not self.percage:
                        kef = 0

                    elif kef == "x" and self.percage:
                        kef = 0.5
                j += 1

            else:
                if listeTab[4][j-1] == "x" and not self.percage:
                    kef = interpolationLineaire(
                        a1, listeTab[j-1], listeTab[j], 0, listeTab[4][j])

                elif listeTab[4][j-1] == "x" and self.percage:
                    kef = interpolationLineaire(
                        a1, listeTab[j-1], listeTab[j], 0.5, listeTab[4][j])

                else:
                    kef = interpolationLineaire(
                        a1, listeTab[j-1], listeTab[j], listeTab[4][j-1], listeTab[4][j])
        self.kef = kef
        return self.kef


    def nef(self, a1_beam1:int, a1_beam2:int) -> tuple:
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
    def pince(self) -> dict:
        """Retourne les différentes pinces minimales pour une pointe en mm."""
        dict_pince = {}
        if self.type_organe == "Tirefond":
            self.d = self.d_vis
        for i, beam in enumerate([self.beam_1, self.beam_2]):
            if self._type_beam[i] not in self.TYPE_BOIS_ASSEMBLAGE:
                continue
            
            alpha = self.alpha[i]
            rho_k = beam.rho_k
            if self.percage:
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


    def prepercage(self, beam: str=["1", "2"], sensible: bool=("False", "True")) -> tuple:
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

