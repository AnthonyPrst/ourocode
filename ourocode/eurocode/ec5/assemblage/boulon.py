# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec5.assemblage.assemblage import Assemblage

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
    def MyRk(self) -> tuple:
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
    def pince(self) -> dict:
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
    
        
    def nef(self, a1_beam1:int, a1_beam2:int) -> tuple:
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
    def Fax_Rk(self, d_int: float=0, d_ext: float=0, filetage_EN1090: bool=("True", "False")) -> tuple:
        """Calcul la résistance axial caractéristique d'un boulon chargé axialement à partir soit de la rondelle soit de la plaque métalique.

        Args:
            d_int (float, optional): diamètre intérieur de la rondelle en mm ou du trou de perçage dans la plaque métallique. Defaults to 0.
            d_ext (float, optional): diamètre extérieur de la rondelle en mm. Defaults to 0.
            filetage_EN1090 (bool, optional): définit si le filetage est conforme à l'EN 1090, soit matricé. Si filetage usiné alors False. Defaults to True.

        Returns:
            FaxRk: la résistance axial d'un boulon en N
        """
        from ourocode.eurocode.ec3.assemblage.tige import Tige
        d_int = d_int * si.mm
        d_ext = d_ext * si.mm

        if self._type_beam[0] == "Métal":
            FtRd = Tige(
                d=self.d.value*10**3, 
                d0=d_int, 
                trou_oblong=False, 
                qualite=self.qualite, 
                tete_fraisee=False, 
                prof_tete_fraisee=0, 
                verif_filetage=True, 
                filetage_EN1090=filetage_EN1090, 
                t=self.beam_1.t.value*10**3, 
                h=self.beam_1.h.value*10**3, 
                classe_acier=self.beam_1.classe_acier, 
                classe_transv=self.beam_1.classe_transv).FtRd
            fc_90_k = float(self.beam_2.caract_meca.loc["fc90k"]) * si.MPa
            if self.nCis == 1:
                d_ext = min(self.beam_1.t*12, 4*self.d, d_ext)
            else:
                d_ext = min(self.beam_1.t*12, 4*self.d)
        else:
            FtRd = Tige(
                d=self.d.value*10**3, 
                d0=d_int, 
                trou_oblong=False, 
                qualite=self.qualite, 
                tete_fraisee=False, 
                prof_tete_fraisee=0, 
                verif_filetage=True, 
                filetage_EN1090=filetage_EN1090, 
                t=0, 
                h=0, 
                classe_acier="S235", 
                classe_transv=3).FtRd
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

