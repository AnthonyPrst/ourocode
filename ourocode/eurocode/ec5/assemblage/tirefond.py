# coding in UTF-8
# by Anthony PARISOT
from copy import deepcopy
import warnings
import math as mt
from math import sin, cos, radians, sqrt, pi

import forallpeople as si
si.environment("structural")
from ourocode.eurocode.core.renderer import handcalc

from ourocode.eurocode.ec5.assemblage.pointe import Pointe
from ourocode.eurocode.ec5.assemblage.boulon import Boulon



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
    def pince_tirefond_axial(self, t: int) -> dict:
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
            warnings.warn("L'épaisseur de bois n'est pas suffisante, il faut un bois de {0} mm minimum !".format(
                12*self.d_vis))
        return {"a1": a1, "a2": a2, "a1CG": a1CG, "a2CG": a2CG}


    @property
    def nefTraction(self) -> tuple:
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


    def faxk(self, l_ef:int, beam:str=["1", "2"]) -> tuple:
        """Calcule la résistance caractéristique à l'arrachement f_ax,k selon EN 1995-1-1 §8.7.2(4).

        Applicable uniquement pour 6 mm ≤ d ≤ 12 mm et 0.6 ≤ d1/d ≤ 0.75.
        Formule : f_ax,k = 0.52 × d⁻°·⁵ × l_ef⁻°·¹ × ρ_k°·⁸ en N/mm².

        Args:
            l_ef (int): Longueur de pénétration de la partie filentée dans l'élément bois, en mm.
            beam (str): Numéro de l'élément bois à considérer ("1" ou "2"). Defaults to "1".

        Returns:
            tuple: (latex_string, f_ax_k) où f_ax_k est la résistance à l'arrachement en N/mm²
                (valeur scalaire, sans unité forallpeople).
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
            warnings.warn(
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
    
    
    def Fax_Rk(self, faxk:float, l_ef:int, alpha:int) -> tuple:
        """Calcule la résistance caractéristique axiale de l'assemblage F_ax,Rk,ass selon EN 1995-1-1 §8.7.2.

        Détermine min(F_ax,alpha,Rk, F_head,Rk) pour chaque tirefond, puis applique n_ef,traction
        pour obtenir la résistance de l'assemblage. Stocke les résultats dans self.FaxRk et self.Fax_Rk_ass.

        Args:
            faxk (float): Résistance caractéristique à l'arrachement perpendiculaire au fil en N/mm².
            l_ef (int): Longueur de pénétration de la partie filentée dans l'élément bois, en mm.
            alpha (int): Angle entre l'axe du tirefond et le fil du bois en ° (doit être ≥ 30°).

        Returns:
            tuple: (latex_string, F_ax_Rk_ass) où F_ax_Rk_ass est la résistance axiale caractéristique
                totale de l'assemblage en N (avec unité si.N).
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


    def FtRk(self) -> tuple:
        """Calcule la résistance caractéristique en traction nette du fil des tirefonds de l'assemblage.

        Formule : F_t,Rk = n_ef,traction × f_tens,k en N.

        Returns:
            tuple: (latex_string, F_t_Rk) où F_t_Rk est la résistance en traction en N (avec unité si.N).
        """
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
    def MyRk(self) -> tuple:
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
    def MyRk(self) -> tuple:
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



