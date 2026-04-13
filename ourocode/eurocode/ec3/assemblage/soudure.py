# coding in UTF-8
import os, sys
from PIL import Image
from math import *
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.ec3.element_droit.plat import Plat


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

