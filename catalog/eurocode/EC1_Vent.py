#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import os
import sys
from PIL import Image

import math as mt
from math import *
from copy import copy

import pandas as pd
import forallpeople as si
from handcalcs.decorator import handcalc

sys.path.append(os.path.join(os.getcwd(),"catalog", "eurocode"))
from A0_Projet import Batiment



si.environment("structural")


def interpolation_lineaire(x, xa, xb, ya, yb):
	"""Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
	y = ya + (x - xa) * ((yb - ya)/(xb - xa))
	return y

def interpolation_logarithmique(x, xa, xb, ya, yb):
	"""Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
	y = (x > xb) * yb + (x < xa) * ya + (x <= xb) * (x >= 1) * (ya - (ya - yb) * mt.log10(x))
	return y


class Vent(Batiment):
	"""Créer une classe qui permet de calculer l'action de vent.
		Cette classe est hérité de la classe Batiment du module A0_Projet.py.
	"""
	RHO_AIR = 1.225 * si.kg/si.m**3
	CPI = [0.2, -0.3]
	VB0 = {"1": 22, "2": 24, "3": 26, "4": 28, 
		   "Guadeloupe": 36, 
		   "Guyane": 17, 
		   "Martinique": 32, 
		   "Mayotte": 30, 
		   "Réunion": 34}

	CAT_TERRAIN = {"0": {"Z0": 0.005, "Zmin": 1},
				   "II": {"Z0": 0.05, "Zmin": 2},
				   "IIIa": {"Z0": 0.2, "Zmin": 5},
				   "IIIb": {"Z0": 0.5, "Zmin": 9},
				   "IV": {"Z0": 1, "Zmin": 15}}

	CAT_ORO = {"Aucun": 1, "Cas 1": "1", "Cas 2": "2"}
	CFR = {"Lisse": 0.01, "Rugueuse": 0.02, "Très rugueuse": 0.04}

	def __init__(self, z:float, terrain:str=CAT_TERRAIN, oro: str=CAT_ORO, CsCd: float= 1, **kwargs):
		"""
		Args:
			z (float): hauteur en m sur le bâtiment ou est étudié le vent (Ze suivant EN 1991-1-4 §7.2.2).
			alt (int): altitude du batiment étudié en m.
			cat_terrain (str): Catégorie de terrain du projet.
			cat_oro (str): Catégorie orographique. Default to "Aucun".
		"""
		super().__init__(**kwargs)
		self.terrain = terrain
		self.oro = oro
		self.z = z *si.m
		self.CsCd = CsCd


	@property
	def cat_terrain(self):
		return __class__.CAT_TERRAIN[self.terrain]

	@property
	def cat_oro(self):
		return __class__.CAT_ORO[self.oro]

	@property
	def zone_vent(self):
		file = "carte_action_region.csv"
		df = self._data_from_csv(file, index_col=1)
		return int(df.loc[str(self.code_INSEE)]["Zone_vent"])

	@property
	def Vb_0(self):
		"""Vb,0 est la valeur de base de la vitesse de référence du vent en m/s

		Returns:
			int: vitesse de référence du vent en m/s
		"""
		return self.VB0[str(self.zone_vent)] * si.m/si.s


	@property
	def Vb(self, Cdir=1, Cseason=1):
		"""Retourne la vitesse de référence du vent en m/s fonction de cette dernière et de la période de l'année
		à une hauteur de 10 m au-dessus d'un sol relevant de la catégorie de terrain II 

		Args:
			Cdir (float): coefficient de direction. Defaults to 1.
			Cseason (float): coefficient de saison. Defaults to 1.

		Returns:
			float: vitesse de référence du vent en m/s
		"""
		V_b_0 = self.Vb_0
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			V_b = V_b_0 * Cdir * Cseason
			return V_b
		return val()


	@property
	def rayon_secteur_angu(self):
		"""Retourne le rayon du secteur angulaire dans lequel la rugosité de terrain est à qualifier fonction de la hauteur du bâtiment étudier

		Returns:
			float: max entre 300m et R
		"""
		z = self.z.value
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			rayon_secteur_angu = max(23 * z**1.2, 300) * si.m
			return rayon_secteur_angu
		return val()


	@property
	def _Kr(self)->float:
		z0II = self.CAT_TERRAIN["II"]["Z0"]
		return 0.19 * (self.cat_terrain["Z0"] / z0II)**0.07


	@property
	def _Cr_z(self):
		"""Retourne le coefficient de rugosité, c r ( z ), tient compte de la variabilité de la vitesse moyenne du vent sur le site de la
		construction 

		Returns:
			float: le coef de rugosité fonction de la hauteur au niveau du sol et de la rugosité de terrain
		"""
		
		if self.cat_terrain["Zmin"] <= self.z.value <= 200:
			Cr_z = self._Kr * mt.log(self.z.value/self.cat_terrain["Z0"])
		else:
			Cr_z = self._Kr * mt.log(self.cat_terrain["Zmin"]/self.cat_terrain["Z0"])

		# print("Kr :", self._Kr)
		# print("Cr_z :", Cr_z)
		return Cr_z


	@property
	def Co_z(self):
		"""Retourne le coeficient orographique du terrain suivant la catégorie orographique et la hauteur z du bâtiment

		Returns:
			float: coef orographique
		"""
		# ATTENTION vérifier le z et zmin si toujours comptatible quand ajout du cas 2 !!!
		if self.z.value < self.cat_terrain["Zmin"]:
			z = self.cat_terrain["Zmin"]
		else:
			z = self.z.value

		match self.cat_oro:
			case "1":
				dict_alti = {"An1": "", "An2": "", "Ae1": "", "Ae2": "",
							 "As1": "", "As2": "", "Ao1": "", "Ao2": ""}

				for cle, value in dict_alti.items():
					dict_alti[cle] = int(input("Altitude {0} en m: ".format(cle)))

				Am = (2*self.alt + dict_alti["An1"] + dict_alti["An2"] + dict_alti["Ae1"] + dict_alti["Ae2"] + 
					 dict_alti["As1"] + dict_alti["As2"] + dict_alti["Ao1"] + dict_alti["Ao2"])/10

				delta_AC = self.alt - Am

				if self.z.value >= 10:
					C0_z = 1 + 0.004 * delta_AC*mt.exp(z-10)
				else:
					C0_z = 1 + 0.004 * delta_AC*mt.exp(0)

			case "2":
				pass
			case _:
				C0_z = 1
		# print("C0_z :", C0_z)
		return max(C0_z, 1)


	@property
	# @handcalc(override="params", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
	def Vm_z(self):
		"""Retourne la vitesse moyenne du vent à une hauteur z au dessus du sol en fonction de la rugosité de terrain,
		   de l'orographie et de la vitesse de référence du vent Vb

		Returns:
			float: vitesse moyenne du vent en m/s
		"""
		C_r_z = self._Cr_z
		C_o_z = self.Co_z
		V_b = self.Vb[1]
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			V_m_z = C_r_z * C_o_z * V_b
			return V_m_z
		return val()
	

	@property
	def _Kl(self):
		"""Retourne le coefficient de turbulence

		Returns:
			float: coef de turbulence
		"""
		match self.cat_oro:
			case "1":
				return self.Co_z * (1 - 2*10**-4 *(mt.log10(self.cat_terrain["Z0"])+3)**6) #4.19NA
			case _:
				return 1 - 2*10**-4 *(mt.log10(self.cat_terrain["Z0"])+3)**6 #4.20NA
		   

	@property
	def _sigma_v(self):
		"""Retourne l'écart type de turbulence

		Returns:
			float: écart type
		"""
		return self._Kr * self.Vb[1] * self._Kl #4.6


	@property
	def _Iv_z(self):
		"""Retourne l'intensité de turbulence à la hauteur z

		Returns:
			float: intensité de la turbulence
		"""
		if self.cat_terrain["Zmin"] <= self.z.value <= 200 or self.z.value < self.cat_terrain["Zmin"] :
			return self._sigma_v / self.Vm_z[1] #4.7
		else:
			print("""Erreur la hauteur max du bâtiment ne peut dépasser 200m, 
					 impossible de calculer l'intensité de turbulence Iv""")


	@property
	def _Qb(self):
		"""Retour la pression dynamique de référence (informatif)

		Returns:
			float: pression en N/m²
		"""
		return 0.5 * self.RHO_AIR * self.Vb[1]**2 #4.10


	@property
	def Qp_z(self):
		"""Retour la pression dynamique de pointe à la hauteur z qui est induite par la vitesse moyenne 
		et les fluctuation rapides de vitesse

		Returns:
			float: pression en N/m²
		"""
		I_v_z = self._Iv_z
		rho_air = self.RHO_AIR
		V_m_z = self.Vm_z[1]
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			Q_p_z = (1 + 7 * I_v_z) * 0.5 * rho_air * V_m_z**2 #4.8
			return Q_p_z
		return val()


	@property
	def _Ce_z(self):
		"""Retourne le coefficient d_bat'exposition à la hauteur z (informatif)

		Returns:
			float: coef d'expo
		"""
		return self.Qp_z[1]/self._Qb #4.9
	

	def We(self, Cpe: float):
		"""Retourne la pression aérodynamique sur les surfaces extérieures

		Args:
			Cpe (float): coefficient de pression extérieure

		Returns:
			float: pression en N/m²
		"""
		Q_p_z = self.Qp_z[1]
		C_pe = Cpe
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			W_e = Q_p_z * C_pe #5.1
			return W_e
		return val()
	

	def Wi(self, Cpi: float):
		"""Retourne la pression aérodynamique sur les surfaces intérieures

		Args:
			Cpi (float): coefficient de pression intérieure

		Returns:
			float: pression en N/m²"""
		Q_p_z = self.Qp_z[1]
		C_pi = Cpi
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			W_e = Q_p_z * C_pi #5.2
			return W_e
		return val()
	

	def Fw_e(self, Aref: float):
		self.CsCd
		pass #5.5

	def Fw_i(self, Aref: float):
		pass #5.6

	
	def Ffr(self, Afr: float, Cfr: str=CFR):
		"""Retourne les forces de frottement sur le bâtiment en N selon l'EN 1991-1-4 §5.7/7.5.

		Args:
			Afr (float): aire de référence en m² correspondant à la surface d'application des forces de frottement.
						Il convient d'appliquer les forces de frottement sur la partie des surfaces extérieures
						parallèle au vent, située au-delà d'une certaine distance des bords au vent ou des angles 
						au vent de la toiture, distance égale à la plus petite valeur de 2 · b ou 4 · h .

			Cfr (str): type de surface de frottement
										- Lisse / exemple: acier, béton lisse
										- Rugueuse / exemple: béton brut, bardeaux bitumés(shingles)
										- Très rugueuse / exemple: ondulation, nervures, pliures

		Returns:
			float: effort en N
		"""
		A_fr = Afr * si.m**2
		C_fr = self.CFR[Cfr]
		Q_p_z = self.Qp_z[1]

		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			F_fr = C_fr * Q_p_z * A_fr
			return F_fr
		return val()
	
	

	def show_Ffr(self):
		"""Affiche l'image du zonage pour un effet d'entrainement en toiture
		"""
		file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Ffr.png")
		image = Image.open(file)
		image.show()


	def K_red_U(self):
		"""Retourne le coefficient de défaut de corrélation entre les pressions aérodynamiques au vent et sous le vent,
		uniquement pour les faces D et E. §7.2.2(3)

		Args:
			h_bat (float): hauteur du bâtiment en m
			d_bat (float): longeur parallèle à la direction du vent en m
		"""
		if self.h_bat/self.d_bat >= 5:
			self.Kred_U = 1
		elif self.h_bat/self.d_bat <= 1:
			self.Kred_U = 0.85
		else:
			self.Kred_U = 0.85 + 0.15 * (((self.h_bat/self.d_bat)-1)/4)
			
		

class Murs_verticaux(Vent):
	def __init__(self, load_area:float, *args, **kwargs):
		"""Créer une classe permettant le calcul des voiles verticaux au vent selon l'EN 1991-1-4 §7.2.2

		Args:
			load_area (float): aire chargée pour le calcul des éléments ou des fixations
		"""
		super().__init__(*args, **kwargs)
		self.load_area = load_area
		self.wind_direction = {"0°":{}, "90°":{}}
		
		for cle in self.wind_direction.keys():
			match cle:
				case "0°":
					self.e = min(self.b_bat, self.h_bat*2)
				case "90°":
					self.e = min(self.d_bat, self.h_bat*2)
					d_bat = copy(self.d_bat)
					b_bat = copy(self.b_bat)
					self.b_bat = d_bat
					self.d_bat = b_bat
		
			self.wind_direction[cle]["geometrie"] = self._geo()
			self.wind_direction[cle]["Cpe"] = self._Cpe(cle)
			

	def _geo(self):
		"""Calcul les dimensions du zonage des parois verticale
		"""
		self._df = self._data_from_csv(os.path.join("vent", "vent_Cpe_mur_verticaux.csv"))
		self._df.reset_index(drop=False, inplace=True)

		if self.e < self.d_bat:  
			geometrie = {"A": self.e/5,
						 "B": 0.8 * self.e,
						 "C": self.d_bat - self.e}

		elif self.e >= 5*self.d_bat:
			geometrie = {"A": self.d_bat}
		else:
			geometrie = {"A": self.e/5,
						 "B": self.d_bat - self.e/5}
		geometrie["D"] = self.b_bat
		geometrie["E"] = self.b_bat

		return geometrie


	def _Cpe(self, cle:str):
		"""Calcul les Cpe des parois verticales

		Args:
			cle (str): sens du vent sur le bâtiment 0° ou 90°
		"""
		h_d = self.h_bat / self.d_bat
		list_CPE_line = []
		borne_inf = None
		if h_d > 0.25 and h_d < 1:
			df_index = [3,4]
			borne_inf = 0.25
			borne_sup = 1
			
		elif h_d > 1 and h_d < 5:
			df_index = [2,3]
			borne_inf = 1
			borne_sup = 5

		elif h_d <= 0.25:
			self._df = self._df[["Zone", "Cpe", "0.25"]]

		elif h_d >= 5:
			self._df = self._df[["Zone", "Cpe", "5"]]

		else:
			self._df = self._df[["Zone", "Cpe", "1"]]


		if borne_inf:
			for i in range(self._df.shape[0]):
				list_CPE_line.append(interpolation_lineaire(h_d, borne_inf, borne_sup, self._df.iloc[i,df_index[1]], self._df.iloc[i,df_index[0]]))
			self._df[round(h_d, 3)] = list_CPE_line
			self._df = self._df[["Zone", "Cpe", round(h_d, 3)]]

		if self.load_area > 1 and self.load_area < 10:
				for i in range(self._df.shape[0]):
					if i % 2 == 0:
						self._df.loc[self._df.shape[0]]= [self._df.iloc[i,0], "CPE "+ str(round(self.load_area, 2)), interpolation_logarithmique(self.load_area, 1, 10, self._df.iloc[i+1,2], self._df.iloc[i,2])]

		self._df.set_index("Zone", inplace=True)
		self._df = self._df.loc[[zone for zone in self.wind_direction[cle]["geometrie"].keys()]]
		return self._df


	def get_wind_dict(self) -> dict:
		return self.wind_direction


	def get_Cpe(self, dir=("0°", "90°")):
		return self.wind_direction[dir]["Cpe"]


	def show_zonage(self):
		"""Affiche l'image du zonage pour les murs verticaux
		"""
		file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cpe_mur_verticaux.png")
		image = Image.open(file)
		image.show()




class Toiture_2pants(Vent):
	"""Créer une classe permetant le calcul d'une toiture à deux versant symétrique au vent selon l'EN 1991-1-4 §7.2.5
	"""
	def __init__(self, load_area: float, *args, **kwargs):
		super().__init__(*args, **kwargs)
		


class Toiture_isolee_1pant(Vent):
	def __init__(self, phi: float, load_area: float, *args, **kwargs):
		"""Créer une classe permetant le calcul d'une toiture isolée à un versant au vent selon l'EN 1991-1-4 §7.3

		Args:
			phi (int | float): le degré d'obstruction sous une toiture isolée entre 0 et 1.
			load_area (int | float): aire chargée pour le calcul des éléments ou des fixations.
		"""
		super().__init__(*args, **kwargs)
		self.phi = phi
		self.load_area = load_area
		self.wind_direction = {"0°":{}, "90°":{}}
		
		for cle in self.wind_direction.keys():
			match cle:
				case "90°":
					d_bat = copy(self.d_bat)
					b_bat = copy(self.b_bat)
					self.b_bat = d_bat
					self.d_bat = b_bat
		
			self.wind_direction[cle]["geometrie"] = self._geo()
			self.wind_direction[cle]["Cp"] = self._Cp()
			

	def _geo(self):
		"""Calcul les surfaces du zonage de la toiture
		"""
		self._df = self._data_from_csv(os.path.join("vent", "vent_Cp_toiture_isolee_1_versant.csv"))
		self._df.reset_index(drop=False, inplace=True)
		 
		geometrie = {"A": {"lenght": [(self.b_bat - self.b_bat/10*2), (self.d_bat - self.d_bat/10*2)], "surface": (self.d_bat - self.d_bat/10*2) * (self.b_bat - self.b_bat/10*2) / mt.cos(mt.radians(self.alpha_toit))},
					"B": {"lenght": [(self.d_bat - self.d_bat/10*2), self.b_bat/10], "surface": ((self.d_bat - self.d_bat/10*2) * self.b_bat/10) / mt.cos(mt.radians(self.alpha_toit))},
					"C": {"lenght": [(self.b_bat - self.b_bat/10*2), self.d_bat/10], "surface": ((self.b_bat - self.b_bat/10*2) * self.d_bat/10) / mt.cos(mt.radians(self.alpha_toit))}}
		return geometrie


	def _Cp(self):
		"""Calcul les Cp d'une toiture isolée

		Args:
			cle (str): sens du vent sur le bâtiment 0° ou 90°
		"""
		
		if self.phi > 0 and self.phi < 1:
			for i in range(0, self._df.shape[0], 3):
				self._df.loc[self._df.shape[0]]= [self._df.iloc[i,0],
												self.phi,
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,2], self._df.iloc[i+2,2]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,3], self._df.iloc[i+2,3]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,4], self._df.iloc[i+2,4]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,5], self._df.iloc[i+2,5])
												]
			self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

		elif self.phi >= 1:
			self._df = self._df[self._df["phi"].isin(["1", "max"])]
		else:
			self._df = self._df[self._df["phi"].isin(["0", "max"])]
			
			
		
		list_alpha_toit= self._df["alpha_toit"].unique()

		if not self.alpha_toit in list_alpha_toit:
			minimum = list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
			maximum = list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
			df_min = self._df[self._df["alpha_toit"]==minimum]
			df_max = self._df[self._df["alpha_toit"]==maximum]

			self._df.reset_index(drop=True, inplace=True)
			for i, phi in enumerate(["max", self.phi]):
				self._df.loc[self._df.shape[0]]= [round(self.alpha_toit, 2),
												phi,
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,2], df_max.iloc[i,2]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,3], df_max.iloc[i,3]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,4], df_max.iloc[i,4]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,5], df_max.iloc[i,5])
												]
			self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

		self._df = self._df[self._df["alpha_toit"]==round(self.alpha_toit, 2)]
		self._df.set_index("alpha_toit", inplace=True)
		return self._df


	def get_wind_dict(self) -> dict:
		return self.wind_direction
	

	def get_geo(self, dir=("0°", "90°")):
		return self.wind_direction[dir]["geometrie"]


	def get_Cp(self, dir=("0°", "90°")):
		return self.wind_direction[dir]["Cp"]


	def show_zonage(self):
		"""Affiche l'image du zonage pour une toiture isolée un versant
		"""
		file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cp_toiture_isolee_1_versant.png")
		image = Image.open(file)
		image.show()



class Toiture_isolee_2pants(Vent):
	def __init__(self, phi: float, load_area: float, *args, **kwargs):
		"""Créer une classe permetant le calcul d_bat'une toiture isolée à deux versants au vent selon l'EN 1991-1-4 §7.3
		ATTENTION : Il ne semble pas y avoir d'inversion de zonage quand le vent est à 0° ou 90° mais une inversion des longueurs de surface A VALIDER !

		Args:
			phi (int | float): le degré d'obstruction sous une toiture isolée entre 0 et 1.
			load_area (int | float): aire chargée pour le calcul des éléments ou des fixations.
		"""
		super().__init__(*args, **kwargs)
		self.phi = phi
		self.load_area = load_area
		self.wind_direction = {"0°":{}, "90°":{}}
		
		for cle in self.wind_direction.keys():
			match cle:
				case "90°":
					d_bat = copy(self.d_bat)
					b_bat = copy(self.b_bat)
					self.b_bat = d_bat
					self.d_bat = b_bat
		
			self.wind_direction[cle]["geometrie"] = self._geo()
			self.wind_direction[cle]["Cp"] = self._Cp()
			

	def _geo(self):
		"""Calcul les surfaces du zonage de la toiture
		"""
		self._df = self._data_from_csv(os.path.join("vent", "vent_Cp_toiture_isolee_2_versants.csv"))
		self._df.reset_index(drop=False, inplace=True)
		 
		geometrie = {"A": {"lenght": [(self.d_bat - self.d_bat/10*2),  self.d_bat - (2*self.d_bat/10) - self.d_bat/5], 
							"surface": ((self.d_bat - self.d_bat/10*2) * (self.d_bat - (2*self.d_bat/10) - self.d_bat/5)) / mt.cos(mt.radians(self.alpha_toit))},
					"B": {"lenght": [(self.d_bat - self.d_bat/10*2),  self.b_bat/10], "surface": ((self.d_bat - self.d_bat/10*2) * self.b_bat/10) / mt.cos(mt.radians(self.alpha_toit))},
					"C": {"lenght": [(self.b_bat - self.b_bat/10*2), self.d_bat/10], "surface": ((self.b_bat - self.b_bat/10*2) * self.d_bat/10) / mt.cos(mt.radians(self.alpha_toit))},
					"D": {"lenght": [(self.b_bat - self.b_bat/10*2)/2, (self.d_bat/5)/2], "surface": (((self.b_bat - self.b_bat/10*2) * self.d_bat/5) / mt.cos(mt.radians(self.alpha_toit))/2)}}
		return geometrie


	def _Cp(self):
		"""Calcul les Cp d'une toiture isolée

		Args:
			cle (str): sens du vent sur le bâtiment 0° ou 90°
		"""
		
		if self.phi > 0 and self.phi < 1:
			for i in range(0, self._df.shape[0], 3):
				self._df.loc[self._df.shape[0]]= [self._df.iloc[i,0],
												self.phi,
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,2], self._df.iloc[i+2,2]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,3], self._df.iloc[i+2,3]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,4], self._df.iloc[i+2,4]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,5], self._df.iloc[i+2,5]),
												interpolation_lineaire(self.phi, 0, 1, self._df.iloc[i+1,6], self._df.iloc[i+2,6])
												]
			self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

		elif self.phi >= 1:
			self._df = self._df[self._df["phi"].isin(["1", "max"])]
		else:
			self._df = self._df[self._df["phi"].isin(["0", "max"])]
			
			
		
		list_alpha_toit= self._df["alpha_toit"].unique()

		if not self.alpha_toit in list_alpha_toit:
			minimum = list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
			maximum = list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
			df_min = self._df[self._df["alpha_toit"]==minimum]
			df_max = self._df[self._df["alpha_toit"]==maximum]

			self._df.reset_index(drop=True, inplace=True)
			for i, phi in enumerate(["max", self.phi]):
				self._df.loc[self._df.shape[0]]= [round(self.alpha_toit, 2),
												phi,
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,2], df_max.iloc[i,2]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,3], df_max.iloc[i,3]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,4], df_max.iloc[i,4]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,5], df_max.iloc[i,5]),
												interpolation_lineaire(self.alpha_toit, minimum, maximum, df_min.iloc[i,6], df_max.iloc[i,6])
												]
			self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

		self._df = self._df[self._df["alpha_toit"]==round(self.alpha_toit, 2)]
		self._df.set_index("alpha_toit", inplace=True)
		return self._df


	def get_wind_dict(self) -> dict:
		return self.wind_direction
	

	def get_Cp(self, dir=("0°", "90°")):
		return self.wind_direction[dir]["Cp"]


	def show_zonage(self):
		"""Affiche l'image du zonage pour une toiture isolée un versant
		"""
		file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cp_toiture_isolee_2_versants.png")
		image = Image.open(file)
		image.show()



if __name__ == "__main__":

	building = Batiment(h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=15, alt=400)
	Action_wind = Vent._from_parent_class(building, terrain="IIIa", oro="Aucun", z=5)
	print(si.environment())
	print(Action_wind.Vb[1])
	qpz= Action_wind.Qp_z
	print(qpz)
	ffr= Action_wind.Ffr(15, "Lisse")
	print(ffr)
	print(Action_wind.rayon_secteur_angu)
	#Action_wind.show_Ffr()
	vertical = Toiture_isolee_2pants._from_parent_class(Action_wind, phi=0.5, load_area=1.9)
	#vertical.show_zonage()
	print(vertical.get_Cp("0°"))