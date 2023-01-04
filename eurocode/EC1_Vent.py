#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import os
import sys
from PIL import Image

import math as mt
from copy import copy

import pandas as pd
sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Building


def interpolation_lineaire(x, xa, xb, ya, yb):
	"""Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
	y = ya + (x - xa) * ((yb - ya)/(xb - xa))
	return y

def interpolation_logarithmique(x, xa, xb, ya, yb):
	"""Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
	y = (x > xb) * yb + (x < xa) * ya + (x <= xb) * (x >= 1) * (ya - (ya - yb) * mt.log10(x))
	return y


class Wind(Building):
	RHO_AIR = 1.225 #kg/m3
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


	def __init__(self, z:float, region_vent:str, terrain:str, oro: str="Aucun", CsCd = 1, **kwargs):
		"""Créer une classe qui permet de calculer l'action de vent

		Args:
			z (float): hauteur en m sur le bâtiment ou est étudié le vent (Ze suivant EN 1991-1-4 §7.2.2).
			alt (int): altitude du batiment étudié en m.
			region_vent (str): Région ou ce situe le projet.
			cat_terrain (str): Catégorie de terrain du projet.
			cat_oro (str): Catégorie orographique. Default to "Aucun".
		"""
		super().__init__(**kwargs)
		self.region_vent = region_vent
		self.terrain = terrain
		self.oro = oro
		self.z = z
		self.CsCd = CsCd


	def _data_from_csv(self, data_file:str):
		"""Retourne un dataframe d_bat'un fichier CSV

		Args:
			data_file (str): Nom du fichier CSV à importer

		Returns:
			data_csv: dataframe pandas du fichier CSV
		"""
		repertory = os.path.join(os.getcwd(), "data", data_file)
		data_csv = pd.read_csv(repertory, sep=';', index_col=0)
		return data_csv


	@property
	def cat_terrain(self):
		return __class__.CAT_TERRAIN[self.terrain]

	@property
	def cat_oro(self):
		return __class__.CAT_ORO[self.oro]

	@property
	def _zone(self):
		file = os.path.join("vent", "zone_vent.csv")
		self.df = self._data_from_csv(file)
		return self.df.loc[self.region_vent][0]


	@property
	def _Vb0(self):
		"""Vb,0 est la valeur de base de la vitesse de référence du vent en m/s

		Returns:
			int: vitesse de référence du vent en m/s
		"""
		return self.VB0[str(self._zone)]


	@property
	def _Vb(self, Cdir=1, Cseason=1):
		"""Retourne la vitesse de référence du vent en m/s fonction de cette dernière et de la période de l'année
		à une hauteur de 10 m au-dessus d'un sol relevant de la catégorie de terrain II 

		Args:
			Cdir (float): coefficient de direction. Defaults to 1.
			Cseason (float): coefficient de saison. Defaults to 1.

		Returns:
			float: vitesse de référence du vent en m/s
		"""
		return self._Vb0* Cdir * Cseason


	@property
	def _rayon_secteur_angu(self):
		"""Retourne le rayon du secteur angulaire dans lequel la rugosité de terrain est à qualifier fonction de la hauteur du bâtiment étudier

		Returns:
			float: max entre 300m et R
		"""
		return max(23 * self.z**1.2, 300)


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
		
		if self.cat_terrain["Zmin"] <= self.z <= 200:
			Cr_z = self._Kr * mt.log(self.z/self.cat_terrain["Z0"])
		else:
			Cr_z = self._Kr * mt.log(self.cat_terrain["Zmin"]/self.cat_terrain["Z0"])

		# print("Kr :", self._Kr)
		# print("Cr_z :", Cr_z)
		return Cr_z


	@property
	def _Co_z(self):
		"""Retourne le coeficient orographique du terrain suivant la catégorie orographique et la hauteur z du bâtiment

		Returns:
			float: coef orographique
		"""
		# ATTENTION vérifier le z et zmin si toujours comptatible quand ajout du cas 2 !!!
		if self.z < self.cat_terrain["Zmin"]:
			z = self.cat_terrain["Zmin"]
		else:
			z = self.z

		match self.cat_oro:
			case "1":
				dict_alti = {"An1": "", "An2": "", "Ae1": "", "Ae2": "",
							 "As1": "", "As2": "", "Ao1": "", "Ao2": ""}

				for cle, value in dict_alti.items():
					dict_alti[cle] = int(input("Altitude {0} en m: ".format(cle)))

				Am = (2*self.alt + dict_alti["An1"] + dict_alti["An2"] + dict_alti["Ae1"] + dict_alti["Ae2"] + 
					 dict_alti["As1"] + dict_alti["As2"] + dict_alti["Ao1"] + dict_alti["Ao2"])/10

				delta_AC = self.alt - Am

				if self.z >= 10:
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
	def _Vm_z(self):
		"""Retourne la vitesse moyenne du vent à une hauteur z au dessus du sol en fonction de la rugosité de terrain,
		   de l'orographie et de la vitesse de référence du vent Vb

		Returns:
			float: vitesse moyenne du vent en m/s
		"""
		return self._Cr_z * self._Co_z * self._Vb
	

	@property
	def _Kl(self):
		"""Retourne le coefficient de turbulence

		Returns:
			float: coef de turbulence
		"""
		match self.cat_oro:
			case "1":
				return self._Co_z * (1 - 2*10**-4 *(mt.log10(self.cat_terrain["Z0"])+3)**6) #4.19NA
			case _:
				return 1 - 2*10**-4 *(mt.log10(self.cat_terrain["Z0"])+3)**6 #4.20NA
		   

	@property
	def _sigma_v(self):
		"""Retourne l'écart type de turbulence

		Returns:
			float: écart type
		"""
		return self._Kr * self._Vb * self._Kl #4.6


	@property
	def _Iv_z(self):
		"""Retourne l'intensité de turbulence à la hauteur z

		Returns:
			float: intensité de la turbulence
		"""

		if self.cat_terrain["Zmin"] <= self.z <= 200 or self.z < self.cat_terrain["Zmin"] :
			return self._sigma_v / self._Vm_z #4.7
		else:
			print("""Erreur la hauteur max du bâtiment ne peut dépasser 200m, 
					 impossible de calculer l'intensité de turbulence Iv""")


	@property
	def _Qb(self):
		"""Retour la pression dynamique de référence (informatif)

		Returns:
			float: pression en N/mm²
		"""
		return 0.5 * self.RHO_AIR * self._Vb**2 #4.10


	@property
	def _Qp_z(self):
		"""Retour la pression dynamique de pointe à la hauteur z qui est induite par la vitesse moyenne 
		et les fluctuation rapides de vitesse

		Returns:
			float: pression en N/mm²
		"""
		return (1 + 7 * self._Iv_z) * 0.5 * self.RHO_AIR * self._Vm_z**2 #4.8


	@property
	def _Ce_z(self):
		"""Retourne le coefficient d_bat'exposition à la hauteur z (informatif)

		Returns:
			float: coef d_bat'expo
		"""
		return self._Qp_z/self._Qb #4.9
	

	def We(self, Cpe):
		"""Retourne la pression aérodynamique sur les surfaces extérieures

		Args:
			Cpe (_type_): coefficient de pression extérieure

		Returns:
			float: pression en N/mm²
		"""
		return self._Qp_z*Cpe #5.1
	

	def Wi(self, Cpi):
		"""Retourne la pression aérodynamique sur les surfaces intérieures

		Args:
			Cpi (_type_): coefficient de pression intérieure

		Returns:
			float: pression en N/mm²"""
		return self._Qp_z*Cpi #5.2

	def Fw_e(self, Aref):
		self.CsCd
		pass #5.5

	def Fw_i(self, Aref):
		pass #5.6
	
	def Ffr(self, Cfr, Afr):
		pass #5.7

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
			
		

class Murs_verticaux(Wind):
	def __init__(self, load_area:int|float, *args, **kwargs):
		"""Créer une classe permettant le calcul des voiles verticaux au vent selon l'EN 1991-1-4 §7.2.2

		Args:
			load_area (int | float): aire chargée pour le calcul des éléments ou des fixations
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
		self.df = self._data_from_csv(os.path.join("vent", "vent_Cpe_mur_verticaux.csv"))
		self.df.reset_index(drop=False, inplace=True)

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
		if h_d > 0.25 and h_d < 1:
			df_index = [3,4]
			borne_inf = 0.25
			borne_sup = 1
			
		elif h_d > 1 and h_d < 5:
			df_index = [2,3]
			borne_inf = 1
			borne_sup = 5

		elif h_d <= 0.25:
			self.df = self.df[["Zone", "Cpe", "0.25"]]

		elif h_d >= 5:
			self.df = self.df[["Zone", "Cpe", "5"]]

		else:
			self.df = self.df[["Zone", "Cpe", "1"]]


		if borne_inf:
			for i in range(self.df.shape[0]):
				list_CPE_line.append(interpolation_lineaire(h_d, borne_inf, borne_sup, self.df.iloc[i,df_index[1]], self.df.iloc[i,df_index[0]]))
			self.df[round(h_d, 3)] = list_CPE_line
			self.df = self.df[["Zone", "Cpe", round(h_d, 3)]]

		if self.load_area > 1 and self.load_area < 10:
				for i in range(self.df.shape[0]):
					if i % 2 == 0:
						self.df.loc[self.df.shape[0]]= [self.df.iloc[i,0], "CPE "+ str(round(self.load_area, 2)), interpolation_logarithmique(self.load_area, 1, 10, self.df.iloc[i+1,2], self.df.iloc[i,2])]

		self.df.set_index("Zone", inplace=True)
		self.df = self.df.loc[[zone for zone in self.wind_direction[cle]["geometrie"].keys()]]
		return self.df


	def get_wind_dict(self) -> dict:
		return self.wind_direction


	def get_Cpe(self, dir="0°"):
		return self.wind_direction[dir]["Cpe"]


	def show_zonage(self):
		"""Affiche l'image du zonage pour les murs verticaux
		"""
		file = os.path.join(os.getcwd(), "data", "vent", "vent_Cpe_mur_verticaux.png")
		image = Image.open(file)
		image.show()




class Toiture_2pants(Wind):
	"""Créer une classe permetant le calcul d'une toiture à deux versant symétrique au vent selon l'EN 1991-1-4 §7.2.5
	"""
	def __init__(self, load_area:int|float, *args, **kwargs):
		super().__init__(*args, **kwargs)
		


class Toiture_isolee_1pant(Wind):
	def __init__(self, phi: int|float, load_area:int|float, *args, **kwargs):
		"""Créer une classe permetant le calcul d'une toiture isolée à un versant au vent selon l'EN 1991-1-4 §7.3

		Args:
			phi (int | float): le degré d_bat'obstruction sous une toiture isolée entre 0 et 1.
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
		self.df = self._data_from_csv(os.path.join("vent", "vent_Cp_toiture_isolee_1_versant.csv"))
		self.df.reset_index(drop=False, inplace=True)
		 
		geometrie = {"A": {"lenght": [(self.d_bat - self.d_bat/10*2), (self.b_bat - self.b_bat/10*2)], "surface": (self.d_bat - self.d_bat/10*2) * (self.b_bat - self.b_bat/10*2) / mt.cos(mt.radians(self.alphatoit))},
					"B": {"lenght": [(self.d_bat - self.d_bat/10*2), self.b_bat/10], "surface": ((self.d_bat - self.d_bat/10*2) * self.b_bat/10) / mt.cos(mt.radians(self.alphatoit))},
					"C": {"lenght": [(self.b_bat - self.b_bat/10*2), self.d_bat/10], "surface": ((self.b_bat - self.b_bat/10*2) * self.d_bat/10) / mt.cos(mt.radians(self.alphatoit))}}
		return geometrie


	def _Cp(self):
		"""Calcul les Cp d'une toiture isolée

		Args:
			cle (str): sens du vent sur le bâtiment 0° ou 90°
		"""
		
		if self.phi > 0 and self.phi < 1:
			for i in range(0, self.df.shape[0], 3):
				self.df.loc[self.df.shape[0]]= [self.df.iloc[i,0],
												self.phi,
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,2], self.df.iloc[i+2,2]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,3], self.df.iloc[i+2,3]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,4], self.df.iloc[i+2,4]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,5], self.df.iloc[i+2,5])
												]
			self.df = self.df[self.df["phi"].isin([self.phi, "max"])]

		elif self.phi >= 1:
			self.df = self.df[self.df["phi"].isin(["1", "max"])]
		else:
			self.df = self.df[self.df["phi"].isin(["0", "max"])]
			
			
		
		list_alphatoit= self.df["alpha_toit"].unique()

		if not self.alphatoit in list_alphatoit:
			minimum = list_alphatoit[list_alphatoit < self.alphatoit].max()
			maximum = list_alphatoit[list_alphatoit > self.alphatoit].min()
			df_min = self.df[self.df["alpha_toit"]==minimum]
			df_max = self.df[self.df["alpha_toit"]==maximum]

			self.df.reset_index(drop=True, inplace=True)
			for i, phi in enumerate(["max", self.phi]):
				self.df.loc[self.df.shape[0]]= [round(self.alphatoit, 2),
												phi,
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,2], df_max.iloc[i,2]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,3], df_max.iloc[i,3]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,4], df_max.iloc[i,4]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,5], df_max.iloc[i,5])
												]
			self.df = self.df[self.df["phi"].isin([self.phi, "max"])]

		self.df = self.df[self.df["alpha_toit"]==round(self.alphatoit, 2)]
		self.df.set_index("alpha_toit", inplace=True)
		return self.df


	def get_wind_dict(self) -> dict:
		return self.wind_direction
	

	def get_Cp(self, dir="0°"):
		return self.wind_direction[dir]["Cp"]


	def show_zonage(self):
		"""Affiche l'image du zonage pour une toiture isolée un versant
		"""
		file = os.path.join(os.getcwd(), "data", "vent", "vent_Cp_toiture_isolee_1_versant.png")
		image = Image.open(file)
		image.show()



class Toiture_isolee_2pants(Wind):
	def __init__(self, phi: int|float, load_area:int|float, *args, **kwargs):
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
		self.df = self._data_from_csv(os.path.join("vent", "vent_Cp_toiture_isolee_2_versants.csv"))
		self.df.reset_index(drop=False, inplace=True)
		 
		geometrie = {"A": {"lenght": [(self.d_bat - self.d_bat/10*2),  self.d_bat - (2*self.d_bat/10) - self.d_bat/5], 
							"surface": ((self.d_bat - self.d_bat/10*2) * (self.d_bat - (2*self.d_bat/10) - self.d_bat/5)) / mt.cos(mt.radians(self.alphatoit))},
					"B": {"lenght": [(self.d_bat - self.d_bat/10*2),  self.b_bat/10], "surface": ((self.d_bat - self.d_bat/10*2) * self.b_bat/10) / mt.cos(mt.radians(self.alphatoit))},
					"C": {"lenght": [(self.b_bat - self.b_bat/10*2), self.d_bat/10], "surface": ((self.b_bat - self.b_bat/10*2) * self.d_bat/10) / mt.cos(mt.radians(self.alphatoit))},
					"D": {"lenght": [(self.b_bat - self.b_bat/10*2)/2, (self.d_bat/5)/2], "surface": (((self.b_bat - self.b_bat/10*2) * self.d_bat/5) / mt.cos(mt.radians(self.alphatoit))/2)}}
		return geometrie


	def _Cp(self):
		"""Calcul les Cp d'une toiture isolée

		Args:
			cle (str): sens du vent sur le bâtiment 0° ou 90°
		"""
		
		if self.phi > 0 and self.phi < 1:
			for i in range(0, self.df.shape[0], 3):
				self.df.loc[self.df.shape[0]]= [self.df.iloc[i,0],
												self.phi,
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,2], self.df.iloc[i+2,2]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,3], self.df.iloc[i+2,3]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,4], self.df.iloc[i+2,4]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,5], self.df.iloc[i+2,5]),
												interpolation_lineaire(self.phi, 0, 1, self.df.iloc[i+1,6], self.df.iloc[i+2,6])
												]
			self.df = self.df[self.df["phi"].isin([self.phi, "max"])]

		elif self.phi >= 1:
			self.df = self.df[self.df["phi"].isin(["1", "max"])]
		else:
			self.df = self.df[self.df["phi"].isin(["0", "max"])]
			
			
		
		list_alphatoit= self.df["alpha_toit"].unique()

		if not self.alphatoit in list_alphatoit:
			minimum = list_alphatoit[list_alphatoit < self.alphatoit].max()
			maximum = list_alphatoit[list_alphatoit > self.alphatoit].min()
			df_min = self.df[self.df["alpha_toit"]==minimum]
			df_max = self.df[self.df["alpha_toit"]==maximum]

			self.df.reset_index(drop=True, inplace=True)
			for i, phi in enumerate(["max", self.phi]):
				self.df.loc[self.df.shape[0]]= [round(self.alphatoit, 2),
												phi,
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,2], df_max.iloc[i,2]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,3], df_max.iloc[i,3]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,4], df_max.iloc[i,4]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,5], df_max.iloc[i,5]),
												interpolation_lineaire(self.alphatoit, minimum, maximum, df_min.iloc[i,6], df_max.iloc[i,6])
												]
			self.df = self.df[self.df["phi"].isin([self.phi, "max"])]

		self.df = self.df[self.df["alpha_toit"]==round(self.alphatoit, 2)]
		self.df.set_index("alpha_toit", inplace=True)
		return self.df


	def get_wind_dict(self) -> dict:
		return self.wind_direction
	

	def get_Cp(self, dir="0°"):
		return self.wind_direction[dir]["Cp"]


	def show_zonage(self):
		"""Affiche l'image du zonage pour une toiture isolée un versant
		"""
		file = os.path.join(os.getcwd(), "data", "vent", "vent_Cp_toiture_isolee_2_versants.png")
		image = Image.open(file)
		image.show()



if __name__ == "__main__":
	
	building = Building(h_bat=5, d_bat=15, b_bat=13.1, alphatoit=15, alt=400)
	Action_wind = Wind._from_parent_class(building, region_vent="88  Vosges", terrain="IIIa", oro="Aucun", z=5)
	vertical = Toiture_isolee_2pants._from_parent_class(Action_wind, phi=0.5, load_area=1.9)
	vertical.show_zonage()
	print(vertical.get_Cp())