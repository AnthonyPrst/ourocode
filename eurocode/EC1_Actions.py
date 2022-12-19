#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import os
import sys

import math
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from objet import Objet

def interpolation_lineaire(x, xa, xb, ya, yb):
	"""Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
	y = ya + (x - xa) * ((yb - ya)/(xb - xa))
	return y

def interpolation_logarithmique(x, ya, yb):
	"""Fait une interpolation linéaire pour trouver un résultat y entre deux valeur xa et xb """
	y = ya - (ya - yb) * math.log10(x)
	return y


class Snow(Objet):
	
	def __init__(self, region:str, alt:int, alphatoit:float, exposition="normal", type_toit="2 versants", blocage=False):
		"""Créer une classe permettant le calcul de l'action variable de neige sur les structure à l'Eurocode 1
		ATTENTION Uniquement pour le versant considéré ! Si plusieur versant alors plusieur classe

		Args:
			region (str): Région ou ce situe le projet suivant "Zone_neige.csv"
			alt (int): Altitude du site du projet en m
			exposition (str, optional): Défini le type d'exposition, normal ou protégé par des batiments
			alentour rendant inexistant le déplacement possible de la neige par le vent (EC1_AN 5.2.7). Defaults to "normal".
			alphatoit (int): Angle du versant en degrés
			type_toit (str, optional): Défini le type de toit, "1 versant", "2 versants". Defaults to "2 versants"
			blocage (bool, optional): Défini si il y a un blocage du glissement de la neige sur la toiture. Defaults to False.
		"""
		super().__init__()
		self.region = region
		self.alt = alt
		self.alphatoit = alphatoit
		self.exposition = exposition
		self.type_toit = type_toit
		self.blocage = blocage
		self.Ct = 1
	

	def __data_from_csv(self, data_file:str):
		"""Retourne un dataframe d'un fichier CSV

		Args:
			data_file (str): Nom du fichier CSV à importer

		Returns:
			data_csv: dataframe pandas du fichier CSV
		"""
		repertory = os.getcwd() + "\Data" + data_file
		data_csv = pd.read_csv(repertory, sep=';', index_col=0)
		return data_csv

	@property
	def _zone(self):
		file = "\zone_neige.csv"
		df = self.__data_from_csv(file)
		return df.loc[self.region][0]
		
	@property
	def _sk200_sAd(self):
		"""Retourne la valeur de neige au sol sk200 et sAd en kN/m2

		Returns:
			dict : dictionnaire avec sk200 et sAd
		"""
		snow_load_carac = {"A1": (0.45, 0, "A"),
						   "A2": (0.45, 1, "B1"),
						   "B1": (0.55, 1, "B1"),
						   "B2": (0.55, 1.35, "B1"),
						   "C1": (0.65, 0, "A"),
						   "C2": (0.65, 1.35, "B1"),
						   "D": (0.90, 1.8, "B1"),
						   "E": (1.40, 0, "A")}
		
		Sk200 = snow_load_carac[self._zone][0]
		SAd = snow_load_carac[self._zone][1]
		
		return {"Sk200": Sk200, "SAd": SAd}
		
	@property
	def _Ce(self):
		"""Détermine le coefficient d'exposition
		"""
		match self.exposition:
			case "normal":
				ce = 1
			case _:
				ce = 1.25
		return ce
				
	@property
	def _mu(self):
		"""Retourne le mu d'une toiture simple et double pan

		Returns:
			float : mu1
		"""
		if 0 <= self.alphatoit <= 30:
			mu1 = 0.8
		elif 30 < self.alphatoit < 60:
			mu1 = 0.8 * (60 - self.alphatoit) / 30
		else:
			mu1 = 0
			
		if self.blocage:
			mu1 = 0.8
		
		return mu1
		
	@property   
	def _sk_altitude(self):
		"""Retourne la charge de neige caractéristique horizontal en kN/m² à une altitude donnée

		Returns:
			float: la charge en kN/m²
		"""
		if self.alt <= 200:
			Sk_alt = self._sk200_sAd["Sk200"]
		elif 200 < self.alt <= 500:
			if self._zone == "E":
				deltaS2 = 0.15 * ((self.alt - 200) /100)
				Sk_alt = self._sk200_sAd["Sk200"] + deltaS2
			else:
				deltaS1 = 0.10 * ((self.alt - 200) /100)
				Sk_alt = self._sk200_sAd["Sk200"] + deltaS1
				
		elif 500 < self.alt <= 1000:
			if self._zone == "E":
				deltaS2 = 0.45 + 0.35 * ((self.alt - 500) /100)
				Sk_alt = self._sk200_sAd["Sk200"] + deltaS2
			else:
				deltaS1 = 0.30 + 0.15  * ((self.alt - 500) /100)
				Sk_alt = self._sk200_sAd["Sk200"] + deltaS1
				
		elif 1000 < self.alt <= 2000:
			if self._zone == "E":
				deltaS2 = 2.20+ 0.70  * ((self.alt - 1000) /100)
				Sk_alt = self._sk200_sAd["Sk200"] + deltaS2
			else:
				deltaS1 = 1.05 + 0.35  * ((self.alt - 1000) /100)
				Sk_alt = self._sk200_sAd["Sk200"] + deltaS1
				
		return Sk_alt

	@property
	def Sn(self):
		"""Retourne la charge de neige normal horizontal en kN/m² fonction du cas de chargement (non distribué/distribué)

		Returns:
			dict: le dictionnaire des charges
		"""
		match self.type_toit:
			case "1 versant":
				cas1 = self._Ce * self.Ct * self._sk_altitude * self._mu
				cas2 = None
			case "2 versants":
				cas1 = self._Ce * self.Ct * self._sk_altitude * self._mu
				cas2 = self._Ce * self.Ct * self._sk_altitude * (self._mu*0.5)
		# Majoration faible pente de 0.2 kN/m²
		if self.alphatoit < 1.718:
			cas1 = cas1 + 0.2
			cas2 = cas2 + 0.2
		Sn = {"cas 1": cas1, "cas 2": cas2}
		return Sn

	
	@property
	def Sx(self):
		"""Retourne la charge de neige accidentelle horizontal en kN/m² fonction du cas de chargement (non distribué/distribué)
		
		Returns:
			dict: le dictionnaire des charges
		"""
		match self.type_toit:
			case "1 versant":
				cas1 = self._Ce * self.Ct * self._sk200_sAd["SAd"] * self._mu
				cas2 = None
			case "2 versants":
				cas1 = self._Ce * self.Ct * self._sk200_sAd["SAd"] * self._mu
				cas2 = self._Ce * self.Ct * self._sk200_sAd["SAd"] * (self._mu*0.5)
			
		Sx = {"cas 1": cas1, "cas 2": cas2}
		return Sx

	@property
	def Se(self):
		"""Retourne la charge linéaire en kN/m à appliquer en débord de toiture si l'altitude est > à 900m

		Returns:
			float: la charge verticale en débord de toiture en kN/m
		"""
		if self.alt > 900:
			gam = 3 # poid de reférence de la neige en kN/m3
			d = self.Sn["cas 1"] / gam
			k = 3 / d
			
			if k > (d * gam):
				k = d * gam
			Se = k * self.Sn["cas 1"]**2 / gam
		else:
			Se = 0
			
		return Se
	
	
	def fs(self, entraxe:int):
		"""Retourne la charge de neige en kN/m sur un obstacle ou une barre à neige

		Args:
			entraxe (int): entraxe horizontal entre deux obstacle à la neige en m

		Returns:
			float: La charge en kN/m
		"""
		self.Fs = self.Sn["cas 1"] * entraxe * math.sin(math.radians(self.alphatoit))
		return self.Fs


#################################### VENT #################################################
class Wind(object):
	RHO_AIR = 1.225 #kg/m3
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

	def __init__(self, region:str, cat_terrain:str, cat_oro:str, z:float, alt:int):
		self.region = region
		self.cat_terrain = self.CAT_TERRAIN[cat_terrain]
		self.cat_oro = self.CAT_ORO[cat_oro]
		self.z = z
		self.alt = alt
		self.CsCd = 1

	def _data_from_csv(self, data_file:str):
		"""Retourne un dataframe d'un fichier CSV

		Args:
			data_file (str): Nom du fichier CSV à importer

		Returns:
			data_csv: dataframe pandas du fichier CSV
		"""
		repertory = os.getcwd() + "\Data" + data_file
		data_csv = pd.read_csv(repertory, sep=';', index_col=0)
		return data_csv


	@property
	def _zone(self):
		file = "\zone_vent.csv"
		df = self._data_from_csv(file)
		return df.loc[self.region][0]


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
			Cr_z = self._Kr * math.log(self.z/self.cat_terrain["Z0"])
		else:
			Cr_z = self._Kr * math.log(self.cat_terrain["Zmin"]/self.cat_terrain["Z0"])

		print("Kr :", self._Kr)
		print("Cr_z :", Cr_z)
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
					C0_z = 1 + 0.004 * delta_AC*math.exp(z-10)
				else:
					C0_z = 1 + 0.004 * delta_AC*math.exp(0)

			case "2":
				pass
			case _:
				C0_z = 1
		print("C0_z :", C0_z)
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
				return self._Co_z * (1 - 2*10**-4 *(math.log10(self.cat_terrain["Z0"])+3)**6) #4.19NA
			case _:
				return 1 - 2*10**-4 *(math.log10(self.cat_terrain["Z0"])+3)**6 #4.20NA
		   

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
		"""Retourne le coefficient d'exposition à la hauteur z (informatif)

		Returns:
			float: coef d'expo
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

	def K_red_U(self, h, d):
		"""Retourne le coefficient de défaut de corrélation entre les pressions aérodynamiques au vent et sous le vent,
		uniquement pour les faces D et E. §7.2.2(3)

		Args:
			h (float): hauteur du bâtiment en m
			d (float): longeur parallèle à la direction du vent en m
		"""
		if h/d >= 5:
			self.Kred_U = 1
		elif h/d <= 1:
			self.Kred_U = 0.85
		else:
			self.Kred_U = 0.85 + 0.15 * (((h/d)-1)/4)
			
			
class Building(Wind):
	def __init__(self, region: str, cat_terrain: str, cat_oro: str, z: float, alt: int, h: float, d: float, b: float, alphatoit:float):
		super().__init__(region, cat_terrain, cat_oro, z, alt)
		self.h = h
		self.d = d
		self.b = b # coté perpendiculaire au vent longpant
		self.alphatoit = alphatoit
  
	def height_ref(self):
		pass # 7.2.1 bâtiment de grande hauteur
		
	
class Gabled_building(Building):
	def __init__(self, region: str, cat_terrain: str, cat_oro: str, z: float, alt: int, h: float, d: float, b: float, alphatoit:float):
		super().__init__(region, cat_terrain, cat_oro, z, alt, h, d, b, alphatoit)
  
		self.wind_direction = {"0°":"", "90°":""}
		
		for cle in self.wind_direction.keys():
			match cle:
				case "0°":
					self.e = min(self.b, self.h*2)
				case "90°":
					self.e = min(self.d, self.h*2)
					self.b = d
					self.d = b
	 
			self.vertical_wall(cle)
			
		
	def vertical_wall(self, cle):
		#print(self.e, self.d, cle)
		df = self._data_from_csv("/vent_CPE_mur_verticaux.csv")
		print(df)
		if self.e < self.d:  
			geometrie = {"A": self.e/5,
						 "B": 0.8 * self.e,
						 "C": self.d - self.e}

		elif self.e >= 5*self.d:
			geometrie = {"A": self.d}
		else:
			geometrie = {"A": self.e/5,
						 "B": self.d - self.e/5}
			
		self.wind_direction[cle] = {"geometrie": geometrie}

		

		
			
					
	
	
if __name__ == "__main__":
	
	Action_snow = Snow("79  Deux-Sèvres", 1000, 30, blocage=True)

	Action_wind = Wind("73  Savoie", "IIIa", "Aucun", 2, 250)
	print(Action_wind._Qp_z)
	building = Gabled_building("73  Savoie", "IIIa", "Aucun", 2, 250, 5, 8, 10, 30)
	print(building.wind_direction)
