#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import os
import sys

import math
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Batiment


class Neige(Batiment):
	REGION_NEIGE = list(Batiment._data_from_csv(Batiment, "zone_neige.csv").index)
	EXPOSITION = ("Normal", "Protégé")
	TYPE_TOIT = ("1 versant", "2 versants")
	def __init__(self, region_neige: str=REGION_NEIGE, exposition: str=EXPOSITION, type_toit :str=TYPE_TOIT, blocage: bool=("False", "True"), **kwargs):
		"""Créer une classe permettant le calcul de l'action variable de neige sur les structure à l'Eurocode 1
		ATTENTION Uniquement pour le versant considéré ! Si plusieur versant alors plusieur classe

		Args:
			region_neige (str): Région ou ce situe le projet suivant "Zone_neige.csv"
			exposition (str, optional): Défini le type d'exposition, normal ou protégé par des batiments
			alentour rendant inexistant le déplacement possible de la neige par le vent (EC1_AN 5.2.7). Defaults to "normal".
			type_toit (str, optional): Défini le type de toit, "1 versant", "2 versants". Defaults to "2 versants"
			blocage (bool, optional): Défini si il y a un blocage du glissement de la neige sur la toiture. Defaults to False.
		"""
		super().__init__(**kwargs)
		self.region_neige = region_neige
		self.exposition = exposition
		self.type_toit = type_toit
		self.blocage = blocage
		self.Ct = 1
	
	@property
	def _zone(self):
		file = "zone_neige.csv"
		self.df = self._data_from_csv(file)
		return self.df.loc[self.region_neige][0]
		
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
			case "Normal":
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
	
	
	def fs(self, entraxe:float):
		"""Retourne la charge de neige en kN/m sur un obstacle ou une barre à neige

		Args:
			entraxe (int): entraxe horizontal entre deux obstacle à la neige en m

		Returns:
			float: La charge en kN/m
		"""
		self.Fs = self.Sn["cas 1"] * entraxe * math.sin(math.radians(self.alphatoit))
		return self.Fs



if __name__ == "__main__":
	Action_snow = Neige("79  Deux-Sèvres", exposition="Normal", type_toit="1 versant", blocage=True, h_bat=5, d_bat=15, b_bat=13.1, alphatoit=15)
	print(Action_snow.Sn)