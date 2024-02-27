#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT

import os
import sys

from math import sin, radians
import pandas as pd

import forallpeople as si
from handcalcs.decorator import handcalc

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Batiment



class Neige(Batiment):
	EXPOSITION = ("Normal", "Protégé")
	TYPE_TOIT = ("1 versant", "2 versants")
	SNOW_LOAD_CARAC = {"A1": (0.45, 0, "A"),
						   "A2": (0.45, 1, "B1"),
						   "B1": (0.55, 1, "B1"),
						   "B2": (0.55, 1.35, "B1"),
						   "C1": (0.65, 0, "A"),
						   "C2": (0.65, 1.35, "B1"),
						   "D": (0.90, 1.8, "B1"),
						   "E": (1.40, 0, "A")}
	def __init__(self, exposition: str=EXPOSITION, type_toit :str=TYPE_TOIT, blocage: bool=("False", "True"), **kwargs):
		"""Créer une classe permettant le calcul de l'action variable de neige sur les structure à l'Eurocode 1
		ATTENTION Uniquement pour le versant considéré ! Si plusieur versant alors plusieur classe

		Args:
			exposition (str, optional): Défini le type d'exposition, normal ou protégé par des batiments
			alentour rendant inexistant le déplacement possible de la neige par le vent (EC1_AN 5.2.7). Defaults to "normal".
			type_toit (str, optional): Défini le type de toit, "1 versant", "2 versants". Defaults to "2 versants"
			blocage (bool, optional): Défini si il y a un blocage du glissement de la neige sur la toiture. Defaults to False.
		"""
		super().__init__(**kwargs)
		self.exposition = exposition
		self.type_toit = type_toit
		self.blocage = blocage
		self.C_t = 1

	
	@property
	def _zone(self):
		file = "carte_action_region.csv"
		df = self._data_from_csv(file, index_col=1)
		return df.loc[str(self.code_INSEE)]["Zone_neige"]
	

	@property
	def S_k_200(self):
		"""Retourne la valeur de neige au sol pour une situation normale Sk200 en kN/m2

		Returns:
			float : Sk200 en kN/m2
		"""
		return self.SNOW_LOAD_CARAC[self._zone][0] * si.kN/si.m**2
	

	@property
	def S_Ad(self):
		"""Retourne la valeur de neige au sol pour une situation accidentelle SAd en kN/m2

		Returns:
			float : SAd en kN/m2
		"""
		return self.SNOW_LOAD_CARAC[self._zone][1] * si.kN/si.m**2
	
		
	@property
	def C_e(self):
		"""Détermine le coefficient d'exposition
		"""
		match self.exposition:
			case "Normal":
				return 1
			case _:
				return 1.25

				
	@property
	def mu(self):
		"""Retourne le mu d'une toiture simple et double pan et a versants multiples selon l'EN 1991-1-3 §5.3
		Ne prend pas en compte les constructions proches ou plus élevées, ni les toitures cylindrique.

		Returns:
			float : mu1, mu2, mu3
		"""
		alpha_toit = self.alpha_toit
		if 0 <= self.alpha_toit <= 30:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
			def val():
				mu_3 = 0.8 + 0.8 * alpha_toit / 30
				return {"mu1": 0.8, "mu2": 0.8, "mu3": mu_3}
		elif 30 < self.alpha_toit < 60:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
			def val():
				mu_1 = 0.8 * (60 - alpha_toit) / 30
				mu_2 = 0.8 * (60 - alpha_toit) / 30
				return {"mu1": mu_1, "mu2": mu_2, "mu3": 1.6}
		elif self.alpha_toit >= 60:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
			def val():
				return {"mu1": 0, "mu2": 0, "mu3": None}
		value = val()
		if self.blocage:
			value[1]["mu1"] = 0.8
		return value
		

	@property   
	def Sk_alt(self):
		"""Retourne la charge de neige caractéristique horizontal en kN/m² à une altitude donnée

		Returns:
			float: la charge en kN/m²
		"""
		
		altitude = int(self.alt.value)
		S_k_200 = self.S_k_200
		
		if altitude <= 200:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
			def val():
				S_k_alt = S_k_200 + 0
				return S_k_alt
		
		elif 200 < altitude <= 500:
			if self._zone == "E":
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					delta_S_2 = 0.15 * (altitude - 200) /100
					S_k_alt = S_k_200 + delta_S_2 * si.kPa
					return S_k_alt
			else:
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					delta_S_1 = 0.10 * (altitude - 200) /100
					S_k_alt = S_k_200 + delta_S_1 * si.kPa
					return S_k_alt
				
		elif 500 < altitude <= 1000:
			if self._zone == "E":
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					delta_S_2 = 0.45 + 0.35 * (altitude - 500) /100
					S_k_alt = S_k_200 + delta_S_2 * si.kPa
					return S_k_alt
			else:
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					delta_S_1 = 0.30 + 0.15  * (altitude - 500) /100
					S_k_alt = S_k_200 + delta_S_1 * si.kPa
					return S_k_alt
				
		elif 1000 < altitude <= 2000:
			if self._zone == "E":
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					delta_S_2 = 2.20 + 0.70  * (altitude - 1000) /100
					S_k_alt = S_k_200 + delta_S_2 * si.kPa
					return S_k_alt
			else:
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					delta_S_1 = 1.05 + 0.35  * (altitude - 1000) /100
					S_k_alt = S_k_200 + delta_S_1 * si.kPa
					return S_k_alt
		return val()
	

	@property
	def Sn(self):
		"""Retourne la charge de neige normal horizontal en kN/m² fonction du cas de chargement (non distribué/distribué)

		Returns:
			dict: le dictionnaire des charges
		"""
		C_e = self.C_e
		C_t = self.C_t
		S_k_alt = self.Sk_alt[1]
		match self.type_toit:
			case "1 versant":
				mu_1 = self.mu[1]["mu1"]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					S_n = C_e * C_t * S_k_alt * mu_1
					return S_n
			case "2 versants":
				mu_2 = self.mu[1]["mu2"]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					cas_1 = C_e * C_t * S_k_alt * mu_2
					cas_2 = C_e * C_t * S_k_alt * mu_2 * 0.5
					return {"cas 1": cas_1, "cas 2": cas_2}	
		value = val()

		# Majoration faible pente de 0.2 kN/m²
		if self.alpha_toit < 1.718:
			match self.type_toit:
				case "1 versant":
					S_n = value[1]
					@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
					def val():
						S_n = S_n + 0.2 # Majoration faible pente
						return S_n
					value2 = val()
					value = (value[0]+value2[0], value2[1])

				case "2 versants":
					S_n1 = value[1]["cas 1"]
					S_n2 = value[1]["cas 2"]
					@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
					def val():
						S_n1 = S_n1 + 0.2 # Majoration faible pente
						S_n2 = S_n2 + 0.2
						return {"cas 1": S_n1, "cas 2": S_n2}
					value2 = val()
					value = (value[0]+value2[0], value2[1])	
		return value

	
	@property
	def Sx(self):
		"""Retourne la charge de neige accidentelle horizontal en kN/m² fonction du cas de chargement (non distribué/distribué)
		
		Returns:
			dict: le dictionnaire des charges
		"""
		C_e = self.C_e
		C_t = self.C_t
		S_Ad = self.S_Ad

		match self.type_toit:
			case "1 versant":
				mu_1 = self.mu[1]["mu1"]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					S_x = C_e * C_t * S_Ad * mu_1
					return S_x
				
			case "2 versants":
				mu_2 = self.mu[1]["mu2"]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
				def val():
					cas_1 = C_e * C_t * S_Ad * mu_2
					return {"cas 1": cas_1}	
		return val()
	

	@property
	def Se(self):
		"""Retourne la charge linéaire en kN/m à appliquer en débord de toiture si l'altitude est > à 900m

		Returns:
			float: la charge verticale en débord de toiture en kN/m
		"""
		match self.type_toit:
			case "1 versant":
				S_n = self.Sn[1]
			case "2 versants":
				S_n = self.Sn[1]["cas 1"]

		if self.alt.value > 900:
			gam = 3 # poid de reférence de la neige en kN/m3
			d = S_n / gam
			k = 3 / d
			
			if k > (d * gam):
				k = d * gam
			Se = k * S_n**2 / gam
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
		entraxe = entraxe * si.m
		alpha_toit = self.alpha_toit
		match self.type_toit:
			case "1 versant":
				S_n = self.Sn[1]
			case "2 versants":
				S_n = self.Sn[1]["cas 1"]
		
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\[", right="\]")
		def val():
			F_s = S_n * entraxe * sin(radians(alpha_toit))
			return F_s
		value = val()
		self.F_s = value[1]
		return value



if __name__ == "__main__":
	Action_snow = Neige("79  Deux-Sèvres", exposition="Normal", type_toit="1 versant", blocage=True, h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=45, alt=0)
	print(Action_snow.Sk_alt)