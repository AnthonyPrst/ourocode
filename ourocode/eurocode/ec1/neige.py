#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
from math import sin, radians
import pandas as pd
import warnings

import forallpeople as si
from handcalcs.decorator import handcalc

from ourocode.eurocode.core.batiment import Batiment

class Neige(Batiment):
	EXPOSITION = ("Normal", "Protégé")
	TYPE_TOIT = ("1 versant", "2 versants", "Versants multiples")
	SNOW_LOAD_CARAC = {"A1": (0.45, 0, "A"),
						   "A2": (0.45, 1, "B1"),
						   "B1": (0.55, 1, "B1"),
						   "B2": (0.55, 1.35, "B1"),
						   "C1": (0.65, 0, "A"),
						   "C2": (0.65, 1.35, "B1"),
						   "D": (0.90, 1.8, "B1"),
						   "E": (1.40, 0, "A")}
	def __init__(self, exposition: str=EXPOSITION, type_toit :str=TYPE_TOIT, blocage: bool=("False", "True"), **kwargs):
		"""Initialise le calcul de l'action variable de neige selon l'EN 1991-1-3 et son Annexe Nationale française.

		Hérite de Batiment. La zone de neige et l'altitude sont issues des attributs du projet
		(code_INSEE → région de neige, alt → correction altimétrique).

		Args:
			exposition (str, optional): Type d'exposition du bâtiment selon EN 1991-1-3 AN §5.2.7 :
				- "Normal" : terrain sans protection particulière (C_e = 1.0).
				- "Protégé" : protégé par bâtiments alentour, déplacement de neige inexistant (C_e = 1.25).
				Defaults to "Normal".
			type_toit (str, optional): Forme du toit pour le calcul des coefficients de forme :
				"1 versant", "2 versants" ou "Versants multiples". Defaults to "2 versants".
			blocage (bool, optional): True si un dispositif bloque le glissement de la neige sur la toiture
				(mu1 = 0.8 quel que soit l'angle). Defaults to False.
			**kwargs: Arguments transmis à la classe parent Batiment.
		"""
		super().__init__(**kwargs)
		self.exposition = exposition
		self.type_toit = type_toit
		self.blocage = blocage
		self.C_t = 1
	
	@property
	def region_neige(self):
		file = "carte_action_region.csv"
		df = self._data_from_csv(file, index_col=1)
		return df.loc[str(self.code_INSEE)]["Zone_neige"]
	

	@property
	def S_k_200(self):
		"""Retourne la valeur caractéristique de la neige au sol à 200 m d'altitude S_k,200 en kN/m².

		Valeur issue du tableau des zones de neige de l'AN français (EN 1991-1-3 NA)
		en fonction de la région de neige déterminée par le code_INSEE du projet.

		Returns:
			forallpeople.Physical: S_k,200 en kN/m².
		"""
		return self.SNOW_LOAD_CARAC[self.region_neige][0] * si.kN/si.m**2
	

	@property
	def S_Ad(self):
		"""Retourne la valeur de neige au sol pour la situation accidentelle S_Ad en kN/m².

		Valeur issue du tableau des zones de neige de l'AN français. Utilisée pour le calcul
		de la charge de neige accidentelle S_x (EN 1991-1-3 NA §4.3).

		Returns:
			forallpeople.Physical: S_Ad en kN/m² (0 si la zone ne prévoit pas de situation accidentelle).
		"""
		return self.SNOW_LOAD_CARAC[self.region_neige][1] * si.kN/si.m**2
	
		
	@property
	def C_e(self):
		"""Retourne le coefficient d'exposition C_e selon EN 1991-1-3 AN §5.2.7.

		- "Normal" : C_e = 1.0.
		- "Protégé" : C_e = 1.25 (bâtiment protégé par obstacles environnants).

		Returns:
			float: Coefficient d'exposition (1.0 ou 1.25).
		"""
		match self.exposition:
			case "Normal":
				return 1
			case _:
				return 1.25

				
	@property
	def mu(self):
		"""Retourne les coefficients de forme de neige μ selon l'EN 1991-1-3 §5.3.

		Ne prend pas en compte les constructions proches ou plus élevées, ni les toitures cylindriques.
		Les valeurs dépendent de l'angle du versant alpha_toit (et alpha_toit2 pour 2 versants).

		Returns:
			dict: Selon le type de toit :
				- "1 versant" : {"mu1": valeur ou tuple (latex, valeur)}.
				- "2 versants" : {"mu2 versant 1": ..., "mu2 versant 2": ...}.
				- "Versants multiples" : {"mu2 versant 1": ..., "mu2 versant 2": ..., "mu3": ...}.
		"""
		
		match self.type_toit:
			case "1 versant":
				alpha_toit = self.alpha_toit
				if self.blocage:
					dict_mu = {"mu1": 0.8}
				
				if 0 <= alpha_toit <= 30:
					dict_mu = {"mu1": 0.8}
				elif 30 < alpha_toit < 60:
					@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
					def val():
						mu_1 = 0.8 * (60 - alpha_toit) / 30
						return mu_1
					dict_mu = {"mu1": val()}
				elif alpha_toit >= 60:
					dict_mu = {"mu1": 0}

					
			case "2 versants" | "Versants multiples":
				alpha_versant_1 = self.alpha_toit
				alpha_versant_2 = self.alpha_toit2
				dict_mu = {"mu2 versant 1": None, "mu2 versant 2": None}
				
				def calculate_mu2(alpha, jupyter_display):
					if 0 <= alpha <= 30:
						return 0.8
					elif 30 < alpha < 60:
						@handcalc(override="short", precision=2, jupyter_display=jupyter_display, left="\\[", right="\\]")
						def val():
							mu_2 = 0.8 * (60 - alpha) / 30
							return mu_2
						return val()[1]
					elif alpha >= 60:
						return 0

				for i, alpha in enumerate((alpha_versant_1, alpha_versant_2)):
					value = calculate_mu2(alpha, self.JUPYTER_DISPLAY)
					if i:
						dict_mu["mu2 versant 2"] = value
					else:
						dict_mu["mu2 versant 1"] = value


				if self.type_toit == "Versants multiples":
					alpha_mean = (alpha_versant_1 + alpha_versant_2) / 2
					if 0 <= alpha_mean <= 30:
						@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
						def val():
							alpha_m = (alpha_versant_1 + alpha_versant_2) / 2
							mu_3 = 0.8 + 0.8 * alpha_m / 30
							return  mu_3
						dict_mu["mu3"] = val()
					elif 30 < alpha_mean < 60:
						dict_mu["mu3"] = 1.6
					elif alpha_mean >= 60:
						dict_mu["mu3"] = None
		return dict_mu
		

	@property   
	def Sk_alt(self):
		"""Retourne la charge de neige caractéristique au sol S_k à l'altitude du projet en kN/m².

		Appliquer la correction altimétrique selon l'AN français de l'EN 1991-1-3 §4.1(1) :
		≤ 200 m : S_k = S_k,200 ; > 200 m : incrément croissant selon la zone (A, B1, E).

		Returns:
			tuple: (latex_string, S_k_alt) où S_k_alt est la charge caractéristique au sol en kN/m²
				(avec unité si.kN/si.m²).
		"""
		
		altitude = int(self.alt.value)
		S_k_200 = self.S_k_200
		
		if altitude <= 200:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
			def val():
				S_k_alt = S_k_200 + 0
				return S_k_alt
		
		elif 200 < altitude <= 500:
			if self.region_neige == "E":
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					delta_S_2 = 0.15 * (altitude - 200) /100
					S_k_alt = S_k_200 + delta_S_2 * si.kPa
					return S_k_alt
			else:
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					delta_S_1 = 0.10 * (altitude - 200) /100
					S_k_alt = S_k_200 + delta_S_1 * si.kPa
					return S_k_alt
				
		elif 500 < altitude <= 1000:
			if self.region_neige == "E":
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					delta_S_2 = 0.45 + 0.35 * (altitude - 500) /100
					S_k_alt = S_k_200 + delta_S_2 * si.kPa
					return S_k_alt
			else:
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					delta_S_1 = 0.30 + 0.15  * (altitude - 500) /100
					S_k_alt = S_k_200 + delta_S_1 * si.kPa
					return S_k_alt
				
		elif 1000 < altitude <= 2000:
			if self.region_neige == "E":
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					delta_S_2 = 2.20 + 0.70  * (altitude - 1000) /100
					S_k_alt = S_k_200 + delta_S_2 * si.kPa
					return S_k_alt
			else:
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					delta_S_1 = 1.05 + 0.35  * (altitude - 1000) /100
					S_k_alt = S_k_200 + delta_S_1 * si.kPa
					return S_k_alt
		return val()
	

	@property
	def Sn(self):
		"""Retourne les charges de neige normale sur toiture S_n en kN/m² selon EN 1991-1-3 §5.3.

		Formule : S_n = C_e × C_t × S_k,alt × μ.
		Majoration de 0.2 kN/m² appliquée automatiquement pour les faibles pentes (alpha < 1.718°).

		Returns:
			tuple: (latex_string, valeur) où valeur dépend du type de toit :
				- "1 versant" : float (charge unique en kN/m²).
				- "2 versants" : dict {"Versant 1": {"cas 1": ..., "cas 2": ...}, "Versant 2": ...}.
				- "Versants multiples" : dict {"cas 1": {"Versant 1": ..., "Versant 2": ...}, "cas 2": ...}.
		"""
		C_e = self.C_e
		C_t = self.C_t
		S_k_alt = self.Sk_alt[1]

		def get_mu_2():
			mu_2_versant_1 = self.mu["mu2 versant 1"]
			mu_2_versant_2 = self.mu["mu2 versant 2"]
			for i, mu in enumerate((mu_2_versant_1, mu_2_versant_2)):
				if isinstance(mu, tuple):
					if i:
						mu_2_versant_2 = mu[1]
					else:
						mu_2_versant_1 = mu[1]
			return mu_2_versant_1, mu_2_versant_2

		match self.type_toit:
			case "1 versant":
				mu_1 = self.mu["mu1"]
				if isinstance(mu_1, tuple):
					mu_1 = mu_1[1]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					S_n = C_e * C_t * S_k_alt * mu_1
					return S_n
				value = val()
				
			case "2 versants":
				mu_2_versant_1, mu_2_versant_2 = get_mu_2()
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					v1_cas_1 = C_e * C_t * S_k_alt * mu_2_versant_1
					v1_cas_2 = C_e * C_t * S_k_alt * mu_2_versant_1 * 0.5
					v2_cas_1 = C_e * C_t * S_k_alt * mu_2_versant_2
					v2_cas_2 = C_e * C_t * S_k_alt * mu_2_versant_2 * 0.5
					return v1_cas_1, v1_cas_2, v2_cas_1, v2_cas_2
				
				latex_sn, value = val()
				v1_cas_1, v1_cas_2, v2_cas_1, v2_cas_2 = value
				value = (latex_sn, {"Versant 1": {"cas 1": v1_cas_1, "cas 2": v1_cas_2},
			 			"Versant 2": {"cas 1": v2_cas_1, "cas 2": v2_cas_2}})
				
			case "Versants multiples":
				mu_2_versant_1, mu_2_versant_2 = get_mu_2()
				mu_3 = self.mu["mu3"]
				if isinstance(mu_3, tuple):
					mu_3 = mu_3[1]

				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					v1_cas_1 = C_e * C_t * S_k_alt * mu_2_versant_1
					v2_cas_1 = C_e * C_t * S_k_alt * mu_2_versant_2
					cas_2 = C_e * C_t * S_k_alt * mu_3
					return v1_cas_1, v2_cas_1, cas_2
				latex_sn, value = val()
				v1_cas_1, v2_cas_1, cas_2 = value
				value = (latex_sn, {"cas 1": {"Versant 1": v1_cas_1, "Versant 2": v2_cas_1,}, "cas 2": cas_2})

		# Majoration faible pente de 0.2 kN/m²
		if self.alpha_toit < 1.718:
			match self.type_toit:
				case "1 versant":
					S_nb = value[1]
					@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
					def val():
						S_n = S_nb + 0.2 * si.kN/si.m**2 # Majoration faible pente
						return S_n
					value2 = val()
					value = (value[0]+value2[0], value2[1])
		return value

	
	@property
	def Sx(self):
		"""Retourne les charges de neige accidentelle sur toiture S_x en kN/m² selon EN 1991-1-3 AN §4.3.

		Formule : S_x = C_e × C_t × S_Ad × μ.
		Utilisé uniquement pour les zones prévoyant une situation accidentelle (S_Ad > 0).

		Returns:
			tuple: (latex_string, valeur) où valeur dépend du type de toit :
				- "1 versant" : float (charge en kN/m²).
				- "2 versants" / "Versants multiples" : dict {"Versant 1 cas 1": ..., "Versant 2 cas 1": ...}.
		"""
		C_e = self.C_e
		C_t = self.C_t
		S_Ad = self.S_Ad

		match self.type_toit:
			case "1 versant":
				mu_1 = self.mu["mu1"]
				if isinstance(mu_1, tuple):
					mu_1 = mu_1[1]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					S_x = C_e * C_t * S_Ad * mu_1
					return S_x
				
			case "2 versants" | "Versants multiples":
				mu_2_versant_1 = self.mu["mu2 versant 1"]
				mu_2_versant_2 = self.mu["mu2 versant 2"]
				for i, mu in enumerate((mu_2_versant_1, mu_2_versant_2)):
					if isinstance(mu, tuple):
						if i:
							mu_2_versant_2 = mu[1]
						else:
							mu_2_versant_1 = mu[1]
				@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
				def val():
					v1_cas_1 = C_e * C_t * S_Ad * mu_2_versant_1
					v2_cas_1 = C_e * C_t * S_Ad * mu_2_versant_2
					return {"Versant 1 cas 1": v1_cas_1, "Versant 2 cas 1": v2_cas_1}
		return val()
	

	@property
	def Se(self):
		"""Retourne la charge linéaire de neige en débord de toiture S_e en kN/m.

		Applicable uniquement pour altitude > 900 m selon EN 1991-1-3 AN §5.7.
		En dessous de 900 m, retourne 0.

		Returns:
			float | dict: Charge en débord de toiture en kN/m :
				- "1 versant" : float.
				- "2 versants" : dict {"Versant 1": ..., "Versant 2": ...}.
		"""
		def calcul_Se(S_n:float):
			"""
			Args:
				S_n (float): charge de neige normal
			"""
			if self.alt.value > 900:
				gam = 3 * si.kN/si.m**3 # poid de reférence de la neige en kN/m3
				d = S_n / gam
				k = 3 / d.value
				
				if (k * si.kN / si.m**2) > (d * gam):
					k = d * gam
				Se = k * S_n**2 / gam
			else:
				Se = 0
			return Se

		match self.type_toit:
			case "1 versant":
				S_n = self.Sn[1]
				return calcul_Se(S_n)
			case "2 versants":
				dict_Se = {"Versant 1": None, "Versant 2": None}
				S_n_v1 = self.Sn[1]["Versant 1"]["cas 1"]
				S_n_v2 = self.Sn[1]["Versant 2"]["cas 1"]

				for i, S_n in enumerate((S_n_v1, S_n_v2)):
					Se = calcul_Se(S_n)
					if i:
						dict_Se["Versant 2"] = Se
					else:
						dict_Se["Versant 1"] = Se
				return dict_Se
	
	def mu2_construction_attenante(self, b1: si.m, b2: si.m, h: si.m, alpha_toit_attenant: float):
		"""Calcule le coefficient de forme μ2 pour une construction attenante selon EN 1991-1-3 §5.3.4.

		Détermine mu_W (neige apportée par le vent), mu_S (neige glissante depuis le toit attenant),
		mu2 = mu_W + mu_S, et la longueur d'application ls.

		Args:
			b1 (float): Longueur du bâtiment attenant (dans la direction perpendiculaire à la façade) en m.
			b2 (float): Longueur du bâtiment considéré (dans la même direction) en m.
			h (float): Différence de hauteur entre les toitures des deux bâtiments en m.
			alpha_toit_attenant (float): Angle de la toiture de la construction attenante en °.
				Si > 15°, mu_S ne peut être calculé directement (avertissement émis).

		Returns:
			dict: {"mu_W": ..., "mu_S": ..., "mu2": ..., "ls": ...} avec ls la longueur d'application en m.
		"""
		b_1 = b1
		b_2 = b2
		gam = 2 * si.kN/si.m**2
		S_k_alt = self.Sk_alt[1]
		muW = (b_1 + b_2) / (2 * h)
		if muW <= gam * h / S_k_alt:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
			def muW():
				mu_W = (b_1 + b_2) / (2 * h)
				return mu_W
		else:
			@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
			def muW():
				mu_W = gam * h / S_k_alt
				return mu_W
		dict_mu = {"mu_W": muW()[1]}
		if alpha_toit_attenant <= 15:
			dict_mu["mu_S"] = 0
		else:
			warnings.warn("Attention mu_S ne peut pas être calculé car l'angle de toiture de la construction attenante est supérieur à 15°.\n \
			Le coeff. est donc calculé via une charge aditionnelle égale à la moitié de la charge maximale totale sur le versant adjacent de la toiture supérieure calculée selon §5.3.3")
			dict_mu["mu_S"] = None
		dict_mu["mu2"] = 0 + dict_mu["mu_W"]
		if 5 * si.m <= 2*h <= 15 * si.m:
			dict_mu["ls"] = 2*h
		elif 2*h > 15 * si.m:
			dict_mu["ls"] = 15 * si.m
		else:
			dict_mu["ls"] = 5 * si.m
		return dict_mu
	
	def fs(self, Sn: si.kN/si.m**2, entraxe:si.m, alpha: float):
		"""Retourne la charge de neige en kN/m sur un obstacle ou une barre à neige

		Args:
			Sn (float): la charge de neige en kN/m²
			entraxe (int): entraxe horizontal entre deux obstacle à la neige en m
			alpha (float): angle du versant considéré.

		Returns:
			float: La charge en kN/m
		"""
		S_n = Sn * si.kN/si.m**2
		entraxe = entraxe * si.m
		@handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
		def val():
			F_s = S_n * entraxe * sin(radians(alpha))
			return F_s
		return val()



if __name__ == "__main__":
	from A0_Projet import Projet
	projet = Projet(code_INSEE=73215, alt=400)
	bat = Batiment._from_parent_class(projet, h_bat=5,b_bat=10,d_bat=10, alpha_toit=28, alpha_toit2=12)
	Action_snow = Neige._from_parent_class(bat, exposition="Normal", type_toit="Versants multiples", blocage=True)
	print(Action_snow.Sn)