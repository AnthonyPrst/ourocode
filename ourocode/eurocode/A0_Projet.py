#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os, sys
import forallpeople as si
from handcalcs.decorator import handcalc

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# from objet import Objet
from ourocode.eurocode.objet import Objet


class Projet(Objet):
	JUPYTER_DISPLAY = False
	def __init__(self, ingenieur: str=None, num_project: str=None, name: str=None, adresse: str=None, 
				code_INSEE: int=None, pays: str="France", alt: si.m=0,**kwargs):
		"""Créer une classe Projet hérité de la classe Objet du fichier objet.py. Cette classe défini le projet, d'ou découle l'ensemble des objets du catalogue.

		Args:
			ingenieur (str, optional): nom de l'ingénieur Defaults to None.
			num_project (str, optional): numéro du projet. Defaults to None.
			name (str, optional): nom du projet. Defaults to None.
			adresse (str, optional): adresse du projet. Defaults to None.
			region (int, optional): numéro INSEE départementale du projet en 5 chiffres. Defaults to None.
			pays (str, optional): pays ou ce situe le projet. Defaults to "France".
			alt (int, optional): altitude du projet en m. Defaults to 0.
		"""
		super().__init__()
		for key, val in kwargs.items():
			setattr(self, key, val)
		self.ingenieur = ingenieur
		self.num_project = num_project
		self.name = name
		self.adresse = adresse
		self.code_INSEE = code_INSEE
		self.pays = pays
		self.alt = alt * si.m
		
	def __str__(self) -> str:
		return "Créer une classe Projet hérité de la classe Objet du fichier objet.py. Cette classe défini le projet, d'ou découle l'ensemble des objets du catalogue."


	def __repr__(self) -> str:
		return super().__repr__()


class Batiment(Projet):
	def __init__(self, h_bat: float, d_bat: float, b_bat: float, alpha_toit: float, alpha_toit2: float=0, *args, **kwargs):
		"""Créer une classe Batiment héritée de Wind, cette classe défini les dimension du bâtiment

		Args:
			h_bat (float): hauteur du bâtiment en m.
			d_bat (float): largeur du bâtiment en m.
			b_bat (float): longueur du bâtiment en m.
			alpha_toit (float): angle de toiture en ° du versant 1.
			alpha_toit2 (float): angle de toiture en ° du versant 2 si il existe sinon 0.
		"""
		super().__init__(*args, **kwargs)
		self.h_bat = h_bat
		self.d_bat = d_bat
		self.b_bat = b_bat # coté perpendiculaire au vent longpant
		self.alpha_toit = alpha_toit
		self.alpha_toit2 = alpha_toit2 
  
	def height_ref(self):
		pass # 7.2.1 bâtiment de grande hauteur


if __name__ == "__main__":
	projet = Projet(num_project="6006.0", commmentaire="c'est mon premier projet")
	building = Batiment(h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=15)
	print(building.h_bat)
	