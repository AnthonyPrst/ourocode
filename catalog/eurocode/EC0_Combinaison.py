# coding in UTF-8 
# by Anthony PARISOT

import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Projet


class Chargement(Projet):
	ACTION = ("Permanente G", 
			"Exploitation Q",
			"Neige normale Sn",
			"Vent pression W+",
			"Vent dépression W-",
			"Neige accidentelle Sx",
			"Sismique Ae")
	
	def __init__(self, **kwargs):
		"""Génère une classe Chargement qui permet de définir le chargement sur la barre.
		Cette classe est hérité de la classe Projet du module A0_Projet.py
		"""
		super().__init__(**kwargs)
		self._list_init_loads = []
	
	@property
	def init_loads(self):
		"""Retourne la liste des charges définis initialement.
		"""
		return self._list_init_loads

	def create_load(self, name: str, load: int, pos: str, action: str=ACTION, direction: str=("Z", "X", "Z local")):
		"""Ajoute une charge dans la liste de chargement de la classe chargement

		Args:
			name (str): nom de la charge.
			load (int): effort en daN ou daN/m suivant le "type_load" sélectionné.
			pos (str): position de la charge sur la barre en mm. 
							- Pour une charge nodale juste la position, ex: "1000".
							- Pour une charge linéique, il faut donner l'intervalle ex: "500/1000".
			action (str): type d'action de l'effort.
			type_load (str): type de chargement nodale ou linéique.
			direction (str): sens de l'effort sur la barre.
		"""	
		type_load = "Linéique"
		if pos.find("/") == -1:
			pos = int(pos)
			type_load = "Nodale"
		load = (len(self._list_init_loads)+1, name, action, type_load, load, pos, direction)
		self._list_init_loads.append(load)
		return load

	
	def create_load_by_list(self, _list_init_loads: list):
		"""Ajoute les charges d'une liste pré-défini dans la liste de chargement

		Args:
			_list_init_loads (list): liste de charge.
		"""
		for load in _list_init_loads:
			self._list_init_loads.append(load)
		return self._list_init_loads

	def del_load(self, index_load: int):
		"""Supprime une charge de l'attribut _list_init_loads par son index

		Args:
			index_load (int): index de la charge à supprimer.
		"""
		return self._list_init_loads.pop(index_load-1)
		


class Combinaison(Chargement):
	
	CAT_TYPE = ("Aucune",
				"Cat A : habitation",
				"Cat B : bureaux", 
				"Cat C : lieu de réunion", 
				"Cat D : zones commerciales", 
				"Cat E : stockage", 
				"Cat F : véhicule <= 30kN",
				"Cat G : véhicule <= 160kN",
				"Cat H : toits"
				)
	DICO_COMBI_ACTION = {"Permanente G": "G", 
						"Exploitation Q": "Q",
						"Neige normale Sn": "Sn",
						"Vent pression W+": "W+",
						"Vent dépression W-": "W-",
						"Neige accidentelle Sx": "Sx",
						"Sismique Ae": "Ae"
						}
	COEF_G = (1, 1.35) # Ginf, Gsup
	COEF_Q= (1.5)  # Qsup

	def __init__(self, ELU_STR: bool=("True", "False"), ELU_STR_ACC: bool=("True", "False"), ELS_C: bool=("True", "False"), ELS_QP: bool=("True", "False"), cat: str=CAT_TYPE, kdef: float=None, **kwargs):
		"""Créer une classe Combinaison qui génère les combinaisons d'action suivant les actions données.
		Cette classe est hérité de la classe Chargement du module EC0_Combinaison.py

		Args:
			ELU_STR (bool): combiner les charges à l'ELU structure, si vrai -> True, sinon False.
			ELU_STR_ACC (bool): combiner les charges à l'ELU accidentel, si vrai -> True, sinon False.
			ELS_C (bool): combiner les charges à l'ELS carctéristique, si vrai -> True, sinon False.
			ELS_QP (bool): combiner les charges à l'ELS quasi permananent, si vrai -> True, sinon False.
			cat (str): catégorie d'exploitation de la zone considéré. Defaults to "Aucune".
			kdef (float): Coefficient permettant de prendre en compte le fluage du bois en fonction de sa classe de service.
				Si le matériaux est autre que du bois laisser vide.
		"""
		super(Chargement, self).__init__(**kwargs)
		self.elu_STR = ELU_STR
		self.elu_STR_ACC = ELU_STR_ACC
		self.els_C = ELS_C
		self.els_QP = ELS_QP
		self.cat = cat
		self.name_combination = []
		self._generate_combinaison()
		if kdef:
			self._els_fleche_bois(kdef)
		

	def _generate_combinaison(self):
		"""Génère les combinaisons de chargement si une liste de chargement à été défini précédement.
		"""
		self.combiActionVariable = [0]*7 
		for load in self._list_init_loads:
			if load[2] == "Permanente G":
				self.combiActionVariable[0] = "G"

			elif load[2] == "Exploitation Q":
				self.combiActionVariable[1] = "Q"

			elif load[2] == "Neige normale Sn":
				self.combiActionVariable[2] = "Sn"

			elif load[2] == "Vent pression W+":
				self.combiActionVariable[3] = "W+"

			elif load[2] == "Vent dépression W-":
				self.combiActionVariable[4] = "W-"

			elif load[2] == "Neige accidentelle Sx":
				self.combiActionVariable[5] = "Sx"
			else:
				self.combiActionVariable[6] = "Ae"

		if self.elu_STR:
			self._ELU_STR()
		if self.elu_STR_ACC:
			self._ELU_STR_ACC()
		if self.els_C:
			self._ELS_C()
		if self.els_QP:
			self._ELS_QP()
		
	@property
	def coef_psy(self):
		""" Retourne les caractéristiques psy sous forme de dictionnaire"""
		dict_psy = {self.cat: {}, "Vent": {}, "Température": {}}
		if self.alt.value > 1000:
			dict_psy["Neige > 1000m"] = {}
		else:
			dict_psy["Neige <= 1000m"] = {}

		data_csv_psy = self._data_from_csv("coeff_psy.csv")
		psy_columns = data_csv_psy.columns.to_list()

		for psy_i in psy_columns:
			for key in dict_psy.keys():
				dict_psy[key][psy_i] = data_csv_psy.loc[key].loc[psy_i]
		return dict_psy


	def _key_action_psy(self, action_variable):
		if action_variable == "Q":
			index = self.cat
		elif action_variable == "Sn":
			if self.alt.value > 1000:
				index = "Neige > 1000m"
			else:
				index = "Neige <= 1000m"

		elif action_variable == "W+" or action_variable == "W-" :
			index = "Vent"
		else:
			index = "Température"
		return index

	
	def _create_array_load(self, name, value, load_list, array):
		load = np.array([name, load_list[0], load_list[1], load_list[2], load_list[3], value, load_list[5], load_list[6]], dtype=object)
		array= np.append(array,[load],axis= 0)
		return array

	def _create_dataframe_load(self, array_load):
		df_load = pd.DataFrame(array_load, columns=['Combinaison', 'Index', 'Nom', 'Action', 'Type', 'Valeur', 'Position', 'Axe'])
		#print(df_load)
		return(df_load)
	

	def _create_list_combination(self, nameCombi: str):
		test = 1
		for i in self.name_combination:
			if i == nameCombi:
				test = 0
				break
		if test:
			self.name_combination.append(nameCombi)
	
		
			
	
	def _ELU_STR(self):
		""" Combinaison à l'ELU STR """
		array_load = np.empty((0, 8))

		for load in self._list_init_loads:
			if self.DICO_COMBI_ACTION[load[2]] not in self.combiActionVariable[5:7]:
				name = "ELU_STR " + self.DICO_COMBI_ACTION[load[2]]
				self._create_list_combination(name)
				value = load[4]
				array_load = self._create_array_load(name, value, load, array_load)
			
			for action in self.combiActionVariable:
				if not action in ("Sx", "Ae"):
					if load[2]== "Permanente G" and action == "G":
							name = "ELU_STR " + str(self.COEF_G[1])+"G"
							self._create_list_combination(name)
							value = self.COEF_G[1] * load[4]
							array_load = self._create_array_load(name, value, load, array_load)
					
					elif action != 0 and action != "G":
						if action != "W-":
							name = "ELU_STR 1.35G + " + str(self.COEF_Q) + action
							self._create_list_combination(name)
						
						if load[2] != "Permanente G" and self.DICO_COMBI_ACTION[load[2]] == action:
							if action == "W-":
								name = "ELU_STR 1G" + " + " + str(self.COEF_Q) + action
								self._create_list_combination(name)
							value = self.COEF_Q * load[4]
							array_load = self._create_array_load(name, value, load, array_load)

						elif load[2] == "Permanente G":
							if action == "W-":
								name = "ELU_STR 1G" + " + " + str(self.COEF_Q) + action
								self._create_list_combination(name)
								value = self.COEF_G[0] * load[4]
							else:
								value = self.COEF_G[1] * load[4]
							array_load = self._create_array_load(name, value, load, array_load)
						
						for index in range(1,4):
							if action != self.combiActionVariable[index] and action != "W-":
								if self.combiActionVariable[index] != 0:
									key_psy = self._key_action_psy(self.combiActionVariable[index])
									name = "ELU_STR 1.35G + " + str(self.COEF_Q) + action + " + " + str(round(self.coef_psy[key_psy]["psy0"] * self.COEF_Q,2)) + self.combiActionVariable[index]
									self._create_list_combination(name)
									
									if load[2] == "Permanente G":
										value = self.COEF_G[1] * load[4]
										array_load = self._create_array_load(name, value, load, array_load)

									elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
										value = float(self.coef_psy[key_psy]["psy0"]) * self.COEF_Q  * load[4]
										array_load = self._create_array_load(name, value, load, array_load)

									elif self.DICO_COMBI_ACTION[load[2]] == action:
										value = self.COEF_Q * load[4]
										array_load = self._create_array_load(name, value, load, array_load)
									
									for index2 in range(1,4):
										if action != self.combiActionVariable[index2]:
											if self.combiActionVariable[index] != self.combiActionVariable[index2]:
												if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
													if self.combiActionVariable[index2] != 0:
														key_psy2 = self._key_action_psy(self.combiActionVariable[index2])
														name = "ELU_STR 1.35G + " + str(self.COEF_Q) + action + " + " + str(round(self.coef_psy[key_psy]["psy0"] * self.COEF_Q,2)) + self.combiActionVariable[index] + " + " + str(round(self.coef_psy[key_psy2]["psy0"] * self.COEF_Q,2)) + self.combiActionVariable[index2]
														self._create_list_combination(name)
														
														if load[2] == "Permanente G":
															value = self.COEF_G[1] * load[4]
															array_load = self._create_array_load(name, value, load, array_load)
														
														elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index2]:
															value = float(self.coef_psy[key_psy2]["psy0"]) * self.COEF_Q  * load[4]
															array_load = self._create_array_load(name, value, load, array_load)

														elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
															value = float(self.coef_psy[key_psy]["psy0"]) * self.COEF_Q  * load[4]
															array_load = self._create_array_load(name, value, load, array_load)

														elif self.DICO_COMBI_ACTION[load[2]] == action:
															value = self.COEF_Q * load[4]
															array_load = self._create_array_load(name, value, load, array_load)
		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELUSTR = self._create_dataframe_load(array_load)
		
	
	def _return_combi_ELUSTR(self, combi):
		return self.df_load_ELUSTR.loc[self.df_load_ELUSTR['Combinaison']==combi]
	

	def _ELU_STR_ACC(self):
		""" Combinaison à l'ELU STR ACC"""
		array_load = np.empty((0, 8))
		acc_action = ("Sx", "Ae")

		for load in self._list_init_loads:
			if self.DICO_COMBI_ACTION[load[2]] not in self.combiActionVariable[1:5]:
				name = "ELU_STR_ACC " + self.DICO_COMBI_ACTION[load[2]]
				self._create_list_combination(name)
				value = load[4]
				array_load = self._create_array_load(name, value, load, array_load)
			
			for action in self.combiActionVariable:
				# if action not in ("W-",):
					if action not in (0, "G"):
						key_psy_action = self._key_action_psy(action)
						if action in acc_action:
							name = "ELU_STR_ACC G + " + action
						elif self.coef_psy[key_psy_action]["psy1"] != 0:
							name = "ELU_STR_ACC_FIRE G + " + str(self.coef_psy[key_psy_action]["psy1"]) + action
						
						if name:
							self._create_list_combination(name)
							if load[2] != "Permanente G" and self.DICO_COMBI_ACTION[load[2]] == action:
								if action in acc_action:
									array_load = self._create_array_load(name, load[4], load, array_load)
								else:
									value = float(self.coef_psy[key_psy_action]["psy1"]) * load[4]
									array_load = self._create_array_load(name, value, load, array_load)
									
							elif load[2] == "Permanente G":
								array_load = self._create_array_load(name, load[4], load, array_load)
						
						for index in range(1,4):
							if action != self.combiActionVariable[index]:
								if self.combiActionVariable[index] != 0 and action not in ("W-",):
									key_psy = self._key_action_psy(self.combiActionVariable[index])
									if self.coef_psy[key_psy]["psy2"] != 0:
										if action in acc_action:
											name_2 = "ELU_STR_ACC G + " + action + " + " + str(self.coef_psy[key_psy]["psy2"]) + self.combiActionVariable[index]
										else:
											name_2 = "ELU_STR_ACC_FIRE G + " + str(self.coef_psy[key_psy_action]["psy1"]) + action + " + " + str(self.coef_psy[key_psy]["psy2"]) + self.combiActionVariable[index]
										
										try:
											self._create_list_combination(name_2)
											if load[2] == "Permanente G":
												array_load = self._create_array_load(name_2, load[4], load, array_load)

											elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
												value = float(self.coef_psy[key_psy]["psy2"]) * load[4]
												array_load = self._create_array_load(name_2, value, load, array_load)

											elif self.DICO_COMBI_ACTION[load[2]] == action:
												if action in acc_action:
													array_load = self._create_array_load(name_2, load[4], load, array_load)
												else:
													value = float(self.coef_psy[key_psy_action]["psy1"]) * load[4]
													array_load = self._create_array_load(name_2, value, load, array_load)
										except:
											pass
										
										for index2 in range(1,4):
											if action != self.combiActionVariable[index2] and action not in ("W-",):
												if self.combiActionVariable[index] != self.combiActionVariable[index2]:
													if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
														if self.combiActionVariable[index2] != 0:
															key_psy2 = self._key_action_psy(self.combiActionVariable[index2])
															if self.coef_psy[key_psy2]["psy2"] != 0:
																if action in acc_action:
																	name_3 = "ELU_STR_ACC G + " + action + " + " + str(self.coef_psy[key_psy]["psy2"]) + self.combiActionVariable[index] + " + " + str(self.coef_psy[key_psy2]["psy2"]) + self.combiActionVariable[index2]
																else:
																	name_3 = "ELU_STR_ACC_FIRE G + " + str(self.coef_psy[key_psy_action]["psy1"]) + action + " + " + str(self.coef_psy[key_psy]["psy2"]) + self.combiActionVariable[index] + " + " + str(self.coef_psy[key_psy2]["psy2"]) + self.combiActionVariable[index2]
																
																try:
																	self._create_list_combination(name_3)

																	if load[2] == "Permanente G":
																		array_load = self._create_array_load(name_3, load[4], load, array_load)
																	
																	elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index2]:
																		value = float(self.coef_psy[key_psy2]["psy2"]) * load[4]
																		array_load = self._create_array_load(name_3, value, load, array_load)

																	elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
																		value = float(self.coef_psy[key_psy]["psy2"]) * load[4]
																		array_load = self._create_array_load(name_3, value, load, array_load)

																	elif self.DICO_COMBI_ACTION[load[2]] == action:
																		if action in acc_action:
																			array_load = self._create_array_load(name_3, load[4], load, array_load)
																		else:
																			value = float(self.coef_psy[key_psy_action]["psy1"]) * load[4]
																			array_load = self._create_array_load(name_3, value, load, array_load)
																except:
																	pass
		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELU_STR_ACC = self._create_dataframe_load(array_load)
		
	
	def _return_combi_ELU_STR_ACC(self, combi):
		return self.df_load_ELU_STR_ACC.loc[self.df_load_ELU_STR_ACC['Combinaison']==combi]
	

	def _ELS_C(self):
		array_load = np.empty((0, 8))

		for load in self._list_init_loads:
			for action in self.combiActionVariable:

				if load[2]== "Permanente G" and action == "G":
					name = "ELS_C G"
					self._create_list_combination(name)
					value = load[4]
					array_load = self._create_array_load(name, value, load, array_load)
				
				elif action != 0 and action != "G":
					name = "ELS_C G + " + action
					self._create_list_combination(name)

					if load[2] != "Permanente G" and self.DICO_COMBI_ACTION[load[2]] == action:
						value = load[4]
						array_load = self._create_array_load(name, value, load, array_load)

					elif load[2] == "Permanente G":
						value = load[4]
						array_load = self._create_array_load(name, value, load, array_load)
				
					for index in range(1,4):
							if action != self.combiActionVariable[index] and action != "W-":
								key_psy = self._key_action_psy(self.combiActionVariable[index])
								
								if self.combiActionVariable[index] != 0 and float(self.coef_psy[key_psy]["psy0"]) != 0:
									name = "ELS_C G + " + action + " + " + str(round(self.coef_psy[key_psy]["psy0"],2)) + self.combiActionVariable[index]
									self._create_list_combination(name)

									if load[2] == "Permanente G":
										value = load[4]
										array_load = self._create_array_load(name, value, load, array_load)

									elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
										value = float(self.coef_psy[key_psy]["psy0"]) * load[4]
										array_load = self._create_array_load(name, value, load, array_load)

									elif self.DICO_COMBI_ACTION[load[2]] == action:
										value = load[4]
										array_load = self._create_array_load(name, value, load, array_load)
									
									for index2 in range(1,4):
										if action != self.combiActionVariable[index2]:
											if self.combiActionVariable[index] != self.combiActionVariable[index2]:
												if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
													if self.combiActionVariable[index2] != 0:
														key_psy2 = self._key_action_psy(self.combiActionVariable[index2])
														name = "ELS_C G + " + action + " + " + str(round(self.coef_psy[key_psy]["psy0"],2)) + self.combiActionVariable[index] + " + " + str(round(self.coef_psy[key_psy2]["psy0"],2)) + self.combiActionVariable[index2]
														self._create_list_combination(name)

														if load[2] == "Permanente G":
															value = load[4]
															array_load = self._create_array_load(name, value, load, array_load)
														
														elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index2]:
															value = float(self.coef_psy[key_psy2]["psy0"]) * load[4]
															array_load = self._create_array_load(name, value, load, array_load)

														elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
															value = float(self.coef_psy[key_psy]["psy0"]) * load[4]
															array_load = self._create_array_load(name, value, load, array_load)

														elif self.DICO_COMBI_ACTION[load[2]] == action:
															value = load[4]
															array_load = self._create_array_load(name, value, load, array_load)

		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELScarac = self._create_dataframe_load(array_load)
	

	def _return_combi_ELScarac(self, combi):
		return self.df_load_ELScarac.loc[self.df_load_ELScarac['Combinaison']==combi]


	def _ELS_QP(self):
		array_load = np.empty((0, 8))

		for load in self._list_init_loads:
			for action in self.combiActionVariable:

				if load[2]== "Permanente G" and action == "G":
					name = "ELS_QP G"
					self._create_list_combination(name)
					value = load[4]
					array_load = self._create_array_load(name, value, load, array_load)
				
				elif action != 0 and action != "G":

					key_psy = self._key_action_psy(action)
					if float(self.coef_psy[key_psy]["psy2"]) != 0:
						name = "ELS_QP G + " +  str(round(self.coef_psy[key_psy]["psy2"],2)) + action
						self._create_list_combination(name)

						if load[2] != "Permanente G" and self.DICO_COMBI_ACTION[load[2]] == action:
							value = float(self.coef_psy[key_psy]["psy2"]) * load[4]
							array_load = self._create_array_load(name, value, load, array_load)

						elif load[2] == "Permanente G":
							value = load[4]
							array_load = self._create_array_load(name, value, load, array_load)
				
						for index in range(1,4):
							if action != self.combiActionVariable[index] and action != "W-":
								if self.combiActionVariable.index(action) < self.combiActionVariable.index(self.combiActionVariable[index]):
									if self.combiActionVariable[index] != 0:
										key_psy1 = self._key_action_psy(self.combiActionVariable[index])
										if float(self.coef_psy[key_psy1]["psy2"]) != 0:
											name = "ELS_QP G + " +  str(round(self.coef_psy[key_psy]["psy2"],2)) + action + " + " + str(round(self.coef_psy[key_psy1]["psy2"],2)) + self.combiActionVariable[index]
											self._create_list_combination(name)
											if load[2] == "Permanente G":
												value = load[4]
												array_load = self._create_array_load(name, value, load, array_load)

											elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
												value = float(self.coef_psy[key_psy1]["psy2"]) * load[4]
												array_load = self._create_array_load(name, value, load, array_load)

											elif self.DICO_COMBI_ACTION[load[2]] == action:
												value = float(self.coef_psy[key_psy]["psy2"]) * load[4]
												array_load = self._create_array_load(name, value, load, array_load)
											
											for index2 in range(1,4):
												if action != self.combiActionVariable[index2]:
													if self.combiActionVariable[index] != self.combiActionVariable[index2]:
														if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
															if self.combiActionVariable[index2] != 0:
																key_psy2 = self._key_action_psy(self.combiActionVariable[index2])
																if float(self.coef_psy[key_psy2]["psy2"]) != 0:
																	name = "ELS_QP G + " + action + " + " + str(round(self.coef_psy[key_psy1]["psy2"],2)) + self.combiActionVariable[index] + " + " + str(round(self.coef_psy[key_psy2]["psy2"],2)) + self.combiActionVariable[index2]
																	self._create_list_combination(name)

																	if load[2] == "Permanente G":
																		value = load[4]
																		array_load = self._create_array_load(name, value, load, array_load)
																	
																	elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index2]:
																		value = float(self.coef_psy[key_psy2]["psy2"]) * load[4]
																		array_load = self._create_array_load(name, value, load, array_load)

																	elif self.DICO_COMBI_ACTION[load[2]] == self.combiActionVariable[index]:
																		value = float(self.coef_psy[key_psy1]["psy2"]) * load[4]
																		array_load = self._create_array_load(name, value, load, array_load)

																	elif self.DICO_COMBI_ACTION[load[2]] == action:
																		value = float(self.coef_psy[key_psy]["psy2"]) * load[4]
																		array_load = self._create_array_load(name, value, load, array_load)

		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELSqp = self._create_dataframe_load(array_load)


	def _return_combi_ELSqp(self, combi):
		return self.df_load_ELSqp.loc[self.df_load_ELSqp['Combinaison']==combi]
	


	def _els_fleche_bois(self, kdef: float):
		"""Génère l'association des combinaisons ELS caractéristique + quasi permanente (fluage, intégrant le Kdef) pour le calcul d'un élément en bois.
		"""

		if self.els_C:
			# On détermine W_inst(Q), pour cela on enlève W_inst_G
			self.df_W_inst_Q = self.df_load_ELScarac[self.df_load_ELScarac["Action"] != "Permanente G"]
			for index in range(self.df_W_inst_Q.shape[0]):
				name = self.df_W_inst_Q.iloc[index,0]
				name_combi = "W_inst " + name[10:]
				self.df_W_inst_Q.iloc[index,0] = name_combi
				self._create_list_combination(name_combi)
			# print(self.df_W_inst_Q)

			if self.els_QP:
				combi_qp = {}
				combi_c = {}
				value_search = ["Q", "S"]

				# On multiplie les combinaisons quasi permanente par Kdef pour trouver Wcreep
				self.df_load_ELSqp["Valeur"] = self.df_load_ELSqp["Valeur"] * kdef

				for name in self.name_combination:
					if name[0:6] == "ELS_QP":
						typeQP = [0]*2
						for i in range(len(value_search)):
							try:
								if name[8:].index(value_search[i]):
									typeQP[i] = value_search[i]
							except ValueError:
								pass
						combi_qp[name] = typeQP

					elif name[0:5] == "ELS_C":
						typeC = [0]*2
						for i in range(len(value_search)):
							try:
								if name[7:].index(value_search[i]):
									if value_search[i] == "Q" and self.cat == "Cat H : toits":
										pass
									elif value_search[i] == "S" and self.alt.value <= 1000:
										pass
									else:
										typeC[i] = value_search[i]
							except ValueError:
								pass
						combi_c[name] = typeC

				list_of_key = list(combi_qp.keys())
				list_of_value = list(combi_qp.values())

				array_load = np.empty((0, 8))
				for name_combi_C, val in combi_c.items():
					position = list_of_value.index(val)
					name_combi_QP = list_of_key[position]
					name_combi = "W_net_fin "+ name_combi_C + " & " + name_combi_QP

					df_c = self.df_load_ELScarac[self.df_load_ELScarac["Combinaison"] == name_combi_C]
					df_qp = self.df_load_ELSqp[self.df_load_ELSqp["Combinaison"] == name_combi_QP]
					for index in range(df_c.shape[0]):
						index_c = df_c.iloc[index, 1]
						name_load = df_c.iloc[index, 2]
						action_load = df_c.iloc[index, 3]
						type_load = df_c.iloc[index, 4]
						valeur_c = df_c.iloc[index, 5]
						position = df_c.iloc[index, 6]
						axe = df_c.iloc[index, 7]
						if len(df_qp[df_qp["Index"] == index_c]):
							valeur_qp = df_qp[df_qp["Index"] == index_c].iloc[0,5]
							valeur = valeur_c + valeur_qp
						else:
							valeur = valeur_c
						load = np.array([name_combi, index_c, name_load, action_load, type_load, valeur, position, axe], dtype=object)
						array_load = np.append(array_load,[load],axis= 0)
						self._create_list_combination(name_combi)

				array_load = array_load[array_load[:, 0].argsort()]
				self.df_W_net_fin = self._create_dataframe_load(array_load)
				# print(self.df_W_net_fin)
 

	def _return_combi_W_inst_Q(self, combi):
		return self.df_W_inst_Q.loc[self.df_W_inst_Q['Combinaison']==combi]
	
	def _return_combi_W_net_fin(self, combi):
		return self.df_W_net_fin.loc[self.df_W_net_fin['Combinaison']==combi]

	@property
	def list_combination(self):
		"""Retourne un data frame avec toute les combinaison créer
		"""
		self.name_combination.sort()
		return pd.DataFrame(self.name_combination, columns=['Combinaison'])
	

	def get_list_combination(self):
		return self.name_combination
	
	
	def _choice_combi_df(self):
		shape = len(self.list_combination)
		dict_load_combi = {}

		for i in range(shape):
			combi = self.list_combination.iloc[i,0]

			if self.elu_STR_ACC and combi[0:11] == 'ELU_STR_ACC':
				df_combi = self._return_combi_ELU_STR_ACC(combi)

			elif self.elu_STR and combi[0:7] == 'ELU_STR':
				df_combi = self._return_combi_ELUSTR(combi)

			elif self.els_C and combi[0:5] == 'ELS_C':
				df_combi = self._return_combi_ELScarac(combi)

			elif self.els_QP and combi[0:6] == 'ELS_QP':
				df_combi = self._return_combi_ELSqp(combi)

			elif self.els_C and combi[0:6] == 'W_inst':
				df_combi = self._return_combi_W_inst_Q(combi)

			elif self.els_C and self.els_QP and combi[0:9] == 'W_net_fin':
				df_combi = self._return_combi_W_net_fin(combi)

			df_combi = df_combi.drop(labels=['Combinaison'], axis=1)
			load_list_combi = df_combi.values.tolist()
			dict_load_combi[combi] =  load_list_combi
		
		return dict_load_combi


	def get_combi_list_load(self, nom: str):
		"""Retourne la liste des charges combinées pour la combinaison sélectionné

		Args:
			nom (str): nom des la combinaison à récupérer. Defaults to "Sélectionner tout".
		"""
		self.list_loads = self._choice_combi_df()[nom]
		return self.list_loads

	
	def min_type_load(self, name_combi: str) -> str:
		"""Retourne le type de chargement de plus courte durée

		Args:
			name_combi (str): Combinaison à analyser
		"""
  
		dictName = {"G": "Permanente",
					"Q": "Moyen terme",
					"Sn": "Court terme",
					"W+": "Instantanee", 
					"W-": "Instantanee",
					"Sx": "Instantanee",
					"Ae": "Instantanee"
					}
		for action in self.combiActionVariable:
			if action:
				indexAction = name_combi[8:].find(action)
				if indexAction > -1 :
					if action == "Sn" and self.alt.value >= 1000 :
						name_load_type = dictName["Q"]
					elif action == "Sn" and self.cat == "Cat H : toits":
						name_load_type = dictName["Sn"]
					else:
						name_load_type = dictName[action]
		return name_load_type

	

if __name__== "__main__":

	_list_init_loads = [[1, 'G', 'Permanente G', 'Linéique', -10, '0/100', 'Z'],
				 [0, 'Ae', "Sismique Ae", 'Linéique', -50, '0/100', 'Z'],
				 [2, '', 'Neige normale Sn', 'Linéique', -200, '0/100', 'Z'],
				 [3, '', 'Exploitation Q', 'Linéique', -100, '0/100', 'Z'],
				 [4, '', "Vent dépression W-", 'Linéique', 300, '0/100', 'Z'],
				 [5, '', "Vent pression W+", 'Linéique', -300, '0/100', 'Z'],
				 ]
	projet = Projet("AP", "6018.0", "", "", 73215, "France", 1200)
	chargement = Chargement._from_parent_class(projet)
	chargement.create_load_by_list(_list_init_loads)
	c1 = Combinaison._from_parent_class(chargement, ELU_STR=False, ELU_STR_ACC=True, ELS_C=False, ELS_QP=False, cat="Cat A : habitation", kdef=0.6)
	rcombi = "ELS_QP G + 0.3Q"
	# print(c1._return_combi_ELUSTR(rcombi))
	print(pd.DataFrame(c1.coef_psy))
	# print(c1.df_load_ELScarac)
	# print(c1.df_load_ELSqp)
	# print(c1.df_W_inst_Q)
	# print(c1.df_load_ELU_STR_ACC)

	print(c1.list_combination)

	
	