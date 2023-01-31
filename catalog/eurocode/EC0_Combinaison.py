# coding in UTF-8 
# by Anthony PARISOT

import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Projet

class Combinaison(Projet):
	"""  """

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
	def __init__(self, list_load:list, cat: str=CAT_TYPE, **kwargs):
		"""Créer une classe Combinaison qui génère les combinaisons d'action suivant les actions données.
		Cette classe est hérité de la classe Projet du module A0_Projet.py

		Args:
			list_load (_type_, optional): liste des charges provenant du module MEF_. Defaults to list.
			cat (str, optional): catégorie d'exploitation de la zone considéré. Defaults to "Categorie A".
		"""
		super().__init__(**kwargs)
		self.list_load = list_load
		self.cat = cat
		self.coefg = (1, 1.35) # Ginf, Gsup
		self.coefq = (1.5)  # Qsup
		self.name_combination = []
		self.dicoCombiAct = {"Permanente G": "G", 
							 "Exploitation Q": "Q",
							 "Neige normale Sn": "Sn",
							 "Vent pression W+": "W+",
							 "Vent dépression W-": "W-",
							 "Neige accidentelle Sx": "Sx",
							 "Sismique Ae": "Ae"
							 }

		self.combiActionVariable = [0]*7 
		for nload in range(len(self.list_load)):
			if self.list_load[nload][2] == "Permanente G":
				self.combiActionVariable[0] = "G"

			elif self.list_load[nload][2] == "Exploitation Q":
				self.combiActionVariable[1] = "Q"

			elif self.list_load[nload][2] == "Neige normale Sn":
				self.combiActionVariable[2] = "Sn"

			elif self.list_load[nload][2] == "Vent pression W+":
				self.combiActionVariable[3] = "W+"

			elif self.list_load[nload][2] == "Vent dépression W+":
				self.combiActionVariable[4] = "W-"

			elif self.list_load[nload][2] == "Neige accidentelle Sx":
				self.combiActionVariable[5] = "Sx"
			else:
				self.combiActionVariable[6] = "Ae"


		self._elu_STR()
		self._els_C()
		self._els_QP()


	@property
	def coef_psy(self):
		""" Retourne les caractéristiques psy sous forme de liste python """
		list_psy = []
		listval = []
		data_csv_psy = self._data_from_csv("coeff_psy.csv")
		if self.cat != 'Aucune':
			for i in range(3):
				listval.append(data_csv_psy.loc[self.cat].loc["psy"+str(i)])
			list_psy.append(listval)
			listval = []
		else: 
			list_psy.append([None, None, None])
		
		for i in range(3):
			if self.alt>1000: 
				listval.append(data_csv_psy.loc["Neige > 1000m"].loc["psy"+str(i)])
			else:
				listval.append(data_csv_psy.loc["Neige < 1000m"].loc["psy"+str(i)])
		list_psy.append(listval)
		listval = []
		for i in range(3):      
			listval.append(data_csv_psy.loc["Vent"].loc["psy"+str(i)])
		list_psy.append(listval)
		
		return list_psy


	def _index_action_psy(self, action_variable):
		if action_variable == "Q":
			index = 0
		elif action_variable == "Sn":
			index = 1
		elif action_variable == "W+" or action_variable == "W-" :
			index = 2
		else:
			index = 0
		return index

	
	def _create_array_load(self, name, value, i, array):
		load = np.array([name, self.list_load[i][0], self.list_load[i][1], self.list_load[i][2], self.list_load[i][3], value, self.list_load[i][5], self.list_load[i][6]], dtype=object)
		array= np.append(array,[load],axis= 0)
		return array

	def _create_dataframe_load(self, array_load):
		df_load = pd.DataFrame(array_load, columns=['Combinaison', 'Index', 'Nom', 'Action', 'Type', 'Valeur', 'Position', 'Axe'])
		#print(df_load)
		return(df_load)
	

	def _create_list_combination(self, nameCombi):
		test = 1
		for i in self.name_combination:
			if i == nameCombi:
				test = 0
				break
		if test:
			self.name_combination.append(nameCombi)
			
	
	def _elu_STR(self):
		""" Combinaison à l'ELU STR """
		lenListLoad = len(self.list_load)
		array_load = np.empty((0, 8))

		for nload in range(lenListLoad):
			name = "ELU_STR " + self.dicoCombiAct[self.list_load[nload][2]]
			self._create_list_combination(name)
			value = self.list_load[nload][4]
			array_load = self._create_array_load(name, value, nload, array_load)
			
			for nAct in range(5):
			
				if self.list_load[nload][2]== "Permanente G" and self.combiActionVariable[nAct] == "G":
					# if self.list_load[nload][2]== "Vent dépression W-" and self.combiActionVariable[nAct] == "W-":
					#     name = "ELU_STR " + str(self.coefg[0])+"G" + str(self.coefq) + self.combiActionVariable[nAct]
					#     self._create_list_combination(name)
					#     value = self.coefg[0] * self.list_load[nload][4]
					#     array_load = self._create_array_load(name, value, nload, array_load)
					# else :
						name = "ELU_STR " + str(self.coefg[1])+"G"
						self._create_list_combination(name)
						value = self.coefg[1] * self.list_load[nload][4]
						array_load = self._create_array_load(name, value, nload, array_load)
				
				elif self.combiActionVariable[nAct] != 0 and self.combiActionVariable[nAct] != "G":
					if self.combiActionVariable[nAct] != "W-":
						name = "ELU_STR 1.35G + " + str(self.coefq) + self.combiActionVariable[nAct]
						self._create_list_combination(name)
					
					if self.list_load[nload][2] != "Permanente G" and self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
						if self.combiActionVariable[nAct] == "W-":
							name = "ELU_STR 1G" + " + " + str(self.coefq) + self.combiActionVariable[nAct]
							self._create_list_combination(name)
						
						value = self.coefq * self.list_load[nload][4]
						array_load = self._create_array_load(name, value, nload, array_load)

					elif self.list_load[nload][2] == "Permanente G":
						if self.combiActionVariable[nAct] == "W-":
							name = "ELU_STR 1G" + " + " + str(self.coefq) + self.combiActionVariable[nAct]
							self._create_list_combination(name)
							value = self.coefg[0] * self.list_load[nload][4]
						else:
							value = self.coefg[1] * self.list_load[nload][4]
						array_load = self._create_array_load(name, value, nload, array_load)
					
					for index in range(1,4):
						if self.combiActionVariable[nAct] != self.combiActionVariable[index] and self.combiActionVariable[nAct] != "W-":
							if self.combiActionVariable[index] != 0:
								typePsy = self._index_action_psy(self.combiActionVariable[index])
								name = "ELU_STR 1.35G + " + str(self.coefq) + self.combiActionVariable[nAct] + " + " + str(round(self.coef_psy[typePsy][0] * self.coefq,2)) + self.combiActionVariable[index]
								self._create_list_combination(name)
								
								if self.list_load[nload][2] == "Permanente G":
									value = self.coefg[1] * self.list_load[nload][4]
									array_load = self._create_array_load(name, value, nload, array_load)

								elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
									value = float(self.coef_psy[typePsy][0]) * self.coefq  * self.list_load[nload][4]
									array_load = self._create_array_load(name, value, nload, array_load)

								elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
									value = self.coefq * self.list_load[nload][4]
									array_load = self._create_array_load(name, value, nload, array_load)
								
								for index2 in range(1,4):
									if self.combiActionVariable[nAct] != self.combiActionVariable[index2]:
										if self.combiActionVariable[index] != self.combiActionVariable[index2]:
											if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
												if self.combiActionVariable[index2] != 0:
													typePsy2 = self._index_action_psy(self.combiActionVariable[index2])
													name = "ELU_STR 1.35G + " + str(self.coefq) + self.combiActionVariable[nAct] + " + " + str(round(self.coef_psy[typePsy][0] * self.coefq,2)) + self.combiActionVariable[index] + " + " + str(round(self.coef_psy[typePsy2][0] * self.coefq,2)) + self.combiActionVariable[index2]
													self._create_list_combination(name)
													
													if self.list_load[nload][2] == "Permanente G":
														value = self.coefg[1] * self.list_load[nload][4]
														array_load = self._create_array_load(name, value, nload, array_load)
													
													elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index2]:
														value = float(self.coef_psy[typePsy2][0]) * self.coefq  * self.list_load[nload][4]
														array_load = self._create_array_load(name, value, nload, array_load)

													elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
														value = float(self.coef_psy[typePsy][0]) * self.coefq  * self.list_load[nload][4]
														array_load = self._create_array_load(name, value, nload, array_load)

													elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
														value = self.coefq * self.list_load[nload][4]
														array_load = self._create_array_load(name, value, nload, array_load)


		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELUSTR = self._create_dataframe_load(array_load)
		
	
	def _return_combi_ELUSTR(self, combi):
		return self.df_load_ELUSTR.loc[self.df_load_ELUSTR['Combinaison']==combi]
	

	def _els_C(self):
		lenListLoad = len(self.list_load)
		array_load = np.empty((0, 8))

		for nload in range(lenListLoad):
			for nAct in range(5):

				if self.list_load[nload][2]== "Permanente G" and self.combiActionVariable[nAct] == "G":
					name = "ELS_C G"
					self._create_list_combination(name)
					value = self.list_load[nload][4]
					array_load = self._create_array_load(name, value, nload, array_load)
				
				elif self.combiActionVariable[nAct] != 0 and self.combiActionVariable[nAct] != "G":
					name = "ELS_C G + " + self.combiActionVariable[nAct]
					self._create_list_combination(name)

					if self.list_load[nload][2] != "Permanente G" and self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
						value = self.list_load[nload][4]
						array_load = self._create_array_load(name, value, nload, array_load)

					elif self.list_load[nload][2] == "Permanente G":
						value = self.list_load[nload][4]
						array_load = self._create_array_load(name, value, nload, array_load)
				
					for index in range(1,4):
							if self.combiActionVariable[nAct] != self.combiActionVariable[index] and self.combiActionVariable[nAct] != "W-":
								typePsy = self._index_action_psy(self.combiActionVariable[index])
								
								if self.combiActionVariable[index] != 0 and float(self.coef_psy[typePsy][0]) != 0:
									name = "ELS_C G + " + self.combiActionVariable[nAct] + " + " + str(round(self.coef_psy[typePsy][0],2)) + self.combiActionVariable[index]
									self._create_list_combination(name)

									if self.list_load[nload][2] == "Permanente G":
										value = self.list_load[nload][4]
										array_load = self._create_array_load(name, value, nload, array_load)

									elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
										value = float(self.coef_psy[typePsy][0]) * self.list_load[nload][4]
										array_load = self._create_array_load(name, value, nload, array_load)

									elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
										value = self.list_load[nload][4]
										array_load = self._create_array_load(name, value, nload, array_load)
									
									for index2 in range(1,4):
										if self.combiActionVariable[nAct] != self.combiActionVariable[index2]:
											if self.combiActionVariable[index] != self.combiActionVariable[index2]:
												if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
													if self.combiActionVariable[index2] != 0:
														typePsy2 = self._index_action_psy(self.combiActionVariable[index2])
														name = "ELS_C G + " + self.combiActionVariable[nAct] + " + " + str(round(self.coef_psy[typePsy][0],2)) + self.combiActionVariable[index] + " + " + str(round(self.coef_psy[typePsy2][0],2)) + self.combiActionVariable[index2]
														self._create_list_combination(name)

														if self.list_load[nload][2] == "Permanente G":
															value = self.list_load[nload][4]
															array_load = self._create_array_load(name, value, nload, array_load)
														
														elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index2]:
															value = float(self.coef_psy[typePsy2][0]) * self.list_load[nload][4]
															array_load = self._create_array_load(name, value, nload, array_load)

														elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
															value = float(self.coef_psy[typePsy][0]) * self.list_load[nload][4]
															array_load = self._create_array_load(name, value, nload, array_load)

														elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
															value = self.list_load[nload][4]
															array_load = self._create_array_load(name, value, nload, array_load)

		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELScarac = self._create_dataframe_load(array_load)
	

	def _return_combi_ELScarac(self, combi):
		return self.df_load_ELScarac.loc[self.df_load_ELScarac['Combinaison']==combi]


	def _els_QP(self):
		lenListLoad = len(self.list_load)
		array_load = np.empty((0, 8))

		for nload in range(lenListLoad):
			for nAct in range(5):

				if self.list_load[nload][2]== "Permanente G" and self.combiActionVariable[nAct] == "G":
					name = "ELS_QP G"
					self._create_list_combination(name)
					value = self.list_load[nload][4]
					array_load = self._create_array_load(name, value, nload, array_load)
				
				elif self.combiActionVariable[nAct] != 0 and self.combiActionVariable[nAct] != "G":

					typePsy = self._index_action_psy(self.combiActionVariable[nAct])
					if float(self.coef_psy[typePsy][2]) != 0:
						name = "ELS_QP G + " +  str(round(self.coef_psy[typePsy][2],2)) + self.combiActionVariable[nAct]
						self._create_list_combination(name)

						if self.list_load[nload][2] != "Permanente G" and self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
							value = float(self.coef_psy[typePsy][2]) * self.list_load[nload][4]
							array_load = self._create_array_load(name, value, nload, array_load)

						elif self.list_load[nload][2] == "Permanente G":
							value = self.list_load[nload][4]
							array_load = self._create_array_load(name, value, nload, array_load)
				
						for index in range(1,4):
								if self.combiActionVariable[nAct] != self.combiActionVariable[index] and self.combiActionVariable[nAct] != "W-":
									if self.combiActionVariable.index(self.combiActionVariable[nAct]) < self.combiActionVariable.index(self.combiActionVariable[index]):
										if self.combiActionVariable[index] != 0:
											typePsy1 = self._index_action_psy(self.combiActionVariable[index])
											if float(self.coef_psy[typePsy1][2]) != 0:
												name = "ELS_QP G + " +  str(round(self.coef_psy[typePsy][2],2)) + self.combiActionVariable[nAct] + " + " + str(round(self.coef_psy[typePsy1][2],2)) + self.combiActionVariable[index]
												self._create_list_combination(name)
												if self.list_load[nload][2] == "Permanente G":
													value = self.list_load[nload][4]
													array_load = self._create_array_load(name, value, nload, array_load)

												elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
													value = float(self.coef_psy[typePsy1][2]) * self.list_load[nload][4]
													array_load = self._create_array_load(name, value, nload, array_load)

												elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
													value = float(self.coef_psy[typePsy][2]) * self.list_load[nload][4]
													array_load = self._create_array_load(name, value, nload, array_load)
												
												for index2 in range(1,4):
													if self.combiActionVariable[nAct] != self.combiActionVariable[index2]:
														if self.combiActionVariable[index] != self.combiActionVariable[index2]:
															if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
																if self.combiActionVariable[index2] != 0:
																	typePsy2 = self._index_action_psy(self.combiActionVariable[index2])
																	if float(self.coef_psy[typePsy2][2]) != 0:
																		name = "ELS_QP G + " + self.combiActionVariable[nAct] + " + " + str(round(self.coef_psy[typePsy1][2],2)) + self.combiActionVariable[index] + " + " + str(round(self.coef_psy[typePsy2][2],2)) + self.combiActionVariable[index2]
																		self._create_list_combination(name)

																		if self.list_load[nload][2] == "Permanente G":
																			value = self.list_load[nload][4]
																			array_load = self._create_array_load(name, value, nload, array_load)
																		
																		elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index2]:
																			value = float(self.coef_psy[typePsy2][2]) * self.list_load[nload][4]
																			array_load = self._create_array_load(name, value, nload, array_load)

																		elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
																			value = float(self.coef_psy[typePsy1][2]) * self.list_load[nload][4]
																			array_load = self._create_array_load(name, value, nload, array_load)

																		elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
																			value = float(self.coef_psy[typePsy][2]) * self.list_load[nload][4]
																			array_load = self._create_array_load(name, value, nload, array_load)

		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELSqp = self._create_dataframe_load(array_load)


	def _return_combi_ELSqp(self, combi):
		return self.df_load_ELSqp.loc[self.df_load_ELSqp['Combinaison']==combi]


	@property
	def list_combination(self):
		"""Retourne un data frame avec toute les combinaison créer

		Returns:
			_type_: _description_
		"""
		self.name_combination.sort()
		return pd.DataFrame(self.name_combination, columns=['Combinaison'])
	
	
	def choice_combi_df(self, elu: bool=("True", "False"), _els_C: bool=("True", "False"), _els_QP: bool=("True", "False")):
		shape = len(self.list_combination)
		dict_load_combi = {}

		for i in range(shape):
			combi = self.list_combination.iloc[i,0]

			if elu and combi[0:7] == 'ELU_STR':
				df_combi = self._return_combi_ELUSTR(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi

			if _els_C and combi[0:5] == 'ELS_C':
				df_combi = self._return_combi_ELScarac(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi
			
			if _els_QP and combi[0:6] == 'ELS_QP':
				df_combi = self._return_combi_ELSqp(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi
		
		return dict_load_combi


	def choice_combi_list(self, list_combination:list, elu: bool=("True", "False"), _els_C: bool=("True", "False"), _els_QP: bool=("True", "False")):
		shape = len(list_combination)
		dict_load_combi = {}

		for i in range(shape):
			combi = list_combination[i]

			if elu and combi[0:7] == 'ELU_STR':
				df_combi = self._return_combi_ELUSTR(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi

			if _els_C and combi[0:5] == 'ELS_C':
				df_combi = self._return_combi_ELScarac(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi
			
			if _els_QP and combi[0:6] == 'ELS_QP':
				df_combi = self._return_combi_ELSqp(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi

		return dict_load_combi


class Calcul_EC0(Combinaison):
	def __init__(self, elu_str: bool = ("True", "False"), els_c: bool = ("True", "False"), els_qp: bool = ("True", "False"), **kwargs):
		super().__init__(**kwargs)
		self.dict_combi = self.choice_combi_df(elu_str, els_c, els_qp)

	def __str__(self) -> str:
		return self.dict_combi

	def min_type_load(self, name_combi:str) -> str:
		"""Retourne le type de chargment de plus courte durée

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
					if action == "Sn" and self.alt >= 1000 :
						name_load_type = dictName["Q"]
					elif action == "Sn" and self.cat == "Cat H : toits":
						name_load_type = dictName["Sn"]
					else:
						name_load_type = dictName[action]
		return name_load_type


		

if __name__== "__main__":

	list_load = [[1, '', 'Permanente G', 'Linéique', -10, '0/100', 'Z'],
				 [0, 'Poids propre', 'Permanente G', 'Linéique', -36, '0/100', 'Z'],
				 [2, '', 'Neige normale Sn', 'Linéique', -200, '0/100', 'Z'],
				 [3, '', 'Exploitation Q', 'Linéique', -150, '0/100', 'Z']]
	c1 = Calcul_EC0(list_load=list_load, alt=1001, cat="Cat A : habitation", h_bat=5, d_bat=15, b_bat=13.1, alphatoit=15)
	rcombi = "ELU_STR 1.35G + 1.5Sn"
	print(c1._return_combi_ELUSTR(rcombi))
	print(c1.coef_psy)
	print(c1.df_load_ELUSTR)


	
	