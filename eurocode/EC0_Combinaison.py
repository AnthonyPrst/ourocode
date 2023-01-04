# coding in UTF-8 
# by Anthony PARISOT

import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.getcwd(), "eurocode"))
from A0_Projet import Project

class Combinaison(Project):
	"""  """
	def __init__(self, list_load:list, cat="Categorie A", **kwargs):
		"""Créer une classe Combinaison qui génère les combinaisons d'action suivant les actions données.
		Cette classe est hérité de la classe Project du module A0_Projet.py

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
		self.dicoCombiAct = {"Permanente G": "G", "Exploitation Q": "Q",
							 "Neige S": "S", "Vent pression W+": "W+",
							 "Vent dépression W-": "W-"
							 }

		self.combiActionVariable = [0]*5 
		for nload in range(len(self.list_load)):
			if self.list_load[nload][2] == "Permanente G":
				self.combiActionVariable[0] = "G"

			elif self.list_load[nload][2] == "Exploitation Q":
				self.combiActionVariable[1] = "Q"

			elif self.list_load[nload][2] == "Neige S":
				self.combiActionVariable[2] = "S"

			elif self.list_load[nload][2] == "Vent pression W+":
				self.combiActionVariable[3] = "W+"

			else:
				self.combiActionVariable[4] = "W-"


		self.listPsy = self.coefPsy()
		self.eluSTR()
		self.elsC()
		self.elsQP()
		
	
	def data_from_csv(self, data_file:str):
		""" Retourne un dataframe d'un fichier CSV """
		repertory = os.getcwd() + "/data/" + data_file
		data_csv = pd.read_csv(repertory, sep=';', index_col=0)
		return data_csv   


	def coefPsy(self):
		""" Retourne les caractéristiques psy sous forme de liste python """
		listPsy = []
		listval = []
		data_csv_psy = self.data_from_csv("coeff_psy.csv")
		if self.cat != 'Aucune':
			for i in range(3):
				listval.append(data_csv_psy.loc[self.cat].loc["psy"+str(i)])
			listPsy.append(listval)
			listval = []
		else: 
			listPsy.append([None, None, None])
		
		for i in range(3):
			if self.alt>1000: 
				listval.append(data_csv_psy.loc["Neige > 1000m"].loc["psy"+str(i)])
			else:
				listval.append(data_csv_psy.loc["Neige < 1000m"].loc["psy"+str(i)])
		listPsy.append(listval)
		listval = []
		for i in range(3):      
			listval.append(data_csv_psy.loc["Vent"].loc["psy"+str(i)])
		listPsy.append(listval)
		
		return listPsy


	def index_action_psy(self ,action_variable):
		if action_variable == "Q":
			index = 0
		elif action_variable == "S":
			index = 1
		elif action_variable == "W+" or action_variable == "W-" :
			index = 2
		else:
			index = 0
		return index

	
	def create_array_load(self, name, value, i, array):
		load = np.array([name, self.list_load[i][0], self.list_load[i][1], self.list_load[i][2], self.list_load[i][3], value, self.list_load[i][5], self.list_load[i][6]], dtype=object)
		array= np.append(array,[load],axis= 0)
		return array

	def create_dataframe_load(self, array_load):
		df_load = pd.DataFrame(array_load, columns=['Combinaison', 'Index', 'Nom', 'Action', 'Type', 'Valeur', 'Position', 'Axe'])
		#print(df_load)
		return(df_load)
	

	def create_list_combination(self, nameCombi):
		test = 1
		for i in self.name_combination:
			if i == nameCombi:
				test = 0
				break
		if test:
			self.name_combination.append(nameCombi)
			
	
	def eluSTR(self):
		""" Combinaison à l'ELU STR """
		lenListLoad = len(self.list_load)
		array_load = np.empty((0, 8))

		for nload in range(lenListLoad):
			name = "ELU_STR " + self.dicoCombiAct[self.list_load[nload][2]]
			self.create_list_combination(name)
			value = self.list_load[nload][4]
			array_load = self.create_array_load(name, value, nload, array_load)
			
			for nAct in range(5):
			
				if self.list_load[nload][2]== "Permanente G" and self.combiActionVariable[nAct] == "G":
					# if self.list_load[nload][2]== "Vent dépression W-" and self.combiActionVariable[nAct] == "W-":
					#     name = "ELU_STR " + str(self.coefg[0])+"G" + str(self.coefq) + self.combiActionVariable[nAct]
					#     self.create_list_combination(name)
					#     value = self.coefg[0] * self.list_load[nload][4]
					#     array_load = self.create_array_load(name, value, nload, array_load)
					# else :
						name = "ELU_STR " + str(self.coefg[1])+"G"
						self.create_list_combination(name)
						value = self.coefg[1] * self.list_load[nload][4]
						array_load = self.create_array_load(name, value, nload, array_load)
				
				elif self.combiActionVariable[nAct] != 0 and self.combiActionVariable[nAct] != "G":
					if self.combiActionVariable[nAct] != "W-":
						name = "ELU_STR 1.35G + " + str(self.coefq) + self.combiActionVariable[nAct]
						self.create_list_combination(name)
					
					if self.list_load[nload][2] != "Permanente G" and self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
						if self.combiActionVariable[nAct] == "W-":
							name = "ELU_STR 1G" + " + " + str(self.coefq) + self.combiActionVariable[nAct]
							self.create_list_combination(name)
						
						value = self.coefq * self.list_load[nload][4]
						array_load = self.create_array_load(name, value, nload, array_load)

					elif self.list_load[nload][2] == "Permanente G":
						if self.combiActionVariable[nAct] == "W-":
							name = "ELU_STR 1G" + " + " + str(self.coefq) + self.combiActionVariable[nAct]
							self.create_list_combination(name)
							value = self.coefg[0] * self.list_load[nload][4]
						else:
							value = self.coefg[1] * self.list_load[nload][4]
						array_load = self.create_array_load(name, value, nload, array_load)
					
					for index in range(1,4):
						if self.combiActionVariable[nAct] != self.combiActionVariable[index] and self.combiActionVariable[nAct] != "W-":
							if self.combiActionVariable[index] != 0:
								typePsy = self.index_action_psy(self.combiActionVariable[index])
								name = "ELU_STR 1.35G + " + str(self.coefq) + self.combiActionVariable[nAct] + " + " + str(round(self.listPsy[typePsy][0] * self.coefq,2)) + self.combiActionVariable[index]
								self.create_list_combination(name)
								
								if self.list_load[nload][2] == "Permanente G":
									value = self.coefg[1] * self.list_load[nload][4]
									array_load = self.create_array_load(name, value, nload, array_load)

								elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
									value = float(self.listPsy[typePsy][0]) * self.coefq  * self.list_load[nload][4]
									array_load = self.create_array_load(name, value, nload, array_load)

								elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
									value = self.coefq * self.list_load[nload][4]
									array_load = self.create_array_load(name, value, nload, array_load)
								
								for index2 in range(1,4):
									if self.combiActionVariable[nAct] != self.combiActionVariable[index2]:
										if self.combiActionVariable[index] != self.combiActionVariable[index2]:
											if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
												if self.combiActionVariable[index2] != 0:
													typePsy2 = self.index_action_psy(self.combiActionVariable[index2])
													name = "ELU_STR 1.35G + " + str(self.coefq) + self.combiActionVariable[nAct] + " + " + str(round(self.listPsy[typePsy][0] * self.coefq,2)) + self.combiActionVariable[index] + " + " + str(round(self.listPsy[typePsy2][0] * self.coefq,2)) + self.combiActionVariable[index2]
													self.create_list_combination(name)
													
													if self.list_load[nload][2] == "Permanente G":
														value = self.coefg[1] * self.list_load[nload][4]
														array_load = self.create_array_load(name, value, nload, array_load)
													
													elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index2]:
														value = float(self.listPsy[typePsy2][0]) * self.coefq  * self.list_load[nload][4]
														array_load = self.create_array_load(name, value, nload, array_load)

													elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
														value = float(self.listPsy[typePsy][0]) * self.coefq  * self.list_load[nload][4]
														array_load = self.create_array_load(name, value, nload, array_load)

													elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
														value = self.coefq * self.list_load[nload][4]
														array_load = self.create_array_load(name, value, nload, array_load)


		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELUSTR = self.create_dataframe_load(array_load)
		
	
	def _return_combi_ELUSTR(self, combi):
		return self.df_load_ELUSTR.loc[self.df_load_ELUSTR['Combinaison']==combi]
	

	def elsC(self):
		lenListLoad = len(self.list_load)
		array_load = np.empty((0, 8))

		for nload in range(lenListLoad):
			for nAct in range(5):

				if self.list_load[nload][2]== "Permanente G" and self.combiActionVariable[nAct] == "G":
					name = "ELS_C G"
					self.create_list_combination(name)
					value = self.list_load[nload][4]
					array_load = self.create_array_load(name, value, nload, array_load)
				
				elif self.combiActionVariable[nAct] != 0 and self.combiActionVariable[nAct] != "G":
					name = "ELS_C G + " + self.combiActionVariable[nAct]
					self.create_list_combination(name)

					if self.list_load[nload][2] != "Permanente G" and self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
						value = self.list_load[nload][4]
						array_load = self.create_array_load(name, value, nload, array_load)

					elif self.list_load[nload][2] == "Permanente G":
						value = self.list_load[nload][4]
						array_load = self.create_array_load(name, value, nload, array_load)
				
					for index in range(1,4):
							if self.combiActionVariable[nAct] != self.combiActionVariable[index] and self.combiActionVariable[nAct] != "W-":
								typePsy = self.index_action_psy(self.combiActionVariable[index])
								
								if self.combiActionVariable[index] != 0 and float(self.listPsy[typePsy][0]) != 0:
									name = "ELS_C G + " + self.combiActionVariable[nAct] + " + " + str(round(self.listPsy[typePsy][0],2)) + self.combiActionVariable[index]
									self.create_list_combination(name)

									if self.list_load[nload][2] == "Permanente G":
										value = self.list_load[nload][4]
										array_load = self.create_array_load(name, value, nload, array_load)

									elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
										value = float(self.listPsy[typePsy][0]) * self.list_load[nload][4]
										array_load = self.create_array_load(name, value, nload, array_load)

									elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
										value = self.list_load[nload][4]
										array_load = self.create_array_load(name, value, nload, array_load)
									
									for index2 in range(1,4):
										if self.combiActionVariable[nAct] != self.combiActionVariable[index2]:
											if self.combiActionVariable[index] != self.combiActionVariable[index2]:
												if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
													if self.combiActionVariable[index2] != 0:
														typePsy2 = self.index_action_psy(self.combiActionVariable[index2])
														name = "ELS_C G + " + self.combiActionVariable[nAct] + " + " + str(round(self.listPsy[typePsy][0],2)) + self.combiActionVariable[index] + " + " + str(round(self.listPsy[typePsy2][0],2)) + self.combiActionVariable[index2]
														self.create_list_combination(name)

														if self.list_load[nload][2] == "Permanente G":
															value = self.list_load[nload][4]
															array_load = self.create_array_load(name, value, nload, array_load)
														
														elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index2]:
															value = float(self.listPsy[typePsy2][0]) * self.list_load[nload][4]
															array_load = self.create_array_load(name, value, nload, array_load)

														elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
															value = float(self.listPsy[typePsy][0]) * self.list_load[nload][4]
															array_load = self.create_array_load(name, value, nload, array_load)

														elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
															value = self.list_load[nload][4]
															array_load = self.create_array_load(name, value, nload, array_load)

		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELScarac = self.create_dataframe_load(array_load)
	

	def _return_combi_ELScarac(self, combi):
		return self.df_load_ELScarac.loc[self.df_load_ELScarac['Combinaison']==combi]


	def elsQP(self):
		lenListLoad = len(self.list_load)
		array_load = np.empty((0, 8))

		for nload in range(lenListLoad):
			for nAct in range(5):

				if self.list_load[nload][2]== "Permanente G" and self.combiActionVariable[nAct] == "G":
					name = "ELS_QP G"
					self.create_list_combination(name)
					value = self.list_load[nload][4]
					array_load = self.create_array_load(name, value, nload, array_load)
				
				elif self.combiActionVariable[nAct] != 0 and self.combiActionVariable[nAct] != "G":

					typePsy = self.index_action_psy(self.combiActionVariable[nAct])
					if float(self.listPsy[typePsy][2]) != 0:
						name = "ELS_QP G + " +  str(round(self.listPsy[typePsy][2],2)) + self.combiActionVariable[nAct]
						self.create_list_combination(name)

						if self.list_load[nload][2] != "Permanente G" and self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
							value = float(self.listPsy[typePsy][2]) * self.list_load[nload][4]
							array_load = self.create_array_load(name, value, nload, array_load)

						elif self.list_load[nload][2] == "Permanente G":
							value = self.list_load[nload][4]
							array_load = self.create_array_load(name, value, nload, array_load)
				
						for index in range(1,4):
								if self.combiActionVariable[nAct] != self.combiActionVariable[index] and self.combiActionVariable[nAct] != "W-":
									if self.combiActionVariable.index(self.combiActionVariable[nAct]) < self.combiActionVariable.index(self.combiActionVariable[index]):
										if self.combiActionVariable[index] != 0:
											typePsy1 = self.index_action_psy(self.combiActionVariable[index])
											if float(self.listPsy[typePsy1][2]) != 0:
												name = "ELS_QP G + " +  str(round(self.listPsy[typePsy][2],2)) + self.combiActionVariable[nAct] + " + " + str(round(self.listPsy[typePsy1][2],2)) + self.combiActionVariable[index]
												self.create_list_combination(name)
												if self.list_load[nload][2] == "Permanente G":
													value = self.list_load[nload][4]
													array_load = self.create_array_load(name, value, nload, array_load)

												elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
													value = float(self.listPsy[typePsy1][2]) * self.list_load[nload][4]
													array_load = self.create_array_load(name, value, nload, array_load)

												elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
													value = float(self.listPsy[typePsy][2]) * self.list_load[nload][4]
													array_load = self.create_array_load(name, value, nload, array_load)
												
												for index2 in range(1,4):
													if self.combiActionVariable[nAct] != self.combiActionVariable[index2]:
														if self.combiActionVariable[index] != self.combiActionVariable[index2]:
															if self.combiActionVariable.index(self.combiActionVariable[index]) < self.combiActionVariable.index(self.combiActionVariable[index2]):
																if self.combiActionVariable[index2] != 0:
																	typePsy2 = self.index_action_psy(self.combiActionVariable[index2])
																	if float(self.listPsy[typePsy2][2]) != 0:
																		name = "ELS_QP G + " + self.combiActionVariable[nAct] + " + " + str(round(self.listPsy[typePsy1][2],2)) + self.combiActionVariable[index] + " + " + str(round(self.listPsy[typePsy2][2],2)) + self.combiActionVariable[index2]
																		self.create_list_combination(name)

																		if self.list_load[nload][2] == "Permanente G":
																			value = self.list_load[nload][4]
																			array_load = self.create_array_load(name, value, nload, array_load)
																		
																		elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index2]:
																			value = float(self.listPsy[typePsy2][2]) * self.list_load[nload][4]
																			array_load = self.create_array_load(name, value, nload, array_load)

																		elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[index]:
																			value = float(self.listPsy[typePsy1][2]) * self.list_load[nload][4]
																			array_load = self.create_array_load(name, value, nload, array_load)

																		elif self.dicoCombiAct[self.list_load[nload][2]] == self.combiActionVariable[nAct]:
																			value = float(self.listPsy[typePsy][2]) * self.list_load[nload][4]
																			array_load = self.create_array_load(name, value, nload, array_load)

		array_load = array_load[array_load[:, 0].argsort()]
		self.df_load_ELSqp = self.create_dataframe_load(array_load)


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
	
	
	def choice_combi_df(self, elu=True, elsc=True, elsqp=True):
		shape = len(self.list_combination)
		dict_load_combi = {}

		for i in range(shape):
			combi = self.list_combination.iloc[i,0]

			if elu and combi[0:7] == 'ELU_STR':
				df_combi = self._return_combi_ELUSTR(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi

			if elsc and combi[0:5] == 'ELS_C':
				df_combi = self._return_combi_ELScarac(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi
			
			if elsqp and combi[0:6] == 'ELS_QP':
				df_combi = self._return_combi_ELSqp(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi
		
		return dict_load_combi


	def choice_combi_list(self, list_combination, elu=True, elsc=True, elsqp=True):
		shape = len(list_combination)
		dict_load_combi = {}

		for i in range(shape):
			combi = list_combination[i]

			if elu and combi[0:7] == 'ELU_STR':
				df_combi = self._return_combi_ELUSTR(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi

			if elsc and combi[0:5] == 'ELS_C':
				df_combi = self._return_combi_ELScarac(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi
			
			if elsqp and combi[0:6] == 'ELS_QP':
				df_combi = self._return_combi_ELSqp(combi)
				df_combi = df_combi.drop('Combinaison', 1)
				load_list_combi = df_combi.values.tolist()
				dict_load_combi[combi] =  load_list_combi

		return dict_load_combi


class Calcul_EC0(Combinaison):
	def __init__(self, elu_str: bool = True, els_c: bool = True, els_qp: bool = True, **kwargs):
		super().__init__(**kwargs)
		self.dict_combi = self.choice_combi_df(elu_str, els_c, els_qp)


	def min_type_load(self, name_combi):
		dictName = {"G": "Permanente",
					"Q": "Moyen terme",
					"S": "Court terme",
					"W+": "Instantanee", 
					"W-": "Instantanee"}

		for action in self.combiActionVariable:
			if action:
				indexAction = name_combi[8:].find(action)
				if indexAction > -1 :
					if action == "S" and self.alt >= 1000 :
						name_load_type = dictName["Q"]
					else:
						name_load_type = dictName[action]
		return name_load_type


		

if __name__== "__main__":

	list_load = [[1, '', 'Permanente G', 'Linéique', -10, '0/100', 'Z'],
				 [0, 'Poids propre', 'Permanente G', 'Linéique', -36, '0/100', 'Z'],
				 [2, '', 'Neige S', 'Linéique', -200, '0/100', 'Z'],
				 [3, '', 'Exploitation Q', 'Linéique', -150, '0/100', 'Z']]
	c1 = Calcul_EC0(list_load=list_load, alt=1001, cat="Cat A : habitation", h_bat=5, d_bat=15, b_bat=13.1, alphatoit=15)
	rcombi = "ELU_STR 1.35G + 1.5S"
	print(c1._return_combi_ELUSTR(rcombi))

	
	