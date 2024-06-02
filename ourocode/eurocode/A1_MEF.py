#! env/scripts/python.exe
#! env/scripts/python.exe
#coding in UTF8

import os, sys

import time
import numpy as np
from matplotlib import pyplot as plt
import math as mt
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

sys.path.append(os.path.join(os.getcwd(), "ourocode"))
from eurocode.A0_Projet import Bar_generator
from eurocode.EC0_Combinaison import Combinaison

# from ourocode.eurocode.EC0_Combinaison import Combinaison
# from ourocode.eurocode.A0_Projet import Bar_generator

class _Base_graph(object): 
    """ Retourne un diagramme de base """
    def __init__(self):
        self.name = "name_combi"
        self.title = ""
        self.color = "darkcyan"
        self.x_label = "Longueur (mm)"
        self.y_label = "Effort (kN)"
        self.unit = " kN"
        self.x_values = []
        self.y_values = []
        self.max_XY_values = (0,0)
        self.min_XY_values = (0,0)
        self.excent_X = 1
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = []
        self.save_path = None
    
    def _base_graphique(self):
        """ Retourne le diagramme """
        # plt.clf()  # Effacer le graphique précédent
        
        plt.figure(self.name, figsize=(12,4))
        plt.gcf().subplots_adjust(left = 0.1, bottom = 0.25, right = 0.9,
                                top = 0.75, wspace = 0, hspace = 0.95)

        #manager = plt.get_current_fig_manager()
        #manager.resize(*manager.window.maxsize())
        
        plt.plot(self.x_values, self.y_values, color=self.color)
        plt.title(self.title + self.name, color=self.color)
        plt.ylabel(self.y_label)
        plt.xlabel (self.x_label)
        
        if self.fill_between:
            plt.fill_between(self.fill_between_X, self.y_values, 0, color=self.color, alpha = 0.2)
            
        plt.annotate('Max: '+ str(round(self.max_XY_values[1],2)) + self.unit, xy=(self.max_XY_values[0], self.max_XY_values[1]),
                    xytext=(self.max_XY_values[0]-self.excent_X, self.max_XY_values[1]-self.excent_Y),
                    arrowprops=dict(facecolor='black',
                    arrowstyle='->'),
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=self.color, lw=1, alpha=0.7))
        plt.annotate('Min: '+ str(round(self.min_XY_values[1],2)) + self.unit, xy=(self.min_XY_values[0], self.min_XY_values[1]),
                    xytext= (self.min_XY_values[0]-self.excent_X, self.min_XY_values[1]+self.excent_Y),
                    arrowprops= dict(facecolor='black',
                    arrowstyle= '->'),
                    bbox= dict(boxstyle="round,pad=0.3", fc="white", ec=self.color, lw=1, alpha=0.7))
        plt.grid()
        
        plt.savefig(self.save_path)


class _Triplets(object):
    def __init__(self) -> None:
        self.data = ([], ([], []))
    
    def __str__(self) -> str:
        return str(self.data)
    
    def append(self, I, J, val):
        self.data[0].append(val)
        self.data[1][0].append(I)
        self.data[1][1].append(J)
    
    
        
class MEF(Bar_generator, _Base_graph):
    SAVE_PATH = os.path.join(Combinaison.PATH_CATALOG, "data","screenshot")
    def __init__(self, combinaison: Combinaison, **kwargs):
        """Classe permettant de créer des poutres MEF et de les calculer. Cette classe est hérité de l'objet Combinaison du module EC0_Combinaison.py
        Args:
            A (float): Section en mm²
        """
        super(Bar_generator, self).__init__(**kwargs)
        _Base_graph.__init__(self)
        self.combinaison = combinaison

    def get_lenght_of_element(self, element_index: int):
        """Retourne la longueur d'un élément en mm à partir de ces coordonnées globale

        Args:
            element_index (int): index de l'élément à récupérer
        """
        element = self.element_list[element_index]
        print(element)
        node_1, node_2 = element
        x1, y1 = self.node_coor[int(node_1)-1]
        x2, y2 = self.node_coor[int(node_2)-1]
        v1 = (x2-x1, y2-y1)
        print(x1, y1, x2, y2)
        return abs(mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2))
    
    
    def _transformation_matrix(self, angle: float):
        """Matrice de transformation pour un élément barre 3D dans un plan perpendiculaire à un axe.
        Si la barre est réellement dans un espace 3D il faut alors calcul les produits scalaire de chaque axe dans une matrice lambda.
        Voir la video youtube https://www.youtube.com/watch?v=J2MwWFcY3NI&list=PLmw2x4fCMxg6Yw-7FqwBNSobXZR6fQ3Sm&index=16
        Args:
            angle (float): angle en degré de la barre par rapport à X
        """
        Assembleur_COO_T = _Triplets()
        rad_angle = mt.radians(angle)
        # lambda_matrix = np.array([[mt.cos(rad_angle), 1, mt.sin(rad_angle)],
        #                         [-mt.sin(rad_angle), 1, mt.cos(rad_angle)],
        #                         [0, 0, 1]])
        lambda_matrix = np.array([[mt.cos(rad_angle), 0, mt.sin(rad_angle)],
                                [0, 1, 0],
                                [-mt.sin(rad_angle), 0, mt.cos(rad_angle)]])
        transformation_matrix = np.zeros((12,12))
        index = 0
        for _ in range(0,4):
            for i in range(0,3):
                for j in range(0,3):
                    Assembleur_COO_T.append(index+i, index+j, lambda_matrix[i,j])
            index += 3

        transformation_matrix = (coo_matrix(Assembleur_COO_T.data, shape=(12, 12))).tocsr()
        # print("Matrice de transformation locale/globale: ",transformation_matrix)
        return transformation_matrix
    
    def _apply_relaxation_to_k(self, k, relaxation, position):
        indices = {
            "start": [0, 1, 2, 3, 4, 5],
            "end": [6, 7, 8, 9, 10, 11]
        }
        
        dof_map = {
            "u": 0,
            "v": 1,
            "w": 2,
            "teta_x": 3,
            "teta_y": 4,
            "teta_z": 5
        }
        
        if position == "start":
            dof_indices = indices["start"]
        elif position == "end":
            dof_indices = indices["end"]
        else:
            raise ValueError("Position must be 'start' or 'end'")
        
        for key, rel in relaxation.items():
            if key != "position" and rel:
                index = dof_map[key]
                k[dof_indices[index], :] = 0
                k[:, dof_indices[index]] = 0
        return k
    

    def _apply_relaxation_to_U(self, element: list, relaxation):
        dof_map = {
            "u": 0,
            "v": 1,
            "w": 2,
            "teta_x": 3,
            "teta_y": 4,
            "teta_z": 5
        }

        if relaxation["position"] == "start":
            index_U = (int(element[0])-1) * 6
        elif relaxation["position"] == "end":
            index_U = (int(element[1])-1) * 6
        else:
            raise ValueError("Position must be 'start' or 'end'")
        print(element)
        print(index_U)
        for key, rel in relaxation.items():
            if key != "position" and rel:
                self.matriceU[index_U+dof_map[key], 0] = True
    

    def _init_matrix_K(self, index: int|str) -> np.array:
        """Initialise la matrice K suivant la longueur de l'élément l en mm
        Args:
            l (int): longueur de l'élément MEF en mm
        Returns:
            numpy.array: matrice de rigidité K
        """
        for key, bar in self.bar_info.items():
            if index in bar["elements"]:
                angle = bar["angle"]
                l = bar["length"] / len(bar["elements"])
                E = bar["material"]["E"]
                EA = E * bar["section"] 
                EIy = E * bar["material"]["Iy"]
                EIz = E * bar["material"]["Iz"]
                GJ = bar["material"]["G"] * bar["material"]["J"]
                print(f"J'ai retrouvé la poutre {key} qui est associé à l'élément {index}, son angle est de {angle} degrés.")
                break
        
        k = np.array([[(EA)/l, 0, 0, 0, 0, 0, -(EA)/l, 0, 0, 0, 0, 0],
                    [0, (12*EIz)/l**3, 0, 0, 0, (6*EIz)/l**2, 0, -(12*EIz)/l**3, 0, 0, 0, (6*EIz)/l**2],
                    [0, 0, (12*EIy)/l**3, 0, -(6*EIy)/l**2, 0, 0, 0, -(12*EIy)/l**3, 0, -(6*EIy)/l**2, 0],
                    [0, 0, 0, GJ/l, 0, 0, 0, 0, 0, -GJ/l, 0, 0],
                    [0, 0, -(6*EIy)/l**2, 0, (4*EIy)/l, 0, 0, 0, (6*EIy)/l**2, 0, (2*EIy)/l, 0],
                    [0, (6*EIz)/l**2, 0, 0, 0, (4*EIz)/l, 0, -(6*EIz)/l**2, 0, 0, 0, (2*EIz)/l],
                    [-(EA)/l, 0, 0, 0, 0, 0, (EA)/l, 0, 0, 0, 0, 0],
                    [0, -(12*EIz)/l**3, 0, 0, 0, -(6*EIz)/l**2, 0, (12*EIz)/l**3, 0, 0, 0, -(6*EIz)/l**2],
                    [0, 0, -(12*EIy)/l**3, 0, (6*EIy)/l**2, 0, 0, 0, (12*EIy)/l**3, 0, (6*EIy)/l**2, 0],
                    [0, 0, 0, -GJ/l, 0, 0, 0, 0, 0, GJ/l, 0, 0],
                    [0, 0, -(6*EIy)/l**2, 0, (2*EIy)/l, 0, 0, 0, (6*EIy)/l**2, 0, (4*EIy)/l,0],
                    [0, (6*EIz)/l**2, 0, 0, 0, (2*EIz)/l, 0, -(6*EIz)/l**2, 0, 0, 0, (4*EIz)/l]])
        
        if bar.get("relaxation"):
            change_matrix_k = False
            for relaxation in bar["relaxation"]:
                match relaxation["position"]:
                    case "start":
                        if index == bar["elements"][0]:
                            change_matrix_k = True
                    case "end":
                        if index == bar["elements"][-1]:
                            change_matrix_k = True

                if change_matrix_k:       
                    k = self._apply_relaxation_to_k(k, relaxation, relaxation["position"])
                    # self._apply_relaxation_to_U(self.element_list[index], relaxation)
                    print("Matrice k local modifiée avec la relaxation: ", k)
                    break 


        T_matrix = self._transformation_matrix(angle)

        self.k_local["ele"+str(index)] = {"k_local": k, "T_global": T_matrix}
        return np.dot(np.dot(T_matrix.T.toarray(), k), T_matrix.toarray())


    def _agregation_matrix_K(self, assembleur_COO_K: list, k: np.array, n1: int, n2: int):
        # Agrégation de la matrice global
        for i in range(-6,0):
            for j in range(-6,0):
                assembleur_COO_K.append(n1+i, n1+j, k[i+6,j+6]+0.00001) # ajout de bruit pour éviter les singularité (+0.001)
                assembleur_COO_K.append(n2+i, n2+j, k[i+12,j+12]+0.00001) # ajout de bruit pour éviter les singularité (+0.001)
                if k[i+6,j+12]:
                    assembleur_COO_K.append(n1+i, n2+j, k[i+6,j+12])
                if k[i+12,j+6]:
                    assembleur_COO_K.append(n2+i, n1+j, k[i+12,j+6]+0.00001) # ajout de bruit pour éviter les singularité (+0.001)
        

    
    def _matrice_K_1D(self):
        """Méthodequi créer la matrice de rigidité K assemblée
        Returns:
            numpy.array: la matrice de rigidité assemblée
        """
        Assembleur_COO_K = _Triplets()
        self.k_local = {}

        for i in range(len(self.element_list)):
            # matrice de rigidité K
            k = self._init_matrix_K(i)
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6
            
            self._agregation_matrix_K(Assembleur_COO_K , k, n1, n2)
        
        self.matriceK = (coo_matrix(Assembleur_COO_K.data, shape=((len(self.element_list)+1)*6, (len(self.element_list)+1)*6))).tocsr() 
             
        # print("\n Matrice de rigidité K  : \n", self.matriceK.toarray(), "\n")
        print("\n Matrice de déplacement U après relaxation : \n",self.matriceU, "\n")
        return self.matriceK
    

    def local_to_global(self, x_local, bar):
        """Convert local x coordinate on a bar to global coordinates (X, Y)."""
        x1, y1 = bar['x1'], bar['y1']
        angle = np.radians(bar['angle'])  # Convert angle to radians

        X = x1 + x_local * np.cos(angle)
        Y = y1 + x_local * np.sin(angle)

        return np.array([X, Y])
    

    def global_to_local(self, global_coor, bar):
        """Convert global coordinates (X, Y) to local x coordinate on a bar."""
        x1, y1 = bar['x1'], bar['y1']
        x2, y2 = bar['x2'], bar['y2']
        
        bar_vector = np.array([x2 - x1, y2 - y1])
        point_vector = global_coor - np.array([x1, y1])
        
        bar_length = np.linalg.norm(bar_vector)
        local_x = np.dot(point_vector, bar_vector) / bar_length

        return local_x
    

    def _nearest_node(self, x_local, bar_id) -> list:
        """Find the nearest node in global coordinates from a local x distance on a given bar."""
        bar = self.bar_info[bar_id]
        value = self.local_to_global(x_local, bar)
        
        # Calculate the Euclidean distance between the given value and all node coordinates
        difference_array = np.linalg.norm(self.node_coor - value, axis=1)
        
        # Find the index of the minimum element from the array
        index = difference_array.argmin()
        nearest_node_coor = self.node_coor[index]
        
        # Convert nearest node's global coordinates to local x coordinate
        nearest_node_local_x = self.global_to_local(nearest_node_coor, bar)
        print("local x: ",nearest_node_local_x)
        
        # Compare local x coordinates
        if np.isclose(nearest_node_local_x, x_local, atol=1e-6):
            position = "exactly at"
        elif nearest_node_local_x < x_local:
            position = "before"
        else:
            position = "after"
        
        index_I_matrice = ((index+1) * 6) - 6
        element = self.element_list[index-1]

        print({"index": index, "nearest node": nearest_node_coor, "nearest node local x": nearest_node_local_x, "index matrice": index_I_matrice, "element": element, "position": position})
        return {"index": index, 
                "nearest node": nearest_node_coor, 
                "nearest node local x": nearest_node_local_x, 
                "index matrice": index_I_matrice, 
                "element": element, 
                "position": position}


    def _matrice_U_1D(self, dict_supports: dict):
        """ Donner en entrée le dictionnaire type des appuis extérieurs, 
        ressort une matrice acceptable pour le calcul MEF
        
        Args:
            dict_supports (dict): dictionnaire des appuis extérieures"""
        
        print(dict_supports)
        self.matriceU = np.ones(((len(self.element_list)+1)*6,1))
        
        for support in dict_supports.values():
            bar_id = support['N° barre']
            nearest_node = self._nearest_node(support["Position de l'appui"], bar_id)
           
            if support["Type d'appui"] == 'Rotule':
                for i in range(0,4):
                    self.matriceU[nearest_node['index matrice'] + i,0] = False
                    
            elif support["Type d'appui"] == 'Encastrement':
                for i in range(0,6):
                    self.matriceU[nearest_node['index matrice'] + i,0] = False

            elif support["Type d'appui"] == 'Simple Z':
                self.matriceU[nearest_node['index matrice'] + 2,0] = False
                self.matriceU[nearest_node['index matrice'] + 3,0] = False

            elif support["Type d'appui"] == 'Simple Y':
                self.matriceU[nearest_node['index matrice'] + 1,0] = False
                self.matriceU[nearest_node['index matrice'] + 3,0] = False

            elif support["Type d'appui"] == 'Simple X':
                self.matriceU[nearest_node['index matrice'],0] = False
                self.matriceU[nearest_node['index matrice'] + 3,0] = False

        print("\n Matrice de déplacement U : \n",self.matriceU, "\n")
        return self.matriceU


    def _matrice_Fext_1D(self, dict_F_ext: dict):
        """ Donner en entrée le dictionnaire type des efforts extérieurs en daN, 
        ressort une matrice acceptable pour le calcul MEF 
        
        Args:
            dict_F_ext (dict): dictionnaire des charges extérieures"""
        print(dict_F_ext)
        Assembleur_COO_F = _Triplets()
        
        self.alpha_R_tetaY = 0
        

        for index_F_ext, F_ext in dict_F_ext.items():
            bar_id = F_ext['N° barre']
            l = self.get_lenght_of_element(self.bar_info[bar_id]["elements"][0])
            print(l)
            angle = mt.radians(self.bar_info[bar_id]["angle"])
            base_global_R_tetaY = np.array([[np.cos(angle), np.sin(angle), 0],
                                            [-np.sin(angle), np.cos(angle), 0],
                                            [0, 0, 1]])
            listFlocal = np.zeros((3,1))
            listFlocal2 = np.zeros((3,1))
            
            if F_ext['Type de charge'] == 'Linéique':
            
                slash_pos = F_ext['Position'].index("/")
                start_pos = int(F_ext['Position'][0:slash_pos])
                start_index = self._nearest_node(start_pos, bar_id)
                
                end_pos = int(F_ext['Position'][slash_pos+1:])
                end_index = self._nearest_node(end_pos, bar_id)
                nearest_node = [[start_index, start_pos], [end_index, end_pos]]

                # print("Longueur élément:", l, " mm")
                for index in range(start_index[2], end_index[2]+6, 6):
                    listFlocal = np.zeros((3,1))
                    listFlocal2 = np.zeros((3,1))

                    if F_ext['Axe'] == 'Z' or F_ext['Axe'] == 'Z local' :
                        force = F_ext['Charge'] * (l/1000) * 10

                        listFlocal[1,0] = F_ext['Charge'] * 10 * (l/1000) / 2 # Rza et Rzb
                        listFlocal2[1,0] = listFlocal[1,0]
                        listFlocal[2,0] = F_ext['Charge'] * 10 * (l/1000) **2 / 12 # Mya et Myb
                        listFlocal2[2,0] = listFlocal[2,0]

                        if index in (start_index[2], end_index[2]):
                            calc_formulaire = False
                            if index != end_index[2]:
                                ind = nearest_node[0]

                                # Cas dans le quelle la force est avant le noeud le plus proche
                                if ind[1] < ind[0][1]:
                                    a = ind[0][1] - ind[1]
                                    b = ind[1] - self.node_coor[ind[0][0]-1][0]
                                    Mya = ((force * b**2)/l) * (4*(b/l) + 3*(b/l)**2)
                                    Myb = -((force * b**2)/12) * (6-8*(b/l) + 3*(b/l)**2)
                                    force = force * b
                                    b = b/2 
                                    a = a + b 
                                    Rza = (force * b**2 * (3 * a + b))/l**3
                                    Rzb = (force * a**2 * (3 * b + a))/l**3
                                    
                                    listFlocal[1,0] = Rza
                                    listFlocal[2,0] = Mya
                                    listFlocal2[1,0] = Rzb
                                    listFlocal2[2,0] = Myb
                                    calc_formulaire = True
                            
                            else:
                                ind = nearest_node[1]
                                # Formule tirée de l'aide mémoire page 113 cas 3 et cas 1 pour RA et RB
                                # Cas dans lequel la force est aprés le noeud le plus proche
                                if ind[1] > ind[0][1]:
                                    a = ind[1] - ind[0][1]
                                    b = self.node_coor[ind[0][0]+1][0]-ind[1]
                                    Mya = ((force * a**2)/12) * (6-8*(a/l) + 3*(a/l)**2)
                                    Myb = -((force * a**2)/l) * (4*(a/l) + 3*(a/l)**2)
                                    force = force * b
                                    a = a/2 
                                    b = a + b 
                                    Rza = (force * b**2 * (3 * a + b))/l**3
                                    Rzb = (force * a**2 * (3 * b + a))/l**3
                                    
                                    listFlocal[1,0] = Rza
                                    listFlocal[2,0] = Mya
                                    listFlocal2[1,0] = Rzb
                                    listFlocal2[2,0] = Myb
                                    calc_formulaire = True

                            if calc_formulaire:
                                if F_ext['Axe'] == 'Z':
                                    listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                                    listFlocal2= np.dot(base_global_R_tetaY, listFlocal2)

                                for i_index in range(3):
                                    j_index = i_index * 2
                                    Assembleur_COO_F.append(index+j_index, 0, listFlocal[i_index,0])
                                    Assembleur_COO_F.append(index+6+j_index, 0, listFlocal2[i_index,0])
                                continue

                        if F_ext['Axe'] == 'Z':
                            listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                            listFlocal2= np.dot(base_global_R_tetaY, listFlocal2)
                            
                    elif F_ext['Axe'] == 'X':
                        listFlocal[0,0] = F_ext['Charge'] * (l/1000) * 10 / 2
                        listFlocal2[0,0] = listFlocal[0,0]
                        listFlocal= np.dot(base_global_R_tetaY, listFlocal)

                    if index != end_index[2]:
                        for i_index in range(3):
                            j_index = i_index * 2
                            Assembleur_COO_F.append(index+j_index, 0, listFlocal[i_index,0])
                            Assembleur_COO_F.append(index+6+j_index, 0, listFlocal2[i_index,0])

            
            elif F_ext['Type de charge'] == 'Nodale':
                nearest_node = self._nearest_node(int(F_ext['Position']), bar_id)
                force = F_ext['Charge'] * 10
                
                if F_ext['Axe'] == 'Z' or F_ext['Axe'] == 'Z local' :                    
                    if nearest_node['position'] == "exactly at":
                        listFlocal[1,0] = force
                    
                    # Cas dans lequel la force est aprés le noeud le plus proche
                    elif nearest_node['position'] == "before":
                        
                        a = F_ext['Position'] - nearest_node['nearest node local x']
                        after_node_local_x = self.global_to_local(self.node_coor[nearest_node['index']+1], self.bar_info[bar_id])
                        b = after_node_local_x - F_ext['Position']
                        Mya = (force * a * b**2)/l**2
                        Myb = -(force * b * a**2)/l**2
                        Rza = (force * b**2 * (3 * a + b))/l**3
                        Rzb = (force * a**2 * (3 * b + a))/l**3
                        
                        listFlocal[1,0] = Rza
                        listFlocal[2,0] = Mya
                        
                        Assembleur_COO_F.append(nearest_node['index matrice']+6+2, 0, Rzb)
                        Assembleur_COO_F.append(nearest_node['index matrice']+6+4, 0, Myb)
                        
                        #print("La force est situé après le noeud le plus proche à: ", a,b)
                    
                    # Cas dans lequel la force est avant le noeud le plus proche
                    else:
                        b = nearest_node['nearest node local x'] - F_ext['Position']
                        before_node_local_x = self.global_to_local(self.node_coor[nearest_node['index']-1], self.bar_info[bar_id])
                        a = F_ext['Position'] - before_node_local_x
                        Mya = (force * a * b**2)/l**2
                        Myb = -(force * b * a**2)/l**2
                        Rza = (force * b**2 * (3 * a + b))/l**3
                        Rzb = (force * a**2 * (3 * b + a))/l**3
                        
                        listFlocal[1,0] = Rzb
                        listFlocal[2,0] = Myb
                        
                        Assembleur_COO_F.append(nearest_node['index matrice']-6+2, 0, Rza)
                        Assembleur_COO_F.append(nearest_node['index matrice']-6+4, 0, Mya)
                        #print("La force est situé avant le noeud le plus proche à: ", a,b)
                    
                elif F_ext['Axe'] == 'X':
                    listFlocal[0,0] = F_ext['Charge'] * 10
                    
                if F_ext['Axe'] in ('Z local', 'X local'):
                    listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                    
                
                Assembleur_COO_F.append(nearest_node['index matrice'], 0, listFlocal[0,0])
                Assembleur_COO_F.append(nearest_node['index matrice']+2, 0, listFlocal[1,0])
                Assembleur_COO_F.append(nearest_node['index matrice']+4, 0, listFlocal[2,0])
                
        self.matriceF = coo_matrix(Assembleur_COO_F.data, shape=((len(self.element_list)+1)*6, 1)).tocsr()
        print("\n Matrice des forces extérieurs Fext : \n",self.matriceF, self.matriceF.shape,"\n")
        return self.matriceF


    def _condition_limite(self):
        self.matriceK_CL = self.matriceK.toarray()
        self.matriceF_CL = self.matriceF.toarray()
     
        j = 0
        for i in range(len(self.matriceU)):
            if not self.matriceU[i]:
                self.matriceK_CL = np.delete(self.matriceK_CL,j,0)
                self.matriceK_CL = np.delete(self.matriceK_CL,j,1)
                self.matriceF_CL = np.delete(self.matriceF_CL,j,0)
            else:
                j += 1

        
        # plt.figure()
        # plt.pcolormesh(self.matriceK_CL)
        # plt.gca().invert_yaxis()
        # plt.colorbar()
        # plt.show()
        self.matriceK_CL = csr_matrix(self.matriceK_CL)
        self.matriceF_CL = csr_matrix(self.matriceF_CL)
        
        # print("\n Matrice des forces extérieure Fext au condition limite: \n", self.matriceF_CL,"\n")
        # print("\n Matrice de rigidité K au condition limite: \n", self.matriceK_CL,"\n")
        


    def _equa_deplacement(self):
        ui = spsolve(self.matriceK_CL,self.matriceF_CL)
        j = 0
        for i in range(len(self.matriceU)):
            if self.matriceU[i]:
                self.matriceU[i] = ui[j]
                j += 1
        # print("\n Solution Déplacement : \n", self.matriceU, "\n")
        return self.matriceU


    def _equa_reaction(self):
        self.ri = np.dot(self.matriceK.toarray(), self.matriceU) - self.matriceF.toarray()
        # print("\n Solution Réaction : \n", self.ri, "\n")
        return self.ri


    def reaction(self):
        """ Reaction en kN ou kN.m le long de l'élément"""
        rx = []
        ry = []
        rz = []
        rmx = []
        rmy = []
        rmz = []
        
        for i in range(len(self.element_list)):
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6

            r1 = self.ri[n1-6:n1,0]
            r2 = self.ri[n2-6:n2,0]
            
            if n1 == 6:
                rx.append(r1[0]/10**3)
                rz.append(r1[2]/10**3)
                rmy.append(r1[4]/10**6)

            rx.append(r2[0]/10**3)
            rz.append(r2[2]/10**3)
            rmy.append(r2[4]/10**6)

        self.react_coor = [rx, rz, rmy]
        return self.react_coor

    
    def reaction_max(self):
        """ Retourne les efforts min et max d'une liste d'effort interne """

        rx_max = max(self.react_coor[0])
        rx_max_index = self.react_coor[0].index(rx_max)
        rx_max_coor = [self.node_coor[rx_max_index,0], rx_max]

        rx_min = min(self.react_coor[0])
        rx_min_index = self.react_coor[0].index(rx_min)
        rx_min_coor = [self.node_coor[rx_min_index,0], rx_min]

        rz_max = max(self.react_coor[1])
        rz_max_index = self.react_coor[1].index(rz_max)
        rz_max_coor = [self.node_coor[rz_max_index,0], rz_max]

        rz_min = min(self.react_coor[1])
        rz_min_index = self.react_coor[1].index(rz_min)
        rz_min_coor = [self.node_coor[rz_min_index,0], rz_min]

        rmy_max = max(self.react_coor[2])
        rmy_max_index = self.react_coor[2].index(rmy_max)
        rmy_max_coor = [self.node_coor[rmy_max_index,0], rmy_max]

        rmy_min = min(self.react_coor[2])
        rmy_min_index = self.react_coor[2].index(rmy_min)
        rmy_min_coor = [self.node_coor[rmy_min_index,0], rmy_min]

        dictMinMax = {"Rx_max": rx_max_coor, "Rx_min": rx_min_coor, 
                      "Rz_max": rz_max_coor, "Rz_min": rz_min_coor, 
                      "RMy_max": rmy_max_coor, "RMy_min": rmy_min_coor}
        return dictMinMax

    
    def deplacement(self):
        uX = []
        uY = []
        uZ = []
        tetaY = []

        for i in range(len(self.element_list)):
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6

            u1 = self.matriceU[n1-6:n1,0]
            u2 = self.matriceU[n2-6:n2,0]

            if n1 == 6:

                uX.append(u1[0])
                uX.append(u2[0])
                uY.append(u1[1])
                uY.append(u2[1])
                uZ.append(u1[2])
                uZ.append(u2[2])
                tetaY.append(u1[4])
                tetaY.append(u2[4])
            else:
                uX.append(u2[0])
                uY.append(u2[1])
                uZ.append(u2[2])
                tetaY.append(u2[4])

        self.u_coor = {"uX": uX, "uY": uY, "uZ": uZ, "teta_Y": tetaY}
        return self.u_coor
        

    def deplacement_max(self):
        depl_min_max = {}
        for depl_type in self.u_coor.keys():

            depl_max = max(self.u_coor[depl_type])
            depl_max_index = self.u_coor[0].index(depl_max)
            depl_max_coor = [self.node_coor[depl_max_index,0], depl_max]

            depl_min = min(self.u_coor[depl_type])
            depl_min_index = self.u_coor[0].index(depl_min)
            depl_min_coor = [self.node_coor[depl_min_index,0], depl_min]

            depl_min_max[depl_type + "_max"] = depl_max_coor
            depl_min_max[depl_type + "_min"] = depl_min_coor
        return depl_min_max


    def effort_interne(self):

        ni, vyi, vzi, myi, mzi = [], [], [], [], []
        
        for i in range(len(self.element_list)):
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6

            # print("noeud 1: ", n1)
            # print("noeud 2: ", n2, "\n")
            
            u_global = self.matriceU[n1-6:n2,0]

            matriceF = self.matriceF.toarray()
            f_0 = matriceF[n1-6:n2,0]

            # effortinterneN1 = np.dot(self.k_local["eleALL"]["k11_local"],u1) + np.dot(self.k_local["eleALL"]["k12_local"], u2) - f1
            # effortinterneN2 = -(np.dot(self.k_local["eleALL"]["k21_local"],u1) + np.dot(self.k_local["eleALL"]["k22_local"], u2) - f2)
            T_matrix = self.k_local["ele"+str(i)]["T_global"]
            internal_forces = np.dot(np.dot(self.k_local["ele"+str(i)]["k_local"], T_matrix.toarray()), u_global) - np.dot(T_matrix.toarray(), f_0)
            # print("effort interne: ", internal_forces)

            if n1 == 6:
                ni.append(internal_forces[0]/10**3)
                vyi.append(internal_forces[1]/10**3)
                vzi.append(internal_forces[2]/10**3)
                myi.append(internal_forces[4]/10**6)
                mzi.append(internal_forces[5]/10**6)
            
            # else:
            #     ni[-1] = ni[-1] + internal_forces[0]/10**3
            #     vi[-1] = vi[-1] + internal_forces[2]/10**3
            #     myi[-1] = myi[-1] + internal_forces[4]/10**6
            
            ni.append(-internal_forces[6]/10**3)
            vyi.append(-internal_forces[7]/10**3)
            vzi.append(-internal_forces[8]/10**3)
            myi.append(-internal_forces[10]/10**6)
            mzi.append(-internal_forces[11]/10**6)
                
        self.ei_coor = {"Nx": ni, "Vy": vyi, "Vz": vzi, "My": myi, "Mz": mzi}
        return self.ei_coor
    
    

    def effort_interne_position(self, position: int):
        position_index = self._nearest_node(position)
        nx = self.ei_coor["Nx"][position_index[0]]
        vy = self.ei_coor["Vy"][position_index[0]]
        vz = self.ei_coor["Vz"][position_index[0]]
        my = self.ei_coor["My"][position_index[0]]
        mz = self.ei_coor["Mz"][position_index[0]]
        print(f"la position la plus proche est {position_index[1]} mm", {"Nx": nx, "Vy": vy, "Vz": vz, "My": my, "Mz": mz})
        return {"Nx": nx, "Vy": vy, "Vz": vz, "My": my, "Mz": mz}
    

    
    def effort_interne_max(self):
        """ Retourne les efforts min et max d'une liste d'effort interne """
        ei_min_max = {}
        for force_type in self.ei_coor.keys():
            force_max = max(self.ei_coor[force_type])
            force_max_index = self.ei_coor[force_type].index(force_max)
            force_max_coor = {"Position": self.node_coor[force_max_index,0], "Effort": force_max}
            
            force_min = min(self.ei_coor[force_type])
            force_min_index = self.ei_coor[force_type].index(force_min)
            force_min_coor = {"Position": self.node_coor[force_min_index,0], "Effort": force_min}

            ei_min_max[force_type + "_max"] = force_max_coor
            ei_min_max[force_type + "_min"] = force_min_coor
        return ei_min_max


    def _graphique_reaction_Z(self, name_combi: str):
        """ Retourne le diagramme des réactions d'appuis Z """

        self.name = name_combi
        self.title = 'Réaction aux appuis z : '
        self.color = "darkcyan"
        self.y_label = "Effort (kN)"
        self.unit = " kN"
        self.x_values = self.node_coor
        self.y_values = self.react_coor[1]
        
        dictRMinMax = self.reaction_max()
        self.max_XY_values = (dictRMinMax["Rz_max"][0], dictRMinMax["Rz_max"][1])
        self.min_XY_values = (dictRMinMax["Rz_min"][0], dictRMinMax["Rz_min"][1])
        self.excent_X = 100
        self.excent_Y = -0.5
        self.fill_between = False
        
        self.save_path = os.path.join(self.SAVE_PATH, "reactionZ.png")
        self._base_graphique()
        
    
    def show_graphique_reaction_Z(self, name_combi: str):
        self._graphique_reaction_Z(name_combi)
        plt.show()


    def _graphique_reaction_X(self, name_combi: str):
        """ Retourne le diagramme des réactions d'appuis X """
        
        self.name = name_combi
        self.title = 'Réaction aux appuis x : '
        self.color = "darkcyan"
        self.y_label = "Effort (kN)"
        self.unit = " kN"
        self.x_values = self.node_coor
        self.y_values = self.react_coor[0]
        
        dictRMinMax = self.reaction_max()
        self.max_XY_values = (dictRMinMax["Rx_max"][0], dictRMinMax["Rx_max"][1])
        self.min_XY_values = (dictRMinMax["Rx_min"][0], dictRMinMax["Rx_min"][1])
        self.excent_X = 100
        self.excent_Y = -0.5
        self.fill_between = False
        
        self.save_path = os.path.join(self.SAVE_PATH, "reactionX.png")
        self._base_graphique()


    def show_graphique_reaction_X(self, name_combi):
        self._graphique_reaction_X(name_combi)
        plt.show()


    def _graphique_fleche(self, name_combi:str):

        self.name = name_combi
        self.title = 'Flèche z : '
        self.color = "g"
        self.y_label = "Déplacement \n(mm)"
        self.unit = " mm"
        self.x_values = self.node_coor
        self.y_values = self.u_coor["uz"]
        
        dictUMinMax = self.deplacement_max()
        self.max_XY_values = (dictUMinMax["Uz_max"][0], dictUMinMax["Uz_max"][1])
        self.min_XY_values = (dictUMinMax["Uz_min"][0], dictUMinMax["Uz_min"][1])
        self.excent_X = 1
        self.excent_Y = 0.05
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, "fleche.png")
        self._base_graphique()
    
    
    def show_graphique_fleche(self, name_combi:str):
        self._graphique_fleche(name_combi)
        plt.show()


    def _graphique_Nx(self, name_combi:str):
        """ Retourne le diagramme des efforts normaux """
        
        self.name = name_combi
        self.title = 'Efforts normaux Nx : '
        self.color = "orange"
        self.x_values = self.node_coor
        self.y_values = self.ei_coor["Nx"]
        
        dictEiMinMax = self.effort_interne_max()
        self.max_XY_values = (dictEiMinMax["Nx_max"]["Position"], dictEiMinMax["Nx_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax["Nx_min"]["Position"], dictEiMinMax["Nx_min"]["Effort"])
        self.excent_X = 100
        self.excent_Y = 0.05
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, "Nx.png")
        self._base_graphique()
    

    def show_graphique_Nx(self, name_combi:str):
        self._graphique_Nx(name_combi)
        plt.show()


    def _graphique_V(self, name_combi:str, axe: str=("y", "z")):
        """ Retourne le diagramme des efforts tranchant suivant l'axe donnée y ou z local"""

        self.name = name_combi
        self.title = f'Efforts tranchants V{axe} : '
        self.color = 'b'
        self.x_values = self.node_coor
        self.y_values = self.ei_coor[f"V{axe}"]
        
        dictEiMinMax = self.effort_interne_max()
        self.max_XY_values = (dictEiMinMax[f"V{axe}_max"]["Position"], dictEiMinMax[f"V{axe}_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax[f"V{axe}_min"]["Position"], dictEiMinMax[f"V{axe}_min"]["Effort"])
        self.excent_X = 50
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, f"V{axe}.png")
        self._base_graphique()
    

    def show_graphique_V(self, name_combi:str, axe: str=("y", "z")):
        self._graphique_V(name_combi, axe)
        plt.show()

    
    def _graphique_M(self, name_combi:str, axe: str=("y", "z")):
        """ Retourne le diagramme des efforts de moment autour de M suivant l'axe y ou z local """

        self.name = name_combi
        self.title = f'Moments M{axe} : '
        self.color = 'r'
        self.y_label = "Effort (kN.m)"
        self.unit = " kN.m"
        self.x_values = self.node_coor
        self.y_values = self.ei_coor[f"M{axe}"]
        
        dictEiMinMax = self.effort_interne_max()
        self.max_XY_values = (dictEiMinMax[f"M{axe}_max"]["Position"], dictEiMinMax[f"M{axe}_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax[f"M{axe}_min"]["Position"], dictEiMinMax[f"M{axe}_min"]["Effort"])
        self.excent_X = 1
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, f"M{axe}.png")
        self._base_graphique()

    def show_graphique_M(self, name_combi, axe: str=("y", "z")):
        self._graphique_M(name_combi, axe)
        plt.show()


    def calcul_1D(self):
        """Calcul une poutre MEF 1D

        Args:
            list_loads (list): liste des charges combinées sur la poutre MEF.
        """
        start = time.time()
        loads = self.combinaison.get_all_loads()
        supports = self.get_supports()
        self._matrice_Fext_1D(loads)
        self._matrice_U_1D(supports)
        self._matrice_K_1D()

        self._condition_limite()
        self._equa_deplacement()
        self._equa_reaction()

        self.reaction()
        self.effort_interne()
        self.deplacement()
        
        end = time.time()
        elapsed = end - start
        print(f'Temps d\'exécution : {elapsed:.2}s')

    
    def graphique(self, name_combi:str):
        """Affiche le graphique des efforts internes, des réactions d'appuis et de la flèche,
            correspondant à la poutre MEF.

        Args:
            name_combi (str): défini le type de graphique à afficher, ELU ou ELS 
                              et écrit le nom de la combinaison dans le graphique.
        """
        if name_combi[0:3] == "ELU":
            self.show_graphique_reaction_Z(name_combi)
            self.show_graphique_reaction_X(name_combi)
            self.show_graphique_Nx(name_combi)
            self.show_graphique_V(name_combi, "z")
            self.show_graphique_M(name_combi, "y")
        else:
            self.show_graphique_fleche(name_combi)
        

if __name__ == '__main__':
    from EC0_Combinaison import Chargement
    # _list_loads = [[1, '', 'Permanente G', 'Linéique', -100, "0/6000", 'Z'],
    #              [2, '', "Neige normale Sn", 'Linéique', -200, "0/6000", 'Z']]
    _list_loads = [[1, '', 'Permanente G', -500, 2000, 'X'], [2, '', 'Permanente G', 500, 2000, 'X'],]
    chargement = Chargement(pays="Japon")
    chargement.create_load_by_list(_list_loads)
    c1 = Combinaison._from_parent_class(chargement, cat="Cat A : habitation", kdef=0.6)
    # print(c1.list_combination)
    rcombi = "ELU_STR G"
    print(c1.get_combi_list_load(rcombi))


    b = 60
    h = 100
    a = b*h
    iy = (b*h**3)/12
    iz = (h*b**3)/12

    beam_gen = Bar_generator()
    beam_gen.add_bar(0,0,2828.42,2828.42,a)
    beam_gen.add_bar(2828.42,2828.42,5656.84,0,a)
    
    beam_gen.add_relaxation(1, "end")
    # beam_gen.add_relaxation(2, "start")
    beam_gen.add_material_by_class(1, iy, iz, "C24")
    beam_gen.add_material_by_class(2, iy, iz, "C24")

    listdeplacement = [[1, "Rotule", 0, 0], [2, "Rotule", beam_gen.bar_info[2]["length"], 0], [2, "Simple Y", 0, 0]]
    beam_gen.create_supports_by_list(listdeplacement)
   
    
    mef = MEF._from_parent_class(beam_gen, combinaison=c1)
    
    mef.calcul_1D()
    mef.show_graph_loads_and_supports()
    # print(mef.reaction_max())
    # mef.graphique(rcombi)
    # mef.show_graphique_fleche(rcombi)
    # mef.show_graph_loads_and_supports()
