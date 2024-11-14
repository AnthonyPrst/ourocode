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

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# from eurocode.A0_Projet import Bar_generator
# from eurocode.EC0_Combinaison import Combinaison

from ourocode.eurocode.EC0_Combinaison import Combinaison
from ourocode.eurocode.A0_Projet import Bar_generator

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

    def get_local_coor_value_in_bar(self, list_value: list, bar_id: int):
            x_coords = []
            y_coords = []
            list_ele = [self.element_list[element] for element in self.bar_info[bar_id]["elements"]]

            for index_ele, ele in enumerate(list_ele):
                node = [int(ele[0])-1, int(ele[1])-1]
                x1, y1, z1 = self.node_coor[node[0]]
                x2, y2, z2 = self.node_coor[node[1]]
                v1 = (x2-x1, y2-y1, z2-z1)
                length = mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2 + abs(v1[2])**2)
                if not index_ele:
                    x_coords.append(0)
                    y_coords.append(list_value[node[0]])
                x_coords.append(x_coords[-1]+length)
                y_coords.append(list_value[node[1]])
            return x_coords, y_coords

    def get_local_internal_forces_coor_value_in_bar(self, list_value: list, internal_type: str, bar_id: int):
        x_coords = []
        y_coords = []
        list_ele = [self.element_list[element] for element in self.bar_info[bar_id]["elements"]]

        # for ele in list_ele:
        for index_ele in self.bar_info[bar_id]["elements"]:
        #     # On récupère l'index de l'élément dans la liste numpy 
        #     print("np xhere", np.where(self.element_list == ele))
        #     index_ele = np.where(self.element_list == ele)[0][0]
            ele = self.element_list[index_ele]

            node = [int(ele[0])-1, int(ele[1])-1]
            # x1, y1 = self.node_coor[node[0]]
            # x2, y2 = self.node_coor[node[1]]
            # v1 = (x2-x1, y2-y1)
            # length = mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2)
            for index_node, node in enumerate(node):
                x_local = self.global_to_local(self.node_coor[node], self.bar_info[bar_id])
                y_local = list_value[index_ele][internal_type][index_node]
                x_coords.append(x_local)
                y_coords.append(y_local)
        print(x_coords, y_coords)
        return x_coords, y_coords
    
    
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

    def remove_index(self, index: int):
        del self.data[0][index]
        del self.data[1][0][index]
        del self.data[1][1][index]

    def to_array(self):
        # Maximum Y and X coordinates
        xmax = max(self.data[1][0])
        ymax = max(self.data[1][1])
        # Target array
        target = np.zeros((xmax+1, ymax+1))
        target[self.data[1][1], self.data[1][0]] = self.data[0]
        print(target)
        return target

    
    
        
class FEM(Bar_generator, _Base_graph):
    SAVE_PATH = os.path.join(Combinaison.PATH_CATALOG, "data","screenshot")
    def __init__(self, combinaison: Combinaison, **kwargs):
        """Classe permettant de créer des poutres MEF et de les calculer. Cette classe est hérité de l'objet Combinaison du module EC0_Combinaison.py
        Args:
            A (float): Section en mm²
        """
        super(Bar_generator, self).__init__(**kwargs)
        _Base_graph.__init__(self)
        self.combinaison = combinaison
        self._fem_info = {}

    def get_lenght_of_element(self, element_index: int):
        """Retourne la longueur d'un élément en mm à partir de ces coordonnées globale

        Args:
            element_index (int): index de l'élément à récupérer
        """
        element = self.element_list[element_index]
        print("élement: ", element)
        node_1, node_2 = element
        x1, y1, z1 = self.node_coor[int(node_1)-1]
        x2, y2, z2 = self.node_coor[int(node_2)-1]
        v1 = (x2-x1, y2-y1, z2-z1)
        # print(x1, y1, x2, y2)
        return abs(mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2 + abs(v1[2])**2))
    
    
    def _transformation_matrix(self, bar: dict):
        """Matrice de transformation pour un élément barre 3D dans un plan perpendiculaire à un axe.
        Si la barre est réellement dans un espace 3D il faut alors calculer les produits scalaire de chaque axe dans une matrice lambda.
        Voir la video youtube https://www.youtube.com/watch?v=J2MwWFcY3NI&list=PLmw2x4fCMxg6Yw-7FqwBNSobXZR6fQ3Sm&index=16
        Args:
            angle (float): angle en degré de la barre par rapport à X
        """
        Assembleur_COO_T = _Triplets()

        # Extraire les coordonnées des extrémités de la barre
        x1, y1, z1 = bar['x1'], bar['y1'], bar['z1']
        x2, y2, z2 = bar['x2'], bar['y2'], bar['z2']
        
        # Vecteur directionnel de la barre
        L = bar["length"]

        x_local_vector_unit = np.array([(x2 - x1) / L, (y2 - y1) / L, (z2 - z1) / L])
        # Vecteur global Z
        Z_global_vector_unit = np.array([0,0,1])

        # Calculer le produit vectoriel pour obtenir le vecteur perpendiculaire y à Z et x
        v = np.cross(Z_global_vector_unit, x_local_vector_unit)  # Axe transverse
        if np.linalg.norm(v) == 0:
                # Si L est parallèle à l'axe de référence, on choisit un autre vecteur de référence (par exemple, l'axe Y)
                y_local_vector_unit = np.array([0, 1, 0])
        else:
            Lv = np.linalg.norm(v) # distance Euclidienne
            y_local_vector_unit = v / Lv
        
        # Calculer le produit vectoriel pour obtenir le vecteur perpendiculaire z à x et y
        w = np.cross(x_local_vector_unit, y_local_vector_unit)  # Axe transverse
        Lw = np.linalg.norm(w) # distance Euclidienne
        z_local_vector_unit = w / Lw
        # print(x_local_vector_unit, y_local_vector_unit, z_local_vector_unit)
        lambda_matrix = np.array([
                                x_local_vector_unit,
                                y_local_vector_unit,
                                z_local_vector_unit
                                ])
        index = 0
        for _ in range(0,4):
            for i in range(0,3):
                for j in range(0,3):
                    Assembleur_COO_T.append(index+i, index+j, lambda_matrix[i,j])
            index += 3

        transformation_matrix = (coo_matrix(Assembleur_COO_T.data, shape=(12, 12))).tocsr()
        # print("Matrice de transformation locale/globale: \n",transformation_matrix.toarray())
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


    # def _apply_relaxation_to_F(self, element: list, relaxation):
    #     dof_map = {
    #         "u": 0,
    #         "v": 1,
    #         "w": 2,
    #         "teta_x": 3,
    #         "teta_y": 4,
    #         "teta_z": 5
    #     }

    #     if relaxation["position"] == "start":
    #         index_F = (int(element[0])-1) * 6
    #     elif relaxation["position"] == "end":
    #         index_F = (int(element[1])-1) * 6
    #     else:
    #         raise ValueError("Position must be 'start' or 'end'")
    #     # print(element)
    #     # print(index_F)
    #     for key, rel in relaxation.items():
    #         if key != "position" and rel:
    #             row = index_F+dof_map[key]
    #             indexes_to_remove = [i for i,x in enumerate(self.matriceF.data[1][0]) if x==row]
    #             for i, index in enumerate(indexes_to_remove):
    #                 index -= i
    #                 self.matriceF.remove_index(index)


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
                # print(f"J'ai retrouvé la poutre {key} qui est associé à l'élément {index}, son angle est de {angle} degrés.")
                break
        k = np.array([
                    [(EA)/l, 0, 0, 0, 0, 0, -(EA)/l, 0, 0, 0, 0, 0],
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
                    [0, (6*EIz)/l**2, 0, 0, 0, (2*EIz)/l, 0, -(6*EIz)/l**2, 0, 0, 0, (4*EIz)/l]
                    ])
        
        if bar.get("relaxation"):
            # change_matrix_k = False
            for relaxation in bar["relaxation"]:
                if relaxation["position"] == "start" and index == bar["elements"][0]:
                    k = self._apply_relaxation_to_k(k, relaxation, relaxation["position"])
                    break
                elif relaxation["position"] == "end" and index == bar["elements"][-1]:
                    k = self._apply_relaxation_to_k(k, relaxation, relaxation["position"])
                    break

        T_matrix = self._transformation_matrix(bar)

        self._fem_info["ele"+str(index)] = {"k_local": k, "T_global": T_matrix}
        # print(index, "matrice k local: \n", k, "\n")
        return np.dot(np.dot(T_matrix.T.toarray(), k), T_matrix.toarray())
    


    def _agregation_matrix_K(self, assembleur_COO_K: list, k: np.array, n1: int, n2: int):
        for i in range(-6,0):
            for j in range(-6,0):
                # print(n1+i, n1+j, n2+i, n2+j, f'i et j: {[i, j]}', k[i+6,j+6], k[i+12,j+12 ], k[i+6,j+12], k[i+12,j+6])
                assembleur_COO_K.append(n1+i, n1+j, k[i+6,j+6]) # ajout de bruit pour éviter les singularité (+0.001)
                assembleur_COO_K.append(n2+i, n2+j, k[i+12,j+12]) # ajout de bruit pour éviter les singularité (+0.001)
                if k[i+6,j+12]:
                    assembleur_COO_K.append(n1+i, n2+j, k[i+6,j+12])
                if k[i+12,j+6]:
                    assembleur_COO_K.append(n2+i, n1+j, k[i+12,j+6]) # ajout de bruit pour éviter les singularité (+0.001)
                    
    
    def _matrice_K_1D(self):
        """Méthodequi créer la matrice de rigidité K assemblée
        Returns:
            numpy.array: la matrice de rigidité assemblée
        """
        Assembleur_COO_K = _Triplets()

        for i in range(len(self.element_list)):
            # matrice de rigidité K
            k = self._init_matrix_K(i)
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6
            
            self._agregation_matrix_K(Assembleur_COO_K , k, n1, n2)

        self.matriceF_concat = coo_matrix(self.matriceF.data, shape=((len(self.element_list)+1)*6, 1)).tocsr()
        self.matriceK = (coo_matrix(Assembleur_COO_K.data, shape=((len(self.element_list)+1)*6, (len(self.element_list)+1)*6))).tocsr()
        np.savetxt("Matrice de rigidité K.csv", self.matriceK.toarray(), delimiter=";", fmt='%f')
        # print("\n Matrice de rigidité K  : \n", self.matriceK.toarray(), "\n")
        print("\n Matrice des forces extérieur  : \n", self.matriceF_concat.toarray(), "\n")
        return self.matriceK
    

    def local_to_global(self, x_local, bar):
        """
        Convertit la coordonnée locale x sur une barre en coordonnées globales (X, Y, Z).

        Args:
            x_local (float): La position locale le long de la barre.
            bar (dict): Dictionnaire contenant les informations sur la barre, notamment les coordonnées globales des extrémités (x1, y1, z1 et x2, y2, z2) et les angles.

        Returns:
            numpy.array: Coordonnées globales (X, Y, Z) correspondant à la position locale.
        """
        # Coordonnées des extrémités de la barre
        x1, y1, z1 = bar['x1'], bar['y1'], bar['z1']
        x2, y2, z2 = bar['x2'], bar['y2'], bar['z2']
        
        # Vecteur directionnel de la barre en 3D
        length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        direction = np.array([x2 - x1, y2 - y1, z2 - z1]) / length  # Vecteur unitaire de direction

        # Coordonnées locales projetées sur la barre (translation le long du vecteur directionnel)
        X = x1 + x_local * direction[0]
        Y = y1 + x_local * direction[1]
        Z = z1 + x_local * direction[2]
        return np.array([X, Y, Z])
    

    def global_to_local(self, global_coor, bar):
        """Convert global coordinates (X, Y) to local x coordinate on a bar."""
        x1, y1, z1 = bar['x1'], bar['y1'], bar['z1']
        x2, y2, z2 = bar['x2'], bar['y2'], bar['z2']
        
        bar_vector = np.array([x2-x1, y2-y1, z2-z1])
        point_vector = global_coor - np.array([x1, y1, z1])
        
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
        # print("local x: ",nearest_node_local_x)
        
        # Compare local x coordinates
        if np.isclose(nearest_node_local_x, x_local, atol=1e-6):
            position = "exactly at"
        elif nearest_node_local_x < x_local:
            position = "before"
        else:
            position = "after"
        
        index_I_matrice = ((index+1) * 6) - 6

        # Recherche de l'élément correspondant dans la poutre
        # On parcourt les éléments de la barre pour vérifier lequel correspond au nœud le plus proche
        for element_id in bar['elements']:
            element_nodes = self.element_list[element_id]
            node1 = self.node_coor[int(element_nodes[0]) - 1]  # Premier nœud de l'élément
            node2 = self.node_coor[int(element_nodes[1]) - 1]  # Deuxième nœud de l'élément
            
            # Vérification si le nœud appartient à l'élément (entre node1 et node2)
            if np.allclose(nearest_node_coor, node1, atol=1e-6) or np.allclose(nearest_node_coor, node2, atol=1e-6):
                # Si le nœud est l'un des nœuds de l'élément, on retourne cet élément
                return {
                    "index": index,
                    "nearest node": nearest_node_coor,
                    "nearest node local x": nearest_node_local_x,
                    "index matrice": index_I_matrice,
                    "element": element_id,  # ID correct de l'élément
                    "position": position
                }

        # Si aucun élément n'est trouvé (ce qui ne devrait pas arriver)
        raise ValueError(f"Aucun élément correspondant trouvé pour le nœud le plus proche à {nearest_node_coor}")


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
                self.matriceU[nearest_node['index matrice'] + 1,0] = False
                self.matriceU[nearest_node['index matrice'] + 3,0] = False

            elif support["Type d'appui"] == 'Simple Y':
                self.matriceU[nearest_node['index matrice'] + 2,0] = False
                self.matriceU[nearest_node['index matrice'] + 3,0] = False

            elif support["Type d'appui"] == 'Simple X':
                self.matriceU[nearest_node['index matrice'],0] = False
                self.matriceU[nearest_node['index matrice'] + 3,0] = False

        # print("\n Matrice de déplacement U : \n",self.matriceU, "\n")
        return self.matriceU
    
    
    def _counter_unique_in_array(self, array: np.ndarray):
        unique, counts = np.unique(array, return_counts=True)
        return dict(zip(unique, counts))
    

    def _create_local_force_arrays(self):
        """Créer des matrices de forces locales initialisées à zéro."""
        return np.zeros((6, 1)), np.zeros((6, 1))
    

    def _get_start_end_indices(self, F_ext, bar_id):
        """Retourne les indices de départ et de fin pour unne force linéique donné."""
        start_pos, end_pos = F_ext['Position'].split("/")
        start_index = self._nearest_node(int(start_pos), bar_id)
        end_index = self._nearest_node(int(end_pos), bar_id)
        return start_index, end_index, [[start_index, int(start_pos)], [end_index, int(end_pos)]]
    

    def _calculate_forces(self, F_ext, l):
        """Calcul des forces locales en fonction de l'axe de la charge."""
        force = F_ext['Charge'] * 10**-2 * l
        listFlocal, listFlocal2 = self._create_local_force_arrays()

        if F_ext['Axe'] in ('Z', 'Z local'):
            listFlocal[1, 0] = listFlocal2[1, 0] = force / 2  # Rza et Rzb
        elif F_ext['Axe'] in ('X', 'X local'):
            listFlocal[0, 0] = listFlocal2[0, 0] = force / 2  # Rxa et Rxb
            
        listFlocal[5, 0] = force * l / 12  # Mya
        listFlocal2[5, 0] = -listFlocal[5, 0]  # Myb

        return listFlocal, listFlocal2
    
    def _adjust_forces_for_nearest_node(self, nearest_node_index, F_ext, l, index, bar_id, listFlocal, listFlocal2):
        """Ajuste les forces pour les nœuds proches en fonction de la position de la charge."""
        force = F_ext['Charge'] * 10**-2 * l
        if index != nearest_node_index[1][0]['index matrice']:
            ind = nearest_node_index[0]
            
            node_local_x = ind[0]['nearest node local x']
            # Cas dans le quelle la force est avant le noeud le plus proche
            # print(ind[1], node_local_x, type(ind[1]), type(node_local_x))
            if ind[1] < node_local_x:
                a = node_local_x - ind[1]
                before_node_local_x = self.global_to_local(self.node_coor[ind[0]['index']-1], self.bar_info[bar_id])
                b = ind[1] - before_node_local_x
                Mya = ((force * b**2)/l) * (4*(b/l) + 3*(b/l)**2)
                Myb = -((force * b**2)/12) * (6-8*(b/l) + 3*(b/l)**2)
                force = force * b
                b = b/2 
                a = a + b 
                Rza = (force * b**2 * (3 * a + b))/l**3
                Rzb = (force * a**2 * (3 * b + a))/l**3
                
                if F_ext['Axe'] in ('Z', 'Z local'):
                    listFlocal[1,0] = Rza
                    listFlocal2[1,0] = Rzb
                elif F_ext['Axe'] in ('X', 'X local'):
                    listFlocal[0,0] = Rza
                    listFlocal2[0,0] = Rzb

                listFlocal[5,0] = Mya
                listFlocal2[5,0] = Myb
        
        else:
            ind = nearest_node_index[1]
            node_local_x = ind[0]['nearest node local x']
            # Formule tirée de l'aide mémoire page 113 cas 3 et cas 1 pour RA et RB
            # Cas dans lequel la force est aprés le noeud le plus proche
            if ind[1] > node_local_x:
                a = ind[1] - node_local_x
                after_node_local_x = self.global_to_local(self.node_coor[ind[0]['index']+1], self.bar_info[bar_id])
                b = after_node_local_x - ind[1]
                Mya = ((force * a**2)/12) * (6-8*(a/l) + 3*(a/l)**2)
                Myb = -((force * a**2)/l) * (4*(a/l) + 3*(a/l)**2)
                force = force * b
                a = a/2 
                b = a + b 
                Rza = (force * b**2 * (3 * a + b))/l**3
                Rzb = (force * a**2 * (3 * b + a))/l**3
                
                if F_ext['Axe'] in ('Z', 'Z local'):
                    listFlocal[1,0] = Rza
                    listFlocal2[1,0] = Rzb
                elif F_ext['Axe'] in ('X', 'X local'):
                    listFlocal[0,0] = Rza
                    listFlocal2[0,0] = Rzb

                listFlocal[5,0] = Mya
                listFlocal2[5,0] = Myb


    def _append_forces_to_assembler(self, index1, index2, listFlocal, listFlocal2, Assembleur_COO_F):
        """Ajoute les forces dans l'assembleur COO."""
        print(index1, index2)
        for i in range(len(listFlocal)):
            Assembleur_COO_F.append(index1 + i, 0, listFlocal[i, 0])
            Assembleur_COO_F.append(index2 + i, 0, listFlocal2[i, 0])


    def _store_forces_by_element(self, counter, element_index_travel, index, bar_id, listFlocal, listFlocal2):
        """Stocke les forces calculées dans un dictionnaire par élément."""
        if counter:
            current_index = int(index/6)
            x_local_of_current_index = self.global_to_local(self.node_coor[current_index], self.bar_info[bar_id])
            print(self._nearest_node(x_local_of_current_index, bar_id))
            current_element = self._nearest_node(x_local_of_current_index, bar_id)['element']

            if current_element not in element_index_travel:
                if not self.fext_by_element.get(current_element):
                    self.fext_by_element[current_element] = {
                            'forces_local_start': np.copy(listFlocal),
                            'forces_local_end': np.copy(listFlocal2)
                    }
                else:
                    self.fext_by_element[current_element]['forces_local_start'] += np.copy(listFlocal)
                    self.fext_by_element[current_element]['forces_local_end'] += np.copy(listFlocal2)
                element_index_travel.append(current_element)


    def _get_index_between_start_and_end_of_bar(self, start_index, end_index, bar_id):
        """
        Récupère les index de la matrice de force entre start_index et end_index en tenant compte des éléments connectés.
        """
        indices = []
        
        # Obtenir la liste des éléments associés à la barre en question
        elements_in_bar = self.bar_info[bar_id]["elements"]

        # Parcourir les éléments pour ajouter les index entre start et end
        for element in elements_in_bar:
            ele = self.element_list[element]
            node1 = ((ele[0]) * 6) - 6
            node2 = ((ele[1]) * 6) - 6
            node = [node1, node2]
            print(node)
            while node:
                node_pop = node.pop(0)
                # if start_index['index matrice'] <= node_pop <= end_index['index matrice'] or start_index['index matrice'] >= node_pop >= end_index['index matrice']:
                if node_pop not in indices:
                    indices.append(node_pop)
        print(indices)
        return indices


    def _matrice_Fext_1D(self, dict_F_ext: dict):
        """Transforme un dictionnaire de charges extérieures en une matrice pour le calcul MEF."""
        print(dict_F_ext)
        Assembleur_COO_F = _Triplets()
        self.fext_by_element = {}
        
        for _ , F_ext in dict_F_ext.items():
            bar_id = F_ext['N° barre']
            print("\n", bar_id, self.bar_info[bar_id])

            l = self.get_lenght_of_element(self.bar_info[bar_id]["elements"][0])

            # On récupère la matrice de transformation des efforts locaux à globaux
            T_matrix = self._transformation_matrix(self.bar_info[bar_id]).toarray()
            T_matrix = T_matrix[0:6,0:6]

            # Créer les listes de forces locales
            listFlocal, listFlocal2 = self._create_local_force_arrays()
            
            if F_ext['Type de charge'] == 'Linéique':
                element_index_travel = []
                start_index, end_index, nearest_node_index = self._get_start_end_indices(F_ext, bar_id)
                print(nearest_node_index)

                sens_inverse = start_index['index matrice'] > end_index['index matrice']
                range_index = range(start_index['index matrice'], end_index['index matrice'] + (6 if not sens_inverse else -6), (6 if not sens_inverse else -6))

                indexs_in_bar = self._get_index_between_start_and_end_of_bar(start_index, end_index, bar_id)
                for counter, index in enumerate(indexs_in_bar):
                    # Calcul des forces et moments
                    listFlocal, listFlocal2 = self._calculate_forces(F_ext, l)

                    if index in (start_index['index matrice'], end_index['index matrice']):
                        self._adjust_forces_for_nearest_node(nearest_node_index, F_ext, l, index, bar_id, listFlocal, listFlocal2)
                    
                    # Appliquer la rotation si nécessaire
                    if F_ext['Axe'] in ('Z local', 'X local'):
                        listFlocal= np.dot(T_matrix, listFlocal)
                        listFlocal2= np.dot(T_matrix, listFlocal2)

                    if sens_inverse:
                        temp = np.copy(listFlocal2)
                        listFlocal2 = listFlocal
                        listFlocal = temp

                    # Ajouter les forces dans l'assembleur COO
                    if index != max(end_index['index matrice'], start_index['index matrice']):
                        if not sens_inverse:
                            self._append_forces_to_assembler(index, indexs_in_bar[counter+1], listFlocal, listFlocal2, Assembleur_COO_F)
                        else:
                            self._append_forces_to_assembler(index, indexs_in_bar[counter-1], listFlocal, listFlocal2, Assembleur_COO_F)

                    # Stocker les forces par élément
                    # print("counter: ", counter, 
                    #       "element_index_travel: ", element_index_travel, 
                    #       "index: ", index, "\n", 
                    #       listFlocal, listFlocal2, '\n')
                    self._store_forces_by_element(counter, element_index_travel, index, bar_id, listFlocal, listFlocal2)

            
            elif F_ext['Type de charge'] == 'Nodale':
                node_pos = int(F_ext['Position'])  # Position du nœud pour une force nodale
                nearest_node = self._nearest_node(node_pos, bar_id)  # Trouver le nœud le plus proche
                force = F_ext['Charge'] * 10 # On passe la charge nodale de daN à N
                                   
                if nearest_node['position'] == "exactly at":
                    if F_ext['Axe'] in ('Z', 'Z local'):
                        listFlocal[1,0] = force
                    elif F_ext['Axe'] in ('X', 'X local'):
                        listFlocal[0,0] = force
                
                # Cas dans lequel la force est aprés le noeud le plus proche
                elif nearest_node['position'] == "before":
                    a = F_ext['Position'] - nearest_node['nearest node local x']
                    after_node_local_x = self.global_to_local(self.node_coor[nearest_node['index']+1], self.bar_info[bar_id])
                    b = after_node_local_x - F_ext['Position']
                    Mya = (force * a * b**2)/l**2
                    Myb = -(force * b * a**2)/l**2
                    Rza = (force * b**2 * (3 * a + b))/l**3
                    Rzb = (force * a**2 * (3 * b + a))/l**3
                    
                    if F_ext['Axe'] in ('Z', 'Z local'):
                        listFlocal[1,0] = Rza
                        Assembleur_COO_F.append(nearest_node['index matrice']+6+2, 0, Rzb)
                    elif F_ext['Axe'] in ('X', 'X local'):
                        listFlocal[0,0] = Rza
                        Assembleur_COO_F.append(nearest_node['index matrice']+6, 0, Rzb)

                    listFlocal[5,0] = Mya
                    Assembleur_COO_F.append(nearest_node['index matrice']+6+6, 0, Myb)
                    
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
                    
                    if F_ext['Axe'] in ('Z', 'Z local'):
                        listFlocal[1,0] = Rzb
                        Assembleur_COO_F.append(nearest_node['index matrice']-6+2, 0, Rza)
                    elif F_ext['Axe'] in ('X', 'X local'):
                        listFlocal[0,0] = Rzb
                        Assembleur_COO_F.append(nearest_node['index matrice']-6, 0, Rza)
                    listFlocal[5,0] = Myb
                    Assembleur_COO_F.append(nearest_node['index matrice']-6+6, 0, Mya)
                    #print("La force est situé avant le noeud le plus proche à: ", a,b)

                    
                if F_ext['Axe'] in ('Z local', 'X local'):
                    listFlocal= np.dot(T_matrix, listFlocal)
                
                for i_index in range(len(listFlocal)):
                    Assembleur_COO_F.append(nearest_node['index matrice']+i_index, 0, listFlocal[i_index,0])
                
                current_element = nearest_node["element"]
                if nearest_node['position'] == "exactly at":
                    applying_node = nearest_node["index"]+1
                    if applying_node == self.element_list[current_element][0]:
                        start_F_local = np.copy(listFlocal)
                        end_F_local = np.zeros((6,1))
                    else:
                        start_F_local = np.zeros((6,1))
                        end_F_local = np.copy(listFlocal)
                    
                    unique_element_counter = self._counter_unique_in_array(self.element_list)
                    if unique_element_counter[applying_node] > 1:
                        start_F_local = start_F_local / unique_element_counter[applying_node]
                        end_F_local = end_F_local / unique_element_counter[applying_node]

                        for i, element in enumerate(self.element_list):
                            if i != current_element and applying_node in element and i in self.bar_info[bar_id]["elements"]:
                                if not self.fext_by_element.get(i):
                                    self.fext_by_element[i] = {
                                        'forces_local_start': end_F_local,
                                        'forces_local_end': start_F_local
                                        }
                                else:
                                    self.fext_by_element[i]['forces_local_start'] += end_F_local
                                    self.fext_by_element[i]['forces_local_end'] += start_F_local
                                

                if not self.fext_by_element.get(current_element):
                    self.fext_by_element[current_element] = {
                        'forces_local_start': start_F_local,
                        'forces_local_end': end_F_local
                        }
                else:
                    self.fext_by_element[current_element]['forces_local_start'] += start_F_local
                    self.fext_by_element[current_element]['forces_local_end'] += end_F_local
        
        self.matriceF = Assembleur_COO_F
        print("\n Matrice des forces extérieurs Fext : \n",self.matriceF,"\n")
        # print("fext_by_element:\n", self.fext_by_element)
        return self.matriceF


    def _condition_limite(self):
        self.matriceK_CL = self.matriceK.toarray()
        self.matriceF_CL = self.matriceF_concat.toarray()
     
        j = 0
        for i in range(len(self.matriceU)):
            if not self.matriceU[i]:
                self.matriceK_CL = np.delete(self.matriceK_CL,j,0)
                self.matriceK_CL = np.delete(self.matriceK_CL,j,1)
                self.matriceF_CL = np.delete(self.matriceF_CL,j,0)
            else:
                j += 1

        
        plt.figure()
        plt.pcolormesh(self.matriceK_CL)
        plt.gca().invert_yaxis()
        plt.colorbar()
        plt.show()
        self.matriceK_CL = csr_matrix(self.matriceK_CL)
        self.matriceF_CL = csr_matrix(self.matriceF_CL)
        
        # print("\n Matrice des forces extérieure Fext au condition limite: \n", self.matriceF_CL.toarray(),"\n")
        # print("\n Matrice de rigidité K au condition limite: \n", self.matriceK_CL.toarray(),"\n")
        
    

    def _equa_deplacement(self):
        ui = np.linalg.lstsq(self.matriceK_CL.toarray(),self.matriceF_CL.toarray())[0]
        # ui = spsolve(self.matriceK_CL,self.matriceF_CL)
        j = 0
        for i in range(len(self.matriceU)):
            if self.matriceU[i]:
                self.matriceU[i] = ui[j]
                j += 1
        # print("\n Solution Déplacement : \n", self.matriceU, "\n")
        return self.matriceU


    def _equa_reaction(self):
        self._ri = np.dot(self.matriceK.toarray(), self.matriceU) - self.matriceF_concat.toarray()
        # print("\n Solution Réaction : \n", self._ri, "\n")
        return self._ri


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

            r1 = self._ri[n1-6:n1,0]
            r2 = self._ri[n2-6:n2,0]
            
            if n1 == 6:
                rx.append(r1[0]/10**3)
                rz.append(r1[1]/10**3)
                rmy.append(r1[5]/10**6)

            rx.append(r2[0]/10**3)
            rz.append(r2[1]/10**3)
            rmy.append(r2[5]/10**6)

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
        uX, uY, uZ, tetaX, tetaY, tetaZ = [], [], [], [], [], []
        ux, uy, uz, tetax, tetay, tetaz = [], [], [], [], [], []

        for i in range(len(self.element_list)):
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6

            u1 = self.matriceU[n1-6:n1,0]
            u2 = self.matriceU[n2-6:n2,0]

            u_global = np.concatenate([self.matriceU[n1-6:n1,0], self.matriceU[n2-6:n2,0]])
            # u_global = self.matriceU[n1-6:n2,0]
            T_matrix = self._fem_info["ele"+str(i)]["T_global"]
            u_local = np.dot(T_matrix.toarray(), u_global)

            if n1 == 6:

                uX.append(u_global[0])
                ux.append(u_local[0])

                uZ.append(u_global[1])
                uz.append(u_local[1])

                uY.append(u_global[2])
                uy.append(u_local[2])

                tetaX.append(u_global[3])
                tetax.append(u_local[3])

                tetaY.append(u_global[4])
                tetay.append(u_local[4])

                tetaZ.append(u_global[5])
                tetaz.append(u_local[5])

            uX.append(u_global[6])
            ux.append(u_local[6])

            uZ.append(u_global[7])
            uz.append(u_local[7])

            uY.append(u_global[8])
            uy.append(u_local[8])

            tetaX.append(u_global[9])
            tetax.append(u_local[9])

            tetaY.append(u_global[10])
            tetay.append(u_local[10])

            tetaZ.append(u_global[11])
            tetaz.append(u_local[11])

        self.u_coor = {"uX": uX, "uY": uY, "uZ": uZ, "teta_X": tetaX, "teta_Y": tetaY, "teta_Z": tetaZ,
                       "ux": ux, "uy": uy, "uz": uz, "teta_x": tetax, "teta_y": tetay, "teta_z": tetaz}
        return self.u_coor
        

    def deplacement_max(self):
        depl_min_max = {}
        for depl_type in self.u_coor.keys():

            depl_max = max(self.u_coor[depl_type])
            depl_max_index = self.u_coor[depl_type].index(depl_max)
            depl_max_coor = {"Position": self.node_coor[depl_max_index,0], "Deplacement": depl_max}

            depl_min = min(self.u_coor[depl_type])
            depl_min_index = self.u_coor[depl_type].index(depl_min)
            depl_min_coor = {"Position": self.node_coor[depl_min_index,0], "Deplacement": depl_min}

            depl_min_max[depl_type + "_max"] = depl_max_coor
            depl_min_max[depl_type + "_min"] = depl_min_coor
        return depl_min_max


    def _effort_interne(self):
        # on boucle sur les différentes barres
        for id_bar, bar in self.bar_info.items():
            bar["internals forces"] = {}
            # On récupère chaque élément composants la barre
            for i_ele in bar["elements"]:
                bar["internals forces"][i_ele]= {"Nx": [], "Vy": [], "Vz": [], "Mx": [], "My": [], "Mz": []}
                n1 = int(self.element_list[i_ele,0])*6
                n2 = int(self.element_list[i_ele,1])*6
                convention_sign = [(-1, 1), (1, -1), (-1, 1)]
                if n1 > n2:
                    convention_sign = [(-1, 1), (1, -1), (1, -1)]
                
                print("Barre: ", id_bar, " Elément : ", self.element_list[i_ele], f"{i_ele+1}/{len(self.element_list)}")
                print("noeud 1: ", n1)
                print("noeud 2: ", n2, "\n")
                
                u_global = np.concatenate([self.matriceU[n1-6:n1,0], self.matriceU[n2-6:n2,0]])
                try:
                    f_0 = np.concatenate([self.fext_by_element[i_ele]["forces_local_start"] ,self.fext_by_element[i_ele]["forces_local_end"]], axis=0).flatten()
                except KeyError:
                    f_0 = np.zeros((12,))
                T_matrix = self._fem_info["ele"+str(i_ele)]["T_global"]

                # print("déplacement:\n", u_global)
                # print("\n force ext.:\n", f_0)
                # print("\n Transform matrix:\n", T_matrix.toarray())
                # print("\n K local matrix:\n", self._fem_info["ele"+str(i_ele)]["k_local"])

                internal_forces = np.dot(self._fem_info["ele"+str(i_ele)]["k_local"], np.dot(T_matrix.toarray(), u_global)) - np.dot(T_matrix.toarray(), f_0)
                # print("déplacement local:\n", np.dot(T_matrix.toarray(), u_global))
                # print("F local total:\n", np.dot(self._fem_info["ele"+str(i_ele)]["k_local"], np.dot(T_matrix.toarray(), u_global)), "\n\n")
                # print("effort interne:\n", internal_forces, "\n\n")
                    
                bar["internals forces"][i_ele]["Nx"].append(convention_sign[0][0] * internal_forces[0]/10**3)
                bar["internals forces"][i_ele]["Nx"].append(convention_sign[0][1] * internal_forces[6]/10**3)

                bar["internals forces"][i_ele]["Vz"].append(convention_sign[1][0] * internal_forces[1]/10**3)
                bar["internals forces"][i_ele]["Vz"].append(convention_sign[1][1] * internal_forces[7]/10**3)

                bar["internals forces"][i_ele]["Vy"].append(convention_sign[1][0] * internal_forces[2]/10**3)
                bar["internals forces"][i_ele]["Vy"].append(convention_sign[1][1] * internal_forces[8]/10**3)

                bar["internals forces"][i_ele]["Mx"].append(convention_sign[2][0] * internal_forces[3]/10**6)
                bar["internals forces"][i_ele]["Mx"].append(convention_sign[2][1] * internal_forces[9]/10**6)

                bar["internals forces"][i_ele]["Mz"].append(convention_sign[2][0] * internal_forces[4]/10**6)
                bar["internals forces"][i_ele]["Mz"].append(convention_sign[2][1] * internal_forces[10]/10**6)

                bar["internals forces"][i_ele]["My"].append(convention_sign[2][0] * internal_forces[5]/10**6)
                bar["internals forces"][i_ele]["My"].append(convention_sign[2][1] * internal_forces[11]/10**6)

    
    def get_internal_forces(self, bar_id: int):
        """Retourne les efforts internes sous forme de dictionnaire avec comme clé l'id de l'élément de la barre.
        On récupère ensuite les efforts par typologie avec un tuple pour chaque élément composant la barre avec le noeud 1 puis le noeud 2.

        Returns:
            dict[dict[tuple]]: Dictionnaire des efforts internes en kN et kN.m par barre.
        """
        return self.bar_info[bar_id]["internals forces"]
    

    def effort_interne_position(self, position: int):
        position_index = self._nearest_node(position)
        nx = self._ei_coor["Nx"][position_index[0]]
        vy = self._ei_coor["Vy"][position_index[0]]
        vz = self._ei_coor["Vz"][position_index[0]]
        my = self._ei_coor["My"][position_index[0]]
        mz = self._ei_coor["Mz"][position_index[0]]
        print(f"la position la plus proche est {position_index[1]} mm", {"Nx": nx, "Vy": vy, "Vz": vz, "My": my, "Mz": mz})
        return {"Nx": nx, "Vy": vy, "Vz": vz, "My": my, "Mz": mz}
    

    def max_internal_forces(self, bar_id: int):
        """ Retourne les efforts min et max d'une barre spécifique"""
        ei_min_max = {}
        force_type_init = {}

        # On itère sur chaque élément de la barre et l'on récupère les efforts internes
        for i_ele, forces in  self.bar_info[bar_id]["internals forces"].items():
            node = [int(self.element_list[i_ele][0])-1, int(self.element_list[i_ele][1])-1]

            # On itère sur les efforts internes par type de force.
            for force_type, values in forces.items():
                # Initialisation du dictionnaire quand un nouveau type de force est présent
                if not force_type_init.get(force_type):
                    force_type_init[force_type] = {"force_max": 0, "force_min": 0}
                    ei_min_max[force_type + "_max"] = {"Position": 0, "Effort": 0}
                    ei_min_max[force_type + "_min"] = {"Position": 0, "Effort": 0}
                force_max = force_type_init[force_type]["force_max"]
                force_min = force_type_init[force_type]["force_min"]
                max_value = max(values)
                min_value = min(values)

                # On test si la valeur est supérieur ou inférieur à la valeur précédente
                if max_value >= force_max:
                    force_type_init[force_type]["force_max"] = max_value
                    force_max_index = values.index(max_value)
                    x_local = self.global_to_local(self.node_coor[node[force_max_index]], self.bar_info[bar_id])
                    force_max_coor = {"Position": x_local, "Effort": max_value}
                    ei_min_max[force_type + "_max"] = force_max_coor
                if min_value <= force_min:
                    force_type_init[force_type]["force_min"] = min_value
                    force_min_index = values.index(min_value)
                    x_local = self.global_to_local(self.node_coor[node[force_min_index]], self.bar_info[bar_id])
                    force_min_coor = {"Position": x_local, "Effort": min_value}
                    ei_min_max[force_type + "_min"] = force_min_coor
        return ei_min_max


    def show_graphique_reaction_Z(self, name_combi: str):
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
        plt.show()


    def show_graphique_reaction_X(self, name_combi: str):
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
        plt.show()


    def show_graphique_fleche(self, bar_id: int, name_combi:str, axe: str=("y", "z")):

        self.name = name_combi
        self.title = f'Barre {bar_id}: Flèche {axe}: '
        self.color = "g"
        self.y_label = "Déplacement \n(mm)"
        self.unit = " mm"
        self.x_values, self.y_values = self.get_local_coor_value_in_bar(self.u_coor[f"u{axe}"], bar_id)
        
        dictUMinMax = self.deplacement_max()
        self.max_XY_values = (dictUMinMax[f"u{axe}_max"]["Position"], dictUMinMax[f"u{axe}_max"]["Deplacement"])
        self.min_XY_values = (dictUMinMax[f"u{axe}_min"]["Position"], dictUMinMax[f"u{axe}_min"]["Deplacement"])
        self.excent_X = 1
        self.excent_Y = 0.05
        self.fill_between = True
        self.fill_between_X = self.x_values
        
        self.save_path = os.path.join(self.SAVE_PATH, f"fleche{axe}.png")
        self._base_graphique()
        plt.show()


    def show_graphique_Nx(self, bar_id: int, name_combi:str):
        """ Retourne le diagramme des efforts normaux """
        
        self.name = name_combi
        self.title = f'Barre {bar_id}: Efforts normaux Nx : '
        self.color = "orange"
        self.x_values, self.y_values = self.get_local_internal_forces_coor_value_in_bar(self.bar_info[bar_id]["internals forces"], "Nx", bar_id)
        
        dictEiMinMax = self.max_internal_forces(bar_id)
        self.max_XY_values = (dictEiMinMax["Nx_max"]["Position"], dictEiMinMax["Nx_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax["Nx_min"]["Position"], dictEiMinMax["Nx_min"]["Effort"])
        self.excent_X = 100
        self.excent_Y = 0.05
        self.fill_between = True
        self.fill_between_X = self.x_values
        
        self.save_path = os.path.join(self.SAVE_PATH, "Nx.png")
        self._base_graphique()
        plt.show()


    def show_graphique_V(self, bar_id: int, name_combi:str, axe: str=("y", "z")):
        """ Retourne le diagramme des efforts tranchant suivant l'axe donnée y ou z local"""

        self.name = name_combi
        self.title = f'Barre {bar_id}: Efforts tranchants V{axe}: '
        self.color = 'b'
        self.x_values, self.y_values = self.get_local_internal_forces_coor_value_in_bar(self.bar_info[bar_id]["internals forces"], f"V{axe}", bar_id)
        
        dictEiMinMax = self.max_internal_forces(bar_id)
        self.max_XY_values = (dictEiMinMax[f"V{axe}_max"]["Position"], dictEiMinMax[f"V{axe}_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax[f"V{axe}_min"]["Position"], dictEiMinMax[f"V{axe}_min"]["Effort"])
        self.excent_X = 50
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.x_values
        
        self.save_path = os.path.join(self.SAVE_PATH, f"V{axe}.png")
        self._base_graphique()
        plt.show()

    
    def show_graphique_M(self, bar_id: int, name_combi:str, axe: str=("y", "z")):
        """ Retourne le diagramme des efforts de moment autour de M suivant l'axe y ou z local """

        self.name = name_combi
        self.title = f'Barre {bar_id}: Moments M{axe} : '
        self.color = 'r'
        self.y_label = "Effort (kN.m)"
        self.unit = " kN.m"
        self.x_values, self.y_values = self.get_local_internal_forces_coor_value_in_bar(self.bar_info[bar_id]["internals forces"], f"M{axe}", bar_id)
        
        dictEiMinMax = self.max_internal_forces(bar_id)
        self.max_XY_values = (dictEiMinMax[f"M{axe}_max"]["Position"], dictEiMinMax[f"M{axe}_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax[f"M{axe}_min"]["Position"], dictEiMinMax[f"M{axe}_min"]["Effort"])
        self.excent_X = 1
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.x_values
        
        self.save_path = os.path.join(self.SAVE_PATH, f"M{axe}.png")
        self._base_graphique()
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
        self._effort_interne()
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
    import json
    b = 100
    h = 100
    a = b*h
    Iy = (b*h**3)/12
    Iz = (h*b**3)/12

    # _list_loads = [
    #             [2, '', -83333.3333333, "0/end", 'Permanente G', 'Z'],
    #             ]

    _list_loads = []
    for i in range(3,7):
        print(i)
        _list_loads.append([i, '', -1000, "0/end", 'Permanente G', 'Z'])

    chargement = Chargement(pays="Japon")

    # chargement.add_bar(0,0,0,360,360,0,a)
    # # chargement.add_bar(840,360,0,360,360,0,a)
    # chargement.add_bar(360,360,0,840,360,0,a)

    chargement.add_bar(0,0,0,4000,0,0,a)
    chargement.add_bar(4000,0,0,8000,0,0,a)

    chargement.add_bar(0,0,0,2000,2000,0,a)
    chargement.add_bar(2000,2000,0,4000,4000,0,a)
    chargement.add_bar(4000,4000,0,6000,2000,0,a)
    chargement.add_bar(6000,2000,0,8000,0,0,a)



    chargement.add_bar(4000,0,0,2000,2000,0,a)
    chargement.add_bar(4000,0,0,4000,4000,0,a)
    chargement.add_bar(6000,2000,0,4000,0,0,a)

    chargement.create_load_by_list(_list_loads)
    c1 = Combinaison._from_parent_class(chargement, cat="Cat A : habitation", kdef=0.6)
    # print(c1.list_combination)
    rcombi = "ELU_STR G"
    print(c1.get_combi_list_load(rcombi))


    
    chargement.add_relaxation(7, "start")
    chargement.add_relaxation(7, "end")
    chargement.add_relaxation(8, "start")
    chargement.add_relaxation(9, "start")
    chargement.add_relaxation(9, "end")

    chargement.add_relaxation(4, "end")
    chargement.add_relaxation(5, "start")
    
    chargement.add_relaxation(1, "start")
    chargement.add_relaxation(2, "end")

    for i in range(9):
        chargement.add_material_by_class(i+1, Iy, Iz, "C24")

    # chargement.add_material_by_mechanical_properties(1, 30*10**6, 0, 0, 1000, 1000)
    # chargement.add_material_by_mechanical_properties(2, 30*10**6, 0, 0, 1000, 1000)

    listdeplacement = [
                    [1, "Rotule", "start", 0],
                    [2, "Rotule", "end", 0],
                    [4, "Simple Y", "end", 0],
                    ]
    chargement.create_supports_by_list(listdeplacement)
   
    
    mef = FEM._from_parent_class(chargement, combinaison=c1)
    
    mef.calcul_1D()
    with open("bars_infos.json", "w") as file:
        json.dump(chargement.bar_info, file, indent=4)
    mef.show_graph_loads_and_supports(scale_internal_forces=40, scale_deplacement=15)
    bar = 2
    mef.show_graphique_reaction_Z("elu")
    mef.show_graphique_Nx(bar, "ELU_STR G")
    mef.show_graphique_V(bar, "ELU_STR G", "z")
    mef.show_graphique_M(bar, "ELU_STR G", "y")
    mef.show_graphique_fleche(bar, "ELU_STR G", "z")
    # print(mef.reaction_max())
    # mef.graphique(rcombi)
    # mef.show_graphique_fleche(rcombi)
    # mef.show_graph_loads_and_supports()
