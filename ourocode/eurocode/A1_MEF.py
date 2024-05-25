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
# from eurocode.EC0_Combinaison import Combinaison

from ourocode.eurocode.EC0_Combinaison import Combinaison

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
    
    
        
class MEF(Combinaison, _Base_graph):
    SAVE_PATH = os.path.join(Combinaison.PATH_CATALOG, "data","screenshot")
    def __init__(self, long: int, E: int, A: float, G: float, J: float, Iy: float, Iz: float, ele: int=500, alphaZ: float=0, alphaY: float=0, alphaX: float=0,**kwargs):
        """Classe permettant de créer des poutres MEF et de les calculer. Cette classe est hérité de l'objet Combinaison du module EC0_Combinaison.py
        Args:
            combinaison (Combinaison): L'objet créer à partir de la classe Combinaison du module EC0_Combinaison.py
            long (int): Longueur de l'élément en mm
            E (int): Module de young en MPa, ce E est le E,mean. Il ne faut absolument pas donner le E,mean,fin sous peine de réaliser le calcul EC5 §2.3.2.2 equ2.7 deux fois !
            A (float): Section en mm²
            G (float): Module de cisaillement en MPa
            J (float): Module de torsion
            Iy (float): Inertie quadratique autour de y en mm4
            Iz (float): Inertie quadratique autour de z en mm4
            ele (int): Nombre d'élément fini (mesh)
            alphaZ (float): Angle d'application des charges entre un repère local et global autour de Z
            alphaY (float): Angle d'application des charges entre un repère local et global autour de y
            alphaX (float): Angle d'application des charges entre un repère local et global autour de x
        """
        super(Combinaison, self).__init__(**kwargs)
        _Base_graph.__init__(self)
        self.long = long
        self.A = A
        self.Iy = Iy
        self.Iz = Iz

        self.ea = E*self.A
        self.eiy = E*self.Iy
        self.eiz = E*self.Iz
        self.gj = G*J

        self.beams = {}
        self.list_supports = {}
        self.bi_connected = 0


    def _get_angle_of_beam(self, v1: tuple):
        """Retourne un angle entre deux points celon le cercle trigo

        Args:
            v1 (tuple): vecteur dans le cercle trigonométrique
        """
        ang1 = np.arctan2(*v1[::-1])
        return np.rad2deg(ang1 % (2 * np.pi))
    

    def add_beam(self, x1: int, y1: int, x2: int, y2: int):
        """Ajoute une poutre au model MEF

        Args:
            x1 (int): position de départ en x en mm
            y1 (int): position de départ en y en mm
            x2 (int): position de fin en x en mm
            y2 (int): position de fin en y en mm
        """
        i_beam = len(self.beams) + 1
        v1 = (x2-x1, y2-y1)
        print(v1)
        length = mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2)
        angle = self._get_angle_of_beam(v1)
        self.beams[str(i_beam)] = {"elements": [], "length": length,
                                    "angle": angle, 
                                    "x1": x1, "y1": y1, "x2": x2, "y2": y2}
        self._create_elements_and_nodes(self.beams[str(i_beam)])
        print("Poutre crée: ", self.beams)
        
        
    def _create_elements_and_nodes(self, beam: dict):
        if beam["length"] < 5000:
            nb_ele = int(mt.ceil(beam["length"]/1000))
        else:
            nb_ele = int(mt.ceil(beam["length"]/1000))

        shape = (0,0)
        if hasattr(self, "element_list"):
            shape = (self.element_list.shape[0] + 1 - self.bi_connected, self.element_list.shape[1])

        # On crée une  node list et une coordonée list et on itère
        nodeCoor = np.zeros((nb_ele+1, 2))
        nodeList = np.zeros((nb_ele,2))
        
        node_coor_to_del = []
        bi_connect = False
        j = 0
        for i in range(nb_ele+1):
            length = beam["length"] / nb_ele * i
            x = round(beam["x1"] + length * mt.cos(mt.radians(beam["angle"])),3)
            y = round(beam["y1"] + length * mt.sin(mt.radians(beam["angle"])),3)
            if not x:
                x = 0
            if not y:
                y = 0
            if shape[0]>1:
                find_array = np.where((self.node_coor == [x,y]).all(axis=1))[0]
                print(find_array, [x,y], i, j, shape[0])
                if find_array.shape[0]:
                    # print(self.element_list[find_array[0]])
                    val = 0
                    if find_array[0] != 0:
                        find_array[0] = find_array[0]-1
                        val = 1
                    if i < nb_ele:
                        nodeList[i,0] = self.element_list[find_array[0]][val]
                        nodeList[i,1] = shape[0]+i+1-j
                    else:
                        nodeList[i-1,0] = shape[0]+i-j
                        nodeList[i-1,1] = self.element_list[find_array[0]][val] + self.bi_connected
                        # on test si on a déja rattacher un premier noeud de la barre à un noeud existant si oui on incrémente le compteur
                        if j == 1:
                            self.bi_connected += 1
                            bi_connect = True

                    j += 1
                    node_coor_to_del.append(i)
                    if not bi_connect:
                        beam["elements"].append(shape[0]+i-j+self.bi_connected)
                    print("salu",shape[0],i,j)
                    
                    continue
            nodeCoor[i,0] = x
            nodeCoor[i,1] = y
            
            if i < nb_ele:
                nodeList[i,0] = shape[0]+i+1 - j
                nodeList[i,1] = shape[0]+i+2 - j
                print(shape[0],i,j)
                beam["elements"].append(shape[0]+i-j+self.bi_connected)

        # print("noeud à supprimer", node_coor_to_del, nodeCoor, nodeCoor.shape)
        j = 0
        for i_del in node_coor_to_del:
            nodeCoor = np.delete(nodeCoor, i_del - j, 0)
            j += 1
        # for i in range(shape[0], shape[0]+nb_ele):
        #     nodeList[i-shape[0],0] = i+1
        #     nodeList[i-shape[0],1] = i+2
        #     beam["elements"].append(i)
        
        
      
        if shape[0]:
            self.node_coor = np.concatenate((self.node_coor, nodeCoor), axis=0)
            self.element_list = np.concatenate((self.element_list, nodeList), axis=0)
        else:
            self.element_list = nodeList
            self.node_coor = nodeCoor
        print(self.element_list, "\n", self.node_coor)
        

    # def node_coord(self):
    #     nodeCoor = np.zeros((self.ele+1, 1))
    #     for i in range(self.ele + 1):
    #         # l0 = 6
    #         # if i == 1 or i == self.ele-1:
    #         #     if i == 1:
    #         #         nodeCoor[i,0]= l0
    #         #     else:
    #         #         nodeCoor[i,0]= round(self.long - l0)
    #         # elif i == 0:
    #         #     nodeCoor[i,0]= 0
    #         # elif i == self.ele:
    #         #     nodeCoor[i,0]= round(self.long)
    #         # else:
    #         #     nodeCoor[i,0]= round(((self.long- l0) / (self.ele-2)) + nodeCoor[i-1,0])
    #         nodeCoor[i,0]= round(self.long/ self.ele * i)

    #     # print("\n Liste des coordonnées des éléments : \n", nodeCoor, "\n")
    #     return nodeCoor


    # def element_list(self):
    #     nodeList = np.zeros((self.ele,2))
    #     for i in range(self.ele):
    #         nodeList[i,0] = i+1
    #         nodeList[i,1] = i+2
        
    #     # print("\n Liste des éléments : \n",nodeList, "\n")
    #     return nodeList
    

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
    

    def _init_matrix_K(self, index: int|str) -> np.array:
        """Initialise la matrice K suivant la longueur de l'élément l en mm
        Args:
            l (int): longueur de l'élément MEF en mm
        Returns:
            numpy.array: matrice de rigidité K
        """
        for key, beam in self.beams.items():
            if index in beam["elements"]:
                angle = beam["angle"]
                l = beam["length"] / len(beam["elements"])
                print(l)
                print(f"J'ai retrouvé la poutre {key} qui est associé à l'élément {index}, son angle est de {angle} degrés.")
                break
        
        k = np.array([[(self.ea)/l, 0, 0, 0, 0, 0, -(self.ea)/l, 0, 0, 0, 0, 0],
                    [0, (12*self.eiz)/l**3, 0, 0, 0, (6*self.eiz)/l**2, 0, -(12*self.eiz)/l**3, 0, 0, 0, (6*self.eiz)/l**2],
                    [0, 0, (12*self.eiy)/l**3, 0, -(6*self.eiy)/l**2, 0, 0, 0, -(12*self.eiy)/l**3, 0, -(6*self.eiy)/l**2, 0],
                    [0, 0, 0, self.gj/l, 0, 0, 0, 0, 0, -self.gj/l, 0, 0],
                    [0, 0, -(6*self.eiy)/l**2, 0, (4*self.eiy)/l, 0, 0, 0, (6*self.eiy)/l**2, 0, (2*self.eiy)/l, 0],
                    [0, (6*self.eiz)/l**2, 0, 0, 0, (4*self.eiz)/l, 0, -(6*self.eiz)/l**2, 0, 0, 0, (2*self.eiz)/l],
                    [-(self.ea)/l, 0, 0, 0, 0, 0, (self.ea)/l, 0, 0, 0, 0, 0],
                    [0, -(12*self.eiz)/l**3, 0, 0, 0, -(6*self.eiz)/l**2, 0, (12*self.eiz)/l**3, 0, 0, 0, -(6*self.eiz)/l**2],
                    [0, 0, -(12*self.eiy)/l**3, 0, (6*self.eiy)/l**2, 0, 0, 0, (12*self.eiy)/l**3, 0, (6*self.eiy)/l**2, 0],
                    [0, 0, 0, -self.gj/l, 0, 0, 0, 0, 0, self.gj/l, 0, 0],
                    [0, 0, -(6*self.eiy)/l**2, 0, (2*self.eiy)/l, 0, 0, 0, (6*self.eiy)/l**2, 0, (4*self.eiy)/l,0],
                    [0, (6*self.eiz)/l**2, 0, 0, 0, (2*self.eiz)/l, 0, -(6*self.eiz)/l**2, 0, 0, 0, (4*self.eiz)/l]])
        
        # self.T_matrix = self._transformation_matrix(angle)
        T_matrix = self._transformation_matrix(angle)

        self.k_local["ele"+str(index)] = {"k_local": k, "T_global": T_matrix}
        return np.dot(np.dot(T_matrix.T.toarray(), k), T_matrix.toarray())


    def _agregation_matrix_K(self, assembleur_COO_K: list, k: np.array, n1: int, n2: int):
        # Agrégation de la matrice global
        for i in range(-6,0):
            for j in range(-6,0):
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
        self.k_local = {}

        for i in range(len(self.element_list)):
            # matrice de rigidité K
            k = self._init_matrix_K(i)
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6
            
            self._agregation_matrix_K(Assembleur_COO_K , k, n1, n2)
        
        self.matriceK = (coo_matrix(Assembleur_COO_K.data, shape=((len(self.element_list)+1)*6, (len(self.element_list)+1)*6))).tocsr() 
             
        # print("\n Matrice de rigidité K  : \n", self.matriceK.toarray(), "\n")
        return self.matriceK 


    def _nearest_node(self, value: int) -> list:
        """Find the nearest node of the beam in function of the value"""
        # element to which nearest value is to be found
        #print("Value to which nearest element is to be found: ", value)
        
        # calculate the difference array
        difference_array = np.absolute(self.node_coor-value)
        
        # find the index of minimum element from the array
        index = difference_array.argmin()
        print("nearest", index)
        nearest_distance = self.node_coor[index][0]
        index_I_matrice = ((index+1) * 6) - 6
        element_coor = self.element_list[index-1]
        
        #print("Nearest element to the given values is : ", nearest_distance)
        #print("Index of nearest value is : ", index)
        return index, nearest_distance, index_I_matrice, element_coor


    def _matrice_U_1D(self, listDepla):
        """ Donner en entrée la liste type des appuis extérieurs, 
        ressort une matrice acceptable pour le calcul MEF """
        self.matriceU = np.ones(((len(self.element_list)+1)*6,1))
        
        for i in range(len(listDepla)):
            index_min = self._nearest_node(int(listDepla[i][2]))
           
            if listDepla[i][1] == 'Rotule':
                for i in range(0,4):
                    self.matriceU[index_min[2] + i,0] = False
                    
            elif listDepla[i][1] == 'Encastrement':
                for i in range(0,6):
                    self.matriceU[index_min[2] + i,0] = False
            else:
                self.matriceU[index_min[2] + 2,0] = False
                self.matriceU[index_min[2] + 3,0] = False

        # print("\n Matrice de déplacement U : \n",self.matriceU, "\n")
        return self.matriceU


    def _matrice_Fext_1D(self, listFext):
        """ Donner en entrée la liste type des efforts extérieurs en daN, 
        ressort une matrice acceptable pour le calcul MEF """
        Assembleur_COO_F = _Triplets()
        
        self.alpha_R_tetaY = 0
        base_global_R_tetaY = np.array([[np.cos(self.alpha_R_tetaY), np.sin(self.alpha_R_tetaY), 0],
                    [-np.sin(self.alpha_R_tetaY), np.cos(self.alpha_R_tetaY), 0],
                    [0, 0, 1]])
        
        # l = abs(int(self.node_coor[1,0]))
        l = abs((self.long/ len(self.element_list)))

        for i in range(len(listFext)):
            listFlocal = np.zeros((3,1))
            listFlocal2 = np.zeros((3,1))
            
            if listFext[i][3] == 'Linéique':
            
                slash_pos = listFext[i][5].index("/")
                start_pos = int(listFext[i][5][0:slash_pos])
                start_index = self._nearest_node(start_pos)
                
                end_pos = int(listFext[i][5][slash_pos+1:])
                end_index = self._nearest_node(end_pos)
                index_min = [[start_index, start_pos], [end_index, end_pos]]

                # print("Longueur élément:", l, " mm")
                for index in range(start_index[2], end_index[2]+6, 6):
                    listFlocal = np.zeros((3,1))
                    listFlocal2 = np.zeros((3,1))

                    if listFext[i][6] == 'Z' or listFext[i][6] == 'Z local' :
                        force = listFext[i][4] * (l/1000) * 10

                        listFlocal[1,0] = listFext[i][4] * 10 * (l/1000) / 2 # Rza et Rzb
                        listFlocal2[1,0] = listFlocal[1,0]
                        listFlocal[2,0] = listFext[i][4] * 10 * (l/1000) **2 / 12 # Mya et Myb
                        listFlocal2[2,0] = listFlocal[2,0]

                        if index in (start_index[2], end_index[2]):
                            calc_formulaire = False
                            if index != end_index[2]:
                                ind = index_min[0]

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
                                ind = index_min[1]
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
                                if listFext[i][6] == 'Z':
                                    listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                                    listFlocal2= np.dot(base_global_R_tetaY, listFlocal2)

                                for i_index in range(3):
                                    j_index = i_index * 2
                                    Assembleur_COO_F.append(index+j_index, 0, listFlocal[i_index,0])
                                    Assembleur_COO_F.append(index+6+j_index, 0, listFlocal2[i_index,0])
                                continue

                        if listFext[i][6] == 'Z':
                            listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                            listFlocal2= np.dot(base_global_R_tetaY, listFlocal2)
                            
                    elif listFext[i][6] == 'X':
                        listFlocal[0,0] = listFext[i][4] * (l/1000) * 10 / 2
                        listFlocal2[0,0] = listFlocal[0,0]
                        listFlocal= np.dot(base_global_R_tetaY, listFlocal)

                    if index != end_index[2]:
                        for i_index in range(3):
                            j_index = i_index * 2
                            Assembleur_COO_F.append(index+j_index, 0, listFlocal[i_index,0])
                            Assembleur_COO_F.append(index+6+j_index, 0, listFlocal2[i_index,0])

            
            elif listFext[i][3] == 'Nodale':
                index_min = self._nearest_node(int(listFext[i][5]))
                force = listFext[i][4] * 10
                
                if listFext[i][6] == 'Z' or listFext[i][6] == 'Z local' :
                    if listFext[i][5] == index_min[1]:
                        listFlocal[1,0] = force
                    
                    # Cas dans lequel la force est aprés le noeud le plus proche
                    elif listFext[i][5] > index_min[1]:
                        
                        a = listFext[i][5] - index_min[1]
                        b = self.node_coor[index_min[0]+1][0]-listFext[i][5]
                        Mya = (force * a * b**2)/l**2
                        Myb = -(force * b * a**2)/l**2
                        Rza = (force * b**2 * (3 * a + b))/l**3
                        Rzb = (force * a**2 * (3 * b + a))/l**3
                        
                        listFlocal[1,0] = Rza
                        listFlocal[2,0] = Mya
                        
                        Assembleur_COO_F.append(index_min[2]+6+2, 0, Rzb)
                        Assembleur_COO_F.append(index_min[2]+6+4, 0, Myb)
                        
                        #print("La force est situé après le noeud le plus proche à: ", a,b)
                    
                    # Cas dans lequel la force est avant le noeud le plus proche
                    else:
                        b = index_min[1] - listFext[i][5]
                        a = listFext[i][5] - self.node_coor[index_min[0]-1][0]
                        Mya = (force * a * b**2)/l**2
                        Myb = -(force * b * a**2)/l**2
                        Rza = (force * b**2 * (3 * a + b))/l**3
                        Rzb = (force * a**2 * (3 * b + a))/l**3
                        
                        listFlocal[1,0] = Rzb
                        listFlocal[2,0] = Myb
                        
                        Assembleur_COO_F.append(index_min[2]-6+2, 0, Rza)
                        Assembleur_COO_F.append(index_min[2]-6+4, 0, Mya)
                        #print("La force est situé avant le noeud le plus proche à: ", a,b)
                    
                elif listFext[i][6] == 'X':
                    listFlocal[0,0] = listFext[i][4] * 10
                    
                if listFext[i][6] == 'Z local':
                    pass
                else:
                    listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                
                Assembleur_COO_F.append(index_min[2], 0, listFlocal[0,0])
                Assembleur_COO_F.append(index_min[2]+2, 0, listFlocal[1,0])
                Assembleur_COO_F.append(index_min[2]+4, 0, listFlocal[2,0])
                
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
        print("\n Solution Déplacement : \n", self.matriceU, "\n")
        return self.matriceU


    def _equa_reaction(self):
        self.ri = np.dot(self.matriceK.toarray(), self.matriceU) - self.matriceF.toarray()
        print("\n Solution Réaction : \n", self.ri, "\n")
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
        ux = []
        uz = []
        tetay = []

        for i in range(len(self.element_list)):
            n1 = int(self.element_list[i,0])*6
            n2 = int(self.element_list[i,1])*6

            u1 = self.matriceU[n1-6:n1,0]
            u2 = self.matriceU[n2-6:n2,0]

            if n1 == 6:

                ux.append(u1[0])
                ux.append(u2[0])
                uz.append(u1[2])
                uz.append(u2[2])
                tetay.append(u1[4])
                tetay.append(u2[4])
            else:
                ux.append(u2[0])
                uz.append(u2[2])
                tetay.append(u2[4])

        self.u_coor = [ux, uz, tetay]
        return self.u_coor
        

    def deplacement_max(self):

        ux_max = max(self.u_coor[0])
        ux_max_index = self.u_coor[0].index(ux_max)
        ux_max_coor = [self.node_coor[ux_max_index,0], ux_max]

        ux_min = min(self.u_coor[0])
        ux_min_index = self.u_coor[0].index(ux_min)
        ux_min_coor = [self.node_coor[ux_min_index,0], ux_min]

        uz_max = max(self.u_coor[1])
        uz_max_index = self.u_coor[1].index(uz_max)
        uz_max_coor = [self.node_coor[uz_max_index,0], uz_max]

        uz_min = min(self.u_coor[1])
        uz_min_index = self.u_coor[1].index(uz_min)
        uz_min_coor = [self.node_coor[uz_min_index,0], uz_min]

        tetay_max = max(self.u_coor[2])
        tetay_max_index = self.u_coor[2].index(tetay_max)
        tetay_max_coor = [self.node_coor[tetay_max_index,0], tetay_max]

        tetay_min = min(self.u_coor[2])
        tetay_min_index = self.u_coor[2].index(tetay_min)
        tetay_min_coor = [self.node_coor[tetay_min_index,0], tetay_min]

        listMinMax = {"Ux_max": ux_max_coor, "Ux_min": ux_min_coor, 
                      "Uz_max": uz_max_coor, "Uz_min": uz_min_coor, 
                      "tetay_max": tetay_max_coor, "tetay_min": tetay_min_coor}

        return listMinMax


    def effort_interne(self):

        ni = []
        vi = []
        myi = []
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
            internal_forces = np.dot(np.dot(self.k_local["eleALL"], self.T_matrix.toarray()), u_global) - np.dot(self.T_matrix.toarray(), f_0)
            # print("effort interne: ", internal_forces)

            if n1 == 6:
                ni.append(internal_forces[0]/10**3)
                vi.append(internal_forces[2]/10**3)
                myi.append(internal_forces[4]/10**6)
            
            # else:
            #     ni[-1] = ni[-1] + internal_forces[0]/10**3
            #     vi[-1] = vi[-1] + internal_forces[2]/10**3
            #     myi[-1] = myi[-1] + internal_forces[4]/10**6
            
            ni.append(-internal_forces[6]/10**3)
            vi.append(-internal_forces[8]/10**3)
            myi.append(-internal_forces[10]/10**6)
                
        self.ei_coor = [ni, vi, myi]
        return self.ei_coor
    
    

    def effort_interne_position(self, position: int):
        position_index = self._nearest_node(position)
        nx = self.ei_coor[0][position_index[0]]
        vz = self.ei_coor[1][position_index[0]]
        my = self.ei_coor[2][position_index[0]]
        print(f"la position la plus proche est {position_index[1]} mm", {"Nx": nx, "Vz": vz, "My": my})
        return {"Nx": nx, "Vz": vz, "My": my}
    

    
    def effort_interne_max(self):
        """ Retourne les efforts min et max d'une liste d'effort interne """

        # def max_min_list_load(get_max: bool=True, nested_list=[]):
        #     item_value = 0
        #     index_item_value = 0
        #     for index, item in enumerate(nested_list):
        #         if isinstance(item, list):
        #             if get_max:
        #                 item = max(item)
        #             else:
        #                 item = min(item)
        #         if get_max and (item > item_value):
        #             item_max = item
        #             index_item_value = index
        #         elif not get_max and (item < item_value):
        #             item_max = item
        #             index_item_value = index
        #     return item_value, index_item_value
                

        # nx_max, nx_max_index = max_min_list_load(get_max=True, nested_list=self.ei_coor[0])
        nx_max = max(self.ei_coor[0])
        nx_max_index = self.ei_coor[0].index(nx_max)
        nx_max_coor = {"Position": self.node_coor[nx_max_index,0], "Effort": nx_max}

        # nx_min, nx_min_index = max_min_list_load(get_max=False, nested_list=self.ei_coor[0])
        nx_min = min(self.ei_coor[0])
        nx_min_index = self.ei_coor[0].index(nx_min)
        nx_min_coor = {"Position": self.node_coor[nx_min_index,0], "Effort": nx_min}

        # vz_max, vz_max_index = max_min_list_load(get_max=True, nested_list=self.ei_coor[1])
        vz_max = max(self.ei_coor[1])
        vz_max_index = self.ei_coor[1].index(vz_max)
        vz_max_coor = {"Position": self.node_coor[vz_max_index,0], "Effort": vz_max}

        # vz_min, vz_min_index = max_min_list_load(get_max=False, nested_list=self.ei_coor[1])
        vz_min = min(self.ei_coor[1])
        vz_min_index = self.ei_coor[1].index(vz_min)
        vz_min_coor = {"Position": self.node_coor[vz_min_index,0], "Effort": vz_min}

        # my_max, my_max_index = max_min_list_load(get_max=True, nested_list=self.ei_coor[2])
        my_max = max(self.ei_coor[2])
        my_max_index = self.ei_coor[2].index(my_max)
        my_max_coor = {"Position": self.node_coor[my_max_index,0], "Effort": my_max}

        # my_min, my_min_index = max_min_list_load(get_max=False, nested_list=self.ei_coor[2])
        my_min = min(self.ei_coor[2])
        my_min_index = self.ei_coor[2].index(my_min)
        my_min_coor = {"Position": self.node_coor[my_min_index,0], "Effort": my_min}

        return {"Nx_max": nx_max_coor,"Nx_min": nx_min_coor, 
                "Vz_max": vz_max_coor, "Vz_min": vz_min_coor, 
                "My_max": my_max_coor, "My_min": my_min_coor}



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
        self.y_values = self.u_coor[1]
        
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
        self.y_values = self.ei_coor[0]
        
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


    def _graphique_Vz(self, name_combi:str):
        """ Retourne le diagramme des efforts tranchant en Z """

        self.name = name_combi
        self.title = 'Efforts tranchants Vz : '
        self.color = 'b'
        # self.x_values = []
        self.x_values = self.node_coor
        self.y_values = self.ei_coor[1]
        # self.y_values = []
        # for index, value in enumerate(self.ei_coor[1]):
        #     if isinstance(value, list):
        #         for i in range(len(value)):
        #             self.x_values.append(self.node_coor[index,0])
        #             self.y_values.append(value[i])
        #     else:
        #         self.x_values.append(self.node_coor[index,0])
        #         self.y_values.append(value[i])
        
        dictEiMinMax = self.effort_interne_max()
        self.max_XY_values = (dictEiMinMax["Vz_max"]["Position"], dictEiMinMax["Vz_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax["Vz_min"]["Position"], dictEiMinMax["Vz_min"]["Effort"])
        self.excent_X = 50
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, "Vz.png")
        self._base_graphique()
    

    def show_graphique_Vz(self, name_combi:str):
        self._graphique_Vz(name_combi)
        plt.show()

    
    def _graphique_My(self, name_combi:str):
        """ Retourne le diagramme des efforts de moment autour de My """

        self.name = name_combi
        self.title = 'Moments My : '
        self.color = 'r'
        self.y_label = "Effort (kN.m)"
        self.unit = " kN.m"
        self.x_values = self.node_coor
        self.y_values = self.ei_coor[2]
        
        dictEiMinMax = self.effort_interne_max()
        self.max_XY_values = (dictEiMinMax["My_max"]["Position"], dictEiMinMax["My_max"]["Effort"])
        self.min_XY_values = (dictEiMinMax["My_min"]["Position"], dictEiMinMax["My_min"]["Effort"])
        self.excent_X = 1
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, "My.png")
        self._base_graphique()

    def show_graphique_My(self, name_combi):
        self._graphique_My(name_combi)
        plt.show()

    
    def show_graph_loads_and_supports(self):
        """Affiche le graphique des charges et des appuis du MEF
        """
        # Conversion des positions en valeurs numériques
        def parse_position(position, charge):
            position = str(position)
            charge = -charge
            parts = position.split("/")
            if len(parts) == 1:
                return [int(parts[0])]*2, [0,charge]
            elif len(parts) == 2:
                return [pos for pos in range(int(parts[0]), int(parts[1]), 1)], [charge for _ in range(int(parts[0]), int(parts[1]), 1)]
            else:
                return None

        # Création du graphique
        plt.figure(figsize=(12, 4))

        plt.arrow(0, 0, 0, 350, width=15, color="blue")
        plt.text(150, -100, "X", ha='right')
        plt.arrow(0, 0, 350, 0, width=15, color="red")
        plt.text(-50, 150, "Z", ha='right')
        for i, load in enumerate(self.list_loads):
            charge = load[4]
            parser = parse_position(load[5], charge)
            charge = round(-charge,2)
            nom = " / ".join([load[1], load[2], load[3], load[6]])
            if len(parser[1]) != 2:
                unit_load = "daN/m"
                plt.plot(parser[0], parser[1], label=nom)
                plt.fill_between(parser[0], parser[1], alpha=0.3)
            else:
                unit_load = "daN"
                plt.plot(parser[0], parser[1], marker="X",label=nom)
            plt.text(parser[0][1]+1000, charge+2, f'{charge} {unit_load}', ha='right')

        for key, beam in self.beams.items():
            x, y = [], []
            for element in beam["elements"]:
                for i, node in enumerate(self.element_list[element]):
                    coor = self.node_coor[int(node)-1]
                    x.append(coor[0])
                    y.append(coor[1])
            plt.plot(x, y, label=f"Poutre N°{key}")
            plt.plot(x[0], y[0], marker='o', mec='gray', mfc="gray")
            plt.plot(x[-1], y[-1], marker='o', mec='gray', mfc="gray")

        for key, support in self.list_supports.items():
            if support["Type d'appui"] == "Rotule":
                support_type = "o"
            elif support["Type d'appui"] == "Encastrement":
                support_type = "s"
            else:
                support_type = "^"
            x = mt.cos(mt.radians(self.beams[str(support["N° poutre"])]["angle"])) * support["Position de l'appui"] + self.beams[str(support["N° poutre"])]["x1"]
            y = mt.sin(mt.radians(self.beams[str(support["N° poutre"])]["angle"])) * support["Position de l'appui"] + self.beams[str(support["N° poutre"])]["y1"]
            plt.plot(x, y, marker=support_type, markersize=10, color="red", label=f"Appui {key} / {support["Type d'appui"]}")


        plt.title('Schématisation de la structure et des charges')
        plt.xlabel('Longueur (mm)')
        plt.ylabel('Hauteur (mm')
        plt.legend()
        plt.grid(True)
        plt.show()


    def get_supports(self):
        """Retourne la liste des appuis définis.
        """
        return self.list_supports

    def create_support(self, beam_number: int, type_appuis: str=("Simple", 'Rotule', 'Encastrement'), pos: int=0, l_appuis: int=0):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            type_appuis (str, optional): type d'appuis à créer. Defaults to ("Simple", 'Rotule', 'Encastrement').
            pos (int, optional): position de l'appuis sur la poutre en mm. Defaults to 0.
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        self.list_supports[str(len(self.list_supports)+1)] = {"N° poutre": beam_number, 
                                                              "Type d'appui": type_appuis, 
                                                              "Position de l'appui": pos, 
                                                              "Longueur d'appui": l_appuis}
        return self.list_supports[str(len(self.list_supports)+1)]

    
    def create_supports_by_list(self, list_supports: list):
        """Ajoute les charges d'une liste pré-défini dans la liste de chargement

        Args:
            list_supports (list): liste de charge.
        """
        for support in list_supports:
            self.list_supports[str(len(self.list_supports)+1)] = {"N° poutre": support[0], 
                                                                  "Type d'appui": support[1], 
                                                                  "Position de l'appui": support[2], 
                                                                  "Longueur d'appui": support[3]}
        return self.list_supports


    def del_support(self, index_load: int):
        """Supprime un appui de l'attribut list_supports par son index

        Args:
            index_load (int): index de l'appui à supprimer.
        """
        return f"L'appui ci-joint à été supprimé: {self.list_supports.pop(str(index_load))}"
    


    def calcul_1D(self):
        """Calcul une poutre MEF 1D

        Args:
            list_loads (list): liste des charges combinées sur la poutre MEF.
        """
        start = time.time()
        self._matrice_K_1D()
        self._matrice_Fext_1D(self.list_loads)
        self._matrice_U_1D(self.list_supports)

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
            self.show_graphique_Vz(name_combi)
            self.show_graphique_My(name_combi)
        else:
            self.show_graphique_fleche(name_combi)
        

if __name__ == '__main__':
    from EC0_Combinaison import Chargement
    # _list_loads = [[1, '', 'Permanente G', 'Linéique', -100, "0/6000", 'Z'],
    #              [2, '', "Neige normale Sn", 'Linéique', -200, "0/6000", 'Z']]
    _list_loads = [[1, '', 'Permanente G', "Nodale", -500, 2500, 'X'],
                [2, '', "Neige normale Sn", 'Linéique', -200, "0/5000", 'Z']]
    chargement = Chargement(pays="Japon")
    chargement.create_load_by_list(_list_loads)
    c1 = Combinaison._from_parent_class(chargement, cat="Cat A : habitation", kdef=0.6)
    # print(c1.list_combination)
    rcombi = "ELU_STR G"
    print(c1.get_combi_list_load(rcombi))
    long = 5000
    ele = 6
    
    b = 60
    h = 100
    a = b*h
    iy = (b*h**3)/12
    iz = (h*b**3)/12
    test = MEF._from_parent_class(c1, long=long,E=11000,A=a, G=690, J=690, Iy=iy, Iz=iz, ele=ele, alphaZ=0, alphaY=0, alphaX=0)

    listdeplacement = [[1, "Rotule", 0, 0], [1, "Rotule", 3000, 0]]
    test.create_supports_by_list(listdeplacement)
    
    test.add_beam(0,0,5000,0)
    test.add_beam(0,0,2500,4000)
    test.add_beam(5000,0,2500,4000)
    test.show_graph_loads_and_supports()
    # test.calcul_1D()
    # print(test.reaction_max())
    # test.graphique(rcombi)
    # test.show_graphique_fleche(rcombi)
    # test.show_graph_loads_and_supports()
