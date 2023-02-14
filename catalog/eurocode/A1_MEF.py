#! env/scripts/python.exe
#! env/scripts/python.exe
#coding in UTF8

import os

import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

from EC0_Combinaison import Combinaison

class _Base_graph(object): 
    """ Retourne un diagramme de base """
    SAVE_PATH = os.path.join(os.getcwd(), "catalog", "data","screenshot")

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
    def __init__(self, long: int, E: int, A: float, G: float, J: float, Iy: float, Iz: float, ele: int=500, alphaZ: float=0, alphaY: float=0, alphaX: float=0,**kwargs):
        """Classe permettant de créer des poutres MEF et de les calculer. Cette classe est hérité de l'objet Combinaison du module EC0_Combinaison.py
        Args:
            long (int): Longueur de l'élément en mm
            E (int): Module de young en MPa
            A (float): Section en mm²
            G (float): Module de cisaillement en MPa
            J (float): Module de torsion
            Iy (float): Inertie quadratique autour de y en mm4
            Iz (float): Inertie quadratique autour de z en mm4
            ele (int): Nombre d'élément MEF
            alphaZ (float): Angle d'application des charges entre un repère local et global autour de Z
            alphaY (float): Angle d'application des charges entre un repère local et global autour de y
            alphaX (float): Angle d'application des charges entre un repère local et global autour de x
        """
        super(Combinaison, self).__init__( **kwargs)
        _Base_graph.__init__(self)
        self.long = long
        self.A = A
        self.Iy = Iy
        self.Iz = Iz

        self.ea = E*self.A
        self.eiy = E*self.Iy
        self.eiz = E*self.Iz
        self.gj = G*J
        self.ele = ele
        self.alpha_R_tetaZ = np.radians(alphaZ)
        self.alpha_R_tetaY = np.radians(alphaY)
        self.alpha_R_tetaX = np.radians(alphaX)
        self.node_coor = self.node_coord()
        self.elementList = self.element_list()
        self.list_supports = []


    def node_coord(self):
        nodeCoor = np.zeros((self.ele+1, 1))
        for i in range(self.ele + 1) :
            # l0 = 6
            # if i == 1 or i == self.ele-1:
            #     if i == 1:
            #         nodeCoor[i,0]= l0
            #     else:
            #         nodeCoor[i,0]= round(self.long - l0)
            # elif i == 0:
            #     nodeCoor[i,0]= 0
            # elif i == self.ele:
            #     nodeCoor[i,0]= round(self.long)
            # else:
            #     nodeCoor[i,0]= round(((self.long- l0) / (self.ele-2)) + nodeCoor[i-1,0])
            nodeCoor[i,0]= round(self.long/ self.ele * i)

        print("\n Liste des coordonnées des éléments : \n", nodeCoor, "\n")
        return nodeCoor


    def element_list(self):
        nodeList = np.zeros((self.ele,2))
        for i in range(self.ele):
            nodeList[i,0] = i+1
            nodeList[i,1] = i+2
        
        print("\n Liste des éléments : \n",nodeList, "\n")
        return nodeList
    

    def _init_matrix_K(self, index: int|str, l: int) -> np.array:
        """Initialise la matrice K suivant la longueur de l'élément l en mm
        Args:
            l (int): longueur de l'élément MEF en mm
        Returns:
            numpy.array: matrice de rigidité K
        """
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

        self.k_local["ele"+str(index)] = {"k11_local": k[0:6,0:6], 
                                        "k12_local": k[0:6,6:12],
                                        "k21_local": k[6:12,0:6],
                                        "k22_local": k[6:12,6:12]}
        return k

    
    def _agregation_matrix_K(self, assembleur_COO_K: list, k: np.array, n1: int, n2: int):
        # Agrégation de la matrice global
        for i in range(-6,0):
                for j in range(-6,0):
                    assembleur_COO_K.append(n1+i, n1+j, k[i+6,j+6]+0.001) # ajout de bruit pour éviter les singularité (+0.001)
                    assembleur_COO_K.append(n2+i, n2+j, k[i+12,j+12]+0.001) # ajout de bruit pour éviter les singularité (+0.001)
                    if k[i+6,j+12]:
                        assembleur_COO_K.append(n1+i, n2+j, k[i+6,j+12])
                        
                    if k[i+12,j+6]:
                        assembleur_COO_K.append(n2+i, n1+j, k[i+12,j+6]+0.001) # ajout de bruit pour éviter les singularité (+0.001)
        

    
    def _matrice_K_1D(self):
        """Méthodequi créer la matrice de rigidité K assemblée
        Returns:
            numpy.array: la matrice de rigidité assemblée
        """
        self.k_local = {}
        l = abs(int(self.node_coor[1,0]))
        
        # matrice de rigidité K
        k = self._init_matrix_K("ALL", l)
        
        # Agrégation de la matrice global
        
        Assembleur_COO_K = _Triplets()
        
        for i in range(self.ele):
            n1 = int(self.elementList[i,0])*6
            n2 = int(self.elementList[i,1])*6
            
            self._agregation_matrix_K(Assembleur_COO_K , k, n1, n2)
        
        self.matriceK = (coo_matrix(Assembleur_COO_K.data, shape=((self.ele+1)*6, (self.ele+1)*6))).tocsr() 
             
        #print("\n Matrice de rigidité K  : \n", self.matriceK, "\n")
        return self.matriceK 


    def _nearest_node(self, value: int) -> list:
        """Find the nearest node of the beam in function of the value"""
        # element to which nearest value is to be found
        #print("Value to which nearest element is to be found: ", value)
        
        # calculate the difference array
        difference_array = np.absolute(self.node_coor-value)
        
        # find the index of minimum element from the array
        index = difference_array.argmin()
        nearest_distance = self.node_coor[index][0]
        index_I_matrice = ((index+1) * 6) - 6
        element_coor = self.elementList[index-1]
        
        #print("Nearest element to the given values is : ", nearest_distance)
        #print("Index of nearest value is : ", index)
        return index, nearest_distance, index_I_matrice, element_coor


    def _matrice_U_1D(self, listDepla):
        """ Donner en entrée la liste type des appuis extérieurs, 
        ressort une matrice acceptable pour le calcul MEF """
        self.matriceU = np.ones(((self.ele+1)*6,1))
        
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

        #print("\n Matrice de déplacement U : \n",self.matriceU, "\n")
        return self.matriceU



    def _matrice_Fext_1D(self, listFext):
        """ Donner en entrée la liste type des efforts extérieurs en daN, 
        ressort une matrice acceptable pour le calcul MEF """
        Assembleur_COO_F = _Triplets()
        
        base_global_R_tetaY = np.array([[np.cos(self.alpha_R_tetaY), np.sin(self.alpha_R_tetaY), 0],
                    [-np.sin(self.alpha_R_tetaY), np.cos(self.alpha_R_tetaY), 0],
                    [0, 0, 1]])
        
        base_global_R_tetaX = np.array([[np.cos(self.alpha_R_tetaX), np.sin(self.alpha_R_tetaX), 0],
                    [-np.sin(self.alpha_R_tetaX), np.cos(self.alpha_R_tetaX), 0],
                    [0, 0, 1]])
        
        l = abs(int(self.node_coor[1,0]))
        
        for i in range(len(listFext)):
            listFlocal = np.zeros((3,1))
            
            if listFext[i][3] == 'Linéique':
            
                slash_pos = listFext[i][5].index("/")
                start_pos = int(listFext[i][5][0:slash_pos])
                start_index = self._nearest_node(start_pos)
                
                end_pos = int(listFext[i][5][slash_pos+1:])
                end_index = self._nearest_node(end_pos)
                index_min = [[start_index, start_pos], [end_index, end_pos]]
                
                
                if listFext[i][6] == 'Z' or listFext[i][6] == 'Z local' :
                    j = 0
                    listFlocal[1,0] = listFext[i][4] * (l/1000) * 10
                    for ind in index_min:
                        force = listFext[i][4] * (l/1000) * 10
                        if ind[1] == ind[0][1]:
                            continue
                        
                        # Formule tirée de l'aide mémoire page 113 cas 3 et cas 1 pour RA et RB
                        # Cas dans le quel la force est aprés le noeud le plus proche
                        elif ind[1] > ind[0][1]:
                            
                            if j:
                            
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
                                
                                Assembleur_COO_F.append(ind[0][2]+2, 0, listFlocal[1,0])
                                Assembleur_COO_F.append(ind[0][2]+4, 0, listFlocal[2,0])
                                Assembleur_COO_F.append(ind[0][2]+6+2, 0, Rzb)
                                Assembleur_COO_F.append(ind[0][2]+6+4, 0, Myb)
                                
                                
                            
                        # Cas dans le quelle la force est avant le noeud le plus proche
                        else:
                            if not j:
                                
                                a = ind[0][1] - ind[1]
                                b = ind[1] - self.node_coor[ind[0][0]-1][0]
                                Mya = ((force * b**2)/l) * (4*(b/l) + 3*(b/l)**2)
                                Myb = -((force * b**2)/12) * (6-8*(b/l) + 3*(b/l)**2)
                                force = force * b
                                b = b/2 
                                a = a + b 
                                Rza = (force * b**2 * (3 * a + b))/l**3
                                Rzb = (force * a**2 * (3 * b + a))/l**3
                                
                                listFlocal[1,0] = Rzb
                                listFlocal[2,0] = Myb
                                
                                Assembleur_COO_F.append(ind[0][2]+2, 0, listFlocal[1,0])
                                Assembleur_COO_F.append(ind[0][2]+4, 0, listFlocal[2,0])
                                Assembleur_COO_F.append(ind[0][2]-6+2, 0, Rza)
                                Assembleur_COO_F.append(ind[0][2]-6+4, 0, Mya)
                        
                        j += 1   
                    
                elif listFext[i][6] == 'X':
                    listFlocal[0,0] = listFext[i][4] * (l/1000) * 10
                
                if listFext[i][6] == 'Z local':
                    pass
                else:
                    listFlocal= np.dot(base_global_R_tetaY, listFlocal)
                
                for index in range(start_index[2], end_index[2]+6, 6):
                    Assembleur_COO_F.append(index, 0, listFlocal[0,0])
                    if index != start_index[2] or index != end_index[2]+6:
                       Assembleur_COO_F.append(index+2, 0, listFlocal[1,0])
            
            
            elif listFext[i][3] == 'Nodale':
                index_min = self._nearest_node(int(listFext[i][5]))
                print(index_min)
                force = listFext[i][4] * 10
                
                if listFext[i][6] == 'Z' or listFext[i][6] == 'Z local' :
                    if listFext[i][5] == index_min[1]:
                        listFlocal[1,0] = force
                    
                    # Cas dans le quel la force est aprés le noeud le plus proche
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
                    
                    # Cas dans le quelle la force est avant le noeud le plus proche
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
                
        print(len(Assembleur_COO_F.data[0])) 
        self.matriceF = coo_matrix(Assembleur_COO_F.data, shape=((self.ele+1)*6, 1)).tocsr()
        # print("\n Matrice des forces extérieurs Fext : \n",self.matriceF, self.matriceF.shape, "\n")
        
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
        
        #print("\n Matrice des forces extérieure Fext au condition limite: \n", self.matriceF_CL,"\n")
        # print("\n Matrice de rigidité K au condition limite: \n", self.matriceK_CL,"\n")
        


    def _equa_deplacement(self):
        ui = spsolve(self.matriceK_CL,self.matriceF_CL)
        j = 0
        for i in range(len(self.matriceU)):
            if self.matriceU[i]:
                self.matriceU[i] = ui[j]
                j += 1
        #print("\n Solution Déplacement : \n", self.matriceU.shape, "\n")
        return self.matriceU


    def _equa_reaction(self):
        self.ri = np.dot(self.matriceK.toarray(), self.matriceU) - self.matriceF.toarray()
        #print("\n Solution Réaction : \n", self.ri, "\n")
        return self.ri


    def reaction(self):
        """ Reaction en kN ou kN.m le long de l'élément"""
        rx = []
        ry = []
        rz = []
        rmx = []
        rmy = []
        rmz = []
        
        for i in range(self.ele):
            n1 = int(self.elementList[i,0])*6
            n2 = int(self.elementList[i,1])*6

            r1 = self.ri[n1-6:n1,0]
            r2 = self.ri[n2-6:n2,0]
            
            if n1 == 6:
                
                rx.append(r1[0]/10**3)
                rx.append(r2[0]/10**3)
                rz.append(r1[2]/10**3)
                rz.append(r2[2]/10**3)
                rmy.append(r1[4]/10**6)
                rmy.append(r2[4]/10**6)
            else:
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

        listMinMax = {"Rx_max": rx_max_coor, "Rx_min": rx_min_coor, 
                      "Rz_max": rz_max_coor, "Rz_min": rz_min_coor, 
                      "RMy_max": rmy_max_coor, "RMy_min": rmy_min_coor}
        return listMinMax

    
    def deplacement(self):
        ux = []
        uz = []
        tetay = []

        for i in range(self.ele):
            n1 = int(self.elementList[i,0])*6
            n2 = int(self.elementList[i,1])*6

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

        for i in range(self.ele):
            n1 = int(self.elementList[i,0])*6
            n2 = int(self.elementList[i,1])*6
            
            u1 = self.matriceU[n1-6:n1,0]
            u2 = self.matriceU[n2-6:n2,0]

            efforinterneN1 = np.dot(self.k_local["eleALL"]["k11_local"],u1) + np.dot(self.k_local["eleALL"]["k12_local"], u2)
            efforinterneN2 = np.dot(self.k_local["eleALL"]["k21_local"],u1) + np.dot(self.k_local["eleALL"]["k22_local"], u2)

            if n1 == 6:
                ni.append(efforinterneN2[0]/10**3)
                ni.append(efforinterneN2[0]/10**3)
                vi.append(efforinterneN2[2]/10**3)
                vi.append(efforinterneN2[2]/10**3)
                myi.append(efforinterneN1[4]/10**6)
                myi.append(efforinterneN2[4]/10**6)
            else:
                ni.append(efforinterneN2[0]/10**3)
                vi.append(efforinterneN2[2]/10**3)
                myi.append(efforinterneN2[4]/10**6)
                
        self.ei_coor = [ni, vi, myi]
        return self.ei_coor


    def effort_interne_max(self):
        """ Retourne les efforts min et max d'une liste d'effort interne """

        nx_max = max(self.ei_coor[0])
        nx_max_index = self.ei_coor[0].index(nx_max)
        nx_max_coor = [self.node_coor[nx_max_index,0], nx_max]

        nx_min = min(self.ei_coor[0])
        nx_min_index = self.ei_coor[0].index(nx_min)
        nx_min_coor = [self.node_coor[nx_min_index,0], nx_min]

        vz_max = max(self.ei_coor[1])
        vz_max_index = self.ei_coor[1].index(vz_max)
        vz_max_coor = [self.node_coor[vz_max_index,0], vz_max]

        vz_min = min(self.ei_coor[1])
        vz_min_index = self.ei_coor[1].index(vz_min)
        vz_min_coor = [self.node_coor[vz_min_index,0], vz_min]

        my_max = max(self.ei_coor[2])
        my_max_index = self.ei_coor[2].index(my_max)
        my_max_coor = [self.node_coor[my_max_index,0], my_max]

        my_min = min(self.ei_coor[2])
        my_min_index = self.ei_coor[2].index(my_min)
        my_min_coor = [self.node_coor[my_min_index,0], my_min]

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
        
        listRMinMax = self.reaction_max()
        self.max_XY_values = (listRMinMax[2][0], listRMinMax[2][1])
        self.min_XY_values = (listRMinMax[3][0], listRMinMax[3][1])
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
        
        listRMinMax = self.reaction_max()
        self.max_XY_values = (listRMinMax[0][0], listRMinMax[0][1])
        self.min_XY_values = (listRMinMax[1][0], listRMinMax[1][1])
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
        
        listUMinMax = self.deplacement_max()
        self.max_XY_values = (listUMinMax[2][0], listUMinMax[2][1])
        self.min_XY_values = (listUMinMax[3][0], listUMinMax[3][1])
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
        
        listEiMinMax = self.effort_interne_max()
        self.max_XY_values = (listEiMinMax[0][0], listEiMinMax[0][1])
        self.min_XY_values = (listEiMinMax[1][0], listEiMinMax[1][1])
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
        self.x_values = self.node_coor
        self.y_values = self.ei_coor[1]
        
        listEiMinMax = self.effort_interne_max()
        self.max_XY_values = (listEiMinMax[2][0], listEiMinMax[2][1])
        self.min_XY_values = (listEiMinMax[3][0], listEiMinMax[3][1])
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
        
        listEiMinMax = self.effort_interne_max()
        self.max_XY_values = (listEiMinMax[4][0], listEiMinMax[4][1])
        self.min_XY_values = (listEiMinMax[5][0], listEiMinMax[5][1])
        self.excent_X = 1
        self.excent_Y = 1
        self.fill_between = True
        self.fill_between_X = self.node_coor.flatten().tolist()
        
        self.save_path = os.path.join(self.SAVE_PATH, "My.png")
        self._base_graphique()

    def show_graphique_My(self, name_combi):
        self._graphique_My(name_combi)
        plt.show()


    def get_supports(self):
        """Retourne la liste des appuis définis.
        """
        return self.list_supports

    def create_support(self, type_appuis: str=("Simple", 'Rotule', 'Encastrement'), pos: int=0, l_appuis: int=0):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            type_appuis (str, optional): type d'appuis à créer. Defaults to ("Simple", 'Rotule', 'Encastrement').
            pos (int, optional): position de l'appuis sur la poutre en mm. Defaults to 0.
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        load = (len(self.list_supports)+1, type_appuis, pos, l_appuis)
        self.list_supports.append(load)
        return load

    
    def create_supports_by_list(self, list_supports: list):
        """Ajoute les charges d'une liste pré-défini dans la liste de chargement

        Args:
            list_supports (list): liste de charge.
        """
        for support in list_supports:
            self.list_supports.append(support)
        return self.list_supports

    def del_support(self, index_load: int):
        """Supprime une charge de l'attribut list_supports par son index

        Args:
            index_load (int): index de la charge à supprimer.
        """
        return self.list_supports.pop(index_load-1)


    def calcul_1D(self):
        """Calcul une poutre MEF 1D

        Args:
            list_loads (list): liste des charges combinées sur la poutre MEF.
        """
        #start = time.time()
        self._matrice_K_1D()
        self._matrice_Fext_1D(self.list_loads)
        self._matrice_U_1D(self.list_supports)

        self._condition_limite()
        self._equa_deplacement()
        self._equa_reaction()

        self.reaction()
        self.effort_interne()
        self.deplacement()
        
        # end = time.time()
        # elapsed = end - start
        # print(f'Temps d\'exécution : {elapsed:.2}s')
    
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
    _list_loads = [[1, '', 'Permanente G', 'Linéique', -10, '0/2000', 'Z'],
                 [0, 'Poids propre', 'Permanente G', 'Linéique', -36, '0/2000', 'Z'],
                 [2, '', 'Neige normale Sn', 'Linéique', -200, '0/2000', 'Z'],
                 [3, '', 'Exploitation Q', 'Linéique', -150, '0/2000', 'Z']]
    chargement = Chargement(pays="Japon")
    chargement.create_load_by_list(_list_loads)
    c1 = Combinaison._from_parent_class(chargement, cat="Cat A : habitation")
    print(c1.list_combination)
    rcombi = "ELU_STR 1.35G + 1.5Sn"
    c1.get_combi_list_load(rcombi)
    long = 2000
    node = int(round(long/100))
    
    b = 140
    h = 600
    a = b*h
    iy = (b*h**3)/12
    iz = (h*b**3)/12
    test = MEF._from_parent_class(c1, long=long,E=11000,A=a, G=350, J=650, Iy=iy, Iz=iz, ele=node, alphaZ=0, alphaY=0, alphaX=0)

    listdeplacement = [[1, "Rotule", 0, 40], [2, "Rotule", 2000, 40]]
    test.create_supports_by_list(listdeplacement)
    
    test.calcul_1D()
    test.graphique("ELU")