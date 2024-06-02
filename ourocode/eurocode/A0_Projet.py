#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os, sys
import math as mt
import numpy as np
from matplotlib import pyplot as plt
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
        """Créer une classe Batiment héritée de Projet, cette classe défini les dimension du bâtiment

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



class Bar_generator(Projet):
    LIST_SECTION = ["Rectangulaire","Circulaire"]
    CLASSE_WOOD = tuple(Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").index)[2:]
    def __init__(self, *args, **kwargs):
        """Créer une classe héritée de Projet, permettant de générer des barres pour générer tout d'abords des charges 
        puis une modélisation MEF après la combinaison des dites charges à l'EC0.

        Voici les étapes de modélisation:
            1) Créer les barres avec la méthode "add_beam"
            2) Attribuer un matériaux aux barres avec la méthode "add_material"
            3) Créer les appuis avec la méthode "create_support" ou "create_support_by_list"
            4) Créer les charges sur les barres à partir de la classe Chargement dans le module EC0_Combinaison
            5) Créer la classe MEF et lancer le calcul avec la méthode "calcul_1D"
            6) Et voilà, plus cas siroter un thé glacé le temps du calcul ! :)
        """
        super().__init__(*args, **kwargs)
        self.bar_info = {}
        self._dict_supports = {}
        self.bi_connected = 0


    def aire(self, b: int, h: int, section: str=LIST_SECTION):
        if section == self.LIST_SECTION[0]:
            return b * h
        else:
            return mt.pi * (b/2)**2
        
    
    def inertie(self, b: int, h: int, section: str=LIST_SECTION):
        """ Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
            b ou d : Largeur ou diamètre de la poutre en mm
            h : Hauteur de la poutre en mm """
        if section == "Rectangulaire":
            I_y = (b * h**3)/12
            I_z = (h * b**3)/12
        else:
            I_y = (mt.pi * b ** 4) / 64
            I_z = I_y
        return [I_y, I_z]
        
        
    def _get_angle_of_bar(self, v1: tuple):
        """Retourne un angle entre deux points celon le cercle trigo

        Args:
            v1 (tuple): vecteur dans le cercle trigonométrique
        """
        ang1 = np.arctan2(*v1[::-1])
        return np.rad2deg(ang1 % (2 * np.pi))
    
    
    def __add_element_to_beam(self, bar: dict):
        """On récupère le dernier élément crée et on l'incrémente pour ajouter un nouvel élément"""
        elements = self.bar_info[list(self.bar_info)[-1]]["elements"]
        # On récupère le dernière index
        if elements:
            last_element = elements[-1]
        elif len(self.bar_info) > 1:
            elements = self.bar_info[list(self.bar_info)[-2]]["elements"]
            last_element = elements[-1]
        else:
            last_element = -1
        bar["elements"].append(last_element + 1)

        
    def _create_elements_and_nodes(self, bar: dict):
        if bar["length"] < 5000:
            nb_ele = int(mt.ceil(bar["length"]/100))
        else:
            nb_ele = int(mt.ceil(bar["length"]/100))

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
            length = bar["length"] / nb_ele * i
            x = round(bar["x1"] + length * mt.cos(mt.radians(bar["angle"])),3)
            y = round(bar["y1"] + length * mt.sin(mt.radians(bar["angle"])),3)
            if not x:
                x = 0
            if not y:
                y = 0

            if shape[0]>1:
                find_array = np.where((self.node_coor == [x,y]).all(axis=1))[0]
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
                        print(self.element_list[find_array[0]][val], self.bi_connected)
                        nodeList[i-1,1] = self.element_list[find_array[0]][val]
                        # on test si on a déja rattacher un premier noeud de la barre à un noeud existant si oui on incrémente le compteur
                        if j == 1:
                            self.bi_connected += 1
                            bi_connect = True
                    # On incrémente le compteur de soustraction de noeud car on a récupérer un noeud déjà existant
                    j += 1
                    node_coor_to_del.append(i)
                    if not bi_connect and i < nb_ele:
                        self.__add_element_to_beam(bar)
                    continue
            nodeCoor[i,0] = x
            nodeCoor[i,1] = y
            
            if i < nb_ele:
                nodeList[i,0] = shape[0]+i+1 - j
                nodeList[i,1] = shape[0]+i+2 - j
                self.__add_element_to_beam(bar)
        j = 0
        for i_del in node_coor_to_del:
            nodeCoor = np.delete(nodeCoor, i_del - j, 0)
            j += 1
      
        if shape[0]:
            self.node_coor = np.concatenate((self.node_coor, nodeCoor), axis=0)
            self.element_list = np.concatenate((self.element_list, nodeList), axis=0)
        else:
            self.element_list = nodeList
            self.node_coor = nodeCoor
        print(self.element_list, "\n", self.node_coor)


    def add_bar(self, x1: int, y1: int, x2: int, y2: int, aire: float):
        """Ajoute une poutre au model MEF

        Args:
            x1 (int): position de départ en x en mm
            y1 (int): position de départ en y en mm
            x2 (int): position de fin en x en mm
            y2 (int): position de fin en y en mm
            aire (float): aire en mm²
        """
        bar_id = len(self.bar_info) + 1
        v1 = (x2-x1, y2-y1)
        length = mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2)
        angle = self._get_angle_of_bar(v1)
        self.bar_info[bar_id] = {"elements": [], 
                                   "length": length,
                                   "section": aire,
                                   "angle": angle, 
                                   "x1": x1, "y1": y1, "x2": x2, "y2": y2,
                                   "relaxation": []}
        self._create_elements_and_nodes(self.bar_info[bar_id])
        print("Poutre crée: ", self.bar_info)


    def add_material_by_class(self, bar_id: int, Iy: float, Iz: float, classe: str=CLASSE_WOOD):
        data_csv_meca = self._data_from_csv("caracteristique_meca_bois.csv")
        material_properties = data_csv_meca.loc[classe]
        E = int(material_properties.loc["E0mean"])
        G = int(material_properties.loc["Gmoy"])
        J = 0
        self.bar_info[bar_id]["material"] = {"classe": classe, "E": E, "G": G, "J": J, "Iy": Iy, "Iz": Iz}


    def add_material_by_mechanical_properties(self, bar_id: int, E: int, G: float, J: float, Iy: float, Iz: float):
        """Ajoute un matériau à une barre par ces caractéristiques mécaniques.
        
        Args:
            bar_id (int): Le numéro de la barre à laquelle ajouter le matériaux
            E (int): Module de young en MPa, ce E est le E,mean. Il ne faut absolument pas donner le E,mean,fin sous peine de réaliser le calcul EC5 §2.3.2.2 equ2.7 deux fois !
            A (float): Section en mm²
            G (float): Module de cisaillement en MPa
            J (float): Module de torsion
            Iy (float): Inertie quadratique autour de y en mm4
            Iz (float): Inertie quadratique autour de z en mm4
        """
        self.bar_info[bar_id]["material"] = {"classe": "Manuel", "E": E, "G": G, "J": J, "Iy": Iy, "Iz": Iz}

    
    def add_relaxation(self, bar_id: int, position: str=("start", "end"), u: bool=False, v: bool=False, w: bool=False, teta_x: bool=False, teta_y: bool=True, teta_z: bool=True):
        """Ajoute une relaxation dans notre barre soit au début de la barre, soit à la fin.
        Ceci est considéré dans les MEF par une matrice de rigidité spécifique avec les éléments relaché égale à 0.

        Args:
            bar_id (int): numéro de la barre à relacher
            position (str, optional): position du relachement sur la barre soit au début soit à la fin. Defaults to ("start", "end").
            u (bool, optional): relachement de l'axe x local, si oui alors True. Defaults to False.
            v (bool, optional): relachement de l'axe y local, si oui alors True. Defaults to False.
            w (bool, optional): relachement de l'axe z local, si oui alors True. Defaults to False.
            teta_x (bool, optional): relachement de l'axe de rotation x local, si oui alors True. Attention de base cette rotation doit toujours être fixé. Defaults to False.
            teta_y (bool, optional): relachement de l'axe de rotation y local, si oui alors True. Defaults to True.
            teta_z (bool, optional): relachement de l'axe de rotation z local, si oui alors True. Defaults to True.
        """
        self.bar_info[bar_id]["relaxation"].append({"position": position, 
                                               "u": u, 
                                               "v": v, 
                                               "w": w, 
                                               "teta_x": teta_x, 
                                               "teta_y": teta_y, 
                                               "teta_z": teta_z
                                               })


    def create_support(self, bar_id: int, type_appuis: str=("Simple", 'Rotule', 'Encastrement'), pos: int=0, l_appuis: int=0):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            bar_id (int): Numéro de la barre sur laquelle positionner l'appuis.
            type_appuis (str, optional): type d'appuis à créer. Defaults to ("Simple", 'Rotule', 'Encastrement').
            pos (int, optional): position de l'appuis sur la barre en mm. Defaults to 0.
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        self._dict_supports[len(self._dict_supports)+1] = {"N° barre": bar_id, 
                                                              "Type d'appui": type_appuis, 
                                                              "Position de l'appui": pos, 
                                                              "Longueur d'appui": l_appuis}
        return self._dict_supports[len(self._dict_supports)+1]

    
    def create_supports_by_list(self, list_supports: list):
        """Ajoute les charges d'une liste pré-défini dans la liste de chargement

        Args:
            list_supports (list): liste de charge.
        """
        for support in list_supports:
            self._dict_supports[len(self._dict_supports)+1] = {"N° barre": support[0], 
                                                                  "Type d'appui": support[1], 
                                                                  "Position de l'appui": support[2], 
                                                                  "Longueur d'appui": support[3]}
        return self._dict_supports


    def del_support(self, index_support: int):
        """Supprime un appui de l'attribut list_supports par son index

        Args:
            index_load (int): index de l'appui à supprimer.
        """
        return f"L'appui ci-joint à été supprimé: {self._dict_supports.pop(index_support)}"
    

    def get_supports(self):
        """Retourne la liste des appuis définis.
        """
        return self._dict_supports
    

    def show_graph_loads_and_supports(self, scale_internal_forces: int=100, scale_deplacement: int=10):
        """Affiche le graphique des charges et des appuis du MEF"""

        # Récupération des efforts internes et conversion du locale au globale
        def get_local_to_global_list(list_value: list, list_ele: list, angle=None):
            def get_local_to_global(y_local, angle, x, y):
                angle = np.radians(angle)
                X = x + y_local * -np.sin(angle)
                Y = y + y_local * np.cos(angle)
                return X, Y

            # Initialisation des listes pour stocker les coordonnées X et Y
            X_coords = []
            Y_coords = []

            for index_ele, ele in enumerate(list_ele):
                node1 = int(ele[0])-1
                node2 = int(ele[1])-1
                x1, y1 = self.node_coor[node1]
                x2, y2 = self.node_coor[node2]
                v1 = (x2-x1, y2-y1)
                if not angle:
                    angle = self._get_angle_of_bar(v1)
                for node in (node1, node2):
                    y_local = list_value[node] * scale_internal_forces
                    x, y = self.node_coor[node]
                    X, Y = get_local_to_global(y_local, angle, x, y)
                    # Ajouter les coordonnées X et Y à la liste correspondante
                    X_coords.append(X)
                    Y_coords.append(Y)
            return X_coords, Y_coords
        
        n_rows, n_cols = 2, 3
        axis_coor = [[row, col] for row in range(n_rows) for col in range(n_cols)]
        print(axis_coor)
        # Calcul des efforts internes
        ei_coor = self.effort_interne()
        
        # Création de la figure et des sous-graphique
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(21, 10))
        
        # Tracé de la structure globale sur le premier sous-graphe
        axs[0, 0].arrow(0, 0, 0, 350, width=15, color="blue")
        axs[0, 0].text(150, -100, "X", ha='right')
        axs[0, 0].arrow(0, 0, 350, 0, width=15, color="red")
        axs[0, 0].text(-50, 150, "Z", ha='right')

        for key, beam in self.bar_info.items():
            x, y = [], []
            ux, uz = [], []
            list_ele = [self.element_list[element] for element in beam["elements"]]
            for element in beam["elements"]:
                for i, node in enumerate(self.element_list[element]):
                    coor = self.node_coor[int(node)-1]
                    x.append(coor[0])
                    y.append(coor[1])
                    ux.append(coor[0] + self.u_coor["uX"][int(node)-1]*scale_deplacement)
                    uz.append(coor[1] + self.u_coor["uZ"][int(node)-1]*scale_deplacement)

            middle_x = x[int(round((len(x) - 1)/2,0))]
            middle_y = y[int(round((len(y) - 1)/2,0))]
            for subplot in axis_coor:
                axs[subplot[0], subplot[1]].plot(x, y, color="black", linestyle='dashed')
                axs[subplot[0], subplot[1]].text(middle_x+100, middle_y+100, f"B{key}", ha="right", color="black")
                axs[subplot[0], subplot[1]].plot(x[0], y[0], marker='o', mec='gray', mfc="gray")
                axs[subplot[0], subplot[1]].plot(x[-1], y[-1], marker='o', mec='gray', mfc="gray")
            axs[0, 1].plot(ux, uz, color='green')

            
            Nx_coor = get_local_to_global_list(ei_coor["Nx"], list_ele)
            axs[0, 2].plot(Nx_coor[0], Nx_coor[1], color='orange')
            
            Vz_coor = get_local_to_global_list(ei_coor["Vz"], list_ele)
            axs[1, 0].plot(Vz_coor[0], Vz_coor[1], color='blue')

            My_coor = get_local_to_global_list(ei_coor["My"], list_ele)
            axs[1, 1].plot(My_coor[0], My_coor[1], color='red')

            Mz_coor = get_local_to_global_list(ei_coor["Mz"], list_ele)
            axs[1, 2].plot(Mz_coor[0], Mz_coor[1], color='purple')

            # U_coor = get_local_to_global_list(self.u_coor["uX"], list_ele, 0)
            # axs[1, 2].plot(U_coor[0], U_coor[1], label="Flèche (u global)", color='green')
        
        for key, support in self._dict_supports.items():
            if support["Type d'appui"] == "Rotule":
                support_type = "o"
            elif support["Type d'appui"] == "Encastrement":
                support_type = "s"
            else:
                support_type = "^"
            x = mt.cos(mt.radians(self.bar_info[support["N° barre"]]["angle"])) * support["Position de l'appui"] + self.bar_info[support["N° barre"]]["x1"]
            y = mt.sin(mt.radians(self.bar_info[support["N° barre"]]["angle"])) * support["Position de l'appui"] + self.bar_info[support["N° barre"]]["y1"]
            for subplot in axis_coor:
                axs[subplot[0], subplot[1]].plot(x, y, marker=support_type, markersize=10, color="red", label=f"Appui {key} / {support['Type d\'appui']}")
                axs[subplot[0], subplot[1]].text(x+100, y+100, f"A{key}", ha="right", color="red")

        titles = ('Schématisation de la structure et des charges', 'Flèche (u global XZ)', 
                  'Effort Normal (Nx)', 'Effort Tranchant (Vz)', 'Moment Fléchissant (My)', 
                  'Moment Fléchissant (Mz)')
        
        for index, subplot in enumerate(axis_coor):
            axs[subplot[0], subplot[1]].set_title(titles[index])
            axs[subplot[0], subplot[1]].set_xlabel('Longueur (mm)')
            axs[subplot[0], subplot[1]].set_ylabel('Hauteur (mm)')
            axs[subplot[0], subplot[1]].legend()
            axs[subplot[0], subplot[1]].grid(True)

        plt.tight_layout()
        plt.show()
        

    # def show_graph_loads_and_supports(self):
    #     """Affiche le graphique des charges et des appuis du MEF
    #     """
    #     # Conversion des positions en valeurs numériques
    #     def parse_position(position, charge):
    #         position = str(position)
    #         charge = -charge
    #         parts = position.split("/")
    #         if len(parts) == 1:
    #             return [int(parts[0])]*2, [0,charge]
    #         elif len(parts) == 2:
    #             return [pos for pos in range(int(parts[0]), int(parts[1]), 1)], [charge for _ in range(int(parts[0]), int(parts[1]), 1)]
    #         else:
    #             return None

    #     # Création du graphique
    #     plt.figure(figsize=(12, 4))

    #     plt.arrow(0, 0, 0, 350, width=15, color="blue")
    #     plt.text(150, -100, "X", ha='right')
    #     plt.arrow(0, 0, 350, 0, width=15, color="red")
    #     plt.text(-50, 150, "Z", ha='right')

    #     # On dessine les barres
    #     for key, beam in self.bar_info.items():
    #         x, y = [], []
    #         for element in beam["elements"]:
    #             for i, node in enumerate(self.element_list[element]):
    #                 coor = self.node_coor[int(node)-1]
    #                 x.append(coor[0])
    #                 y.append(coor[1])
    #         plt.plot(x, y, label=f"Barre N°{key}", marker="o")
    #         plt.plot(x[0], y[0], marker='o', mec='gray', mfc="gray")
    #         plt.plot(x[-1], y[-1], marker='o', mec='gray', mfc="gray")

    #     # # On génère les charges
    #     # for i, load in enumerate(self.list_loads):
    #     #     charge = load[4]
    #     #     parser = parse_position(load[5], charge)
    #     #     charge = round(-charge,2)
    #     #     nom = " / ".join([load[1], load[2], load[3], load[6]])
    #     #     if len(parser[1]) != 2:
    #     #         unit_load = "daN/m"
    #     #         plt.plot(parser[0], parser[1], label=nom)
    #     #         plt.fill_between(parser[0], parser[1], alpha=0.3)
    #     #     else:
    #     #         unit_load = "daN"
    #     #         plt.plot(parser[0], parser[1], marker="X",label=nom)
    #     #     plt.text(parser[0][1]+1000, charge+2, f'{charge} {unit_load}', ha='right')

    #     # On génère les appuis
    #     for key, support in self._dict_supports.items():
    #         if support["Type d'appui"] == "Rotule":
    #             support_type = "o"
    #         elif support["Type d'appui"] == "Encastrement":
    #             support_type = "s"
    #         else:
    #             support_type = "^"
    #         x = mt.cos(mt.radians(self.bar_info[support["N° barre"]]["angle"])) * support["Position de l'appui"] + self.bar_info[support["N° barre"]]["x1"]
    #         y = mt.sin(mt.radians(self.bar_info[support["N° barre"]]["angle"])) * support["Position de l'appui"] + self.bar_info[support["N° barre"]]["y1"]
    #         plt.plot(x, y, marker=support_type, markersize=10, color="red", label=f"Appui {key} / {support["Type d'appui"]}")

    #     plt.title('Schématisation de la structure et des charges')
    #     plt.xlabel('Longueur (mm)')
    #     plt.ylabel('Hauteur (mm')
    #     plt.legend()
    #     plt.grid(True)
    #     plt.show()


if __name__ == "__main__":
    # projet = Projet(num_project="6006.0", commmentaire="c'est mon premier projet")
    # building = Batiment(h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=15)
    # print(building.h_bat)
    beam_gen = Bar_generator()
    beam_gen.add_bar(0,0,2500,4000,20000)
    beam_gen.add_material_by_class(1, 60000, 60000, "C24")
    print(beam_gen.bar_info[1])

    listdeplacement = [[1, "Rotule", 0, 0], [1, "Rotule", beam_gen.bar_info[1]["length"], 0]]
    beam_gen.create_supports_by_list(listdeplacement)

    beam_gen.show_graph_loads_and_supports()