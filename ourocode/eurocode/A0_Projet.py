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
        """Retourne les angles de la barre par rapport aux plans XY, XZ et YZ.

        Args:
            v1 (tuple): vecteur 3D représentant la barre (x2-x1, y2-y1, z2-z1)

        Returns:
            dict: Dictionnaire contenant les angles par rapport aux plans XY, XZ et YZ en degrés.
        """
        def _get_local_angle_of_bar(v1: tuple):
            """
            Retourne l'angle dans le plan local de la barre, autour de l'axe longitudinal.
            
            Args:
                v1 (tuple): vecteur 3D représentant la barre (x2-x1, y2-y1, z2-z1)

            Returns:
                float: L'angle de rotation autour de l'axe longitudinal en degrés.
            """
            # Calcul de l'axe longitudinal de la barre (normalisation du vecteur)
            L = np.array(v1) / np.linalg.norm(v1)
            
            # Choisir un vecteur de référence global, par exemple (1, 0, 0) pour l'axe X
            ref_vec = np.array([1, 0, 0])
            
            # Calculer le produit vectoriel pour obtenir un vecteur perpendiculaire à L
            T = np.cross(L, ref_vec)  # Axe transverse
            if np.linalg.norm(T) == 0:
                # Si L est parallèle à l'axe de référence, on choisit un autre vecteur de référence (par exemple, l'axe Y)
                ref_vec = np.array([0, 1, 0])
                T = np.cross(L, ref_vec)
            
            # Normaliser l'axe transverse
            T = T / np.linalg.norm(T)
            
            # Calculer l'angle entre l'axe de référence projeté et l'axe transverse
            angle_local = np.arctan2(np.linalg.norm(np.cross(ref_vec, T)), np.dot(ref_vec, T))
            
            return np.rad2deg(angle_local)
        
        dx, dy, dz = v1

        # Angle dans le plan XY
        angle_xy = np.arctan2(dy, dx)
        
        # Angle dans le plan XZ
        angle_xz = np.arctan2(dz, dx)
        
        # Angle dans le plan YZ
        angle_yz = np.arctan2(dz, dy)
        
        # Convertir les angles en degrés
        angles = {
            "XY": np.rad2deg(angle_xy % (2 * np.pi)),
            "XZ": np.rad2deg(angle_xz % (2 * np.pi)),
            "YZ": np.rad2deg(angle_yz % (2 * np.pi)),
            "local": _get_local_angle_of_bar(v1)
        }
        
        return angles
        # ang1 = np.arctan2(*v1[::-1])
        # return np.rad2deg(ang1 % (2 * np.pi))
    
    
    def __add_element_to_bar(self, bar: dict, elements):
        """On récupère le dernier élément crée et on l'incrémente pour ajouter un nouvel élément"""
        initial_len_of_element_list = len(self.element_list)
        for i in range(len(elements)):
            bar["elements"].append(initial_len_of_element_list + i)


    def __generate_mesh_number(self, bar_lenght: int|float):
        """Calcul le nombre d'élément à créer en fonction de la longueur de la barre

        Args:
            bar_lenght (int | float): Longueur de la barre en mm
        """
        if bar_lenght < 4000:
            # nb_ele = int(mt.ceil(bar_lenght/800))
            nb_ele = 20
        else:
            # nb_ele = int(mt.ceil(bar_lenght/100))
            nb_ele = 20
        return nb_ele

        
    def _create_elements_and_nodes(self, bar: dict):
        """Crée les éléments et les nœuds associés à une barre en utilisant numpy."""
        nb_ele = self.__generate_mesh_number(bar["length"])

        # Si c'est la première barre, on initialise les tableaux
        if not hasattr(self, "node_coor"):
            self.node_coor = np.empty((0, 3), dtype=float)  # Coordonnées des nœuds (X, Y, Z)
            self.element_list = np.empty((0, 2), dtype=int)  # Liste des éléments (nœud 1, nœud 2)

        # Calcul de la direction de la barre dans l'espace 3D (vecteur direction unitaire)
        direction_vector = np.array([bar["x2"] - bar["x1"], bar["y2"] - bar["y1"], bar["z2"] - bar["z1"]])
        unit_vector = direction_vector / bar["length"]  # Vecteur direction unitaire (pour projeter les nœuds)

        # Fonction de génération des coordonnées d'un nœud
        def generate_node(i):
            length_ratio = bar["length"] / nb_ele * i
            node_position = np.array([bar["x1"], bar["y1"], bar["z1"]]) + length_ratio * unit_vector
            return np.round(node_position, 3)  # Arrondi pour éviter les imprécisions numériques

        # Vérifier si un nœud existe déjà et le retourner, sinon l'ajouter
        def get_or_add_node(coor):
            # On vérifie si le nœud existe déjà dans la liste
            idx = np.where((self.node_coor == coor).all(axis=1))[0]
            if len(idx) > 0:
                return idx[0]  # Retourner l'index du nœud existant
            else:
                # Ajouter le nouveau nœud
                self.node_coor = np.vstack([self.node_coor, coor])
                return self.node_coor.shape[0] - 1  # Retourner l'index du nouveau nœud

        # Créer les nœuds et les éléments
        nodes = np.zeros(nb_ele + 1, dtype=int)  # Liste des nœuds
        for i in range(nb_ele + 1):
            node_coor = generate_node(i)
            nodes[i] = get_or_add_node(node_coor)

        # Créer les éléments entre les nœuds
        elements = np.vstack([nodes[:-1]+1, nodes[1:]+1]).T  # Création des éléments en reliant les nœuds
        # Ajouter les éléments à la liste globale
        self.__add_element_to_bar(bar, elements)
        self.element_list = np.vstack([self.element_list, elements])

        print("Elements list: \n", self.element_list)
        print("Node coordinates: \n", self.node_coor)


    def add_bar(self, x1: int, y1: int, z1: int, x2: int, y2: int, z2: int, aire: float):
        """Ajoute une poutre au model MEF

        Args:
            x1 (int): position de départ en x en mm
            y1 (int): position de départ en y en mm
            z1 (int): position de départ en z en mm
            x2 (int): position de fin en x en mm
            y2 (int): position de fin en y en mm
            z2 (int): position de fin en z en mm
            aire (float): aire en mm²
        """
        bar_id = len(self.bar_info) + 1
        v1 = (x2-x1, y2-y1, z2-z1)
        length = mt.sqrt(abs(v1[0])**2 + abs(v1[1])**2 + abs(v1[2])**2)
        angle = self._get_angle_of_bar(v1)
        self.bar_info[bar_id] = {"elements": [], 
                                   "length": length,
                                   "section": aire,
                                   "angle": angle, 
                                   "x1": x1, "y1": y1, "z1": z1, "x2": x2, "y2": y2, "z2": z2,
                                   "relaxation": []}
        print("Poutre crée: ", self.bar_info)
        self._create_elements_and_nodes(self.bar_info[bar_id])


    def add_material_by_class(self, bar_id: int, Iy: float, Iz: float, classe: str=CLASSE_WOOD):
        data_csv_meca = self._data_from_csv("caracteristique_meca_bois.csv")
        material_properties = data_csv_meca.loc[classe]
        E = int(material_properties.loc["E0mean"])
        G = int(material_properties.loc["Gmoy"])
        J = 0
        # ATTENTION: Iy / Iz est différent des Eurocodes ! y est verticale et est donc dans le sens théorique eurocode de Z ! z est donc dans le sens théorique de Y !
        # Il faut donc faire attention au dimenssion de l'élément ! Pour une poutre de section b: 100mm h: 200mm -> Iy = 200 * 100^3 / 12 et Iz = 100 * 200^3 / 12.
        self.bar_info[bar_id]["material"] = {"classe": classe, "E": E, "G": G, "J": J, "Iy": Iz, "Iz": Iy}


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
        # ATTENTION: Iy / Iz est différent des Eurocodes ! y est verticale et est donc dans le sens théorique eurocode de Z ! z est donc dans le sens théorique de Y !
        # Il faut donc faire attention au dimenssion de l'élément ! Pour une poutre de section b: 100mm h: 200mm -> Iy = 200 * 100^3 / 12 et Iz = 100 * 200^3 / 12.
        self.bar_info[bar_id]["material"] = {"classe": "Manuel", "E": E, "G": G, "J": J, "Iy": Iz, "Iz": Iy}

    
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


    def create_support(self, bar_id: int, type_appuis: str=("Simple X", "Simple Y", "Simple Z", "Rotule", "Encastrement"), pos: int=0, l_appuis: int=0):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            bar_id (int): Numéro de la barre sur laquelle positionner l'appuis.
            type_appuis (str, optional): type d'appuis à créer. Defaults to ("Simple", 'Rotule', 'Encastrement').
            pos (int, str, optional): position de l'appuis sur la barre en mm ou chaine de caractère représentant le début de la barre: "start", la fin de la barre: "end", le milieu: "middle ou encore un pourcentage de la longueur de la barre: ex. "40%". Defaults to 0.
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        match pos:
            case "start":
                pos = 0
            case "end":
                pos = self.bar_info[bar_id]["length"]
            case "middle":
                pos = int(round(self.bar_info[bar_id]["length"] / 2, 0))
            case str(x) if '%' in x:
                pos = pos.split("%")[0]
                pos.replace(" ", "")
                pos = int(round(self.bar_info[bar_id]["length"] *  int(pos) / 100, 0))
            case _:
                pass

        self._dict_supports[len(self._dict_supports)+1] = {"N° barre": bar_id, 
                                                              "Type d'appui": type_appuis, 
                                                              "Position de l'appui": pos, 
                                                              "Longueur d'appui": l_appuis}
        return self._dict_supports[len(self._dict_supports)]

    
    def create_supports_by_list(self, list_supports: list):
        """Ajoute les charges d'une liste pré-défini dans la liste de chargement

        Args:
            list_supports (list): liste de charge.
        """
        for support in list_supports:
            print(support)
            self.create_support(*support)
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
    

    def show_graph_loads_and_supports(self, scale_internal_forces: int=20, scale_deplacement: int=10):
        """Affiche le graphique des charges et des appuis du MEF"""

        # Récupération des efforts internes et conversion du locale au globale
        def get_local_to_global_list(list_value: list, internal_force_type: str, list_ele: list, angle=None):
            print(list_value, internal_force_type, list_ele)
            def get_local_to_global(y_local, angle, x, y):
                angle_XY = np.radians(angle["XY"])
                X = x + y_local * -np.sin(angle_XY)
                Y = y + y_local * np.cos(angle_XY)
                return X, Y

            # Initialisation des listes pour stocker les coordonnées X et Y
            X_coords = []
            Y_coords = []

            for ele in list_ele:
                node1 = int(ele[0])-1
                node2 = int(ele[1])-1

                # On récupère l'index de l'élément dans la liste numpy
                print(ele, type(ele), np.where(np.all(self.element_list == ele, axis=1)))
                index_ele = np.where(np.all(self.element_list == ele, axis=1))[0][0]

                x1, y1, z1 = self.node_coor[node1]
                x2, y2, z2 = self.node_coor[node2]
                v1 = (x2-x1, y2-y1, z2-z1)
                if not angle:
                    angle = self._get_angle_of_bar(v1)
                for index_node, node in enumerate((node1, node2)):
                    y_local = list_value[index_ele][internal_force_type][index_node] * scale_internal_forces
                    # print("node: ", node, "if_type: ", internal_force_type, " list_value: ", list_value[index_ele][internal_force_type][index_node])
                    x, y, z = self.node_coor[node]
                    X, Y = get_local_to_global(y_local, angle, x, y)
                    # Ajouter les coordonnées X et Y à la liste correspondante
                    X_coords.append(X)
                    Y_coords.append(Y)
            return X_coords, Y_coords
        
        n_rows, n_cols = 3, 3
        axis_coor = [[row, col] for row in range(n_rows) for col in range(n_cols)]
        
        # Création de la figure et des sous-graphique
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 12))
        
        # Tracé de la structure globale sur le premier sous-graphe
        axs[0, 0].arrow(0, 0, 0, 350, width=15, color="blue")
        axs[0, 0].text(150, -100, "X", ha='right')
        axs[0, 0].arrow(0, 0, 350, 0, width=15, color="red")
        axs[0, 0].text(-50, 150, "Z", ha='right')

        for id_bar, beam in self.bar_info.items():
            x, y = [], []
            ux, uy, uz = [], [], []
            list_ele = [self.element_list[element] for element in beam["elements"]]
            for element in beam["elements"]:
                for i, node in enumerate(self.element_list[element]):
                    coor = self.node_coor[int(node)-1]
                    x.append(coor[0])
                    y.append(coor[1])
                    ux.append(coor[0] + self.u_coor["uX"][int(node)-1]*scale_deplacement)
                    uy.append(coor[0] + self.u_coor["uY"][int(node)-1]*scale_deplacement)
                    uz.append(coor[1] + self.u_coor["uZ"][int(node)-1]*scale_deplacement)

            middle_x = x[int(round((len(x) - 1)/2,0))]
            middle_y = y[int(round((len(y) - 1)/2,0))]
            for subplot in axis_coor:
                axs[subplot[0], subplot[1]].plot(x, y, color="black", linestyle='dashed')
                axs[subplot[0], subplot[1]].text(middle_x+100, middle_y+100, f"B{id_bar}", ha="right", color="black")
                axs[subplot[0], subplot[1]].plot(x[0], y[0], marker='o', mec='gray', mfc="gray")
                axs[subplot[0], subplot[1]].plot(x[-1], y[-1], marker='o', mec='gray', mfc="gray")


            axs[0, 1].plot(ux, uz, color='green')
            axs[0, 2].plot(uy, uz, color='green')
            internal_forces_params = {"Nx": ["orange", [1, 0]],
                                      "Vy": ["maroon", [1, 1]],
                                      "Vz": ["blue", [1, 2]],
                                      "My": ["red", [2, 0]],
                                      "Mz": ["purple", [2, 1]],
                                      }
            
            for internal_f, value in internal_forces_params.items():
                If_coor = get_local_to_global_list(beam["internals forces"], internal_f, list_ele)
                axs[value[1][0], value[1][1]].plot(If_coor[0], If_coor[1], color=value[0])

                if 0 <= beam['angle']["XY"] < 90 or 180 <= beam['angle']["XY"] < 270:
                    axs[value[1][0], value[1][1]].fill_between(x, If_coor[1], y, alpha=0.3, interpolate=True, color=value[0])
                else:
                    axs[value[1][0], value[1][1]].fill_betweenx(y, If_coor[0], x, alpha=0.3, interpolate=True, color=value[0])
                
        
        for key, support in self._dict_supports.items():
            if support["Type d'appui"] == "Rotule":
                support_type = "o"
            elif support["Type d'appui"] == "Encastrement":
                support_type = "s"
            else:
                support_type = "^"
            x = mt.cos(mt.radians(self.bar_info[support["N° barre"]]["angle"]["XY"])) * support["Position de l'appui"] + self.bar_info[support["N° barre"]]["x1"]
            y = mt.sin(mt.radians(self.bar_info[support["N° barre"]]["angle"]["XY"])) * support["Position de l'appui"] + self.bar_info[support["N° barre"]]["y1"]
            for subplot in axis_coor:
                axs[subplot[0], subplot[1]].plot(x, y, marker=support_type, markersize=10, color="red", label=f"A{key} / {support['Type d\'appui']}")
                axs[subplot[0], subplot[1]].text(x+100, y+100, f"A{key}", ha="right", color="red")


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
        
        # # On génère les charges
        # for i, load in enumerate(self.list_loads):
        #     charge = load[4]
        #     parser = parse_position(load[5], charge)
        #     charge = round(-charge,2)
        #     nom = " / ".join([load[1], load[2], load[3], load[6]])
        #     if len(parser[1]) != 2:
        #         unit_load = "daN/m"
        #         plt.plot(parser[0], parser[1], label=nom)
        #         plt.fill_between(parser[0], parser[1], alpha=0.3)
        #     else:
        #         unit_load = "daN"
        #         plt.plot(parser[0], parser[1], marker="X",label=nom)
        #     plt.text(parser[0][1]+1000, charge+2, f'{charge} {unit_load}', ha='right')


        titles = ('Schématisation de la structure et des charges', 'Flèche (u global XZ)', 'Flèche (u global YZ)',
                  'Effort Normal (Nx)', 'Effort Tranchant (Vy)', 'Effort Tranchant (Vz)', 
                  'Moment Fléchissant (My)',  'Moment Fléchissant (Mz)', "")
        
        for index, subplot in enumerate(axis_coor):
            axs[subplot[0], subplot[1]].set_title(titles[index])
            axs[subplot[0], subplot[1]].set_xlabel('Longueur (mm)')
            axs[subplot[0], subplot[1]].set_ylabel('Hauteur (mm)')
            axs[subplot[0], subplot[1]].grid(True)
            if not index:
                axs[subplot[0], subplot[1]].legend()

        plt.tight_layout()
        plt.show()


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