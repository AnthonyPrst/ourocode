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



class Beam_generator(Projet):
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
    
    
    def __add_element_to_beam(self, beam: dict):
        """On récupère le dernier élément crée et on l'incrémente pour ajouter un nouvel élément"""
        elements = self.beams[str(list(self.beams)[-1])]["elements"]
        # On récupère le dernière index
        if elements:
            last_element = elements[-1]
        elif len(self.beams) > 1:
            elements = self.beams[str(list(self.beams)[-2])]["elements"]
            last_element = elements[-1]
        else:
            last_element = -1
        beam["elements"].append(last_element + 1)

        
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
                        self.__add_element_to_beam(beam)
                    continue
            nodeCoor[i,0] = x
            nodeCoor[i,1] = y
            
            if i < nb_ele:
                nodeList[i,0] = shape[0]+i+1 - j
                nodeList[i,1] = shape[0]+i+2 - j
                self.__add_element_to_beam(beam)
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
        self.beams[str(i_beam)] = {"elements": [], 
                                   "length": length,
                                   "angle": angle, 
                                   "x1": x1, "y1": y1, "x2": x2, "y2": y2}
        self._create_elements_and_nodes(self.beams[str(i_beam)])
        print("Poutre crée: ", self.beams)


    def add_material_by_class(self, beam_number: int, classe: str=CLASSE_WOOD):
        pass

    def add_material_by_mechanical_properties(self, beam_number: int, E: int, G: float, J: float, Iy: float, Iz: float):
        """Ajoute un matériau à une barre par ces caractéristiques mécaniques.
        
        Args:
            beam_number (int): Le numéro de la barre à laquelle ajouter le matériaux
            E (int): Module de young en MPa, ce E est le E,mean. Il ne faut absolument pas donner le E,mean,fin sous peine de réaliser le calcul EC5 §2.3.2.2 equ2.7 deux fois !
            A (float): Section en mm²
            G (float): Module de cisaillement en MPa
            J (float): Module de torsion
            Iy (float): Inertie quadratique autour de y en mm4
            Iz (float): Inertie quadratique autour de z en mm4
        """
        self.beams[str(beam_number)]["material"] = {"classe": "Manuel", "E": E, "G": G, "J": J, "Iy": Iy, "Iz": Iz}


    def create_support(self, beam_number: int, type_appuis: str=("Simple", 'Rotule', 'Encastrement'), pos: int=0, l_appuis: int=0):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            type_appuis (str, optional): type d'appuis à créer. Defaults to ("Simple", 'Rotule', 'Encastrement').
            pos (int, optional): position de l'appuis sur la barre en mm. Defaults to 0.
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        self.list_supports[str(len(self.list_supports)+1)] = {"N° barre": beam_number, 
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
            self.list_supports[str(len(self.list_supports)+1)] = {"N° barre": support[0], 
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
    

    def get_supports(self):
        """Retourne la liste des appuis définis.
        """
        return self.list_supports
    

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

        # On dessine les barres
        for key, beam in self.beams.items():
            x, y = [], []
            for element in beam["elements"]:
                for i, node in enumerate(self.element_list[element]):
                    coor = self.node_coor[int(node)-1]
                    x.append(coor[0])
                    y.append(coor[1])
            plt.plot(x, y, label=f"Barre N°{key}")
            plt.plot(x[0], y[0], marker='o', mec='gray', mfc="gray")
            plt.plot(x[-1], y[-1], marker='o', mec='gray', mfc="gray")

        # On génère les charges
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

        # On génère les appuis
        for key, support in self.list_supports.items():
            if support["Type d'appui"] == "Rotule":
                support_type = "o"
            elif support["Type d'appui"] == "Encastrement":
                support_type = "s"
            else:
                support_type = "^"
            x = mt.cos(mt.radians(self.beams[str(support["N° barre"])]["angle"])) * support["Position de l'appui"] + self.beams[str(support["N° barre"])]["x1"]
            y = mt.sin(mt.radians(self.beams[str(support["N° barre"])]["angle"])) * support["Position de l'appui"] + self.beams[str(support["N° barre"])]["y1"]
            plt.plot(x, y, marker=support_type, markersize=10, color="red", label=f"Appui {key} / {support["Type d'appui"]}")

        plt.title('Schématisation de la structure et des charges')
        plt.xlabel('Longueur (mm)')
        plt.ylabel('Hauteur (mm')
        plt.legend()
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    # projet = Projet(num_project="6006.0", commmentaire="c'est mon premier projet")
    # building = Batiment(h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=15)
    # print(building.h_bat)
    beam_gen = Beam_generator()
    beam_gen.add_beam(0,0,0,5000)
    beam_gen.add_beam(0,5000,5000,5000)
    beam_gen.add_beam(5000,0,5000,5000)

    listdeplacement = [[1, "Rotule", 0, 0], [1, "Rotule", 3000, 0]]
    beam_gen.create_supports_by_list(listdeplacement)

    beam_gen.show_graph_loads_and_supports()