#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os, sys
import math as mt
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image

from PySide6.QtWidgets import QFileDialog, QMessageBox
import forallpeople as si
from handcalcs.decorator import handcalc
from Pynite import FEModel3D
from Pynite.Rendering import Renderer

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# from eurocode.objet import Objet

from ourocode.eurocode.objet import Objet


class Projet(Objet):
    JUPYTER_DISPLAY = False
    DICO_COMBI_ACTION = {
        "Permanente G": "G",
        "Exploitation Q": "Q",
        "Neige normale Sn": "Sn",
        "Vent pression W+": "W+",
        "Vent dépression W-": "W-",
        "Neige accidentelle Sx": "Sx",
        "Sismique Ae": "Ae",
    }

    def __init__(
        self,
        ingenieur: str = None,
        num_project: str = None,
        name: str = None,
        adresse: str = None,
        code_INSEE: int = None,
        pays: str = "France",
        alt: si.m = 0,
        **kwargs,
    ):
        """Créer une classe Projet hérité de la classe Objet du fichier objet.py. Cette classe définit le projet, d'ou découle l'ensemble des objets du catalogue.

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
        return "Créer une classe Projet hérité de la classe Objet du fichier objet.py. Cette classe définit le projet, d'ou découle l'ensemble des objets du catalogue."

    def __repr__(self) -> str:
        return super().__repr__()


class Batiment(Projet):
    def __init__(
        self,
        h_bat: float,
        d_bat: float,
        b_bat: float,
        alpha_toit: float,
        alpha_toit2: float = 0,
        *args,
        **kwargs,
    ):
        """Créer une classe Batiment héritée de Projet, cette classe définit les dimension du bâtiment

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
        self.b_bat = b_bat  # coté perpendiculaire au vent longpant
        self.alpha_toit = alpha_toit
        self.alpha_toit2 = alpha_toit2

    def height_ref(self):
        pass  # 7.2.1 bâtiment de grande hauteur


class Model_generator(Projet):
    ACTION = (
        "Permanente G",
        "Exploitation Q",
        "Neige normale Sn",
        "Vent pression W+",
        "Vent dépression W-",
        "Neige accidentelle Sx",
        "Sismique Ae",
    )
    LIST_SECTION = ["Rectangulaire", "Circulaire"]
    CLASSE_WOOD = tuple(
        Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").index
    )[2:]

    def __init__(self, *args, **kwargs):
        """Créer une classe héritée de Projet, permettant de générer des barres pour générer tout d'abords des charges
        puis une modélisation MEF après la combinaison des dites charges à l'EC0.

        Voici les étapes de modélisation:
            1) Créer les sections avec la méthode "add_section"
            2) Créer les matériaux avec la méthode "add_material"
            3) Créer les barres avec la méthode "add_member"
            4) Créer les appuis avec la méthode "create_support"
            5) Créer les charges sur les barres avec les méthodes "create_dist_load" et "create_point_load"
            7) Transmettre cette classe dans les argument de la classe Combinaison du module EC0_Combinaison pour le calcul des combinaisons.
            8) Exécuté la méthode "generate_model" avec la classe Combinaison.
            9) Récupérer les efforts internes , les déformations et afficher les graphiques associés avec les méthodes correspondante.
        """
        super().__init__(*args, **kwargs)
        self._data = {
            "nodes": {},
            "sections": {},
            "materials": {},
            "members": {},
            "supports": {},
            "loads": {},
        }
        self._model = None

    def get_all_data(self) -> dict:
        """Retourne l'ensemble des données du model"""
        return self._data

    def export_data(self):
        """Export les données du model au format JSON

        Args:
            filename (str, optional): nom du fichier à créer. Defaults to "model.json".
        """
        import json
        from PySide6.QtWidgets import QFileDialog, QApplication
        from PySide6.QtCore import Qt

        # QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
        app = QApplication(sys.argv)
        filename = QFileDialog.getSaveFileName(
            filter="'JSON' (*.json)",
            selectedFilter=".json",
        )[0]
        with open(filename, "w") as f:
            json.dump(self._data, f, indent=4)

    def show_sign_convention(self):
        """Affiche l'image de la convention de signe d'une barre'"""
        file = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "sign_convention.png"
        )
        image = Image.open(file)
        image.show()

    ################## Noeud ##################

    def add_node(self, X: si.mm, Y: si.mm, Z: si.mm, comment: str = None) -> str:
        """Ajoute un noeud au model MEF

        Args:
            X (int): position en X global en mm
            Y (int): position en Y global en mm
            Z (int): position en Z global en mm
        """
        node_id = "N" + str(len(self._data["nodes"]) + 1)
        self._data["nodes"][node_id] = {
            "X": X * si.mm,
            "Y": Y * si.mm,
            "Z": Z * si.mm,
            "Commentaire": comment,
        }
        return node_id

    def _add_node_to_model(self, node_id: str):
        """Ajoute un noeud au model MEF

        Args:
            node_id (str): id du noeud à ajouter
        """
        node = self._data["nodes"][node_id]
        self._model.add_node(
            node_id,
            node["X"].value * 10**3,
            node["Y"].value * 10**3,
            node["Z"].value * 10**3,
        )
        return node_id

    def _add_nodes_to_model(self):
        """Ajoute tous les noeuds au model MEF"""
        for node_id in self._data["nodes"].keys():
            self._add_node_to_model(node_id)
        return self._model.nodes

    def get_node(self, node_id: str) -> dict:
        """Retourne les coordonnée du noeud

        Args:
            node_id (str): id du noeud à récupérer

        Returns:
            dict: dictionnaire contenant les coordonnées du noeud
        """
        return self._data["nodes"][node_id]

    def get_all_nodes(self) -> dict:
        return self._data["nodes"]

    ################## Barre ##################

    def _get_angle_of_bar(self, vector: tuple):
        """Retourne les angles de la barre par rapport aux plans XY, XZ et YZ.

        Args:
            vector (tuple): vecteur 3D représentant la barre (x2-x1, y2-y1, z2-z1)

        Returns:
            dict: Dictionnaire contenant les angles par rapport aux plans XY, XZ et YZ en degrés.
        """

        def _get_local_angle_of_bar(vector: tuple):
            """
            Retourne l'angle dans le plan local de la barre, autour de l'axe longitudinal.

            Args:
                vector (tuple): vecteur 3D représentant la barre (x2-x1, y2-y1, z2-z1)

            Returns:
                float: L'angle de rotation autour de l'axe longitudinal en degrés.
            """
            # Calcul de l'axe longitudinal de la barre (normalisation du vecteur)
            L = np.array(vector) / np.linalg.norm(vector)

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
            angle_local = np.arctan2(
                np.linalg.norm(np.cross(ref_vec, T)), np.dot(ref_vec, T)
            )

            return np.rad2deg(angle_local)

        dx, dy, dz = vector

        # Angle dans le plan XY
        angle_xy = np.arctan2(dy.value(), dx.value())

        # Angle dans le plan XZ
        angle_xz = np.arctan2(dz.value(), dx.value())

        # Angle dans le plan YZ
        angle_yz = np.arctan2(dz.value(), dy.value())

        # Convertir les angles en degrés
        angles = {
            "XY": np.rad2deg(angle_xy % (2 * np.pi)),
            "XZ": np.rad2deg(angle_xz % (2 * np.pi)),
            "YZ": np.rad2deg(angle_yz % (2 * np.pi)),
            "local": _get_local_angle_of_bar(vector),
        }
        # ang1 = np.arctan2(*v1[::-1])
        # return np.rad2deg(ang1 % (2 * np.pi))
        return angles

    def add_member(
        self,
        node1: str,
        node2: str,
        material: str,
        section: str,
        rotation: float = 0,
        tension_only: bool = ("False", "True"),
        compression_only: bool = ("False", "True"),
        name: str = None,
        comment: str = None,
    ):
        """Ajoute une poutre au model MEF

        Args:
            node1 (str): id du noeud 1
            node2 (str): id du noeud 2
            material (str): id du matériaux
            section (str): id de la section
            rotation (float, optional): angle de rotation de la poutre en °. Defaults to 0.
            tension_only (bool, optional): si True, la poutre ne peut que subir des efforts de traction. Defaults to False.
            compression_only (bool, optional): si True, la poutre ne peut que subir des efforts de compression. Defaults to False.
            name (str, optional): Nom de la membrure, si vous remplissez cet argument alors c'est le nom de la membrure,
                                  sinon la fonction en génère un automatique.
            comment (str, optional): commentaire sur la poutre. Defaults to None.
        """
        if name:
            member_id = name
        else:
            member_id = "M" + str(len(self._data["members"]) + 1)
        node_coor_1 = self.get_node(node1)
        node_coor_2 = self.get_node(node2)
        x1, y1, z1 = (
            node_coor_1["X"].value * 10**3,
            node_coor_1["Y"].value * 10**3,
            node_coor_1["Z"].value * 10**3,
        )
        x2, y2, z2 = (
            node_coor_2["X"].value * 10**3,
            node_coor_2["Y"].value * 10**3,
            node_coor_2["Z"].value * 10**3,
        )
        vector = (x2 - x1, y2 - y1, z2 - z1)
        length = (
            mt.sqrt(abs(vector[0]) ** 2 + abs(vector[1]) ** 2 + abs(vector[2]) ** 2)
            * si.mm
        )

        self._data["members"][member_id] = {
            "Noeuds": [node1, node2],
            "Longueur": length,
            "Section": section,
            "Matériaux": material,
            "Rotation": rotation,
            "Relaxation": {"start": None, "end": None},
            "Commentaire": comment,
        }
        if tension_only:
            self._data["members"][member_id]["Tension uniquement"] = True
        elif compression_only:
            self._data["members"][member_id]["Compression uniquement"] = True

        return member_id

    def _add_member_to_model(self, member_id: str):
        """Ajoute une poutre au model MEF

        Args:
            member_id (str): id de la poutre à ajouter
        """
        node1, node2 = self._data["members"][member_id]["Noeuds"]
        tension_only = False
        compression_only = False
        if self._data["members"][member_id].get("Tension uniquement"):
            tension_only = True
        elif self._data["members"][member_id].get("Compression uniquement"):
            compression_only = True
        self._model.add_member(
            member_id,
            node1,
            node2,
            self._data["members"][member_id]["Matériaux"],
            self._data["members"][member_id]["Section"],
            self._data["members"][member_id]["Rotation"],
            tension_only,
            compression_only,
        )

        # Ajout des releases
        list_releases = []
        has_release = False
        for pos, release in self._data["members"][member_id]["Relaxation"].items():
            if release:
                has_release = True
                list_releases.append([val for val in release.values()])
            else:
                list_releases.append([False] * 6)

        if has_release:
            list_releases = list_releases[0] + list_releases[1]
            print(list_releases)
            self._model.def_releases(
                member_id,
                *list_releases,
            )
        return member_id

    def _add_members_to_model(self):
        """Ajoute toutes les poutres au model MEF"""
        for member_id in self._data["members"].keys():
            self._add_member_to_model(member_id)
        return self._model.members

    def get_member(self, member_id: str) -> dict:
        """Retourne une membrure par son id

        Args:
            member_id (str): id de la membrure à récupérer
        """
        return self._data["members"][member_id]
    
    def get_member_length(self, member_id: str):
        """Retourne la longueur du membrure par son id

        Args:
            member_id (str): id de la membrure à récupérer
        """
        return self._data["members"][member_id]["Longueur"]

    def get_all_members(self) -> dict:
        return self._data["members"]

    ################## Matériaux ##################

    def add_material_by_class(self, classe: str = CLASSE_WOOD) -> str:
        if not self._data["materials"].get(classe):
            data_csv_meca = self._data_from_csv("caracteristique_meca_bois.csv")
            material_properties = data_csv_meca.loc[classe]
            E = int(material_properties.loc["E0mean"])
            G = int(material_properties.loc["Gmoy"])
            nu = (E / (2 * G)) - 1  # Coefficient de Poisson
            self._data["materials"][classe] = {
                "classe": classe,
                "E": E * si.MPa,
                "G": G * si.MPa,
                "nu": nu,
                "rho": int(material_properties.loc["rhomean"]) * si.kg / si.m**3,
            }
        return classe

    def add_material_by_mechanical_properties(
        self,
        name: str,
        E: si.MPa,
        G: si.MPa,
        nu: float,
        rho: float,
    ):
        """Ajoute un matériau à une barre par ces caractéristiques mécaniques.

        Args:
            name (str): nom du matériau
            E (int): Module de young en MPa, ce E est le E,mean. Il ne faut absolument pas donner le E,mean,fin sous peine de réaliser le calcul EC5 §2.3.2.2 equ2.7 deux fois !
            G (float): Module de cisaillement en MPa
            nu (float): Coefficient de Poisson
            rho (float): Masse volumique en kg/m3
        """
        self._data["materials"][name] = {
            "classe": "Manuel",
            "E": E * si.MPa,
            "G": G * si.MPa,
            "nu": nu,
            "rho": rho * si.kg / si.m**3,
        }
        return name

    def _add_material_to_model(self, material_id: str):
        """Ajoute un matériau au model MEF

        Args:
            material_id (str): id du matériau à ajouter
        """
        material = self._data["materials"][material_id]
        self._model.add_material(
            material_id,
            material["E"].value * 10**-6,
            material["G"].value * 10**-6,
            material["nu"],
            material["rho"].value,
        )
        return material_id

    def _add_materials_to_model(self):
        for mat_id in self._data["materials"].keys():
            self._add_material_to_model(mat_id)
        return self._model.materials

    def get_material(self, material_id: str) -> dict:
        return self._data["materials"][material_id]

    def get_all_materials(self) -> dict:
        return self._data["materials"]

    ################## Section ##################

    def aire(self, b: si.mm, h: si.mm, section: str = LIST_SECTION):
        b = b * si.mm
        h = h * si.mm
        if section == self.LIST_SECTION[0]:
            return b * h
        else:
            return mt.pi * (b / 2) ** 2

    def inertie(self, b: si.mm, h: si.mm, section: str = LIST_SECTION):
        """Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
        b ou d : Largeur ou diamètre de la poutre en mm
        h : Hauteur de la poutre en mm"""
        b = b * si.mm
        h = h * si.mm
        if section == "Rectangulaire":
            I_z = (b * h**3) / 12
            I_y = (h * b**3) / 12
        else:
            I_y = (mt.pi * b**4) / 64
            I_z = I_y
        return {"Iy": I_y, "Iz": I_z}

    def add_section(self, b: si.mm, h: si.mm, J: si.mm**4, section: str = LIST_SECTION):
        """Ajoute une section à la liste de section

        Args:
            b (int): largeur de la section en mm
            h (int): hauteur de la section en mm
            J (float): Module de torsion en mm4
            section (str): type de section "Circulaire" ou "Rectangulaire".
        """
        if section not in self.LIST_SECTION:
            raise ValueError(
                f"Le type de section {section} n'est pas reconnu. Les types de sections disponibles sont: {self.LIST_SECTION}"
            )
        match section:
            case "Rectangulaire":
                name = "".join(["R", str(b), "X", str(h)])
            case "Circulaire":
                name = "".join(["C", str(b)])
        inertie = self.inertie(b, h, section)
        self._data["sections"][name] = {
            "Section": section,
            "b": b * si.mm,
            "h": h * si.mm,
            "Aire": self.aire(b, h, section),
            "Iy": inertie["Iy"],
            "Iz": inertie["Iz"],
            "J": J * si.mm**4,
        }
        return name

    def add_section_by_property(
        self, name: str, aire: si.mm**2, Iy: si.mm**4, Iz: si.mm**4, J: si.mm**4
    ):
        """Ajoute une section à la liste de section

        Args:
            name (str): nom de la section
            aire (float): aire de la section en mm²

            ATTENTION: Iy est la petite inertie et Iz est la grande inertie !
            Iy (float): Inertie quadratique autour de y en mm4
            Iz (float): Inertie quadratique autour de z en mm4

            J (float): Module de torsion en mm4
        """
        self._data["sections"][name] = {
            "Section": "Manuel",
            "Aire": aire * si.mm**2,
            "Iy": Iy * si.mm**4,
            "Iz": Iz * si.mm**4,
            "J": J * si.mm**4,
        }
        return self._data["sections"][name]

    def _add_section_to_model(self, section_id: str):
        """Ajoute une section au model MEF

        Args:
            section_id (str): id de la section à ajouter
        """
        section = self._data["sections"][section_id]
        self._model.add_section(
            section_id,
            section["Aire"].value * 10**6,
            section["Iy"].value * 10**12,
            section["Iz"].value * 10**12,
            section["J"].value * 10**12,
        )
        return section_id

    def _add_sections_to_model(self):
        for section_id in self._data["sections"].keys():
            self._add_section_to_model(section_id)
        return self._model.sections

    def get_section(self, section_id: str) -> dict:
        return self._data["sections"][section_id]

    def get_all_sections(self) -> dict:
        return self._data["sections"]

    ################## Relachement ##################

    def add_release(
        self,
        member_id: str,
        position: str = ("start", "end"),
        u: bool = ("False", "True"),
        v: bool = ("False", "True"),
        w: bool = ("False", "True"),
        teta_x: bool = ("False", "True"),
        teta_y: bool = ("False", "True"),
        teta_z: bool = ("False", "True"),
    ):
        """Ajoute une relaxation sur une membrure soit au début, soit à la fin.
        Ceci est considéré dans les MEF par une matrice de rigidité spécifique avec les éléments relaché égale à 0.

        Args:
            member_id (int): numéro de la membrure à relacher
            position (str, optional): position du relachement sur la barre soit au début soit à la fin. Defaults to ("start", "end").
            u (bool, optional): relachement de l'axe x local, si oui alors True.
            v (bool, optional): relachement de l'axe y local, si oui alors True.
            w (bool, optional): relachement de l'axe z local, si oui alors True.
            teta_x (bool, optional): relachement de l'axe de rotation x local, si oui alors True. Attention de base cette rotation doit toujours être fixé.
            teta_y (bool, optional): relachement de l'axe de rotation y local, si oui alors True.
            teta_z (bool, optional): relachement de l'axe de rotation z local, si oui alors True.
        """
        self._data["members"][member_id]["Relaxation"][position] = {
            "u": u,
            "v": v,
            "w": w,
            "teta_x": teta_x,
            "teta_y": teta_y,
            "teta_z": teta_z,
        }
        return self._data["members"][member_id]["Relaxation"][position]

    ################## Appuis ##################

    def add_support(
        self,
        node_id: str,
        DX: bool = ("True", "False"),
        DY: bool = ("True", "False"),
        DZ: bool = ("True", "False"),
        RX: bool = ("True", "False"),
        RY: bool = ("True", "False"),
        RZ: bool = ("True", "False"),
        l_appuis: int = 0,
    ):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            node_id (int): Numéro du noeud sur lequel positionner l'appuis.
            DX (bool, optional): Blocage en translation de l'axe X global, si oui alors True.
            DY (bool, optional): Blocage en translation de l'axe Y global, si oui alors True.
            DZ (bool, optional): Blocage en translation de l'axe Z global, si oui alors True.
            RX (bool, optional): Blocage en rotation de l'axe X global, si oui alors True.
            RY (bool, optional): Blocage en translation de l'axe X global, si oui alors True.
            RZ (bool, optional): Blocage en translation de l'axe X global, si oui alors True.
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        support_id = "S" + str(len(self._data["supports"]) + 1)
        self._data["supports"][support_id] = {
            "Noeud": node_id,
            "DX": DX,
            "DY": DY,
            "DZ": DZ,
            "RX": RX,
            "RY": RY,
            "RZ": RZ,
            "Longueur d'appui": l_appuis,
        }
        return self._data["supports"][support_id]

    def create_supports_by_list(self, list_supports: list):
        """Ajoute les charges d'une liste pré-définit dans la liste de chargement

        Args:
            list_supports (list): liste de charge.
        """
        for support in list_supports:
            self.add_support(*support)

    def del_support(self, support_id: str):
        """Supprime un appui de l'attribut list_supports par son index

        Args:
            support_id (int): id de l'appuis à supprimer.
        """
        return f"L'appui à été supprimé: {self._data["supports"].pop(support_id)}"

    def _add_support_to_model(self, support_id: str):
        """Ajoute un appui au model MEF

        Args:
            support_id (str): id de l'appui à ajouter
        """
        support = self._data["supports"][support_id]
        self._model.def_support(
            support["Noeud"],
            support["DX"],
            support["DY"],
            support["DZ"],
            support["RX"],
            support["RY"],
            support["RZ"],
        )
        return support_id

    def _add_supports_to_model(self):
        for support_id in self._data["supports"].keys():
            self._add_support_to_model(support_id)
        return "Appuis ajoutés"

    def get_all_supports(self) -> dict:
        return self._data["supports"]

    ################## Chargements ##################

    def _convert_pos(self, pos_index: int | str, member_id: str) -> int:
        """Converti la position en valeur recevable par la fonction create_load

        Args:
                pos_index (int | str): position de la charge sur la barre

        Returns:
                int: la position numérique sur la barre
        """
        match pos_index:
            case "start":
                pos_index = 0 * si.mm
            case "end":
                pos_index = self._data["members"][member_id]["Longueur"]
            case "middle":
                pos_index = round(self._data["members"][member_id]["Longueur"] / 2, 0)
            case str(x) if "%" in x:
                pos_index = pos_index.split("%")[0]
                pos_index.replace(" ", "")
                pos_index = round(
                    self._data["members"][member_id]["Longueur"] * int(pos_index) / 100,
                    0,
                )
            case _:
                pos_index = pos_index * si.mm
        return pos_index

    def create_dist_load(
        self,
        member_id: str,
        name: str,
        start_load: float,
        end_load: float,
        start_pos: str = None,
        end_pos: str = None,
        action: str = ACTION,
        direction: str = ("Fx", "Fy", "Fz", "FX", "FY", "FZ"),
        comment: str = None,
    ):
        """Ajoute une charge dans la liste de chargement de la classe chargement

        Args:
                member_id (str): id de la barre à charger.
                name (str): nom de la charge.
                start_load (int): effort de départ en kN/m.
                end_load (int): effort de fin en kN/m.
                start_pos (str, optional): position de début de la charge sur la barre en mm. En complément il est possible de mettre "start", "middle"
                                           ou un pourcentage pour définitr la position de la charge.
                end_pos (str, optional): position de début de la charge sur la barre en mm. En complément il est possible de mettre "end", "middle"
                                         ou un pourcentage pour définitr la position de la charge.
                action (str): type d'action de l'effort.
                direction (str): sens de l'effort sur la barre.
                comment (str, optional): commentaire sur la charge.
        """
        load_id = "L" + str(len(self._data["loads"]) + 1)
        type_load = "Distribuée"
        if not start_pos:
            start_pos = "start"
        if not end_pos:
            end_pos = "end"
        start_pos = self._convert_pos(start_pos, member_id)
        end_pos = self._convert_pos(end_pos, member_id)

        self._data["loads"][load_id] = {
            "N° barre": member_id,
            "Nom": name,
            "Action": action,
            "Type de charge": type_load,
            "Charge": {
                "start": start_load * si.kN / si.m,
                "end": end_load * si.kN / si.m,
            },
            "Position": {"start": start_pos, "end": end_pos},
            "Axe": direction,
            "Commentaire": comment,
        }
        return self._data["loads"][load_id]

    def create_point_load(
        self,
        member_id: str,
        name: str,
        load: int,
        pos: str = None,
        action: str = ACTION,
        direction: str = (
            "Fx",
            "Fy",
            "Fz",
            "Mx",
            "My",
            "Mz",
            "FX",
            "FY",
            "FZ",
            "MX",
            "MY",
            "MZ",
        ),
        comment: str = None,
    ):
        """Ajoute une charge dans la liste de chargement de la classe chargement

        Args:
                member_id (str): id de la barre à charger.
                name (str): nom de la charge.
                load (int): effort de départ en kN ou kN.m.
                pos (str, optional): position de la charge sur la barre en mm. En complément il est possible de mettre "start", "middle", "end"
                                     ou un pourcentage pour définitr la position de la charge.
                action (str): type d'action de l'effort.
                direction (str): sens de l'effort sur la barre.
                comment (str, optional): commentaire sur la charge.
        """
        load_id = "L" + str(len(self._data["loads"]) + 1)
        if "F" in direction:
            type_load = "Concentrée"
            load = load * si.kN
        else:
            type_load = "Moment"
            load = load * si.kN * si.m

        pos = self._convert_pos(pos, member_id)

        self._data["loads"][load_id] = {
            "N° barre": member_id,
            "Nom": name,
            "Action": action,
            "Type de charge": type_load,
            "Charge": load,
            "Position": pos,
            "Axe": direction,
            "Commentaire": comment,
        }
        return self._data["loads"][load_id]

    def create_load_by_list(
        self, list_loads: list, type_load: str = ("Distribuée", "Autre")
    ):
        """Ajoute les charges d'une liste pré-définit dans la liste de chargement

        Args:
                list_loads (list): liste de charge.
        """
        for load in list_loads:
            if type_load == "Distribuée":
                self.create_dist_load(*load)
            else:
                self.create_point_load(*load)

    def del_load(self, load_id: str):
        """Supprime une charge de l'attribut _data["loads"] par son index

        Args:
                index_load (int): index de la charge à supprimer.
        """
        self._data["loads"].pop(load_id)

    def _add_load_to_model(self, load_id: str):
        load = self._data["loads"][load_id]
        case = self.DICO_COMBI_ACTION[load["Action"]]
        match load["Type de charge"]:
            case "Distribuée":
                self._model.add_member_dist_load(
                    load["N° barre"],
                    load["Axe"],
                    load["Charge"]["start"].value * 10**-3,
                    load["Charge"]["end"].value * 10**-3,
                    load["Position"]["start"].value * 10**3,
                    load["Position"]["end"].value * 10**3,
                    case,
                )
            case _:
                member = self._data["members"][load["N° barre"]]
                long = member["Longueur"]
                local_axes = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")
                if load["Position"] == 0 and load["Axe"] not in local_axes:
                    node1 = member["Noeuds"][0]
                    self._model.add_node_load(
                        node1, load["Axe"], load["Charge"].value, case
                    )
                elif load["Position"] == long and load["Axe"] not in local_axes:
                    node2 = member["Noeuds"][1]
                    self._model.add_node_load(
                        node2, load["Axe"], load["Charge"].value, case
                    )
                else:
                    self._model.add_member_pt_load(
                        load["N° barre"],
                        load["Axe"],
                        load["Charge"].value,
                        load["Position"].value * 10**3,
                        case,
                    )
        return load_id

    def _add_loads_to_model(self):
        for load_id in self._data["loads"].keys():
            self._add_load_to_model(load_id)
        return self._model.load_cases

    def get_all_loads(self):
        """Retourne la liste des charges définits initialement."""
        return self._data["loads"]

    def get_member_loads(self, member_id: str) -> list:
        """Retourne la liste des charges définits initialement."""
        return [
            load
            for load in self._data["loads"].values()
            if load["N° barre"] == member_id
        ]

    def generate_model(self):
        self._model = FEModel3D()
        self._add_nodes_to_model()
        self._add_materials_to_model()
        self._add_sections_to_model()
        self._add_members_to_model()
        self._add_supports_to_model()
        self._add_loads_to_model()

    def _add_load_combos_to_model(self, combos: dict, tag: str):
        for combo, factor in combos.items():
            self._model.add_load_combo(combo, factor, tag)


class Model_result(Projet):
    ANALYZE_TYPE = ("Général", "Linéaire", "Second ordre")

    def __init__(
        self,
        model_generator: object,
        analyze_type: str = ANALYZE_TYPE,
        check_stability: bool = ("False", "True"),
        *args,
        **kwargs,
    ):
        """
        Cette classe permet de lancer l'analyse aux éléments finis et de récupérer et d'afficher les résultats du modèle.
        Cette classe est héritée de Projet dans le module A0_Projet.py.
        Elle nécessite toutefois la création d'un modèle de calcul avant d'être instanciée.
        Le modèle de calcul est créé par la classe Model_generator dans le module A0_Projet.py.

        Args:
                model_generator (Model_generator): le modèle de calcul à utiliser pour la génération des combinaisons d'action.
                analize_type (str): Définit le type d'analyse à réaliser
                check_stability (bool, optional): Définit si vous voulez vérifier la stabilité du modèle.
                    Ceci ralentit le calcul, à utiliser donc quand cas de débuguage. Defaults to False.
        """
        super().__init__(*args, **kwargs)
        self._model_generator = model_generator
        self.analyze_type = analyze_type
        self.check_stability = check_stability
        self._analyze()

    def _base_graph(
        self,
        title: str,
        combo_name: str,
        x_values,
        y_values,
        x_label: str,
        y_label: str,
        color: str,
        fill_between: bool = True,
        savefig: bool = False,
        filepath: str=None
    ):
        """Retourne un diagramme"""
        plt.clf()  # Effacer le graphique précédent
        plt.figure(self.name, figsize=(12, 4))
        plt.gcf().subplots_adjust(
            left=0.1, bottom=0.25, right=0.9, top=0.75, wspace=0, hspace=0.95
        )

        # manager = plt.get_current_fig_manager()
        # manager.resize(*manager.window.maxsize())

        plt.plot(x_values, y_values, color=color)
        plt.title(f"{title}\n{combo_name}", color=color)
        plt.ylabel(y_label)
        plt.xlabel(x_label)

        if fill_between:
            plt.fill_between(x_values, y_values, 0, color=color, alpha=0.2)
        plt.grid()
        if savefig:
            if not filepath:
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
            plt.savefig(filepath)
            return filepath
        else:
            plt.show()
        

    def _analyze(self):
        """Lance l'analyse du modèle aux éléments finis."""
        orphaned_nodes = self._model_generator._model.orphaned_nodes()
        if orphaned_nodes:
            return f"Les noeuds suivants ne sont pas connectés: {orphaned_nodes}"
        else:
            if self.analyze_type == self.ANALYZE_TYPE[0]:
                self._model_generator._model.analyze(
                    check_stability=self.check_stability
                )
            elif self.analyze_type == self.ANALYZE_TYPE[1]:
                self._model_generator._model.analyze_linear(
                    check_stability=self.check_stability
                )
            else:
                self._model_generator._model.analyze_PDelta(
                    check_stability=self.check_stability
                )

    def get_internal_force(
        self,
        member_id: str,
        combination: str,
        type: str = ("Nx", "Vy", "Vz", "Mx", "My", "Mz"),
        n_points: int = 20,
    ) -> np.array:
        """Retourne une table des efforts internes d'une membrure pour le type d'effort donné.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
            type (str): Le type d'effort interne à retourner. Defaults to ("Nx", "Vy", "Vz", "Mx", "My", "Mz").
            n_points (int, optional): le nombre de valeur à retrouner le long de la membrure. Defaults to 20.
        """
        match type:
            case "Nx":
                return self._model_generator._model.members[member_id].axial_array(
                    n_points=n_points, combo_name=combination
                )
            case "Vy":
                return self._model_generator._model.members[member_id].shear_array(
                    "Fy", n_points=n_points, combo_name=combination
                )
            case "Vz":
                return self._model_generator._model.members[member_id].shear_array(
                    "Fz", n_points=n_points, combo_name=combination
                )
            case "Mx":
                return self._model_generator._model.members[member_id].torque_array(
                    n_points=n_points, combo_name=combination
                )
            case "My":
                return self._model_generator._model.members[member_id].moment_array(
                    "My", n_points=n_points, combo_name=combination
                )
            case "Mz":
                return self._model_generator._model.members[member_id].moment_array(
                    "Mz", n_points=n_points, combo_name=combination
                )

    def get_min_max_internal_force(self, member_id: str, combination: str) -> np.array:
        """Retourne le maximum et minimum des efforts internes d'une membrure donnée.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
        """
        dict_internal_forces = {}
        for type in ("Nx", "Vy", "Vz", "Mx", "My", "Mz"):
            if type == "Nx":
                max = self._model_generator._model.members[member_id].max_axial(
                    combo_name=combination
                )
                min = self._model_generator._model.members[member_id].min_axial(
                    combo_name=combination
                )
            elif type == "Vy":
                max = self._model_generator._model.members[member_id].max_shear(
                    "Fy", combo_name=combination
                )
                min = self._model_generator._model.members[member_id].min_shear(
                    "Fy", combo_name=combination
                )
            elif type == "Vz":
                max = self._model_generator._model.members[member_id].max_shear(
                    "Fz", combo_name=combination
                )
                min = self._model_generator._model.members[member_id].min_shear(
                    "Fz", combo_name=combination
                )
            elif type == "Mx":
                max = self._model_generator._model.members[member_id].max_torque(
                    combo_name=combination
                )
                min = self._model_generator._model.members[member_id].min_torque(
                    combo_name=combination
                )
            elif type == "My":
                max = self._model_generator._model.members[member_id].max_moment(
                    "My", combo_name=combination
                )
                min = self._model_generator._model.members[member_id].min_moment(
                    "My", combo_name=combination
                )
            elif type == "Mz":
                max = self._model_generator._model.members[member_id].max_moment(
                    "Mz", combo_name=combination
                )
                min = self._model_generator._model.members[member_id].min_moment(
                    "Mz", combo_name=combination
                )
            if "M" in type:
                si_unit = si.N * si.mm
            else:
                si_unit = si.N
            dict_internal_forces[type] = {"Min": min * si_unit, "Max": max * si_unit}
        return dict_internal_forces
    
    def get_absolute_internal_force(self, member_id: str, combination: str, type: str = ("Nx", "Vy", "Vz", "Mx", "My", "Mz")):
        """Retourne la valeur de d'effort absolue pour le type d'effort donné.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
            type (str): Le type d'effort interne à retourner. Defaults to ("Nx", "Vy", "Vz", "Mx", "My", "Mz").
        """
        ei = self.get_min_max_internal_force(member_id, combination)
        return max(abs(ei[type]["Min"]), ei[type]["Max"])

    def show_internal_force_of_member(
        self,
        member_id: str,
        combination: str,
        type: str = ("Nx", "Vy", "Vz", "Mx", "My", "Mz"),
        n_points: int = 20,
        screenshot: bool = ("False", "True"),
        filepath: str=None,
    ):
        """Retourne un graphique des efforts internes d'une membrure suivant le type d'effort et la combinaison choisit.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
            type (str): Le type d'effort interne à retourner. Defaults to ("Nx", "Vy", "Vz", "Mx", "My", "Mz").
            n_points (int, optional): le nombre de valeur à retrouner le long de la membrure. Defaults to 20.
            screenshot (bool, optional): Définit si l'on souhaite enregistrer un screenshot du graph, si oui alors True. Defaults to False
            filepath (str, optional): Si le screenshot est souhaité, on peut rentrer un chemin de sauvegarde automatique, ce qui permettra je générer dynamiquement un rapport.
                                      Sinon le logiciel vous demande à chaque fois le chemin de sauvegarde.
        """
        x_label = "Longueur (mm)"
        x_value, y_value = self.get_internal_force(
            member_id, combination, type, n_points
        )
        if type.startswith("N"):
            title = f"Barre {member_id}: Effort normal {type}"
            color = "orange"
            y_label = "Effort (kN)"
            y_value = y_value * 10**-3
        elif type.startswith("V"):
            title = f"Barre {member_id}: Cisaillement {type}"
            color = "b"
            y_label = "Effort (kN)"
            y_value = y_value * 10**-3
        else:
            title = f"Barre {member_id}: Moments {type}"
            color = "r"
            y_label = "Effort (kN.m)"
            y_value = y_value * 10**-6
        return self._base_graph(
            title,
            combination,
            x_value,
            y_value,
            x_label,
            y_label,
            color=color,
            savefig=screenshot,
            filepath=filepath
        )

    def get_deflection(
        self,
        member_id: str,
        combination: str,
        direction: str = ("dx", "dy", "dz"),
        n_points: int = 20,
    ) -> np.array:
        """Retourne une table des déformation d'une membrure pour la direction locale donnée.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
            direction (str): La direction locale à retourner. Defaults to ("dx", "dy", "dz").
            n_points (int, optional): le nombre de valeur à retrouner le long de la membrure. Defaults to 20.
        """
        return self._model_generator._model.members[member_id].deflection_array(
            direction, n_points=n_points, combo_name=combination
        )

    def get_min_max_deflection(self, member_id: str, combination: str) -> np.array:
        """Retourne le maximum et minimum des efforts internes d'une membrure donnée.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
        """
        dict_deflection = {}
        for type in ("dx", "dy", "dz"):
            max = self._model_generator._model.members[member_id].max_deflection(
                type, combo_name=combination
            )
            min = self._model_generator._model.members[member_id].min_deflection(
                type, combo_name=combination
            )
            dict_deflection[type] = {"Min": min * si.mm, "Max": max * si.mm}
        return dict_deflection
    
    def get_absolute_max_deflection(self, member_id: str, combination: str, direction: str = ("dx", "dy", "dz")):
        """Retourne la valeur de déplacement absolue pour la direction de la flèche donnée.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
            direction (str): La direction locale à retourner. Defaults to ("dx", "dy", "dz").
        """
        deflection = self.get_min_max_deflection(member_id, combination)
        return max(abs(deflection[direction]["Min"]), deflection[direction]["Max"])
        

    def show_deflection_of_member(
        self,
        member_id: str,
        combination: str,
        direction: str = ("dx", "dy", "dz"),
        n_points: int = 20,
        screenshot: bool = ("False", "True"),
        filepath: str=None,
    ):
        """Retourne un graphique des efforts internes d'une membrure suivant le type d'effort et la combinaison choisit.

        Args:
            member_id (str): Le nom de la membrure à analyser
            combination (str): Le nom de la combinaison à récupérer
            direction (str): La direction locale à retourner. Defaults to ("dx", "dy", "dz").
            n_points (int, optional): le nombre de valeur à retrouner le long de la membrure. Defaults to 20.
            screenshot (bool, optional): Définit si l'on souhaite enregistrer un screenshot du graph, si oui alors True. Defaults to False
            filepath (str, optional): Si le screenshot est souhaité, on peut rentrer un chemin de sauvegarde automatique, ce qui permettra je générer dynamiquement un rapport.
                                      Sinon le logiciel vous demande à chaque fois le chemin de sauvegarde.
        """
        title = f'Barre {member_id}: Flèche {direction}'
        x_label = "Longueur (mm)"
        y_label = "Déplacement\n(mm)"
        color = "g"
        x_value, y_value = self.get_deflection(
            member_id, combination, direction, n_points
        )
        return self._base_graph(
            title,
            combination,
            x_value,
            y_value,
            x_label,
            y_label,
            color=color,
            savefig=screenshot,
            filepath=filepath
        )


    def show_model(
        self,
        combination: str,
        annotation_size: int = 70,
        deformed_shape: bool = True,
        deformed_scale: int = 20,
        screenshot: bool = ("False", "True"),
        filepath: str=None,
    ):
        """Retourne le model dans un graphique 3D.

        Args:
            combination (str): Nom de la combinaison à afficher.
            annotation_size (int, optional): Tailles des annotations. Defaults to 70.
            deformed_shape (bool, optional): Définit si l'on doit afficher les déformations globales. Defaults to True.
            deformed_scale (int, optional): Définit l'effet d'échelle des déformations. Defaults to 20.
            screenshot (bool, optional): Définit si l'on souhaite enregistrer un screenshot du graph, si oui alors True. Defaults to False
            filepath (str, optional): Si le screenshot est souhaité, on peut rentrer un chemin de sauvegarde automatique, ce qui permettra je générer dynamiquement un rapport.
                                      Sinon le logiciel vous demande à chaque fois le chemin de sauvegarde.
        """
        renderer = Renderer(self._model_generator._model)
        renderer.combo_name = combination
        renderer.annotation_size = annotation_size
        renderer.deformed_shape = deformed_shape
        renderer.deformed_scale = deformed_scale
        if screenshot:
            if not filepath:
                interaction = True
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
                QMessageBox.information(
                    None,
                    "Screenshot",
                    "Vous pouvez bouger le modèle pour prendre le screenshot.\nUne fois prêt, cliquer sur Q pour faire le screenshot.",
                )
            else:
                interaction = False
            renderer.screenshot(filepath, interact=interaction)
            return filepath
        else:
            renderer.render_model()

    def get_global_displacement_of_node(self, node_id: str, combination: str):
        """Retourne un dictionnaire des coordonnées de déplacment global d'un noeud dans le model pour la combinaison définit.

        Args:
            node_id (str): id du noeud à récupérer.
            combination (str): Nom de la combinaison à analyser.
        """
        DX = self._model_generator._model.nodes[node_id].DX[combination] * si.mm
        DY = self._model_generator._model.nodes[node_id].DY[combination] * si.mm
        DZ = self._model_generator._model.nodes[node_id].DZ[combination] * si.mm
        RX = self._model_generator._model.nodes[node_id].RX[combination]
        RY = self._model_generator._model.nodes[node_id].RY[combination]
        RZ = self._model_generator._model.nodes[node_id].RZ[combination]
        return {"DX": DX, "DY": DY, "DZ": DZ, "RX": RX, "RY": RY, "RZ": RZ}

    def _get_node_reaction(self, node_id: str, combination: str):
        """Retourne un dictionnaire des réaction d'un noeud.

        Args:
            node_id (str): id du noeud à récupérer.
            combination (str): Nom de la combinaison à analyser.
        """
        FX = self._model_generator._model.nodes[node_id].RxnFX[combination] * si.N
        FY = self._model_generator._model.nodes[node_id].RxnFY[combination] * si.N
        FZ = self._model_generator._model.nodes[node_id].RxnFZ[combination] * si.N
        MX = (
            self._model_generator._model.nodes[node_id].RxnMX[combination]
            * si.N
            * si.mm
        )
        MY = (
            self._model_generator._model.nodes[node_id].RxnMY[combination]
            * si.N
            * si.mm
        )
        MZ = (
            self._model_generator._model.nodes[node_id].RxnMZ[combination]
            * si.N
            * si.mm
        )
        return {"FX": FX, "FY": FY, "FZ": FZ, "MX": MX, "MY": MY, "MZ": MZ}

    def get_supports_reactions(self, combination: str):
        """Retourne un dictionnaire des réaction aux appuis.

        Args:
            combination (str): Nom de la combinaison à analyser.
        """
        reaction = {}
        supports = self._model_generator.get_all_supports()
        for support_id, support in supports.items():
            node_id = support["Noeud"]
            reaction[support_id] = self._get_node_reaction(node_id, combination)
        return reaction


if __name__ == "__main__":
    # projet = Projet(num_project="6006.0", commmentaire="c'est mon premier projet")
    # building = Batiment(h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=15)
    # print(building.h_bat)
    beam_gen = Model_generator()
    node1 = beam_gen.add_node(0, 0, 0)
    node2 = beam_gen.add_node(5000, 0, 0)
    inertie = beam_gen.inertie(100, 200, "Rectangulaire")
    aire = beam_gen.aire(100, 200, "Rectangulaire")
    memb1 = beam_gen.add_member(node1, node2, aire, inertie["Iy"], inertie["Iz"])
    beam_gen.add_material_by_class(memb1, "C24")

    listdeplacement = [
        ["N1", "Rotule YZ", 0],
        ["N2", "Rotule YZ", 0],
    ]
    beam_gen.create_supports_by_list(listdeplacement)
    beam_gen.export_data()
    print(beam_gen._data)
