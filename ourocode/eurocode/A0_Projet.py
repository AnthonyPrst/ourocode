#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os, sys
import math as mt
import numpy as np
from matplotlib import pyplot as plt
import forallpeople as si
from handcalcs.decorator import handcalc

sys.path.append(os.path.join(os.getcwd(), "ourocode"))
from eurocode.objet import Objet

# from ourocode.eurocode.objet import Objet


class Projet(Objet):
    JUPYTER_DISPLAY = False

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
        self.b_bat = b_bat  # coté perpendiculaire au vent longpant
        self.alpha_toit = alpha_toit
        self.alpha_toit2 = alpha_toit2

    def height_ref(self):
        pass  # 7.2.1 bâtiment de grande hauteur


class Model_generator(Projet):
    LIST_SECTION = ["Rectangulaire", "Circulaire"]
    CLASSE_WOOD = tuple(
        Projet._data_from_csv(Projet, "caracteristique_meca_bois.csv").index
    )[2:]

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
        self._data = {"nodes": {}, "members": {}, "supports": {}}

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

    def aire(self, b: int, h: int, section: str = LIST_SECTION):
        if section == self.LIST_SECTION[0]:
            return b * h
        else:
            return mt.pi * (b / 2) ** 2

    def inertie(self, b: int, h: int, section: str = LIST_SECTION):
        """Retourne le moment quadratique d'une section rectangulaire en mm4 avec pour argument :
        b ou d : Largeur ou diamètre de la poutre en mm
        h : Hauteur de la poutre en mm"""
        if section == "Rectangulaire":
            I_y = (b * h**3) / 12
            I_z = (h * b**3) / 12
        else:
            I_y = (mt.pi * b**4) / 64
            I_z = I_y
        return {"Iy": I_y, "Iz": I_z}

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
            angle_local = np.arctan2(
                np.linalg.norm(np.cross(ref_vec, T)), np.dot(ref_vec, T)
            )

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
            "local": _get_local_angle_of_bar(v1),
        }

        return angles
        # ang1 = np.arctan2(*v1[::-1])
        # return np.rad2deg(ang1 % (2 * np.pi))

    def add_node(self, X: int, Y: int, Z: int, comment: str = None) -> str:
        """Ajoute un noeud au model MEF

        Args:
            X (int): position en X global en mm
            Y (int): position en Y global en mm
            Z (int): position en Z global en mm
        """
        node_id = "N" + str(len(self._data["nodes"]) + 1)
        self._data["nodes"][node_id] = {"X": X, "Y": Y, "Z": Z, "Commentaire": comment}
        return node_id

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

    def get_all_members(self) -> dict:
        return self._data["members"]

    def get_all_supports(self) -> dict:
        return self._data["supports"]

    def add_member(
        self,
        node1: str,
        node2: str,
        aire: float,
        Iy: float,
        Iz: float,
        comment: str = None,
    ):
        """Ajoute une poutre au model MEF

        Args:
            node1 (str): id du noeud 1
            node2 (str): id du noeud 2
            aire (float): airede la section en mm²

            ATTENTION: Iy est la grande inertie et Iz est la petite inertie !
            Iy (float): Inertie quadratique autour de y en mm4
            Iz (float): Inertie quadratique autour de z en mm4
        """
        member_id = "M" + str(len(self._data["members"]) + 1)
        node_coor_1 = self.get_node(node1)
        node_coor_2 = self.get_node(node2)
        x1, y1, z1 = node_coor_1["X"], node_coor_1["Y"], node_coor_1["Z"]
        x2, y2, z2 = node_coor_2["X"], node_coor_2["Y"], node_coor_2["Z"]
        vector = (x2 - x1, y2 - y1, z2 - z1)
        length = mt.sqrt(
            abs(vector[0]) ** 2 + abs(vector[1]) ** 2 + abs(vector[2]) ** 2
        )
        angle = self._get_angle_of_bar(vector)

        self._data["members"][member_id] = {
            "Noeuds": [node1, node2],
            "Longueur": length,
            "Section": aire,
            "Iy": Iz,
            "Iz": Iy,
            "Angle": angle,
            "Relaxation": {"start": None, "end": None},
            "Commentaire": comment,
        }
        return member_id

    def get_member(self, member_id: str) -> dict:
        """Retourne une membrure par son id

        Args:
            node_id (str): id de la membrure à récupérer

        Returns:
            dict: dictionnaire de la membrure
        """
        return self._data["members"][member_id]

    def add_material_by_class(self, member_id: int, classe: str = CLASSE_WOOD) -> str:
        data_csv_meca = self._data_from_csv("caracteristique_meca_bois.csv")
        material_properties = data_csv_meca.loc[classe]
        E = int(material_properties.loc["E0mean"])
        G = int(material_properties.loc["Gmoy"])
        J = 0
        nu = (E / (2 * G)) - 1  # Coefficient de Poisson
        self._data["members"][member_id]["Matériaux"] = {
            "classe": classe,
            "E": E,
            "G": G,
            "J": J,
            "nu": nu,
        }
        return self._data["members"][member_id]["Matériaux"]

    def add_material_by_mechanical_properties(
        self, member_id: int, E: int, G: float, J: float, nu: float
    ):
        """Ajoute un matériau à une barre par ces caractéristiques mécaniques.

        Args:
            member_id (int): Le numéro de la barre à laquelle ajouter le matériaux
            E (int): Module de young en MPa, ce E est le E,mean. Il ne faut absolument pas donner le E,mean,fin sous peine de réaliser le calcul EC5 §2.3.2.2 equ2.7 deux fois !
            A (float): Section en mm²
            G (float): Module de cisaillement en MPa
            J (float): Module de torsion en mm4
            nu (float): Coefficient de Poisson
        """
        # ATTENTION: Iy / Iz est différent des Eurocodes ! y est verticale et est donc dans le sens théorique eurocode de Z ! z est donc dans le sens théorique de Y !
        # Il faut donc faire attention au dimenssion de l'élément ! Pour une poutre de section b: 100mm h: 200mm -> Iy = 200 * 100^3 / 12 et Iz = 100 * 200^3 / 12.
        self._data["members"][member_id]["Matériaux"] = {
            "classe": "Manuel",
            "E": E,
            "G": G,
            "J": J,
            "nu": nu,
        }
        return self._data["members"][member_id]["Matériaux"]

    def add_release(
        self,
        member_id: int,
        position: str = ("start", "end"),
        u: bool = False,
        v: bool = False,
        w: bool = False,
        teta_x: bool = False,
        teta_y: bool = True,
        teta_z: bool = True,
    ):
        """Ajoute une relaxation sur une membrure soit au début, soit à la fin.
        Ceci est considéré dans les MEF par une matrice de rigidité spécifique avec les éléments relaché égale à 0.

        Args:
            member_id (int): numéro de la membrure à relacher
            position (str, optional): position du relachement sur la barre soit au début soit à la fin. Defaults to ("start", "end").
            u (bool, optional): relachement de l'axe x local, si oui alors True. Defaults to False.
            v (bool, optional): relachement de l'axe y local, si oui alors True. Defaults to False.
            w (bool, optional): relachement de l'axe z local, si oui alors True. Defaults to False.
            teta_x (bool, optional): relachement de l'axe de rotation x local, si oui alors True. Attention de base cette rotation doit toujours être fixé. Defaults to False.
            teta_y (bool, optional): relachement de l'axe de rotation y local, si oui alors True. Defaults to True.
            teta_z (bool, optional): relachement de l'axe de rotation z local, si oui alors True. Defaults to True.
        """
        self._data[member_id]["relaxation"][position](
            {
                "position": position,
                "u": u,
                "v": v,
                "w": w,
                "teta_x": teta_x,
                "teta_y": teta_y,
                "teta_z": teta_z,
            }
        )
        return self._data[member_id]["relaxation"][-1]

    def create_support(
        self,
        node_id: int,
        type_appuis: str = (
            "Simple X",
            "Simple Y",
            "Simple Z",
            "Rotule YZ",
            "Rotule XY",
            "Rotule XZ",
            "Encastrement",
        ),
        l_appuis: int = 0,
    ):
        """Ajoute un appuis dans la liste d'appuis de la classe MEF

        Args:
            node_id (int): Numéro du noeud sur lequel positionner l'appuis.
            type_appuis (str, optional): type d'appuis à créer. Defaults to ("Simple", 'Rotule', 'Encastrement').
            l_appuis (int, optional): longueur d'appuis sur la poutre en mm. Defaults to 0.
        """
        support_id = "S" + str(len(self._data["supports"]) + 1)
        if type_appuis not in (
            "Simple X",
            "Simple Y",
            "Simple Z",
            "Rotule YZ",
            "Rotule XY",
            "Rotule XZ",
            "Encastrement",
        ):
            raise ValueError(
                f"Le type d'appuis {type_appuis} n'est pas reconnu. Les types d'appuis disponibles sont: {self._data['supports']}"
            )

        self._data["supports"][support_id] = {
            "Noeud": node_id,
            "Type d'appui": type_appuis,
            "Longueur d'appui": l_appuis,
        }
        return self._data["supports"][support_id]

    def create_supports_by_list(self, list_supports: list):
        """Ajoute les charges d'une liste pré-défini dans la liste de chargement

        Args:
            list_supports (list): liste de charge.
        """
        for support in list_supports:
            self.create_support(*support)
        return self._data["supports"]

    def del_support(self, support_id: int):
        """Supprime un appui de l'attribut list_supports par son index

        Args:
            support_id (int): id de l'appuis à supprimer.
        """
        return f"L'appui à été supprimé: {self._data["supports"].pop(support_id)}"

    def get_supports(self):
        """Retourne la liste des appuis définis."""
        return self._data["supports"]


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
