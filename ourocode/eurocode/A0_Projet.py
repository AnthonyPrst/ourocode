#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import os, sys
import math as mt
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image

import forallpeople as si
from handcalcs.decorator import handcalc

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# from eurocode.objet import Objet

from ourocode.eurocode.objet import Objet


class Projet(Objet):
    """Classe définissant les informations générales d'un projet de structure.

    Cette classe est la racine de la hiérarchie des objets ourocode.
    Elle contient les informations administratives et géographiques du projet,
    nécessaires pour l'application des annexes nationales françaises des Eurocodes.
    """

    DICO_COMBI_ACTION = {
        "Permanente G": "G",
        "Exploitation Q": "Q",
        "Neige normale Sn": "Sn",
        "Vent pression W+": "W+",
        "Vent dépression W-": "W-",
        "Neige accidentelle Sx": "Sx",
        "Sismique Ae": "Ae",
    }
    CAT_TYPE = (
        "Aucune",
        "Cat A : habitation",
        "Cat B : bureaux",
        "Cat C : lieu de réunion",
        "Cat D : zones commerciales",
        "Cat E : stockage",
        "Cat F : véhicule <= 30kN",
        "Cat G : véhicule <= 160kN",
        "Cat H : toits",
    )
    PAYS = ("France")

    def __init__(
        self,
        ingenieur: str = None,
        num_project: str = None,
        name: str = None,
        adresse: str = None,
        code_INSEE: int = None,
        pays: str = PAYS,
        alt: si.m = 0,
        **kwargs,
    ):
        """Initialise un projet avec ses informations générales.

        Cette classe définit le contexte du projet, indispensable pour
        l'application correcte des annexes nationales françaises des Eurocodes.
        Tous les objets de calcul (éléments de structure, charges, etc.)
        héritent indirectement de cette classe.

        Args:
            ingenieur (str, optional): Nom de l'ingénieur responsable.
                Defaults to None.
            num_project (str, optional): Numéro de référence du projet.
                Defaults to None.
            name (str, optional): Nom ou désignation du projet.
                Defaults to None.
            adresse (str, optional): Adresse géographique du chantier.
                Defaults to None.
            code_INSEE (int, optional): Code INSEE à 5 chiffres du département
                ou de la commune pour les données climatiques. Defaults to None.
            pays (str, optional): Pays où se situe le projet. Defaults to "France".
                Attention : ce package intègre uniquement les annexes nationales
                françaises des Eurocodes.
            alt (si.m, optional): Altitude du projet en mètres, utilisée pour
                le calcul de la neige et du vent. Defaults to 0.
            **kwargs: Arguments supplémentaires ajoutés dynamiquement.
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

class Batiment(Projet):
    """Classe définissant la géométrie d'un bâtiment pour les calculs de structure.

    Cette classe décrit les dimensions principales du bâtiment et les
    caractéristiques de sa toiture, nécessaires pour le calcul des
    charges climatiques (neige, vent) et sismiques.

    Attributs de classe:
        ETAGE (tuple): Liste des niveaux courants (RDC à Toiture).
    """

    ETAGE = ("RDC", "R+1", "R+2", "R+3", "R+4", "Toiture")

    def __init__(
        self,
        h_bat: si.m,
        d_bat: si.m,
        b_bat: si.m,
        alpha_toit: float,
        alpha_toit2: float = 0,
        *args,
        **kwargs,
    ):
        """Initialise les dimensions du bâtiment.

        Args:
            h_bat (si.m): Hauteur totale du bâtiment en mètres, mesurée depuis
                le soubassement rigide ou les fondations (référence sismique).
            d_bat (si.m): Largeur du bâtiment en mètres (dimension perpendiculaire
                au vent dominant pour le calcul du vent).
            b_bat (si.m): Longueur du bâtiment en mètres.
            alpha_toit (float): Pente du premier versant de toiture en degrés.
                0° pour un toit plat, valeur positive pour un versant.
            alpha_toit2 (float, optional): Pente du second versant pour les
                toits à deux pans (0° si toit à un seul versant ou plat).
                Defaults to 0.
            *args: Arguments positionnels transmis à la classe parent Projet.
            **kwargs: Arguments nommés transmis à la classe parent Projet
                (ingenieur, code_INSEE, alt, etc.).
        """
        super().__init__(*args, **kwargs)
        self.h_bat = h_bat * si.m
        self.d_bat = d_bat * si.m
        self.b_bat = b_bat * si.m  # coté perpendiculaire au vent longpant
        self.alpha_toit = alpha_toit
        self.alpha_toit2 = alpha_toit2


class Model_generator(Projet):
    """Générateur de modèles de calcul par éléments finis (MEF).

    Cette classe permet de construire un modèle de structure complet pour
    le calcul par éléments finis avec Pynite. Elle gère la définition des
    nœuds, barres, sections, matériaux, appuis et chargements.

    Le flux de travail typique est :
        1. Définir les sections avec `add_section`
        2. Définir les matériaux avec `add_material_by_class` ou `add_material_by_mechanical_properties`
        3. Créer les nœuds avec `add_node`
        4. Créer les barres avec `add_member`
        5. Définir les appuis avec `add_support` ou `add_support_spring`
        6. Appliquer les charges avec `create_dist_load` et `create_point_load`
        7. Lancer les combinaisons via la classe `Combinaison`
        8. Analyser les résultats avec `Model_result`
    """

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
        """Initialise un générateur de modèle MEF vide.

        Args:
            *args: Arguments positionnels transmis à Projet.
            **kwargs: Arguments nommés transmis à Projet.
        """
        super().__init__(*args, **kwargs)
        self._data = {
            "nodes": {},
            "sections": {},
            "materials": {},
            "members": {},
            "supports": {"classic": {}, "spring": {}},
            "loads": {},
        }
        self._model = None

    def get_all_data(self) -> dict:
        """Retourne l'ensemble des données du modèle MEF.

        Returns:
            dict: Dictionnaire contenant toutes les données structurées :
                - "nodes": dictionnaire des nœuds
                - "sections": dictionnaire des sections
                - "materials": dictionnaire des matériaux
                - "members": dictionnaire des barres
                - "supports": dictionnaire des appuis (classiques et ressorts)
                - "loads": dictionnaire des chargements
        """
        return self._data

    def export_data(self):
        """Exporte les données du modèle au format JSON via boîte de dialogue.

        Ouvre une boîte de dialogue Qt pour choisir l'emplacement du fichier JSON.
        Le fichier contient l'intégralité des données du modèle (nœuds, barres,
        matériaux, sections, appuis et charges).

        Note:
            L'export est utile pour sauvegarder un modèle ou le transférer
            vers une autre application.
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
        """Affiche l'image de la convention de signe pour les efforts sur barre.

        Ouvre une fenêtre affichant le schéma de la convention de signes
        utilisée par Pynite pour les efforts internes (efforts normaux,
        tranchants, moments fléchissants) et les déplacements.

        Note:
            La convention de signe suit le repère local de chaque barre.
        """
        file = os.path.join(
            self.PATH_CATALOG, "data", "screenshot", "sign_convention.png"
        )
        image = Image.open(file)
        image.show()

    ################## Noeud ##################

    def add_node(self, X: si.mm, Y: si.mm, Z: si.mm, comment: str = None) -> str:
        """Ajoute un nœud au modèle MEF.

        Crée un nœud dans le système de coordonnées globales (X, Y, Z).
        L'identifiant est généré automatiquement sous la forme "N1", "N2", etc.

        Args:
            X (si.mm): Coordonnée X dans le repère global, en millimètres.
            Y (si.mm): Coordonnée Y dans le repère global, en millimètres.
            Z (si.mm): Coordonnée Z dans le repère global, en millimètres.
            comment (str, optional): Commentaire descriptif pour le nœud.

        Returns:
            str: Identifiant unique du nœud créé (ex: "N1").

        Exemple:
            >>> model = Model_generator()
            >>> n1 = model.add_node(0, 0, 0, comment="Appui gauche")
            >>> n2 = model.add_node(5000, 0, 0, comment="Appui droit")
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
        """Retourne les coordonnées et informations d'un nœud.

        Args:
            node_id (str): Identifiant du nœud à récupérer (ex: "N1").

        Returns:
            dict: Dictionnaire contenant :
                - "X", "Y", "Z": coordonnées avec unités (si.mm)
                - "Commentaire": texte descriptif (ou None)

        Raises:
            KeyError: Si le nœud n'existe pas.
        """
        return self._data["nodes"][node_id]

    def get_all_nodes(self) -> dict:
        """Retourne l'ensemble des nœuds du modèle.

        Returns:
            dict: Dictionnaire de tous les nœuds avec leurs coordonnées,
                indexé par les identifiants de nœuds.
        """
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
        poids_propre: bool = ("True", "False"),
        rotation: float = 0,
        tension_only: bool = ("False", "True"),
        compression_only: bool = ("False", "True"),
        name: str = None,
        comment: str = None,
    ):
        """Ajoute une barre (poutre ou colonne) au modèle MEF.

        Crée une barre entre deux nœuds existants, avec un matériau et une section
        prédéfinis. La longueur est calculée automatiquement d'après les coordonnées.

        Args:
            node1 (str): Identifiant du nœud de départ (ex: "N1").
            node2 (str): Identifiant du nœud d'arrivée (ex: "N2").
            material (str): Identifiant du matériau (créé via add_material_by_class).
            section (str): Identifiant de la section (créé via add_section).
            poids_propre (bool, optional): Si True, génère automatiquement une
                charge répartie correspondant au poids propre de la barre.
                Defaults to True.
            rotation (float, optional): Angle de rotation de la section en degrés
                autour de l'axe longitudinal. Defaults to 0.
            tension_only (bool, optional): Si True, la barre travaille uniquement
                en traction (ex: tirant). Defaults to False.
            compression_only (bool, optional): Si True, la barre travaille uniquement
                en compression (ex: étai). Defaults to False.
            name (str, optional): Nom personnalisé de la barre (doit être unique).
                Si None, un identifiant automatique "M1", "M2"... est généré.
            comment (str, optional): Commentaire descriptif pour la barre.

        Returns:
            str: Identifiant unique de la barre créée.

        Raises:
            KeyError: Si les nœuds, le matériau ou la section n'existent pas.

        Note:
            Les arguments tension_only et compression_only sont mutuellement exclusifs.
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
            "Poids propre": poids_propre,
            "Rotation": rotation,
            "Relaxation": {"start": None, "end": None},
            "Commentaire": comment,
        }
        if poids_propre:
            rho = self._data["materials"][material]["rho"].value # en kg/m^3
            aire = self._data["sections"][section]["Aire"].value
            pp = rho * aire * 10**-2
            # Important: appeler explicitement la méthode de la classe de base pour éviter
            # le dispatch vers Wood_beam_model.create_dist_load (signature différente)
            Model_generator.create_dist_load(
                self,
                member_id,
                f"PP_{member_id}",
                -pp,
                -pp,
                "start",
                "end",
                action="Permanente G",
                direction="FY",
                comment=f"Poids propre {member_id}")
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
        """Retourne les informations d'une barre par son identifiant.

        Args:
            member_id (str): Identifiant de la barre à récupérer (ex: "M1").

        Returns:
            dict: Dictionnaire contenant les propriétés de la barre :
                - "Noeuds": liste [node1_id, node2_id]
                - "Longueur": longueur calculée en mm
                - "Section": identifiant de la section
                - "Matériaux": identifiant du matériau
                - "Rotation": angle de rotation en degrés
                - "Relaxation": dictionnaire des relâchements

        Raises:
            KeyError: Si la barre n'existe pas.
        """
        return self._data["members"][member_id]
    
    def get_member_length(self, member_id: str) -> si.mm:
        """Retourne la longueur d'une barre par son identifiant.

        Args:
            member_id (str): Identifiant de la barre (ex: "M1").

        Returns:
            si.mm: Longueur de la barre avec unité.
        """
        return self._data["members"][member_id]["Longueur"]

    def get_all_members(self) -> dict:
        """Retourne l'ensemble des barres du modèle.

        Returns:
            dict: Dictionnaire de toutes les barres, indexé par leurs identifiants.
        """
        return self._data["members"]

    ################## Matériaux ##################

    def add_material_by_class(self, classe: str = CLASSE_WOOD) -> str:
        """Ajoute un matériau bois au modèle par sa classe de résistance.

        Charge les caractéristiques mécaniques depuis les données normatives
        (caracteristique_meca_bois.csv) et crée le matériau avec :
        - Module de Young E (E0mean)
        - Module de cisaillement G (Gmoy)
        - Coefficient de Poisson nu (calculé)
        - Masse volumique rho (rhomean)

        Args:
            classe (str): Classe de résistance du bois selon l'EC5.
                Valeurs courantes : "C14", "C16", "C18", "C20", "C24", "C27", "C30",
                "D30", "D35", "D40", "D50", "D60", "D70",
                "GL20h", "GL24h", "GL28h", "GL32h", etc.

        Returns:
            str: Identifiant du matériau créé (égal à la classe fournie).

        Raises:
            KeyError: Si la classe n'existe pas dans la base de données.

        Note:
            Si le matériau existe déjà, la méthode retourne simplement son identifiant
            sans recréer les propriétés.
        """
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
        """Ajoute un matériau défini manuellement par ses propriétés mécaniques.

        Permet de créer un matériau personnalisé en spécifiant directement
        les caractéristiques mécaniques, utile pour les matériaux non standard
        ou les matériaux autres que le bois.

        Args:
            name (str): Nom unique du matériau.
            E (si.MPa): Module de Young (module d'élasticité longitudinal).
                Important : fournir E0mean, pas E0mean,fin (le facteur de fluage
                sera appliqué automatiquement si nécessaire par ailleurs).
            G (si.MPa): Module de cisaillement (module d'élasticité transversal).
            nu (float): Coefficient de Poisson (rapport des déformations transversale
                et longitudinale).
            rho (float): Masse volumique en kg/m³.

        Returns:
            str: Identifiant du matériau créé (égal au nom fourni).

        Note:
            Pour le bois, il est recommandé d'utiliser `add_material_by_class`
            qui garantit la cohérence avec les valeurs normatives.
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
        """Retourne les propriétés d'un matériau par son identifiant.

        Args:
            material_id (str): Identifiant du matériau (ex: "C24").

        Returns:
            dict: Dictionnaire contenant :
                - "classe": type de matériau (nom de classe ou "Manuel")
                - "E": Module de Young avec unité (si.MPa)
                - "G": Module de cisaillement avec unité (si.MPa)
                - "nu": Coefficient de Poisson
                - "rho": Masse volumique avec unité (si.kg/si.m**3)
        """
        return self._data["materials"][material_id]

    def get_all_materials(self) -> dict:
        """Retourne l'ensemble des matériaux du modèle.

        Returns:
            dict: Dictionnaire de tous les matériaux, indexé par leurs identifiants.
        """
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
        """Ajoute une section transversale au modèle.

        Crée une section rectangulaire ou circulaire avec calcul automatique
        de l'aire et des moments quadratiques d'inertie.

        Args:
            b (si.mm): Largeur de la section (ou diamètre pour section circulaire)
                en millimètres.
            h (si.mm): Hauteur de la section en millimètres.
                Ignoré pour les sections circulaires.
            J (si.mm**4): Module de torsion (constante de torsion) en mm⁴.
                Pour une section rectangulaire pleine : J ≈ k * b * h³ où k ≈ 0.33
                Pour une section circulaire pleine : J = π * d⁴ / 32
            section (str): Type de section. Valeurs acceptées :
                - "Rectangulaire" : section rectangulaire pleine
                - "Circulaire" : section circulaire pleine (b = diamètre)

        Returns:
            str: Identifiant unique de la section créé :
                - "R{b}X{h}" pour une section rectangulaire (ex: "R100X200")
                - "C{b}" pour une section circulaire (ex: "C300")

        Raises:
            ValueError: Si le type de section n'est pas reconnu.

        Note:
            L'identifiant est généré automatiquement à partir des dimensions.
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
        """Ajoute une section personnalisée définie par ses propriétés mécaniques.

        Permet de créer une section sans géométrie prédéfinie (IPE, HEA, etc.)
        en spécifiant directement l'aire et les moments d'inertie.

        Args:
            name (str): Nom unique de la section (ex: "IPE200", "HEA140").
            aire (si.mm**2): Aire de la section en mm².
            Iy (si.mm**4): Moment quadratique d'inertie autour de l'axe Y (faible inertie pour une section rectangulaire dans le logiciel).
                Note : Pour les sections rectangulaires, Iy = h × b³ / 12.
            Iz (si.mm**4): Moment quadratique d'inertie autour de l'axe Z (forte inertie pour une section rectangulaire dans le logiciel).
                Note : Pour les sections rectangulaires, Iz = b × h³ / 12.
            J (si.mm**4): Module de torsion (constante de torsion) en mm⁴.

        Returns:
            dict: Dictionnaire des propriétés de la section créée.
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
        """Retourne les propriétés d'une section par son identifiant.

        Args:
            section_id (str): Identifiant de la section (ex: "R100X200", "C300", "IPE200").

        Returns:
            dict: Dictionnaire contenant les propriétés de la section :
                - "Section": type de section ("Rectangulaire", "Circulaire", "Manuel")
                - "b", "h": dimensions pour les sections prédéfinies
                - "Aire": aire avec unité (si.mm**2)
                - "Iy", "Iz": inerties avec unités (si.mm**4)
                - "J": module de torsion avec unité (si.mm**4)
        """
        return self._data["sections"][section_id]

    def get_all_sections(self) -> dict:
        """Retourne l'ensemble des sections du modèle.

        Returns:
            dict: Dictionnaire de toutes les sections, indexé par leurs identifiants.
        """
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
        """Ajoute un relâchement (libération de degrés de liberté) sur une barre.

        Définit des conditions de dégagement aux extrémités d'une barre pour
        créer des rotules, des articulations ou des glissements partiels.
        La matrice de rigidité locale est modifiée en conséquence dans Pynite.

        Args:
            member_id (str): Identifiant de la barre à relâcher (ex: "M1").
            position (str, optional): Position du relâchement sur la barre.
                "start" = nœud de départ (i), "end" = nœud d'arrivée (j).
                Defaults to "start".
            u (bool, optional): Relâchement de la translation selon l'axe x local
                (longitudinal). Defaults to False.
            v (bool, optional): Relâchement de la translation selon l'axe y local.
                Defaults to False.
            w (bool, optional): Relâchement de la translation selon l'axe z local.
                Defaults to False.
            teta_x (bool, optional): Relâchement de la rotation selon l'axe x local
                (torsion). Par défaut toujours bloquée pour éviter les modes rigides.
                Defaults to False.
            teta_y (bool, optional): Relâchement de la rotation selon l'axe y local
                (moment fléchissant My). Defaults to False.
            teta_z (bool, optional): Relâchement de la rotation selon l'axe z local
                (moment fléchissant Mz). Defaults to False.

        Returns:
            dict: Configuration du relâchement créé pour la position spécifiée.

        Note:
            Un relâchement complet de la rotation (teta_y=True, teta_z=True)
            crée une articulation parfaite. Toutes les translations relâchées
            créent un glisseur.
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
        """Ajoute un appui classique (encastrement, rotule, glisseur) sur un nœud.

        Définit les conditions de déplacement (translation et rotation) bloquées
        ou libres en chaque nœud selon le repère global (X, Y, Z).

        Args:
            node_id (str): Identifiant du nœud d'appui (ex: "N1").
            DX (bool, optional): Bloque la translation selon l'axe X global.
                True = bloqué, False = libre. Defaults to True.
            DY (bool, optional): Bloque la translation selon l'axe Y global.
                Defaults to True.
            DZ (bool, optional): Bloque la translation selon l'axe Z global.
                Defaults to True.
            RX (bool, optional): Bloque la rotation selon l'axe X global (torsion).
                Defaults to True.
            RY (bool, optional): Bloque la rotation selon l'axe Y global.
                Defaults to True.
            RZ (bool, optional): Bloque la rotation selon l'axe Z global.
                Defaults to True.
            l_appuis (int, optional): Longueur d'appui sur la poutre en mm,
                utilisée pour la vérification à la compression perpendiculaire.
                Defaults to 0.

        Returns:
            dict: Configuration de l'appui créé avec son identifiant généré.
        """
        support_id = "S" + str(len(self._data["supports"]["classic"]) + 1)
        self._data["supports"]["classic"][support_id] = {
            "Noeud": node_id,
            "DX": DX,
            "DY": DY,
            "DZ": DZ,
            "RX": RX,
            "RY": RY,
            "RZ": RZ,
            "Longueur d'appui": l_appuis,
        }
        return self._data["supports"]["classic"][support_id]

    def add_support_spring(self,
        node_id: str,
        dof: str = ("DX", "DY", "DZ", "RX", "RY", "RZ"),
        stiffness: si.kN/si.m = 0,
        limit_direction: str = ("Aucune limitation", "Tension uniquement", "Compression uniquement")
    ):
        """Ajoute un appui élastique (ressort) sur un nœud dans une direction donnée.

        Modélise un appui avec raideur finie (sol élastique, appui flexible,
        suspension) ou un tirant/compression unidirectionnel.

        Args:
            node_id (str): Identifiant du nœud d'appui (ex: "N1").
            dof (str, optional): Degré de liberté sur lequel appliquer la raideur.
                "DX", "DY", "DZ" pour les translations, "RX", "RY", "RZ" pour les rotations.
                Defaults to "DX".
            stiffness (si.kN/si.m): Raideur de l'appui.
                - Pour DX/DY/DZ : kN/m (translation)
                - Pour RX/RY/RZ : kN·m/rad (rotation)
                Defaults to 0.
            limit_direction (str, optional): Limite le comportement du ressort
                à un seul sens de sollicitation.
                - "Aucune limitation" : ressort bidirectionnel
                - "Tension uniquement" : ne travaille qu'en traction (ex: tirant)
                - "Compression uniquement" : ne travaille qu'en compression (ex: étai)
                Defaults to "Aucune limitation".

        Returns:
            dict: Configuration de l'appui ressort créé.
        """
        support_id = "SSpring" + str(len(self._data["supports"]["spring"]) + 1)
        self._data["supports"]["spring"][support_id] = {
            "Noeud": node_id,
            "Dof": dof,
            "Raideur": stiffness * si.kN/si.m,
            "Limite de direction": limit_direction,
        }
        return self._data["supports"]["spring"][support_id]

    def create_supports_by_list(self, list_supports: list):
        """Ajoute les support d'une liste pré-définit.

        Args:
            list_supports (list): liste de support.
        """
        for support in list_supports:
            self.add_support(*support)

    def del_support(self, support_id: str):
        """Supprime un appui classique par son identifiant.

        Args:
            support_id (str): Identifiant de l'appui à supprimer (ex: "S1").

        Returns:
            str: Message de confirmation avec les détails de l'appui supprimé.
        """
        return f"L'appui à été supprimé: {self._data["supports"]["classic"].pop(support_id)}"

    def _add_support_to_model(self, support_id: str, support_type: str):
        """Ajoute un appui au model MEF

        Args:
            support_id (str): id de l'appui à ajouter
            support_type (str): type de support classic or spring
        """
        support = self._data["supports"][support_type].get(support_id, None)
        if support is None:
            raise ValueError(f"Appui {support_id} non trouvé. Le support n'a pas été ajouter au modèle.")
        if support_type == "classic":
            self._model.def_support(
                support["Noeud"],
                support["DX"],
                support["DY"],
                support["DZ"],
                support["RX"],
                support["RY"],
                support["RZ"],
            )
        else:
            limit_val = {"Aucune limitation": None, "Tension uniquement": "+", "Compression uniquement": "-"} 
            self._model.def_support_spring(
                support["Noeud"],
                support["Dof"],
                support["Raideur"].value * 10**3,
                limit_val[support["Limite de direction"]],
            )
        return support_id

    def _add_supports_to_model(self):
        for support_type, supports in self._data["supports"].items():
            for support_id, support in supports.items():
                self._add_support_to_model(support_id, support_type)
        return "Appuis ajoutés"

    def get_all_supports(self) -> dict:
        """Retourne l'ensemble des appuis du modèle.

        Returns:
            dict: Dictionnaire contenant les appuis classiques et ressorts :
                - "classic": dictionnaire des appuis classiques
                - "spring": dictionnaire des appuis ressorts
        """
        return self._data["supports"]

    ################## Chargements ##################

    def _convert_pos(self, pos_index: int | str, member_id: str) -> int:
        """Convertit une position textuelle ou relative en valeur numérique absolue.

        Transforme les notations symboliques en position réelle sur la barre
        pour le calcul des charges distribuées et ponctuelles.

        Args:
            pos_index (int | str): Position de la charge. Formats acceptés :
                - "start" : début de la barre (0 mm)
                - "end" : fin de la barre (longueur totale)
                - "middle" : milieu de la barre (L/2)
                - "XX%" : pourcentage de la longueur (ex: "25%")
                - int : position directe en millimètres
            member_id (str): Identifiant de la barre concernée.

        Returns:
            int: Position numérique absolue en millimètres sur la barre.

        Exemples:
            >>> _convert_pos("start", "M1") -> 0
            >>> _convert_pos("50%", "M1") -> 2500 (sur une barre de 5000 mm)
            >>> _convert_pos(1000, "M1") -> 1000
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
        """Ajoute une charge répartie linéairement sur une barre.

        Crée une charge distribuée (poids propre, neige, vent, etc.) appliquée
        sur une portion ou la totalité d'une barre. La charge peut être
        uniforme (start_load = end_load) ou trapézoïdale.

        Args:
            member_id (str): Identifiant de la barre à charger (ex: "M1").
            name (str): Nom descriptif de la charge (ex: "Neige", "PP poutre").
            start_load (float): Intensité de charge au point de départ en kN/m.
                Valeur positive selon le sens de l'axe choisi.
            end_load (float): Intensité de charge au point d'arrivée en kN/m.
                Pour une charge uniforme, end_load = start_load.
            start_pos (str, optional): Position de début de la charge.
                Formats acceptés : "start" (début), "middle" (milieu),
                "XX%" (pourcentage), ou valeur numérique en mm.
                Defaults to "start" (toute la barre).
            end_pos (str, optional): Position de fin de charge. Mêmes formats.
                Defaults to "end" (toute la barre).
            action (str): Type d'action selon DICO_COMBI_ACTION.
                Ex: "Permanente G", "Neige normale Sn", "Vent pression W+".
            direction (str): Direction et sens de l'effort dans le repère local.
                Forces : "Fx", "Fy", "Fz" (local) ou "FX", "FY", "FZ" (global).
                Defaults to "Fy" (vertical vers le bas en local).
            comment (str, optional): Commentaire descriptif sur la charge.

        Returns:
            dict: Configuration de la charge créée avec son identifiant généré.

        Note:
            La convention de signe suit Pynite : Fy négatif = charge vers le bas
            pour une barre horizontale.
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
        """Ajoute une charge ponctuelle (force ou moment) sur une barre.

        Crée une force concentrée ou un moment concentré appliqué à une position
        précise sur une barre (charge en trémie, point d'application d'une poutre,
        moment d'encastrement équivalent, etc.).

        Args:
            member_id (str): Identifiant de la barre à charger (ex: "M1").
            name (str): Nom descriptif de la charge (ex: "Charge ponctuelle P1").
            load (float): Valeur de l'effort.
                - Force : en kN si direction commence par "F"
                - Moment : en kN·m si direction commence par "M"
            pos (str, optional): Position de la charge sur la barre.
                Formats acceptés : "start" (début), "end" (fin), "middle" (milieu),
                "XX%" (pourcentage), ou valeur numérique en mm.
                Defaults to "middle".
            action (str): Type d'action selon DICO_COMBI_ACTION.
                Ex: "Permanente G", "Exploitation Q".
            direction (str): Direction et type d'effort.
                Forces : "Fx", "Fy", "Fz" (local) ou "FX", "FY", "FZ" (global).
                Moments : "Mx", "My", "Mz" (local) ou "MX", "MY", "MZ" (global).
                Defaults to "Fy".
            comment (str, optional): Commentaire descriptif sur la charge.

        Returns:
            dict: Configuration de la charge ponctuelle créée.
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
        """Ajoute plusieurs charges en batch depuis une liste.

        Permet de créer rapidement plusieurs charges similaires à partir d'une
        liste de paramètres, utile pour les modèles répétitifs ou les imports.

        Args:
            list_loads (list): Liste de tuples contenant les arguments de charge.
                Chaque tuple doit correspondre aux arguments de create_dist_load
                ou create_point_load selon le type_load.
                Ex: [("M1", "PP", -0.5, -0.5, "start", "end", "Permanente G", "FY"), ...]
            type_load (str, optional): Type de charges dans la liste.
                "Distribuée" : utilise create_dist_load pour chaque élément.
                "Autre" : utilise create_point_load pour chaque élément.
                Defaults to "Distribuée".

        Returns:
            None: Les charges sont ajoutées directement au modèle.

        Exemple:
            >>> charges = [
            ...     ("M1", "Neige", -1.5, -1.5, "start", "end", "Neige normale Sn", "FY"),
            ...     ("M2", "Neige", -1.5, -1.5, "start", "end", "Neige normale Sn", "FY"),
            ... ]
            >>> model.create_load_by_list(charges, "Distribuée")
        """
        for load in list_loads:
            if type_load == "Distribuée":
                self.create_dist_load(*load)
            else:
                self.create_point_load(*load)

    def del_load(self, load_id: str):
        """Supprime une charge par son identifiant.

        Args:
            load_id (str): Identifiant de la charge à supprimer (ex: "L1").

        Returns:
            dict: Dictionnaire de la charge supprimée.
        """
        return self._data["loads"].pop(load_id)

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

    def get_all_loads(self) -> dict:
        """Retourne l'ensemble des charges définies dans le modèle.

        Returns:
            dict: Dictionnaire de toutes les charges, indexé par leurs identifiants
                ("L1", "L2", etc.). Chaque entrée contient les propriétés complètes
                de la charge (type, intensité, position, barre concernée).
        """
        return self._data["loads"]

    def get_member_loads(self, member_id: str) -> list:
        """Retourne les charges appliquées sur une barre spécifique.

        Filtre l'ensemble des charges du modèle pour ne retourner que celles
        appliquées sur la barre identifiée.

        Args:
            member_id (str): Identifiant de la barre (ex: "M1").

        Returns:
            list: Liste des dictionnaires de charges appliquées à cette barre.

        Exemple:
            >>> model.get_member_loads("M1")
            [{'N° barre': 'M1', 'Nom': 'PP', 'Action': 'Permanente G', ...}, ...]
        """
        return [
            load
            for load in self._data["loads"].values()
            if load["N° barre"] == member_id
        ]

    def generate_model(self):
        """Génère et assemble le modèle MEF complet dans Pynite.

        Cette méthode interne (appelée automatiquement par Combinaison) :
        1. Crée l'objet FEModel3D de Pynite
        2. Ajoute tous les nœuds, matériaux, sections
        3. Ajoute toutes les barres avec leurs relâchements
        4. Ajoute tous les appuis (classiques et ressorts)
        5. Ajoute tous les chargements

        Note:
            Cette méthode est généralement appelée automatiquement par la
            classe Combinaison. Elle n'a pas besoin d'être invoquée manuellement
            dans un flux de travail standard.
        """
        from Pynite import FEModel3D
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

class Wood_beam_model(Model_generator):
    """Générateur de modèle MEF simplifié pour poutres bois continues.

    Cette classe spécialisée crée automatiquement un modèle de poutre
    sur appuis multiples avec les hypothèses suivantes :
    - Barre isostatique
    - Distance entre appuis égale (pas de porte-à-faux)
    - Possibilité d'inclinaison et de dévers

    Elle est particulièrement adaptée pour :
    - Poutres, pannes, chevrons
    - Poteaux verticaux (inclinaison = 90°)
    - Pannes à dévers (dévers ≠ 0)

    Flux de travail :
        1. Créer le modèle avec Wood_beam_model
        2. Ajouter les charges avec create_dist_load / create_point_load
        3. Passer à la classe Combinaison pour les combinaisons
        4. Analyser les résultats avec Model_result

    Note:
        Cette classe surcharge certaines méthodes de Model_generator pour
        gérer automatiquement le poids propre des barres multiples.
    """

    def __init__(
        self,
        longueur: si.mm,
        b: si.mm,
        h: si.mm,
        section: str=Model_generator.LIST_SECTION,
        classe: str = Model_generator.CLASSE_WOOD,
        poids_propre: bool = ("True", "False"),
        nbr_appuis: int = 2,
        l_appuis: float = 0,
        devers: float = 0,
        inclinaison: float = 0,
        *args,
        **kwargs,
    ):
        """Initialise un modèle de poutre bois simplifié.

        Crée automatiquement les nœuds, barres et appuis répartis uniformément
        sur la longueur totale. Génère autant de barres que de travées
        (nbr_appuis - 1).

        Args:
            longueur (si.mm): Longueur totale de la barre en millimètres.
            b (si.mm): Épaisseur de la section en millimètres (largeur).
            h (si.mm): Hauteur de la section en millimètres.
            section (str, optional): Type de section transversale.
                "Rectangulaire" ou "Circulaire". Defaults to "Rectangulaire".
            classe (str, optional): Classe de résistance du bois selon l'EC5.
                Ex: "C24", "GL28h". Defaults to "C24".
            poids_propre (bool, optional): Générer automatiquement le poids
                propre de chaque barre comme charge répartie. Defaults to True.
            nbr_appuis (int, optional): Nombre d'appuis rotulés (≥ 2).
                2 = poutre simple sur deux appuis, 3 = poutre continue sur 3 appuis.
                Defaults to 2.
            l_appuis (float, optional): Longueur d'appui en mm pour la vérification
                à la compression perpendiculaire. Defaults to 0.
            devers (float, optional): Angle de dévers en degrés (rotation autour
                de l'axe longitudinal X). Utilisé pour les pannes inclinées
                transversalement. Defaults to 0.
            inclinaison (float, optional): Angle d'inclinaison en degrés (rotation
                autour de l'axe Z). 0° = poutre horizontale, 90° = poteau vertical.
                Utile pour les chevrons ou poteaux inclinés. Defaults to 0.
        """
        super().__init__(*args, **kwargs)
        self.longueur = longueur * si.mm
        self.b = b * si.mm
        self.h = h * si.mm
        self.section = section
        self.classe = classe
        self.poids_propre = poids_propre
        self.nbr_appuis = nbr_appuis
        self.l_appuis = l_appuis
        self.devers = devers
        self.inclinaison = inclinaison

        material = self.add_material_by_class(classe)
        section = self.add_section(b, h, 0, section)

        d_appuis = longueur / (nbr_appuis - 1)
        for i in range(nbr_appuis):
            node = self.add_node(int(float(i * d_appuis * mt.cos(mt.radians(inclinaison)))), int(float(i * d_appuis * mt.sin(mt.radians(inclinaison)))), 0)
            if not i:
                self.add_support(node, DX=True, DY=True, DZ=True, RX=True, RY=False, RZ=False, l_appuis=l_appuis)
            else:
                self.add_support(node, DX=False, DY=True, DZ=True, RX=True, RY=False, RZ=False, l_appuis=l_appuis)

        for i in range(nbr_appuis-1):
            self.add_member(f"N{i+1}", f"N{i+2}", material, section, poids_propre, rotation=devers, tension_only=False, compression_only=False)

    def _convert_beam_pos(self, pos_index: int | str) -> int:
        """Converti la position en valeur recevable par la fonction create_load

        Args:
                pos_index (int | str): position de la charge sur la barre

        Returns:
                int: la position numérique sur la barre
        """
        match pos_index:
            case "start":
                pos_index = 0 * si.mm
                member_id = "M1"
            case "end":
                member_id = "M" + str(len(self._data["members"]))
                pos_index = self._data["members"][member_id]["Longueur"]
            case "middle":
                dist = self.longueur / 2
                sum_dist = 0
                for member_id in self._data["members"].keys():
                    sum_dist = sum_dist + self._data["members"][member_id]["Longueur"]
                    if sum_dist >= dist:
                        delta_end_dist = sum_dist - dist
                        pos_index = self._data["members"][member_id]["Longueur"] - delta_end_dist
                        break
                pos_index = round(pos_index, 0)
            case str(x) if "%" in x:
                pos_index = pos_index.split("%")[0]
                pos_index.replace(" ", "")
                dist = round(
                    self.longueur * int(pos_index) / 100,
                    0,
                )
                sum_dist = 0
                for member_id in self._data["members"].keys():
                    sum_dist = sum_dist + self._data["members"][member_id]["Longueur"]
                    if sum_dist >= dist:
                        delta_end_dist = sum_dist - dist
                        pos_index = self._data["members"][member_id]["Longueur"] - delta_end_dist
                        break
                pos_index = round(pos_index, 0)
            case _:
                dist = float(pos_index) * si.mm
                sum_dist = 0
                for member_id in self._data["members"].keys():
                    sum_dist = sum_dist + self._data["members"][member_id]["Longueur"]
                    if sum_dist >= dist:
                        delta_end_dist = sum_dist - dist
                        pos_index = self._data["members"][member_id]["Longueur"] - delta_end_dist
                        break
                pos_index = round(pos_index, 0)
        return pos_index.value*10**3, member_id

    def create_dist_load(
        self,
        name: str,
        start_load: float,
        end_load: float,
        start_pos: str = None,
        end_pos: str = None,
        action: str = Model_generator.ACTION,
        direction: str = ("Fx", "Fy", "Fz", "FX", "FY", "FZ"),
        comment: str = None,
    ):
        """Ajoute une charge distribuée sur la barre

        Args:
            name (str): nom de la charge.
            start_load (int): effort de départ en kN/m.
            end_load (int): effort de fin en kN/m.
            start_pos (str, optional): position de début de la charge sur la barre en mm. En complément il est possible de mettre "start", "middle"
                                        ou un pourcentage pour définir la position de la charge.
            end_pos (str, optional): position de début de la charge sur la barre en mm. En complément il est possible de mettre "end", "middle"
                                        ou un pourcentage pour définir la position de la charge.
            action (str): type d'action de l'effort.
            direction (str): sens de l'effort sur la barre.
            comment (str, optional): commentaire sur la charge.
        """
        if not start_pos:
            start_pos = "start"
        if not end_pos:
            end_pos = "end"
        start_pos, member_id_start = self._convert_beam_pos(start_pos)
        end_pos, member_id_end = self._convert_beam_pos(end_pos)
        if member_id_start != member_id_end:
            id_start = int(member_id_start.split("M")[1])
            id_end = int(member_id_end.split("M")[1])
            super().create_dist_load(f"M{id_start}", name, start_load, end_load, start_pos, 'end', action, direction, comment)
            if id_start != id_end-1:
                for i in range(id_start, id_end-1):
                    super().create_dist_load(f"M{i+1}", name, start_load, end_load, 'start', 'end', action, direction, comment)
            super().create_dist_load(f"M{id_end}", name, start_load, end_load, 'start', end_pos, action, direction, comment)
        else:
            super().create_dist_load(member_id_start, name, start_load, end_load, start_pos, end_pos, action, direction, comment)

    def create_point_load(
        self,
        name: str,
        load: float,
        pos: str = None,
        action: str = Model_generator.ACTION, 
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
            comment: str = None
    ):
        """Ajoute une charge nodale sur la barre

        Args:
            name (str): nom de la charge.
            load (int): effort de départ en kN ou kN.m.
            pos (str, optional): position de la charge sur la barre en mm. En complément il est possible de mettre "start", "middle", "end"
                                    ou un pourcentage pour définir la position de la charge.
            action (str): type d'action de l'effort.
            direction (str): sens de l'effort sur la barre.
            comment (str, optional): commentaire sur la charge.
        """
        pos, member_id = self._convert_beam_pos(pos)
        return super().create_point_load(member_id, name, load, pos, action, direction, comment)
        


class Model_result(Projet):
    """Classe d'analyse et d'extraction des résultats MEF.

    Cette classe permet de lancer l'analyse aux éléments finis et de récupérer
    les résultats (efforts, déformées, réactions) du modèle généré par
    Model_generator et résolu par Combinaison.

    Elle hérite de Projet pour maintenir le contexte du projet.
    L'analyse est lancée automatiquement lors de l'instanciation.

    Flux de travail :
        1. Créer un modèle avec Model_generator
        2. Définir les combinaisons avec Combinaison
        3. Instancier Model_result pour analyser et extraire les résultats
        4. Utiliser les méthodes get_* ou show_* pour visualiser

    Types d'analyse disponibles:
        - "Général": Analyse complète (linéaire + P-Delta si nécessaire)
        - "Linéaire": Analyse linéaire élastique uniquement
        - "Second ordre": Analyse P-Delta (effets du second ordre)
    """

    ANALYZE_TYPE = ("Général", "Linéaire", "Second ordre")

    def __init__(
        self,
        model_generator: object,
        analyze_type: str = ANALYZE_TYPE,
        check_stability: bool = ("False", "True"),
        *args,
        **kwargs,
    ):
        """Initialise l'analyseur de résultats MEF.

        Lance automatiquement l'analyse du modèle après initialisation.
        Vérifie que tous les nœuds sont connectés avant de résoudre.

        Args:
            model_generator (Model_generator): Instance du modèle généré
                contenant les nœuds, barres, matériaux et chargements.
            analyze_type (str): Type d'analyse à réaliser.
                "Général", "Linéaire" ou "Second ordre". Defaults to "Général".
            check_stability (bool, optional): Active la vérification de stabilité
                du modèle. Ralentit le calcul, utile pour le débogage.
                Defaults to False.
            *args: Arguments transmis à la classe parent Projet.
            **kwargs: Arguments nommés transmis à Projet.

        Raises:
            ValueError: Si des nœuds orphelins (non connectés) sont détectés.
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
        """Génère un graphique matplotlib de base pour les résultats.

        Méthode interne utilisée par show_internal_force_of_member et
        show_deflection_of_member pour créer des diagrammes cohérents.

        Args:
            title (str): Titre principal du graphique.
            combo_name (str): Nom de la combinaison affichée (sous-titre).
            x_values (array): Valeurs pour l'axe horizontal (position le long de la barre).
            y_values (array): Valeurs pour l'axe vertical (effort ou déplacement).
            x_label (str): Label de l'axe X.
            y_label (str): Label de l'axe Y avec unité.
            color (str): Couleur matplotlib du tracé (ex: "r", "b", "g", "orange").
            fill_between (bool, optional): Remplit l'aire sous la courbe si True.
                Defaults to True.
            savefig (bool, optional): Sauvegarde automatique si True, affichage interactif sinon.
                Defaults to False.
            filepath (str, optional): Chemin de sauvegarde si savefig=True.
                Ouvre une boîte de dialogue si None.

        Returns:
            str: Chemin du fichier sauvegardé si savefig=True.
            None: Affiche le graphique si savefig=False.

        Note:
            Cette méthode est interne et ne devrait pas être appelée directement.
        """
        # plt.clf()  # Effacer le graphique précédent
        plt.figure(self.name, figsize=(11, 4))
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
                from PySide6.QtWidgets import QFileDialog
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
            plt.savefig(filepath)
            return filepath
        else:
            plt.show()
        

    def _analyze(self):
        """Lance l'analyse du modèle aux éléments finis.

        Méthode interne appelée automatiquement lors de l'instanciation.
        Détecte les nœuds orphelins et lance le solveur adapté au type d'analyse.

        Raises:
            RuntimeError: Si des nœuds orphelins sont détectés dans le modèle.
        """
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

    def get_member_length(self, member_id: str) -> float:
        """Retourne la longueur d'une barre par son identifiant.

        Args:
            member_id (str): Identifiant de la barre (ex: "M1").

        Returns:
            float: Longueur de la barre en millimètres avec unité (si.mm).
        """
        return self._model_generator._data["members"][member_id]["Longueur"]

    def get_internal_force(
        self,
        member_id: str,
        combination: str,
        type: str = ("Nx", "Vy", "Vz", "Mx", "My", "Mz"),
        n_points: int = 20,
    ) -> np.array:
        """Retourne les efforts internes le long d'une barre pour une combinaison donnée.

        Extrait les valeurs d'efforts (effort normal, cisaillement, moment)
        répartis le long de la barre avec un nombre de points défini.

        Args:
            member_id (str): Identifiant de la barre (ex: "M1").
            combination (str): Nom de la combinaison de charges à analyser
                (ex: "ELU_STR", "ELS_QP").
            type (str): Type d'effort interne à récupérer.
                "Nx" = effort normal, "Vy"/"Vz" = effort tranchant,
                "Mx" = moment de torsion, "My"/"Mz" = moment fléchissant.
                Defaults to "Nx".
            n_points (int, optional): Nombre de points de discrétisation
                le long de la barre. Plus ce nombre est élevé, plus la courbe
                est lisse. Defaults to 20.

        Returns:
            tuple: (positions, valeurs) où :
                - positions (array): Abscisses le long de la barre en mm
                - valeurs (array): Valeurs de l'effort correspondant

        Note:
            Les valeurs sont retournées dans les unités de base de Pynite (N ou N·mm).
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

    def get_min_max_internal_force(self, member_id: str|list, combination: str|list) -> dict:
        """Retourne les valeurs minimales et maximales des efforts internes.

        Analyse tous les types d'efforts (Nx, Vy, Vz, Mx, My, Mz) et retourne
        les valeurs extrêmes avec la combinaison correspondante. Supporte les
        barres continues (liste de barres) et les tags de combinaisons.

        Args:
            member_id (str | list): Identifiant de la barre ou liste de barres
                pour une barre continue (ex: "M1" ou ["M1", "M2", "M3"]).
            combination (str | list): Nom de la combinaison ou liste de tags
                de combinaisons (ex: "ELU_STR" ou ["ELU_STR", "ELU_ACC"]).
                Si des tags sont fournis, le max/min est cherché parmi toutes
                les combinaisons correspondantes.

        Returns:
            dict: Dictionnaire structuré par type d'effort :
                {
                    "Nx": {"Min": (valeur, combinaison), "Max": (valeur, combinaison)},
                    "Vy": {"Min": (...), "Max": (...)}, ...
                }
                Les valeurs incluent les unités (N pour efforts, N·mm pour moments).

        Raises:
            ValueError: Si une barre de la liste n'existe pas dans le modèle.

        """
        dict_internal_forces = {}
        for type in ("Nx", "Vy", "Vz", "Mx", "My", "Mz"):
            if "[" in member_id:
                import json
                member_id = json.loads(member_id.replace("'", "\""))
            elif isinstance(member_id, str):
                member_id = [member_id]
            max_value = 0
            min_value = 0
            combi_max = None
            combi_min = None
            for member in member_id:
                if member not in self._model_generator._model.members:
                    raise ValueError(f"La membrure {member} n'est pas dans le model MEF")
                
                if type == "Nx":
                    max = self._model_generator._model.members[member].max_axial(
                        combo_tags=combination
                    )
                    min = self._model_generator._model.members[member].min_axial(
                        combo_tags=combination
                    )
                elif type == "Vy":
                    max = self._model_generator._model.members[member].max_shear(
                        "Fy", combo_tags=combination
                    )
                    min = self._model_generator._model.members[member].min_shear(
                        "Fy", combo_tags=combination
                    )
                elif type == "Vz":
                    max = self._model_generator._model.members[member].max_shear(
                        "Fz", combo_tags=combination
                    )
                    min = self._model_generator._model.members[member].min_shear(
                        "Fz", combo_tags=combination
                    )
                elif type == "Mx":
                    max = self._model_generator._model.members[member].max_torque(
                        combo_tags=combination
                    )
                    min = self._model_generator._model.members[member].min_torque(
                        combo_tags=combination
                    )
                elif type == "My":
                    max = self._model_generator._model.members[member].max_moment(
                        "My", combo_tags=combination
                    )
                    min = self._model_generator._model.members[member].min_moment(
                        "My", combo_tags=combination
                    )
                elif type == "Mz":
                    max = self._model_generator._model.members[member].max_moment(
                        "Mz", combo_tags=combination
                    )
                    min = self._model_generator._model.members[member].min_moment(
                        "Mz", combo_tags=combination
                    )
                if isinstance(max, tuple):
                    if max[0] > max_value:
                        max_value = max[0]
                        combi_max = max[1]
                    if min[0] < min_value:
                        min_value = min[0]
                        combi_min = min[1]
                else:
                    if max > max_value:
                        max_value = max
                        combi_max = combination
                    if min < min_value:
                        min_value = min
                        combi_min = combination
            
            if "M" in type:
                si_unit = si.N * si.mm
            else:
                si_unit = si.N
            dict_internal_forces[type] = {"Min": (min_value * si_unit, combi_min), "Max": (max_value * si_unit, combi_max)}
        return dict_internal_forces
    
    def get_absolute_internal_force(self, member_id: str|list, combination: str|list, type: str = ("Nx", "Vy", "Vz", "Mx", "My", "Mz"), get_combo_name: bool = ("False", "True")):
        """Retourne la valeur absolue maximale d'un type d'effort spécifique.

        Détermine automatiquement si le maximum absolu est la valeur minimale
        (la plus négative) ou maximale (la plus positive) et retourne sa
        valeur absolue. Utile pour les vérifications de dimensionnement.

        Args:
            member_id (str | list): Identifiant de la barre ou liste de barres.
            combination (str | list): Nom de la combinaison ou tags de combinaisons.
            type (str): Type d'effort à analyser ("Nx", "Vy", "Vz", "Mx", "My", "Mz").
                Defaults to "Nx".
            get_combo_name (bool, optional): Si True, retourne également le nom
                de la combinaison associée. Defaults to False.

        Returns:
            float | dict: Valeur absolue de l'effort maximal.
                Si get_combo_name=True, retourne un dict :
                {"Effort": valeur, "Combinaison": nom_combinaison}
        """
        ei = self.get_min_max_internal_force(member_id, combination)
        max = "Max"
        if abs(ei[type]["Min"][0]) > ei[type]["Max"][0]:
            max = "Min"
        if get_combo_name:
            return {"Effort": abs(ei[type][max][0]), "Combinaison": ei[type][max][1]}
        else:
            return abs(ei[type][max][0])
    

    def show_internal_force_of_member(
        self,
        member_id: str|list,
        combination: str,
        type: str = ("Nx", "Vy", "Vz", "Mx", "My", "Mz"),
        n_points: int = 20,
        screenshot: bool = ("False", "True"),
        filepath: str=None,
    ):
        """Affiche le diagramme des efforts internes d'une barre.

        Génère un graphique matplotlib montrant la distribution de l'effort
        interne choisi le long de la barre (ou des barres continues).
        Les valeurs sont automatiquement converties en kN ou kN·m.

        Args:
            member_id (str | list): Identifiant de la barre ou liste de barres
                pour une barre continue.
            combination (str): Nom de la combinaison à afficher (ex: "ELU_STR").
            type (str): Type d'effort à diagrammer.
                - "Nx": Effort normal (orange)
                - "Vy", "Vz": Effort tranchant (bleu)
                - "Mx", "My", "Mz": Moments (rouge)
                Defaults to "Nx".
            n_points (int, optional): Nombre de points pour le tracé.
                Defaults to 20.
            screenshot (bool, optional): Si True, sauvegarde l'image sans
                afficher la fenêtre interactive. Defaults to False.
            filepath (str, optional): Chemin de sauvegarde si screenshot=True.
                Ouvre une boîte de dialogue si None.

        Returns:
            str: Chemin du fichier sauvegardé si screenshot=True.
            None: Affiche le graphique interactif si screenshot=False.

        Note:
            Les barres sont concaténées pour former un diagramme continu
            dans le cas d'une barre continue (ex: ["M1", "M2", "M3"]).
        """
        x_label = "Longueur (mm)"
        if "[" in member_id:
            import json
            member_id = json.loads(member_id.replace("'", "\""))
        elif isinstance(member_id, str):
            member_id = [member_id]
        x_value = []
        y_value = []
        for member in member_id:
            if member not in self._model_generator._model.members:
                raise ValueError(f"La membrure {member} n'est pas dans le model MEF")
            x_local, y_local = self.get_internal_force(
                member, combination, type, n_points
            )
            if len(x_value) == 0:
                x_value = x_local
                y_value = y_local
            else:
                x_local = x_local + x_value[-1]
                x_value = np.concatenate((x_value, x_local))
                y_value = np.concatenate((y_value, y_local))
            
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
        """Retourne les déplacements/déformées le long d'une barre.

        Extrait les valeurs de déplacement aux nœuds intermédiaires de la barre
        pour une direction donnée et une combinaison spécifique.

        Args:
            member_id (str): Identifiant de la barre (ex: "M1").
            combination (str): Nom de la combinaison à analyser
                (ex: "ELS_QP", "W_inst_Q").
            direction (str): Direction du déplacement.
                "dx" = selon l'axe longitudinal (allongement/raccourcissement)
                "dy" = flèche dans le plan vertical local
                "dz" = flèche dans le plan horizontal local
                Defaults to "dy" (flèche verticale usuelle).
            n_points (int, optional): Nombre de points de discrétisation.
                Defaults to 20.

        Returns:
            tuple: (positions, déplacements) où :
                - positions (array): Abscisses le long de la barre en mm
                - déplacements (array): Valeurs de déplacement en mm

        Note:
            Les déplacements sont relatifs aux conditions d'appui (la ligne
            élastique est calculée par rapport à la position déformée des appuis).
        """
        return self._model_generator._model.members[member_id].deflection_array(
            direction, n_points=n_points, combo_name=combination
        )

    def get_min_max_deflection(self, member_id: str|list, combination: str|list) -> dict:
        """Retourne les déplacements minimaux et maximaux d'une barre.

        Analyse les trois directions de déplacement (dx, dy, dz) et retourne
        les valeurs extrêmes avec la combinaison correspondante. Supporte les
        barres continues et les tags de combinaisons ELS.

        Args:
            member_id (str | list): Identifiant de la barre ou liste de barres
                pour une barre continue.
            combination (str | list): Nom de la combinaison ou tags ELS
                (ex: "ELS_QP", ["W_inst_Q", "W_net_fin"]).

        Returns:
            dict: Dictionnaire structuré par direction :
                {
                    "dx": {"Min": (valeur, combinaison), "Max": (valeur, combinaison)},
                    "dy": {"Min": (...), "Max": (...)},
                    "dz": {"Min": (...), "Max": (...)}
                }
                Les valeurs incluent l'unité (si.mm).

        Note:
            Les valeurs "Min" peuvent être négatives (déplacement dans le sens
            opposé à l'axe). La flèche maximale est généralement la valeur absolue
            la plus grande en valeur absolue dans la direction verticale (dy).
        """
        dict_deflection = {}
        for type in ("dx", "dy", "dz"):
            if "[" in member_id:
                import json
                member_id = json.loads(member_id.replace("'", "\""))
            elif isinstance(member_id, str):
                member_id = [member_id]
            max_value = 0
            min_value = 0
            combi_max = None
            combi_min = None
            for member in member_id:
                if member not in self._model_generator._model.members:
                    raise ValueError(f"La membrure {member} n'est pas dans le model MEF")
                max = self._model_generator._model.members[member].max_deflection(
                    type, combo_tags=combination
                )

                min = self._model_generator._model.members[member].min_deflection(
                    type, combo_tags=combination
                )
                if isinstance(max, tuple):
                    if max[0] > max_value:
                        max_value = max[0]
                        combi_max = max[1]
                    if min[0] < min_value:
                        min_value = min[0]
                        combi_min = min[1]
                else:
                    if max > max_value:
                        max_value = max
                        combi_max = combination
                    if min < min_value:
                        min_value = min
                        combi_min = combination

            dict_deflection[type] = {"Min": (min_value * si.mm, combi_min), "Max": (max_value * si.mm, combi_max)}
        return dict_deflection
    
    def get_absolute_max_deflection(self, member_id: str|list, combination: str|list, direction: str = ("dx", "dy", "dz"), get_combo_name: bool=("False", "True")):
        """Retourne la valeur absolue maximale de déplacement pour une direction donnée.

        Détermine automatiquement si le maximum absolu est la valeur minimale
        ou maximale et retourne sa valeur absolue. Méthode standard pour
        récupérer la flèche maximale à vérifier.

        Args:
            member_id (str | list): Identifiant de la barre ou liste de barres.
            combination (str | list): Nom de la combinaison ou tags de combinaisons.
                (ex: "ELS_QP", ["W_inst_Q", "W_net_fin"]).
            direction (str): Direction à analyser ("dx", "dy", "dz").
                Defaults to "dy" (flèche verticale).
            get_combo_name (bool, optional): Si True, retourne également le nom
                de la combinaison associée. Defaults to False.

        Returns:
            float | dict: Valeur absolue du déplacement maximal en mm.
                Si get_combo_name=True, retourne un dict :
                {"Flèche": valeur, "Combinaison": nom_combinaison}
        """
        deflection = self.get_min_max_deflection(member_id, combination)
        max = "Max"
        if abs(deflection[direction]["Min"][0]) > deflection[direction]["Max"][0]:
            max = "Min"
        if get_combo_name:
            return {"Flèche": abs(deflection[direction][max][0]), "Combinaison": deflection[direction][max][1]}
        else:
            return abs(deflection[direction][max][0])
        

    def show_deflection_of_member(
        self,
        member_id: str|list,
        combination: str,
        direction: str = ("dx", "dy", "dz"),
        n_points: int = 20,
        screenshot: bool = ("False", "True"),
        filepath: str=None,
    ):
        """Affiche le diagramme de déformée (flèche) d'une barre.

        Génère un graphique matplotlib montrant la ligne élastique de la barre
        pour la direction et combinaison spécifiées. Utile pour visualiser
        la déformée et identifier les zones de flèche maximale.

        Args:
            member_id (str | list): Identifiant de la barre ou liste de barres
                pour une barre continue.
            combination (str): Nom de la combinaison à afficher (ex: "W_inst_Q").
            direction (str): Direction de la déformée à afficher.
                "dx" = allongement, "dy" = flèche verticale, "dz" = flèche horizontale.
                Defaults to "dy".
            n_points (int, optional): Nombre de points pour le tracé.
                Defaults to 20.
            screenshot (bool, optional): Si True, sauvegarde l'image.
                Defaults to False.
            filepath (str, optional): Chemin de sauvegarde si screenshot=True.

        Returns:
            str: Chemin du fichier sauvegardé si screenshot=True.
            None: Affiche le graphique interactif si screenshot=False.

        Note:
            La déformée est tracée en vert avec remplissage pour visualiser
            l'amplitude des déplacements.
        """
        title = f'Barre {member_id}: Flèche {direction}'
        x_label = "Longueur (mm)"
        y_label = "Déplacement\n(mm)"
        color = "g"
        if "[" in member_id:
            import json
            member_id = json.loads(member_id.replace("'", "\""))
        elif isinstance(member_id, str):
            member_id = [member_id]
        x_value = []
        y_value = []
        for member in member_id:
            if member not in self._model_generator._model.members:
                raise ValueError(f"La membrure {member} n'est pas dans le model MEF")
            x_local, y_local = self.get_deflection(
                member, combination, direction, n_points
            )
            if len(x_value) == 0:
                x_value = x_local
                y_value = y_local
            else:
                x_local = x_local + x_value[-1]
                x_value = np.concatenate((x_value, x_local))
                y_value = np.concatenate((y_value, y_local))

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
        show_loads: bool = ("True", "False"),
        diagrams: str = ("Aucun", "Fx", "Fy", "Fz", "My", "Mz", "Tx", "Flèche"),
        scale: int = 1000,
        screenshot: bool = ("False", "True"),
        filepath: str=None,
    ):
        """Affiche une visualisation 3D interactive du modèle MEF.

        Rendu 3D complet du modèle avec possibilité d'afficher les charges,
        les efforts internes colorés, ou la déformée. Utilise le renderer
        Pynite avec interaction utilisateur.

        Args:
            combination (str): Nom de la combinaison à afficher.
            annotation_size (int, optional): Taille des annotations de nœuds.
                Defaults to 70.
            show_loads (bool, optional): Affiche les charges appliquées.
                Defaults to True.
            diagrams (str, optional): Type de diagramme à superposer.
                - "Aucun": Modèle fil de fer seul
                - "Fx", "Fy", "Fz", "My", "Mz", "Tx": Diagramme d'efforts coloré
                - "Flèche": Déformée amplifiée
                Defaults to "Aucun".
            scale (int, optional): Facteur d'échelle pour les diagrammes et
                la déformée. Defaults to 1000.
            screenshot (bool, optional): Si True, capture l'image sans interaction.
                Defaults to False.
            filepath (str, optional): Chemin de sauvegarde. Ouvre une boîte de
                dialogue si None et screenshot=True.

        Returns:
            str: Chemin du fichier sauvegardé si screenshot=True.
            None: Affiche la fenêtre interactive si screenshot=False.

        Note:
            En mode interactif (screenshot=False), utilisez la souris pour
            tourner le modèle (clic gauche), zoomer (molette), paner (clic droit).
            Appuyez sur Q pour fermer.
        """
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        from Pynite.Rendering import Renderer
        renderer = Renderer(self._model_generator._model)
        renderer.combo_name = combination
        renderer.annotation_size = annotation_size
        renderer.render_loads = show_loads
        if diagrams == "Flèche":
            renderer.deformed_shape = True
            renderer.deformed_scale = scale
        elif diagrams != "Aucun":
            renderer.member_diagrams = diagrams
            renderer.diagram_scale = scale
            
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
            renderer.screenshot(filepath, interact=interaction, reset_camera=True)
            return filepath
        else:
            renderer.render_model()

    def get_global_displacement_of_node(self, node_id: str, combination: str) -> dict:
        """Retourne les déplacements/rotations globaux d'un nœud.

        Extrait les valeurs de déplacement et rotation dans le repère global
        (X, Y, Z) pour un nœud spécifique et une combinaison donnée.

        Args:
            node_id (str): Identifiant du nœud (ex: "N1").
            combination (str): Nom de la combinaison à analyser.

        Returns:
            dict: Déplacements et rotations du nœud :
                {
                    "DX", "DY", "DZ": translations en mm (si.mm)
                    "RX", "RY", "RZ": rotations en radians (sans unité)
                }

        Note:
            Les déplacements sont relatifs à la position initiale du nœud.
            Les rotations positives suivent la convention de la main droite.
        """
        DX = self._model_generator._model.nodes[node_id].DX[combination] * si.mm
        DY = self._model_generator._model.nodes[node_id].DY[combination] * si.mm
        DZ = self._model_generator._model.nodes[node_id].DZ[combination] * si.mm
        RX = self._model_generator._model.nodes[node_id].RX[combination]
        RY = self._model_generator._model.nodes[node_id].RY[combination]
        RZ = self._model_generator._model.nodes[node_id].RZ[combination]
        return {"DX": DX, "DY": DY, "DZ": DZ, "RX": RX, "RY": RY, "RZ": RZ}

    def _get_node_reaction(self, node_id: str, combination: str) -> dict:
        """Retourne les réactions d'appui en un nœud (méthode interne).

        Extrait les forces et moments de réaction dans le repère global
        pour un nœud d'appui. Cette méthode est principalement utilisée
        par get_supports_reactions.

        Args:
            node_id (str): Identifiant du nœud d'appui (ex: "N1").
            combination (str): Nom de la combinaison à analyser.

        Returns:
            dict: Réactions au nœud :
                {
                    "FX", "FY", "FZ": forces en N (si.N)
                    "MX", "MY", "MZ": moments en N·mm (si.N*si.mm)
                }

        Note:
            Cette méthode est interne. Pour récupérer toutes les réactions
            d'appui, utilisez plutôt get_supports_reactions().
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

    def get_supports_reactions(self, combination: str) -> dict:
        """Retourne les réactions d'appui pour tous les appuis du modèle.

        Parcourt tous les appuis définis dans le modèle et récupère les
        forces et moments de réaction pour la combinaison spécifiée.

        Args:
            combination (str): Nom de la combinaison à analyser (ex: "ELU_STR").

        Returns:
            dict: Dictionnaire indexé par identifiant d'appui :
                {
                    "S1": {"FX": ..., "FY": ..., "FZ": ..., "MX": ..., ...},
                    "S2": {...},
                    ...
                }
                Forces en N, moments en N·mm.

        Note:
            Un appui rotulé retournera des moments nuls (MX=MY=MZ=0).
            Un glisseur retournera une force nulle dans la direction libre.
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
