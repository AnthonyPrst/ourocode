# coding in UTF-8
# by Anthony PARISOT
import os
import sys
import math as mt
import numpy as np
from PIL import Image

import forallpeople as si

from ourocode.eurocode.core.projet import Projet


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
            "structural_members": {},
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

    ################## Barre structurale (regroupement de barres FEM continues) ##################

    def _member_unit_vector(self, member_id: str, from_node: str) -> np.ndarray:
        """Retourne le vecteur unitaire d'une barre, orienté depuis ``from_node``.

        Args:
            member_id (str): identifiant de la barre FEM.
            from_node (str): noeud de départ (doit être une extrémité de la barre).

        Returns:
            np.ndarray: vecteur unitaire 3D (sans dimension).
        """
        n1, n2 = self._data["members"][member_id]["Noeuds"]
        other = n2 if n1 == from_node else n1
        c_from = self._data["nodes"][from_node]
        c_to = self._data["nodes"][other]
        v = np.array(
            [
                (c_to["X"] - c_from["X"]).value,
                (c_to["Y"] - c_from["Y"]).value,
                (c_to["Z"] - c_from["Z"]).value,
            ]
        )
        norm = np.linalg.norm(v)
        return v / norm if norm > 0 else v

    def _end_at_node(self, member_id: str, node_id: str) -> str:
        """Retourne ``"start"`` ou ``"end"`` selon l'extrémité de la barre au noeud donné."""
        n1, _ = self._data["members"][member_id]["Noeuds"]
        return "start" if node_id == n1 else "end"

    def _has_flexural_release_at(self, member_id: str, end: str) -> bool:
        """Indique si la barre a un relâchement en rotation (rotule fléchie) à l'extrémité.

        Un relâchement ``teta_y`` ou ``teta_z`` rompt la continuité de la flexion
        et donc la barre structurale.
        """
        rel = self._data["members"][member_id]["Relaxation"].get(end)
        if not rel:
            return False
        return bool(rel.get("teta_y") or rel.get("teta_z"))

    def detect_continuous_members(
        self, angle_tol_deg: float = 1.0
    ) -> list[list[str]]:
        """Détecte automatiquement les chaînes de barres FEM formant une barre continue.

        Deux barres FEM sont considérées comme continues au travers d'un noeud partagé si :

        1. Elles utilisent le même matériau et la même section.
        2. Aucune des deux n'a de relâchement en rotation (``teta_y`` ou ``teta_z``)
           à l'extrémité partagée — une rotule rompt la continuité en flexion.
        3. Elles sont colinéaires au noeud partagé à ``angle_tol_deg`` près :
           les deux vecteurs sortants du noeud partagé sont quasi-opposés
           (``cos(angle) ≈ -1``).

        Le degré du noeud n'est pas contraint : un noeud T ou Y (contrefiche fixée
        en milieu d'arbalétrier par ex.) ne rompt pas la continuité de l'arbalétrier
        tant que les barres de l'arbalétrier sont colinéaires et sans relâchement.
        La présence d'un appui intermédiaire ne rompt pas non plus la continuité.

        Args:
            angle_tol_deg (float): Tolérance angulaire en degrés pour la colinéarité.
                Valeurs typiques : 0.5° (strict) à 5° (tolérant sur CAO imprécise).
                Defaults to 1.0.

        Returns:
            list[list[str]]: Liste de chaînes ordonnées de ``member_id``.
                Chaque chaîne est orientée d'une extrémité libre vers l'autre.
                Une barre isolée (non continue) est retournée comme une chaîne
                d'un seul élément.

        Exemple:
            >>> chains = model.detect_continuous_members()
            >>> # [["M1"], ["M2", "M3", "M4"], ["M5"]]
            >>> # M2→M3→M4 forment une solive continue sur 4 appuis par exemple.

        Note:
            La méthode ne modifie pas le modèle. Pour créer effectivement les
            barres structurales à partir du résultat, utilisez ``group_members``
            ou ``auto_group_continuous_members``.
        """
        from collections import defaultdict

        # 1. Index noeud -> [(member_id, end)]
        node_to_ends: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for mid, m in self._data["members"].items():
            n1, n2 = m["Noeuds"]
            node_to_ends[n1].append((mid, "start"))
            node_to_ends[n2].append((mid, "end"))

        cos_threshold = -np.cos(np.radians(angle_tol_deg))

        def can_connect(m1: str, m2: str, shared: str) -> bool:
            mem1 = self._data["members"][m1]
            mem2 = self._data["members"][m2]
            if mem1["Matériaux"] != mem2["Matériaux"]:
                return False
            if mem1["Section"] != mem2["Section"]:
                return False
            # Ne pas exiger degré 2 au noeud : une contrefiche (T/Y) ne rompt
            # pas la continuité de l'arbalétrier si les deux barres sont
            # colinéaires et sans relâchement.
            if self._has_flexural_release_at(
                m1, self._end_at_node(m1, shared)
            ) or self._has_flexural_release_at(m2, self._end_at_node(m2, shared)):
                return False
            v1 = self._member_unit_vector(m1, shared)
            v2 = self._member_unit_vector(m2, shared)
            # Continuité : v1 et v2 sortent du noeud partagé en directions opposées
            return float(np.dot(v1, v2)) <= cos_threshold

        # 2. Construire la table des voisins : member_id -> {node_partage: voisin}
        # Au noeud partagé il peut y avoir N barres (T, Y…) : on cherche parmi
        # toutes celles qui passent le test de continuité (même matériau, même
        # section, colinéaires, sans relâchement). Au plus une sera retenue.
        neighbors: dict[str, dict[str, str]] = {}
        for mid, m in self._data["members"].items():
            n1, n2 = m["Noeuds"]
            nb: dict[str, str] = {}
            for shared in (n1, n2):
                ends = node_to_ends[shared]
                other = next(
                    (om for om, _ in ends if om != mid and can_connect(mid, om, shared)),
                    None,
                )
                if other is not None:
                    nb[shared] = other
            neighbors[mid] = nb

        # 3. Parcourir les chaînes (graphe localement linéaire)
        visited: set[str] = set()
        chains: list[list[str]] = []
        for start_mid in self._data["members"].keys():
            if start_mid in visited:
                continue
            n1, n2 = self._data["members"][start_mid]["Noeuds"]

            def walk(direction_node: str) -> list[str]:
                chain_ids: list[str] = []
                cur_mid = start_mid
                out_node = direction_node
                while True:
                    nxt = neighbors[cur_mid].get(out_node)
                    if nxt is None or nxt == start_mid or nxt in chain_ids:
                        break
                    chain_ids.append(nxt)
                    nxt_nodes = self._data["members"][nxt]["Noeuds"]
                    out_node = (
                        nxt_nodes[0] if nxt_nodes[1] == out_node else nxt_nodes[1]
                    )
                    cur_mid = nxt
                return chain_ids

            forward = walk(n2)
            backward = walk(n1)
            chain = list(reversed(backward)) + [start_mid] + forward
            visited.update(chain)
            chains.append(chain)
        return chains

    def group_members(
        self,
        name: str,
        member_ids: list[str],
        role: str = None,
        design_params: dict = None,
        comment: str = None,
    ) -> str:
        """Crée une barre structurale (regroupement de barres FEM continues).

        Regroupe une ou plusieurs barres FEM en un unique élément structurel,
        destiné à être vérifié comme un tout selon un Eurocode (EC2, EC3, EC5, etc.).
        La continuité n'est pas vérifiée ici : utilisez ``detect_continuous_members``
        en amont si vous voulez une vérification automatique.

        Args:
            name (str): Nom unique de la barre structurale (ex: "Solive_S1").
            member_ids (list[str]): Liste ordonnée des identifiants de barres FEM
                constituant la barre structurale (ex: ["M12", "M13", "M14"]).
                Pour une barre simple non continue : ``[member_id]``.
            role (str, optional): Rôle structurel (ex: "Solive", "Panne", "Poteau",
                "Arbalétrier", "Moise", "Poutre", "Voile"). Utilisé par les modules
                de vérification.
            design_params (dict, optional): Paramètres de design spécifiques à l'Eurocode
                utilisé. Le contenu dépend du matériau (EC2=béton, EC3=acier, EC5=bois).
                Exemple EC5: {"classe_bois": "C24", "cs": 1, "Hi": 12, "Hf": 12,
                "effet_systeme": False, "type_element_fleche": "Solives"}.
                Exemple EC3: {"classe_acier": "S355", "classe_section": 1}.
            lo_rel_y (float, optional): Longueur de flambement/déversement autour
                de l'axe y, en mm. Si None, à déterminer automatiquement au moment
                de la vérification.
            lo_rel_z (float, optional): Longueur de flambement/déversement autour
                de l'axe z, en mm.
            comment (str, optional): Commentaire descriptif.

        Returns:
            str: Nom de la barre structurale créée (= ``name``).

        Raises:
            ValueError: Si une barre FEM listée n'existe pas, si ``member_ids`` est
                vide, ou si ``name`` existe déjà.

        Exemple:
            >>> chains = model.detect_continuous_members()
            >>> # supposons chains[0] == ["M1", "M2", "M3"]
            >>> model.group_members("Solive_S1", chains[0],
            ...                     role="Solive",
            ...                     design_params={"classe_bois": "C24",})
        """
        if not member_ids:
            raise ValueError("`member_ids` ne peut pas être vide.")
        if name in self._data["structural_members"]:
            raise ValueError(
                f"La barre structurale '{name}' existe déjà. "
                f"Utilisez `del_structural_member` d'abord pour la remplacer."
            )
        for mid in member_ids:
            if mid not in self._data["members"]:
                raise ValueError(f"La barre FEM '{mid}' n'existe pas dans le modèle.")

        total_length = sum(
            (self._data["members"][mid]["Longueur"] for mid in member_ids),
            start=0 * si.mm,
        )

        design = dict(design_params) if design_params else {}

        # Paramètres de flambement (EC5 §6.3.2) — longueurs et types d'appui par axe
        comp_params = self._determine_compression_params(member_ids)
        design["lo_flamb_y"] = comp_params["lo_flamb_y"]
        design["lo_flamb_z"] = comp_params["lo_flamb_z"]
        design["type_appuis_y"] = comp_params["type_appuis_y"]
        design["type_appuis_z"] = comp_params["type_appuis_z"]

        # Paramètres de déversement en flexion (EC5 §6.3.3) — longueurs et coeflef par axe
        flexion_params = self._determine_flexion_params(member_ids)
        design["lo_rel_y"] = flexion_params["lo_rel_y"]
        design["lo_rel_z"] = flexion_params["lo_rel_z"]
        design["coeflef_y"] = flexion_params["coeflef_y"]
        design["coeflef_z"] = flexion_params["coeflef_z"]

        self._data["structural_members"][name] = {
            "Barres FEM": list(member_ids),
            "Rôle": role,
            "Longueur totale": total_length,
            "Design": design,
            "Commentaire": comment,
        }
        return name

    def _node_is_rotationally_fixed(
        self, node_id: str, member_id: str, end: str, axis: str
    ) -> bool:
        """Indique si une extrémité de barre est effectivement encastrée en rotation.

        Un encastrement en rotation sur ``axis`` requiert **deux conditions**
        simultanées :

        1. Le nœud bloque la rotation via un appui classique (``R{axis}=True``).
        2. La barre elle-même n'a pas de relâchement ``teta_{axis}`` à cette
           extrémité — sinon la rotule de barre annule l'appui en rotation.

        Si le nœud n'a pas d'appui classique (nœud de ferme interne, faîtage…),
        la rotation est libre quelle que soit la continuité de la barre.

        Args:
            node_id: Identifiant du nœud à analyser.
            member_id: Barre FEM dont l'extrémité est au nœud.
            end: ``"start"`` ou ``"end"``.
            axis: ``"y"`` ou ``"z"``.

        Returns:
            bool: ``True`` si le nœud constitue un encastrement en rotation
                pour la barre analysée.
        """
        rot_key = "R" + axis.upper()
        support = self._get_support_conditions_at_node(node_id)
        if not support[rot_key]:
            return False
        return not self._has_rotation_release_at(member_id, end, axis)

    def _type_appuis_between_nodes(
        self,
        node_start: str, node_end: str,
        member_start: str, member_end: str,
        axis: str,
        end_str_start: str, end_str_end: str,
    ) -> str:
        """Retourne le type d'appui en flambement entre deux nœuds consécutifs.

        Args:
            node_start: nœud de début du segment.
            node_end: nœud de fin du segment.
            member_start: barre FEM dont l'extrémité est en ``node_start``.
            member_end: barre FEM dont l'extrémité est en ``node_end``.
            axis: ``"y"`` ou ``"z"``.
            end_str_start: ``"start"`` ou ``"end"`` pour ``member_start`` à ``node_start``.
            end_str_end: ``"start"`` ou ``"end"`` pour ``member_end`` à ``node_end``.

        Returns:
            str: Type d'appui selon ``Compression.COEF_LF``.

        Note:
            L'encastrement en rotation est déterminé par ``_node_is_rotationally_fixed`` :
            il requiert à la fois un appui classique bloquant ``R{axis}`` **et** l'absence
            de relâchement ``teta_{axis}`` sur la barre analysée à cette extrémité.
            Un nœud sans appui classique (nœud de ferme, faîtage…) est toujours une rotule,
            même si la barre n'a pas de relâchement explicite.
        """
        start_fixed = self._node_is_rotationally_fixed(
            node_start, member_start, end_str_start, axis
        )
        end_fixed = self._node_is_rotationally_fixed(
            node_end, member_end, end_str_end, axis
        )
        if start_fixed and end_fixed:
            return "Encastré - Encastré"
        if start_fixed or end_fixed:
            # Un seul côté est encastré (rotation bloquée sans relâchement).
            # L'autre extrémité doit être libre (pas de translation latérale bloquée)
            # pour que ce soit un porte-à-faux, sinon c'est Encastré-Rotule.
            lat_dof = "DZ" if axis == "y" else "DY"
            s_lat = self._get_support_conditions_at_node(node_start)
            e_lat = self._get_support_conditions_at_node(node_end)
            other_lat_free = (start_fixed and not e_lat[lat_dof]) or (end_fixed and not s_lat[lat_dof])
            if other_lat_free:
                return "Encastré 1 côté"
            return "Encastré - Rotule"
        return "Rotule - Rotule"

    def _determine_type_appuis(self, member_ids: list[str], axis: str) -> str:
        """Détermine le type d'appui aux **extrémités** de la barre structurale
        pour le calcul de flambement.

        Utilise ``_ordered_chain_nodes`` pour respecter la topologie réelle
        (barres éventuellement à sens alternés).

        Args:
            member_ids (list[str]): Liste des identifiants de barres FEM.
            axis (str): Axe de flambement ("y" ou "z").

        Returns:
            str: Type d'appui selon Compression.COEF_LF.
        """
        chain_nodes = self._ordered_chain_nodes(member_ids)
        start_node = chain_nodes[0][0]
        end_node = chain_nodes[-1][0]
        end_str_start = self._end_at_node(member_ids[0], start_node)
        end_str_end = self._end_at_node(member_ids[-1], end_node)
        return self._type_appuis_between_nodes(
            start_node, end_node,
            member_ids[0], member_ids[-1],
            axis,
            end_str_start, end_str_end,
        )

    def _determine_compression_params(
        self, member_ids: list[str]
    ) -> dict:
        """Calcule les paramètres de flambement (EC5 §6.3.2) par axe.

        Pour chaque axe (y, z) :

        - Les appuis en **translation** latérale (``DZ`` pour axe y, ``DY`` pour
          axe z) découpent la barre structurale en sous-portées.
        - La longueur de flambement ``lo_flamb`` retenue est la longueur de la
          sous-portée la plus longue (cas le plus défavorable).
        - Le ``type_appuis`` est déterminé aux extrémités de cette sous-portée
          par analyse des conditions de rotation (``RY``/``RZ``) et des
          relâchements éventuels.

        La convention de DDL est :

        =========  =================  ==================
        Axe         Translation lat.   Rotation bloquante
        =========  =================  ==================
        y (plan XZ) DZ                 RY
        z (plan XY) DY                 RZ
        =========  =================  ==================

        Args:
            member_ids (list[str]): Liste ordonnée de barres FEM.

        Returns:
            dict: Clés ``lo_flamb_y``, ``lo_flamb_z``,
                ``type_appuis_y``, ``type_appuis_z`` (flottants mm et str).
        """
        chain_nodes = self._ordered_chain_nodes(member_ids)  # [(node, abs), …]
        node_at_abs: dict[float, str] = {pos: nid for nid, pos in chain_nodes}

        def _worst_span_for_axis(lat_dof: str, axis: str) -> tuple[float, str]:
            rot_dof = "R" + axis.upper()
            spans = self._compute_spans(member_ids, lat_dof, rot_dof)
            worst_lo = 0.0
            worst_type = "Rotule - Rotule"
            for x0, x1, span_len in spans:
                node_s = node_at_abs[x0]
                node_e = node_at_abs[x1]
                # Trouve les barres FEM dont les extrémités tombent sur ces nœuds
                mid_s = next(
                    m for m in member_ids
                    if node_s in self._data["members"][m]["Noeuds"]
                )
                mid_e = next(
                    m for m in reversed(member_ids)
                    if node_e in self._data["members"][m]["Noeuds"]
                )
                end_s = self._end_at_node(mid_s, node_s)
                end_e = self._end_at_node(mid_e, node_e)
                t = self._type_appuis_between_nodes(
                    node_s, node_e, mid_s, mid_e, axis, end_s, end_e
                )
                if span_len >= worst_lo:
                    worst_lo = span_len
                    worst_type = t
            return worst_lo, worst_type

        lo_flamb_y, type_appuis_y = _worst_span_for_axis("DY", "z") #inversé pour convertir le repère locale à l'EUROCODE 5
        lo_flamb_z, type_appuis_z = _worst_span_for_axis("DZ", "y") #inversé pour convertir le repère locale à l'EUROCODE 5
        return {
            "lo_flamb_y": lo_flamb_y,
            "lo_flamb_z": lo_flamb_z,
            "type_appuis_y": type_appuis_y,
            "type_appuis_z": type_appuis_z,
        }

    def _ordered_chain_nodes(
        self, member_ids: list[str]
    ) -> list[tuple[str, float]]:
        """Retourne la liste ordonnée ``[(node_id, abscisse_mm), …]`` pour une chaîne
        de barres FEM consécutives, en respectant la topologie réelle (sens des
        barres potentiellement alternés).

        Args:
            member_ids (list[str]): Barres FEM de la chaîne, dans l'ordre structural.

        Returns:
            list[tuple[str, float]]: Du premier nœud au dernier, avec abscisse
                cumulée en mm.
        """
        if not member_ids:
            return []

        # Pour la première barre on choisit Noeuds[0] comme nœud de départ.
        first = self._data["members"][member_ids[0]]["Noeuds"]
        result: list[tuple[str, float]] = [(first[0], 0.0)]
        prev_node = first[0]
        cumul = 0.0

        for mid in member_ids:
            n1, n2 = self._data["members"][mid]["Noeuds"]
            length_mm = self._data["members"][mid]["Longueur"].value * 1000.0
            # Le nœud de sortie est celui qui n'est pas le nœud d'entrée
            next_node = n2 if n1 == prev_node else n1
            cumul += length_mm
            result.append((next_node, cumul))
            prev_node = next_node

        return result

    def _compute_spans(
        self,
        member_ids: list[str],
        lateral_dof: str,
        rot_dof: str,
    ) -> list[tuple[float, float, float]]:
        """Retourne la liste des sous-portées entre appuis latéraux consécutifs.

        Parcourt la séquence ordonnée de barres FEM et identifie les nœuds
        qui bloquent le déversement. Chaque sous-portée est décrite par ses
        abscisses de début et de fin ainsi que sa longueur.

        Args:
            member_ids (list[str]): Liste ordonnée des identifiants de barres FEM.
            lateral_dof (str): DDL en translation bloquant le déversement
                (ex: ``"DZ"`` pour déversement selon Y).
            rot_dof (str): DDL en rotation bloquant le déversement
                (ex: ``"RY"`` pour déversement selon Y).

        Returns:
            list[tuple[float, float, float]]: Liste de ``(x_debut, x_fin, longueur)``
                en mm, du premier au dernier segment entre appuis latéraux.
                S'il n'y a aucun appui latéral, retourne un seul segment couvrant
                toute la longueur.

        Note:
            Un nœud est considéré appui latéral s'il bloque ``lateral_dof``
            ou ``rot_dof``.
        """
        nodes_with_pos = self._ordered_chain_nodes(member_ids)
        cumul = nodes_with_pos[-1][1] if nodes_with_pos else 0.0

        bracing_positions: list[float] = []
        for node_id, pos in nodes_with_pos:
            support = self._get_support_conditions_at_node(node_id)
            if support[lateral_dof] or support[rot_dof]:
                bracing_positions.append(pos)

        all_positions = sorted(set([0.0] + bracing_positions + [cumul]))
        return [
            (all_positions[i], all_positions[i + 1], all_positions[i + 1] - all_positions[i])
            for i in range(len(all_positions) - 1)
        ]

    def _compute_lo_rel(
        self,
        member_ids: list[str],
        lateral_dof: str,
        rot_dof: str,
    ) -> float:
        """Calcule la longueur de déversement maximale entre appuis latéraux.

        Délègue à ``_compute_spans`` et retourne la longueur du segment le plus long.

        Args:
            member_ids (list[str]): Liste ordonnée des identifiants de barres FEM.
            lateral_dof (str): DDL en translation bloquant le déversement.
            rot_dof (str): DDL en rotation bloquant le déversement.

        Returns:
            float: Longueur de déversement maximale en mm.
        """
        spans = self._compute_spans(member_ids, lateral_dof, rot_dof)
        return max(s[2] for s in spans) if spans else 0.0

    def _determine_flexion_params(self, member_ids: list[str]) -> dict:
        """Détermine les paramètres de déversement en flexion selon l'EC5.

        Analyse les charges appliquées et les conditions d'appui pour déterminer
        les longueurs efficaces de déversement et les coefficients associés.

        Args:
            member_ids (list[str]): Liste des identifiants de barres FEM.

        Returns:
            dict: Paramètres de flexion avec clés :
                - "lo_rel_y" (float): Longueur de déversement axe y en mm
                - "lo_rel_z" (float): Longueur de déversement axe z en mm
                - "coeflef_y" (float): Coefficient de longueur efficace axe y
                - "coeflef_z" (float): Coefficient de longueur efficace axe z
                - "pos_charge" (str): Position de la charge verticale

        Note:
            Les coefficients coeflef dépendent du type de chargement :
            - 1.0 : moment constant
            - 0.9 : charge répartie (défaut)
            - 0.8 : charge concentrée centrale
            - 0.5 : porte-à-faux charge répartie
            - 0.8 : porte-à-faux charge concentrée bout

            La position de la charge est déterminée par la direction des charges
            verticales appliquées.
        """
        # Vérifie si c'est un porte-à-faux (un côté encastré, l'autre libre)
        # _ordered_chain_nodes respecte la topologie réelle des barres enchaînées.
        chain_nodes = self._ordered_chain_nodes(member_ids)
        start_node = chain_nodes[0][0]
        end_node = chain_nodes[-1][0]

        start_support = self._get_support_conditions_at_node(start_node)
        end_support = self._get_support_conditions_at_node(end_node)

        start_blocked = start_support["DX"] or start_support["DY"] or start_support["DZ"]
        end_blocked = end_support["DX"] or end_support["DY"] or end_support["DZ"]

        is_cantilever = (start_blocked and not end_blocked) or (end_blocked and not start_blocked)

        # Abscisses cumulées des barres (calculé une seule fois, capturé par closure)
        cumul_map: dict[str, tuple[float, float]] = {}
        x_cur = 0.0
        for _mid in member_ids:
            _len = self._data["members"][_mid]["Longueur"].value * 1000.0
            cumul_map[_mid] = (x_cur, x_cur + _len)
            x_cur += _len

        def _coeflef_for_span(
            x0: float, x1: float, is_cant: bool
        ) -> float:
            """Retourne le coeflef le plus défavorable pour un segment [x0, x1].

            Pour chaque charge verticale portant sur une barre dont l'abscisse
            locale chevauche le segment, détermine la contribution (distribuée,
            ponctuelle centrale, ponctuelle en bout). Le coeflef le plus élevé
            (le plus défavorable) est retenu ; 1.0 est utilisé si aucune charge
            n'est trouvée (moment constant, hypothèse conservative).
            """
            span_len = x1 - x0
            local_dist = False
            local_center = False
            local_end = False

            for load in self._data["loads"].values():
                if load["N° barre"] not in member_ids:
                    continue
                if load["Axe"] not in ("Fy", "Fz", "FY", "FZ"):
                    continue
                mid = load["N° barre"]
                bar_x0, bar_x1 = cumul_map[mid]
                # La barre doit chevaucher le segment
                if bar_x1 <= x0 or bar_x0 >= x1:
                    continue
                if load["Type de charge"] == "Distribuée":
                    local_dist = True
                elif load["Type de charge"] == "Ponctuelle":
                    pos_rel = load.get("Position", 0.5)
                    # Abscisse absolue de la charge
                    x_load = bar_x0 + pos_rel * (bar_x1 - bar_x0)
                    if x0 <= x_load <= x1:
                        # Position relative dans le segment
                        pos_in_span = (x_load - x0) / span_len if span_len > 0 else 0.5
                        if 0.4 <= pos_in_span <= 0.6:
                            local_center = True
                        else:
                            local_end = True

            if is_cant:
                return 0.8 if local_end else 0.5
            else:
                if local_center:
                    return 0.8
                if local_dist:
                    return 0.9
                return 1.0  # moment constant — hypothèse conservative

        def _worst_case_axis(
            lateral_dof: str, rot_dof: str
        ) -> tuple[float, float]:
            """Retourne ``(lo_rel, coeflef)`` les plus défavorables pour un axe.

            Pour chaque sous-portée entre appuis latéraux, calcule le coeflef
            associé. La paire ``(lo_rel × coeflef)`` maximale détermine la
            longueur efficace de déversement la plus défavorable.
            """
            spans = self._compute_spans(member_ids, lateral_dof, rot_dof)
            worst_lo = 0.0
            worst_coeflef = 0.9  # valeur par défaut
            for x0, x1, span_len in spans:
                c = _coeflef_for_span(x0, x1, is_cantilever)
                if span_len * c >= worst_lo * worst_coeflef:
                    worst_lo = span_len
                    worst_coeflef = c
            return worst_lo, worst_coeflef

        lo_rel_y, coeflef_y = _worst_case_axis("DZ", "RY") #inversé pour convertir le repère locale à l'EUROCODE 5
        lo_rel_z, coeflef_z = _worst_case_axis("DY", "RZ") #inversé pour convertir le repère locale à l'EUROCODE 5

        return {
            "lo_rel_y": lo_rel_y,
            "lo_rel_z": lo_rel_z,
            "coeflef_y": coeflef_y,
            "coeflef_z": coeflef_z,
        }

    def auto_group_continuous_members(
        self,
        angle_tol_deg: float = 1.0,
        role: str = None,
        name_prefix: str = "SM",
        only_continuous: bool = ("False", "True"),
        design_params: dict = None,
    ) -> list[str]:
        """Détecte et regroupe automatiquement toutes les barres continues du modèle.

        Combine ``detect_continuous_members`` et ``group_members`` : chaque chaîne
        détectée devient une barre structurale nommée ``{name_prefix}{i}``.

        Args:
            angle_tol_deg (float): Tolérance angulaire pour la détection.
                Defaults to 1.0.
            role (str, optional): Rôle appliqué à toutes les barres créées.
            name_prefix (str): Préfixe pour les noms auto-générés.
                Defaults to "SM".
            only_continuous (bool): Si True, ignore les chaînes d'une seule barre FEM
                (barres isolées non continues). Defaults to False.
            design_params (dict, optional): Paramètres de design spécifiques à l'Eurocode.
                Transmis à ``group_members`` pour chaque barre créée.

        Returns:
            list[str]: Noms des barres structurales créées, dans l'ordre de détection.

        Exemple:
            >>> names = model.auto_group_continuous_members(
            ...     role="Solive", design_params={"classe_bois": "C24", "cs": 1})
            >>> # ["SM1", "SM2", ...]
        """
        chains = self.detect_continuous_members(angle_tol_deg=angle_tol_deg)
        created: list[str] = []
        idx = 1
        for chain in chains:
            if only_continuous and len(chain) < 2:
                continue
            name = f"{name_prefix}{idx}"
            while name in self._data["structural_members"]:
                idx += 1
                name = f"{name_prefix}{idx}"
            self.group_members(
                name, chain, role=role,
                design_params=design_params,
            )
            created.append(name)
            idx += 1
        return created

    def get_structural_member(self, name: str) -> dict:
        """Retourne les données d'une barre structurale par son nom.

        Args:
            name (str): Nom de la barre structurale (ex: "Solive_S1").

        Returns:
            dict: Dictionnaire contenant :

                - ``"Barres FEM"`` : liste ordonnée des ``member_id`` FEM
                - ``"Rôle"`` : rôle structurel
                - ``"Longueur totale"`` : longueur cumulée avec unité (``si.mm``)
                - ``"Design"`` : paramètres de design spécifiques à l'Eurocode
                - ``"Commentaire"`` : texte libre

        Raises:
            KeyError: Si le nom n'existe pas.
        """
        return self._data["structural_members"][name]

    def get_all_structural_members(self) -> dict:
        """Retourne toutes les barres structurales du modèle.

        Returns:
            dict: Dictionnaire ``{name: données}`` de toutes les barres structurales.
        """
        return self._data["structural_members"]

    def del_structural_member(self, name: str) -> dict:
        """Supprime une barre structurale par son nom.

        Args:
            name (str): Nom de la barre structurale à supprimer.

        Returns:
            dict: Données de la barre structurale supprimée.

        Raises:
            KeyError: Si le nom n'existe pas.
        """
        return self._data["structural_members"].pop(name)

    def _iter_structural_members(self):
        """Itère sur toutes les barres structurales du modèle.

        Yields:
            tuple[str, dict]: Couple ``(name, data)`` pour chaque barre structurale.

        Exemple:
            >>> for name, sm in model.iter_structural_members():
            ...     print(name, sm["Rôle"], sm["Longueur totale"])
        """
        yield from self._data["structural_members"].items()

    def debug_node_compression(self, structural_member_name: str, axis: str = "y") -> None:
        """Affiche les informations de debug pour le calcul du type d'appui en compression.

        Affiche les conditions d'appui et de relâchement aux nœuds d'extrémité de la
        barre structurale, ainsi que le résultat final du type d'appui.

        Args:
            structural_member_name (str): Nom de la barre structurale à analyser.
            axis (str): Axe de flambement ("y" ou "z"). Defaults to "y".
        """
        sm = self._data["structural_members"][structural_member_name]
        member_ids = sm["Barres FEM"]
        chain_nodes = self._ordered_chain_nodes(member_ids)
        start_node = chain_nodes[0][0]
        end_node = chain_nodes[-1][0]
        mid_s = member_ids[0]
        mid_e = member_ids[-1]
        end_s = self._end_at_node(mid_s, start_node)
        end_e = self._end_at_node(mid_e, end_node)
        s_sup = self._get_support_conditions_at_node(start_node)
        e_sup = self._get_support_conditions_at_node(end_node)
        s_rel = self._has_rotation_release_at(mid_s, end_s, axis)
        e_rel = self._has_rotation_release_at(mid_e, end_e, axis)
        rot_key = "R" + axis.upper()
        start_fixed = s_sup[rot_key] and not s_rel
        end_fixed = e_sup[rot_key] and not e_rel
        print(f"=== debug_node_compression : '{structural_member_name}' axe={axis} ===")
        print(f"  Barres FEM      : {member_ids}")
        print(f"  Nœud début      : {start_node}  (barre {mid_s}, extrémité '{end_s}')")
        print(f"    Appui classique : {s_sup}")
        print(f"    Relâchement teta_{axis} : {s_rel}")
        print(f"    → start_fixed = {start_fixed}")
        print(f"  Nœud fin        : {end_node}  (barre {mid_e}, extrémité '{end_e}')")
        print(f"    Appui classique : {e_sup}")
        print(f"    Relâchement teta_{axis} : {e_rel}")
        print(f"    → end_fixed = {end_fixed}")
        result = self._determine_type_appuis(member_ids, axis)
        print(f"  → type_appuis   : '{result}'")

    def get_type_appuis_for_compression(
        self,
        structural_member_name: str,
        axis: str = "y",
    ) -> str:
        """Détermine le type d'appui pour le calcul de flambement en compression.

        Analyse les conditions d'appui aux extrémités d'une barre structurale
        (membre continu ou simple) pour déterminer le coefficient β de longueur
        efficace selon l'EC5.

        Args:
            structural_member_name (str): Nom de la barre structurale à analyser
                (ex: "Poteau_P1", "SM1").
            axis (str): Axe de flambement à considérer ("y" ou "z").
                Defaults to "y".

        Returns:
            str: Type d'appui selon Compression.COEF_LF :
                - "Encastré 1 côté" : β = 2.0 (console)
                - "Rotule - Rotule" : β = 1.0 (articulé-articulé)
                - "Encastré - Rotule" : β = 0.7
                - "Encastré - Encastré" : β = 0.5
                - "Encastré - Rouleau" : β = 1.0 (encastré-glissière)

        Raises:
            KeyError: Si la barre structurale n'existe pas.

        Note:
            La méthode analyse :
            1. Les appuis classiques (DX, DY, DZ, RX, RY, RZ) aux nœuds d'extrémité
            2. Les relâchements (releases) en rotation aux extrémités des barres FEM
            3. La continuité du membre (barre continue ou simple)

            Pour un axe donné (y ou z), une rotation bloquée sur cet axe
            correspond à un encastrement, une rotation libre à une rotule.
        """
        sm = self._data["structural_members"][structural_member_name]
        member_ids = sm["Barres FEM"]

        return self._determine_type_appuis(member_ids, axis)

    def _get_support_conditions_at_node(self, node_id: str) -> dict:
        """Retourne les conditions d'appui au nœud spécifié.

        Args:
            node_id (str): Identifiant du nœud.

        Returns:
            dict: Conditions d'appui avec clés DX, DY, DZ, RX, RY, RZ.
                True = bloqué, False = libre.
                Si aucun appui, retourne toutes les directions libres.
        """
        # Cherche un appui classique sur ce nœud
        for support_id, support in self._data["supports"]["classic"].items():
            if support["Noeud"] == node_id:
                return {
                    "DX": support["DX"],
                    "DY": support["DY"],
                    "DZ": support["DZ"],
                    "RX": support["RX"],
                    "RY": support["RY"],
                    "RZ": support["RZ"],
                }
        # Pas d'appui = libre en translation et rotation
        return {"DX": False, "DY": False, "DZ": False, "RX": False, "RY": False, "RZ": False}

    def _has_rotation_release_at(self, member_id: str, end: str, axis: str) -> bool:
        """Indique si la barre a un relâchement en rotation à l'extrémité.

        Args:
            member_id (str): Identifiant de la barre FEM.
            end (str): Extrémité ("start" ou "end").
            axis (str): Axe de rotation ("y" ou "z").

        Returns:
            bool: True si relâchement en rotation sur l'axe, False sinon.
        """
        rel = self._data["members"][member_id]["Relaxation"].get(end)
        if not rel:
            return False
        rot_key = "teta_" + axis
        return bool(rel.get(rot_key, False))

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
