# coding in UTF-8
# by Anthony PARISOT
import os, sys
import math as mt
import numpy as np
from matplotlib import pyplot as plt

import forallpeople as si

from ourocode.eurocode.core.projet import Projet


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

