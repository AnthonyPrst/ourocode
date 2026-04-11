# coding in UTF-8 
# by Anthony PARISOT
import os
import json
import math as mt
import csv
import importlib.resources as pkg_resources
from collections.abc import Mapping, Iterable
from PIL import Image
import pandas as pd
import pickle
import inspect
from IPython.display import display, Latex

import forallpeople as si
si.environment("structural")


def get_package_path(package):
    # Obtenir le chemin du fichier principal du module
    package_path = os.path.dirname(package.__file__)
    return package_path


class Objet(object):
    """Classe permettant la sauvegarde et l'ouverture d'objets dans des fichiers .oco (format pickle).

    Cette classe de base fournit les fonctionnalités essentielles pour :
    - La sérialisation/désérialisation des objets avec gestion des unités physiques
    - L'accès aux données normatives (CSV, JSON)
    - Les opérations mathématiques et de synthèse de résultats
    - La conversion des unités SI

    Toutes les classes du package ourocode héritent de cette classe.
    """
    JUPYTER_DISPLAY = False
    OPERATOR = ("+", "-", "x", "/")
    _csv_cache: dict[tuple, pd.DataFrame] = {}
    _json_cache: dict[str, dict] = {}
    try:
        import ourocode
        PATH_CATALOG = os.path.join(get_package_path(ourocode))
    except ImportError:
        PATH_CATALOG = os.path.join(os.getcwd(), "ourocode")

    def _data_from_csv(self, data_file: str, index_col=0):
        """Charge et retourne les données normatives depuis un fichier CSV.

        Cette méthode permet d'accéder aux tables de données Eurocode stockées
        dans le répertoire data/ du package (caractéristiques mécaniques, kmod, etc.)
        Les résultats sont mis en cache au niveau de la classe pour éviter
        de relire le disque à chaque appel.

        Args:
            data_file (str): Nom du fichier CSV à charger (ex: "kmod.csv", "gammaM.csv")
            index_col (int, optional): Colonne à utiliser comme index. Defaults to 0.

        Returns:
            pd.DataFrame: DataFrame pandas contenant les données normatives.
        """
        cache_key = (data_file, index_col)
        if cache_key not in Objet._csv_cache:
            repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
            Objet._csv_cache[cache_key] = pd.read_csv(repertory, sep=';', header=0, index_col=index_col)
        return Objet._csv_cache[cache_key].copy()
    
    def _data_from_json(self, data_file: str):
        """Charge et retourne les données depuis un fichier JSON sous forme de DataFrame.

        Args:
            data_file (str): Nom du fichier JSON à charger.

        Returns:
            pd.DataFrame: DataFrame pandas contenant les données.
        """
        repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
        data_json = pd.read_json(repertory)
        return data_json

    def _load_json(self, data_file: str):
        """Charge et retourne les données depuis un fichier JSON sous forme de dictionnaire.

        Les résultats sont mis en cache au niveau de la classe.

        Args:
            data_file (str): Nom du fichier JSON à charger.

        Returns:
            dict: Dictionnaire contenant les données du fichier.
        """
        if data_file not in Objet._json_cache:
            repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
            with open(repertory, "r", encoding="utf-8") as json_file:
                Objet._json_cache[data_file] = json.load(json_file)
        return Objet._json_cache[data_file]

    def _assign_handcalcs_value(self, handcalc_value: tuple, args: list[str]):
        """Assigne les valeurs calculées par handcalcs aux attributs de l'objet.

        Cette méthode extrait la valeur numérique du résultat handcalcs (qui peut
        contenir du LaTeX) et l'assigne aux attributs spécifiés.

        Args:
            handcalc_value (tuple): Tuple retourné par le décorateur @handcalc
                contenant (latex_string, numeric_result) ou (result,) selon le mode
            args (list[str]): Liste des noms d'attributs à assigner, dans le même
                ordre que les valeurs retournées par le calcul.
        """
        if isinstance(handcalc_value, tuple):
            if self.JUPYTER_DISPLAY:
                for i, value in enumerate(handcalc_value):
                    print(args[i], value)
                    setattr(self, args[i], value)
            else:
                for i, value in enumerate(handcalc_value[1]):
                    setattr(self, args[i], value)
    
    @property
    def objet(self):
        """Retourne l'objet lui même.
        """
        return self

    def _physical_to_dict(self, obj):
        """Convertit un objet Physical (forallpeople) en dictionnaire sérialisable.

        Cette méthode permet de sérialiser les objets contenant des unités physiques
        (m, mm, MPa, kN, etc.) pour la sauvegarde JSON ou pickle.
        La conversion préserve la valeur numérique et l'unité sous forme de chaîne.

        Args:
            obj: Objet à convertir (dict, list, tuple, ou Physical)

        Returns:
            Objet converti avec les valeurs Physical transformées en dictionnaires
            contenant '_physical_value' et '_physical_unit'.
        """
        if isinstance(obj, dict):
            return {k: self._physical_to_dict(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._physical_to_dict(x) for x in obj]
        elif isinstance(obj, si.Physical):
            split_vals = obj.split(base_value=True)
            return {
                '_physical_value': float(split_vals[0]),
                '_physical_unit': str(split_vals[1])
            }
        return obj

    @staticmethod
    def _resolve_unit_expr(unit_expr: str):
        """Résout une expression d'unité de manière sécurisée sans eval().

        Gère les unités simples (m, mm, kN), les puissances (mm**2, m**4)
        et les produits (kN*m, N*mm**2) par parsing explicite.

        Args:
            unit_expr (str): Expression d'unité normalisée.

        Returns:
            Physical: Objet unité forallpeople correspondant.

        Raises:
            ValueError: Si l'unité est inconnue.
        """
        si_ns = vars(si)
        # Produit : kN*m, N*mm**2, etc.
        if "*" in unit_expr and "**" not in unit_expr.split("*")[0]:
            parts = unit_expr.split("*", 1)
            return Objet._resolve_unit_expr(parts[0]) * Objet._resolve_unit_expr(parts[1])
        # Puissance : mm**2, m**4, etc.
        if "**" in unit_expr:
            base, exp_str = unit_expr.split("**", 1)
            # Gérer le cas kN*m**2 déjà splité à gauche
            if not exp_str.lstrip("-").isdigit():
                raise ValueError(f"Exposant non numérique: {exp_str}")
            base_unit = si_ns.get(base)
            if base_unit is None:
                raise ValueError(f"Unité de base inconnue: {base}")
            return base_unit ** int(exp_str)
        # Unité simple : m, mm, kN, MPa, etc.
        unit_obj = si_ns.get(unit_expr)
        if unit_obj is None:
            raise ValueError(f"Unité inconnue: {unit_expr}")
        return unit_obj

    def _dict_to_physical(self, data):
        """Reconstruit un objet Physical (forallpeople) à partir d'un dictionnaire.

        Cette méthode est l'inverse de _physical_to_dict. Elle restaure les
        objets Physical à partir de leur représentation sérialisée.
        Gère les conversions d'unités avec différentes notations (², ³, ·, etc.)

        Args:
            data: Données sérialisées (dict, list, ou valeur simple)

        Returns:
            Objet avec les valeurs physiques restaurées.
        """
        if isinstance(data, dict):
            if '_physical_value' in data:
                # Reconstruire l'objet Physical
                value = data['_physical_value']
                unit_str = str(data['_physical_unit'].split()[-1])
                # Normaliser quelques notations possibles
                unit_expr = (
                    unit_str
                    .replace('^', '**')
                    .replace('²', '**2')
                    .replace('³', '**3')
                    .replace('⁴', '**4')
                    .replace('·', '*')
                )
                if "/" in unit_expr:
                    unit_expr = unit_expr.replace("/", "_")
                try:
                    unit_obj = self._resolve_unit_expr(unit_expr)
                except ValueError as e:
                    raise ValueError(f"Unité inconnue ou non prise en charge: {unit_str}") from e
                return value * unit_obj
            return {k: self._dict_to_physical(v) for k, v in data.items()}
        elif isinstance(data, (list, tuple)):
            return [self._dict_to_physical(x) for x in data]
        return data

    def get_value(self, value: dict|list|str, index: int=None, key: str=None, get_keys: bool=("False", "True"),):
        """Extrait et retourne une valeur depuis une structure de données complexe.

        Cette méthode utilitaire permet de naviguer dans les dictionnaires,
        listes, ou chaînes JSON pour extraire des valeurs spécifiques.

        Args:
            value (dict|list|str): La structure de données source.
            index (int, optional): Index à extraire dans une liste Python.
                Note : sous Python, le premier élément est à l'index 0.
            key (str, optional): Clé à extraire d'un dictionnaire Python,
                ou clé JSON à extraire d'une chaîne.
            get_keys (bool, optional): Si True, retourne la liste des clés
                d'un dictionnaire au lieu d'une valeur.

        Returns:
            La valeur extraite selon les critères spécifiés.

        Exemple:
            >>> obj.get_value({"a": 1, "b": 2}, key="a")
            1
            >>> obj.get_value([10, 20, 30], index=1)
            20
            >>> obj.get_value({"x": 1, "y": 2}, get_keys=True)
            ["x", "y"]
        """
        if index and isinstance(value, list):
            value = value[index]
        elif index and isinstance(value, str):
            value = list(value)[index]
        elif get_keys and isinstance(value, dict):
            value = list(value.keys())
        elif key and isinstance(value, dict):
            value = value[key]
        elif key and isinstance(value, str):
            value = json.loads(f"{value}".replace("'", "\"")).get(key)
        return value
        
    def operation_between_values(self, value1: float, value2: float, operator: str=OPERATOR):
        """Effectue une opération arithmétique entre deux valeurs.

        Opérateurs supportés : "+" (addition), "-" (soustraction),
        "x" (multiplication), "/" (division).

        Args:
            value1 (float): Première opérande.
            value2 (float): Deuxième opérande.
            operator (str): Opérateur à appliquer. Defaults to "+".

        Returns:
            float: Résultat de l'opération.

        Raises:
            ValueError: Si l'opérateur n'est pas reconnu.
        """
        if operator not in self.OPERATOR:
            raise ValueError(f"Invalid operator: {operator}")
        match operator:
            case "+":
                result = value1 + value2
            case "-":
                result = value1 - value2
            case "x":
                result = value1 * value2
            case "/":
                result = value1 / value2
        return result
    
    def abs_value(self, value: float):
        """Retourne la valeur absolue (module) d'un nombre.

        Args:
            value (float): Valeur dont on veut la valeur absolue.

        Returns:
            float: Valeur absolue (toujours positive ou nulle).
        """
        return abs(float(value))

    
    def max(self, value1: float, value2: float):
        """Retourne la valeur maximale entre deux nombres.

        Args:
            value1 (float): Première valeur à comparer.
            value2 (float): Deuxième valeur à comparer.

        Returns:
            float: La plus grande des deux valeurs.
        """
        return max(float(value1), float(value2))

    def min(self, value1: float, value2: float):
        """Retourne la valeur minimale entre deux nombres.

        Args:
            value1 (float): Première valeur à comparer.
            value2 (float): Deuxième valeur à comparer.

        Returns:
            float: La plus petite des deux valeurs.
        """
        return min(float(value1), float(value2))
    
    def _extract_numbers(self, value, absolute: bool=True):
        """Extrait récursivement tous les nombres d'une structure de données imbriquée.

        Parcourt les dictionnaires, listes et tuples pour extraire toutes les
        valeurs numériques, y compris celles encapsulées dans des objets Physical.

        Args:
            value: Structure de données à parcourir (dict, list, tuple, ou valeur simple)
            absolute (bool, optional): Si True, retourne les valeurs absolues.
                Defaults to True.

        Returns:
            list: Liste de toutes les valeurs numériques extraites.
        """
        numbers = []
        if isinstance(value, (int, float, Physical)):
            if isinstance(value, Physical):
                value = value.split(base_value=True)[0]
            if absolute:
                value = abs(value)
            numbers.append(float(value))
        elif isinstance(value, (list, tuple)):
            for item in value:
                numbers.extend(self._extract_numbers(item, absolute))
        elif isinstance(value, dict):
            for item in value.values():
                numbers.extend(self._extract_numbers(item, absolute))
        return numbers

    def max_list(self, iterable: dict|list|tuple, absolute: bool=("False", "True")):
        """Retourne la valeur maximale d'une liste ou d'un dictionnaire.

        Cette méthode parcourt récursivement la structure pour extraire
        toutes les valeurs numériques et retourne la maximum.

        Args:
            iterable (dict|list|tuple): Structure de données contenant des valeurs numériques.
            absolute (bool, optional): Si True, compare les valeurs absolues.
                Defaults to False.

        Returns:
            float: La valeur maximale trouvée.

        Exemple:
            >>> obj.max_list([1, -5, 3], absolute=True)
            5
            >>> obj.max_list({"a": 10, "b": [20, 5]})
            20
        """
        result = self.get_value(iterable, get_keys=False)
        result = self._extract_numbers(result, absolute)  
        return max(result)
    

    def min_list(self, iterable: dict|list|tuple|str, absolute: bool=("False", "True")):
        """Retourne la valeur minimale d'une liste ou d'un dictionnaire.

        Cette méthode parcourt récursivement la structure pour extraire
        toutes les valeurs numériques et retourne le minimum.

        Args:
            iterable (dict|list|tuple|str): Structure de données contenant des valeurs numériques.
            absolute (bool, optional): Si True, compare les valeurs absolues.
                Defaults to False.

        Returns:
            float: La valeur minimale trouvée.

        Exemple:
            >>> obj.min_list([1, -5, 3], absolute=True)
            1
            >>> obj.min_list({"a": 10, "b": [20, 5]})
            5
        """
        result = self.get_value(iterable, get_keys=False)
        result = self._extract_numbers(result, absolute)  
        return min(result)

    def get_trigonometric_value(self, value: float, operator: str=("COS", "SIN", "TAN", "ACOS", "ASIN", "ATAN")):
        """Calcule une fonction trigonométrique avec un angle exprimé en degrés.

        Convertit automatiquement l'angle de degrés en radians avant le calcul.

        Args:
            value (float): Angle en degrés.
            operator (str): Fonction trigonométrique à appliquer.
                Choix possibles : "COS", "SIN", "TAN", "ACOS", "ASIN", "ATAN"

        Returns:
            float: Résultat de la fonction trigonométrique.

        Raises:
            ValueError: Si la fonction trigonométrique n'est pas reconnue.

        Exemple:
            >>> obj.get_trigonometric_value(90, "SIN")
            1.0
            >>> obj.get_trigonometric_value(0, "COS")
            1.0
        """
        if operator not in ("COS", "SIN", "TAN", "ACOS", "ASIN", "ATAN"):
            raise ValueError(f"La fonction trigonométrique {operator} n'est pas reconnue.")
        match operator:
            case "COS":
                result = mt.cos(mt.radians(float(value)))
            case "SIN":
                result = mt.sin(mt.radians(float(value)))
            case "TAN":
                result = mt.tan(mt.radians(float(value)))
            case "ACOS":
                result = mt.acos(mt.radians(float(value)))
            case "ASIN":
                result = mt.asin(mt.radians(float(value)))
            case "ATAN":
                result = mt.atan(mt.radians(float(value)))
        return result
    
    def _add_synthese_taux_travail(
        self,
        lignes: Iterable,
        col_normale: str="Situation normale",
        col_incendie: str="Situation d'incendie",
        arrondi: int|None=3
    ) -> pd.DataFrame:
        """Construit un tableau pandas de synthèse des taux de travail.

        Cette méthode agrège les résultats de vérification (taux de travail)
        pour la situation normale et la situation d'incendie dans un DataFrame
        partagé entre toutes les méthodes de vérification.

        Args:
            lignes (Iterable): Éléments de synthèse sous forme de :
                - dict avec clés : "effort"/"Effort", "SN"/"situation_normale",
                  "SI"/"situation_incendie"
                - tuple/list : (effort, taux_situation_normale, taux_situation_incendie)
            col_normale (str, optional): Nom de la colonne pour la situation normale.
                Defaults to "Situation normale".
            col_incendie (str, optional): Nom de la colonne pour la situation d'incendie.
                Defaults to "Situation d'incendie".
            arrondi (int|None, optional): Nombre de décimales à conserver.
                None pour ne pas arrondir. Defaults to 3.

        Returns:
            pd.DataFrame: Tableau de synthèse cumulatif avec les taux de travail.
        """
        colonnes = ["Effort agissant dimensionnant", col_normale, col_incendie]
        # DataFrame commun stocké sur l'instance
        dataframe = getattr(self, "_synthese_taux_df", None)
        lignes_preparees = []

        def _arrondi(val):
            if val is None or (isinstance(val, float) and mt.isnan(val)):
                return None
            try:
                val = float(val)
                return round(val, arrondi) if arrondi is not None else val
            except Exception:
                return val

        for ligne in lignes:
            if isinstance(ligne, Mapping):
                effort = ligne.get("effort") or ligne.get("Effort")
                sn = (
                    ligne.get("SN")
                    or ligne.get("situation_normale")
                    or ligne.get(col_normale)
                    or ligne.get("normale")
                )
                si_val = (
                    ligne.get("SI")
                    or ligne.get("situation_incendie")
                    or ligne.get(col_incendie)
                    or ligne.get("incendie")
                )
            elif isinstance(ligne, (list, tuple)):
                effort, sn, si_val = (list(ligne) + [None] * 3)[:3]
            else:
                raise TypeError("Chaque ligne doit être un dict, une liste ou un tuple.")

            lignes_preparees.append({
                "Effort agissant dimensionnant": effort,
                col_normale: _arrondi(sn),
                col_incendie: _arrondi(si_val),
            })

        df_nouveau = pd.DataFrame(lignes_preparees, columns=colonnes)

        if dataframe is not None:
            dataframe = dataframe.copy()
            df_final = pd.concat([dataframe, df_nouveau], join="inner", ignore_index=True)
        else:
            df_final = df_nouveau

        # On mémorise le tableau commun sur l'instance pour les appels suivants
        self._synthese_taux_df = df_final
        return df_final
    
    def synthese_taux_travail(self):
        """Construit un tableau pandas de synthèse des taux de travail.

        Cette méthode agrège les résultats de vérification (taux de travail)
        pour la situation normale et la situation d'incendie dans un DataFrame
        partagé entre toutes les méthodes de vérification.

        Returns:
            pd.DataFrame: Tableau de synthèse final avec la ligne des maximums.

        Note:
            Le DataFrame est partagé entre toutes les méthodes de vérification
            d'une même instance via l'attribut _synthese_taux_df.
        """
        df = getattr(self, "_synthese_taux_df", None)
        if df is None or df.empty:
            return df
        # On suppose que les 2 dernières colonnes sont les taux (situation normale / incendie)
        col_normale, col_incendie = df.columns[-2], df.columns[-1]
        max_normale = pd.to_numeric(df[col_normale], errors="coerce").max(skipna=True)
        max_incendie = pd.to_numeric(df[col_incendie], errors="coerce").max(skipna=True)
        ligne_max = {
            df.columns[0]: "Max",
            col_normale: max_normale,
            col_incendie: max_incendie,
        }
        self._synthese_taux_df = pd.concat([df, pd.DataFrame([ligne_max], columns=df.columns)], ignore_index=True)
        return self._synthese_taux_df

    
    def save_data(self, data: dict, type_data: str=("JSON", "CSV"), path: str=None):
        """Sauvegarde les données dans un fichier JSON ou CSV via boîte de dialogue.

        Convertit automatiquement les unités physiques (Physical) en valeurs
        sérialisables avant la sauvegarde.

        Args:
            data (dict): Données à sauvegarder sous forme de dictionnaire.
                Les valeurs Physical sont automatiquement converties.
            type_data (str): Format de sauvegarde : "JSON" ou "CSV".
            path (str, optional): Chemin du fichier à créer. Si None, une boîte
                de dialogue Qt s'ouvre pour choisir l'emplacement.

        Note:
            Pour JSON : sauvegarde indentée avec encodage UTF-8.
            Pour CSV : format simple avec une ligne d'en-tête.
        """
        from PySide6.QtWidgets import QFileDialog
        data = self._physical_to_dict(data)
        if type_data == "JSON":
            save_file_path = path if path else QFileDialog.getSaveFileName(
                filter="JSON (*.json)",
                selectedFilter=".json",
            )[0]
            with open(save_file_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=4)
        elif type_data == "CSV":
            save_file_path = path if path else QFileDialog.getSaveFileName(
                filter="CSV (*.csv)",
                selectedFilter=".csv",
            )[0]
            with open(save_file_path, "w", newline="") as f:
                w = csv.DictWriter(f, data.keys())
                w.writeheader()
                w.writerow(data)

    def load_data(self, type_data: str=("JSON", "CSV"), path: str=None):
        """Charge les données depuis un fichier JSON ou CSV via boîte de dialogue.

        Reconstruit automatiquement les unités physiques (Physical) à partir
        des données sérialisées.

        Args:
            type_data (str): Format du fichier : "JSON" ou "CSV".
            path (str, optional): Chemin du fichier à charger. Si None, une boîte
                de dialogue Qt s'ouvre pour sélectionner le fichier.

        Returns:
            dict: Dictionnaire contenant les données chargées avec les unités
                physiques restaurées.
        """
        from PySide6.QtWidgets import QFileDialog
        if type_data == "JSON":
            file_path = path if path else QFileDialog.getOpenFileName(
                filter="JSON (*.json)",
                selectedFilter=".json",
            )[0]
            with open(file_path, "r", encoding="utf-8") as f:
                return self._dict_to_physical(json.load(f))
        elif type_data == "CSV":
            file_path = path if path else QFileDialog.getOpenFileName(
                filter="CSV (*.csv)",
                selectedFilter=".csv",
            )[0]
            with open(file_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                return self._dict_to_physical({row['key']: row['value'] for row in reader})

    @classmethod
    def _convert_unit_physical(cls, value: int|float, si_unit: si.Physical, unit_to_convert: si.Physical):
            """Convertit une valeur d'une unité SI vers une autre unité.

            Cette méthode gère les conversions entre les unités courantes en
            calcul de structure : m/mm, N/kN/daN, Pa/MPa, et les unités
            combinées (moments, pressions, etc.).

            Args:
                value (int|float): Valeur numérique à convertir.
                si_unit (si.Physical): Unité source (ex: si.m, si.N).
                unit_to_convert (si.Physical): Unité cible (ex: si.mm, si.kN).

            Returns:
                float: Valeur convertie dans l'unité cible.

            Exemple:
                >>> Objet._convert_unit_physical(1.5, si.m, si.mm)
                1500.0
                >>> Objet._convert_unit_physical(5000, si.N, si.kN)
                5.0
            """
            si_unit, unit_to_convert = str(si_unit), str(unit_to_convert)
            if si_unit != unit_to_convert:
                if si_unit == str(si.m):
                    if unit_to_convert == str(si.mm):
                        return value * 10**3
                    elif unit_to_convert == str(si.cm):
                        return value * 10**2
                    elif unit_to_convert == str(si.km):
                        return value * 10**-3
                elif si_unit == str(si.m**2):
                    if unit_to_convert == str(si.mm**2):
                        return value * 10**6
                    elif unit_to_convert == str(si.cm**2):
                        return value * 10**4
                elif si_unit == str(si.m**3):
                    if unit_to_convert == str(si.mm**3):
                        return value * 10**9
                    elif unit_to_convert == str(si.cm**3):
                        return value * 10**6
                elif si_unit == str(si.m**4):
                    if unit_to_convert == str(si.mm**4):
                        return value * 10**12
                    elif unit_to_convert == str(si.cm**4):
                        return value * 10**8
                elif si_unit == str(si.N):
                    if unit_to_convert == str(si.kN):
                        return value * 10**-3
                    elif unit_to_convert == str(si.daN):
                        return value * 10**-1
                elif si_unit == str(si.N*si.m):
                    if unit_to_convert == str(si.kN*si.m):
                        return value * 10**-3
                    elif unit_to_convert == str(si.daN*si.m):
                        return value * 10**-1
                    elif unit_to_convert == str(si.N*si.mm):
                        return value * 10**3
                elif si_unit == str(si.N/si.m):
                    if unit_to_convert == str(si.kN/si.m):
                        return value * 10**-3
                    elif unit_to_convert == str(si.daN/si.m):
                        return value * 10**-1
                    elif unit_to_convert == str(si.N/si.mm):
                        return value * 10**-3
                elif si_unit == str(si.N/si.m**2):
                    if unit_to_convert == str(si.kN/si.m**2):
                        return value * 10**-3
                    elif unit_to_convert == str(si.daN/si.m**2):
                        return value * 10**-1
                    elif unit_to_convert == str(si.N/si.mm**2):
                        return value * 10**-6
                elif si_unit == str(si.Pa):
                    if unit_to_convert == str(si.kPa):
                        return value * 10**-3
                    elif unit_to_convert == str(si.MPa):
                        return value * 10**-6
            return value
    
    @classmethod           
    def _reset_physical_dictionnary(cls, objet: object, dictionnary: dict) -> dict:
        """Réinitialise les valeurs physiques d'un dictionnaire d'arguments.

        Lors du chaînage de classes via _from_parent_class, les unités physiques
        peuvent être multipliées par elles-mêmes. Cette méthode réinitialise
        les valeurs en extrayant la partie numérique et en réappliquant
        l'unité attendue par l'annotation de type de la classe parent.

        Args:
            objet (object): Objet dont on veut réinitialiser les valeurs physiques.
            dictionnary (dict): Dictionnaire d'arguments contenant potentiellement
                des valeurs Physical malformées.

        Returns:
            dict: Dictionnaire avec les valeurs physiques corrigées.
        """
        dict_physical = {}
        # Si un argument utilise forallpeople on récupère que la valeur pour ne pas multiplier l'unité par elle même
        for key, val in dictionnary.items():
            if isinstance(val, si.Physical):
                physical = val.split(base_value=True)
                # On test si l'objet est une classe ou une instance de classe
                if isinstance(objet, type):
                    mro = objet.mro()
                else:
                    mro = type(objet).mro()
                for objt in mro:
                    spec = inspect.getfullargspec(objt.__init__).annotations
                    if spec.get(key):
                        unit = spec[key]
                        value = cls._convert_unit_physical(physical[0], physical[1], unit)
                        dict_physical[key] = value
                        break
        return dict_physical
    
    @classmethod           
    def _reset_physical_object(cls, objet: object):
        """Réinitialise les valeurs physiques d'un objet en utilisant son __dict__.

        Wrapper de _reset_physical_dictionnary qui extrait automatiquement
        le dictionnaire d'attributs de l'objet.

        Args:
            objet (object): Objet à réinitialiser.

        Returns:
            dict: Dictionnaire des attributs avec valeurs physiques corrigées.
        """
        dictionnary = objet.__dict__
        return cls._reset_physical_dictionnary(objet, dictionnary)
    
    

    @classmethod
    def _from_dict(cls, dictionary:dict):
        """Instancie une classe à partir d'un dictionnaire d'arguments.

        Cette méthode de classe permet de créer une nouvelle instance
        à partir d'un dictionnaire de paramètres.

        Args:
            dictionary (dict): Dictionnaire contenant les arguments nommés
                pour le constructeur de la classe.

        Returns:
            instance: Nouvelle instance de la classe initialisée avec les
                arguments du dictionnaire.
        """
        return cls(**dictionary)
    
    @classmethod
    def _from_parent_class(cls, objet: list|object, **kwargs):
        """Instancie une classe en héritant des attributs d'objets existants.

        Cette méthode de classe est le mécanisme central d'enchaînement des
        vérifications Eurocode. Elle permet de créer une nouvelle instance
        en copiant les attributs d'objets parent(s), avec possibilité de
        surcharge via les kwargs.

        Pattern : héritage par composition pour les vérifications en cascade.

        Args:
            objet (object|list): Objet ou liste d'objets sources dont on veut
                hériter les attributs. Les derniers objets écrasent les précédents.
            **kwargs: Arguments nommés qui écrasent les attributs hérités.

        Returns:
            instance: Nouvelle instance de la classe avec les attributs fusionnés.

        Exemple:
            >>> barre = Barre(b=100, h=200, classe="C24")
            >>> flexion = Flexion._from_parent_class(barre, lo_rel_y=3000)
            # Flexion hérite de b, h, classe depuis Barre
        """
        dict_objet = {}
        
        # Récupération des attributs des objets sources
        if isinstance(objet, list):
            # Pour une liste d'objets, chaque objet écrase les précédents
            for obj in objet:
                if hasattr(obj, "__dict__"):
                    dict_objet.update(obj.__dict__)
                    dict_objet.update(cls._reset_physical_object(obj))
        elif hasattr(objet, "__dict__"):
            # Pour un seul objet
            dict_objet.update(objet.__dict__)
            dict_objet.update(cls._reset_physical_object(objet))
        
        # On met à jour avec les kwargs qui écrase tout
        dict_objet.update(kwargs)
        return cls(**dict_objet)

    
    def _save_muliple_objects(self, object: list):
        from PySide6.QtWidgets import QFileDialog
        save_file_path = QFileDialog.getSaveFileName(
                filter="Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                selectedFilter=".oco",
            )[0]
        # with filedialog.asksaveasfile('wb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
        with open(save_file_path, "wb") as f:
            for ligne in object:
                pickle.dump(ligne, f)
    
    
    def save_object(self):
        from PySide6.QtWidgets import QFileDialog
        save_file_path = QFileDialog.getSaveFileName(
                filter="Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                selectedFilter=".oco",
            )[0]
        with open(save_file_path, "wb") as f:
            pickle.dump(self, f)


    
    def _show_element(self, picture: str):
        """Affiche l'image des caractéristiques d'une entaille au cisaillement
        """
        file = os.path.join(self.PATH_CATALOG, "data", "screenshot", picture)
        image = Image.open(file)
        image.show()
            
            
    @classmethod
    def _open_multiple_objects(cls):
        from PySide6.QtWidgets import QFileDialog
        data = []
        # with filedialog.askopenfile('rb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
        file_path = QFileDialog.getOpenFileName(
                    filter="Ourea catalog object (*.oco);;'Text Document' (*.txt)", selectedFilter=".oco"
                )[0]
        with open(file_path, "rb") as f:
            while True:
                try:
                    data.append(pickle.load(f))
                except EOFError:
                    break
            return data
    
    @classmethod
    def _open_object(cls, path: str=None):
        from PySide6.QtWidgets import QFileDialog
        if not path:
            file_path = QFileDialog.getOpenFileName(
                    filter="Ourea catalog object (*.oco);;'Text Document' (*.txt)", selectedFilter=".oco"
                )[0]
            with open(file_path, "rb") as f:
                return pickle.load(f)
        else:
            with open(path, mode="rb") as f:
                return pickle.load(f)
