# coding in UTF-8
# by Anthony PARISOT
import math as mt
import inspect

import forallpeople as si
si.environment("structural")
from forallpeople import Physical


class MathUtilsMixin:
    """Mixin fournissant les utilitaires mathématiques et de conversion d'unités.

    Méthodes :
    - abs_value, max, min : wrappers arithmétiques
    - _extract_numbers, max_list, min_list : extraction récursive de valeurs numériques
    - get_trigonometric_value : fonctions trigonométriques en degrés
    - _convert_unit_physical : conversion entre unités SI
    - _reset_physical_dictionnary / _reset_physical_object : réinitialisation des Physical
    """

    OPERATOR = ("+", "-", "x", "/")
    # Registre des unites physiques supportees par set_value.
    # L'ordre du tuple SET_VALUE_UNITS sert de valeurs du combobox dans l'UI.
    # Les cles sont en ASCII (pas d'exposants) pour la saisie utilisateur.
    _PHYSICAL_UNITS = {
        # Longueurs
        "m": si.m, "mm": si.mm,
        # Surfaces
        "m2": si.m**2, "mm2": si.mm**2,
        # Volumes
        "m3": si.m**3, "mm3": si.mm**3,
        # Inerties
        "m4": si.m**4, "mm4": si.mm**4,
        # Forces
        "N": si.N, "kN": si.kN,
        # Moments
        "N.m": si.N * si.m, "kN.m": si.kN * si.m,
        "N.mm": si.N * si.mm,
        # Forces lineiques
        "N/m": si.N / si.m, "kN/m": si.kN / si.m,
        "N/mm": si.N / si.mm,
        # Pressions / contraintes
        "Pa": si.Pa, "kPa": si.kPa, "MPa": si.MPa,
        "N/m2": si.N / si.m**2, "kN/m2": si.kN / si.m**2,
        "N/mm2": si.N / si.mm**2,
    }
    # Tuple d'unites propose a l'utilisateur (sert de widget combobox).
    # "none" = valeur brute retournee telle quelle (int/float/str).
    _SET_VALUE_UNITS = ("Aucune",) + tuple(_PHYSICAL_UNITS.keys())
    
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

    @staticmethod
    def interpolation_lineaire(x, xa, xb, ya, yb):
        """Fait une interpolation linéaire pour trouver un résultat y entre deux valeurs xa et xb.

        Args:
            x: Valeur pour laquelle on cherche l'interpolation.
            xa: Première borne inférieure.
            xb: Première borne supérieure.
            ya: Valeur correspondante à xa.
            yb: Valeur correspondante à xb.

        Returns:
            float: Valeur interpolée y.
        """
        y = ya + (x - xa) * ((yb - ya) / (xb - xa))
        return y

    @staticmethod
    def interpolation_logarithmique(x, xa, xb, ya, yb):
        """Fait une interpolation logarithmique pour trouver un résultat y entre deux valeurs xa et xb.

        Args:
            x: Valeur pour laquelle on cherche l'interpolation.
            xa: Première borne inférieure.
            xb: Première borne supérieure.
            ya: Valeur correspondante à xa.
            yb: Valeur correspondante à xb.

        Returns:
            float: Valeur interpolée y.
        """
        y = (
            (x > xb) * yb
            + (x < xa) * ya
            + (x <= xb) * (x >= 1) * (ya - (ya - yb) * mt.log10(x))
        )
        return y
