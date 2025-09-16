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
from PySide6.QtWidgets import QFileDialog

import forallpeople as si
si.environment("structural")


def get_package_path(package):
    # Obtenir le chemin du fichier principal du module
    package_path = os.path.dirname(package.__file__)
    return package_path


class Objet(object):
    """Classe permetant la sauvegarde ou l'ouverture d'un objet ou de plusieur sous un fichier .ec
    """
    JUPYTER_DISPLAY = False
    OPERATOR = ("+", "-", "x", "/")
    try:
        import ourocode
        PATH_CATALOG = os.path.join(get_package_path(ourocode))
    except:
        PATH_CATALOG = os.path.join(os.getcwd(), "ourocode")

    def _data_from_csv(self, data_file: str, index_col=0):
        """ Retourne un dataframe d'un fichier CSV """
        repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
        data_csv = pd.read_csv(repertory, sep=';', header=0, index_col=index_col)
        return data_csv
    
    def _data_from_json(self, data_file: str):
        """ Retourne un dataframe d'un fichier JSON """
        repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
        data_json = pd.read_json(repertory)
        return data_json

    def _load_json(self, data_file: str):
        """ Retourne un dict d'un fichier JSON """
        repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
        with open(repertory, "r", encoding="utf-8") as json_file:
            data = json.load(json_file)
        return data

    def _assign_handcalcs_value(self, handcalc_value: tuple, args: list[str]):
        """Assigne les valeurs des calculs handcalc aux arguments de l'objet.
        Les arguments doivent être dans le même ordre que les valeurs des calculs handcalc."""
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
        """Convertit un objet Physical en dictionnaire sérialisable."""
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

    def _dict_to_physical(self, data):
        """Reconstruit un objet Physical à partir d'un dictionnaire."""
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
                # Évaluer l'expression d'unité dans l'espace de noms de forallpeople
                try:
                    unit_obj = eval(unit_expr, {"__builtins__": {}}, vars(si))
                except Exception:
                    # En dernier recours, tenter un getattr direct si c'est un symbole simple
                    try:
                        unit_obj = getattr(si, unit_expr)
                    except Exception as e:
                        raise ValueError(f"Unité inconnue ou non prise en charge: {unit_str}") from e
                return value * unit_obj
            return {k: self._dict_to_physical(v) for k, v in data.items()}
        elif isinstance(data, (list, tuple)):
            return [self._dict_to_physical(x) for x in data]
        return data

    def get_value(self, value: dict|list|str, index: int=None, key: str=None, get_keys: bool=("False", "True"),):
        """Retourne l'argument transmit.

        Args:
            value (dict|list|str): la valeur à retourner.
            index (int, optional): index à retourner dans une liste python. 
                Attention sous pyhon le premier élément d'une liste ce trouve à l'index 0.
            key (str, optional): clé à renvoyer dans un dictionnaire python.
            get_keys (bool, optional): permet de retourner les clés d'un dictionnaire python.
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
        """Retourne l'opération donnée entre la valeur 1 et la valeur 2.
        Pour les calculs trigonométrique, la valeur 2 est en degré."""
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
        """Retourne la valeur absolue.
        """
        return abs(float(value))

    
    def max(self, value1: float, value2: float):
        """Retourne la valeur max entre la valeur 1 et valeur 2.
        """
        return max(float(value1), float(value2))

    def min(self, value1: float, value2: float):
        """Retourne la valeur min entre la valeur 1 et valeur 2.
        """
        return min(float(value1), float(value2))
    
    def _extract_numbers(self, value, absolute: bool=True):
        """Fonction récursive pour extraire tous les nombres d'une structure de données imbriquée."""
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
        """Retourne la valeur max d'une liste, ou d'un dictionnaire.
        
        Args:
            iterable (dict|list|tuple): la liste ou le dictionnaire à parcourir.
            absolute (bool, optional): permet de retourner la valeur absolue. Defaults to False.
        """
        result = self.get_value(iterable, get_keys=False)
        result = self._extract_numbers(result, absolute)  
        return max(result)
    

    def min_list(self, iterable: dict|list|tuple|str, absolute: bool=("False", "True")):
        """Retourne la valeur min d'une liste, ou d'un dictionnaire.
        
        Args:
            iterable (dict|list|tuple|str): la liste ou le dictionnaire à parcourir.
            absolute (bool, optional): permet de retourner la valeur absolue. Defaults to False.
        """
        result = self.get_value(iterable, get_keys=False)
        result = self._extract_numbers(result, absolute)  
        return min(result)

    def get_trigonometric_value(self, value: float, operator: str=("COS", "SIN", "TAN", "ACOS", "ASIN", "ATAN")):
        """Retourne la valeur trigonométrique donnée en degré."""
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
    
    def save_data(self, data: dict, type_data: str=("JSON", "CSV"), path: str=None):
        """Sauvegarde les données dans un fichier JSON ou CSV.

        Args:
            type_data (str): le type de données à sauvegarder (JSON, CSV).
            data (dict): les données à sauvegarder sous forme de dictionnaire.
            path (str, optional): Chemin du fichier à créer, s'il n'est pas fourni, une boite de dialogue s'ouvre pour choisir le fichier. 
            Defaults to None.
        """
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
        """Charge les données depuis un fichier JSON ou CSV.

        Args:
            type_data (str): le type de données à charger (JSON, CSV).
            path (str, optional): Chemin du fichier à charger, s'il n'est pas fourni, une boite de dialogue s'ouvre pour choisir le fichier. 
            Defaults to None.
        """
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
            """Convertie l'unité de base dans l'unité nécessaire à l'instanciation de la classe parente

            Args:
                value (int|float): valeur à convertir
                si_unit (si.Physical): unité si de base
                unit_to_convert (si.Physical): unité dans la quelle convertir
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
        """Class méthode permetant de réinitialiser les valeurs physiques d'un dictionnaire d'argument d'une classe parent.

        Args:
            objet (object): l'objet à réinitailiser
            dictionnary (dict): le dictionnaire d'argument de la classe parent

        Returns:
            dict: le dictionnaire d'argument de la classe parent avec les valeurs physiques réinitialisées
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
        """Class méthode permetant de réinitialiser les valeurs physiques d'un objet d'une classe parent.
        """
        dictionnary = objet.__dict__
        return cls._reset_physical_dictionnary(objet, dictionnary)
    
    

    @classmethod
    def _from_dict(cls, dictionary:dict):
        """Class méthode permetant l'intanciation des classe hérité de la classe parent, par une classe déjà instanciée.

        Args:
            object (class object): l'objet Element déjà créer par l'utilisateur
        """ 
        return cls(**dictionary)
    
    @classmethod
    def _from_parent_class(cls, objet: list|object, **kwargs):
        """Class méthode permetant l'intanciation des classes héritées de la classe parent, par une classe déjà instanciée.
        
        Les clés dans kwargs écrasent les clés existantes dans dict_objet.
        
        Args:
            objet (object|list): L'objet ou la liste d'objets à partir desquels créer la nouvelle instance
            **kwargs: Arguments additionnels qui écraseront les attributs des objets sources
            
        Returns:
            Une nouvelle instance de la classe avec les attributs des objets sources et des kwargs
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
        save_file_path = QFileDialog.getSaveFileName(
                filter="Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                selectedFilter=".oco",
            )[0]
        # with filedialog.asksaveasfile('wb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
        with open(save_file_path, "wb") as f:
            for ligne in object:
                pickle.dump(ligne, f)
    
    
    def save_object(self):
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
        if not path:
            file_path = QFileDialog.getOpenFileName(
                    filter="Ourea catalog object (*.oco);;'Text Document' (*.txt)", selectedFilter=".oco"
                )[0]
            with open(file_path, "rb") as f:
                return pickle.load(f)
        else:
            with open(path, mode="rb") as f:
                return pickle.load(f)
