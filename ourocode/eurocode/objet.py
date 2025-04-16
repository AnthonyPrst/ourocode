# coding in UTF-8 
# by Anthony PARISOT

import os
import json
import math as mt
import importlib.resources as pkg_resources
from PIL import Image
import pandas as pd
import pickle
import inspect
from PySide6.QtWidgets import QFileDialog
# from tkinter import filedialog

import forallpeople as si
si.environment("structural")


def get_package_path(package):
    # Obtenir le chemin du fichier principal du module
    package_path = os.path.dirname(package.__file__)
    return package_path


class Objet(object):
    """Classe permetant la sauvegarde ou l'ouverture d'un objet ou de plusieur sous un fichier .ec
    """
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
    
    
    @property
    def objet(self):
        """Retourne l'objet lui même.
        """
        return self

    def get_value(self, value: str, index: int=None, key: str=None):
        """Retourne l'argument transmit.

        Args:
            index (int, optional): index à retourner dans une liste python. 
                Attention sous pyhon le premier élément d'une liste ce trouve à l'index 0.
            key (str, optional): clé à renvoyer dans un dictionnaire python.
        """
        if index:
            value = list(value)[index]
        elif key:
            value = f"{value.replace("'", '"')}"
            value = json.loads(value)[key]
        return str(value)
        
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
    
    def cos(self, value: float):
        """Retourne le cosinus de la valeur donnée en degré."""
        return mt.cos(mt.radians(float(value)))
    
    def sin(self, value: float):
        """Retourne le sinus de la valeur donnée en degré."""
        return mt.sin(mt.radians(float(value)))
    
    def tan(self, value: float):
        """Retourne la tangente de la valeur donnée en degré."""
        return mt.tan(mt.radians(float(value)))
    

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
                elif si_unit == str(si.m**2):
                    if unit_to_convert == str(si.mm**2):
                        return value * 10**6
                elif si_unit == str(si.m**3):
                    if unit_to_convert == str(si.mm**3):
                        return value * 10**9
                elif si_unit == str(si.m**4):
                    if unit_to_convert == str(si.mm**4):
                        return value * 10**12
                elif si_unit == str(si.N):
                    if unit_to_convert == str(si.kN):
                        return value * 10**-3
                elif si_unit == str(si.Pa):
                    if unit_to_convert == str(si.MPa):
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
        """Class méthode permetant l'intanciation des classe hérité de la classe parent, par une classe déjà instanciée.

        Args:
            object (class object): l'objet Element déjà créer par l'utilisateur
        """
        dict_objet = {}
        if isinstance(objet, list):
            for obj in objet:
                dict_objet.update(obj.__dict__)
                dict_objet.update(Objet._reset_physical_object(obj))
        else:
            dict_objet = objet.__dict__
            dict_objet.update(Objet._reset_physical_object(objet))
        return cls(**dict_objet, **kwargs)

    
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
