# coding in UTF-8 
# by Anthony PARISOT

import os
from PIL import Image
import pandas as pd
import pickle
from tkinter import filedialog

import forallpeople as si
si.environment("structural")


class Objet(object):
    """Classe permetant la sauvegarde ou l'ouverture d'un objet ou de plusieur sous un fichier .ec
    """
    try:
        from api.base.constants import EXECUTABLE_PATH
        PATH_CATALOG = os.path.join(EXECUTABLE_PATH, "catalog")
    except:
        PATH_CATALOG = os.path.join(os.getcwd(), "catalog")
    def __init__(self) -> None:
        pass

    def _data_from_csv(self, data_file: str, index_col=0):
            """ Retourne un dataframe d'un fichier CSV """
            repertory = os.path.join(__class__.PATH_CATALOG, "data", data_file)
            data_csv = pd.read_csv(repertory, sep=';', header=0, index_col=index_col)
            return data_csv
    
    @property
    def objet(self):
        """Retourne l'objet lui même.
        """
        return self

    def return_value(self, value: str):
        """Retourne l'argument donnée.
        """
        print(value)
        return value
    

    @classmethod
    def _from_dict(cls, dictionary:dict):
        """Class méthode permetant l'intanciation des classe hérité de la classe parent, par une classe déjà instanciée.

        Args:
            object (class object): l'objet Element déjà créer par l'utilisateur
        """ 
        # dict_physical = {}
        # # Si on argument utilise forallpeople on récupère que la valeur pour ne pas multiplier l'unité par elle même
        # for key, val in dictionary.items():
        #     if isinstance(val, si.Physical):
        #         physical = val.split(base_value=False)
        #         dict_physical[key] = physical[0]

        # dictionary.update(dict_physical)
        return cls(**dictionary)
    
    @classmethod
    def _from_parent_class(cls, objet: list|object, **kwargs):
        """Class méthode permetant l'intanciation des classe hérité de la classe parent, par une classe déjà instanciée.

        Args:
            object (class object): l'objet Element déjà créer par l'utilisateur
        """
        def reset_physical(dictionnary):
            dict_physical = {}
            # Si un argument utilise forallpeople on récupère que la valeur pour ne pas multiplier l'unité par elle même
            for key, val in dictionnary.items():
                if isinstance(val, si.Physical):
                    physical = val.split(base_value=False)
                    dict_physical[key] = physical[0]

            dict_objet.update(dict_physical)


        dict_objet = {}
        if isinstance(objet, list):
            for obj in objet:
                dict_objet.update(obj.__dict__)
                reset_physical(obj.__dict__)
        else:
            dict_objet = objet.__dict__
            reset_physical(dict_objet)
            
        # print("dict_objet :", dict_objet)
        return cls(**dict_objet, **kwargs)

    
    def _save_muliple_objects(self, object: list):
        with filedialog.asksaveasfile('wb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
            for ligne in object:
                pickle.dump(ligne, f)
    
    def save_object(self):
        with filedialog.asksaveasfile('wb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
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
        with filedialog.askopenfile('rb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
            while True:
                try:
                    data.append(pickle.load(f))
                except EOFError:
                    break
            return data
    
    @classmethod
    def _open_object(cls):
        with filedialog.askopenfile('rb', filetypes=(("Ourea catalog object", "*.oco"), ('Text Document', '*.txt')), defaultextension='.oco') as f:
            return pickle.load(f)
