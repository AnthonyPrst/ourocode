# coding in UTF-8 
# by Anthony PARISOT
import os
import json
from PIL import Image
try:
    from IPython.display import display, Latex
except ImportError:
    display = None
    Latex = None

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.mixins.serialization import SerializationMixin
from ourocode.eurocode.mixins.data_loader import DataLoaderMixin
from ourocode.eurocode.mixins.math_utils import MathUtilsMixin
from ourocode.eurocode.mixins.synthese import SyntheseMixin


class Objet(SerializationMixin, DataLoaderMixin, MathUtilsMixin, SyntheseMixin):
    """Classe de base du package ourocode.

    Cette classe fournit les fonctionnalités essentielles pour :
    - La sérialisation/désérialisation des objets avec gestion des unités physiques
      (via SerializationMixin)
    - L'accès aux données normatives CSV/JSON (via DataLoaderMixin)
    - Les opérations mathématiques et de conversion d'unités (via MathUtilsMixin)
    - La synthèse des taux de travail (via SyntheseMixin)

    Toutes les classes du package ourocode héritent de cette classe.
    """
    JUPYTER_DISPLAY = False
    

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

    
    def _show_element(self, picture: str):
        """Affiche l'image des caractéristiques d'une entaille au cisaillement
        """
        file = os.path.join(self.PATH_CATALOG, "data", "screenshot", picture)
        image = Image.open(file)
        image.show()
