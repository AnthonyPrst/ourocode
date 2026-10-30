# coding in UTF-8
# by Anthony PARISOT
import os
import json
import inspect
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
        """Retourne l'objet lui même."""
        return self

    @staticmethod
    def _flag(value) -> bool:
        """Normalise un argument booléen au format pseudo-enum en ``bool``.

        La convention du package documente les options booléennes par un tuple
        ``("False", "True")`` ou ``("True", "False")`` dans les signatures, le
        premier élément étant la valeur par défaut. Comme un tuple non vide est
        toujours évalué ``True`` par Python, un défaut ``("False", "True")``
        non renseigné activait le flag à tort. Cette méthode convertit la
        valeur reçue en vrai booléen :

        - ``bool`` → retourné tel quel
        - ``"True"`` / ``"False"`` (insensible à la casse) → ``True`` / ``False``
        - ``tuple`` / ``list`` → le premier élément est normalisé (défaut)
        - ``int`` / ``float`` → ``bool(value)``
        - ``None`` → ``False``

        Args:
            value: Valeur à normaliser.

        Returns:
            bool: Le booléen équivalent.

        Raises:
            ValueError: Si la chaîne n'est pas une valeur booléenne reconnue.
        """
        if isinstance(value, bool):
            return value
        if value is None:
            return False
        if isinstance(value, (list, tuple)):
            if not value:
                return False
            return Objet._flag(value[0])
        if isinstance(value, str):
            normalized = value.strip().lower()
            if normalized in ("true", "vrai", "oui", "yes", "1"):
                return True
            if normalized in ("false", "faux", "non", "no", "0", ""):
                return False
            raise ValueError(f"Valeur booléenne non reconnue : {value!r}. Attendu : 'True' ou 'False'.")
        if isinstance(value, (int, float)):
            return bool(value)
        return bool(value)

    @staticmethod
    def _check_choice(value, choices, name: str):
        """Valide qu'une valeur appartient à une liste de choix autorisés.

        Applique la convention pseudo-enum du package : lorsque la valeur
        reçue est le tuple/liste des options lui-même (défaut de signature
        non renseigné), le premier élément est pris comme valeur par défaut,
        à l'image du comportement des combobox de l'interface graphique.

        Args:
            value: Valeur à valider.
            choices: Itérable des valeurs acceptées.
            name (str): Nom du paramètre pour le message d'erreur.

        Returns:
            La valeur validée (ou le premier choix si ``value`` était le
            conteneur d'options).

        Raises:
            ValueError: Si la valeur n'appartient pas aux choix.
        """
        # Conteneur reçu en entrée = défaut de signature pseudo-enum non
        # renseigné : le premier élément est la valeur par défaut (comportement
        # des combobox de l'interface graphique).
        if isinstance(value, dict):
            value = next(iter(value))
        elif isinstance(value, (list, tuple)) and not isinstance(value, str):
            value = value[0]
        if value not in choices:
            raise ValueError(f"{name} '{value}' invalide. Valeurs acceptées : {list(choices)}")
        return value

    @staticmethod
    def _check_dimension(value, name: str, strict: bool = True):
        """Valide une dimension géométrique (mm) fournie en valeur brute.

        Args:
            value: Valeur numérique de la dimension.
            name (str): Nom du paramètre pour le message d'erreur.
            strict (bool): True → valeur > 0 requise ; False → valeur >= 0.

        Raises:
            ValueError: Si la dimension est négative (ou nulle si strict).
        """
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise ValueError(f"{name} doit être une valeur numérique en mm, reçu : {value!r}")
        if strict and value <= 0:
            raise ValueError(f"{name} doit être strictement positif (> 0 mm), reçu : {value}")
        if not strict and value < 0:
            raise ValueError(f"{name} ne peut pas être négatif, reçu : {value}")

    def set_value(self, value: float | int | str = None, unit: str = MathUtilsMixin._SET_VALUE_UNITS):
        """Retourne une valeur typee avec une unite physique optionnelle.

        Equivalent inverse de ``get_value`` : permet a l'utilisateur de
        saisir une valeur numerique brute et de l'encapsuler dans un objet
        ``Physical`` (forallpeople) avec l'unite choisie parmi la liste
        definie dans ``_PHYSICAL_UNITS`` (m, mm, cm, km, m2, N, kN, daN,
        N.m, kN/m, Pa, MPa, etc.).

        Args:
            value (float|int|str): la valeur numerique a encapsuler. Les
                chaines numeriques ("5", "5,3") sont converties.
            unit (str): unite cible. La valeur
                speciale ``"Aucune"`` (par defaut) retourne la valeur brute
                sans conversion.

        Returns:
            - ``Physical`` si une unite valide est specifiee (ex: ``5 * si.m``)
            - valeur brute (int/float/str) si ``unit == "Aucune"`` ou unite
              inconnue

        Exemple:
            >>> obj.set_value(5, unit="m")
            5.000 m
            >>> obj.set_value("5,3", unit="kN")
            5.300 kN
            >>> obj.set_value(42, unit="Aucune")
            42
        """
        # "unit" arrive sous forme de tuple au premier chargement (valeur par
        # defaut du widget combobox) ; on normalise en chaine.
        if not isinstance(unit, str):
            unit = "Aucune"

        if value is None:
            return None

        if unit == "Aucune" or unit not in self._PHYSICAL_UNITS:
            return value

        # Conversion robuste en float (gere int/float/str numerique)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            numeric = float(value)
        elif isinstance(value, str):
            try:
                numeric = float(value.strip().replace(",", "."))
            except ValueError:
                return value
        else:
            # Valeur non convertible (Physical deja type, objet, etc.) :
            # on retourne tel quel pour ne pas ecraser le sens.
            return value

        return numeric * self._PHYSICAL_UNITS[unit]

    def get_value(
        self,
        value: dict | list | str,
        index: int | None = None,
        key: str = None,
        get_keys: bool = ("False", "True"),
    ):
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
        get_keys = self._flag(get_keys)
        if index is not None and isinstance(value, list):
            value = value[index]
        elif index is not None and isinstance(value, str):
            value = list(value)[index]
        elif get_keys and isinstance(value, dict):
            value = list(value.keys())
        elif key and isinstance(value, dict):
            value = value[key]
        elif key and isinstance(value, str):
            value = json.loads(f"{value}".replace("'", '"')).get(key)
        return value

    @classmethod
    def _from_dict(cls, dictionary: dict):
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
    def _from_parent_class(cls, objet: list | object, **kwargs):
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

        # On ne propage que les attributs qui correspondent à un paramètre
        # attendu par l'un des constructeurs de la classe cible (MRO).
        # Les attributs internes (privés), les états calculés (taux, dataframes,
        # objets modèle) et autres données recalculées par le constructeur sont
        # ainsi exclus de la propagation.
        accepted = cls._accepted_init_params()
        dict_objet = {k: v for k, v in dict_objet.items() if not k.startswith("_") and k in accepted}

        # On met à jour avec les kwargs qui écrase tout
        dict_objet.update(kwargs)
        return cls(**dict_objet)

    @classmethod
    def _accepted_init_params(cls) -> set:
        """Retourne l'ensemble des noms de paramètres nommés acceptés par les
        ``__init__`` de la hiérarchie de la classe (MRO).

        Les ``*args``/``**kwargs`` ne sont pas inclus : seuls les paramètres
        explicitement déclarés dans les signatures sont retenus pour filtrer
        la propagation des attributs par :meth:`_from_parent_class`.
        """
        params = set()
        for klass in cls.__mro__:
            try:
                signature = inspect.signature(klass.__init__)
            except (AttributeError, TypeError, ValueError):
                continue
            for name, param in signature.parameters.items():
                if param.kind in (param.POSITIONAL_OR_KEYWORD, param.KEYWORD_ONLY):
                    params.add(name)
        return params

    def _show_element(self, picture: str):
        """Affiche l'image des caractéristiques d'une entaille au cisaillement"""
        file = os.path.join(self.PATH_CATALOG, "data", "screenshot", picture)
        image = Image.open(file)
        image.show()
