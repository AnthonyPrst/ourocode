# coding in UTF-8
# by Anthony PARISOT
import forallpeople as si

from ourocode.eurocode.core.objet import Objet


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
