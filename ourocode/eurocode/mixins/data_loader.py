# coding in UTF-8
# by Anthony PARISOT
import os
import json

import pandas as pd


def get_package_path(package):
    # Obtenir le chemin du fichier principal du module
    package_path = os.path.dirname(package.__file__)
    return package_path


class DataLoaderMixin:
    """Mixin fournissant le chargement des données normatives (CSV, JSON).

    Attributs de classe :
    - PATH_CATALOG : chemin racine du package ourocode
    - _csv_cache / _json_cache : caches partagés

    Méthodes :
    - _data_from_csv : charge un CSV avec cache
    - _data_from_json : charge un JSON en DataFrame
    - _load_json : charge un JSON en dict avec cache
    """

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
        if cache_key not in DataLoaderMixin._csv_cache:
            repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
            DataLoaderMixin._csv_cache[cache_key] = pd.read_csv(repertory, sep=';', header=0, index_col=index_col)
        return DataLoaderMixin._csv_cache[cache_key].copy()
    
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
        if data_file not in DataLoaderMixin._json_cache:
            repertory = os.path.join(self.PATH_CATALOG, "data", data_file)
            with open(repertory, "r", encoding="utf-8") as json_file:
                DataLoaderMixin._json_cache[data_file] = json.load(json_file)
        return DataLoaderMixin._json_cache[data_file]
