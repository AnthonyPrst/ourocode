# coding in UTF-8
# by Anthony PARISOT
import os
import json
import csv
import pickle

import forallpeople as si
si.environment("structural")


class SerializationMixin:
    """Mixin fournissant la sérialisation/désérialisation avec gestion des unités physiques.

    Méthodes :
    - _physical_to_dict / _dict_to_physical : conversion Physical ↔ dict
    - _resolve_unit_expr : parsing sécurisé d'expressions d'unité
    - save_data / load_data : sauvegarde/chargement JSON/CSV
    - save_object / _save_muliple_objects : pickle .oco
    - _open_object / _open_multiple_objects : lecture pickle .oco
    """

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
            return SerializationMixin._resolve_unit_expr(parts[0]) * SerializationMixin._resolve_unit_expr(parts[1])
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
