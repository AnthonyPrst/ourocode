# coding in UTF-8
# by Anthony PARISOT
import os
import json
import csv
import pickle
import logging
import importlib

import forallpeople as si

si.environment("structural")

logger = logging.getLogger(__name__)


# Correspondance (ancien_module, classe) -> nouveau sous-module ourocode.eurocode.
# Permet d'ouvrir les fichiers .oco créés avec l'arborescence v1.x
# (ex : ourocode.eurocode.EC5_Element_droit.Barre).
_RELOCATED_CLASSES = {
    ("objet", "Objet"): "core.objet",
    ("A0_Projet", "Projet"): "core.projet",
    ("A0_Projet", "Batiment"): "core.batiment",
    ("A0_Projet", "Model_generator"): "core.model_generator",
    ("A0_Projet", "Wood_beam_model"): "core.model_generator",
    ("A0_Projet", "Model_result"): "core.model_result",
    ("EC0_Combinaison", "Combinaison"): "core.combinaison",
    ("EC1_Exploitation", "Exploitation"): "ec1.exploitation",
    ("EC1_Neige", "Neige"): "ec1.neige",
    ("EC1_Vent", "Vent"): "ec1.vent",
    ("EC1_Vent", "Murs_verticaux"): "ec1.vent",
    ("EC1_Vent", "Toiture_terrasse_acrotere"): "ec1.vent",
    ("EC1_Vent", "Toiture_1_pant"): "ec1.vent",
    ("EC1_Vent", "Toiture_2_pants"): "ec1.vent",
    ("EC1_Vent", "Toiture_isolee_1_pant"): "ec1.vent",
    ("EC1_Vent", "Toiture_isolee_2_pants"): "ec1.vent",
    ("EC3_Element_droit", "Plat"): "ec3.element_droit.plat",
    ("EC3_Element_droit", "Traction"): "ec3.element_droit.traction",
    ("EC3_Element_droit", "Compression"): "ec3.element_droit.compression",
    ("EC3_Element_droit", "Cisaillement"): "ec3.element_droit.cisaillement",
    ("EC3_Element_droit", "Flexion"): "ec3.element_droit.flexion",
    ("EC3_Feu", "_TemperatureGaz"): "ec3.feu",
    ("EC3_Feu", "_CoeffFeu"): "ec3.feu",
    ("EC3_Feu", "Feu_acier"): "ec3.feu",
    ("EC3_Assemblage", "Tige"): "ec3.assemblage.tige",
    ("EC3_Assemblage", "_Tige_scellement"): "ec3.assemblage.tige",
    ("EC3_Assemblage", "Tige_scellement_plaque_ancrage"): "ec3.assemblage.tige",
    ("EC3_Assemblage", "Tige_scellement_crochet"): "ec3.assemblage.tige",
    ("EC3_Assemblage", "Soudure"): "ec3.assemblage.soudure",
    ("EC3_Assemblage", "Platine_assise_compression_beton"): "ec3.assemblage.platine",
    ("EC3_Assemblage", "Platine_assise_compression_bois"): "ec3.assemblage.platine",
    ("EC3_Assemblage", "Platine_assise_traction"): "ec3.assemblage.platine",
    ("EC5_Element_droit", "Barre"): "ec5.element_droit.barre",
    ("EC5_Element_droit", "Flexion"): "ec5.element_droit.flexion",
    ("EC5_Element_droit", "Traction"): "ec5.element_droit.traction",
    ("EC5_Element_droit", "Compression"): "ec5.element_droit.compression",
    ("EC5_Element_droit", "Compression_perpendiculaire"): "ec5.element_droit.compression",
    ("EC5_Element_droit", "Compression_inclinees"): "ec5.element_droit.compression",
    ("EC5_Element_droit", "Cisaillement"): "ec5.element_droit.cisaillement",
    ("EC5_Assemblage", "Assemblage"): "ec5.assemblage.assemblage",
    ("EC5_Assemblage", "Embrevement"): "ec5.assemblage.embrevement",
    ("EC5_Assemblage", "Pointe"): "ec5.assemblage.pointe",
    ("EC5_Assemblage", "Agrafe"): "ec5.assemblage.agrafe",
    ("EC5_Assemblage", "Boulon"): "ec5.assemblage.boulon",
    ("EC5_Assemblage", "Broche"): "ec5.assemblage.broche",
    ("EC5_Assemblage", "_Tirefond"): "ec5.assemblage.tirefond",
    ("EC5_Assemblage", "Tirefond_inf_7"): "ec5.assemblage.tirefond",
    ("EC5_Assemblage", "Tirefond_sup_6"): "ec5.assemblage.tirefond",
    ("EC5_BLC", "Poutre_simple_decroissance"): "ec5.blc",
    ("EC5_CVT", "MOB"): "ec5.cvt",
    ("EC5_Feu", "Feu"): "ec5.feu.feu",
    ("EC5_Feu", "Flexion_feu"): "ec5.feu.flexion_feu",
    ("EC5_Feu", "Traction_feu"): "ec5.feu.traction_feu",
    ("EC5_Feu", "Compression_feu"): "ec5.feu.compression_feu",
    ("EC5_Feu", "Cisaillement_feu"): "ec5.feu.cisaillement_feu",
    ("EC8_Sismique", "Sismique"): "ec8.sismique",
}


class _CompatUnpickler(pickle.Unpickler):
    """Unpickler redirigeant les classes de l'ancienne arborescence v1.x.

    Les fichiers ``.oco`` créés avant la restructuration en sous-packages
    (``core``, ``ec1``, ``ec3``, ``ec5``, ``ec8``) référencent d'anciens
    modules comme ``ourocode.eurocode.EC5_Element_droit.Barre``. Cette
    classe redirige ces références vers leur nouvel emplacement.

    Note :
        La redirection est faite par nom de classe ; si des attributs ont
        été renommés entre versions, l'objet restauré peut être incomplet.
    """

    def find_class(self, module, name):
        old_module = module.split(".")[-1]
        new_sub = _RELOCATED_CLASSES.get((old_module, name))
        if new_sub:
            new_module = f"ourocode.eurocode.{new_sub}"
            logger.debug("Redirection pickle : %s.%s -> %s.%s", module, name, new_module, name)
            module = new_module
        return super().find_class(module, name)


def _file_dialog_save(filter_str: str, selected: str) -> str:
    """Ouvre un QFileDialog de sauvegarde et retourne le chemin choisi."""
    from PySide6.QtWidgets import QFileDialog

    return QFileDialog.getSaveFileName(filter=filter_str, selectedFilter=selected)[0]


def _file_dialog_open(filter_str: str, selected: str) -> str:
    """Ouvre un QFileDialog d'ouverture et retourne le chemin choisi."""
    from PySide6.QtWidgets import QFileDialog

    return QFileDialog.getOpenFileName(filter=filter_str, selectedFilter=selected)[0]


class SerializationMixin:
    """Mixin fournissant la sérialisation/désérialisation avec gestion des unités physiques.

    Méthodes :
    - _physical_to_dict / _dict_to_physical : conversion Physical ↔ dict
    - _resolve_unit_expr : parsing sécurisé d'expressions d'unité
    - save_data / load_data : sauvegarde/chargement JSON/CSV
    - save_object / _save_muliple_objects : sauvegarde .oco (JSON versionné)
    - _open_object / _open_multiple_objects : lecture .oco (JSON ou pickle legacy)

    Format .oco :
        L'écriture produit un JSON versionné ``{"_format": "oco/1", "_class": ..., "data": {...}}``.
        La lecture accepte ce format JSON ainsi que les anciens fichiers pickle
        (détection automatique au premier octet). Le pickle n'est plus écrit :
        il reste supporté en lecture seule pour la compatibilité ascendante.
    """

    _OCO_FORMAT = "oco/1"

    @staticmethod
    def _physical_to_dict(obj):
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
            return {k: SerializationMixin._physical_to_dict(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [SerializationMixin._physical_to_dict(x) for x in obj]
        elif isinstance(obj, si.Physical):
            split_vals = obj.split(base_value=True)
            return {"_physical_value": float(split_vals[0]), "_physical_unit": str(split_vals[1])}
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

    @staticmethod
    def _dict_to_physical(data):
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
            if "_physical_value" in data:
                # Reconstruire l'objet Physical
                value = data["_physical_value"]
                unit_str = str(data["_physical_unit"].split()[-1])
                # Normaliser quelques notations possibles
                unit_expr = (
                    unit_str.replace("^", "**")
                    .replace("²", "**2")
                    .replace("³", "**3")
                    .replace("⁴", "**4")
                    .replace("·", "*")
                )
                if "/" in unit_expr:
                    unit_expr = unit_expr.replace("/", "_")
                try:
                    unit_obj = SerializationMixin._resolve_unit_expr(unit_expr)
                except ValueError as e:
                    raise ValueError(f"Unité inconnue ou non prise en charge: {unit_str}") from e
                return value * unit_obj
            return {k: SerializationMixin._dict_to_physical(v) for k, v in data.items()}
        elif isinstance(data, (list, tuple)):
            return [SerializationMixin._dict_to_physical(x) for x in data]
        return data

    def save_data(self, data: dict, type_data: str = ("JSON", "CSV"), path: str = None):
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
            save_file_path = (
                path
                if path
                else QFileDialog.getSaveFileName(
                    filter="JSON (*.json)",
                    selectedFilter=".json",
                )[0]
            )
            with open(save_file_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=4)
        elif type_data == "CSV":
            save_file_path = (
                path
                if path
                else QFileDialog.getSaveFileName(
                    filter="CSV (*.csv)",
                    selectedFilter=".csv",
                )[0]
            )
            with open(save_file_path, "w", newline="") as f:
                w = csv.DictWriter(f, data.keys())
                w.writeheader()
                w.writerow(data)

    def load_data(self, type_data: str = ("JSON", "CSV"), path: str = None):
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
            file_path = (
                path
                if path
                else QFileDialog.getOpenFileName(
                    filter="JSON (*.json)",
                    selectedFilter=".json",
                )[0]
            )
            with open(file_path, "r", encoding="utf-8") as f:
                return self._dict_to_physical(json.load(f))
        elif type_data == "CSV":
            file_path = (
                path
                if path
                else QFileDialog.getOpenFileName(
                    filter="CSV (*.csv)",
                    selectedFilter=".csv",
                )[0]
            )
            with open(file_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                return self._dict_to_physical({row["key"]: row["value"] for row in reader})

    # ------------------------------------------------------------------
    # Persistance .oco — format JSON versionné "oco/1"
    # ------------------------------------------------------------------

    @staticmethod
    def _oco_serialize(obj):
        """Sérialise récursivement un objet pour le format .oco JSON.

        Gère les ``Physical``, les conteneurs et les objets ourocode imbriqués
        (sérialisés par référence de classe + données)."""
        if isinstance(obj, dict):
            return {k: SerializationMixin._oco_serialize(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [SerializationMixin._oco_serialize(x) for x in obj]
        if isinstance(obj, si.Physical):
            return SerializationMixin._physical_to_dict(obj)
        if isinstance(obj, SerializationMixin):
            return {
                "_oco_object": True,
                "_class": f"{type(obj).__module__}.{type(obj).__qualname__}",
                "data": SerializationMixin._oco_serialize(obj.__dict__),
            }
        return obj

    @staticmethod
    def _oco_deserialize(data):
        """Reconstruit les objets ``Physical`` et ourocode imbriqués d'un .oco JSON."""
        if isinstance(data, dict):
            if data.get("_oco_object"):
                cls = SerializationMixin._resolve_oco_class(data["_class"])
                obj = cls.__new__(cls)
                obj.__dict__.update(SerializationMixin._oco_deserialize(data["data"]))
                return obj
            if "_physical_value" in data:
                return SerializationMixin._dict_to_physical(data)
            return {k: SerializationMixin._oco_deserialize(v) for k, v in data.items()}
        if isinstance(data, list):
            return [SerializationMixin._oco_deserialize(x) for x in data]
        return SerializationMixin._dict_to_physical(data)

    @staticmethod
    def _resolve_oco_class(class_path: str):
        """Résout une classe depuis son chemin ``module.QualName``.

        Seules les classes du package ``ourocode`` sont acceptées — un fichier
        .oco ne peut pas instancier de classe arbitraire."""
        module_name, _, class_name = class_path.rpartition(".")
        if not module_name.startswith("ourocode"):
            raise ValueError(f"Classe .oco refusée (hors package ourocode) : {class_path}")
        target = importlib.import_module(module_name)
        for part in class_name.split("."):
            target = getattr(target, part)
        return target

    @classmethod
    def _oco_payload(cls, obj) -> dict:
        """Construit le payload JSON d'un objet ourocode.

        Les attributs non sérialisables en JSON (DataFrames, figures, handles
        Qt…) sont ignorés avec un avertissement plutôt que de faire échouer
        toute la sauvegarde."""
        data = {}
        for key, value in obj.__dict__.items():
            try:
                serialized = cls._oco_serialize(value)
                json.dumps(serialized)
            except (TypeError, ValueError):
                logger.warning(
                    "Attribut %r de %s non sérialisable en JSON : ignoré à la sauvegarde .oco.",
                    key,
                    type(obj).__name__,
                )
                continue
            data[key] = serialized
        return {
            "_format": cls._OCO_FORMAT,
            "_class": f"{type(obj).__module__}.{type(obj).__qualname__}",
            "data": data,
        }

    @classmethod
    def _oco_load(cls, payload: dict):
        """Instancie un objet depuis un payload .oco JSON (sans appeler __init__)."""
        if payload.get("_format") != cls._OCO_FORMAT:
            raise ValueError(f"Format .oco non supporté : {payload.get('_format')!r}")
        cls_obj = cls._resolve_oco_class(payload["_class"])
        obj = cls_obj.__new__(cls_obj)
        obj.__dict__.update(cls._oco_deserialize(payload["data"]))
        return obj

    @staticmethod
    def _oco_is_json(file_path: str) -> bool:
        """Détecte le format du fichier : JSON (``{``) ou pickle legacy."""
        with open(file_path, "rb") as f:
            return f.read(64).lstrip()[:1] == b"{"

    def _save_muliple_objects(self, object: list, path: str = None):
        """Sauvegarde une liste d'objets dans un fichier .oco (JSON versionné).

        Args:
            object (list): Liste des objets ourocode à sérialiser.
            path (str, optional): Chemin du fichier .oco à créer. Si None,
                une boîte de dialogue Qt s'ouvre pour choisir l'emplacement.
        """
        save_file_path = (
            path
            if path
            else _file_dialog_save(
                "Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                ".oco",
            )
        )
        if not save_file_path:
            return
        payload = {
            "_format": self._OCO_FORMAT,
            "_objects": [self._oco_payload(obj) for obj in object],
        }
        with open(save_file_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    def save_object(self, path: str = None):
        """Sauvegarde l'objet courant dans un fichier .oco (JSON versionné).

        Format ``oco/1`` : ``{"_format", "_class", "data"}`` où ``data`` est
        l'état ``__dict__`` sérialisé (Physical → dict valeur/unité, objets
        ourocode imbriqués → référence de classe). Les attributs non
        sérialisables en JSON sont ignorés avec un avertissement.

        Args:
            path (str, optional): Chemin du fichier .oco à créer. Si None,
                une boîte de dialogue Qt s'ouvre pour choisir l'emplacement.
        """
        save_file_path = (
            path
            if path
            else _file_dialog_save(
                "Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                ".oco",
            )
        )
        if not save_file_path:
            return
        with open(save_file_path, "w", encoding="utf-8") as f:
            json.dump(self._oco_payload(self), f, indent=2)

    @classmethod
    def _open_multiple_objects(cls, path: str = None):
        """Charge une liste d'objets depuis un fichier .oco.

        Le format JSON ``oco/1`` est lu nativement ; les fichiers pickle
        (v1.x et pickles .oco) sont automatiquement détectés et chargés via
        ``_CompatUnpickler``.

        Args:
            path (str, optional): Chemin du fichier .oco à charger. Si None,
                une boîte de dialogue Qt s'ouvre pour sélectionner le fichier.

        Returns:
            list: Liste des objets désérialisés.

        Warning:
            Les fichiers pickle hérités sont chargés avec ``pickle.load`` qui
            peut exécuter du code arbitraire. N'ouvrez que des fichiers .oco
            dont vous maîtrisez la provenance.
        """
        file_path = (
            path
            if path
            else _file_dialog_open(
                "Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                ".oco",
            )
        )
        if cls._oco_is_json(file_path):
            with open(file_path, "r", encoding="utf-8") as f:
                payload = json.load(f)
            if "_objects" in payload:
                return [cls._oco_load(p) for p in payload["_objects"]]
            return [cls._oco_load(payload)]
        data = []
        with open(file_path, "rb") as f:
            while True:
                try:
                    data.append(_CompatUnpickler(f).load())
                except EOFError:
                    break
        return data

    @classmethod
    def _open_object(cls, path: str = None):
        """Charge un objet depuis un fichier .oco.

        Le format JSON ``oco/1`` est lu nativement ; les fichiers pickle
        (v1.x et pickles .oco) sont automatiquement détectés et chargés via
        ``_CompatUnpickler``.

        Args:
            path (str, optional): Chemin du fichier .oco à charger. Si None,
                une boîte de dialogue Qt s'ouvre pour sélectionner le fichier.

        Returns:
            object: L'objet désérialisé.

        Warning:
            Les fichiers pickle hérités sont chargés avec ``pickle.load`` qui
            peut exécuter du code arbitraire. N'ouvrez que des fichiers .oco
            dont vous maîtrisez la provenance.
        """
        file_path = (
            path
            if path
            else _file_dialog_open(
                "Ourea catalog object (*.oco);;'Text Document' (*.txt)",
                ".oco",
            )
        )
        if cls._oco_is_json(file_path):
            with open(file_path, "r", encoding="utf-8") as f:
                return cls._oco_load(json.load(f))
        with open(file_path, "rb") as f:
            return _CompatUnpickler(f).load()
