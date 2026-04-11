# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/),
et ce projet adhère au [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

### Changed
- PySide6, PyNiteFEA et pyvista sont désormais des dépendances optionnelles (`pip install "ourocode[full]"`).
- Les messages utilisateur (`print`) convertis en `warnings.warn` dans EC5_Element_droit, EC5_Feu, EC5_Assemblage.
- Import wildcard `from math import *` supprimé dans EC1_Vent.py.

### Fixed
- Suppression des `print` de debug dans EC5_Element_droit, EC1_Vent, EC3_Feu, EC5_Assemblage.
- Correction du shadow de la variable `si` (module forallpeople) dans `objet.py` `_add_synthese_taux_travail`.
- Remplacement du `bare except` par `except ImportError` dans `objet.py`.
- Correction des badges du README (licence, CI, release).

### Added
- CI GitHub Actions : lint (ruff) + tests (pytest) sur push/PR vers main.
- CHANGELOG.md.
- Documentation des dépendances optionnelles dans le README.

## [1.9.0] - 2025

### Added
- Module EC3_Feu : calcul température acier au feu selon EN 1993-1-2.
- Module EC8_Sismique : spectre élastique selon EN 1998-1.
