# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/),
et ce projet adhère au [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

### Added
- Documentation API complète avec MkDocs-Material + mkdocstrings (déploiement GitHub Pages automatique).
- Job `docs` dans le workflow CI : déploiement automatique sur `gh-pages` après passage des tests.
- Badge "Documentation" dans le README avec lien vers GitHub Pages.
- Section "Documentation" dans le README avec liens vers le guide, la référence API et les normes.
- Module `EC5_BLC.py` : paramètres spécifiques au bois lamellé-collé (EN 1995-1-1).
- Dépendances optionnelles `[docs]` dans `pyproject.toml` (`mkdocs-material`, `mkdocstrings[python]`, `mkdocs-include-markdown-plugin`, `mike`).

## [1.10.0] - 2025

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
