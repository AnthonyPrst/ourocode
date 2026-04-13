# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/),
et ce projet adhère au [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

### Changed
- **Refactoring Phase 1 — Restructuration complète du package** :
  - Création des sous-packages `core/`, `ec1/`, `ec3/`, `ec5/`, `ec8/` sous `ourocode/eurocode/`.
  - `A0_Projet.py` éclaté en `core/projet.py`, `core/batiment.py`, `core/model_generator.py`, `core/model_result.py`.
  - `EC0_Combinaison.py` → `core/combinaison.py`.
  - `objet.py` → `core/objet.py`.
  - `EC1_Neige.py`, `EC1_Vent.py`, `EC1_Exploitation.py` → `ec1/neige.py`, `ec1/vent.py`, `ec1/exploitation.py`.
  - `EC3_Element_droit.py` éclaté en `ec3/element_droit/plat.py`, `traction.py`, `compression.py`, `cisaillement.py`, `flexion.py`.
  - `EC3_Assemblage.py` éclaté en `ec3/assemblage/tige.py`, `soudure.py`, `platine.py`.
  - `EC3_Feu.py` → `ec3/feu.py`.
  - `EC5_Element_droit.py` éclaté en `ec5/element_droit/barre.py`, `flexion.py`, `traction.py`, `compression.py`, `cisaillement.py`.
  - `EC5_Assemblage.py` éclaté en `ec5/assemblage/assemblage.py`, `embrevement.py`, `pointe.py`, `agrafe.py`, `boulon.py`, `broche.py`, `tirefond.py`.
  - `EC5_Feu.py` éclaté en `ec5/feu/feu.py`, `flexion_feu.py`, `traction_feu.py`, `compression_feu.py`, `cisaillement_feu.py`.
  - `EC5_BLC.py` → `ec5/blc.py`, `EC5_CVT.py` → `ec5/cvt.py`.
  - `EC8_Sismique.py` → `ec8/sismique.py`.
  - Suppression de tous les anciens fichiers monolithiques.
  - Mise à jour de tous les imports dans les tests et la documentation mkdocs.

### Added
- Documentation API complète avec MkDocs-Material + mkdocstrings (déploiement GitHub Pages automatique).
- Job `docs` dans le workflow CI : déploiement automatique sur `gh-pages` après passage des tests.
- Badge "Documentation" dans le README avec lien vers GitHub Pages.
- Section "Documentation" dans le README avec liens vers le guide, la référence API et les normes.
- Module `EC5_BLC.py` : paramètres spécifiques au bois lamellé-collé (EN 1995-1-1).
- Dépendances optionnelles `[docs]` dans `pyproject.toml` (`mkdocs-material`, `mkdocstrings[python]`, `mkdocs-include-markdown-plugin`, `mike`).

### Fixed
- Import conditionnel de `IPython` dans `core/objet.py` pour éviter une erreur si le module n'est pas installé.
- Correction du paramètre `section` manquant dans les tests CVT.

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
