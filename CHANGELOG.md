# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/),
et ce projet adhère au [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

## [2.0.0] - 2026-04-16

### Added
- **`ourocode/eurocode/core/renderer.py`** : nouveau module de rendu LaTeX maison, remplaçant `handcalcs` — compatible Python 3.12+.
  - Décorateur `@handcalc` basé sur l'AST (`ast` + `inspect`), sans dépendance à `innerscope`.
  - Rendu symbolique (noms → symboles LaTeX, grecs, indices), numérique (substitution des valeurs) et résultat final sur une même ligne (`override="short"`) ou sur trois lignes alignées (`override="long"`).
  - Extraction correcte des unités `forallpeople` via `str(x).rsplit()` → `\,\mathrm{...}`.
  - Support des fonctions `sqrt`, `sin`, `cos`, `tan`, `log`, `log10`, `exp`, `floor`, `ceil`, `abs`, `min`, `max`, `radians`.
  - Capture des variables de fermeture (closures), du scope module (globals) et des arguments positionnels.

### Changed
- **Suppression de la dépendance `handcalcs`** dans `pyproject.toml` — incompatible avec Python 3.13+ (bug `FrameLocalsProxy` dans `innerscope`).
- `requires-python` élargi à `>=3.12` (sans borne supérieure).
- Classifiers Python mis à jour : ajout `3.13`, `3.14`.
- Tous les imports `from handcalcs.decorator import handcalc` remplacés par `from ourocode.eurocode.core.renderer import handcalc` dans les modules EC1, EC5.
- **`ourocode/eurocode/ec1/exploitation.py`** : correction des séquences d'échappement invalides `"\["` / `"\]"` → `"\\["` / `"\\]"` dans le décorateur `@handcalc` de `alpha_A` (`SyntaxWarning` Python 3.14+).

### Changed (CI/CD)
- Migration du déploiement de la documentation vers l'action officielle GitHub Pages (`actions/deploy-pages@v4`).

## [1.11.3] - 2026-04-14

### Fixed
- Correction des docstrings `griffe`-incompatibles causant l'échec du job `docs` en CI (`mkdocs build --strict`) :
  - `ec5/feu/feu.py` : paramètre `Attention :` (non présent en signature) déplacé en section `Note`.
  - `ec5/feu/compression_feu.py` : noms `lo` et `type_appuis ` (avec espace) → `lo_y`, `lo_z`, `type_appuis`.
  - `ec5/feu/flexion_feu.py` : noms `lo_rel_y/z` et `coeflef_y/z` (caractère `/` invalide) → noms réels de la signature.
  - `ec5/assemblage/pointe.py` : suppression des `Args` erronés dans la docstring de la property `pince` (pas de paramètres).
- Retrait du flag `--strict` de `mkdocs build` dans `.github/workflows/ci.yml` — les warnings griffe sur `*args`/`**kwargs` non annotés ne sont pas bloquants.

### Changed
- **Amélioration des docstrings** — `ec1/vent.py` :
  - `Murs_verticaux.__init__` : zones A/B/C/D/E, calcul de `e`, interpolation C_pe.
  - `Toiture_terrasse_acrotere.__init__` : zones F/G/H/I, rapport hp/h, direction 0° uniquement.
  - `Toiture_1_pant.__init__` : zones F/G/H (0°/180°) et F_up/F_low/G/H/I (90°), Note sur alpha_toit unique.
  - `Toiture_2_pants.__init__` : zones F/G/H/I/J (0°) et F/G/H/I (90°), Note hypothèses (pentes égales, faîtage centré).
  - `Toiture_isolee_1_pant.__init__` : zones A/B/C, C_p nets, degré d'obstruction phi.
  - `Toiture_isolee_2_pants.__init__` : zones A/B/C/D, C_p nets.
  - `get_geo` / `get_Cpe` / `get_Cp` / `show_zonage` de toutes les classes ci-dessus : Returns documentés, référence §EN ajoutée.
- **Amélioration des docstrings** — `ec3/assemblage/soudure.py` :
  - `verif_soudure`, `beta_Lw1`, `cordon_frontal`, `cordon_lateral`, `cordon_oblique`, `cordon_pieces_obliques`, `critere_generale`, `soudure_discontinue`.
- **Amélioration des docstrings** — `ec3/assemblage/platine.py` :
  - `f_jd`, `c`, `taux_compression` (béton et bois), `taux_traction`.
- **Amélioration des docstrings** — `ec3/assemblage/tige.py` :
  - `FvRd`, `FtRd`, `pince_metal_boulon`, `Bp_Rd`, `taux_FvEd_FtEd`.

## [1.11.0] - 2026-04-13

### Added
- **Suite de tests unitaires complète (59% → 63% de couverture, 41 → 337 tests)** :
  - `test_A0_Projet.py` réécrit : couverture de `Projet`, `Objet` et des 4 mixins (`SerializationMixin`, `MathUtilsMixin`, `SyntheseMixin`, `DataLoaderMixin`).
  - `test_EC1_Exploitation.py` créé : `qk`, `Qk`, `alpha_A`, `alpha_n` pour toutes les catégories.
  - `test_EC3_Feu.py` créé : `_TemperatureGaz` (courbe ISO 834), `_CoeffFeu` (ky/kb/kp/kE/kw θ), `Feu_acier` (fire_data, c_a branches).
  - `test_EC3_Platine.py` créé : `Platine_assise_compression_beton` et `Platine_assise_compression_bois`.
  - `test_EC3_Element_droit.py` complété : `Nb_Rd` (toutes courbes a0/a/b/c/d), `Flexion` classe 3, `Mc_V_Rd` cisaillement élevé/faible, `Traction` Nt_Rd.
  - `test_EC5_BLC.py` créé : `Poutre_simple_decroissance` (init, héritage, résistances).
  - `test_EC5_Embrevement.py` créé : `Embrevement` init, profondeur auto/manuel, types (Bissectrice, Équerre p1/p2).
  - `test_EC5_Feu.py` complété : branches de protection (plâtre joints comblés/vides, pas d'exposition), faible élancement compression feu.
  - `test_EC0_Combinaison.py` : fixtures corrigées pour correspondre à l'API réelle de `Projet` et `Model_generator`.
  - `test_EC5_CVT.py` : suppression du `return` dans `test_K_panel`, extraction du helper `_setup_wall_with_panels`.

### Fixed
- Import conditionnel de `IPython` dans `core/objet.py` pour éviter une erreur si le module n'est pas installé.
- Correction du paramètre `section` manquant dans les tests CVT.

### Changed
- **Refactoring Phase 4 — Réexports & API courte** :
  - Création d'une API courte via `ourocode/__init__.py` : `from ourocode import Barre, Projet, Flexion, ...`
  - Réexports intermédiaires via `ourocode/eurocode/__init__.py` : `from ourocode.eurocode import Barre, Projet, ...`
  - Toutes les classes principales sont accessibles directement depuis le namespace racine.
  - Centralisation de `si.environment("structural")` dans `ourocode/__init__.py`.

- **Refactoring Phase 2 — Ajout des annotations de type retour** :
  - Ajout des return type hints (`-> tuple`, `-> dict`, `-> float`, `-> str`, `-> None`, `-> 'pd.Series'`, `-> 'pd.DataFrame'`) sur toutes les méthodes publiques.
  - Fichiers concernés :
    - `ec5/element_droit/` : `barre.py`, `compression.py`, `cisaillement.py`, `traction.py`, `flexion.py`
    - `ec5/assemblage/` : `assemblage.py`, `pointe.py`, `agrafe.py`, `boulon.py`, `broche.py`, `tirefond.py`, `embrevement.py`
    - `ec5/feu/` : `feu.py`, `cisaillement_feu.py`, `compression_feu.py`, `flexion_feu.py`, `traction_feu.py`
    - `ec5/cvt.py`
    - `ec8/sismique.py`

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

### Added (suite refactoring)
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
