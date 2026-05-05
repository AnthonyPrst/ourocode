# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/),
et ce projet adhère au [Semantic Versioning](https://semver.org/lang/fr/).

## [Unreleased]

### Fixed
- **`detect_continuous_members`** : correction du bug qui empêchait la détection de la continuité d'un arbalétrier (ou de toute barre) lorsqu'une contrefiche (ou autre barre transversale) était fixée en nœud intermédiaire (jonction T/Y). La contrainte de degré 2 au nœud a été supprimée : la colinéarité seule suffit à distinguer les barres continues des branchements. La recherche du voisin parcourt désormais toutes les barres du nœud pour trouver celle qui passe le test de continuité, et non plus uniquement la première trouvée. Les tests de `TestDetectContinuousMembers` ont été mis à jour en conséquence (2 nouveaux cas de régression : jonction T perpendiculaire et jonction T avec relâchement).
- **`Compression` / `_determine_compression_params` / `verification_EC5`** : correction du calcul des longueurs et conditions d'appui de flambement (EC5 §6.3.2). Trois bugs corrigés simultanément : (1) `Compression.__init__` utilisait un seul `type_appuis`/`coef_lef` appliqué identiquement aux axes y et z — remplacé par `type_appuis_y`/`type_appuis_z` et `coef_lef_y`/`coef_lef_z` indépendants (`coef_lef` conservé pour rétro-compatibilité) ; (2) `_determine_type_appuis` utilisait `Noeuds[0/1]` pour les nœuds d'extrémité (bug topologique) et ignorait les appuis intermédiaires — remplacé par `_determine_compression_params` qui découpe la barre en sous-portées via `_compute_spans` et détermine le `type_appuis` à chaque sous-portée via le nouveau helper `_type_appuis_between_nodes` ; (3) `verification_EC5` passait `type_appuis_z` pour les deux axes à `Compression` — désormais `lo_flamb_y`/`lo_flamb_z` (flambement, stockés séparément de `lo_rel_y`/`lo_rel_z`) et `type_appuis_y`/`type_appuis_z` sont transmis correctement par axe. Rétro-compatibilité assurée via fallback `lo_flamb_y → lo_rel_y` si la clé est absente du Design.
- **`_determine_flexion_params`** : correction du calcul des coefficients `coeflef_y` / `coeflef_z` sur les barres structurales multi-membrures. Un seul `coeflef` global était calculé pour l'ensemble de la barre et appliqué identiquement aux deux axes, sans tenir compte de la découpe en sous-portées par les appuis latéraux ni des différences de chargement par sous-portée. La méthode calcule désormais `(lo_rel, coeflef)` **par axe et par sous-portée** : pour chaque sous-portée délimitée par les appuis latéraux d'un axe, `_coeflef_for_span` détermine le coeflef local en fonction des charges verticales (distribuées, ponctuelle centrale, ponctuelle en bout) effectivement présentes sur ce segment ; la paire `(lo_rel × coeflef)` maximale — c'est-à-dire la plus défavorable — est retenue indépendamment pour Y et Z. Le helper `_compute_spans` a été introduit (refactorisation de `_compute_lo_rel`) pour exposer la liste des segments entre appuis.
- **`_compute_lo_rel` / `_determine_flexion_params`** : correction du bug de parcours topologique sur les barres structurales multi-membrures. Les deux méthodes utilisaient `Noeuds[0]`/`Noeuds[1]` directement, supposant à tort que toutes les barres FEM d'une chaîne sont orientées dans le même sens. Les nœuds d'extrémité et intermédiaires pouvaient ainsi être incorrects (ex : dernier nœud pointant vers le milieu de la chaîne plutôt que vers l'extrémité libre). Correction via le nouveau helper **`_ordered_chain_nodes(member_ids)`** qui parcourt la topologie réelle de la chaîne (en déduisant le sens de chaque barre depuis le nœud précédent), produisant une liste ordonnée `[(node_id, abscisse_mm), …]` fiable quelle que soit l'orientation des barres FEM individuelles.

### Added
- **`ourocode/eurocode/ec5/verification.py`** : nouvelle **classe `Verification_EC5(Projet)`** de vérification EC5 d'un modèle MEF complet (barres structurales continues ou non).
  - Factory :meth:`Verification_EC5.from_combinaison(combinaison, model_result)` qui hérite via `_from_parent_class([model_generator, combinaison])` des attributs Projet (ingénieur, nom, code INSEE, altitude) et attache la combinaison + le résultat MEF.
  - :meth:`verify(name, elu_filter="ELU_ALL", type_bat, pos_charge, coeflef_y, coeflef_z, type_appuis_compression, n_points)` **boucle sur toutes les combinaisons ELU** : pour chaque combo, `kmod` est dérivé via :meth:`Combinaison.min_type_load` et `γM` via :meth:`Combinaison.type_combi`, puis les efforts gouvernants (Nx signé, Vy, My, Mz) sont extraits via :meth:`Model_result.get_internal_force` (ciblage par `combo_name`). Pour chaque vérification (Flexion, Cisaillement, Traction **ou** Compression selon le signe de Nx), le taux maximum rencontré et la combinaison gouvernante associée sont retenus. La flexion inclut l'interaction flexo-compression/flexo-traction via `Flexion.taux_m_d(compression=…, traction=…)`. La flèche ELS (W_inst(Q), W_net,fin) est calculée indépendamment via le max tag-based.
  - :meth:`synthese(**kwargs)` : agrège toutes les barres structurales en un unique DataFrame trié par taux décroissant, avec gestion des erreurs (section manuelle, classe inconnue) en DataFrame secondaire.
  - Helpers privés : `_verify_combo_elu`, `_verify_fleche`, `_efforts_for_combo`, `_deflection`, `_resolve_section`, `_resolve_classe_bois`, `_design_or_default`, `_lo_mm`, `_check_state`, `_taux_to_df`.
  - Exposée depuis `ourocode.eurocode.ec5` via `__init__.py`.
- **`ourocode/eurocode/core/combinaison.py`** : nouvelle méthode `Combinaison.type_combi(name_combi) -> "Fondamentales" | "Accidentelles"` utilisée par `Verification_EC5` pour résoudre dynamiquement le coefficient partiel `γM` selon l'EC5 §2.4.1 (EN 1990 §6.4.3.2/3.3). Lève `ValueError` pour une combinaison non-ELU.
- **`tests/test_EC5_Verification.py`** : 18 tests d'intégration couvrant :
  - `Combinaison.type_combi` (Fondamentales, Accidentelles, ValueError hors ELU).
  - Factory `Verification_EC5.from_combinaison` : instanciation, héritage des attributs Projet depuis le modèle MEF, état valide pour `verify`.
  - Poutre simple 1 travée (DataFrame non vide, Flexion+Cisaillement+Flèche présents, Traction/Compression absents sans charge axiale, ordres de grandeur, combo gouvernante `ELU_STR 1.35G + 1.5Q`).
  - Boucle par combinaison : `min_type_load` retourne les bonnes durées par combo, cohérence entre `ELU_ALL` et `ELU_STR` en l'absence d'accidentelles.
  - Poutre continue 2 travées (auto-groupage + synthèse, tri décroissant, présence de la combo gouvernante).
- **`ourocode/eurocode/core/model_generator.py`** : nouveau concept de **barre structurale** regroupant une ou plusieurs barres FEM continues, pour préparer l'itération et la vérification Eurocode à l'échelle d'une barre physique (solive continue, panne sur plusieurs appuis, etc.).
  - `detect_continuous_members(angle_tol_deg=1.0)` : détection automatique des chaînes colinéaires de barres FEM partageant matériau/section, sans rotule (`teta_y`/`teta_z`) aux extrémités communes et sans T-junction.
  - `group_members(name, member_ids, role, classe_bois, cs, Hi, Hf, effet_systeme, type_element_fleche, lo_rel_y, lo_rel_z, comment)` : regroupement manuel avec configuration de design EC5.
  - `auto_group_continuous_members(angle_tol_deg, role, name_prefix, only_continuous, **design_kwargs)` : détection + regroupement en une passe.
  - `get_structural_member`, `get_all_structural_members`, `del_structural_member`, `iter_structural_members` : accès et itération.
  - Ajout de `"structural_members"` dans `self._data`.
- **`tests/test_model_generator.py`** : 17 tests couvrant la détection (continuité, rupture par section/matériau/rotule/T-junction/angle, tolérance angulaire, orientation inversée des extrémités) et la gestion des barres structurales (regroupement, validations, itération, auto-groupage).

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
