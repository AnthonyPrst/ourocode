# Ourocode — Contexte agent

## Vue d'ensemble du projet

**Ourocode** est une bibliothèque Python (v2.1.8, Python ≥ 3.12) de calcul de structure selon les Eurocodes avec Annexes Nationales Françaises.
Auteur : Anthony PARISOT — [OUREA STRUCTURE](https://ourea-structure.fr)
Licence : Apache 2.0

---

## Architecture du code

```
ourocode/
├── ourocode/
│   ├── __init__.py                  # Réexport API courte : from ourocode import Barre, Flexion...
│   ├── __version__.py
│   ├── data/                        # Données CSV/JSON de référence normative
│   │   ├── caracteristique_meca_acier.csv / _bois.csv / _panel.csv / _beton.csv
│   │   ├── carte_action_region.csv  # Zones neige/vent/sismique par code INSEE
│   │   ├── coeff_psy.csv, kmod.csv, kdef.csv, kfi.csv, gammaM.csv
│   │   ├── limite_fleche.csv, exploitation.csv
│   │   ├── section_boulon.csv, qualite_acier.csv, vis_bois.csv (ROTHOBLAAS)
│   │   ├── sismique/                # Données sismiques (catégories, sols…)
│   │   └── vent/                    # Données vent par zones
│   └── eurocode/
│       ├── core/
│       │   ├── objet.py             # Classe de base Objet
│       │   ├── projet.py            # Projet
│       │   ├── batiment.py          # Batiment
│       │   ├── combinaison.py       # EN 1990 – Combinaisons d'actions
│       │   ├── model_generator.py   # Modèle MEF via PyNiteFEA (+ Wood_beam_model)
│       │   ├── model_result.py      # Exploitation des résultats MEF
│       │   └── _renderer.py         # Décorateur @handcalc interne (remplace handcalcs, compat 3.13+)
│       ├── ec1/
│       │   ├── exploitation.py      # EN 1991-1-1 – Charges d'exploitation
│       │   ├── neige.py             # EN 1991-1-3 – Action neige
│       │   └── vent.py              # EN 1991-1-4 – Action vent
│       ├── ec3/
│       │   ├── element_droit/       # plat, profil, traction, compression, cisaillement, flexion
│       │   ├── assemblage/          # tige, soudure, platine
│       │   └── feu.py               # Vérification feu acier
│       ├── ec5/
│       │   ├── element_droit/       # barre, flexion, traction, compression, cisaillement, verification_EC5
│       │   ├── assemblage/          # assemblage, pointe, agrafe, boulon, broche, tirefond, embrevement
│       │   ├── feu/                 # feu, flexion_feu, traction_feu, compression_feu, cisaillement_feu
│       │   ├── blc.py               # Poutre_simple_decroissance (lamellé-collé)
│       │   └── cvt.py               # Murs à ossature bois (MOB/CVT)
│       ├── ec8/
│       │   └── sismique.py          # EN 1998 – Action sismique
│       └── mixins/
│           ├── data_loader.py       # DataLoaderMixin : CSV/JSON normatifs + caches
│           ├── math_utils.py        # MathUtilsMixin : conversions d'unités, interpolation
│           ├── serialization.py     # SerializationMixin : JSON/CSV + pickle .oco (+ compat v1.x)
│           └── synthese.py          # SyntheseMixin : tableau des taux de travail
└── tests/
    ├── test_A0_Projet.py, test_model_generator.py, test_EC0_Combinaison.py
    ├── test_EC1_Exploitation.py, test_EC1_Neige.py
    ├── test_EC3_Element_droit.py, test_EC3_Assemblage.py, test_EC3_Feu.py, test_EC3_Platine.py
    ├── test_EC5_Element_droit.py, test_EC5_Assemblage.py, test_EC5_BLC.py
    ├── test_EC5_CVT.py, test_EC5_Embrevement.py, test_EC5_Feu.py, test_EC5_Verification.py
    └── test_EC8_Sismique.py
```

---

## Hiérarchie des classes (héritage)

```
Objet  (core/objet.py) = SerializationMixin + DataLoaderMixin + MathUtilsMixin + SyntheseMixin
└── Projet  (core/projet.py)
    ├── Batiment  (core/batiment.py)
    │   ├── Neige   (ec1/neige.py)
    │   ├── Vent → Murs_verticaux, Toiture_*  (ec1/vent.py)
    │   ├── Sismique (ec8/sismique.py)
    │   └── MOB     (ec5/cvt.py)
    ├── Combinaison  (core/combinaison.py)
    ├── Model_generator → Wood_beam_model  (core/model_generator.py)
    │   └── Model_result (core/model_result.py)
    ├── Barre  (ec5/element_droit/barre.py)
    │   ├── Flexion, Traction, Compression (+_perpendiculaire, +_inclinees),
    │   │   Cisaillement, Verification_EC5
    │   └── Feu → Flexion_feu, Traction_feu, Compression_feu, Cisaillement_feu
    ├── Plat / _Profil  (ec3/element_droit/)  ← acier
    │   └── Traction, Compression, Cisaillement, Flexion (EC3), Feu_acier
    └── Assemblage  (ec5/assemblage/assemblage.py)
        └── Pointe, Agrafe, Boulon, Tirefond_inf_7, Tirefond_sup_6,
            Vis_catalogue, Broche, Embrevement
```

---

## Classe de base : `Objet`

Toutes les classes héritent de `Objet` (`core/objet.py`), qui compose 4 mixins :

- **`SerializationMixin`** : `_physical_to_dict`/`_dict_to_physical`, `_resolve_unit_expr` (parsing d'unités sans eval), `save_data`/`load_data` (JSON/CSV), `save_object`/`_open_object`/`_open_multiple_objects` (pickle `.oco`).
  - Toutes ces méthodes acceptent un `path` optionnel → headless ; sinon `QFileDialog` (extra `[gui]`).
  - `_CompatUnpickler` redirige les classes des fichiers `.oco` créés avec l'arborescence v1.x (`EC5_Element_droit`, `A0_Projet`, …) vers les nouveaux modules.
- **`DataLoaderMixin`** : `PATH_CATALOG`, `_data_from_csv` (cache), `_data_from_json`, `_load_json`.
- **`MathUtilsMixin`** : `_convert_unit_physical`, interpolations, `_PHYSICAL_UNITS` / `_SET_VALUE_UNITS`.
- **`SyntheseMixin`** : `synthese_taux_travail`, `_add_synthese_taux_travail`.

Sur `Objet` directement :
- **`_from_parent_class(objet, **kwargs)`** : instanciation d'une sous-classe depuis un parent (pattern central).
- **`_from_dict(dictionary)`** : instanciation depuis un dict.
- **`_assign_handcalcs_value(handcalc_value, args)`** : assigne les résultats `@handcalc` aux attributs.
- **`set_value` / `get_value`** : encapsulation/extraction de valeurs (API GUI).
- **`JUPYTER_DISPLAY`** : booléen de classe, contrôle l'affichage LaTeX dans Jupyter.

---

## Unités : `forallpeople` (si)

Toutes les grandeurs physiques utilisent `forallpeople` avec l'environnement `"structural"` (initialisé dans `ourocode/__init__.py`).

| Entrée utilisateur | Unité interne |
|--------------------|---------------|
| dimensions (b, h, t) | `si.mm` |
| longueurs de projet | `si.m` |
| forces | `si.kN` ou `si.N` |
| contraintes | `si.MPa` ou `si.N/si.mm**2` |

Toujours fournir les valeurs numériques brutes au constructeur (ex : `b=100` pour 100 mm) — la conversion `* si.mm` est faite **dans** le constructeur.

---

## Génération LaTeX avec `@handcalc`

Le décorateur `@handcalc` est une réimplémentation **interne** (`core/_renderer.py`) de handcalcs, compatible Python 3.13+. Les méthodes de vérification retournent `(latex_str, valeur)`.

```python
from IPython.display import display, Latex
latex_taux, taux = panne_flexion.taux_m_d()
display(Latex(latex_taux))
```

---

## Logging

Le package utilise `logging` (`logger = logging.getLogger(__name__)` par module, `NullHandler` sur le logger `"ourocode"`). Ne pas utiliser `print()` dans le code de la bibliothèque (sauf blocs `if __name__ == "__main__":` de dev et méthodes explicitement diagnostiques comme `debug_node_compression`).

---

## Dépendances principales

| Package | Usage |
|---------|-------|
| `forallpeople` | Calcul avec unités physiques SI |
| `pandas` | Lecture des CSV de référence |
| `Pillow` | Affichage des images de schémas |
| `matplotlib` | Graphiques (ex : températures feu) |
| `PyNiteFEA` | Modèle éléments finis 3D |
| `pyvista[all,trame]` | Visualisation 3D |
| `PySide6` | Dialogues fichiers Qt (extra `[gui]`) |

Note : `handcalcs` n'est **plus** une dépendance — remplacé par `core/_renderer.py`.

Extras : `[gui]` (PySide6), `[test]` (pytest, pytest-cov), `[docs]` (mkdocs), `[dev]` (black).

---

## Tests

```bash
pytest tests/ -v --tb=short
pytest --cov=ourocode --cov-report=term-missing   # couverture mesurée : ~66 %
```

État connu : quelques tests échouent sur la branche de dev en cours (profilés EC3 catalogue incomplet, attentes `Vis_catalogue`, format tuple latex `Verification_EC5`) — travaux en cours, non liés à l'infrastructure.

---

## CI/CD (GitHub Actions — `.github/workflows/ci.yml`)

| Job | Dépend de | Condition |
|-----|-----------|----------|
| `lint` | — | Toujours (`ruff check ourocode/ tests/`) |
| `test` | `lint` | Matrice Python **3.12 / 3.13 / 3.14**, `pytest --cov-fail-under=60` |
| `docs` | `lint` | Push sur `main` uniquement → GitHub Pages |

Autres workflows : `pypi.yml` (publication), `release.yml`.

URL de la doc déployée : `https://anthonyprst.github.io/ourocode/`

---

## Documentation (`docs/` + `mkdocs.yml`)

MkDocs-Material + mkdocstrings. `mkdocs serve` pour prévisualiser.

---

## Conventions et patterns importants

0. **CHANGELOG.md à maintenir** : toute modification notable doit être documentée dans `CHANGELOG.md` sous `[Unreleased]` ou `[X.Y.Z]` (Keep a Changelog : `Added`, `Changed`, `Fixed`, `Removed`).

1. **Pattern `_from_parent_class`** :
   ```python
   barre = Barre(b=100, h=200, classe="C24", cs=2, ...)
   flexion = Flexion._from_parent_class(barre, lo_rel_y=5000, ...)
   ```
2. **`**kwargs` remontants** : tous les constructeurs transmettent `**kwargs` à `super().__init__()` → les paramètres `Projet` peuvent être passés directement à `Barre`.
3. **Arguments à choix multiples** : un tuple/liste dans la signature = énumération des valeurs acceptées (ex : `section=["Rectangulaire","Circulaire"]`). C'est une convention documentaire, pas une validation.
4. **Unités en entrée** : valeurs numériques brutes, conversion dans le constructeur (`self.b = b * si.mm`).
5. **Résultats des méthodes de vérification** : tuple `(latex_str, valeur)`.
6. **Persistance `.oco`** : format pickle (risque sécurité si fichier non fiable) ; `_CompatUnpickler` assure la compat v1.x. Toujours passer `path` pour un usage headless.
7. **Convention d'axes MEF ↔ EC5** : les axes PyNiteFEA et EC5 diffèrent (voir CHANGELOG 2.1.8 — inversions historiques). Toute modification dans `model_result.py` ou `verification_EC5.py` doit être validée par les tests de signe existants.
