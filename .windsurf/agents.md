# Ourocode — Contexte agent Cascade

## Vue d'ensemble du projet

**Ourocode** est une bibliothèque Python (v1.10.0, Python ≥ 3.12) de calcul de structure selon les Eurocodes avec Annexes Nationales Françaises.  
Auteur : Anthony PARISOT — [OUREA STRUCTURE](https://ourea-structure.fr)  
Licence : Apache 2.0

---

## Architecture du code

```
ourocode/
├── ourocode/
│   ├── __version__.py
│   ├── data/                          # Données CSV/JSON de référence normative
│   │   ├── caracteristique_meca_acier.csv
│   │   ├── caracteristique_meca_bois.csv
│   │   ├── caracteristique_meca_panel.csv
│   │   ├── caracteristique_meca_beton.csv
│   │   ├── carte_action_region.csv    # Zones neige/vent/sismique par code INSEE
│   │   ├── coeff_psy.csv
│   │   ├── kmod.csv / kdef.csv / kfi.csv / gammaM.csv
│   │   ├── limite_fleche.csv
│   │   ├── exploitation.csv
│   │   ├── section_boulon.csv
│   │   ├── qualite_acier.csv
│   │   ├── sismique/                  # Données sismiques (catégories, sols…)
│   │   └── vent/                      # Données vent par zones
│   └── eurocode/
│       ├── objet.py                   # Classe de base Objet
│       ├── A0_Projet.py               # Classes Projet, Batiment, Model_generator
│       ├── EC0_Combinaison.py         # EN 1990 – Combinaisons d'actions
│       ├── EC1_Exploitation.py        # EN 1991 – Charges d'exploitation
│       ├── EC1_Neige.py               # EN 1991 – Action neige
│       ├── EC1_Vent.py                # EN 1991 – Action vent
│       ├── EC3_Element_droit.py       # EN 1993 – Éléments droits acier
│       ├── EC3_Assemblage.py          # EN 1993 – Assemblages acier
│       ├── EC3_Feu.py                 # EN 1993 – Vérification feu acier
│       ├── EC5_Element_droit.py       # EN 1995 – Éléments droits bois
│       ├── EC5_BLC.py                 # EN 1995 – Paramètres bois lamellé-collé
│       ├── EC5_Assemblage.py          # EN 1995 – Assemblages bois
│       ├── EC5_CVT.py                 # EN 1995 – Murs à ossature bois (MOB/CVT)
│       ├── EC5_Feu.py                 # EN 1995 – Vérification feu bois
│       └── EC8_Sismique.py            # EN 1998 – Action sismique
└── tests/
    ├── test_A0_Projet.py
    ├── test_EC0_Combinaison.py
    ├── test_EC1_Neige.py
    ├── test_EC3_Assemblage.py
    ├── test_EC3_Element_droit.py
    ├── test_EC5_Assemblage.py
    ├── test_EC5_CVT.py
    ├── test_EC5_Element_droit.py
    ├── test_EC5_Feu.py
    └── test_EC8_Sismique.py
```

---

## Hiérarchie des classes (héritage)

```
Objet  (objet.py)
└── Projet  (A0_Projet.py)
    ├── Batiment  (A0_Projet.py)
    │   ├── Neige   (EC1_Neige.py)
    │   ├── Sismique (EC8_Sismique.py)
    │   └── MOB     (EC5_CVT.py)
    ├── Combinaison  (EC0_Combinaison.py)
    ├── Barre  (EC5_Element_droit.py)
    │   ├── Flexion, Traction, Compression, Cisaillement,
    │   │   Compression_perpendiculaire, Compression_inclinees, Deformation…
    │   └── Feu  (EC5_Feu.py)
    ├── Plat  (EC3_Element_droit.py)  ← acier
    └── Assemblage  (EC5_Assemblage.py)
        └── Pointe, Agrafe, Boulon, Tirefond, Broche…
```

---

## Classe de base : `Objet`

Toutes les classes héritent de `Objet` (`objet.py`). Points clés :

- **`PATH_CATALOG`** : résolution automatique du chemin vers `ourocode/data/`.
- **`_data_from_csv(data_file)`** : lecture des CSV de référence via pandas.
- **`_from_parent_class(objet, **kwargs)`** : instanciation d'une sous-classe à partir d'une instance parente existante (pattern central du catalogue).
- **`_from_dict(dictionary)`** : instanciation depuis un dictionnaire.
- **`_assign_handcalcs_value(handcalc_value, args)`** : assigne les résultats `handcalcs` aux attributs d'instance.
- **`_physical_to_dict` / `_dict_to_physical`** : sérialisation/désérialisation des objets `forallpeople.Physical`.
- **`save_data` / `load_data`** : persistence JSON ou CSV (ouvre un `QFileDialog` si pas de chemin fourni).
- **`save_object` / `_open_object`** : sérialisation pickle vers fichier `.oco`.
- **`synthese_taux_travail()` / `_add_synthese_taux_travail()`** : tableau pandas de synthèse des taux de travail (situation normale + incendie).
- **`JUPYTER_DISPLAY`** : booléen de classe, contrôle l'affichage LaTeX dans Jupyter.

---

## Unités : `forallpeople` (si)

Toutes les grandeurs physiques utilisent `forallpeople` avec l'environnement `"structural"`.

| Entrée utilisateur | Unité interne |
|--------------------|---------------|
| dimensions (b, h, t) | `si.mm` |
| longueurs de projet | `si.m` |
| forces | `si.kN` ou `si.N` |
| contraintes | `si.MPa` ou `si.N/si.mm**2` |

La méthode `_convert_unit_physical` gère les conversions entre unités SI.  
Toujours fournir les valeurs numériques brutes au constructeur (ex : `b=100` pour 100 mm) — la conversion `* si.mm` est faite **dans** le constructeur.

---

## Génération LaTeX avec `handcalcs`

Chaque méthode de vérification retourne un tuple `(latex_str, valeur)` via le décorateur `@handcalc`.  
Exemple d'usage dans Jupyter :

```python
from IPython.display import display, Latex
latex_taux, taux = panne_flexion.taux_m_d()
display(Latex(latex_taux))
```

---

## Classes et modules clés

### `A0_Projet.py` — Base projet
- **`Projet`** : définit ingénieur, numéro, adresse, code INSEE, altitude.
- **`Batiment`** : ajoute dimensions (h, d, b), angles de toiture.
- **`Model_generator`** : modèle éléments finis via `PyNiteFEA` (FEModel3D).

### `EC0_Combinaison.py` — EN 1990
- **`Combinaison`** : génère les combinaisons ELU_STR, ELU_STR_ACC, ELS_C, ELS_QP à partir d'un `Model_generator`.
- Actions : G, Q, Sn, W+, W-, Sx, Ae.

### `EC1_Neige.py` — EN 1991-1-3
- **`Neige`** → hérite de `Batiment`.
- Récupère la zone neige via code INSEE (`carte_action_region.csv`).

### `EC1_Vent.py` — EN 1991-1-4
- Classes vent héritant de `Batiment`.

### `EC3_Element_droit.py` — EN 1993-1-1
- **`Plat`** : plaque acier (t, h, b, classe_acier S235/S355…, classe transversale 1-3).
- Classe 4 non implémentée → `ValueError`.

### `EC3_Assemblage.py` — EN 1993 assemblages acier

### `EC3_Feu.py` — EN 1993 feu

### `EC5_Element_droit.py` — EN 1995-1-1
- **`Barre`** : pièce bois (b, h, section, classe C24/GL24h/OSB…, cs 1-3, humidité Hi/Hf).
- Sous-classes de vérification : `Flexion`, `Traction`, `Compression`, `Cisaillement`, `Compression_perpendiculaire`, `Compression_inclinees`, `Deformation`.

### `EC5_Assemblage.py` — EN 1995 assemblages bois
- **`Assemblage`** : bois/bois ou bois/métal, prend deux objets `Barre` ou `Plat`.
- Sous-classes : `Pointe`, `Agrafe`, `Boulon`, `Tirefond`, `Broche`.

### `EC5_CVT.py` — EN 1995 §9.2.4 Murs à ossature bois
- **`MOB`** → hérite de `Batiment`.
- Méthode A : systèmes de murs, panneaux, connecteurs, calcul efforts CVT + déformations.

### `EC5_Feu.py` — EN 1995-1-2
- **`Feu`** → hérite de `Barre`.
- Vitesse de carbonisation, section résiduelle, vérifications ELU incendie.

### `EC8_Sismique.py` — EN 1998
- **`Sismique`** → hérite de `Batiment`.
- Méthode des forces latérales, classes de ductilité (DCL/DCM/DCH), spectre de réponse.

---

## Données de référence (`ourocode/data/`)

| Fichier | Contenu |
|---------|---------|
| `caracteristique_meca_bois.csv` | fm, ft, fc, E… pour C14→GL32h, LVL, OSB |
| `caracteristique_meca_panel.csv` | Caractéristiques panneaux OSB/CP |
| `caracteristique_meca_acier.csv` | fy, fu pour S235, S275, S355, S420, S460 |
| `carte_action_region.csv` | Zone neige/vent/sismique par code INSEE (5 chiffres) |
| `kmod.csv` | Coefficients kmod par classe de service et durée de chargement |
| `kdef.csv` | Coefficients de fluage kdef |
| `gammaM.csv` | Coefficients partiels γM par type de bois |
| `limite_fleche.csv` | Limites de flèche par type d'élément |
| `sismique/` | Catégories d'importance, classes de sol |

---

## Dépendances principales

| Package | Usage | Optionnel |
|---------|-------|----------|
| `forallpeople` | Calcul avec unités physiques SI | Non |
| `handcalcs` | Génération de LaTeX depuis du Python | Non |
| `pandas` | Lecture des CSV de référence | Non |
| `Pillow` | Affichage des images de schémas | Non |
| `PyNiteFEA` | Modèle éléments finis 3D | Oui (`[mef]`) |
| `pyvista` | Visualisation 3D | Oui (`[mef]`) |
| `PySide6` | Dialogues fichiers Qt | Oui (`[gui]`) |
| `mkdocs-material` | Génération de la documentation | Oui (`[docs]`) |
| `mkdocstrings[python]` | Extraction des docstrings | Oui (`[docs]`) |
| `ruff` | Linting (CI) | Oui (`[dev]`) |
| `pytest` / `pytest-cov` | Tests unitaires | Oui (`[test]`) |

Installation des extras :

```bash
pip install "ourocode[full]"    # gui + mef
pip install "ourocode[docs]"   # documentation
pip install "ourocode[test]"   # tests
pip install "ourocode[dev]"    # linting
```

---

## Tests

```bash
pytest tests/ -v --tb=short
pytest --cov=ourocode --cov-report html
```

Fichiers de test dans `tests/`, nommés `test_<module>.py`.  
Couvrent les modules principaux. Pas de test dédié à `EC1_Vent.py` ni `EC3_Feu.py` pour l'instant.

---

## CI/CD (GitHub Actions — `.github/workflows/ci.yml`)

Trois jobs séquentiels déclenchés sur push/PR vers `main` :

| Job | Dépend de | Condition |
|-----|-----------|----------|
| `lint` | — | Toujours |
| `test` | `lint` | Toujours |
| `docs` | `test` | Push sur `main` uniquement |

- **`lint`** : `ruff check ourocode/ tests/`
- **`test`** : `pip install -e ".[test]"` + `pytest tests/ -v --tb=short`
- **`docs`** : `pip install -e ".[docs]"` + `mkdocs gh-deploy --force --clean` → déploie sur la branche `gh-pages`

URL de la doc déployée : `https://anthonyprst.github.io/ourocode/`

---

## Documentation (`docs/` + `mkdocs.yml`)

Outil : **MkDocs-Material** + **mkdocstrings[python]**.

Structure :

```
docs/
├── index.md          # Page d'accueil
├── guide.md          # Guide de démarrage
├── normes.md         # Référence EN 1990/1991/1993/1995/1998
├── changelog.md      # Inclusion automatique du CHANGELOG.md
├── assets/
│   └── logo.png
└── api/
    ├── objet.md
    ├── A0_Projet.md
    ├── EC0_Combinaison.md
    ├── EC1_Exploitation.md
    ├── EC1_Neige.md
    ├── EC1_Vent.md
    ├── EC3_Element_droit.md
    ├── EC3_Feu.md
    ├── EC3_Assemblage.md
    ├── EC5_Element_droit.md
    ├── EC5_BLC.md
    ├── EC5_Feu.md
    ├── EC5_Assemblage.md
    ├── EC5_CVT.md
    └── EC8_Sismique.md
```

Commandes locales :

```bash
pip install "ourocode[docs]"
mkdocs serve    # prévisualisation sur http://127.0.0.1:8000
mkdocs build    # génération statique dans site/
```

---

## Conventions et patterns importants

0. **CHANGELOG.md à maintenir** : toute modification notable (ajout de module, correction de bug, nouvelle fonctionnalité, changement de dépendance, mise à jour de doc ou CI) **doit être documentée dans `CHANGELOG.md`** sous la section `[Unreleased]` ou dans une nouvelle version `[X.Y.Z]`. Format [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) avec les sous-sections `Added`, `Changed`, `Fixed`, `Removed`.

1. **Pattern `_from_parent_class`** : pour enchaîner les vérifications, on passe un objet déjà instancié à une sous-classe :
   ```python
   barre = Barre(b=100, h=200, classe="C24", cs=2, ...)
   flexion = Flexion._from_parent_class(barre, lo_rel_y=5000, ...)
   ```
2. **`**kwargs` remontants** : tous les constructeurs acceptent `**kwargs` transmis à `super().__init__()`, ce qui permet de passer les paramètres `Projet` directement lors de la création d'une `Barre`.
3. **Arguments à choix multiples** : les paramètres dont la valeur est un tuple/liste dans la signature sont des énumérations de valeurs acceptées (ex : `section=["Rectangulaire","Circulaire"]`).
4. **Unités en entrée** : toujours des valeurs numériques brutes, la conversion est faite dans le constructeur (`self.b = b * si.mm`).
5. **Résultats des méthodes de vérification** : tuple `(latex_str, valeur_numérique)`.
6. **`QFileDialog`** est requis pour `save_data`/`load_data` sans argument `path` → nécessite un environnement Qt actif.
