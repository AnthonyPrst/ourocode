# Guide de démarrage

## Prérequis

- Python ≥ 3.12
- `pip install ourocode`

---

## Concepts clés

### Hiérarchie des classes

Toutes les classes du package héritent de la même classe de base :

```
Objet
└── Projet          (A0_Projet)
    ├── EC0 — Combinaisons de charges
    ├── EC1 — Actions (neige, vent, exploitation)
    ├── EC3 — Structures acier
    │   ├── Element_droit, Feu, Assemblage
    └── EC5 — Structures bois
        ├── Barre, Flexion, Traction, Compression, Cisaillement…
        └── Feu, Assemblage, CVT
```

### Pattern `_from_parent_class`

La plupart des classes de vérification s'instancient à partir d'un objet parent (ex. `Barre`) via la méthode de classe `_from_parent_class`. Cela permet de chaîner les vérifications en réutilisant les propriétés déjà calculées.

```python
barre = Barre(...)                              # objet géométrique + matériau
flexion = Flexion._from_parent_class(barre, ...) # hérite des propriétés de barre
```

### Unités physiques

Les valeurs sont manipulées avec la bibliothèque `forallpeople` (système SI). Les entrées en **mm** et **mm²** sont généralement attendues ; les résultats sont exprimés en unités normalisées (MPa, kN, etc.).

### Génération LaTeX

Chaque méthode de calcul renvoie un tuple `(latex_str, valeur)`. La chaîne LaTeX peut être affichée dans Jupyter via :

```python
from IPython.display import display, Latex
display(Latex(latex_str))
```

---

## Exemples par module

### EN 1995 — Flexion d'une panne bois (EC5_Element_droit)

```python
from ourocode.eurocode.EC5_Element_droit import Barre, Flexion

panne = Barre(
    b=100, h=200,
    section="Rectangulaire",
    classe="C24",
    cs=2,          # Classe de service 2
    Hi=12, Hf=12   # Hauteur du bâtiment (m)
)

panne_flexion = Flexion._from_parent_class(
    panne,
    lo_rel_y=5000, lo_rel_z=5000,
    coeflef=0.9,
    pos="Charge sur fibre comprimée"
)

latex_fmd, fmd = panne_flexion.f_m_d("Moyen terme", "Fondamentales")
latex_sigma, sigma = panne_flexion.sigma_m_d(20, 0)  # My=20 kN.m, Mz=0
latex_taux, taux = panne_flexion.taux_m_d()

display(Latex(latex_taux))
```

### EN 1995 — Assemblage boulonné (EC5_Assemblage)

```python
from ourocode.eurocode.EC5_Assemblage import Assemblage

# Consulter la référence API pour les paramètres complets
```

### EN 1993 — Vérification acier (EC3_Element_droit)

```python
from ourocode.eurocode.EC3_Element_droit import Plat

# Consulter la référence API pour les paramètres complets
```

### EN 1991 — Charges de neige (EC1_Neige)

```python
from ourocode.eurocode.EC1_Neige import Neige

# Consulter la référence API pour les paramètres complets
```

### EN 1990 — Combinaisons de charges (EC0_Combinaison)

```python
from ourocode.eurocode.EC0_Combinaison import Combinaison

# Consulter la référence API pour les paramètres complets
```

---

## Génération de la documentation localement

```bash
pip install "ourocode[docs]"
mkdocs serve
```

Ouvrir ensuite [http://127.0.0.1:8000](http://127.0.0.1:8000) dans votre navigateur.

---

## Exécution des tests

```bash
pip install "ourocode[test]"
pytest tests/ -v --tb=short
```

Avec rapport de couverture :

```bash
pytest --cov=ourocode --cov-report=html
```
