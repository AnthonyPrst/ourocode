# Ourocode

**Bibliothèque Python pour le calcul de structure selon les Eurocodes (AN Françaises)**

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/AnthonyPrst/ourocode/blob/main/LICENSE)
[![CI](https://github.com/AnthonyPrst/ourocode/actions/workflows/ci.yml/badge.svg)](https://github.com/AnthonyPrst/ourocode/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/AnthonyPrst/ourocode)](https://github.com/AnthonyPrst/ourocode/releases)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://anthonyprst.github.io/ourocode/)

---

## Présentation

**Ourocode** est une bibliothèque Python regroupant les formules normalisées issues des différentes parties des Eurocodes aux Annexes Nationales Françaises :

| Norme | Domaine |
|-------|---------|
| **EN 1990** | Bases de calcul des structures |
| **EN 1991** | Actions sur les structures (neige, vent, exploitation) |
| **EN 1993** | Calcul des structures en acier |
| **EN 1995** | Calcul des structures en bois |
| **EN 1998** | Calcul des actions sismiques |

Elle permet de construire un **catalogue de fonctions réutilisables** pour générer rapidement des **notes de calcul personnalisées**, intégrables dans des applications de vérification ou de génération de rapports.

---

## Installation rapide

```bash
pip install ourocode
```

Avec les dépendances optionnelles :

```bash
pip install "ourocode[full]"    # GUI (PySide6) + MEF (PyNiteFEA + pyvista)
pip install "ourocode[gui]"     # PySide6 uniquement
```

---

## Exemple minimal

```python
from IPython.display import display, Latex
from ourocode.eurocode.ec5.element_droit import Barre, Flexion

panne = Barre(b=100, h=200, section="Rectangulaire", classe="C24", cs=2, Hi=12, Hf=12)
panne_flexion = Flexion._from_parent_class(panne, lo_rel_y=5000, lo_rel_z=5000, coeflef=0.9, pos="Charge sur fibre comprimée")
latex_taux, taux = panne_flexion.taux_m_d()
display(Latex(latex_taux))
```

---

## Navigation

- [Guide de démarrage](guide.md) — Premiers pas et exemples d'utilisation
- [Référence API](api/objet.md) — Documentation complète des classes et méthodes
- [Référence des normes](normes.md) — EN 1990, EN 1991, EN 1993, EN 1995, EN 1998
- [Changelog](changelog.md) — Historique des versions

---

## Auteur

Développé par **Anthony PARISOT**, ingénieur structure bois & développeur, dans le cadre du projet open source [OUREA STRUCTURE](https://ourea-structure.fr).
