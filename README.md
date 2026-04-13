
# 📐 Ourocode - Bibliothèque Python pour le calcul de structure selon les Eurocodes

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![CI](https://github.com/AnthonyPrst/ourocode/actions/workflows/ci.yml/badge.svg)](https://github.com/AnthonyPrst/ourocode/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/AnthonyPrst/ourocode)](https://github.com/AnthonyPrst/ourocode/releases)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://anthonyprst.github.io/ourocode/)

---

## 🔍 Description

**Ourocode** est une bibliothèque Python regroupant les formules normalisées issues des différentes parties des Eurocodes aux Annexes Nationales Française :

- **EN 1990** : Bases de calcul
- **EN 1991** : Actions sur les structures
- **EN 1993** : Calcul des structures en acier
- **EN 1995** : Calcul des structures en bois
- **EN 1998** : Calcul des actions sismiques

Elle permet de construire un **catalogue de fonctions réutilisables** pour générer rapidement des **notes de calcul personnalisées**, intégrables dans des applications de vérification ou de génération de rapports.

---

## 🚀 Objectifs

- Offrir une **implémentation fiable et vérifiée** des formules Eurocode.
- Permettre une **utilisation modulaire** dans des interfaces No-Code, des scripts ou des applications.
- Fournir une **base open source transparente**, adaptée aux ingénieurs, bureaux d’études, enseignants ou développeurs.

---

## 📦 Installation

> **Note :** Python 3.12 requis. Les dépendances tierces ne sont pas encore compatibles avec Python 3.13+.

```bash
pip install ourocode
```

> Installez directement depuis le dépôt :
```bash
pip install git+https://github.com/AnthonyPrst/ourocode.git
```

> Avec les dépendances optionnelles :
```bash
pip install "ourocode[full]"    # GUI (PySide6) + MEF (PyNiteFEA + pyvista)
pip install "ourocode[gui]"     # PySide6 uniquement
pip install "ourocode[mef]"     # PyNiteFEA + pyvista uniquement
pip install "ourocode[docs]"    # Génération de la documentation
pip install "ourocode[test]"    # Outils de test (pytest)
```

---

## ✨ Fonctionnalités

- Calculs normalisés : flexion, traction, cisaillement, flambement, assemblage EC5 et EC3, flèche…
- Support des classes de résistance bois massif (C14→C40), lamellé-collé (GL20h→GL32h) et acier (S235, S355…)
- Prise en compte des effets de feu (EN 1993-1-2, EN 1995-1-2)
- Calcul des actions sismiques et spectre de réponse (EN 1998-1)
- Vérification des murs à ossature bois — MOB/CVT (EN 1995-1-1 §9)
- Intégration avec `handcalcs` pour génération LaTeX des formules
- Compatible avec le logiciel [OUREA-NoCode](https://ourea-structure.fr/telechargement)
- Organisé en modules clairs par norme et fonction (EN 1990, EN 1991, EN 1993, EN 1995, EN 1998)

---

## 🛠 Exemple d'utilisation dans Jupyter

```python
from IPython.display import display, Latex
from ourocode.eurocode.ec5.element_droit import Barre, Flexion

panne = Barre(b=100, h=200, section="Rectangulaire", classe="C24", cs=2, Hi=12, Hf=12)
panne_flexion = Flexion._from_parent_class(panne, lo_rel_y=5000, lo_rel_z=5000, coeflef=0.9, pos="Charge sur fibre comprimée")
latex_fmd, fmd = panne_flexion.f_m_d("Moyen terme", "Fondamentales")
latex_sigmamd, sigmamd = panne_flexion.sigma_m_d(20, 0)
latex_taux, taux = panne_flexion.taux_m_d()
display(Latex(latex_taux))
```

## ✅ Tests & couverture

```bash
pytest tests/ -v --tb=short
pytest --cov=ourocode --cov-report html
```

Les tests couvrent les modules principaux.

---

## 📚 Documentation

La documentation complète est disponible sur **[GitHub Pages](https://anthonyprst.github.io/ourocode/)** :

- [Guide de démarrage](https://anthonyprst.github.io/ourocode/guide/) — Premiers pas et exemples
- [Référence API](https://anthonyprst.github.io/ourocode/api/objet/) — Classes et méthodes de chaque module
- [Référence des normes](https://anthonyprst.github.io/ourocode/normes/) — EN 1990, EN 1991, EN 1993, EN 1995, EN 1998

Générer la documentation localement :

```bash
pip install "ourocode[docs]"
mkdocs serve
```

---

## 🤝 Contribuer

Les contributions sont les bienvenues ! Pour proposer une amélioration ou corriger un bug n'hésité pas à fork le dépôt !

---

## 📄 Licence

Distribué sous licence **Apache License 2.0**.

---

## 👷 Auteur

Développé par **Anthony PARISOT**, ingénieur structure bois & développeur, dans le cadre du projet open source [OUREA STRUCTURE](https://ourea-structure.fr).

---

## ⭐ Si vous trouvez ce projet utile...

N'hésitez pas à [⭐️ le repo GitHub](https://github.com/AnthonyPrst/ourocode) pour le soutenir !
