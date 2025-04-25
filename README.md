
# 📐 Eurocode - Bibliothèque Python pour le calcul structurel selon les Eurocodes

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/eurocode-calc.svg)](https://pypi.org/project/eurocode-calc/)
[![Tests](https://img.shields.io/github/actions/workflow/status/ton-org/eurocode/python-app.yml?branch=main)](https://github.com/ton-org/eurocode/actions)
[![Coverage](https://img.shields.io/codecov/c/github/ton-org/eurocode)](https://codecov.io/gh/ton-org/eurocode)

---

## 🔍 Description

**Eurocode** est une bibliothèque Python regroupant les formules normalisées issues des différentes parties des Eurocodes :

- **EN 1990** : Bases de calcul
- **EN 1991** : Actions sur les structures
- **EN 1993** : Calcul des structures en acier
- **EN 1995** : Calcul des structures en bois

Elle permet de construire un **catalogue de fonctions réutilisables** pour générer rapidement des **notes de calcul personnalisées**, intégrables dans des applications de vérification ou de génération de rapports.

---

## 🚀 Objectifs

- Offrir une **implémentation fiable et vérifiée** des formules Eurocode.
- Permettre une **utilisation modulaire** dans des interfaces No-Code, des scripts ou des applications.
- Fournir une **base open source transparente**, adaptée aux ingénieurs, bureaux d’études, enseignants ou développeurs.

---

## 📦 Installation

```bash
pip install eurocode-calc
```

> Ou installez directement depuis le dépôt :
```bash
pip install git+https://github.com/ton-org/eurocode.git
```

---

## ✨ Fonctionnalités

- Calculs normalisés : flexion, traction, cisaillement, flambement, stabilité, flèche…
- Support des classes de résistance bois (C24, GL24h, etc.) et acier (S235, S355…)
- Prise en compte des effets de feu (EN 1995-1-2)
- Intégration avec `handcalcs` pour génération LaTeX des formules
- Compatible `PySide6` pour des interfaces graphiques
- Organisé en modules clairs par norme (EN1990, EN1991, etc.)

---

## 🛠 Exemple d'utilisation

```python
from eurocode.EC5_Element_droit import Barre

panne = Barre(b=100, h=200, section="Rectangulaire", classe="C24", cs=2, Hi=12, Hf=12)
latex, fmd = panne._f_type_d("fm0k", "Moyen terme", "Fondamentales")

print(f"Résistance de flexion : {fmd:.2f} MPa")
```

---

## 📚 Documentation

Une documentation détaillée avec exemples d’utilisation, formules LaTeX et visualisations est en cours de rédaction.

📘 Accès (prochainement) : [eurocode.readthedocs.io](https://eurocode.readthedocs.io)

---

## ✅ Tests & couverture

```bash
pytest --cov=eurocode --cov-report=html
```

Les tests couvrent les modules principaux : `Barre`, `Flexion`, `Cisaillement`, `Feu`, `Traction`, etc.

---

## 🤝 Contribuer

Les contributions sont les bienvenues ! Pour proposer une amélioration ou corriger un bug :

1. Fork le dépôt
2. Crée une branche (`git checkout -b feature/ta-fonction`)
3. Commits (`git commit -am "feat: ajout nouvelle vérif"`),
4. Pull request 📥

---

## 📄 Licence

Distribué sous licence **MIT** – libre d’usage, même commercial, avec attribution.

---

## 👷 Auteur

Développé par **Anthony PARISOT**, ingénieur structure bois & développeur, dans le cadre du projet open source [OUREA STRUCTURE](https://ourea-structure.fr).

---

## ⭐ Si vous trouvez ce projet utile...

N'hésitez pas à [⭐️ le repo GitHub](https://github.com/ton-org/eurocode) pour le soutenir !
