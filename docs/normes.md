# Référence des normes Eurocodes

Cette page récapitule les normes européennes (Eurocodes) implémentées dans Ourocode, ainsi que les Annexes Nationales Françaises (AN) appliquées.

---

## EN 1990 — Bases de calcul des structures

**Objet :** Établit les principes et les exigences fondamentales pour la sécurité, la durabilité et la fonctionnalité des structures.

**Contenu principal :**

- Situations de calcul (persistante, transitoire, accidentelle, sismique)
- Combinaisons d'actions aux états limites ultimes (ELU) et de service (ELS)
- Coefficients partiels γ et coefficients de combinaison ψ
- Annexe A1 : Bâtiments

**Module Ourocode :** [`core.combinaison`](api/EC0_Combinaison.md)

---

## EN 1991 — Actions sur les structures

### EN 1991-1-1 : Charges permanentes et charges d'exploitation

**Objet :** Densités volumiques des matériaux, charges de planchers, toitures, balcons, etc.

**Module Ourocode :** [`ec1.exploitation`](api/EC1_Exploitation.md)

---

### EN 1991-1-3 : Charges de neige

**Objet :** Détermination des charges de neige sur les toitures selon la zone géographique, l'altitude et la forme de toiture.

**Principaux paramètres :**

- Valeur caractéristique `sk` selon la zone de neige (AN France : zones A1 à E)
- Coefficient de forme `μ1, μ2, μ3` selon la pente du toit et son type
- Coefficient d'exposition `Ce` et coefficient thermique `Ct`

**Module Ourocode :** [`ec1.neige`](api/EC1_Neige.md)

---

### EN 1991-1-4 : Actions du vent

**Objet :** Calcul des pressions et forces dues au vent sur les structures et éléments de structure.

**Principaux paramètres :**

- Vitesse de référence du vent `vb,0` selon la zone de vent (AN France : zones 1 à 4)
- Profil de vitesse selon la catégorie de terrain (0 à IV)
- Coefficients de pression extérieure `Cpe` selon la géométrie (murs, toitures 1 et 2 versants, toiture terrasse)
- Pression dynamique de pointe `qp(z)`

**Module Ourocode :** [`ec1.vent`](api/EC1_Vent.md)

---

## EN 1993 — Calcul des structures en acier

### EN 1993-1-1 : Règles générales

**Objet :** Vérification des sections transversales et des barres en acier sous sollicitations composées.

**Vérifications couvertes :**

- Traction, compression, flexion, cisaillement, torsion
- Flambement par flexion (courbes a0, a, b, c, d)
- Déversement

**Module Ourocode :** [`ec3.element_droit`](api/EC3_Element_droit.md)

---

### EN 1993-1-2 : Calcul en situation d'incendie

**Objet :** Comportement des structures acier en situation d'incendie, évolution de la température.

**Module Ourocode :** [`ec3.feu`](api/EC3_Feu.md)

---

### EN 1993-1-8 : Calcul des assemblages

**Objet :** Assemblages boulonnés et soudés en acier.

**Module Ourocode :** [`ec3.assemblage`](api/EC3_Assemblage.md)

---

## EN 1995 — Calcul des structures en bois

### EN 1995-1-1 : Règles générales

**Objet :** Vérification des éléments droits en bois massif, lamellé-collé et dérivés du bois.

**Vérifications couvertes :**

- Traction parallèle au fil
- Compression parallèle et perpendiculaire au fil
- Flexion simple et déviée
- Cisaillement
- Flambement (`kc,y`, `kc,z`)
- Déversement
- Flèche (classe de déformation, facteur `kdef`)
- Vérification du contreventement des MOBs.

**Classes de résistance bois massif :** C14, C16, C18, C20, C22, C24, C27, C30, C35, C40  
**Classes de résistance lamellé-collé :** GL20h, GL22h, GL24h, GL26h, GL28h, GL30h, GL32h

**Modules Ourocode :** [`ec5.element_droit`](api/EC5_Element_droit.md), [`ec5.blc`](api/EC5_BLC.md), [`ec5.cvt`](api/EC5_CVT.md)

---

### EN 1995-1-2 : Calcul en situation d'incendie

**Objet :** Comportement des structures bois en situation d'incendie, vitesse de carbonisation.

**Paramètres clés :**

- Vitesse de carbonisation nominale `β0` et de calcul `βn` selon l'essence
- Section résiduelle efficace après carbonisation

**Module Ourocode :** [`ec5.feu`](api/EC5_Feu.md)

---

### EN 1995-1-1 §8 : Assemblages mécaniques

**Objet :** Assemblages avec organes mécaniques (boulons, broches, pointes, vis, connecteurs).

**Module Ourocode :** [`ec5.assemblage`](api/EC5_Assemblage.md)

---

## EN 1998 — Calcul des structures pour leur résistance aux séismes

### EN 1998-1 : Règles générales

**Objet :** Spectre de réponse élastique et de calcul selon la zone sismique et le type de sol.

**Paramètres :**

- Zone sismique (AN France : zones 1 à 5)
- Classe d'importance du bâtiment
- Type de spectre (type 1 ou type 2)

**Module Ourocode :** [`ec8.sismique`](api/EC8_Sismique.md)

---

## Annexes Nationales Françaises

Toutes les implémentations s'appuient sur les valeurs recommandées des **Annexes Nationales Françaises (NF EN)** publiées par l'AFNOR.  
Les paramètres propres à chaque AN (coefficients ψ, valeurs de référence du vent, zones de neige, zones sismiques, etc.) sont pré-intégrés dans les fichiers CSV du répertoire `ourocode/data/`.
