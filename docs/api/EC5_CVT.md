# EC5_CVT — EN 1995-1-1 §9

Calcul des éléments de contreventement des MOB (Murs Ossature Bois) selon EN 1995-1-1 §9.

::: ourocode.eurocode.ec5.cvt
    options:
      members: true
      show_root_heading: true
      show_source: false

## Import depuis Cadwork 3D

1. Copier `ourocode/eurocode/ec5/cadwork/ourocode_export_mob.py` dans
   `<userprofil>/3d/API.x64/ourocode_export_mob/ourocode_export_mob.py`.
2. Dans Cadwork, sélectionner les parois ossature, panneaux et ouvertures à exporter
   (sans sélection : tous les éléments visibles), puis lancer le plugin `ourocode_export_mob`.
3. Dans ourocode :

```python
mob = MOB._from_parent_class(batiment, connecteurs_int=pointe_int, connecteurs_ext=pointe_ext)
mob.load_cadwork_data("export_cadwork.json", couturage_ext_rive=150, couturage_ext_inter=300)
mob.save_walls_data("murs.json")  # rechargeable ensuite avec load_walls_data
```

Les options du plugin (`DEFAULT_OPTIONS` : mots-clés des panneaux, épaisseur maximale, tolérance)
se modifient en tête du fichier. Les panneaux écartés et les avertissements sont listés dans
les clés `ignored` et `warnings` du JSON.

::: ourocode.eurocode.ec5.cadwork.ourocode_export_mob
    options:
      members: true
      show_root_heading: true
      show_source: false
