"""
Ourocode Eurocode modules - Niveau intermédiaire
Permet : from ourocode.eurocode import Barre, Projet, ...
"""

# Sub-packages
from ourocode.eurocode import core
from ourocode.eurocode import ec1
from ourocode.eurocode import ec3
from ourocode.eurocode import ec5
from ourocode.eurocode import ec8

# Core classes
from ourocode.eurocode.core import (
    Objet,
    Projet,
    Batiment,
    Model_generator,
    Wood_beam_model,
    Model_result,
    Combinaison,
)

# EC1 - Actions
from ourocode.eurocode.ec1 import Neige, Vent, Exploitation

# EC3 - Acier
from ourocode.eurocode.ec3 import (
    Plat, Traction as Traction_EC3, Compression as Compression_EC3,
    Cisaillement as Cisaillement_EC3, Flexion as Flexion_EC3,
    Tige, Tige_scellement_plaque_ancrage, Tige_scellement_crochet,
    Soudure,
    Platine_assise_compression_beton, Platine_assise_compression_bois, Platine_assise_traction,
)

# EC5 - Bois (classes principales)
from ourocode.eurocode.ec5 import (
    Barre, Flexion, Traction, Compression, Compression_perpendiculaire,
    Compression_inclinees, Cisaillement,
    Assemblage, Embrevement, Pointe, Agrafe, Boulon, Broche,
    Tirefond_inf_7, Tirefond_sup_6,
    Feu, Flexion_feu, Traction_feu, Compression_feu, Cisaillement_feu,
    Poutre_simple_decroissance, MOB,
)

# EC8 - Sismique
from ourocode.eurocode.ec8 import Sismique

__all__ = [
    # Sub-packages
    "core", "ec1", "ec3", "ec5", "ec8",
    # Core
    "Objet", "Projet", "Batiment", "Model_generator", "Wood_beam_model",
    "Model_result", "Combinaison",
    # EC1
    "Neige", "Vent", "Exploitation",
    # EC3
    "Plat", "Traction_EC3", "Compression_EC3", "Cisaillement_EC3", "Flexion_EC3",
    "Tige", "Tige_scellement_plaque_ancrage", "Tige_scellement_crochet",
    "Soudure",
    "Platine_assise_compression_beton", "Platine_assise_compression_bois", "Platine_assise_traction",
    # EC5
    "Barre", "Flexion", "Traction", "Compression", "Compression_perpendiculaire",
    "Compression_inclinees", "Cisaillement",
    "Assemblage", "Embrevement", "Pointe", "Agrafe", "Boulon", "Broche",
    "Tirefond_inf_7", "Tirefond_sup_6",
    "Feu", "Flexion_feu", "Traction_feu", "Compression_feu", "Cisaillement_feu",
    "Poutre_simple_decroissance", "MOB",
    # EC8
    "Sismique",
]
