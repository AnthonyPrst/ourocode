# coding in UTF-8
# by Anthony PARISOT
import pandas as pd

import forallpeople as si
si.environment("structural")

from ourocode.eurocode.ec3.element_droit.plat import Plat


class _Profil(Plat):
    """Classe pour les profilés laminés acier selon EN 1993-1-1.
    
    Permet de générer des poutres en I ou H (IPE, HEA, HEB, UPN, etc.)
    à partir de désignations normalisées selon EN 10365.
    
    Les caractéristiques géométriques sont chargées depuis le CSV
    `profil_acier.csv` et les propriétés mécaniques (fy, fu) sont
    déterminées selon la classe d'acier et l'épaisseur maximale.
    
    Hérite de Plat pour être compatible avec les classes de vérification
    (Compression, Traction, Flexion, Cisaillement) via _from_parent_class.
    """
    # Chargement des données de profils
    _PROFIL_DATA = Plat._data_from_csv(Plat, "profil_acier.csv", index_col=None)
    TYPES_PROFIL = tuple(_PROFIL_DATA["type"].unique())
    DESIGNATIONS = tuple(_PROFIL_DATA["designation"])
    
    def __init__(self, designation: str, classe_acier: str = "S235", 
                 classe_transv: int = 1, **kwargs):
        """Initialise un profilé laminé acier.
        
        Args:
            designation (str): Désignation du profil (ex: "IPE 300", "HEA 200", "HEB 260")
            classe_acier (str, optional): Classe d'acier selon EN 1993-1-1 Tableau 3.1
                (ex: "S235", "S275", "S355"). Defaults to "S235".
            classe_transv (int, optional): Classe transversale de la section (1, 2 ou 3)
                selon EN 1993-1-1 §5.5. Defaults to 1.
            **kwargs: Arguments transmis à la classe parent Plat.
                
        Raises:
            ValueError: Si designation n'est pas dans DESIGNATIONS.
        """
        # Validation de la désignation avant appel du parent
        if designation not in self.DESIGNATIONS:
            raise ValueError(
                f"Désignation '{designation}' invalide. Valeurs acceptées : {self.DESIGNATIONS}"
            )
        
        self.designation = designation
        
        # Chargement des données géométriques depuis le CSV
        self._load_profil_data()
        
        # Appel du parent Plat avec t=tw (épaisseur âme)
        # h et b ne sont pas passés car ils seraient écrasés par Plat
        super().__init__(t=self.tw.value*10**3, h=self.h.value*10**3, b=self.b.value*10**3,
                        classe_acier=classe_acier, classe_transv=classe_transv,
                        **kwargs)
        
        # Restauration des valeurs de h et b après l'appel du parent
        self.h = self._h_backup
        self.b = self._b_backup
        
        # Surcharge de fy et fu selon l'épaisseur maximale (max(tw, tf))
        self.__fy_fu()
    
    def _load_profil_data(self):
        """Charge les caractéristiques géométriques du profil depuis le CSV."""
        df = self._PROFIL_DATA
        profil = df[df["designation"] == self.designation].iloc[0]
        
        # Dimensions géométriques
        self.h = profil["h"] * si.mm
        self.b = profil["b"] * si.mm
        self.tw = profil["tw"] * si.mm  # épaisseur âme
        self.tf = profil["tf"] * si.mm  # épaisseur semelle
        self.r = profil["r"] * si.mm    # rayon congé
        
        # Sauvegarde pour restauration après appel super().__init__
        self._h_backup = self.h
        self._b_backup = self.b
        
        # Caractéristiques de section
        self.A = profil["A"] * si.mm**2
        self.Av = profil["Av"] * si.mm**2
        self.Avz = profil["Avz"] * si.mm**2 if pd.notna(profil["Avz"]) else None
        
        # Moments quadratiques
        self.Iy = profil["Iy"] * si.mm**4
        self.Iz = profil["Iz"] * si.mm**4
        
        # Modules de flexion
        self.Wely = profil["Wely"] * si.mm**3
        self.Welz = profil["Welz"] * si.mm**3
        self.Wply = profil["Wply"] * si.mm**3
        self.Wplz = profil["Wplz"] * si.mm**3
        
        # Rayons de giration
        self.iy = profil["iy"] * si.mm
        self.iz = profil["iz"] * si.mm
        
        # Masse linéique
        self.masse = profil["masse"] * si.kg / si.m
    
    def __fy_fu(self):
        """Définit fy et fu en fonction de la classe d'acier et de l'épaisseur maximale.
        
        Surcharge la méthode de Plat pour utiliser max(tw, tf) au lieu de t uniquement,
        car les profilés ont deux épaisseurs différentes (âme et semelle).
        """
        t_max = max(self.tw, self.tf)
        
        if t_max <= 40 * si.mm:
            self.fy = self._Plat__classe_acier.loc["t<= 40  fy"] * si.MPa
            self.fu = self._Plat__classe_acier.loc["t<= 40  fu"] * si.MPa
        elif 40 * si.mm < t_max <= 80 * si.mm:
            self.fy = self._Plat__classe_acier.loc["40<t<= 80  fy"] * si.MPa
            self.fu = self._Plat__classe_acier.loc["40<t<= 80  fu"] * si.MPa
    
    @property
    def type_profil(self):
        """Retourne le type de profil (IPE, HEA, HEB, etc.)."""
        df = self._PROFIL_DATA
        return df[df["designation"] == self.designation].iloc[0]["type"]
