# coding in UTF-8
# by Anthony PARISOT
import forallpeople as si

si.environment("structural")

from ourocode.eurocode.ec3.element_droit.plat import Plat


class Profil(Plat):
    """Classe pour les profilés laminés acier selon EN 1993-1-1.

    Permet de définir un profilé IPE, IPN, HEA, HEB, HEM ou UPN à partir de sa
    désignation normalisée (ex : "IPE 300"). Les caractéristiques géométriques sont
    chargées depuis `data/profil_acier.csv` (valeurs en mm, mm², mm³, mm⁴) et les
    caractéristiques mécaniques (fy, fu) sont déterminées selon la classe d'acier
    et l'épaisseur maximale max(tw, tf).

    La désignation est propagée par `_from_parent_class` : les classes de vérification
    (Compression, Traction, Flexion, Cisaillement) récupèrent automatiquement A, Avz,
    Iy, Iz, Wel, Wpl, etc.

    Attributs (avec unités) : h, b, tw, tf, r, A, Avz, Avy, Iy, Iz, Wely, Welz,
    Wply, Wplz (None pour les UPN), iy, iz, masse.
    """

    def __init__(
        self,
        designation: str = Plat.DESIGNATIONS,
        classe_acier: str = Plat.CLASSE_STEEL,
        classe_transv: int = ("1", "2", "3"),
        **kwargs,
    ):
        """Initialise un profilé laminé acier.

        Args:
            designation (str): Désignation du profil (ex: "IPE 300", "HEA 200", "UPN 200").
            classe_acier (str, optional): Classe d'acier selon EN 1993-1-1 Tableau 3.1
                (ex: "S235", "S275", "S355"). Defaults to "S235".
            classe_transv (int, optional): Classe transversale de la section (1, 2 ou 3)
                selon EN 1993-1-1 §5.5. Defaults to 1.
            **kwargs: Arguments transmis à la classe parent Plat.

        Raises:
            ValueError: Si designation n'est pas dans DESIGNATIONS.
            ValueError: Si l'épaisseur maximale dépasse 80 mm (hors tableau 3.1 implémenté).
        """
        if isinstance(designation, (list, tuple)):
            designation = designation[0]
        super().__init__(classe_acier=classe_acier, classe_transv=classe_transv, designation=designation, **kwargs)

    @classmethod
    def liste(cls, type_profil: str = None) -> tuple:
        """Retourne les désignations disponibles, éventuellement filtrées par type (ex: "IPE")."""
        if type_profil is None:
            return cls.DESIGNATIONS
        df = cls._PROFIL_DATA
        return tuple(df.loc[df["type"] == type_profil, "designation"])
