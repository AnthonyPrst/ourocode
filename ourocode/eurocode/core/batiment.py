# coding in UTF-8
# by Anthony PARISOT
import forallpeople as si

from ourocode.eurocode.core.projet import Projet


class Batiment(Projet):
    """Classe définissant la géométrie d'un bâtiment pour les calculs de structure.

    Cette classe décrit les dimensions principales du bâtiment et les
    caractéristiques de sa toiture, nécessaires pour le calcul des
    charges climatiques (neige, vent) et sismiques.

    Attributs de classe:
        ETAGE (tuple): Liste des niveaux courants (RDC à Toiture).
    """

    ETAGE = ("RDC", "R+1", "R+2", "R+3", "R+4", "Toiture")

    def __init__(
        self,
        h_bat: si.m,
        d_bat: si.m,
        b_bat: si.m,
        alpha_toit: float,
        alpha_toit2: float = 0,
        *args,
        **kwargs,
    ):
        """Initialise les dimensions du bâtiment.

        Args:
            h_bat (si.m): Hauteur totale du bâtiment en mètres, mesurée depuis
                le soubassement rigide ou les fondations (référence sismique).
            d_bat (si.m): Largeur du bâtiment en mètres (dimension perpendiculaire
                au vent dominant pour le calcul du vent).
            b_bat (si.m): Longueur du bâtiment en mètres.
            alpha_toit (float): Pente du premier versant de toiture en degrés.
                0° pour un toit plat, valeur positive pour un versant.
            alpha_toit2 (float, optional): Pente du second versant pour les
                toits à deux pans (0° si toit à un seul versant ou plat).
                Defaults to 0.
            *args: Arguments positionnels transmis à la classe parent Projet.
            **kwargs: Arguments nommés transmis à la classe parent Projet
                (ingenieur, code_INSEE, alt, etc.).
        """
        super().__init__(*args, **kwargs)
        self.h_bat = h_bat * si.m
        self.d_bat = d_bat * si.m
        self.b_bat = b_bat * si.m  # coté perpendiculaire au vent longpant
        self.alpha_toit = alpha_toit
        self.alpha_toit2 = alpha_toit2
