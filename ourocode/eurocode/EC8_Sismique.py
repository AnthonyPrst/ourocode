#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc
from ourocode.eurocode.A0_Projet import Batiment

class Sismique(Batiment):
    CAT_IMPORTANCE = {"I": 0.8, "II": 1, "III": 1.2, "IV": 1.4}
    CLASSE_SOL = {"A": 0.5, "B": 0.6, "C": 0.7, "D": 0.8, "E": 0.9}
    def __init__(self, cat_importance: str = CAT_IMPORTANCE, classe_sol: str = CLASSE_SOL, **kwargs):
        """
        Créer une classe qui permet de calculer l'action sismique selon la méthode des forces latérales.
        Attention, tout les bâtiments ne ce prêtent pas à ce type d'étude (voir EN 1998).
        Cette classe est hérité de la classe Batiment du module A0_Projet.py.

        Args:
            cat_importance (str): Catégorie d'importance du bâtiment:
                - I : Bâtiments dans lesquels il n'y a aucune activité humaine nécessitant un séjour de longue durée.
                - II : 
                    - Habitations individuelles.
                    - Établissements recevant du public (ERP) de catégories 4 et 5.
                    - Habitations collectives de hauteur inférieure à 28 m.
                    - Bureaux ou établissements commerciaux non ERP, h ≤ 28 m, max. 300 pers. 
                    - Bâtiments industriels pouvant accueillir au plus 300 personnes.
                    - Parcs de stationnement ouverts au public.
                - III : 
                    - ERP de catégories 1, 2 et 3.
                    - Habitations collectives et bureaux, h > 28 m.
                    - Bâtiments pouvant accueillir plus de 300 personnes. 
                    - Établissements sanitaires et sociaux.
                    - Centres de production collective d'énergie.
                    - Établissements scolaires.
                - IV : 
                    - Bâtiments indispensables à la sécurité civile, la défense nationale et le maintien de l'ordre public.
                    - Bâtiments assurant le maintien des communications, la production et le stockage d'eau potable, la distribution publique de l'énergie.
                    - Bâtiments assurant le contrôle de la sécurité aérienne.
                    - Établissements de santé nécessaires à la gestion de crise.
                    - Centres météorologiques.
            classe_sol (str): Classe du sol:
                - A : 
                - B : 
                - C : 
                - D : 
                - E : 
        """
        super().__init__(**kwargs)
        self.cat_importance = cat_importance
        self.classe_sol = classe_sol
        self.gamma_1 = self.CAT_IMPORTANCE[self.cat_importance]
        if not self._has_to_be_analyzed():
            raise ValueError("Le bâtiment n'est pas à analyser au niveau sismique")

    @property
    def cat_importance_table(self):
        file = "carte_action_region.csv"
        data_cat_imp = self._data_from_csv(file, index_col=1)
        return data_cat_imp

    @property
    def region_sismique(self):
        file = "carte_action_region.csv"
        df = self._data_from_csv(file, index_col=1)
        return df.loc[str(self.code_INSEE)]["Alea_sismique"]
 
    def _has_to_be_analyzed(self):
        """Retourne True si le bâtiment doit être analysé selon la catégorie d'importance"""
        if self.region_sismique == "Zone 1":
            return False
        elif self.cat_importance == "I":
            return False
        elif self.cat_importance == "II" and self.region_sismique == "Zone 2":
            return False
        else:
            return True