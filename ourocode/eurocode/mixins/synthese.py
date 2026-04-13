# coding in UTF-8
# by Anthony PARISOT
import math as mt
from collections.abc import Mapping, Iterable

import pandas as pd


class SyntheseMixin:
    """Mixin fournissant la synthèse des taux de travail.

    Méthodes :
    - _add_synthese_taux_travail : agrège les taux dans un DataFrame cumulatif
    - synthese_taux_travail : finalise le tableau avec la ligne des maximums
    """

    def _add_synthese_taux_travail(
        self,
        lignes: Iterable,
        col_normale: str="Situation normale",
        col_incendie: str="Situation d'incendie",
        arrondi: int|None=3
    ) -> pd.DataFrame:
        """Construit un tableau pandas de synthèse des taux de travail.

        Cette méthode agrège les résultats de vérification (taux de travail)
        pour la situation normale et la situation d'incendie dans un DataFrame
        partagé entre toutes les méthodes de vérification.

        Args:
            lignes (Iterable): Éléments de synthèse sous forme de :
                - dict avec clés : "effort"/"Effort", "SN"/"situation_normale",
                  "SI"/"situation_incendie"
                - tuple/list : (effort, taux_situation_normale, taux_situation_incendie)
            col_normale (str, optional): Nom de la colonne pour la situation normale.
                Defaults to "Situation normale".
            col_incendie (str, optional): Nom de la colonne pour la situation d'incendie.
                Defaults to "Situation d'incendie".
            arrondi (int|None, optional): Nombre de décimales à conserver.
                None pour ne pas arrondir. Defaults to 3.

        Returns:
            pd.DataFrame: Tableau de synthèse cumulatif avec les taux de travail.
        """
        colonnes = ["Effort agissant dimensionnant", col_normale, col_incendie]
        # DataFrame commun stocké sur l'instance
        dataframe = getattr(self, "_synthese_taux_df", None)
        lignes_preparees = []

        def _arrondi(val):
            if val is None or (isinstance(val, float) and mt.isnan(val)):
                return None
            try:
                val = float(val)
                return round(val, arrondi) if arrondi is not None else val
            except Exception:
                return val

        for ligne in lignes:
            if isinstance(ligne, Mapping):
                effort = ligne.get("effort") or ligne.get("Effort")
                sn = (
                    ligne.get("SN")
                    or ligne.get("situation_normale")
                    or ligne.get(col_normale)
                    or ligne.get("normale")
                )
                si_val = (
                    ligne.get("SI")
                    or ligne.get("situation_incendie")
                    or ligne.get(col_incendie)
                    or ligne.get("incendie")
                )
            elif isinstance(ligne, (list, tuple)):
                effort, sn, si_val = (list(ligne) + [None] * 3)[:3]
            else:
                raise TypeError("Chaque ligne doit être un dict, une liste ou un tuple.")

            lignes_preparees.append({
                "Effort agissant dimensionnant": effort,
                col_normale: _arrondi(sn),
                col_incendie: _arrondi(si_val),
            })

        df_nouveau = pd.DataFrame(lignes_preparees, columns=colonnes)

        if dataframe is not None:
            dataframe = dataframe.copy()
            df_final = pd.concat([dataframe, df_nouveau], join="inner", ignore_index=True)
        else:
            df_final = df_nouveau

        # On mémorise le tableau commun sur l'instance pour les appels suivants
        self._synthese_taux_df = df_final
        return df_final
    
    def synthese_taux_travail(self):
        """Construit un tableau pandas de synthèse des taux de travail.

        Cette méthode agrège les résultats de vérification (taux de travail)
        pour la situation normale et la situation d'incendie dans un DataFrame
        partagé entre toutes les méthodes de vérification.

        Returns:
            pd.DataFrame: Tableau de synthèse final avec la ligne des maximums.

        Note:
            Le DataFrame est partagé entre toutes les méthodes de vérification
            d'une même instance via l'attribut _synthese_taux_df.
        """
        df = getattr(self, "_synthese_taux_df", None)
        if df is None or df.empty:
            return df
        # On suppose que les 2 dernières colonnes sont les taux (situation normale / incendie)
        col_normale, col_incendie = df.columns[-2], df.columns[-1]
        max_normale = pd.to_numeric(df[col_normale], errors="coerce").max(skipna=True)
        max_incendie = pd.to_numeric(df[col_incendie], errors="coerce").max(skipna=True)
        ligne_max = {
            df.columns[0]: "Max",
            col_normale: max_normale,
            col_incendie: max_incendie,
        }
        self._synthese_taux_df = pd.concat([df, pd.DataFrame([ligne_max], columns=df.columns)], ignore_index=True)
        return self._synthese_taux_df
