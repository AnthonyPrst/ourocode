# coding in UTF-8
# by Anthony PARISOT
"""Vérification EC5 d'un modèle MEF — classe :class:`Verification_EC5`.

Ce module fournit une classe de haut niveau qui orchestre :

    MEF (Model_generator + Model_result)  +  Combinaisons (Combinaison)
             |
             v
    ** boucle sur toutes les combinaisons ELU **
             |
             v
    Barre / Flexion / Cisaillement / Traction / Compression  (EC5)
             +
    Flèche ELS (W_inst_Q, W_net_fin)

Pour chaque barre structurale (:meth:`Model_generator.group_members` ou
:meth:`auto_group_continuous_members`) :

1. Pour chaque combinaison ELU (ELU_STR + ELU_STR_ACC) :
   - ``kmod`` est déterminé via :meth:`Combinaison.min_type_load` (durée de
     chargement la plus courte présente dans la combinaison, EN 1995-1-1 §2.3.1.2) ;
   - ``typecombi`` (``"Fondamentales"`` ou ``"Accidentelles"``) est
     déterminé via :meth:`Combinaison.type_combi` pour le coefficient ``γM`` ;
   - les efforts gouvernants (Nx signé, Vy, My, Mz) sont extraits pour cette
     combinaison précise via :meth:`Model_result.get_min_max_internal_force` et
     :meth:`Model_result.get_absolute_internal_force` ;
   - les taux de travail Flexion / Cisaillement / Traction ou Compression
     sont calculés.
2. Pour chaque type de vérification, le taux maximal rencontré et la
   combinaison gouvernante associée sont retenus.
3. La flèche ELS est calculée séparément (pas de boucle par combo : on
   prend la flèche max sur tous les ``W_inst_Q`` et ``W_net_fin``).

Les paramètres de vérification (type de bâtiment, position de charge,
coefficients de longueur efficace…) sont configurés **globalement** à
l'initialisation et peuvent être **surchargés individuellement** via
:meth:`verify`.

Utilisation typique :

.. code-block:: python

    combi = Combinaison(model, ELU_STR=True, ELS_C=True, ELS_QP=True, kdef=0.6, ...)
    result = Model_result(model, ...)
    verif = Verification_EC5(
        combi, result,
        type_bat="Bâtiment courant",
    )
    verif.synthese()             # DataFrame agrégé de toutes les barres
    verif.verify("SM1")          # Détail d'une barre structurale
    verif.verify("SM2", coeflef_y=0.7)  # Surcharge individuelle
"""
from __future__ import annotations

import pandas as pd

from ourocode.eurocode.core.projet import Projet
from ourocode.eurocode.ec5.element_droit import (
    Barre,
    Flexion,
    Traction,
    Compression,
    Cisaillement,
)


# ---------------------------------------------------------------------------
# Conversions d'unités
# ---------------------------------------------------------------------------
# Pynite renvoie les efforts en unités cohérentes : N (forces), N·mm (moments)
# car le modèle MEF est alimenté en E [MPa], longueurs [mm], A [mm²], I [mm⁴].
# - Force : N -> kN : /1000
# - Moment : N·mm -> kN·m : /1e6
# - Flèche via Physical ``.value`` : base SI = m -> mm : *1000

def _length_to_mm(physical) -> float:
    return float(physical.value) * 1000.0


def _as_tag_list(tags):
    """Normalise ``tags`` en liste : Pynite considère une chaîne comme itérable."""
    if tags is None:
        return None
    if isinstance(tags, (list, tuple)):
        return list(tags)
    return [tags]


# ---------------------------------------------------------------------------
# Classe principale
# ---------------------------------------------------------------------------

class Verification_EC5(Projet):
    """Vérifie les barres structurales d'un modèle selon l'EC5 (ELU + Flèche ELS).

    Hérite de :class:`~ourocode.eurocode.core.projet.Projet` pour s'insérer
    dans la hiérarchie ourocode et bénéficier des méthodes communes
    (``_from_parent_class``, persistence, etc.). La classe stocke en interne
    la :class:`~ourocode.eurocode.core.combinaison.Combinaison` et le
    :class:`~ourocode.eurocode.core.model_result.Model_result` fournis, et
    expose :meth:`verify` et :meth:`synthese` pour lancer les vérifications.

    Les paramètres de vérification (coefficients de longueur efficace,
    types d'appui pour compression) sont automatiquement récupérés depuis
    le champ ``Design`` de chaque barre structurale (déterminés par
    ``group_members``). Les arguments ``type_bat``, ``type_ele``, ``pos_charge``,
    ``Hi``, ``Hf``, ``effet_systeme``, ``classe_service`` définis à l'initialisation
    s'appliquent à toutes les barres et peuvent être surchargés ponctuellement
    via :meth:`verify`.

    Pour chaque combinaison ELU, :meth:`Combinaison.min_type_load` et
    :meth:`Combinaison.type_combi` sont utilisées pour déterminer
    dynamiquement ``kmod`` et ``γM`` (EN 1995-1-1 §2.3.1.2 et §2.4.1).
    Le taux maximal est ensuite retenu par type de vérification avec la
    combinaison gouvernante associée.
    """

    #: Tag par défaut pour filtrer les combinaisons ELU à vérifier.
    DEFAULT_ELU_FILTER = "ELU_ALL"
    #: Tag ELS pour la flèche instantanée W_inst(Q).
    DEFAULT_W_INST_TAG = "W_inst_Q"
    #: Tag ELS pour la flèche nette finale W_net,fin.
    DEFAULT_W_NET_FIN_TAG = "W_net_fin"
    #: Valeurs acceptées pour ``elu_filter``.
    ELU_FILTERS = ("ELU_ALL", "ELU_STR", "ELU_STR_ACC")

    def __init__(
        self,
        combinaison: object,
        model_result: object,
        elu_filter: str = ELU_FILTERS,
        Hi: int = 12,
        Hf: int = 12,
        effet_systeme: bool = ("False", "True"),
        classe_service: int = Barre.CS,
        type_bat: str = Barre.TYPE_BAT,
        type_ele: str = Barre.TYPE_ELE,
        pos_charge: str = Flexion.LOAD_POS,
        **kwargs,
    ):
        """Initialise une vérification EC5.

        Args:
            combinaison (Combinaison): Instance :class:`~ourocode.eurocode.core.combinaison.Combinaison`
                déjà instanciée et dont les combos ont été générées.
            model_result (Model_result): Instance :class:`~ourocode.eurocode.core.model_result.Model_result`
                déjà analysée.
            elu_filter (str, optional): Tag pour filtrer les combinaisons ELU.
                Valeurs acceptées : ``"ELU_ALL"``, ``"ELU_STR"``, ``"ELU_STR_ACC"``.
                Defaults to ``"ELU_ALL"``.
            Hi (int, optional): Humidité initiale du bois en service (%).
                Defaults to ``12``.
            Hf (int, optional): Humidité finale du bois en service (%).
                Defaults to ``12``.
            effet_systeme (bool, optional): Active l'effet système (EN 1995-1-1 §6.7).
                Defaults to ``False``.
            classe_service (str, optional): Classe de service globale (``Barre.CS``).
                Defaults to 1.
            type_bat (str, optional): Type de bâtiment pour les limites de flèche
                (``Barre.TYPE_BAT``). Defaults to ``Barre.TYPE_BAT``.
            type_ele (str, optional): Type d'élément pour les limites de flèche
                (``Barre.TYPE_ELE``). Peut être surchargé par le champ ``Design``
                de la barre structurale. Defaults to ``Barre.TYPE_ELE``.
            pos_charge (str, optional): Position de la charge verticale par rapport
                à la section (``Flexion.LOAD_POS``). Defaults to ``Flexion.LOAD_POS``.
            **kwargs: Transmis à :class:`~ourocode.eurocode.core.projet.Projet`
                (``ingenieur``, ``name``, ``code_INSEE``, ``alt``…).

        Note:
            Les paramètres de vérification (coefficients de longueur efficace,
            types d'appui pour la compression) sont désormais automatiquement
            récupérés depuis le champ ``Design`` de chaque barre structurale
            (déterminés par :meth:`Model_generator.group_members`).
        """
        super().__init__(**kwargs)
        self._combinaison = combinaison
        self._result = model_result
        if combinaison is not None and not hasattr(self, "_model_generator"):
            self._model_generator = combinaison._model_generator

        # Stockage des paramètres de vérification globaux
        self.elu_filter = elu_filter
        self.type_bat = type_bat
        self.type_ele = type_ele
        self.pos_charge = pos_charge
        self.Hi = Hi
        self.Hf = Hf
        self.effet_systeme = effet_systeme
        self.classe_service = int(classe_service)

        # Cache des objets EC5 par (structural_member_name, combo_name)
        # Peuplé automatiquement lors de chaque appel à verify()
        self._objects_cache: dict[tuple[str, str], dict] = {}


    # ------------------------------------------------------------------ #
    # API publique
    # ------------------------------------------------------------------ #
    def verify(
        self,
        name: str,
        type_ele: str = Barre.TYPE_ELE,
        pos_charge: str = Flexion.LOAD_POS,
        Hi: int = 12,
        Hf: int = 12,
        effet_systeme: bool = ("False", "True"),
        classe_service: int = Barre.CS,
        coeflef_y: float = None,
        coeflef_z: float = None,
        type_appuis_y: str = None,
        type_appuis_z: str = None,
    ) -> dict:
        """Vérifie une barre structurale en bouclant sur toutes les combos ELU.

        Pour chaque combinaison du filtre, les efforts sont extraits puis les
        taux Flexion / Cisaillement / Traction ou Compression sont calculés
        avec le ``kmod`` et la ``typecombi`` appropriés. Le taux maximal sur
        toutes les combos est retenu par vérification. La flèche ELS est
        calculée séparément via les tags ``W_inst_Q`` et ``W_net_fin``.

        Les paramètres ``coeflef_y``, ``coeflef_z`` ainsi que les types d'appui
        pour la compression sont automatiquement récupérés depuis le champ
        ``Design`` du membre structural (déterminés par ``group_members``).
        Les arguments fournis ici permettent de surcharger ces valeurs.

        Args:
            name (str): Nom de la barre structurale (clé dans
                :meth:`Model_generator.get_structural_member`).
            type_ele (str, optional): Surcharge du type d'élément pour les
                limites de flèche. Priorité inférieure au champ ``Design``
                de la barre structurale. Defaults to ``self.type_ele``.
            pos_charge (str, optional): Surcharge de la position de charge.
                Defaults to ``self.pos_charge``.
            Hi (int, optional): Surcharge de l'humidité initiale (%).
                Defaults to 12.
            Hf (int, optional): Surcharge de l'humidité finale (%).
                Defaults to 12.
            effet_systeme (bool, optional): Surcharge de l'effet système.
                Defaults to 1.
            classe_service (int, optional): Surcharge de la classe de service.
                Defaults to 1.
            coeflef_y (float, optional): Surcharge du coef. de longueur efficace y.
                Si ``None`` (défaut), la valeur du champ ``Design`` est utilisée
                (fallback 0.9 si absente).
            coeflef_z (float, optional): Surcharge du coef. de longueur efficace z.
                Même logique que ``coeflef_y``.
            type_appuis_y (str, optional): Surcharge du type d'appui de flambement
                selon y (ex: ``"Rotule - Rotule"``, ``"Encastré - Rotule"``…).
                Si ``None`` (défaut), la valeur du champ ``Design`` est utilisée
                (fallback ``"Rotule - Rotule"`` si absente).
            type_appuis_z (str, optional): Surcharge du type d'appui de flambement
                selon z. Même logique que ``type_appuis_y``.

        Returns:
            dict: ``{"name", "taux", "dataframe"}`` où :

            - ``"name"`` (str) : nom de la barre structurale ;
            - ``"taux"`` (dict) : ``{verif_name: {"taux": float, "combinaison": str}}``
              pour chaque vérification active (``None`` si non applicable) ;
            - ``"dataframe"`` (pd.DataFrame) : tableau synthétique avec colonnes
              ``["Barre structurale", "Vérification", "Taux", "Statut", "Combinaison"]``.

        Raises:
            ValueError: Si la combinaison ou le résultat MEF est absent, ou si
                la section est définie manuellement (``b``/``h`` inconnus).
            KeyError: Si ``name`` n'existe pas dans le modèle.
        """
        self._check_state()

        w_inst_tag = self.DEFAULT_W_INST_TAG
        w_net_fin_tag = self.DEFAULT_W_NET_FIN_TAG
        sm = self._model_generator.get_structural_member(name)
        member_ids = sm["Barres FEM"]

        params = self._resolve_member_params(
            sm, Hi, Hf, effet_systeme, classe_service,
            coeflef_y=coeflef_y, coeflef_z=coeflef_z,
            type_appuis_y=type_appuis_y, type_appuis_z=type_appuis_z,
        )
        L_mm = params["L_mm"]
        lo_y = params["lo_y"]
        lo_z = params["lo_z"]
        lo_flamb_y = params["lo_flamb_y"]
        lo_flamb_z = params["lo_flamb_z"]
        coeflef_y = params["coeflef_y"]
        coeflef_z = params["coeflef_z"]
        type_appuis_y = params["type_appuis_y"]
        type_appuis_z = params["type_appuis_z"]
        barre_kwargs = params["barre_kwargs"]

        # --- Accumulateur de taux ELU sur toutes les combos
        max_taux: dict[str, dict] = {
            "Flexion": {"taux": 0.0, "combinaison": None},
            "Cisaillement": {"taux": 0.0, "combinaison": None},
            "Traction": {"taux": 0.0, "combinaison": None},
            "Compression": {"taux": 0.0, "combinaison": None},
        }

        # --- Boucle sur les combinaisons ELU
        combos = self._combinaison.get_list_combination(self.elu_filter) or []
        for combo_name in combos:
            taux_combo, objects = self._verify_combo_elu(
                combo_name=combo_name,
                member_ids=member_ids,
                barre_kwargs=barre_kwargs,
                lo_y=lo_y, lo_z=lo_z,
                lo_flamb_y=lo_flamb_y, lo_flamb_z=lo_flamb_z,
                coeflef_y=coeflef_y, coeflef_z=coeflef_z,
                pos_charge=pos_charge,
                type_appuis_y=type_appuis_y,
                type_appuis_z=type_appuis_z,
                return_objects=True,
            )
            # Stockage dans le cache pour get_combo_objects()
            self._objects_cache[(name, combo_name)] = objects
            for verif, taux in taux_combo.items():
                if taux is not None and taux > max_taux[verif]["taux"]:
                    max_taux[verif] = {"taux": float(taux), "combinaison": combo_name}

        # --- Flèche ELS (combo-indépendante : on prend le max sur tous les W_*)
        fleche_results = self._verify_fleche(
            member_ids=member_ids,
            barre_kwargs=barre_kwargs,
            L_mm=L_mm,
            type_ele=type_ele,
            type_bat=self.type_bat,
            w_inst_tag=w_inst_tag,
            w_net_fin_tag=w_net_fin_tag,
        )

        taux_final: dict = {}
        for verif, data in max_taux.items():
            taux_final[verif] = data if data["taux"] > 0 else None
        taux_final.update(fleche_results)

        return {
            "name": name,
            "taux": taux_final,
            "dataframe": self._taux_to_df(name, taux_final),
        }

    def get_combo_objects(
        self,
        name: str,
        combo_name: str,
    ) -> dict:
        """Retourne les objets EC5 instanciés pour une barre et une combinaison.

        Lit directement depuis le cache interne peuplé par :meth:`verify`.
        Aucun recalcul n'est effectué.

        Args:
            name (str): Nom de la barre structurale.
            combo_name (str): Nom de la combinaison ELU
                (ex: ``"ELU_STR 1.35G + 1.5Q"``).

        Returns:
            dict: ``{"Flexion", "Cisaillement", "Traction", "Compression"}``
                où chaque valeur est l'objet EC5 instancié et calculé,
                ou ``None`` si la vérification n'est pas applicable
                (ex: pas de compression → ``"Compression"`` est ``None``).

        Raises:
            KeyError: Si ``verify(name)`` n'a pas été appelé au préalable,
                ou si ``combo_name`` n'appartient pas au filtre ELU utilisé.

        Exemple:
            >>> verif.verify("SM1")
            >>> objs = verif.get_combo_objects("SM1", "ELU_STR 1.35G + 1.5Q")
            >>> flexion = objs["Flexion"]
            >>> display(Latex(flexion.taux_m_d()[0]))
        """
        self._check_state()
        key = (name, combo_name)
        if key not in self._objects_cache:
            raise KeyError(
                f"Aucun résultat en cache pour '{name}' / '{combo_name}'. "
                f"Appelez d'abord verify('{name}') pour peupler le cache."
            )
        return self._objects_cache[key]

    def synthese(self) -> pd.DataFrame:
        """Agrège les vérifications EC5 de toutes les barres structurales.

        Parcourt toutes les barres structurales via :meth:`Model_generator._iter_structural_members`
        et concatène les tableaux individuels en un unique DataFrame trié par
        barre structurale, puis par taux décroissant au sein de chaque barre.

        Les paramètres de vérification sont ceux définis à l'initialisation.
        Pour des paramètres spécifiques à une barre, utilisez :meth:`verify` directement.

        Returns:
            pd.DataFrame: Colonnes ``["Barre structurale", "Vérification",
                "Taux", "Statut", "Combinaison"]``, triées par barre structurale
                puis par taux décroissant. En cas d'erreur sur une barre
                (section manuelle, classe inconnue…), un DataFrame
                ``["Barre structurale", "Erreur"]`` est retourné à la place
                si aucune barre n'a pu être vérifiée.
        """
        self._check_state()
        frames: list[pd.DataFrame] = []
        errors: list[tuple[str, str]] = []
        for name, _ in self._model_generator._iter_structural_members():
            try:
                res = self.verify(name,
                    type_ele=self.type_ele,
                    pos_charge=self.pos_charge,
                    Hi=self.Hi,
                    Hf=self.Hf,
                    effet_systeme=self.effet_systeme,
                    classe_service=self.classe_service)
                if not res["dataframe"].empty:
                    frames.append(res["dataframe"])
            except (ValueError, KeyError) as exc:
                errors.append((name, str(exc)))

        if not frames:
            if errors:
                return pd.DataFrame(
                    errors, columns=["Barre structurale", "Erreur"]
                )
            return pd.DataFrame(
                columns=[
                    "Barre structurale", "Vérification",
                    "Taux", "Statut", "Combinaison",
                ]
            )
        df = pd.concat(frames, ignore_index=True)
        return df.sort_values(by=["Barre structurale", "Taux"], ascending=[True, False]).reset_index(drop=True)

    # ------------------------------------------------------------------ #
    # Vérification d'une combinaison ELU (coeur de l'algorithme)
    # ------------------------------------------------------------------ #
    def _verify_combo_elu(
        self,
        combo_name: str,
        member_ids: list[str],
        barre_kwargs: dict,
        lo_y: float, lo_z: float,
        lo_flamb_y: float, lo_flamb_z: float,
        coeflef_y: float, coeflef_z: float,
        pos_charge: str,
        type_appuis_y: str,
        type_appuis_z: str,
        return_objects: bool = False,
    ) -> dict:
        """Retourne les taux ELU (Flexion, Cisaillement, Traction ou Compression)
        pour UNE combinaison spécifique, avec le ``kmod`` et la ``typecombi``
        propres à cette combinaison.

        La clef ``Traction`` et ``Compression`` sont mutuellement exclusives :
        c'est le signe de Nx dominant en valeur absolue sur la combo qui
        tranche. La flexion inclut systématiquement l'interaction
        flexo-compression (eq. 6.19/6.20/6.23/6.24) ou flexo-traction
        (eq. 6.17/6.18) via :meth:`Flexion.taux_m_d`.

        Args:
            return_objects (bool): Si True, retourne aussi les objets EC5
                instanciés sous la clé ``"objects"`` du dictionnaire.
                Defaults to False.
        """
        load_time = self._combinaison.min_type_load(combo_name)
        typecombi = self._combinaison.type_combi(combo_name)

        efforts = self._result.get_min_max_internal_force(member_ids, combo_name)
        Vy_abs = self._result.get_absolute_internal_force(member_ids, combo_name, "Vy", get_combo_name=False).value
        Vz_abs = self._result.get_absolute_internal_force(member_ids, combo_name, "Vz", get_combo_name=False).value
        My_abs = self._result.get_absolute_internal_force(member_ids, combo_name, "My", get_combo_name=False).value
        Mz_abs = self._result.get_absolute_internal_force(member_ids, combo_name, "Mz", get_combo_name=False).value
        # efforts = self._efforts_for_combo(member_ids, combo_name, n_points)

        taux_combo = {
            "Flexion": None,
            "Cisaillement": None,
            "Traction": None,
            "Compression": None,
        }
        objects: dict[str, object] = {
            "Flexion": None,
            "Cisaillement": None,
            "Traction": None,
            "Compression": None,
        }

        barre = Barre(**barre_kwargs)

        # --- Traction / Compression (exclusif, signe-dépendant)
        traction_obj = None
        compression_obj = None
        N_pos = efforts["Nx"]["Max"][0].value  # N, > 0 si compresion
        N_neg = efforts["Nx"]["Min"][0].value   # N, < 0 si traction
        if N_pos > 1e-9 and N_pos >= abs(N_neg):
            compression_obj = Compression._from_parent_class(
                barre,
                lo_y=lo_flamb_y, lo_z=lo_flamb_z,
                type_appuis_y=type_appuis_y, type_appuis_z=type_appuis_z,
            )
            compression_obj.f_c_0_d(loadtype=load_time, typecombi=typecombi)
            compression_obj.sigma_c_0_d(Fc0d=abs(N_pos / 10**3))
            compression_obj.taux_c_0_d()
            taux_combo["Compression"] = max(compression_obj.taux_c_0_rd.values())
            objects["Compression"] = compression_obj
        elif abs(N_neg) > 1e-9:
            traction_obj = Traction._from_parent_class(barre)
            traction_obj.f_t_0_d(loadtype=load_time, typecombi=typecombi)
            traction_obj.sigma_t_0_d(Ft0d=N_neg / 10**3)
            traction_obj.taux_t_0_d()
            taux_combo["Traction"] = max(traction_obj.taux_t_0_rd.values())
            objects["Traction"] = traction_obj

        # --- Flexion (incluant l'interaction)
        flexion_obj = None
        if Mz_abs > 1e-9 or My_abs > 1e-9:
            flexion_obj = Flexion._from_parent_class(
                barre,
                lo_rel_y=lo_y, lo_rel_z=lo_z,
                coeflef_y=coeflef_y, coeflef_z=coeflef_z,
                pos=pos_charge,
            )
            flexion_obj.f_m_d(loadtype=load_time, typecombi=typecombi)
            flexion_obj.sigma_m_d(
                My=My_abs / 10**3, Mz=Mz_abs / 10**3,
            )
            flexion_obj.taux_m_d(
                compression=compression_obj, traction=traction_obj,
            )
            taux_combo["Flexion"] = max(flexion_obj.taux_m_rd.values())
            objects["Flexion"] = flexion_obj

        # --- Cisaillement
        cisaillement_obj = None
        if Vy_abs > 1e-9:
            cisaillement_obj = Cisaillement._from_parent_class(barre)
            cisaillement_obj.f_v_d(loadtype=load_time, typecombi=typecombi)
            cisaillement_obj.tau_d(Vd=Vy_abs / 10**3)
            cisaillement_obj.taux_tau_d()
            taux_combo["Cisaillement"] = max(cisaillement_obj.taux_tau_rd.values())
            objects["Cisaillement"] = cisaillement_obj

        if return_objects:
            return taux_combo, objects
        return taux_combo

    # ------------------------------------------------------------------ #
    # Helpers — résolution des paramètres d'un membre structural
    # ------------------------------------------------------------------ #
    def _resolve_member_params(
        self,
        sm: dict,
        Hi, Hf,
        effet_systeme,
        classe_service,
        coeflef_y: float = None,
        coeflef_z: float = None,
        type_appuis_y: str = None,
        type_appuis_z: str = None,
    ) -> dict:
        """Résout tous les paramètres de calcul d'un membre structural.

        Centralise la lecture de la section, du matériau et des valeurs
        ``Design`` pour éviter la duplication entre :meth:`verify` et
        :meth:`get_combo_objects`.

        Returns:
            dict: Clés ``L_mm``, ``lo_y``, ``lo_z``, ``lo_flamb_y``,
                ``lo_flamb_z``, ``coeflef_y``, ``coeflef_z``,
                ``type_appuis_y``, ``type_appuis_z``, ``barre_kwargs``.

                ``lo_y``/``lo_z`` sont les longueurs de déversement (flexion,
                EC5 §6.3.3) ; ``lo_flamb_y``/``lo_flamb_z`` sont les longueurs
                de flambement (compression, EC5 §6.3.2).
        """
        b_mm, h_mm, section_type = self._resolve_section(sm)
        classe = self._resolve_classe_bois(sm)
        L_mm = _length_to_mm(sm["Longueur totale"])
        lo_y = self._lo_mm(sm, "lo_rel_y", L_mm)
        lo_z = self._lo_mm(sm, "lo_rel_z", L_mm)
        lo_flamb_y = self._lo_mm(sm, "lo_flamb_y", lo_y)
        lo_flamb_z = self._lo_mm(sm, "lo_flamb_z", lo_z)
        coeflef_y = coeflef_y if coeflef_y is not None else self._design_or_default(sm, "coeflef_y", 0.9)
        coeflef_z = coeflef_z if coeflef_z is not None else self._design_or_default(sm, "coeflef_z", 0.9)
        type_appuis_y = type_appuis_y if type_appuis_y is not None else self._design_or_default(sm, "type_appuis_y", "Rotule - Rotule")
        type_appuis_z = type_appuis_z if type_appuis_z is not None else self._design_or_default(sm, "type_appuis_z", "Rotule - Rotule")
        barre_kwargs = dict(
            b=b_mm, h=h_mm, section=section_type,
            Hi=Hi, Hf=Hf, classe=classe, cs=classe_service,
            effet_systeme=effet_systeme,
        )
        return {
            "L_mm": L_mm,
            "lo_y": lo_y,
            "lo_z": lo_z,
            "lo_flamb_y": lo_flamb_y,
            "lo_flamb_z": lo_flamb_z,
            "coeflef_y": coeflef_y,
            "coeflef_z": coeflef_z,
            "type_appuis_y": type_appuis_y,
            "type_appuis_z": type_appuis_z,
            "barre_kwargs": barre_kwargs,
        }

    # ------------------------------------------------------------------ #
    # Flèche ELS (pas de boucle par combo : tag-based)
    # ------------------------------------------------------------------ #
    def _verify_fleche(
        self,
        member_ids: list[str],
        barre_kwargs: dict,
        L_mm: float,
        type_ele,
        type_bat: str,
        w_inst_tag: str,
        w_net_fin_tag: str,
    ) -> dict:
        """Calcule les taux de flèche ELS max via ``Barre.fleche``.

        Retourne ``{"Flèche W_inst(Q)": {...}, "Flèche W_net,fin": {...}}``
        avec des entrées ``None`` si la combinaison ou le type d'élément
        est indisponible.
        """
        out = {
            "Flèche W_inst(Q)": None,
            "Flèche W_net,fin": None,
        }
        if type_ele is None:
            return out
        Ed_WinstQ = self._deflection(member_ids, w_inst_tag)
        Ed_Wnetfin = self._deflection(member_ids, w_net_fin_tag)
        if Ed_WinstQ == 0 and Ed_Wnetfin == 0:
            return out
        barre_fleche = Barre(**barre_kwargs)
        barre_fleche.fleche(
            long=L_mm,
            Ed_WinstQ=_length_to_mm(Ed_WinstQ["Flèche"]),
            Ed_Wnetfin=_length_to_mm(Ed_Wnetfin["Flèche"]),
            Ed_Wfin=0, Ed_W2=0,
            type_ele=type_ele,
            type_bat=type_bat,
        )
        els = barre_fleche.taux_ELS
        if "Winst(Q)" in els:
            out["Flèche W_inst(Q)"] = {
                "taux": float(els["Winst(Q)"]),
                "combinaison": Ed_WinstQ["Combinaison"],
            }
        if "Wnet,fin" in els:
            out["Flèche W_net,fin"] = {
                "taux": float(els["Wnet,fin"]),
                "combinaison": Ed_Wnetfin["Combinaison"],
            }
        return out

    # ------------------------------------------------------------------ #
    # Helpers — extraction efforts par combinaison via Pynite
    # ------------------------------------------------------------------ #
    def _efforts_for_combo(
        self, member_ids: list[str], combo_name: str, n_points: int,
    ) -> dict:
        """Retourne les efforts gouvernants pour une combo en kN / kN·m.

        Utilise :meth:`Model_result.get_internal_force` qui cible une
        combinaison précise par ``combo_name`` (et non par tag). Itère sur
        toutes les barres FEM et tous les points de discrétisation.

        Returns:
            dict: ``{"Nx_min", "Nx_max", "Vy_abs", "My_abs", "Mz_abs"}``
                (kN ou kN·m, ``Nx_min/max`` signés).
        """
        nx_min = 0.0
        nx_max = 0.0
        vy_abs = 0.0
        my_abs = 0.0
        mz_abs = 0.0
        for mid in member_ids:
            _, nx_arr = self._result.get_internal_force(mid, combo_name, "Nx", n_points)
            _, vy_arr = self._result.get_internal_force(mid, combo_name, "Vy", n_points)
            _, my_arr = self._result.get_internal_force(mid, combo_name, "My", n_points)
            _, mz_arr = self._result.get_internal_force(mid, combo_name, "Mz", n_points)
            nx_min = min(nx_min, float(min(nx_arr)))
            nx_max = max(nx_max, float(max(nx_arr)))
            vy_abs = max(vy_abs, float(max(abs(v) for v in vy_arr)))
            my_abs = max(my_abs, float(max(abs(v) for v in my_arr)))
            mz_abs = max(mz_abs, float(max(abs(v) for v in mz_arr)))
        # Pynite : forces en N, moments en N·mm -> kN et kN·m
        return {
            "Nx_min": nx_min / 1_000.0,
            "Nx_max": nx_max / 1_000.0,
            "Vy_abs": vy_abs / 1_000.0,
            "My_abs": my_abs / 1_000_000.0,
            "Mz_abs": mz_abs / 1_000_000.0,
        }

    def _deflection(self, member_ids, tag, direction: str = "dy") -> float:
        """Flèche absolue max en mm pour un tag donné (tous combos du tag)."""
        phys = self._result.get_absolute_max_deflection(
            member_ids, _as_tag_list(tag), direction, get_combo_name=True,
        )
        return phys

    # ------------------------------------------------------------------ #
    # Helpers — résolution de la géométrie / matériau depuis le modèle
    # ------------------------------------------------------------------ #
    def _resolve_section(self, sm: dict) -> tuple[float, float, str]:
        first_mid = sm["Barres FEM"][0]
        section_id = self._model_generator.get_member(first_mid)["Section"]
        section = self._model_generator.get_section(section_id)
        if section["Section"] == "Manuel" or "b" not in section:
            raise ValueError(
                f"Section '{section_id}' définie par propriétés : "
                f"dimensions b, h inconnues. Utilisez `add_section`."
            )
        b_mm = float(section["b"].value) * 1000.0
        h_mm = float(section["h"].value) * 1000.0
        return b_mm, h_mm, section["Section"]

    def _resolve_classe_bois(self, sm: dict) -> str:
        design = sm.get("Design") or {}
        if design.get("classe_bois"):
            return design["classe_bois"]
        first_mid = sm["Barres FEM"][0]
        mat_id = self._model_generator.get_member(first_mid)["Matériaux"]
        mat = self._model_generator.get_material(mat_id)
        if mat["classe"] == "Manuel":
            raise ValueError(
                f"Matériau '{mat_id}' manuel : classe EC5 inconnue. "
                f"Renseignez `classe_bois` ou utilisez `add_material_by_class`."
            )
        return mat_id

    @staticmethod
    def _design_or_default(sm: dict, key: str, default):
        val = (sm.get("Design") or {}).get(key)
        return default if val is None else val

    @classmethod
    def _lo_mm(cls, sm: dict, key: str, default_mm: float) -> float:
        val = (sm.get("Design") or {}).get(key)
        if val is None:
            return default_mm
        if hasattr(val, "value"):
            return _length_to_mm(val)
        return float(val)

    # ------------------------------------------------------------------ #
    # Autres
    # ------------------------------------------------------------------ #
    def _check_state(self) -> None:
        if self._combinaison is None:
            raise ValueError(
                "Verification_EC5 : aucune Combinaison fournie. "
            )
        if self._result is None:
            raise ValueError(
                "Verification_EC5 : aucun Model_result fourni. "
            )

    @staticmethod
    def _taux_to_df(name: str, taux: dict) -> pd.DataFrame:
        rows = []
        for verif, data in taux.items():
            if data is None:
                continue
            t = float(data["taux"])
            rows.append(
                {
                    "Barre structurale": name,
                    "Vérification": verif,
                    "Taux": round(t, 3),
                    "Statut": "OK" if t <= 1.0 else "NOK",
                    "Combinaison": data["combinaison"],
                }
            )
        return pd.DataFrame(
            rows,
            columns=[
                "Barre structurale", "Vérification",
                "Taux", "Statut", "Combinaison",
            ],
        )
