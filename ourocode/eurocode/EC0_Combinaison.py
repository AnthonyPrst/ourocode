# coding in UTF-8
# by Anthony PARISOT
import os, sys
import numpy as np
import pandas as pd

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# # from A0_Projet import Bar_generator
from ourocode.eurocode.A0_Projet import Projet, Model_generator


class Combinaison(Projet):
    CAT_TYPE = (
        "Aucune",
        "Cat A : habitation",
        "Cat B : bureaux",
        "Cat C : lieu de réunion",
        "Cat D : zones commerciales",
        "Cat E : stockage",
        "Cat F : véhicule <= 30kN",
        "Cat G : véhicule <= 160kN",
        "Cat H : toits",
    )
    COEF_G = (1, 1.35)  # Ginf, Gsup
    COEF_Q = 1.5  # Qsup
    ANALYZE_TYPE = ("Général", "Linéaire", "Second ordre")

    def __init__(
        self,
        model_generator: object,
        ELU_STR: bool = ("True", "False"),
        ELU_STR_ACC: bool = ("True", "False"),
        ELS_C: bool = ("True", "False"),
        ELS_QP: bool = ("True", "False"),
        cat: str = CAT_TYPE,
        kdef: float = None,
        type_psy_2: str = ["Court terme", "Moyen terme", "Long terme"],
        **kwargs,
    ):
        """Créer une classe Combinaison qui génère les combinaisons d'action suivant les actions données.
        Cette classe est hérité de la classe Projet du module A0_Projet.py.
        Elle nécessite toutefois la création d'un modèle de calcul avant d'être instanciée.
        Le modèle de calcul est créé par la classe Model_generator dans le module A0_Projet.py.

        Args:
                model_generator (Model_generator): le modèle de calcul à utiliser pour la génération des combinaisons d'action.
                ELU_STR (bool): combiner les charges à l'ELU structure, si vrai -> True, sinon False.
                ELU_STR_ACC (bool): combiner les charges à l'ELU accidentel, si vrai -> True, sinon False.
                ELS_C (bool): combiner les charges à l'ELS carctéristique, si vrai -> True, sinon False.
                ELS_QP (bool): combiner les charges à l'ELS quasi permananent, si vrai -> True, sinon False.
                cat (str): catégorie d'exploitation de la zone considéré. Defaults to "Aucune".
                kdef (float): Coefficient permettant de prendre en compte le fluage du bois en fonction de sa classe de service.
                        Si le matériaux est autre que du bois laisser vide.
                type_psy_2: détermine le type de psy 2 à récupérer notamment pour le calcul de la flèche dans le bois.
                        Court terme = 0, Moyen terme = calcul du psy en fonction de l'action donnant le psy le plus élevé, Long terme = 1.
        """
        super().__init__(**kwargs)
        self._combo_tags = []
        self._model_generator = model_generator
        self._model_generator.generate_model()
        self.elu_STR = ELU_STR
        self.elu_STR_ACC = ELU_STR_ACC
        self.els_C = ELS_C
        self.els_QP = ELS_QP
        self.cat = cat
        self._name_combination = []
        self._generate_combinaison()
        if kdef:
            self.type_psy_2 = type_psy_2
            self._els_fleche_bois(kdef)

    def _generate_combinaison(self):
        """Génère les combinaisons de chargement si une liste de chargement à été définit précédement."""
        self.combiActionVariable = [0] * 7
        for load in self._model_generator.get_all_loads().values():
            action = load["Action"]
            if action == "Permanente G":
                self.combiActionVariable[0] = "G"

            elif action == "Exploitation Q":
                self.combiActionVariable[1] = "Q"

            elif action == "Neige normale Sn":
                self.combiActionVariable[2] = "Sn"

            elif action == "Vent pression W+":
                self.combiActionVariable[3] = "W+"

            elif action == "Vent dépression W-":
                self.combiActionVariable[4] = "W-"

            elif action == "Neige accidentelle Sx":
                self.combiActionVariable[5] = "Sx"

            else:
                self.combiActionVariable[6] = "Ae"

        if self.elu_STR:
            self._ELU_STR()
            self._combo_tags.append("ELU_STR")
        if self.elu_STR_ACC:
            self._ELU_STR_ACC()
            self._combo_tags.append("ELU_STR_ACC")
        if self.els_C:
            self._ELS_C()
            self._combo_tags.append("ELS_C")
        if self.els_QP:
            self._ELS_QP()
            self._combo_tags.append("ELS_QP")

    @property
    def coef_psy(self):
        """Retourne les caractéristiques psy sous forme de dictionnaire"""
        dict_psy = {self.cat: {}, "Vent": {}, "Température": {}}
        if self.alt.value > 1000:
            dict_psy["Neige > 1000m"] = {}
        else:
            dict_psy["Neige <= 1000m"] = {}

        data_csv_psy = self._data_from_csv("coeff_psy.csv")
        psy_columns = data_csv_psy.columns.to_list()

        for psy_i in psy_columns:
            for key in dict_psy.keys():
                dict_psy[key][psy_i] = data_csv_psy.loc[key].loc[psy_i]
        return dict_psy

    def _key_action_psy(self, action_variable):
        if action_variable == "Q":
            index = self.cat
        elif action_variable == "Sn":
            if self.alt.value > 1000:
                index = "Neige > 1000m"
            else:
                index = "Neige <= 1000m"

        elif action_variable == "W+" or action_variable == "W-":
            index = "Vent"
        else:
            index = "Température"
        return index

    def _create_array_load(
        self,
        name: str,
        index: str,
        factor_load: float,
        load: float | dict,
        load_dict: dict,
        array,
    ):
        load = np.array(
            [
                name,
                str(index),
                load_dict["N° barre"],
                load_dict["Nom"],
                load_dict["Action"],
                load_dict["Type de charge"],
                factor_load,
                load,
                load_dict["Position"],
                load_dict["Axe"],
            ],
            dtype=object,
        )
        array = np.append(array, [load], axis=0)
        return array

    def _create_dataframe_load(self, array_load):
        df_load = pd.DataFrame(
            array_load,
            columns=[
                "Combinaison",
                "Index",
                "N° de barre",
                "Nom",
                "Action",
                "Type",
                "Facteur de charge",
                "Charge combinée",
                "Position",
                "Axe",
            ],
        )
        # print(df_load)
        return df_load
    
    def _create_analyze_type(self, tag: str, df_type_load: pd.DataFrame):
        """Génère les load combo dans le modèle FEM

        Args:
            tag (str): combo_tag à donner dans le modèle
            df_type_load (pd.Dataframe): dataframe des combinaisons associé au tag (ex: ELU_STR, ELS_C)
        """
        unique_combinations = df_type_load["Combinaison"].unique()
        combo_dict = {}
        for combo_name in unique_combinations:
            combo_dict[combo_name] = self._get_combi_factor_load(combo_name)
        self._model_generator._add_load_combos_to_model(combo_dict, tag)


    def _calcs_combined_load(self, factor: float, load: float):
        """Convertit la valeur de la charge en fonction du type de charge."""
        if isinstance(load, dict):
            return {"start": load["start"] * factor, "end": load["end"] * factor}
        else:
            return load * factor

    def _create_list_combination(self, nameCombi: str):
        if nameCombi not in self._name_combination:
            self._name_combination.append(nameCombi)

    def _add_combination(self, combi_name: str, index_load: str, factor: float, value: float| dict, load_dict: dict, array_load):
        """Ajoute une combinaison et met à jour array_load."""
        self._create_list_combination(combi_name)
        return self._create_array_load(combi_name, index_load, factor, value, load_dict, array_load)

    def _ELU_STR(self):
        """Combinaison à l'ELU STR"""
        array_load = np.empty((0, 10))

        for index_load, load in self._model_generator.get_all_loads().items():
            load_action = self.DICO_COMBI_ACTION[load["Action"]]
            if (
                load_action
                not in self.combiActionVariable[5:7]
            ):
                name = "ELU_STR " + load_action
                array_load = self._add_combination(name, index_load, 1, load["Charge"], load, array_load)

            for var_action in self.combiActionVariable:
                if not var_action in ("Sx", "Ae"):
                    if load["Action"] == "Permanente G" and var_action == "G":
                        name = "ELU_STR " + str(self.COEF_G[1]) + "G"
                        value = self._calcs_combined_load(
                            self.COEF_G[1], load["Charge"]
                        )
                        array_load = self._add_combination(name, index_load, self.COEF_G[1], value, load, array_load)

                    elif var_action != 0 and var_action != "G":
                        if var_action != "W-":
                            name = "ELU_STR 1.35G + " + str(self.COEF_Q) + var_action
                            self._create_list_combination(name)

                        if (
                            load["Action"] != "Permanente G"
                            and load_action == var_action
                        ):
                            if var_action == "W-":
                                name = "ELU_STR 1G" + " + " + str(self.COEF_Q) + var_action
                                self._create_list_combination(name)
                            value = self._calcs_combined_load(
                                self.COEF_Q, load["Charge"]
                            )
                            array_load = self._create_array_load(
                                name, index_load, self.COEF_Q, value, load, array_load
                            )

                        elif load["Action"] == "Permanente G":
                            if var_action == "W-":
                                name = "ELU_STR 1G" + " + " + str(self.COEF_Q) + var_action
                                self._create_list_combination(name)
                                factor = self.COEF_G[0]
                                value = self._calcs_combined_load(
                                    factor, load["Charge"]
                                )
                            else:
                                factor = self.COEF_G[1]
                                value = self._calcs_combined_load(
                                    factor, load["Charge"]
                                )
                            array_load = self._create_array_load(
                                name, index_load, factor, value, load, array_load
                            )

                        for index in range(1, 4):
                            if (
                                var_action != self.combiActionVariable[index]
                                and var_action != "W-"
                            ):
                                # selon AN EC1 on ne doit pas combiner les charges d'entretien avec les charges de neige et de vent en toiture
                                if (
                                    var_action == "Q"
                                    or self.combiActionVariable[index] == "Q"
                                ) and self.cat == "Cat H : toits":
                                    continue
                                if self.combiActionVariable[index] != 0:
                                    key_psy = self._key_action_psy(
                                        self.combiActionVariable[index]
                                    )
                                    name = (
                                        "ELU_STR 1.35G + "
                                        + str(self.COEF_Q)
                                        + var_action
                                        + " + "
                                        + str(
                                            round(
                                                self.coef_psy[key_psy]["psy0"]
                                                * self.COEF_Q,
                                                2,
                                            )
                                        )
                                        + self.combiActionVariable[index]
                                    )
                                    self._create_list_combination(name)

                                    if load["Action"] == "Permanente G":
                                        value = self._calcs_combined_load(
                                            self.COEF_G[1], load["Charge"]
                                        )
                                        array_load = self._create_array_load(
                                            name,
                                            index_load,
                                            self.COEF_G[1],
                                            value,
                                            load,
                                            array_load,
                                        )

                                    elif (
                                        load_action
                                        == self.combiActionVariable[index]
                                    ):
                                        factor = (
                                            float(self.coef_psy[key_psy]["psy0"])
                                            * self.COEF_Q
                                        )
                                        value = self._calcs_combined_load(
                                            factor, load["Charge"]
                                        )
                                        array_load = self._create_array_load(
                                            name,
                                            index_load,
                                            factor,
                                            value,
                                            load,
                                            array_load,
                                        )

                                    elif (
                                        load_action == var_action
                                    ):
                                        value = self._calcs_combined_load(
                                            self.COEF_Q, load["Charge"]
                                        )
                                        array_load = self._create_array_load(
                                            name,
                                            index_load,
                                            self.COEF_Q,
                                            value,
                                            load,
                                            array_load,
                                        )

                                    for index2 in range(1, 4):
                                        if var_action != self.combiActionVariable[index2]:
                                            if (
                                                self.combiActionVariable[index]
                                                != self.combiActionVariable[index2]
                                            ):
                                                if self.combiActionVariable.index(
                                                    self.combiActionVariable[index]
                                                ) < self.combiActionVariable.index(
                                                    self.combiActionVariable[index2]
                                                ):
                                                    if (
                                                        self.combiActionVariable[index2]
                                                        != 0
                                                    ):
                                                        key_psy2 = self._key_action_psy(
                                                            self.combiActionVariable[
                                                                index2
                                                            ]
                                                        )
                                                        name = (
                                                            "ELU_STR 1.35G + "
                                                            + str(self.COEF_Q)
                                                            + var_action
                                                            + " + "
                                                            + str(
                                                                round(
                                                                    self.coef_psy[
                                                                        key_psy
                                                                    ]["psy0"]
                                                                    * self.COEF_Q,
                                                                    2,
                                                                )
                                                            )
                                                            + self.combiActionVariable[
                                                                index
                                                            ]
                                                            + " + "
                                                            + str(
                                                                round(
                                                                    self.coef_psy[
                                                                        key_psy2
                                                                    ]["psy0"]
                                                                    * self.COEF_Q,
                                                                    2,
                                                                )
                                                            )
                                                            + self.combiActionVariable[
                                                                index2
                                                            ]
                                                        )
                                                        self._create_list_combination(
                                                            name
                                                        )

                                                        if (
                                                            load["Action"]
                                                            == "Permanente G"
                                                        ):
                                                            value = self._calcs_combined_load(
                                                                self.COEF_G[1],
                                                                load["Charge"],
                                                            )
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    self.COEF_G[1],
                                                                    value,
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

                                                        elif (
                                                            self.DICO_COMBI_ACTION[
                                                                load["Action"]
                                                            ]
                                                            == self.combiActionVariable[
                                                                index2
                                                            ]
                                                        ):
                                                            factor = (
                                                                float(
                                                                    self.coef_psy[
                                                                        key_psy2
                                                                    ]["psy0"]
                                                                )
                                                                * self.COEF_Q
                                                            )
                                                            self._calcs_combined_load(
                                                                factor,
                                                                load["Charge"],
                                                            )
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    factor,
                                                                    value,
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

                                                        elif (
                                                            self.DICO_COMBI_ACTION[
                                                                load["Action"]
                                                            ]
                                                            == self.combiActionVariable[
                                                                index
                                                            ]
                                                        ):
                                                            factor = (
                                                                float(
                                                                    self.coef_psy[
                                                                        key_psy
                                                                    ]["psy0"]
                                                                )
                                                                * self.COEF_Q
                                                            )
                                                            value = self._calcs_combined_load(
                                                                factor,
                                                                load["Charge"],
                                                            )
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    factor,
                                                                    value,
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

                                                        elif (
                                                            self.DICO_COMBI_ACTION[
                                                                load["Action"]
                                                            ]
                                                            == var_action
                                                        ):
                                                            value = self._calcs_combined_load(
                                                                self.COEF_Q,
                                                                load["Charge"],
                                                            )
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    self.COEF_Q,
                                                                    value,
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )
        array_load = array_load[array_load[:, 0].argsort()]
        self._df_load_ELU_STR = self._create_dataframe_load(array_load)
        self._create_analyze_type("ELU_STR", self._df_load_ELU_STR)

    def get_ELU_STR_loads(self):
        """Retourne les charges des combinaisons ELU STR"""
        return self._df_load_ELU_STR

    def _return_combi_ELUSTR(self, combi):
        return self._df_load_ELU_STR.loc[self._df_load_ELU_STR["Combinaison"] == combi]

    def _ELU_STR_ACC(self):
        """Combinaison à l'ELU STR ACC"""
        array_load = np.empty((0, 10))
        acc_action = ("Sx", "Ae")

        for index_load, load in self._model_generator.get_all_loads().items():
            if (
                self.DICO_COMBI_ACTION[load["Action"]]
                not in self.combiActionVariable[1:5]
            ):
                name = "ELU_STR_ACC " + self.DICO_COMBI_ACTION[load["Action"]]
                self._create_list_combination(name)
                array_load = self._create_array_load(
                    name, index_load, 1, load["Charge"], load, array_load
                )

            for action in self.combiActionVariable:
                # if action not in ("W-",):
                if action not in (0, "G"):
                    key_psy_action = self._key_action_psy(action)
                    if action in acc_action:
                        name = "ELU_STR_ACC G + " + action
                    elif self.coef_psy[key_psy_action]["psy1"] != 0:
                        name = (
                            "ELU_STR_ACC_FIRE G + "
                            + str(self.coef_psy[key_psy_action]["psy1"])
                            + action
                        )

                    if name:
                        self._create_list_combination(name)
                        if (
                            load["Action"] != "Permanente G"
                            and self.DICO_COMBI_ACTION[load["Action"]] == action
                        ):
                            if action in acc_action:
                                array_load = self._create_array_load(
                                    name,
                                    index_load,
                                    1,
                                    load["Charge"],
                                    load,
                                    array_load,
                                )
                            else:
                                factor = float(self.coef_psy[key_psy_action]["psy1"])
                                value = self._calcs_combined_load(
                                    factor, load["Charge"]
                                )
                                array_load = self._create_array_load(
                                    name, index_load, factor, value, load, array_load
                                )

                        elif load["Action"] == "Permanente G":
                            array_load = self._create_array_load(
                                name, index_load, 1, load["Charge"], load, array_load
                            )

                    for index in range(1, 4):
                        if action != self.combiActionVariable[index]:
                            if self.combiActionVariable[index] != 0 and action not in (
                                "W-",
                            ):
                                # selon AN EC1 on ne doit pas combiner les charges d'entretien avec les charges de neige et de vent en toiture
                                if (
                                    action == "Q"
                                    or self.combiActionVariable[index] == "Q"
                                ) and self.cat == "Cat H : toits":
                                    continue
                                key_psy = self._key_action_psy(
                                    self.combiActionVariable[index]
                                )
                                if self.coef_psy[key_psy]["psy2"] != 0:
                                    if action in acc_action:
                                        name_2 = (
                                            "ELU_STR_ACC G + "
                                            + action
                                            + " + "
                                            + str(self.coef_psy[key_psy]["psy2"])
                                            + self.combiActionVariable[index]
                                        )
                                    else:
                                        name_2 = (
                                            "ELU_STR_ACC_FIRE G + "
                                            + str(self.coef_psy[key_psy_action]["psy1"])
                                            + action
                                            + " + "
                                            + str(self.coef_psy[key_psy]["psy2"])
                                            + self.combiActionVariable[index]
                                        )

                                    try:
                                        self._create_list_combination(name_2)
                                        if load["Action"] == "Permanente G":
                                            array_load = self._create_array_load(
                                                name_2,
                                                index_load,
                                                1,
                                                load["Charge"],
                                                load,
                                                array_load,
                                            )

                                        elif (
                                            self.DICO_COMBI_ACTION[load["Action"]]
                                            == self.combiActionVariable[index]
                                        ):
                                            factor = float(
                                                self.coef_psy[key_psy]["psy2"]
                                            )
                                            value = self._calcs_combined_load(
                                                factor, load["Charge"]
                                            )
                                            array_load = self._create_array_load(
                                                name_2,
                                                index_load,
                                                factor,
                                                value,
                                                load,
                                                array_load,
                                            )

                                        elif (
                                            self.DICO_COMBI_ACTION[load["Action"]]
                                            == action
                                        ):
                                            if action in acc_action:
                                                array_load = self._create_array_load(
                                                    name_2,
                                                    index_load,
                                                    1,
                                                    load["Charge"],
                                                    load,
                                                    array_load,
                                                )
                                            else:
                                                factor = float(
                                                    self.coef_psy[key_psy_action][
                                                        "psy1"
                                                    ]
                                                )
                                                value = self._calcs_combined_load(
                                                    factor, load["Charge"]
                                                )
                                                array_load = self._create_array_load(
                                                    name_2,
                                                    index_load,
                                                    factor,
                                                    value,
                                                    load,
                                                    array_load,
                                                )
                                    except:
                                        pass

                                    for index2 in range(1, 4):
                                        if action != self.combiActionVariable[
                                            index2
                                        ] and action not in ("W-",):
                                            if (
                                                self.combiActionVariable[index]
                                                != self.combiActionVariable[index2]
                                            ):
                                                if self.combiActionVariable.index(
                                                    self.combiActionVariable[index]
                                                ) < self.combiActionVariable.index(
                                                    self.combiActionVariable[index2]
                                                ):
                                                    if (
                                                        self.combiActionVariable[index2]
                                                        != 0
                                                    ):
                                                        key_psy2 = self._key_action_psy(
                                                            self.combiActionVariable[
                                                                index2
                                                            ]
                                                        )
                                                        if (
                                                            self.coef_psy[key_psy2][
                                                                "psy2"
                                                            ]
                                                            != 0
                                                        ):
                                                            if action in acc_action:
                                                                name_3 = (
                                                                    "ELU_STR_ACC G + "
                                                                    + action
                                                                    + " + "
                                                                    + str(
                                                                        self.coef_psy[
                                                                            key_psy
                                                                        ]["psy2"]
                                                                    )
                                                                    + self.combiActionVariable[
                                                                        index
                                                                    ]
                                                                    + " + "
                                                                    + str(
                                                                        self.coef_psy[
                                                                            key_psy2
                                                                        ]["psy2"]
                                                                    )
                                                                    + self.combiActionVariable[
                                                                        index2
                                                                    ]
                                                                )
                                                            else:
                                                                name_3 = (
                                                                    "ELU_STR_ACC_FIRE G + "
                                                                    + str(
                                                                        self.coef_psy[
                                                                            key_psy_action
                                                                        ]["psy1"]
                                                                    )
                                                                    + action
                                                                    + " + "
                                                                    + str(
                                                                        self.coef_psy[
                                                                            key_psy
                                                                        ]["psy2"]
                                                                    )
                                                                    + self.combiActionVariable[
                                                                        index
                                                                    ]
                                                                    + " + "
                                                                    + str(
                                                                        self.coef_psy[
                                                                            key_psy2
                                                                        ]["psy2"]
                                                                    )
                                                                    + self.combiActionVariable[
                                                                        index2
                                                                    ]
                                                                )

                                                            try:
                                                                self._create_list_combination(
                                                                    name_3
                                                                )

                                                                if (
                                                                    load["Action"]
                                                                    == "Permanente G"
                                                                ):
                                                                    array_load = self._create_array_load(
                                                                        name_3,
                                                                        index_load,
                                                                        1,
                                                                        load["Charge"],
                                                                        load,
                                                                        array_load,
                                                                    )

                                                                elif (
                                                                    self.DICO_COMBI_ACTION[
                                                                        load["Action"]
                                                                    ]
                                                                    == self.combiActionVariable[
                                                                        index2
                                                                    ]
                                                                ):
                                                                    factor = float(
                                                                        self.coef_psy[
                                                                            key_psy2
                                                                        ]["psy2"]
                                                                    )
                                                                    value = self._calcs_combined_load(
                                                                        factor,
                                                                        load["Charge"],
                                                                    )
                                                                    array_load = self._create_array_load(
                                                                        name_3,
                                                                        index_load,
                                                                        factor,
                                                                        value,
                                                                        load,
                                                                        array_load,
                                                                    )

                                                                elif (
                                                                    self.DICO_COMBI_ACTION[
                                                                        load["Action"]
                                                                    ]
                                                                    == self.combiActionVariable[
                                                                        index
                                                                    ]
                                                                ):
                                                                    factor = float(
                                                                        self.coef_psy[
                                                                            key_psy
                                                                        ]["psy2"]
                                                                    )
                                                                    value = self._calcs_combined_load(
                                                                        factor,
                                                                        load["Charge"],
                                                                    )
                                                                    array_load = self._create_array_load(
                                                                        name_3,
                                                                        index_load,
                                                                        factor,
                                                                        value,
                                                                        load,
                                                                        array_load,
                                                                    )

                                                                elif (
                                                                    self.DICO_COMBI_ACTION[
                                                                        load["Action"]
                                                                    ]
                                                                    == action
                                                                ):
                                                                    if (
                                                                        action
                                                                        in acc_action
                                                                    ):
                                                                        array_load = self._create_array_load(
                                                                            name_3,
                                                                            index_load,
                                                                            1,
                                                                            load[
                                                                                "Charge"
                                                                            ],
                                                                            load,
                                                                            array_load,
                                                                        )
                                                                    else:
                                                                        factor = float(
                                                                            self.coef_psy[
                                                                                key_psy_action
                                                                            ][
                                                                                "psy1"
                                                                            ]
                                                                        )
                                                                        value = self._calcs_combined_load(
                                                                            factor,
                                                                            load[
                                                                                "Charge"
                                                                            ],
                                                                        )
                                                                        array_load = self._create_array_load(
                                                                            name_3,
                                                                            index_load,
                                                                            factor,
                                                                            value,
                                                                            load,
                                                                            array_load,
                                                                        )
                                                            except:
                                                                pass
        array_load = array_load[array_load[:, 0].argsort()]
        self._df_load_ELU_STR_ACC = self._create_dataframe_load(array_load)
        self._create_analyze_type("ELU_STR_ACC", self._df_load_ELU_STR_ACC)

    def get_ELU_STR_ACC_loads(self):
        """Retourne les charges de la combinaison ELU STR ACC"""
        return self._df_load_ELU_STR_ACC

    def _return_combi_ELU_STR_ACC(self, combi):
        return self._df_load_ELU_STR_ACC.loc[
            self._df_load_ELU_STR_ACC["Combinaison"] == combi
        ]

    def _ELS_C(self):
        array_load = np.empty((0, 10))

        for index_load, load in self._model_generator.get_all_loads().items():
            for action in self.combiActionVariable:
                if not action in ("Sx", "Ae"):

                    if load["Action"] == "Permanente G" and action == "G":
                        name = "ELS_C G"
                        self._create_list_combination(name)
                        array_load = self._create_array_load(
                            name, index_load, 1, load["Charge"], load, array_load
                        )

                    elif action != 0 and action != "G":
                        name = "ELS_C G + " + action
                        self._create_list_combination(name)

                        if (
                            load["Action"] != "Permanente G"
                            and self.DICO_COMBI_ACTION[load["Action"]] == action
                        ):
                            array_load = self._create_array_load(
                                name, index_load, 1, load["Charge"], load, array_load
                            )

                        elif load["Action"] == "Permanente G":
                            array_load = self._create_array_load(
                                name, index_load, 1, load["Charge"], load, array_load
                            )

                        for index in range(1, 4):
                            if (
                                action != self.combiActionVariable[index]
                                and action != "W-"
                            ):
                                # selon AN EC1 on ne doit pas combiner les charges d'entretien avec les charges de neige et de vent en toiture
                                if (
                                    action == "Q"
                                    or self.combiActionVariable[index] == "Q"
                                ) and self.cat == "Cat H : toits":
                                    continue
                                key_psy = self._key_action_psy(
                                    self.combiActionVariable[index]
                                )

                                if (
                                    self.combiActionVariable[index] != 0
                                    and float(self.coef_psy[key_psy]["psy0"]) != 0
                                ):
                                    name = (
                                        "ELS_C G + "
                                        + action
                                        + " + "
                                        + str(round(self.coef_psy[key_psy]["psy0"], 2))
                                        + self.combiActionVariable[index]
                                    )
                                    self._create_list_combination(name)

                                    if load["Action"] == "Permanente G":
                                        array_load = self._create_array_load(
                                            name,
                                            index_load,
                                            1,
                                            load["Charge"],
                                            load,
                                            array_load,
                                        )

                                    elif (
                                        self.DICO_COMBI_ACTION[load["Action"]]
                                        == self.combiActionVariable[index]
                                    ):
                                        factor = float(self.coef_psy[key_psy]["psy0"])
                                        value = self._calcs_combined_load(
                                            factor, load["Charge"]
                                        )
                                        array_load = self._create_array_load(
                                            name,
                                            index_load,
                                            factor,
                                            value,
                                            load,
                                            array_load,
                                        )

                                    elif (
                                        self.DICO_COMBI_ACTION[load["Action"]] == action
                                    ):
                                        array_load = self._create_array_load(
                                            name,
                                            index_load,
                                            1,
                                            load["Charge"],
                                            load,
                                            array_load,
                                        )

                                    for index2 in range(1, 4):
                                        if action != self.combiActionVariable[index2]:
                                            if (
                                                self.combiActionVariable[index]
                                                != self.combiActionVariable[index2]
                                            ):
                                                if self.combiActionVariable.index(
                                                    self.combiActionVariable[index]
                                                ) < self.combiActionVariable.index(
                                                    self.combiActionVariable[index2]
                                                ):
                                                    if (
                                                        self.combiActionVariable[index2]
                                                        != 0
                                                    ):
                                                        key_psy2 = self._key_action_psy(
                                                            self.combiActionVariable[
                                                                index2
                                                            ]
                                                        )
                                                        name = (
                                                            "ELS_C G + "
                                                            + action
                                                            + " + "
                                                            + str(
                                                                round(
                                                                    self.coef_psy[
                                                                        key_psy
                                                                    ]["psy0"],
                                                                    2,
                                                                )
                                                            )
                                                            + self.combiActionVariable[
                                                                index
                                                            ]
                                                            + " + "
                                                            + str(
                                                                round(
                                                                    self.coef_psy[
                                                                        key_psy2
                                                                    ]["psy0"],
                                                                    2,
                                                                )
                                                            )
                                                            + self.combiActionVariable[
                                                                index2
                                                            ]
                                                        )
                                                        self._create_list_combination(
                                                            name
                                                        )

                                                        if (
                                                            load["Action"]
                                                            == "Permanente G"
                                                        ):
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    1,
                                                                    load["Charge"],
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

                                                        elif (
                                                            self.DICO_COMBI_ACTION[
                                                                load["Action"]
                                                            ]
                                                            == self.combiActionVariable[
                                                                index2
                                                            ]
                                                        ):
                                                            factor = float(
                                                                self.coef_psy[key_psy2][
                                                                    "psy0"
                                                                ]
                                                            )
                                                            value = self._calcs_combined_load(
                                                                factor, load["Charge"]
                                                            )
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    factor,
                                                                    value,
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

                                                        elif (
                                                            self.DICO_COMBI_ACTION[
                                                                load["Action"]
                                                            ]
                                                            == self.combiActionVariable[
                                                                index
                                                            ]
                                                        ):
                                                            factor = float(
                                                                self.coef_psy[key_psy][
                                                                    "psy0"
                                                                ]
                                                            )
                                                            value = self._calcs_combined_load(
                                                                factor, load["Charge"]
                                                            )
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    factor,
                                                                    value,
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

                                                        elif (
                                                            self.DICO_COMBI_ACTION[
                                                                load["Action"]
                                                            ]
                                                            == action
                                                        ):
                                                            value = load["Charge"]
                                                            array_load = (
                                                                self._create_array_load(
                                                                    name,
                                                                    index_load,
                                                                    1,
                                                                    load["Charge"],
                                                                    load,
                                                                    array_load,
                                                                )
                                                            )

        array_load = array_load[array_load[:, 0].argsort()]
        self._df_load_ELScarac = self._create_dataframe_load(array_load)
        self._create_analyze_type("ELS_C", self._df_load_ELScarac)

    def get_ELS_C_loads(self):
        """Retourne les charges des combinaisons ELS CARACTERISTIQUE"""
        return self._df_load_ELScarac

    def _return_combi_ELScarac(self, combi):
        return self._df_load_ELScarac.loc[self._df_load_ELScarac["Combinaison"] == combi]

    def _ELS_QP(self):
        array_load = np.empty((0, 10))

        for index_load, load in self._model_generator.get_all_loads().items():
            for action in self.combiActionVariable:

                if load["Action"] == "Permanente G" and action == "G":
                    name = "ELS_QP G"
                    self._create_list_combination(name)
                    array_load = self._create_array_load(
                        name, index_load, 1, load["Charge"], load, array_load
                    )

                elif action != 0 and action != "G":

                    key_psy = self._key_action_psy(action)
                    if float(self.coef_psy[key_psy]["psy2"]) != 0:
                        name = (
                            "ELS_QP G + "
                            + str(round(self.coef_psy[key_psy]["psy2"], 2))
                            + action
                        )
                        self._create_list_combination(name)

                        if (
                            load["Action"] != "Permanente G"
                            and self.DICO_COMBI_ACTION[load["Action"]] == action
                        ):
                            factor = float(self.coef_psy[key_psy]["psy2"])
                            value = self._calcs_combined_load(factor, load["Charge"])
                            array_load = self._create_array_load(
                                name, index_load, factor, value, load, array_load
                            )

                        elif load["Action"] == "Permanente G":
                            array_load = self._create_array_load(
                                name, index_load, 1, load["Charge"], load, array_load
                            )

                        for index in range(1, 4):
                            if (
                                action != self.combiActionVariable[index]
                                and action != "W-"
                            ):
                                if self.combiActionVariable.index(
                                    action
                                ) < self.combiActionVariable.index(
                                    self.combiActionVariable[index]
                                ):
                                    if self.combiActionVariable[index] != 0:
                                        key_psy1 = self._key_action_psy(
                                            self.combiActionVariable[index]
                                        )
                                        if float(self.coef_psy[key_psy1]["psy2"]) != 0:
                                            name = (
                                                "ELS_QP G + "
                                                + str(
                                                    round(
                                                        self.coef_psy[key_psy]["psy2"],
                                                        2,
                                                    )
                                                )
                                                + action
                                                + " + "
                                                + str(
                                                    round(
                                                        self.coef_psy[key_psy1]["psy2"],
                                                        2,
                                                    )
                                                )
                                                + self.combiActionVariable[index]
                                            )
                                            self._create_list_combination(name)
                                            if load["Action"] == "Permanente G":
                                                array_load = self._create_array_load(
                                                    name,
                                                    index_load,
                                                    1,
                                                    load["Charge"],
                                                    load,
                                                    array_load,
                                                )

                                            elif (
                                                self.DICO_COMBI_ACTION[load["Action"]]
                                                == self.combiActionVariable[index]
                                            ):
                                                factor = float(
                                                    self.coef_psy[key_psy1]["psy2"]
                                                )
                                                value = self._calcs_combined_load(
                                                    factor, load["Charge"]
                                                )
                                                array_load = self._create_array_load(
                                                    name,
                                                    index_load,
                                                    factor,
                                                    value,
                                                    load,
                                                    array_load,
                                                )

                                            elif (
                                                self.DICO_COMBI_ACTION[load["Action"]]
                                                == action
                                            ):
                                                factor = float(
                                                    self.coef_psy[key_psy]["psy2"]
                                                )
                                                value = self._calcs_combined_load(
                                                    factor, load["Charge"]
                                                )
                                                array_load = self._create_array_load(
                                                    name,
                                                    index_load,
                                                    factor,
                                                    value,
                                                    load,
                                                    array_load,
                                                )

                                            for index2 in range(1, 4):
                                                if (
                                                    action
                                                    != self.combiActionVariable[index2]
                                                ):
                                                    if (
                                                        self.combiActionVariable[index]
                                                        != self.combiActionVariable[
                                                            index2
                                                        ]
                                                    ):
                                                        if self.combiActionVariable.index(
                                                            self.combiActionVariable[
                                                                index
                                                            ]
                                                        ) < self.combiActionVariable.index(
                                                            self.combiActionVariable[
                                                                index2
                                                            ]
                                                        ):
                                                            if (
                                                                self.combiActionVariable[
                                                                    index2
                                                                ]
                                                                != 0
                                                            ):
                                                                key_psy2 = self._key_action_psy(
                                                                    self.combiActionVariable[
                                                                        index2
                                                                    ]
                                                                )
                                                                if (
                                                                    float(
                                                                        self.coef_psy[
                                                                            key_psy2
                                                                        ]["psy2"]
                                                                    )
                                                                    != 0
                                                                ):
                                                                    name = (
                                                                        "ELS_QP G + "
                                                                        + action
                                                                        + " + "
                                                                        + str(
                                                                            round(
                                                                                self.coef_psy[
                                                                                    key_psy1
                                                                                ][
                                                                                    "psy2"
                                                                                ],
                                                                                2,
                                                                            )
                                                                        )
                                                                        + self.combiActionVariable[
                                                                            index
                                                                        ]
                                                                        + " + "
                                                                        + str(
                                                                            round(
                                                                                self.coef_psy[
                                                                                    key_psy2
                                                                                ][
                                                                                    "psy2"
                                                                                ],
                                                                                2,
                                                                            )
                                                                        )
                                                                        + self.combiActionVariable[
                                                                            index2
                                                                        ]
                                                                    )
                                                                    self._create_list_combination(
                                                                        name
                                                                    )

                                                                    if (
                                                                        load["Action"]
                                                                        == "Permanente G"
                                                                    ):
                                                                        array_load = self._create_array_load(
                                                                            name,
                                                                            index_load,
                                                                            1,
                                                                            load[
                                                                                "Charge"
                                                                            ],
                                                                            load,
                                                                            array_load,
                                                                        )

                                                                    elif (
                                                                        self.DICO_COMBI_ACTION[
                                                                            load[
                                                                                "Action"
                                                                            ]
                                                                        ]
                                                                        == self.combiActionVariable[
                                                                            index2
                                                                        ]
                                                                    ):
                                                                        factor = float(
                                                                            self.coef_psy[
                                                                                key_psy2
                                                                            ][
                                                                                "psy2"
                                                                            ]
                                                                        )
                                                                        value = self._calcs_combined_load(
                                                                            factor,
                                                                            load[
                                                                                "Charge"
                                                                            ],
                                                                        )
                                                                        array_load = self._create_array_load(
                                                                            name,
                                                                            index_load,
                                                                            factor,
                                                                            value,
                                                                            load,
                                                                            array_load,
                                                                        )

                                                                    elif (
                                                                        self.DICO_COMBI_ACTION[
                                                                            load[
                                                                                "Action"
                                                                            ]
                                                                        ]
                                                                        == self.combiActionVariable[
                                                                            index
                                                                        ]
                                                                    ):
                                                                        factor = float(
                                                                            self.coef_psy[
                                                                                key_psy1
                                                                            ][
                                                                                "psy2"
                                                                            ]
                                                                        )
                                                                        value = self._calcs_combined_load(
                                                                            factor,
                                                                            load[
                                                                                "Charge"
                                                                            ],
                                                                        )
                                                                        array_load = self._create_array_load(
                                                                            name,
                                                                            index_load,
                                                                            factor,
                                                                            value,
                                                                            load,
                                                                            array_load,
                                                                        )

                                                                    elif (
                                                                        self.DICO_COMBI_ACTION[
                                                                            load[
                                                                                "Action"
                                                                            ]
                                                                        ]
                                                                        == action
                                                                    ):
                                                                        factor = float(
                                                                            self.coef_psy[
                                                                                key_psy
                                                                            ][
                                                                                "psy2"
                                                                            ]
                                                                        )
                                                                        value = self._calcs_combined_load(
                                                                            factor,
                                                                            load[
                                                                                "Charge"
                                                                            ],
                                                                        )
                                                                        array_load = self._create_array_load(
                                                                            name,
                                                                            index_load,
                                                                            factor,
                                                                            value,
                                                                            load,
                                                                            array_load,
                                                                        )

        array_load = array_load[array_load[:, 0].argsort()]
        self._df_load_ELSqp = self._create_dataframe_load(array_load)
        self._create_analyze_type("ELS_QP", self._df_load_ELSqp)

    def get_ELS_QP_loads(self):
        """Retourne les charges des combinaisons ELS Quasi Permanente"""
        return self._df_load_ELSqp

    def _return_combi_ELSqp(self, combi):
        return self._df_load_ELSqp.loc[self._df_load_ELSqp["Combinaison"] == combi]

    def get_psy_2_by_combination(self, name: str):
        """Récupère le psy 2 en fonction de l'action la plus défavorable soit le psy2 le plus élevé

        Args:
                name (str): nom de la combinaison
        """
        match self.type_psy_2:
            case "Court terme":
                psy_2 = 0
            case "Long terme":
                psy_2 = 1
            case _:
                psy_2 = 0
                if name[0:6] == "ELS_QP":
                    for action in ["Q", "S"]:
                        try:
                            if name[8:].index(action):
                                key_psy = self._key_action_psy(action)
                                psy_2 = max(psy_2, self.coef_psy[key_psy]["psy2"])
                        except ValueError:
                            pass
        return psy_2

    def _els_fleche_bois(self, kdef: float):
        """Génère l'association des combinaisons ELS caractéristique + quasi permanente (fluage, intégrant le Kdef) pour le calcul d'un élément en bois."""

        if self.els_C:
            self._combo_tags.append("W_inst_Q")
            # On détermine W_inst(Q), pour cela on enlève W_inst_G
            self._df_W_inst_Q = self._df_load_ELScarac[
                self._df_load_ELScarac["Action"] != "Permanente G"
            ]
            for index in range(self._df_W_inst_Q.shape[0]):
                name = self._df_W_inst_Q.iloc[index, 0]
                name_combi = "W_inst " + name[10:]
                self._df_W_inst_Q.iloc[index, 0] = name_combi
                self._create_list_combination(name_combi)
            self._create_analyze_type("W_inst_Q", self._df_W_inst_Q)
            # print(self._df_W_inst_Q)

            if self.els_QP:
                self._combo_tags.append("W_net_fin")
                combi_qp = {}
                combi_c = {}
                value_search = ["Q", "S"]

                for name in self._name_combination:
                    if name[0:6] == "ELS_QP":
                        typeQP = [0] * 2
                        for i in range(len(value_search)):
                            try:
                                if name[8:].index(value_search[i]):
                                    typeQP[i] = value_search[i]
                            except ValueError:
                                pass
                        combi_qp[name] = typeQP

                    elif name[0:5] == "ELS_C":
                        typeC = [0] * 2
                        for i in range(len(value_search)):
                            try:
                                if name[7:].index(value_search[i]):
                                    if (
                                        value_search[i] == "Q"
                                        and self.cat == "Cat H : toits"
                                    ):
                                        pass
                                    elif (
                                        value_search[i] == "S"
                                        and self.alt.value <= 1000
                                    ):
                                        pass
                                    else:
                                        typeC[i] = value_search[i]
                            except ValueError:
                                pass
                        combi_c[name] = typeC

                list_of_key = list(combi_qp.keys())
                list_of_value = list(combi_qp.values())

                array_load = np.empty((0, 10))
                for name_combi_C, val in combi_c.items():
                    position = list_of_value.index(val)
                    name_combi_QP = list_of_key[position]
                    name_combi = "W_net_fin " + name_combi_C + " & " + name_combi_QP
                    psy_2 = self.get_psy_2_by_combination(name_combi_QP)

                    df_c = self._df_load_ELScarac[
                        self._df_load_ELScarac["Combinaison"] == name_combi_C
                    ]
                    df_qp = self._df_load_ELSqp[
                        self._df_load_ELSqp["Combinaison"] == name_combi_QP
                    ]
                    for index in range(df_c.shape[0]):
                        index_c = df_c.iloc[index, 1]
                        beam_number = df_c.iloc[index, 2]
                        name_load = df_c.iloc[index, 3]
                        action_load = df_c.iloc[index, 4]
                        type_load = df_c.iloc[index, 5]
                        factor_c = df_c.iloc[index, 6]
                        valeur_c = df_c.iloc[index, 7]
                        position = df_c.iloc[index, 8]
                        axe = df_c.iloc[index, 9]
                        if len(df_qp[df_qp["Index"] == index_c]):
                            # On multiplie les combinaisons quasi permanente par Kdef pour trouver Wcreep et par la formule 2.3.2.2 equ 2.7/2.10 pour modifier artificiellement E,mean en E,mean,fin
                            factor_kdef = (1 + kdef * psy_2) * kdef
                            factor_qp = df_qp[df_qp["Index"] == index_c].iloc[0, 6]
                            total_factor = factor_c + factor_kdef * factor_qp
                            if type_load == "Distribuée":
                                valeur = {"start": 0, "end": 0}
                                for pos, load in df_qp[df_qp["Index"] == index_c].iloc[0, 7].items():
                                    valeur_qp = load * factor_kdef
                                    valeur[pos] += (valeur_qp + valeur_c[pos])
                                    
                            else:
                                valeur_qp = (
                                    df_qp[df_qp["Index"] == index_c].iloc[0, 7]
                                    * factor_kdef
                                )
                                valeur = valeur_c + valeur_qp
                            
                        else:
                            valeur = valeur_c
                            total_factor = factor_c
                        load = np.array(
                            [
                                name_combi,
                                index_c,
                                beam_number,
                                name_load,
                                action_load,
                                type_load,
                                total_factor,
                                valeur,
                                position,
                                axe,
                            ],
                            dtype=object,
                        )
                        array_load = np.append(array_load, [load], axis=0)
                        self._create_list_combination(name_combi)

                array_load = array_load[array_load[:, 0].argsort()]
                self._df_W_net_fin = self._create_dataframe_load(array_load)
                self._create_analyze_type("W_net_fin", self._df_W_net_fin)
                # print(self._df_W_net_fin)

    def get_W_inst_Q_loads(self):
        """Retourne les charges des combinaisons Winst (Q)"""
        if hasattr(self, "_df_W_inst_Q"):
            return self._df_W_inst_Q
        return None

    def _return_combi_W_inst_Q(self, combi):
        return self._df_W_inst_Q.loc[self._df_W_inst_Q["Combinaison"] == combi]

    def get_W_net_fin_loads(self):
        """Retourne les charges des combinaisons Wnet fin"""
        if hasattr(self, "_df_W_net_fin"):
            return self._df_W_net_fin
        return None

    def _return_combi_W_net_fin(self, combi):
        return self._df_W_net_fin.loc[self._df_W_net_fin["Combinaison"] == combi]

    @property
    def list_combination(self):
        """Retourne un data frame avec toute les combinaison créer"""
        self._name_combination.sort()
        return pd.DataFrame(self._name_combination, columns=["Combinaison"])

    def get_list_combination(self, type: str=("Toutes", "ELU_ALL", "ELU_STR", "ELU_STR_ACC", "ELS_ALL", "ELS_C", "ELS_QP", "W")):
        """Retourne la liste des combinaisons pour le type sélectionné.

        Args:
            type (str): Le type de combinaison à retourner. Defaults to ("Toutes", "ELU_ALL", "ELU_STR", "ELU_STR_ACC", "ELS_ALL", "ELS_C", "ELS_QP", "Flèche bois").

        Returns:
            list: Liste des combinaisons
        """
        match type:
            case "Toutes":
                return self._name_combination
            case "ELU_ALL":
                return [name for name in self._name_combination if name.startswith("ELU")]
            case "ELU_STR":
                return [name for name in self._name_combination if name.startswith("ELU_STR ")]
            case "ELU_STR_ACC":
                return [name for name in self._name_combination if name.startswith("ELU_STR_ACC")]
            case "ELS_ALL":
                return [name for name in self._name_combination if name.startswith(("ELS", "W"))]
            case "ELS_C":
                return [name for name in self._name_combination if name.startswith("ELS_C")]
            case "ELS_QP":
                return [name for name in self._name_combination if name.startswith("ELS_QP")]
            case "Flèche bois":
                return [name for name in self._name_combination if name.startswith("W")]

    def _choice_combi_df(self):
        shape = len(self.list_combination)
        dict_load_combi = {}

        for i in range(shape):
            combi = self.list_combination.iloc[i, 0]

            if self.elu_STR_ACC and combi[0:11] == "ELU_STR_ACC":
                df_combi = self._return_combi_ELU_STR_ACC(combi)

            elif self.elu_STR and combi[0:7] == "ELU_STR":
                df_combi = self._return_combi_ELUSTR(combi)

            elif self.els_C and combi[0:5] == "ELS_C":
                df_combi = self._return_combi_ELScarac(combi)

            elif self.els_QP and combi[0:6] == "ELS_QP":
                df_combi = self._return_combi_ELSqp(combi)

            elif self.els_C and combi[0:6] == "W_inst":
                df_combi = self._return_combi_W_inst_Q(combi)

            elif self.els_C and self.els_QP and combi[0:9] == "W_net_fin":
                df_combi = self._return_combi_W_net_fin(combi)

            df_combi = df_combi.drop(labels=["Combinaison"], axis=1)
            load_list_combi = df_combi.values.tolist()
            dict_load_combi[combi] = load_list_combi

        return dict_load_combi

    def get_combi_list_load(self, nom: str):
        """Retourne la liste des charges combinées pour la combinaison sélectionné

        Args:
                nom (str): nom de la combinaison à récupérer. Defaults to "Sélectionner tout".
        """
        self.list_loads = self._choice_combi_df()[nom]
        return self.list_loads
    
    def _get_combi_factor_load(self, nom: str) -> dict:
        """Retourne la liste des facteurs de charge pour la combinaison sélectionné

        Args:
                nom (str): nom des la combinaison à récupérer. Defaults to "Sélectionner tout".
        """
        self.list_loads = self.get_combi_list_load(nom)
        dict_loads = {}
        for load in self.list_loads:
            key_load = self.DICO_COMBI_ACTION[load[3]]
            dict_loads[key_load] = load[5]
        return dict_loads


    def min_type_load(self, name_combi: str) -> str:
        """Retourne le type de chargement de plus courte durée

        Args:
                name_combi (str): Combinaison à analyser
        """

        dictName = {
            "G": "Permanente",
            "Q": "Moyen terme",
            "Sn": "Court terme",
            "W+": "Instantanee",
            "W-": "Instantanee",
            "Sx": "Instantanee",
            "Ae": "Instantanee",
        }
        for action in self.combiActionVariable:
            if action:
                indexAction = name_combi[8:].find(action)
                if indexAction > -1:
                    if action == "Sn" and self.alt.value >= 1000:
                        name_load_type = dictName["Q"]
                    elif action == "Sn" and self.cat == "Cat H : toits":
                        name_load_type = dictName["Sn"]
                    else:
                        name_load_type = dictName[action]
        return name_load_type
    

if __name__ == "__main__":

    _dict_loads = [
        [1, "G", "Permanente G", -10, "0/4000", "Z"],
        [2, "", "Neige accidentelle Sx", -200, "0/6000", "Z"],
        [2, "", "Neige normale Sn", -200, "0/6000", "Z"],
        [3, "", "Exploitation Q", -100, 3610, "Z"],
    ]
    projet = Projet("AP", "6018.0", "", "", 73215, "France", 1200)
    chargement = Model_generator._from_parent_class(projet)
    chargement.create_load_by_list(_dict_loads)
    # "Cat A : habitation"
    # "Cat H : toits"
    c1 = Combinaison._from_parent_class(
        chargement,
        ELU_STR=True,
        ELU_STR_ACC=False,
        ELS_C=True,
        ELS_QP=True,
        cat="Cat H : toits",
        kdef=0.6,
        type_psy_2="Moyen terme",
    )
    rcombi = "ELS_QP G + 0.3Q"
    # print(c1._return_combi_ELUSTR(rcombi))
    print(pd.DataFrame(c1.coef_psy))
    # print(c1._df_load_ELScarac)
    # print(c1._df_load_ELSqp)
    # print(c1._df_W_inst_Q)
    # print(c1._df_load_ELU_STR_ACC)

    print(c1.list_combination)
    print(c1.get_combi_list_load("ELU_STR 1.35G + 1.5Q"))
