#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import os
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

# Ajout du chemin du projet au PYTHONPATH
sys.path.append(str(Path(__file__).parent.parent))

# Import des modules à tester
from ourocode.eurocode.EC0_Combinaison import Combinaison
from ourocode.eurocode.A0_Projet import Model_generator, Projet

# Fixture pour un projet de base
@pytest.fixture
def basic_project():
    return Projet(
        nom="Projet Test",
        auteur="Testeur",
        date="2025-05-27",
        unite_longueur="m",
        unite_force="kN",
        unite_masse="kg",
        unite_temp="°C",
        alt=100,
    )

# Fixture pour un générateur de modèle de base
@pytest.fixture
def basic_model_generator(basic_project):
    model = Model_generator(projet=basic_project)
    # Ajout d'une barre de test
    n1 = model.add_node(X=0, Y=0, Z=0)
    n2 = model.add_node(X=5000, Y=2000, Z=0)
    section = model.add_section(100, 200, 0, "Rectangulaire")
    material = model.add_material_by_class("C24")
    m1 = model.add_member(n1, n2, material, section)
    # Ajout de charges de test
    model.create_dist_load(
        member_id=m1,
        name="Charge permanente",
        start_load=-5,
        end_load=-5,
        start_pos="start",
        end_pos="end",
        action="Permanente G",
        direction="FY"
    )
    model.create_dist_load(
        member_id=m1,
        name="Charge d'exploitation",
        start_load=-3,
        end_load=-3,
        start_pos="start",
        end_pos="end",
        action="Exploitation Q",
        direction="FY"
    )
    model.create_dist_load(
        member_id=m1,
        name="Charge de vent +",
        start_load=-1,
        end_load=-1,
        start_pos="start",
        end_pos="end",
        action="Vent pression W+",
        direction="Fy"
    )
    model.create_dist_load(
        member_id=m1,
        name="Charge de neige",
        start_load=-1,
        end_load=-1,
        start_pos="start",
        end_pos="end",
        action="Neige normale Sn",
        direction="FY"
    )
    return model

class Test_Combinaison:
    def test_initialization(self, basic_model_generator):
        """Teste l'initialisation de la classe Combinaison"""
        combi = Combinaison(
            model_generator=basic_model_generator,
            ELU_STR=True,
            ELU_STR_ACC=False,
            ELS_C=True,
            ELS_QP=True,
            cat="Cat H : toits"
        )
        
        # Vérification des attributs de base
        assert combi.elu_STR is True
        assert combi.elu_STR_ACC is False
        assert combi.els_C is True
        assert combi.els_QP is True
        assert combi.cat == "Cat H : toits"
        assert hasattr(combi, '_combo_tags')
        assert hasattr(combi, 'combiActionVariable')

    def test_coef_psy_property(self, basic_model_generator):
        """Teste la propriété coef_psy"""
        combi = Combinaison(
            model_generator=basic_model_generator,
            ELU_STR=True,
            cat="Cat H : toits"
        )
        
        coefs = combi.coef_psy
        
        # Vérification des clés du dictionnaire
        assert "Cat H : toits" in coefs
        assert "Vent" in coefs
        assert "Température" in coefs
        
        # Vérification des sous-clés
        assert "psy0" in coefs["Cat H : toits"]
        assert "psy1" in coefs["Cat H : toits"]
        assert "psy2" in coefs["Cat H : toits"]

    def test_key_action_psy(self, basic_model_generator):
        """Teste la méthode _key_action_psy"""
        combi = Combinaison(
            model_generator=basic_model_generator,
            ELU_STR=True,
            cat="Cat H : toits"
        )
        
        # Test avec différentes actions
        assert combi._key_action_psy("Q") == "Cat H : toits"
        assert combi._key_action_psy("Sn") == "Neige <= 1000m"
        assert combi._key_action_psy("W+") == "Vent"
        assert combi._key_action_psy("T") == "Température"

    def test_create_dataframe_load(self, basic_model_generator):
        """Teste la création d'un DataFrame de charges"""
        combi = Combinaison(
            model_generator=basic_model_generator,
            ELU_STR=True,
            cat="Cat H : toits"
        )
        
        # Création d'un tableau de test
        array_load = np.empty((0, 10), dtype=object)
        load_dict = {
            "N° barre": 1,
            "Nom": "Charge test",
            "Action": "Permanente G",
            "Type de charge": "Charge répartie",
            "Position": 0,
            "Axe": "Y"
        }
        
        # Ajout d'une charge
        array_load = combi._create_array_load(
            name="COMB_ELU1",
            index=1,
            factor_load=1.35,
            load=-6.75,  # 1.35 * 5 kN/m
            load_dict=load_dict,
            array=array_load
        )
        
        # Création du DataFrame
        df = combi._create_dataframe_load(array_load)
        
        # Vérifications
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert df.iloc[0]["Combinaison"] == "COMB_ELU1"
        assert df.iloc[0]["Facteur de charge"] == 1.35
        assert df.iloc[0]["Charge combinée"] == -6.75
        assert df.iloc[0]["Action"] == "Permanente G"

    def test_generate_combinaison(self, basic_model_generator):
        """Teste la génération des combinaisons"""
        combi = Combinaison(
            model_generator=basic_model_generator,
            ELU_STR=True,
            ELU_STR_ACC=False,
            ELS_C=True,
            ELS_QP=True,
            cat="Cat H : toits"
        )
        
        # Vérification que les tags de combinaisons sont corrects
        assert "ELU_STR" in combi._combo_tags
        assert "ELS_C" in combi._combo_tags
        assert "ELS_QP" in combi._combo_tags
        assert "ELU_STR_ACC" not in combi._combo_tags  # Désactivé dans ce test
        
        # Vérification des actions variables détectées
        assert "G" in combi.combiActionVariable
        assert "Q" in combi.combiActionVariable
        assert "Sn" in combi.combiActionVariable
        assert "W+" in combi.combiActionVariable

