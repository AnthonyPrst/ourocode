#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import pytest
import forallpeople as si
import pandas as pd
import numpy as np

si.environment("structural")
sys.path.insert(1, "./")
from ourocode.eurocode import EC5_CVT as EC5
from ourocode.eurocode.A0_Projet import Batiment

@pytest.fixture
def batiment():
    """Fixture pour créer une instance de Batiment pour les tests."""
    return Batiment(h_bat=3, d_bat=10.0, b_bat=20.0, alpha_toit=10)

class TestMOB:
    """Classe de test pour la classe MOB (Mur Ossature Bois)."""
    
    @pytest.fixture
    def mob_instance(self, batiment):
        """Fixture pour créer une instance de MOB pour les tests."""
        mob = EC5.MOB._from_parent_class(batiment)
        return mob
    
    def test_add_wall(self, mob_instance):
        """Test l'ajout d'un système de mur."""
        sys_name = mob_instance.add_wall(
            h_etage=3.0,
            h_sys_MOB=2.7,
            l_sys_MOB=10.0,
            etage="RDC",
            name="MurPrincipal"
        )
        
        assert sys_name == "MurPrincipal"
        assert "MurPrincipal" in mob_instance._data_MOB
        assert mob_instance._data_MOB["MurPrincipal"]["Etage"] == "RDC"
        assert mob_instance._data_MOB["MurPrincipal"]["Hauteur étage"] == 3.0 * si.m
        assert mob_instance._data_MOB["MurPrincipal"]["Hauteur du système de mur"] == 2.7 * si.m
        assert mob_instance._data_MOB["MurPrincipal"]["Longueur du système de mur"] == 10.0 * si.m
        assert mob_instance._data_MOB["MurPrincipal"]["Murs"] == {}
    
    def test_add_internal_wall(self, mob_instance):
        """Test l'ajout d'un mur interne."""
        # D'abord ajouter un système de mur
        sys_name = mob_instance.add_wall(
            h_etage=3.0,
            h_sys_MOB=2.7,
            l_sys_MOB=10.0,
            etage="RDC"
        )
        
        # Ajouter un mur interne
        wall_name = mob_instance.add_internal_wall(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            Kser_ext=695,
            Kser_int=695,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
            name="MurInterne1"
        )
        
        assert len(mob_instance.data[sys_name]["Murs"]) == 1
        wall = mob_instance.data[sys_name]["Murs"][wall_name]
        assert wall_name == "MurInterne1"
        assert wall["Position en X"] == 1.5 * si.m
        assert wall["Longueur"] == 3.0 * si.m
        assert wall["Kser ext"] == 695 * si.N / si.mm
        assert wall["Kser int"] == 695 * si.N / si.mm
        assert wall["Couturage ext. rive"] == 150 * si.mm
        assert wall["Couturage ext. intermédiaire"] == 300 * si.mm
        assert wall["Couturage int. rive"] == 150 * si.mm
        assert wall["Couturage int. intermédiaire"] == 300 * si.mm
        assert wall["Panneaux"] == {}
    
    def test_add_panel_to_wall(self, mob_instance):
        """Test l'ajout d'un panneau à un mur."""
        # Ajouter un système de mur
        sys_name = mob_instance.add_wall(
            h_etage=3.0,
            h_sys_MOB=2.7,
            l_sys_MOB=10.0,
            etage="RDC"
        )
        
        # Ajouter un mur interne
        wall_name = mob_instance.add_internal_wall(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            Kser_ext=695,
            Kser_int=695,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )
        
        # Ajouter un panneau
        panel_name = mob_instance.add_panel_to_wall(
            sys_name=sys_name,
            wall_name=wall_name,
            number=2,
            e_panel=15,
            h_panel=2700,
            b_panel=1250,
            position_panel="Extérieur",
            type_panel="OSB/3 11-18 mm"
        )
        
        # Vérifier que le panneau a été correctement ajouté
        wall = mob_instance.data[sys_name]["Murs"][wall_name]
        panel = wall["Panneaux"][panel_name]
        assert len(wall["Panneaux"]) == 1
        assert panel["Nombre"] == 2
        assert panel["Epaisseur"] == 15 * si.mm
        assert panel["Hauteur"] == 2700 * si.mm
        assert panel["Largeur"] == 1250 * si.mm
        assert panel["Type"] == "OSB/3 11-18 mm"
    
    def test_save_and_load_wall_data(self, batiment, mob_instance, tmp_path):
        """Test la sauvegarde et le chargement des données des murs."""
        # Ajouter un système de mur avec un mur et un panneau
        sys_name = mob_instance.add_wall(
            h_etage=3.0,
            h_sys_MOB=2.7,
            l_sys_MOB=10.0,
            etage="RDC"
        )
        
        wall_name = mob_instance.add_internal_wall(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            Kser_ext=695,
            Kser_int=695,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )
        
        panel_name = mob_instance.add_panel_to_wall(
            sys_name=sys_name,
            wall_name=wall_name,
            number=2,
            e_panel=15,
            h_panel=2700,
            b_panel=1250,
            position_panel="Extérieur",
            type_panel="OSB/3 11-18 mm"
        )
        
        # Sauvegarder dans un fichier temporaire
        test_file = tmp_path / "test_walls.json"
        mob_instance.save_walls_data(test_file)
        
        # Vérifier que le fichier a été créé
        assert test_file.exists()
        
        # Créer une nouvelle instance et charger les données
        new_mob = EC5.MOB._from_parent_class(batiment)
        new_mob.load_walls_data(test_file)
        print(new_mob.data)
        
        # Vérifier que les données ont été correctement chargées
        assert sys_name in new_mob.data
        assert len(new_mob.data[sys_name]["Murs"]) == 1
        assert len(new_mob.data[sys_name]["Murs"][wall_name]["Panneaux"]) == 1
        assert new_mob.data[sys_name]["Murs"][wall_name]["Kser ext"] == 695 * si.N / si.mm
        assert new_mob.data[sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Nombre"] == 2
        assert new_mob.data[sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Epaisseur"] == 15 * si.mm
        assert new_mob.data[sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Hauteur"] == 2700 * si.mm
        assert new_mob.data[sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Largeur"] == 1250 * si.mm
        assert new_mob.data[sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Type"] == "OSB/3 11-18 mm"
    
    def test_K_panel(self, mob_instance):
        """Test le calcul de la raideur des panneaux."""
        # Ajouter un système de mur
        sys_name = mob_instance.add_wall(
            h_etage=2.8,
            h_sys_MOB=2.5,
            l_sys_MOB=10.0,
            etage="RDC"
        )
        
        # Ajouter un mur avec des propriétés de couturage
        wall_name = mob_instance.add_internal_wall(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            Kser_ext=695,
            Kser_int=695,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
        )
        
        # Ajouter des panneaux de test
        panel1 = mob_instance.add_panel_to_wall(
            sys_name=sys_name,
            wall_name=wall_name,
            number=2,
            e_panel=15,
            h_panel=2500,
            b_panel=1196,
            position_panel="Extérieur",
            type_panel="OSB/3 11-18 mm"
        )
        
        panel2 = mob_instance.add_panel_to_wall(
            sys_name=sys_name,
            wall_name=wall_name,
            number=1,
            e_panel=15,
            h_panel=2500,
            b_panel=1070,
            position_panel="Extérieur",
            type_panel="OSB/3 11-18 mm"
        )
        
        # Appeler la méthode à tester
        df_kp = mob_instance._k_panel()
        
        # Vérifications
        # 1. Vérifier que les DataFrames ne sont pas vides
        # print(df_kp)
        assert not df_kp.empty
        
        # 2. Vérifier les colonnes du DataFrame des raideurs
        expected_columns = ['Hauteur', 'Type', 'Kser_p unit.', 'Coeff. équivalence']
        assert all(col in df_kp.columns for col in expected_columns)
        
        # 3. Vérifier que les raideurs sont positives
        assert (df_kp['Kser_p unit.'] > 0).all()
        
        # 4. Vérifier que le coefficient d'équivalence est compris entre 0 et 1
        assert (df_kp['Coeff. équivalence'] > 0).all()
        assert (df_kp['Coeff. équivalence'] <= 1.0).all()
        
        # 5. Vérifier que le panneau le plus large a un coefficient d'équivalence de 1.0
        max_width_idx = df_kp['Largeur'].idxmax()
        assert np.isclose(df_kp.loc[max_width_idx, 'Coeff. équivalence'], 1.0)
        
        return sys_name, wall_name
    
    def test_k_wall(self, mob_instance):
        """Test le calcul de la raideur des murs et systèmes de murs."""
        # D'abord, configurer les données de test en utilisant test_K_panel
        sys_name, wall_name = self.test_K_panel(mob_instance)
        
        # Ajouter un deuxième mur avec des panneaux différents
        wall_name2 = mob_instance.add_internal_wall(
            sys_name=sys_name,
            position=4.5,
            l_MOB=3.0,
            Kser_ext=695,
            Kser_int=0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
        )
        
        # Ajouter des panneaux au deuxième mur
        mob_instance.add_panel_to_wall(
            sys_name=sys_name,
            wall_name=wall_name2,
            number=3,
            e_panel=15,
            h_panel=2500,
            b_panel=1000,
            position_panel="Extérieur",
            type_panel="OSB/3 11-18 mm"
        )
        
        # Appeler la méthode à tester
        df_kp = mob_instance._k_panel()
        df_ksw, df_kw = mob_instance._k_wall()
        print(df_kp)
        print(df_ksw)
        print(df_kw)
        
        # Vérifications
        # 1. Vérifier que les DataFrames ne sont pas vides
        assert not df_ksw.empty
        assert not df_kw.empty
        
        # 2. Vérifier les colonnes des DataFrames
        expected_ksw_columns = ['Kser système de mur', 'Coeff. équivalence']
        expected_kw_columns = ['Kser mur', 'Coeff. équivalence']
        
        assert all(col in df_ksw.columns for col in expected_ksw_columns)
        assert all(col in df_kw.columns for col in expected_kw_columns)
        
        # 3. Vérifier que les valeurs de Kser sont positives
        assert (df_ksw['Kser système de mur'] > 0).all()
        assert (df_kw['Kser mur'] > 0).all()
        
        # 4. Vérifier que la somme des Kser des murs est égale au Kser du système
        kser_system = df_ksw.loc[sys_name, 'Kser système de mur']
        kser_sum = df_kw['Kser mur'].sum()
        assert kser_system == kser_sum
        