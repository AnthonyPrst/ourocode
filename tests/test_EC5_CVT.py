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
from ourocode.eurocode.ec5 import cvt as EC5
from ourocode.eurocode.core.batiment import Batiment
from ourocode.eurocode.ec5.element_droit.barre import Barre
from ourocode.eurocode.ec5.assemblage import Pointe


@pytest.fixture
def batiment():
    """Fixture pour créer une instance de Batiment pour les tests."""
    return Batiment(h_bat=3, d_bat=10.0, b_bat=20.0, alpha_toit=10)


@pytest.fixture
def connectors():
    """Crée des connecteurs (pointes) et panneaux/montants compatibles avec MOB."""
    # Panneau OSB intérieur/extérieur
    panel = Barre(b=15, h=1196, section="Rectangulaire", classe="OSB/3 11-18 mm", cs=1, effet_systeme=False)
    montant = Barre(b=60, h=120, section="Rectangulaire", classe="C24", cs=1, effet_systeme=False)
    # Connecteurs par pointes (paramètres raisonnables)
    conn_int = Pointe(
        d=2.5,
        dh=5,
        l=55,
        qualite="6.8",
        n=1,
        type_organe="Pointe circulaire lisse",
        percage=False,
        beam_1=panel,
        beam_2=montant,
        nfile=1,
        nCis=1,
    )
    conn_ext = Pointe(
        d=2.5,
        dh=5,
        l=55,
        qualite="6.8",
        n=1,
        type_organe="Pointe circulaire lisse",
        percage=False,
        beam_1=panel,
        beam_2=montant,
        nfile=1,
        nCis=1,
    )
    return conn_int, conn_ext


class TestMOB:
    """Classe de test pour la classe MOB (Mur Ossature Bois)."""

    @pytest.fixture
    def mob_instance(self, batiment, connectors):
        """Fixture pour créer une instance de MOB avec connecteurs pour les tests."""
        conn_int, conn_ext = connectors
        mob = EC5.MOB._from_parent_class(batiment, connecteurs_int=conn_int, connecteurs_ext=conn_ext)
        return mob

    def test_add_wall(self, mob_instance):
        """Test l'ajout d'un système de mur."""
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(
            h_etage=3.0,
            h_sys_MOB=2.7,
            l_sys_MOB=10.0,
            etage=etage,
        )

        assert sys_name == "S1"
        assert etage in mob_instance.data
        assert "S1" in mob_instance.data[etage]
        assert mob_instance.data[etage]["S1"]["Etage"] == "RDC"
        assert mob_instance.data[etage]["S1"]["Hauteur étage"] == 3.0 * si.m
        assert mob_instance.data[etage]["S1"]["Hauteur du système de mur"] == 2.7 * si.m
        assert mob_instance.data[etage]["S1"]["Longueur du système de mur"] == 10.0 * si.m
        assert mob_instance.data[etage]["S1"]["Murs"] == {}

    def test_add_internal_wall(self, mob_instance):
        """Test l'ajout d'un mur interne."""
        # D'abord ajouter un système de mur
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=3.0, h_sys_MOB=2.7, l_sys_MOB=10.0, etage=etage)

        # Ajouter un mur interne
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )

        assert len(mob_instance.data[etage][sys_name]["Murs"]) == 1
        wall = mob_instance.data[etage][sys_name]["Murs"][wall_name]
        assert wall_name == "S1_W1"
        assert wall["Position en X"] == 1.5 * si.m
        assert wall["Longueur"] == 3.0 * si.m
        assert wall["Ff,Rk ext"] == mob_instance._data_connecteurs["Extérieur"]["Ff,Rk"]
        assert wall["Ff,Rk int"] == mob_instance._data_connecteurs["Intérieur"]["Ff,Rk"]
        assert wall["Kser ext"] == mob_instance._data_connecteurs["Extérieur"]["Kser"]
        assert wall["Kser int"] == mob_instance._data_connecteurs["Intérieur"]["Kser"]
        assert wall["Couturage ext. rive"] == 150 * si.mm
        assert wall["Couturage ext. intermédiaire"] == 300 * si.mm
        assert wall["Couturage int. rive"] == 150 * si.mm
        assert wall["Couturage int. intermédiaire"] == 300 * si.mm
        assert wall["Panneaux"] == {}

    def test_add_panel_to_wall(self, mob_instance):
        """Test l'ajout d'un panneau à un mur."""
        # Ajouter un système de mur
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=3.0, h_sys_MOB=2.7, l_sys_MOB=10.0, etage=etage)

        # Ajouter un mur interne
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )

        # Ajouter un panneau
        panel_name = mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            h_panel=2700,
            b_panel=1250,
            number=2,
            position_panel="Extérieur",
        )

        # Vérifier que le panneau a été correctement ajouté
        wall = mob_instance.data[etage][sys_name]["Murs"][wall_name]
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
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=3.0, h_sys_MOB=2.7, l_sys_MOB=10.0, etage=etage)

        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )

        # Capture des valeurs attendues strictes avant sauvegarde
        expected_kser_ext = mob_instance._data_connecteurs["Extérieur"]["Kser"]
        expected_epaisseur = mob_instance._connecteurs["Intérieur"].beam_1.b
        expected_type = mob_instance._connecteurs["Intérieur"].beam_1.classe

        panel_name = mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=2,
            h_panel=2700,
            b_panel=1250,
            position_panel="Extérieur",
        )

        # Sauvegarder dans un fichier temporaire
        test_file = tmp_path / "test_walls.json"
        mob_instance.save_walls_data(test_file)

        # Vérifier que le fichier a été créé
        assert test_file.exists()

        # Créer une nouvelle instance et charger les données
        new_mob = EC5.MOB._from_parent_class(batiment)
        new_mob.load_walls_data(test_file)

        # Vérifier que les données ont été correctement chargées
        assert etage in new_mob.data
        assert sys_name in new_mob.data[etage]
        assert len(new_mob.data[etage][sys_name]["Murs"]) == 1
        assert len(new_mob.data[etage][sys_name]["Murs"][wall_name]["Panneaux"]) == 1
        assert new_mob.data[etage][sys_name]["Murs"][wall_name]["Kser ext"] == expected_kser_ext
        assert new_mob.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Nombre"] == 2
        assert (
            new_mob.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Epaisseur"] == expected_epaisseur
        )
        assert new_mob.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Hauteur"] == 2700 * si.mm
        assert new_mob.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Largeur"] == 1250 * si.mm
        assert new_mob.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]["Type"] == expected_type

    def test_K_panel(self, mob_instance):
        """Test le calcul de la raideur des panneaux."""
        # Ajouter un système de mur
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=2.8, h_sys_MOB=2.5, l_sys_MOB=10.0, etage=etage)

        # Ajouter un mur avec des propriétés de couturage
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )

        # Ajouter des panneaux de test
        mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=2,
            h_panel=2500,
            b_panel=1196,
            position_panel="Extérieur",
        )

        mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=2500,
            b_panel=1070,
            position_panel="Intérieur",
        )

        mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=2500,
            b_panel=1070,
            position_panel="Extérieur",
        )

        # Appeler la méthode à tester
        df_kp = mob_instance._k_panel()

        # Vérifications
        # 1. Vérifier que les DataFrames ne sont pas vides
        # print(df_kp)
        assert not df_kp.empty

        # 2. Vérifier les colonnes du DataFrame des raideurs
        expected_columns = ["Hauteur", "Type", "Kser,p unit.", "Coeff. équivalence"]
        assert all(col in df_kp.columns for col in expected_columns)

        # 3. Vérifier que les raideurs sont positives
        assert (df_kp["Kser,p unit."] > 0).all()

        # 4. Vérifier que le coefficient d'équivalence est compris entre 0 et 1
        assert (df_kp["Coeff. équivalence"] > 0).all()
        assert (df_kp["Coeff. équivalence"] <= 1.0).all()

        # 5. Vérifier que le panneau le plus large a un coefficient d'équivalence de 1.0
        max_width_idx = df_kp["Largeur"].idxmax()
        assert np.isclose(df_kp.loc[max_width_idx, "Coeff. équivalence"], 1.0)

    def _setup_wall_with_panels(self, mob_instance):
        """Helper : crée un système de mur avec 3 panneaux et retourne (sys_name, wall_name)."""
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=2.8, h_sys_MOB=2.5, l_sys_MOB=10.0, etage=etage)
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.5,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )
        mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=2,
            h_panel=2500,
            b_panel=1196,
            position_panel="Extérieur",
        )
        mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=2500,
            b_panel=1070,
            position_panel="Intérieur",
        )
        mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=2500,
            b_panel=1070,
            position_panel="Extérieur",
        )
        return sys_name, wall_name

    def test_k_wall(self, mob_instance):
        """Test le calcul de la raideur des murs et systèmes de murs."""
        # Configurer les données de test via le helper
        sys_name, wall_name = self._setup_wall_with_panels(mob_instance)

        # Ajouter un deuxième mur avec des panneaux différents
        wall_name2 = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=4.5,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )

        # Ajouter des panneaux au deuxième mur
        mob_instance.add_panel_to_wall(
            wall_name=wall_name2,
            number=3,
            h_panel=2500,
            b_panel=1000,
            position_panel="Extérieur",
        )

        # Appeler la méthode à tester
        mob_instance._k_panel()
        # Calcul des raideurs par murs et par systèmes via l'API publique
        mob_instance.get_Kser_wall()
        df_ksw = mob_instance.data_Kser_sys_walls
        df_kw = mob_instance.data_Kser_walls
        # print(df_kp)
        # print(df_ksw)
        # print(df_kw)

        # Vérifications
        # 1. Vérifier que les DataFrames ne sont pas vides
        assert not df_ksw.empty
        assert not df_kw.empty

        # 2. Vérifier les colonnes des DataFrames
        expected_ksw_columns = ["Kser système de mur", "Coeff. équivalence"]
        expected_kw_columns = ["Kser mur", "Coeff. équivalence"]

        assert all(col in df_ksw.columns for col in expected_ksw_columns)
        assert all(col in df_kw.columns for col in expected_kw_columns)

        # 3. Vérifier que les valeurs de Kser sont positives
        assert (df_ksw["Kser système de mur"] > 0).all()
        assert (df_kw["Kser mur"] > 0).all()

        # 4. Vérifier que la somme des Kser des murs est égale au Kser du système
        kser_system = df_ksw.loc[sys_name, "Kser système de mur"]
        kser_sum = df_kw["Kser mur"].sum()
        assert kser_system == kser_sum

    def test_taux_ancrage(self, mob_instance):
        """Test le calcul des taux d'ancrage des murs (soulèvement et cisaillement)."""
        # Créer un système de mur avec charges directes pour remplir data_walls_loads via taux_walls
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(
            h_etage=3.0,
            h_sys_MOB=2.7,
            l_sys_MOB=10.0,
            etage=etage,
            Fv_Ed_ELU=80.0,  # kN
            Fv_Ed_ELS=30.0,  # kN
        )

        # Ajouter deux murs et des panneaux contreventants des deux côtés
        wall_name1 = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.0,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )
        mob_instance.add_panel_to_wall(
            wall_name=wall_name1,
            number=2,
            h_panel=2500,
            b_panel=1200,
            position_panel="Extérieur",
        )
        mob_instance.add_panel_to_wall(
            wall_name=wall_name1,
            number=2,
            h_panel=2500,
            b_panel=1200,
            position_panel="Intérieur",
        )

        wall_name2 = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=5.0,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )
        mob_instance.add_panel_to_wall(
            wall_name=wall_name2,
            number=1,
            h_panel=2500,
            b_panel=1000,
            position_panel="Extérieur",
        )

        # Calculer les taux de murs pour remplir data_walls_loads (via charges directes)
        mob_instance.taux_walls()

        # Définir les barres pour ancrages
        montant = Barre(b=45, h=120, section="Rectangulaire", classe="C24", cs=1, effet_systeme=False)
        lisse_basse = Barre(b=60, h=120, section="Rectangulaire", classe="C24", cs=1, effet_systeme=False)
        lisse_impl = Barre(b=80, h=120, section="Rectangulaire", classe="C24", cs=1, effet_systeme=False)

        # Lancer le calcul des taux d'ancrage avec valeurs raisonnables
        df_taux_anc = mob_instance.taux_ancrage(
            montant=montant,
            lisse_basse=lisse_basse,
            lisse_impl=lisse_impl,
            Ft_Rk_wood=20.0,  # kN (équerre/bois)
            Ft_Rd_concrete=30.0,  # kN (fixation au béton)
            Fv_Rk_anc_lb=6.0,  # kN par ancrage lisse basse
            e_anc_lb=600,  # mm
            Fv_Rk_anc_li=8.0,  # kN par ancrage lisse implantation
            e_anc_li=800,  # mm
        )

        # Vérifications
        print(df_taux_anc)
        assert not df_taux_anc.empty
        # Les deux murs ajoutés doivent être présents
        assert wall_name1 in df_taux_anc.index
        assert wall_name2 in df_taux_anc.index

        expected_columns = [
            "Étage",
            "Longueur mur",
            "Ft,Ed",
            "Ft,Rd",
            "Taux ancr. soulèvement (%)",
            "Fv,anc,Ed",
            "Fv,anc,lb,Rd",
            "Nbr ancr,lb",
            "entraxe anc,lb",
            "Taux ancr. lisse basse (%)",
            "Fv,anc,li,Rd",
            "Nbr ancr,li",
            "entraxe anc,li",
            "Taux ancr. lisse impl. (%)",
        ]
        print(df_taux_anc.columns)
        assert all(col in df_taux_anc.columns for col in expected_columns)

        # Le nombre d'ancrages doit correspondre à floor(longueur / entraxe)
        long_w1 = mob_instance.data_walls[wall_name1]["Longueur"]
        nb_lb_w1 = df_taux_anc.loc[wall_name1, "Nbr ancr,lb"]
        nb_li_w1 = df_taux_anc.loc[wall_name1, "Nbr ancr,li"]
        assert nb_lb_w1 == int(np.floor(long_w1 / (600 * si.mm)))
        assert nb_li_w1 == int(np.floor(long_w1 / (800 * si.mm)))

        # Les entraxes doivent être stockés avec unités mm
        assert df_taux_anc.loc[wall_name1, "entraxe anc,lb"] == 600 * si.mm
        assert df_taux_anc.loc[wall_name1, "entraxe anc,li"] == 800 * si.mm

        # Les taux doivent être des valeurs positives
        assert (df_taux_anc["Taux ancr. soulèvement (%)"] >= 0).all()
        assert (df_taux_anc["Taux ancr. lisse basse (%)"] >= 0).all()
        assert (df_taux_anc["Taux ancr. lisse impl. (%)"] >= 0).all()

    def test_non_contreventant_panel_excluded(self, mob_instance):
        """Un panneau qui ne respecte pas le rapport d'aspect ne doit pas être pris en compte dans la raideur."""
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=3.0, h_sys_MOB=2.7, l_sys_MOB=5.0, etage=etage)
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=0.5,
            l_MOB=2.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=0,
            couturage_int_inter=0,
        )
        # Panneau non contreventant: h/4 > b -> 2500/4 = 625 > 500
        panel_name = mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=2500,
            b_panel=500,
            position_panel="Extérieur",
        )
        # Flag dans les données brutes
        etage_data = mob_instance.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]
        assert etage_data["Contreventement"] is False
        # Le calcul des K doit ignorer ce panneau (raideur de mur nulle ici faute d'autres panneaux)
        mob_instance.get_Kser_wall()
        k_wall = mob_instance.data_Kser_walls.loc[wall_name, "Kser mur"]
        assert k_wall == 0

    def test_stacked_panels_processed(self, mob_instance):
        """Les panneaux empilés (on_top_of) sont traités via la raideur équivalente sans double comptage."""
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=3.0, h_sys_MOB=2.7, l_sys_MOB=6.0, etage=etage)
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=1.0,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=150,
            couturage_int_inter=300,
        )
        bottom = mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=1250,
            b_panel=1196,
            position_panel="Extérieur",
        )
        top = mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=1,
            h_panel=1250,
            b_panel=1196,
            position_panel="Extérieur",
            on_top_of=bottom,
        )
        # Le calcul ne doit pas lever d'erreur et la raideur du mur doit être > 0
        mob_instance.get_Kser_wall()
        k_wall = mob_instance.data_Kser_walls.loc[wall_name, "Kser mur"]
        assert k_wall > 0
        # Les 2 panneaux existent dans les données des panneaux
        assert bottom in mob_instance.data_Kser_panels.index
        assert top in mob_instance.data_Kser_panels.index

    def test_single_sided_panel(self, mob_instance):
        """Configuration avec panneau sur une seule face: pas de panneau double et raideur > 0."""
        etage = "RDC"
        sys_name = mob_instance.add_sys_wall(h_etage=3.0, h_sys_MOB=2.7, l_sys_MOB=6.0, etage=etage)
        wall_name = mob_instance.add_wall_to_sys(
            sys_name=sys_name,
            position=0.8,
            l_MOB=3.0,
            couturage_ext_rive=150,
            couturage_ext_inter=300,
            couturage_int_rive=0,
            couturage_int_inter=0,
        )
        panel_name = mob_instance.add_panel_to_wall(
            wall_name=wall_name,
            number=2,
            h_panel=2500,
            b_panel=1196,
            position_panel="Extérieur",
        )
        # Vérifications d'état
        p = mob_instance.data[etage][sys_name]["Murs"][wall_name]["Panneaux"][panel_name]
        assert p["Panneau double"] == "Non"
        # La raideur doit se calculer sans la face intérieure
        mob_instance.get_Kser_wall()
        k_wall = mob_instance.data_Kser_walls.loc[wall_name, "Kser mur"]
        assert k_wall > 0


class TestMOBDiaphragme:
    """Répartition diaphragme rigide, murs béton et visualisation."""

    @pytest.fixture
    def mob(self, batiment, connectors):
        conn_int, conn_ext = connectors
        return EC5.MOB._from_parent_class(batiment, connecteurs_int=conn_int, connecteurs_ext=conn_ext)

    @staticmethod
    def _add_sys(mob, Xg, Yg, alpha, l_sys=10.0, etage="RDC"):
        sys_name = mob.add_sys_wall(
            h_etage=2.8, h_sys_MOB=2.5, l_sys_MOB=l_sys, etage=etage, Xg=Xg, Yg=Yg, alpha_sw=alpha
        )
        wall = mob.add_wall_to_sys(
            sys_name=sys_name, position=1.0, l_MOB=2.4, couturage_ext_rive=150, couturage_ext_inter=300
        )
        mob.add_panel_to_wall(wall_name=wall, number=2, h_panel=2500, b_panel=1196, position_panel="Extérieur")
        return sys_name

    def _box(self, mob):
        """Bâtiment rectangulaire 20 x 10 m, un système de mur par façade."""
        return [
            self._add_sys(mob, 10, 0, 0, 20),
            self._add_sys(mob, 10, 10, 0, 20),
            self._add_sys(mob, 0, 5, 90),
            self._add_sys(mob, 20, 5, 90),
        ]

    def test_solver_equilibre(self):
        elements = [(0, 0, 0, 2e6), (3, 8, 0, 1e6), (0, 4, np.pi / 2, 3e6), (12, 2, np.pi / 2, 1e6)]
        Fx, Fy, XF, YF = 1e4, -4e3, 7.0, 3.0
        res = EC5.MOB._solve_diaphragm(elements, Fx, Fy, XF, YF)
        x, y, al, _ = np.array(elements).T
        F = res["F"]
        assert np.isclose((F * np.cos(al)).sum(), Fx)
        assert np.isclose((F * np.sin(al)).sum(), Fy)
        # Moment autour de l'origine
        assert np.isclose((F * (x * np.sin(al) - y * np.cos(al))).sum(), XF * Fy - YF * Fx)
        # Force appliquée au centre de raideur : pas de rotation
        res_cr = EC5.MOB._solve_diaphragm(elements, Fx, Fy, res["xO"], res["yO"])
        assert abs(res_cr["theta"]) < 1e-12

    def test_solver_instable(self):
        with pytest.raises(ValueError):
            EC5.MOB._solve_diaphragm([(0, 0, 0, 1e6), (0, 5, 0, 1e6), (0, 9, 0, 1e6)], 1e3, 0, 0, 0)

    def test_repartition_symetrique(self, mob):
        s1, s2, s3, s4 = self._box(mob)
        df = mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=15, XFw=10, YFw=5, alpha_Fw=0, etage="RDC")
        F = {s: df.loc[s, "Fv,Ed,i ELU"].value for s in (s1, s2, s3, s4)}
        assert np.isclose(F[s1], 7.5e3) and np.isclose(F[s2], 7.5e3)
        assert abs(F[s3]) < 1e-6 and abs(F[s4]) < 1e-6
        tors = mob.synthese_torsion().loc["RDC"]
        assert np.isclose(tors["xO"].value, 10) and np.isclose(tors["yO"].value, 5)
        # Les murs et panneaux reçoivent l'effort du système
        assert np.isclose(mob.data_walls_loads.loc[f"{s1}_W1", "Fv,Ed,i ELU"].value, 7.5e3)

    def test_excentricite_et_e_acc(self, mob):
        s1, s2, *_ = self._box(mob)
        df = mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=15, XFw=10, YFw=8, etage="RDC")
        F1, F2 = df.loc[s1, "Fv,Ed,i ELU"].value, df.loc[s2, "Fv,Ed,i ELU"].value
        assert F2 > F1 and np.isclose(F1 + F2, 15e3)
        df_acc = mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=15, XFw=10, YFw=8, etage="RDC", e_acc=0.5)
        assert df_acc.loc[s2, "Fv,Ed,i ELU"].value > F2

    def test_mur_beton(self, mob):
        s1, s2, s3, s4 = self._box(mob)
        b1 = mob.add_concrete_wall(l_wall=5, t_wall=0.2, h_wall=2.8, X=5, Y=2, alpha_w=0, etage="RDC")
        wall = mob.data_concrete_walls[b1]
        E = 0.5 * 31000e6
        I, A = 0.2 * 5**3 / 12, 0.2 * 5
        K = 1 / (2.8**3 / (3 * E * I) + 1.2 * 2.8 / (E / 2.4 * A))
        assert np.isclose(wall["Kser"].value, K)
        mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=15, XFw=10, YFw=5, etage="RDC")
        F_b = mob.data_concrete_loads.loc[b1, "Fv,Ed,i ELU"].value
        assert F_b > 0.9 * 15e3
        assert np.isclose(mob.data_concrete_loads.loc[b1, "M,base,Ed ELU"].value, F_b * 2.8)
        F_mob = sum(mob.data_sys_walls_loads.loc[s, "Fv,Ed,i ELU"].value for s in (s1, s2))
        assert np.isclose(F_b + F_mob, 15e3)
        # Le centre de raideur se déplace vers le refend
        assert 2 <= mob.synthese_torsion().loc["RDC", "yO"].value < 2.1
        assert mob.add_concrete_wall(l_wall=1, t_wall=0.2, h_wall=2.8, K_override=100, etage="RDC") == "B2"
        assert np.isclose(mob.data_concrete_walls["B2"]["Kser"].value, 1e8)

    def test_ku_elu(self, mob):
        """ELU avec Ku = 2/3 Kser pour les MOB : le refend béton reprend une part plus importante qu'à l'ELS."""
        s1, s2, s3, s4 = self._box(mob)
        b1 = mob.add_concrete_wall(l_wall=1, t_wall=0.2, h_wall=2.8, K_override=2, X=5, Y=2, etage="RDC")
        mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=10, XFw=10, YFw=5, etage="RDC")
        K_b = mob.data_concrete_walls[b1]["Kser"].value
        part_ELS = K_b * mob.data_concrete_loads.loc[b1, "delta,v,Ed,i ELS"].value / 10e3
        part_ELU = mob.data_concrete_loads.loc[b1, "Fv,Ed,i ELU"].value / 10e3
        assert part_ELU > part_ELS
        F_mob = sum(mob.data_sys_walls_loads.loc[s, "Fv,Ed,i ELU"].value for s in (s1, s2))
        assert np.isclose(part_ELU * 10e3 + F_mob, 10e3)
        tors = mob.synthese_torsion().loc["RDC"]
        assert tors["yO ELU"].value < tors["yO"].value

    def test_save_load_beton(self, mob, batiment, tmp_path):
        self._box(mob)
        mob.add_concrete_wall(l_wall=5, t_wall=0.2, h_wall=2.8, X=5, Y=5, etage="RDC")
        path = str(tmp_path / "mob.json")
        mob.save_walls_data(path)
        new = EC5.MOB._from_parent_class(batiment)
        new.load_walls_data(path)
        assert list(new.data["RDC"].keys()) == ["S1", "S2", "S3", "S4"]
        assert "B1" in new.data_concrete_walls

    def test_show_plan(self, mob, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        self._box(mob)
        mob.add_concrete_wall(l_wall=5, t_wall=0.2, h_wall=2.8, X=5, Y=5, etage="RDC")
        mob.calculate_loads_and_deformations(Fw_Ed_ELS=10, Fw_Ed_ELU=15, XFw=10, YFw=7, alpha_Fw=30, etage="RDC")
        mob.taux_walls()
        path = tmp_path / "plan.png"
        fig, ax = mob.show_plan(filepath=str(path), show=False)
        assert path.exists()
        assert ax.get_title() == "Contreventement - RDC"
