#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import pytest
import forallpeople as si
si.environment("structural")
# Import des modules à tester
from ourocode.eurocode import EC8_Sismique as EC8
from ourocode.eurocode import A0_Projet as A0

# Fixture pour un projet de base
@pytest.fixture
def project():
    return A0.Projet(
        name="Projet Test",
        code_INSEE="73011",
        pays="France",
        alt=400,
    )

@pytest.fixture
def batiment(project):
    return A0.Batiment._from_parent_class(
        project,
        h_bat=10,
        d_bat=30,
        b_bat=20,
        alpha_toit=15
    )

class Test_Sismique:
    @pytest.fixture
    def sismique(self, batiment):
        sismique = EC8.Sismique._from_parent_class(
            batiment,
            Kbx=15000,
            Kby=10000,
            type_constructif_x = "Panneaux de murs avec diaphragmes cloués, assemblés par clous et boulons", 
            type_constructif_y = "Panneaux de murs avec diaphragmes cloués, assemblés par clous et boulons",
            regulier_plan = True,
            regulier_elevation = True,
            cat_importance="II",
            classe_sol="C",
            type_spectre="Type 2"
        )
        sismique.add_gravity_load(name="LRDC1",load=4.2, surface=330, etage="RDC", z_i=3.5, action="Permanente G", categorie_Q="Aucune", occupations="Autres", comment="Charge plancher G")
        sismique.add_gravity_load(name="LRDC2",load=1.05, surface=296.1, etage="RDC", z_i=3.5, action="Permanente G", categorie_Q="Aucune", occupations="Autres", comment="Charge mur G")
        sismique.add_gravity_load(name="LToiture1",load=2.5, surface=120, etage="Toiture", z_i=7, action="Permanente G", categorie_Q="Aucune", occupations="Autres", comment="Charge toiture G")
        sismique.add_gravity_load(name="LToiture2",load=1.05, surface=154, etage="Toiture", z_i=7, action="Permanente G", categorie_Q="Aucune", occupations="Autres", comment="Charge mur toiture G")
        sismique.add_gravity_load(name="LRDC3",load=3, surface=330, etage="RDC", z_i=3.5, action="Exploitation Q", categorie_Q="Cat C : lieu de réunion", occupations="Étages à occupations indépendantes", comment="Charge plancher Q")
        return sismique

    def test_init(self, sismique):
        assert sismique.region_sismique == "Zone 4"
        assert sismique.gamma_1 == 1
        assert sismique.h_bat == 10
        assert sismique.d_bat == 30
        assert sismique.b_bat == 20
        assert sismique.alpha_toit == 15
        assert sismique.type_constructif == {"x": "Panneaux de murs avec diaphragmes cloués, assemblés par clous et boulons", "y": "Panneaux de murs avec diaphragmes cloués, assemblés par clous et boulons"}
        assert sismique.coeff_comportement == {"x": {'Classe de ductilité': 'DCH', 'q': 3}, "y": {'Classe de ductilité': 'DCH', 'q': 3}}
        assert sismique.K_b == {"Raideur du bâtiment dans la direction x": 15000*si.kN / si.m, "Raideur du bâtiment dans la direction y": 10000*si.kN / si.m}
        assert sismique.regulier_plan == True
        assert sismique.regulier_elevation == True
        assert sismique.cat_importance == "II"
        assert sismique.classe_sol == "C"
        assert sismique.type_spectre == "Type 2"
    
    @pytest.mark.parametrize(
        "cat_importance, expected_output",
        [
            ("II", True),  # Catégories II-IV nécessitent une analyse
            ("III", True),
            ("IV", True),
        ],
    )
    def test_has_to_be_analyzed(self, sismique, cat_importance, expected_output):
        """Teste si une analyse sismique est nécessaire selon la catégorie d'importance"""
        sismique.cat_importance = cat_importance
        assert sismique._has_to_be_analyzed() == expected_output

    def test_add_gravity_load(self, sismique):
        """Teste l'ajout de charges gravitaires"""
        
        # Vérification que la charge a bien été ajoutée
        assert "RDC" in sismique.gravity_loads
        assert "LRDC1" in sismique.gravity_loads["RDC"]
        assert "LRDC2" in sismique.gravity_loads["RDC"]
        assert "LToiture1" in sismique.gravity_loads["Toiture"]
        assert "LToiture2" in sismique.gravity_loads["Toiture"]
        
        # Vérification des valeurs
        charge = sismique.gravity_loads["RDC"]["LRDC1"]
        assert charge["Charge gravitaire"] == 4.2 * si.kN / si.m**2
        assert charge["Surface"] == 330 * si.m**2
        assert charge["Zi"] == 3.5 * si.m
        assert charge["Action"] == "Permanente G"

    def test_periode_fondamentale(self, sismique):
        """Teste le calcul de la période fondamentale"""
        T1 = sismique.T1[1]
        assert T1["x"] == 0.2558960726545056
        assert T1["y"] == 0.31340740259285516

    def test_Fb(self, sismique):
        """Teste le calcul de l'effort tranchant à la base"""
        Fb = sismique.Fb[1]
        assert Fb["x"].value == 500633.0275229359
        assert Fb["y"].value == 500633.0275229359

    def test_Fi(self, sismique):
        """Teste le calcul de l'effort tranchant à la base"""
        Fi = sismique.Fi("Toiture")[1]
        assert Fi["x"].value == 158462.87502152807
        assert Fi["y"].value == 158462.87502152807

    def test_coeff_torsion_accidentelle(self, sismique):
        """Teste le calcul du coefficient de torsion accidentelle"""
        delta = sismique.coeff_torsion_accidentelle(15, 30)[1]
        assert delta == 1.6

    def test_effet_second_ordre(self, sismique):
        """Teste le calcul de l'effet du second ordre"""
        p_delta = sismique.coeff_second_ordre(dr=80, V_tot=342.17, etage="RDC")[1]
        assert p_delta == 1.0

    def test_effort_sismique_final_capacite(self, sismique):
        """Teste le calcul de l'effort sismique final"""
        F_final = sismique.F_sismique_final_capacite(
            etage="RDC",
            gamma_d_x="Rupture fragile",
            gamma_d_y="Rupture ductile",
            Omega_x=1.5,
            Omega_y=1.3,
            eta_torsion_x=1.2,
            eta_torsion_y=1.1,
            P_delta_x=1.05,
            P_delta_y=1.05
        )[1]
        assert F_final["x"].value == 840712.0646959591
        assert F_final["y"].value == 565145.3323789504

    def test_effort_sismique_final_dissipatif(self, sismique):
        """Teste le calcul de l'effort sismique final"""
        F_final = sismique.F_sismique_final_dissipatif(
            etage="RDC",
            eta_torsion_x=1.2,
            eta_torsion_y=1.1,
            P_delta_x=1.05,
            P_delta_y=1.05
        )[1]
        assert F_final["x"].value == 431134.39215177385
        assert F_final["y"].value == 395206.5261391261

    def test_effort_non_structural(self, sismique):
        """Teste le calcul de l'effort sismique pour les éléments non structuraux"""
        Fa = sismique.Fa(150, 0.2, 0.15, 10, "Garde-corps ou ornements", "Autres")[1]
        assert Fa["x"].value == 191.80770522277845
        assert Fa["y"].value == 154.77238988189353

    def test_limitation_dommages(self, sismique):
        """Teste le calcul du taux de limitation des dommages"""
        tx = sismique.taux_limitations_dommages(dr=80, etage="RDC", type_dommages="Ductiles")[1]
        assert tx == 0.42560000000000003