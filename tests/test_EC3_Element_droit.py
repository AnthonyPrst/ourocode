"""Tests unitaires pour le module EC3_Element_droit.py."""
import pytest
import math
import sys
import os

# Ajout du répertoire parent au chemin Python pour les imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import forallpeople as si
si.environment("structural")
from ourocode.eurocode.ec3.element_droit import Plat, Traction, Compression, Cisaillement, Flexion

class TestPlat:
    """Tests pour la classe de base Plat."""
    
    def test_initialisation(self):
        """Test l'initialisation de la classe Plat."""
        element = Plat(t=10, h=100, b=200, classe_acier="S235", classe_transv=1)
        assert element.t == 10.0 * si.mm
        assert element.h == 100.0 * si.mm
        assert element.classe_acier == "S235"
        assert element.classe_transv == 1
        
    def test_fy_fu_values(self):
        """Test les valeurs de fy et fu pour différentes épaisseurs."""
        # Test avec t <= 40mm
        element1 = Plat(t=20, h=100, b=200, classe_acier="S235", classe_transv=1)
        assert element1.fy == 235000000.0  # 235 MPa en Pa
        assert element1.fu == 360000000.0  # 360 MPa en Pa
        
        # Test avec 40mm < t <= 80mm
        element2 = Plat(t=50, h=100, b=200, classe_acier="S235", classe_transv=1)
        assert element2.fy == 215000000.0  # 215 MPa en Pa
        assert element2.fu == 360000000.0  # 360 MPa en Pa


class TestTraction:
    """Tests pour la classe Traction."""
    
    def test_Npl_Rd(self):
        """Test le calcul de la résistance plastique en traction."""
        traction = Traction(A=1000, Anet=0, ass_cat_C=False, t=10, b=200, h=100, classe_acier="S235", classe_transv=1)
        assert traction._Npl_Rd[1] == 235 * si.kN
    
    def test_Nu_Rd(self):
        """Test le calcul de la résistance ultime en traction."""
        traction = Traction(A=1000, Anet=800, t=10, b=200, h=100, ass_cat_C=False, classe_acier="S235", classe_transv=1)
        assert traction._Nu_Rd[1] == 207.36 * si.kN
    
    def test_Nt_Rd_assemblage_cat_C(self):
        """Test le calcul de Nt,Rd pour un assemblage de catégorie C."""
        traction = Traction(A=1000, Anet=800, ass_cat_C=True, t=10, b=200, h=100, classe_acier="S235", classe_transv=1)
        assert traction.Nt_Rd[1] == 188 * si.kN


class TestCompression:
    """Tests pour la classe Compression."""
    
    def test_Nc_Rd(self):
        """Test le calcul de la résistance en compression."""
        compression = Compression(A=1000, Iy=8333333.33, Iz=8333.33, lo_y=2000, lo_z=3000, type_appuis="Encastré 1 côté", t=10, b=200, h=100, classe_acier="S235", classe_transv=3)
        assert compression.Nc_Rd[1] == 235 * si.kN
    
    def test_lamb(self):
        """Test le calcul de l'élancement."""
        compression = Compression(A=1000, Iy=833333.33, Iz=8333.33, lo_y=2000, lo_z=3000, type_appuis="Encastré 1 côté", t=10, b=200, h=100, classe_acier="S235", classe_transv=3)
        # Récupération des valeurs calculées
        lamb_values = compression.lamb
        
        # Vérifications de base
        assert isinstance(lamb_values, dict)
        assert 'y' in lamb_values
        assert 'z' in lamb_values
        
        # Vérification que les valeurs sont positives
        assert lamb_values['y'] > 0
        assert lamb_values['z'] > 0
        assert lamb_values['y'] == 138.5640648826383
        assert lamb_values['z'] == 2078.461384774971


class TestCisaillement:
    """Tests pour la classe Cisaillement."""
    
    def test_Vpl_Rd(self):
        """Test le calcul de la résistance au cisaillement."""
        # Création d'une instance avec des valeurs de test
        cisaillement = Cisaillement(Av=500, t=10, h=100, b=200, classe_acier="S235", classe_transv=3)
        assert cisaillement.Vpl_Rd[1].value == 67838.65662978104

class TestFlexion:
    """Tests pour la classe Flexion."""
    
    def test_Mc_Rd_classe_1_2(self):
        """Test le calcul du moment résistant pour les sections de classe 1 et 2."""
        # Création d'une instance avec des valeurs de test (classe 1)
        flexion = Flexion(W=20000, t=10, h=100, b=200, classe_acier="S235", classe_transv=3)
        assert flexion.Mc_Rd[1] == 4.7 * si.kN * si.m
    
    def test_Mc_V_Rd(self):
        """Test le calcul du moment résistant avec prise en compte du cisaillement."""
        # Création d'une instance avec des valeurs de test
        flexion = Flexion(W=20000, t=10, h=100, b=200, classe_acier="S235", classe_transv=3)
        # Appel de la méthode à tester
        result = flexion.Mc_V_Rd(Av=500, V_Ed=40)
        print(result[0])
        assert result[1].value == 4548.954955674857


class TestCompressionFlambement:
    """Tests pour les branches de flambement de Compression."""

    def test_Nb_Rd(self):
        """Test la résistance au flambement Nb,Rd (branches chi y et z)."""
        comp = Compression(
            A=3000, Iy=6750000, Iz=562500,
            lo_y=3000, lo_z=3000,
            courbe_flamb={"y": "b", "z": "c"},
            type_appuis="Rotule - Rotule",
            t=15, h=150, b=100,
            classe_acier="S235", classe_transv=1,
        )
        latex, val = comp.Nb_Rd
        assert "y" in val and "z" in val
        assert val["y"].value > 0
        assert val["z"].value > 0
        assert val["y"].value > val["z"].value  # courbe b moins pénalisante que c

    def test_lamb_rel_et_chi_positifs(self):
        """Test que _lamb_rel_Axe et _chi sont strictement positifs."""
        comp = Compression(
            A=3000, Iy=6750000, Iz=562500,
            lo_y=3000, lo_z=3000,
            courbe_flamb={"y": "a", "z": "b"},
            type_appuis="Rotule - Rotule",
            t=15, h=150, b=100,
            classe_acier="S235", classe_transv=1,
        )
        lamb_rel = comp._lamb_rel_Axe
        chi = comp._chi
        assert lamb_rel["y"] > 0
        assert lamb_rel["z"] > 0
        assert 0 < chi["y"] <= 1
        assert 0 < chi["z"] <= 1

    def test_courbes_flambement_toutes(self):
        """Test les différentes courbes de flambement (a0, a, b, c, d)."""
        for courbe in ("a0", "a", "b", "c", "d"):
            comp = Compression(
                A=3000, Iy=6750000, Iz=562500,
                lo_y=3000, lo_z=3000,
                courbe_flamb={"y": courbe, "z": courbe},
                type_appuis="Rotule - Rotule",
                t=15, h=150, b=100,
                classe_acier="S235", classe_transv=1,
            )
            _, val = comp.Nb_Rd
            assert val["y"].value > 0


class TestFlexionSupplementaire:
    """Tests complémentaires pour Flexion."""

    def test_Mc_Rd_classe_3(self):
        """Test du moment résistant pour section de classe 3 (module élastique)."""
        flexion = Flexion(W=18000, t=10, h=100, b=200, classe_acier="S235", classe_transv=3)
        latex, val = flexion.Mc_Rd
        assert val > 0
        assert abs(val.value - 18000e-9 * 235e6) < 1

    def test_Mc_V_Rd_cisaillement_eleve(self):
        """Test du moment réduit quand V_Ed/Vpl_Rd > 0.5 (rho actif)."""
        flexion = Flexion(W=20000, t=10, h=100, b=200, classe_acier="S235", classe_transv=1)
        # Vpl_Rd ≈ 67.8 kN -> V_Ed=40 kN donne V/Vpl~0.59 > 0.5
        result_eleve = flexion.Mc_V_Rd(Av=500, V_Ed=40)
        result_faible = flexion.Mc_V_Rd(Av=500, V_Ed=10)
        assert result_eleve[1].value < result_faible[1].value

    def test_Mc_V_Rd_cisaillement_faible(self):
        """Test du moment non réduit quand V_Ed/Vpl_Rd <= 0.5."""
        flexion = Flexion(W=20000, t=10, h=100, b=200, classe_acier="S235", classe_transv=1)
        mc = flexion.Mc_Rd[1]
        result = flexion.Mc_V_Rd(Av=500, V_Ed=5)
        assert abs(result[1].value - mc.value) < 1e-3


class TestTractionSupplementaire:
    """Tests complémentaires pour Traction."""

    def test_Nt_Rd_anet_egal_A(self):
        """Quand Anet=A (section nette = section brute), Nt_Rd = Npl_Rd (car Nu_Rd > Npl_Rd)."""
        traction = Traction(A=2000, Anet=2000, ass_cat_C=False,
                            t=10, b=200, h=100, classe_acier="S235", classe_transv=1)
        _, nt = traction.Nt_Rd
        _, npl = traction._Npl_Rd
        assert nt.value == npl.value

    def test_Nt_Rd_avec_anet(self):
        """Quand Anet>0, retourne min(Npl_Rd, Nu_Rd)."""
        traction = Traction(A=2000, Anet=1500, ass_cat_C=False,
                            t=10, b=200, h=100, classe_acier="S235", classe_transv=1)
        _, nt = traction.Nt_Rd
        _, npl = traction._Npl_Rd
        _, nu = traction._Nu_Rd
        assert abs(nt.value - min(npl.value, nu.value)) < 1e-3


class TestPlatValidations:
    """Tests des validations de Plat."""

    def test_classe_acier_invalide(self):
        """ValueError sur classe d'acier inconnue."""
        with pytest.raises(ValueError, match="S999"):
            Plat(t=10, h=100, b=200, classe_acier="S999", classe_transv=1)

    def test_classe_transv_invalide(self):
        """ValueError sur classe transversale 4 (non développée)."""
        with pytest.raises(ValueError):
            Plat(t=10, h=100, b=200, classe_acier="S235", classe_transv=4)

    def test_fy_fu_epaisseur_superieure_40mm(self):
        """fy réduit pour t > 40mm."""
        p_mince = Plat(t=20, h=100, b=200, classe_acier="S355", classe_transv=1)
        p_epais = Plat(t=50, h=100, b=200, classe_acier="S355", classe_transv=1)
        assert p_epais.fy.value < p_mince.fy.value


if __name__ == "__main__":
    pytest.main(["-v", "test_EC3_Element_droit.py"])
