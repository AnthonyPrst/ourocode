#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import pytest
import forallpeople as si
from copy import deepcopy

si.environment("structural")
sys.path.insert(1, "./")
from ourocode.eurocode import EC5_Assemblage as EC5_Assem
from ourocode.eurocode import EC5_Element_droit as EC5_Elem
from ourocode.eurocode import EC3_Element_droit as EC3_Elem

# Fixtures pour les barres de test
@pytest.fixture
def barre_bois():
    return EC5_Elem.Barre(
        b=100,
        h=200,
        section="Rectangulaire",
        classe="GL24h",
        cs=2,
        effet_systeme=True,
    )

@pytest.fixture
def panneau_bois():
    return EC5_Elem.Barre(
        b=12,
        h=1250,
        section="Rectangulaire",
        classe="OSB/3 11-18 mm",
        cs=2,
        effet_systeme=True,
    )

@pytest.fixture
def barre_metal():
    return EC3_Elem.Element(6, 200, "S275", classe_transv="2")

# Fixture pour un assemblage bois/bois
@pytest.fixture
def assemblage_bois_bois(barre_bois):
    barre1 = barre_bois
    barre2 = deepcopy(barre_bois)
    return EC5_Assem.Assemblage(beam_1=barre1, beam_2=barre2, nfile=2, nCis=1)

# Fixture pour un assemblage bois/métal
@pytest.fixture
def assemblage_bois_metal(barre_bois, barre_metal):
    return EC5_Assem.Assemblage(beam_1=barre_bois, beam_2=barre_metal, nfile=2, nCis=1)

@pytest.fixture
def pointe(barre_bois, panneau_bois):
    return EC5_Assem.Pointe(
        d=2.5,
        dh=5,
        l=55,
        qualite="6.8",
        n=9,
        alpha1=0,
        alpha2=0,
        type_organe="Autres pointes",
        percage=False,
        beam_1=panneau_bois,
        beam_2=barre_bois,
        nfile=1,
        nCis=1
    )
        
class Test_Assemblage:
    def test_init_bois_bois(self, assemblage_bois_bois):
        assert assemblage_bois_bois.type_assemblage == "Bois/Bois"
        assert assemblage_bois_bois.nfile == 2
        assert assemblage_bois_bois.nCis == 1

    def test_init_bois_metal(self, assemblage_bois_metal):
        assert assemblage_bois_metal.type_assemblage == "Bois/Métal"

    def test_rho_mean_ass_bois_bois(self, assemblage_bois_bois, barre_bois):
        # Pour un assemblage bois/bois, rho_mean_ass est la racine carrée du produit des deux rhomean
        rho_mean = (int(barre_bois.caract_meca.loc["rhomean"]) * int(barre_bois.caract_meca.loc["rhomean"])) ** 0.5
        assert assemblage_bois_bois.rho_mean_ass == rho_mean

    def test_rho_mean_ass_bois_metal(self, assemblage_bois_metal, barre_bois):
        # Pour un assemblage bois/métal, rho_mean_ass est le rhomean du bois
        assert assemblage_bois_metal.rho_mean_ass == int(barre_bois.caract_meca.loc["rhomean"])

class Test_Pointe:
    def test_init(self, pointe):
        print(pointe._type_beam)
        assert pointe.type_assemblage == "Bois/Bois"
        assert pointe._type_beam == ["PP/OSB","Bois"]
        assert pointe.d == 2.5 * si.mm
        assert pointe.dh == 5 * si.mm
        assert pointe.l == 55 * si.mm
        assert pointe.qualite == "6.8"
        assert pointe.n == 9
        assert pointe.alpha == [0, 0]
        assert pointe.type_organe == "Autres pointes"
        assert pointe.percage == False

    def test_MyRk(self, pointe):
        assert pointe.MyRk[1].value == 1.9494698713738494
    
    def test_nef(self, pointe):
        pointe.nef(a1_beam1=150, a1_beam2=150)
        assert pointe._nef == 9

    def test_Kser_ass(self, pointe):
        assert pointe.Kser_ass[1].value == 7197869.965862667

    def test_Ku_ass(self, pointe):
        assert pointe.Ku_ass[1].value == 4798579.977241778
    
    def test_prepercage(self, pointe):
        assert pointe.prepercage(beam="1", sensible=True)[1] == 35 * si.mm
    
    def test_FvRk(self, pointe):
        pointe.nef(a1_beam1=150, a1_beam2=150)
        pointe.Fax_Rk()
        Fvrk = pointe.FvRk(effet_corde=True)
        assert Fvrk[1].value == 4653.815456320474


class Test_Agrafe:
    @pytest.fixture
    def agrafe(self, barre_bois, panneau_bois):
        return EC5_Assem.Agrafe(
            d=1.5,
            b_agrafe=13.5,
            l=50,
            qualite="8.8",
            n=9,
            angle_sup_30=True,
            alpha1=0,
            alpha2=0,
            beam_1=panneau_bois,
            beam_2=barre_bois,
            nfile=1,
            nCis=1
        )

    def test_init(self, agrafe):
        assert agrafe.type_assemblage == "Bois/Bois"
        assert agrafe._type_beam == ["PP/OSB","Bois"]
        assert agrafe.d == 1.5 * si.mm
        assert agrafe.b_agrafe == 13.5 * si.mm
        assert agrafe.l == 50 * si.mm
        assert agrafe.qualite == "8.8"
        assert agrafe.n == 9
        assert agrafe.angle_sup_30 == True
        assert agrafe.alpha == [0, 0]
        assert agrafe.type_organe == "Agrafe"
        assert agrafe.percage == False

    def test_MyRk(self, agrafe):
        assert agrafe.MyRk[1].value == 0.50625
    
    def test_nef(self, agrafe):
        # Test du calcul du nombre efficace d'agrafes
        agrafe.nef()
        assert agrafe._nef == 9

    def test_Kser_ass(self, agrafe):
        assert agrafe.Kser_ass[1].value == 14281.908958018748

    def test_Ku_ass(self, agrafe):
        assert agrafe.Ku_ass[1].value == 9521.272638679166

    def test_pince(self, agrafe):
        print(agrafe.pince["barre 2"])
        pince_b1 = agrafe.pince["barre 1"]
        pince_b2 = agrafe.pince["barre 2"]
        assert pince_b1["a1"] == 30 * si.mm
        assert pince_b2["a1"] == 30 * si.mm
        assert pince_b1["a2"] == 22.5 * si.mm
        assert pince_b2["a2"] == 22.5 * si.mm
        assert pince_b1["a3t"] == 30 * si.mm
        assert pince_b2["a3t"] == 30 * si.mm
        assert pince_b1["a3c"] == 22.5 * si.mm
        assert pince_b2["a3c"] == 22.5 * si.mm
        assert pince_b1["a4t"] == 22.5 * si.mm
        assert pince_b2["a4t"] == 22.5 * si.mm
        assert pince_b1["a4c"] == 15 * si.mm
        assert pince_b2["a4c"] == 15 * si.mm

    def test_FvRk(self, agrafe):
        agrafe.nef()
        agrafe.Fax_Rk()
        Fvrk = agrafe.FvRk(effet_corde=True)
        assert Fvrk[1].value == 4862.449519109065


class Test_Boulon:
    @pytest.fixture
    def boulon(self, barre_bois, barre_metal):
        return EC5_Assem.Boulon(
            d=12,
            qualite="4.6",
            n=2,
            alpha1=0,
            alpha2=90,
            beam_1=barre_bois,
            beam_2=barre_metal,
            nfile=2,
            nCis=2
        )

    def test_init(self, boulon):
        assert boulon.type_assemblage == "Bois/Métal"
        assert boulon.d == 12 * si.mm
        assert boulon.qualite == "4.6"
        assert boulon.n == 2
        assert boulon.alpha == [0, 90]
        assert boulon.type_organe == "Boulon"

    def test_MyRk(self, boulon):
        assert boulon.MyRk[1].value == 76.74542328693614
    
    def test_nef(self, boulon):
        boulon.nef(a1_beam1=100, a1_beam2=100)
        assert boulon._nef == 1.6697284479494683

    def test_Kser_ass(self, boulon):
        kser = boulon.Kser_ass
        assert kser[1].value == 71853400.84930278

    def test_Ku_ass(self, boulon):
        ku = boulon.Ku_ass
        assert ku[1].value == 47902267.23286852

    def test_pince(self, boulon):
        pinces_b1 = boulon.pince["barre 1"]
        assert pinces_b1["a1"] == 60 * si.mm
        assert pinces_b1["a2"] == 48 * si.mm
        assert pinces_b1["a3t"] == 84 * si.mm
        assert pinces_b1["a3c"] == 48 * si.mm
        assert pinces_b1["a4t"] == 36 * si.mm
        assert pinces_b1["a4c"] == 36 * si.mm

    def test_FvRk(self, boulon):
        boulon.nef(a1_beam1=100, a1_beam2=100)
        boulon.Fax_Rk(d_int=13, d_ext=24, filetage_EN1090=True)
        Fvrk = boulon.FvRk(effet_corde=True)
        assert Fvrk[1].value == 81704.51133271086

    def test_FbsRk(self, boulon):
        boulon.FvRk(effet_corde=False)
        FbsRk = boulon.FbsRk(dp=12, a1=100, a2=100, a3t=100, num_beam=1)
        assert FbsRk[1].value == 18170552.3689843
        

