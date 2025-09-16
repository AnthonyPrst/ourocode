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
        assert pointe.Kser_ass()[1].value == 7197869.965862667

    def test_Ku_ass(self, pointe):
        assert pointe.Ku_ass()[1].value == 4798579.977241778
    
    def test_prepercage(self, pointe):
        assert pointe.prepercage(beam="1", sensible=True)[1] == 35 * si.mm
    
    def test_FvRk(self, pointe):
        pointe.nef(a1_beam1=150, a1_beam2=150)
        pointe.Fax_Rk()
        Fvrk = pointe.FvRk(effet_corde=True)
        assert Fvrk[1].value == 5537.902373351718


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
        assert agrafe.Kser_ass()[1].value == 14281.908958018748

    def test_Ku_ass(self, agrafe):
        assert agrafe.Ku_ass()[1].value == 9521.272638679166

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
        assert Fvrk[1].value == 5257.055494983215


class Test_Boulon:
    @pytest.fixture
    def boulon(self, barre_bois):
        return EC5_Assem.Boulon(
            d=12,
            qualite="4.6",
            n=2,
            alpha1=45,
            alpha2=0,
            beam_1=barre_bois,
            beam_2=barre_bois,
            nfile=2,
            nCis=2
        )

    def test_init(self, boulon):
        assert boulon.type_assemblage == "Bois/Bois"
        assert boulon.d == 12 * si.mm
        assert boulon.qualite == "4.6"
        assert boulon.n == 2
        assert boulon.alpha == [45, 0]
        assert boulon.type_organe == "Boulon"

    def test_MyRk(self, boulon):
        assert boulon.MyRk[1].value == 76.74542328693614
    
    def test_nef(self, boulon):
        boulon.nef(a1_beam1=100, a1_beam2=100)
        assert boulon._nef == 1.6697284479494683

    def test_Kser_ass(self, boulon):
        kser = boulon.Kser_ass()
        assert kser[1].value == 35926700.42465139

    def test_Ku_ass(self, boulon):
        ku = boulon.Ku_ass()
        assert ku[1].value == 23951133.61643426

    def test_pince(self, boulon):
        pinces_b1 = boulon.pince["barre 1"]
        pinces_b2 = boulon.pince["barre 2"]
        assert pinces_b1["a1"].value == 0.056485281374238576
        assert pinces_b1["a2"] == 48 * si.mm
        assert pinces_b1["a3t"] == 84 * si.mm
        assert pinces_b1["a3c"].value == 0.06291168824543143
        assert pinces_b1["a4t"].value == 0.040970562748477136
        assert pinces_b1["a4c"].value == 0.036000000000000004
        assert pinces_b2["a1"] == 60 * si.mm
        assert pinces_b2["a2"] == 48 * si.mm
        assert pinces_b2["a3t"] == 84 * si.mm
        assert pinces_b2["a3c"] == 48 * si.mm
        assert pinces_b2["a4t"] == 36 * si.mm
        assert pinces_b2["a4c"] == 36 * si.mm

    def test_FvRk(self, boulon):
        boulon.nef(a1_beam1=100, a1_beam2=100)
        boulon.Fax_Rk(d_int=13, d_ext=24, filetage_EN1090=True)
        Fvrk = boulon.FvRk(effet_corde=True)
        assert Fvrk[1].value == 55632.21408197254

    def test_FbsRk(self, boulon):
        boulon.FvRk(effet_corde=False)
        FbsRk = boulon.FbsRk(dp=12, a1=100, a2=100, a3t=100, num_beam=1)
        assert FbsRk[1].value == 253440.0

    def test_F90Rk(self, boulon):
        F90Rk = boulon.F90Rk(b=100, h=200, he=150)
        F90Rd = boulon.F_Rd(F90Rk[1].value*10**-3, loadtype="Moyen terme")
        assert F90Rk[1].value == 34292.85639896449
        assert F90Rd[1].value == 21103.29624551661


class Test_Broche:
    @pytest.fixture
    def broche(self, barre_bois, barre_metal):
        return EC5_Assem.Broche(
            d=16,
            qualite="4.6",
            n=2,
            alpha1=0,
            alpha2=90,
            beam_1=barre_bois,
            beam_2=barre_metal,
            nfile=2,
            nCis=2
        )

    def test_init(self, broche):
        assert broche.type_assemblage == "Bois/Métal"
        assert broche.d == 16 * si.mm
        assert broche.qualite == "4.6"
        assert broche.n == 2
        assert broche.alpha == [0, 90]
        assert broche.type_organe == "Broche"
        assert broche.FaxRk == 0

    def test_pince(self, broche):
        pinces = broche.pince
        assert pinces["barre 1"]["a1"] == 80 * si.mm    
        assert pinces["barre 1"]["a2"] == 48 * si.mm
        assert pinces["barre 1"]["a3t"] == 112 * si.mm
        assert pinces["barre 1"]["a3c"] == 56 * si.mm
        assert pinces["barre 1"]["a4t"] == 48 * si.mm
        assert pinces["barre 1"]["a4c"] == 48 * si.mm
    
    def test_Kser_ass(self, broche):
        kser = broche.Kser_ass()
        assert kser[1].value == 95804534.46573706

    def test_Ku_ass(self, broche):
        ku = broche.Ku_ass()
        assert ku[1].value == 63869689.643824704

    def test_FvRk(self, broche):
        broche.nef(a1_beam1=100, a1_beam2=100)
        Fvrk = broche.FvRk(effet_corde=False)
        assert Fvrk[1].value == 118572.10105692707

    def test_FbsRk(self, broche):
        broche.FvRk(effet_corde=False)
        FbsRk = broche.FbsRk(dp=22, a1=100, a2=100, a3t=100, num_beam=1)
        assert FbsRk[1].value == 21467088.26594395


class Test_Tirefond_inf_7:
    @pytest.fixture
    def tirefond(self, barre_bois):
        return EC5_Assem.Tirefond_inf_7(
            d=6,
            d1=4.5,
            ds=3.5,
            dh=10,
            l=180,
            n=2,
            rho_a=380,
            fhead=15,
            ftensk=10000,
            alpha1=0,
            alpha2=90,
            percage=False,
            beam_1=barre_bois,
            beam_2=barre_bois,
            nfile=2,
            nCis=1
        )

    def test_init(self, tirefond):
        assert tirefond.type_assemblage == "Bois/Bois"
        assert tirefond.d == 4.95 * si.mm
        assert tirefond.qualite == "6.8"
        assert tirefond.n == 2
        assert tirefond.alpha == [0, 90]
        assert tirefond.type_organe == "Tirefond"
        assert tirefond.percage is False

    def test_MyRk(self, tirefond):
        assert tirefond.MyRk[1].value == 11.514524874935903
    
    def test_nef(self, tirefond):
        tirefond.nef(a1_beam1=60, a1_beam2=60)
        assert tirefond._nef == 1.8025009252216604

    def test_nefTraction(self, tirefond):
        assert tirefond.nefTraction[1] == 3.4822022531844965

    def test_pince(self, tirefond):
        pinces = tirefond.pince
        assert pinces["barre 1"]["a1"] == 60 * si.mm
        assert pinces["barre 1"]["a2"] == 30 * si.mm
        assert pinces["barre 1"]["a3t"] == 90 * si.mm
        assert pinces["barre 1"]["a3c"] == 60 * si.mm
        assert pinces["barre 1"]["a4t"] == 30 * si.mm
        assert pinces["barre 1"]["a4c"] == 30 * si.mm
        assert pinces["barre 2"]["a1"] == 30 * si.mm
        assert pinces["barre 2"]["a2"] == 30 * si.mm
        assert pinces["barre 2"]["a3t"] == 60 * si.mm
        assert pinces["barre 2"]["a3c"] == 60 * si.mm
        assert pinces["barre 2"]["a4t"] == 42 * si.mm
        assert pinces["barre 2"]["a4c"] == 30 * si.mm

    def test_FvRk(self, tirefond):
        tirefond.nef(a1_beam1=60, a1_beam2=60)
        tirefond.Fax_Rk(faxk=10, l_ef=80, alpha=90)
        Fvrk = tirefond.FvRk(effet_corde=True)
        assert Fvrk[1].value == 7553.234585878883

    def test_FaxRk(self, tirefond):
        FaxRk = tirefond.Fax_Rk(faxk=10, l_ef=80, alpha=90)
        assert FaxRk[1].value == 5278.213553724563

    def test_FtRk(self, tirefond):
        FtRk = tirefond.FtRk()
        assert FtRk[1].value == 34822.02253184497

    def test_FbsRk(self, tirefond):
        tirefond.nef(a1_beam1=60, a1_beam2=60)
        tirefond.Fax_Rk(faxk=10, l_ef=80, alpha=90)
        tirefond.FvRk(effet_corde=True)
        FbsRk = tirefond.FbsRk(dp=6, a1=60, a2=60, a3t=100, num_beam=1)
        assert FbsRk[1].value == 155520.0

class Test_Tirefond_sup_6:
    @pytest.fixture
    def tirefond(self, barre_bois):
        return EC5_Assem.Tirefond_sup_6(
            d=11,
            d1=9.5,
            ds=7,
            dh=18,
            l=180,
            n=2,
            rho_a=380,
            fhead=15,
            ftensk=10000,
            MyRk=80000,
            alpha1=0,
            alpha2=90,
            beam_1=barre_bois,
            beam_2=barre_bois,
            nfile=1,
            nCis=1
        )

    def test_init(self, tirefond):
        assert tirefond.type_assemblage == "Bois/Bois"
        assert tirefond.d == 10.45 * si.mm
        assert tirefond.qualite == "6.8"
        assert tirefond.n == 2
        assert tirefond.alpha == [0, 90]
        assert tirefond.type_organe == "Tirefond"

    def test_MyRk(self, tirefond):
        assert tirefond.MyRk[1].value == 80.0
    
    def test_nef(self, tirefond):
        tirefond.nef(a1_beam1=100, a1_beam2=100)
        assert tirefond._nef == 1.7064477171757193

    def test_nefTraction(self, tirefond):
        assert tirefond.nefTraction[1] == 1.8660659830736148

    def test_pince(self, tirefond):
        pinces = tirefond.pince
        assert pinces["barre 1"]["a1"] == 55 * si.mm
        assert pinces["barre 1"]["a2"] == 44 * si.mm
        assert pinces["barre 1"]["a3t"] == 80 * si.mm
        assert pinces["barre 1"]["a3c"] == 44 * si.mm
        assert pinces["barre 1"]["a4t"] == 33 * si.mm
        assert pinces["barre 1"]["a4c"] == 33 * si.mm
        assert pinces["barre 2"]["a1"] == 44 * si.mm
        assert pinces["barre 2"]["a2"] == 44 * si.mm
        assert pinces["barre 2"]["a3t"] == 80 * si.mm
        assert pinces["barre 2"]["a3c"] == 77 * si.mm
        assert pinces["barre 2"]["a4t"] == 44 * si.mm
        assert pinces["barre 2"]["a4c"] == 33 * si.mm

    def test_FvRk(self, tirefond):
        tirefond.nef(a1_beam1=60, a1_beam2=60)
        tirefond.Fax_Rk(faxk=10, l_ef=80, alpha=90)
        Fvrk = tirefond.FvRk(effet_corde=True)
        assert Fvrk[1].value == 11295.461411468756

    def test_FaxRk(self, tirefond):
        FaxRk = tirefond.Fax_Rk(faxk=10, l_ef=80, alpha=90)
        assert FaxRk[1].value == 9164.41973069982

    def test_FtRk(self, tirefond):
        FtRk = tirefond.FtRk()
        assert FtRk[1].value == 18660.65983073615

        