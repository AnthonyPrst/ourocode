#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import pytest
import forallpeople as si

si.environment("structural")
sys.path.insert(1, "./")
from ourocode.eurocode import EC5_Element_droit as EC5


@pytest.fixture
def barre():
    return EC5.Barre(
        b=100,
        h=300,
        section="Rectangulaire",
        classe="GL24h",
        cs=2,
        effet_systeme="True",
    )


@pytest.fixture
def traction(barre, load_and_combi_data):
    traction = EC5.Traction._from_parent_class(barre)
    traction.sigma_t_0_d(20)
    traction.f_t_0_d(**load_and_combi_data)
    traction.taux_t_0_d()
    return traction


@pytest.fixture
def compression(barre, load_and_combi_data):
    compression = EC5.Compression._from_parent_class(
        barre, lo_y=7000, lo_z=7000, type_appuis="Rotule - Rotule"
    )
    compression.sigma_c_0_d(50)
    compression.f_c_0_d(**load_and_combi_data)
    compression.taux_c_0_d()
    return compression


@pytest.fixture
def flexion(barre, load_and_combi_data):
    flexion = EC5.Flexion._from_parent_class(
        barre, lo=7000, coeflef=0.8, pos="Charge sur fibre comprimée"
    )
    flexion.sigma_m_d(100, axe="y")
    flexion.f_m_d(**load_and_combi_data)
    flexion.taux_m_d()
    return flexion


@pytest.fixture
def load_and_combi_data():
    return {"loadtype": "Moyen terme", "typecombi": "Fondamentales"}


class Test_Barre(object):
    @pytest.fixture(autouse=True)
    def setup_method(self, barre):
        self.barre = barre

    def test_init(self):
        assert self.barre.b == 100 * si.mm
        assert self.barre.h == 300 * si.mm
        assert self.barre.section == "Rectangulaire"
        assert self.barre.classe == "GL24h"
        assert self.barre.cs == 2
        assert self.barre.effet_systeme == "True"

    def test_Kdef(self):
        assert self.barre.K_def == 0.8

    def test_Emean_fin(self):
        self.barre.Emean_fin(psy_2=0.6)
        assert self.barre.E_mean_fin.value == 7770270270.27027

    def test_fleche(self):
        self.barre.fleche(
            long=5000,  # en mm
            Ed_WinstQ=5,  # flèche instantanée sous charge variable
            Ed_Wnetfin=10,  # flèche nette finale
            Ed_Wfin=0,  # flèche finale (non utilisée ici)
            type_ele="Élément structuraux",
            type_bat="Bâtiments courants",
        )
        assert isinstance(self.barre.taux_ELS, dict)
        assert self.barre.taux_ELS["Winst(Q)"] == 30
        assert self.barre.taux_ELS["Wnet,fin"] == 40
        assert self.barre.taux_ELS["Wfin"] == 0


class Test_Flexion(object):
    @pytest.fixture(autouse=True)
    def setup_method(self, flexion, compression, traction, load_and_combi_data):
        self.flexion = flexion
        self.flexion.taux_m_d(compression, traction)

    def test_init(self):
        assert self.flexion.lo == 7 * si.m
        assert self.flexion.coeflef == 0.8
        assert self.flexion.pos == "Charge sur fibre comprimée"

    def test_Kh(self):
        assert self.flexion.K_h == {"y": 1.0717734625362931, "z": 1.1}

    def test_Km(self):
        assert self.flexion.K_m == 0.7

    def test_sig_m_crit(self):
        assert self.flexion.sigma_m_crit[1].value == 40258064.51612905

    def test_lamb_rel_m(self):
        assert self.flexion.lamb_rel_m[1] == 0.7721099961494126

    def test_sigma_m_d(self):
        assert self.flexion.sigma_m_rd["y"].value == 66666666.66666668
        assert self.flexion.sigma_m_rd["z"].value == 0

    def test_f_m_d(self):
        assert self.flexion.f_type_rd == 16.896 * si.MPa

    def test_taux_m_d(self):
        assert self.flexion.taux_m_rd == {
            "equ6.11": 3.7320951513985103,
            "equ6.12": 2.627652689825585,
            "equ6.33y": 3.7530932632673393,
            "equ6.33z": 0.0,
            "equ6.19": 3.693248628902421,
            "equ6.20": 2.5888061673294955,
            "equ6.35zyz": 15.743161072121572,
            "equ6.35yzz": 5.410545292606226,
            "equ6.35yzy": 3.8840891485379205,
            "equ6.35zyy": 13.755871509131088,
        }


class Test_Traction:
    @pytest.fixture(autouse=True)
    def setup_method(self, traction):
        self.traction = traction

    def test_f_t_0_d(self):
        assert self.traction.f_type_rd == 12.288 * si.MPa

    def test_sigma_t_0_d(self):
        assert self.traction.sigma_t_0_rd.value == 666666.6666666667

    def test_taux_t_0_d(self):
        assert self.traction.taux_t_0_rd == {"equ6.1": 0.05062027948875909}


class Test_Compression:
    @pytest.fixture(autouse=True)
    def setup_method(self, compression, flexion):
        self.compression = compression
        compression.taux_c_0_d(flexion)

    def test_f_c_0_d(self):
        assert self.compression.f_type_rd == 15.360 * si.MPa

    def test_sigma_c_0_d(self):
        assert self.compression.sigma_c_0_rd.value == 1666666.6666666667

    def test_taux_c_0_d(self):
        print(self.compression.taux_c_0_rd)
        assert self.compression.taux_c_0_rd == {
            "equ6.23": 3.8840891485379205,
            "equ6.24": 4.234484439675713,
            "equ6.2": 0.10850694444444445,
        }


class Test_Compression_Perpendiculaire:
    @pytest.fixture(autouse=True)
    def setup_method(self, barre, load_and_combi_data):
        self.compression_perpendiculaire = (
            EC5.Compression_perpendiculaire._from_parent_class(
                barre,
                b_appuis=100,
                l_appuis=50,
                l1d=7000,
                l1g=7000,
                ad=50,
                ag=0,
                type_appuis_90="Appuis discret",
            )
        )
        self.compression_perpendiculaire.sigma_c_90_d(30)
        self.compression_perpendiculaire.f_c_90_d(**load_and_combi_data)
        self.compression_perpendiculaire.taux_c_90_d()

    def test_f_c_90_d(self):
        assert self.compression_perpendiculaire.f_type_rd == 1.600 * si.MPa

    def test_sigma_c_90_d(self):
        assert self.compression_perpendiculaire.sigma_c_90_rd.value == 3750000.0

    def test_taux_c_90_d(self):
        assert self.compression_perpendiculaire.taux_c_90_rd == {
            "equ6.3": 1.3392857142857142
        }

    @pytest.mark.parametrize(
        "type_appuis_90, classe, h_calcul, expected_output",
        [
            ("Appuis discret", "C24", 300, 1.5),
            ("Appuis discret", "GL24h", 300, 1.75),
            ("Appuis continu", "C24", 300, 1.25),
            ("Appuis continu", "GL24h", 300, 1.5),
            ("Appuis discret", "C24", 450, 1.5),
            ("Appuis discret", "GL24h", 450, 1.75),
            ("Appuis continu", "C24", 450, 1.5),
            ("Appuis continu", "GL24h", 450, 1.75),
        ],
    )
    def test_K_c90(self, type_appuis_90, classe, h_calcul, expected_output):
        self.compression_perpendiculaire.type_appuis_90 = type_appuis_90
        self.compression_perpendiculaire.classe = classe
        self.compression_perpendiculaire.h_calcul = h_calcul * si.mm
        assert self.compression_perpendiculaire.K_c90 == expected_output


class Test_Cisaillement:
    @pytest.fixture(autouse=True)
    def setup_method(self, barre, load_and_combi_data):
        self.cisaillement = EC5.Cisaillement._from_parent_class(barre)
        self.cisaillement.Kv(200, 40, 0, ent="Dessous")
        self.cisaillement.tau_d(50)
        self.cisaillement.f_v_d(**load_and_combi_data)
        self.cisaillement.taux_tau_d()

    def test_Kv(self):
        assert self.cisaillement.K_v == 0.6459233802767231

    @pytest.mark.parametrize(
        "cs, classe, h_calcul, expected_output",
        [
            (1, "C24", 200, 0.67),
            (1, "GL24h", 200, 1),
            (2, "C24", 200, 0.67),
            (2, "GL24h", 200, 0.67),
            (3, "GL24h", 200, 0.67)
        ],
    )
    def test_K_cr(self, cs, classe, h_calcul, expected_output):
        self.cisaillement.cs = cs
        self.cisaillement.classe = classe
        self.cisaillement.h_calcul = h_calcul * si.mm
        assert self.cisaillement.K_cr == expected_output

    def test_f_v_d(self):
        assert self.cisaillement.f_type_rd == 2.240 * si.MPa

    def test_tau_d(self):
        assert self.cisaillement.tau_rd.value == 5597014.925373133

    def test_taux_tau_d(self):
        assert self.cisaillement.taux_tau_rd == {
            "equ6.13": 2.49866737739872,
            "equ6.60": 3.868364969740303,
        }


class Test_Compression_Inclinees:
    @pytest.fixture(autouse=True)
    def setup_method(self, compression, load_and_combi_data):
        compression_perpendiculaire = (
            EC5.Compression_perpendiculaire._from_parent_class(
                compression,
                b_appuis=100,
                l_appuis=50,
                l1d=7000,
                l1g=7000,
                ad=0,
                ag=50,
                type_appuis_90="Appuis discret",
            )
        )
        self.compression_inclinees = EC5.Compression_inclinees._from_parent_class(
            compression_perpendiculaire, alpha=45  # Angle d'inclinaison en degrés
        )
        self.compression_inclinees.sigma_c_alpha_d(50)
        self.compression_inclinees.taux_c_alpha_d(**load_and_combi_data)

    def test_sigma_c_alpha_d(self):
        assert self.compression_inclinees.sigma_c_alpha_rd.value == 9999999.999999998

    def test_taux_c_alpha_d(self):
        assert self.compression_inclinees.taux_c_alpha_rd == {
            "equ6.16": 2.111235119047619
        }

    def test_alpha(self):
        assert self.compression_inclinees.alpha == 45
