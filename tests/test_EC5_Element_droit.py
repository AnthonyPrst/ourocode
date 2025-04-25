#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import pytest
import forallpeople as si
si.environment("structural")
sys.path.insert(1, './')
from ourocode.eurocode import EC5_Element_droit as EC5

@pytest.fixture
def barre():
    return EC5.Barre(b=100, h=300, section="Rectangulaire", classe='GL24h', cs=2, effet_systeme="True")

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
        assert self.barre.classe == 'GL24h'
        assert self.barre.cs == 2
        assert self.barre.effet_systeme == "True"

    def test_Kdef(self):
        assert self.barre.K_def == 0.8
    
    def test_Emean_fin(self):
        self.barre.Emean_fin(psy_2=0.6)
        assert self.barre.E_mean_fin.value == 7770270270.27027

    def test_fleche(self):
        self.barre.fleche(
            long=5000,            # en mm
            Ed_WinstQ=5,          # flèche instantanée sous charge variable
            Ed_Wnetfin=10,        # flèche nette finale
            Ed_Wfin=0,            # flèche finale (non utilisée ici)
            type_ele="Élément structuraux",
            type_bat="Bâtiments courants"
        )
        assert isinstance(self.barre.taux_ELS, dict)
        assert self.barre.taux_ELS["Winst(Q)"] == 30
        assert self.barre.taux_ELS["Wnet,fin"] == 40
        assert self.barre.taux_ELS["Wfin"] == 0

class Test_Flexion(object):
    @pytest.fixture(autouse=True)
    def setup_method(self, barre, load_and_combi_data):
        self.flexion = EC5.Flexion._from_parent_class(barre, lo=7000, coeflef=0.8, pos="Charge sur fibre comprimée")
        self.flexion.sigma_m_d(100, axe='y')
        self.flexion.f_m_d(**load_and_combi_data)
        self.flexion.taux_m_d()

    def test_init(self):
        assert self.flexion.lo == 7 * si.m
        assert self.flexion.coeflef == 0.8
        assert self.flexion.pos == 0
        
    def test_Kh(self):
        # print(self.flexion.K_h)
        assert self.flexion.K_h == {'y': 1.0717734625362931, 'z': 1.1}
    
    def test_Km(self):
        assert self.flexion.K_m == 0.7
    
    def test_sig_m_crit(self):
        #print(self.flexion.sigma_m_crit)
        assert self.flexion.sigma_m_crit[1].value == 40258064.51612905
    
    def test_lamb_rel_m(self):
        #print(self.flexion.lamb_rel_m)
        assert self.flexion.lamb_rel_m[1] == 0.7721099961494126
        
    def test_sigma_m_d(self):
        # print(self.flexion.sigma_m_rd)
        assert self.flexion.sigma_m_rd["y"].value ==  66666666.66666668
        assert self.flexion.sigma_m_rd["z"].value ==  0

    def test_f_m_d(self):
        assert self.flexion.f_type_rd == 16.896 * si.MPa

    def test_taux_m_d(self):
        # print(self.flexion.taux_m_rd)
        assert self.flexion.taux_m_rd == {'equ6.11': 3.681474871909751, 'equ6.12': 2.577032410336826, 'equ6.33y': 3.7530932632673393, 'equ6.33z': 0.0}


class Test_Traction:
    @pytest.fixture(autouse=True)
    def setup_method(self, barre, load_and_combi_data):
        self.traction = EC5.Traction._from_parent_class(barre)
        self.traction.sigma_t_0_d(20)
        self.traction.f_t_0_d(**load_and_combi_data)
        self.traction.taux_t_0_d()

    def test_f_t_0_d(self):
        assert self.traction.f_type_rd == 12.288 * si.MPa

    def test_sigma_t_0_d(self):
        assert self.traction.sigma_t_0_rd.value == 666666.6666666667

    def test_taux_t_0_d(self):
        assert self.traction.taux_t_0_rd == {'equ6.1': 0.05062027948875909}