#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import sys
sys.path.insert(1, './')
from eurocode import EC5_Element_droit as EC5


class Test_Flexion(object):
    
    def setup_method(self):
        self.flexion = EC5.Flexion(b=300, h=600, lo=5000, coeflef=0.8, pos=0, classe='GL24h', cs=2)
    
    def test_lo(self):
        assert self.flexion.lo == 5000
    
    def test_coeflef(self):
        assert self.flexion.coeflef == 0.8
        
    def test_pos(self):
        assert self.flexion.pos == 0
        
    def test_Kh(self):
        #print(self.flexion.Kh)
        assert self.flexion.Kh == {'y': 1, 'z': 1.0717734625362931}
    
    def test_Km(self):
        assert self.flexion.Km == 0.7
    
    def test_sig_m_crit(self):
        #print(self.flexion.sig_m_crit)
        assert self.flexion.sig_m_crit == 216.0
    
    def test_lamb_rel_m(self):
        #print(self.flexion.lamb_rel_m)
        assert self.flexion.lamb_rel_m == 0.3333333333333333
        
    def test_sigma_m_d(self):
        self.flexion.sigma_m_d(100, axe='y')
        #print(self.flexion.sigma_m_rd)
        assert self.flexion.sigma_m_rd == {'y': 5.555555555555555, 'z': 0}
    
    def test_taux_m_d(self):
        self.flexion.f_type_d()
        self.flexion.sigma_m_d(100, axe='y')
        self.flexion.taux_m_d(20, 10, 30)
        #print(self.flexion.taux_m_rd)
        assert self.flexion.taux_m_rd == {'equ6.11': 0.48225308641975306, 'equ6.12': 0.33757716049382713, 'equ6.33y': 0.48225308641975306, 'equ6.33z': 0.0, 'equ6.35zyz': 10.232568039361379, 'equ6.35yzz': 10.482253086419753, 'equ6.35yzy': 20.482253086419753, 'equ6.35zyy': 20.23256803936138, 'equ6.17': 30.482253086419753, 'equ6.18': 30.337577160493826}