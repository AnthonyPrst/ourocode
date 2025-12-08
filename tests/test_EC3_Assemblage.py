#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import pytest
import forallpeople as si

si.environment("structural")
sys.path.insert(1, "./")

from ourocode.eurocode import EC3_Assemblage as EC3_Assem
from ourocode.eurocode import EC3_Element_droit as EC3_Elem


@pytest.fixture
def plat_acier():
    return EC3_Elem.Plat(t=10, h=200, classe_acier="S235", classe_transv=1)


@pytest.fixture
def tige(plat_acier):
    # Utilise la méthode _from_parent_class comme dans le bloc __main__ du module
    return EC3_Assem.Tige._from_parent_class(
        plat_acier,
        d=16,
        d0=18,
        qualite="4.6",
        verif_filetage=True,
        filetage_EN1090=True,
    )


@pytest.fixture
def soudure(plat_acier):
    # Paramètres choisis pour respecter les vérifications de verif_soudure
    return EC3_Assem.Soudure._from_parent_class(
        plat_acier,
        t2=10,
        gorge=5,
        l=60,
        retour_soudure=True,
        alpha=90,
    )


@pytest.fixture
def tige_scellement_plaque(tige):
    # Tige de scellement avec plaque circulaire, paramètres choisis pour respecter les vérifications
    return EC3_Assem.Tige_scellement_plaque_ancrage._from_parent_class(
        tige,
        tr=9,
        dr=50,
        l_ancr=200,
        d1=100,
        e_tiges=120,
        n=4,
        n_traction=2,
        classe_beton="C25/30"
    )


@pytest.fixture
def tige_scellement_crochet(tige):
    # Tige de scellement avec crochet, paramètres choisis pour satisfaire les contraintes géométriques
    return EC3_Assem.Tige_scellement_crochet._from_parent_class(
        tige,
        l1=200,
        l2=30,
        rayon_crochet=60,
        n=4,
        n_traction=2,
        classe_beton="C25/30"
    )


class Test_Tige:
    def test_init(self, tige):
        assert tige.d == 16 * si.mm
        assert tige.d0 == 18 * si.mm
        assert tige.qualite == "4.6"
        assert tige.verif_filetage is True
        assert tige.percage_surdimensionne is False
        assert tige.fyb > 0 * si.MPa
        assert tige.fub > 0 * si.MPa
        assert tige.As > 0 * si.mm**2
        assert tige.An > 0 * si.mm**2

    def test_pince_metal_boulon(self, tige):
        pince = tige.pince_metal_boulon(trous_oblongs=False, corrosion=False)
        # Vérifie la présence des clés attendues
        assert "e1" in pince
        assert "e2" in pince
        assert "p1" in pince
        assert "p2" in pince
        # Pour d0 = 18 mm, les minimums doivent suivre l'EN 1993-1-8 (formules du code)
        assert pince["e1"]["e1_min"] == pytest.approx(21.6)
        assert pince["p1"]["p1_min"] == pytest.approx(39.6)

    def test_FvRd_FtRd(self, tige):
        FvRd = tige.FvRd
        FtRd = tige.FtRd
        # FvRd est un dict avec deux plans de cisaillement
        assert "filetage" in FvRd
        assert "lisse" in FvRd
        assert FvRd["filetage"][1].value > 0
        assert FvRd["lisse"][1].value > 0
        # FtRd est un tuple (latex, valeur)
        assert FtRd[1].value > 0

    def test_taux_FvEd_FtEd(self, tige):
        latex, taux = tige.taux_FvEd_FtEd(FvEd=10, FtEd=5)
        assert isinstance(latex, str)
        assert "taux_t_d" in taux
        assert "taux_v_d_lisse" in taux
        assert taux["taux_t_d"] >= 0
        assert taux["taux_v_d_lisse"] >= 0

    def test_Bp_Rd(self, tige):
        # Test d'une valeur simple sans vérifier le résultat exact
        Bp_Rd = tige.Bp_Rd(d_ecrou=30, d_head_bl=30)
        assert Bp_Rd[1].value > 0

    def test_Fb_Rd_et_taux_Fb_d(self, tige):
        Fb = tige.Fb_Rd(e1=40, e2=40, p1=80, p2=80)
        assert "F_bx_Rd" in Fb[1]
        assert "F_by_Rd" in Fb[1]
        assert Fb[1]["F_bx_Rd"].value > 0
        assert Fb[1]["F_by_Rd"].value > 0

        res = tige.taux_Fb_d(FvEd=50, alpha=30, FbxRd=Fb[1]["F_bx_Rd"].value * 10**-3, FbyRd=Fb[1]["F_by_Rd"].value * 10**-3)
        assert isinstance(res[1], tuple)
        assert res[1][0] >= 0
        assert res[1][1] >= 0
        assert res[1][2] >= 0

    def test_Veff_Rd_et_taux_Veff_d(self, tige):
        Veff = tige.Veff_Rd(Lnt=200, Lvt=200, effort="Centré")
        assert Veff[1].value > 0

        taux = tige.taux_Veff_d(N_Ed=50, N_Veff_Rd=Veff[1].value * 10**-3, V_Ed=25, V_Veff_Rd=Veff[1].value * 10**-3)
        assert taux[1] >= 0

    def test_axe_articulation(self, tige):
        res = tige.axe_articulation(FvEd=100, jeu=2, t_flasque=10, t_ame=8)
        assert len(res[1]) == 5
        assert tige.taux_axe_articulation["taux_v_d"] == res[1][0]
        assert tige.taux_axe_articulation["taux_m_d"] == res[1][1]
        assert tige.taux_axe_articulation["taux_m_v_d"] == res[1][2]
        assert tige.taux_axe_articulation["taux_f_b_flasque_d"] == res[1][3]
        assert tige.taux_axe_articulation["taux_f_b_ame_d"] == res[1][4]


class Test_Tige_Erreurs:
    def test_trou_surdimensionne_declenche_erro(self, plat_acier):
        # Pour d = 16 mm, un trou normal est d0 <= d + 3 ou 4 mm (selon le cas).
        # On force un trou très surdimensionné pour déclencher l'exception.
        with pytest.raises(ValueError):
            EC3_Assem.Tige._from_parent_class(
                plat_acier,
                d=16,
                d0=30,
                qualite="4.6",
                verif_filetage=False,
                filetage_EN1090=True,
            )


class Test_Soudure:
    def test_init(self, soudure):
        assert soudure.t == 10 * si.mm
        assert soudure.t2 == 10 * si.mm
        assert soudure.gorge == 5 * si.mm
        assert soudure.l == 60 * si.mm
        assert soudure.retour_soudure is True
        assert soudure.alpha == 90

    def test_beta_w_et_lef(self, soudure):
        assert soudure.beta_w > 0
        # Avec retour de soudure, lef = l
        assert soudure.lef == 60 * si.mm

    def test_cordon_frontal_lateral(self, soudure):
        cf = soudure.cordon_frontal(N_Ed=100)
        cl = soudure.cordon_lateral(V_Ed=50)
        assert cf[1] > 0
        assert cl[1] > 0

    def test_cordon_oblique_et_pieces_obliques(self, soudure):
        co = soudure.cordon_oblique(alpha_cordon=30, N_Ed=80)
        cp = soudure.cordon_pieces_obliques(N_Ed=80)
        assert co[1] > 0
        assert cp[1] > 0

    def test_critere_generale(self, soudure):
        taux = soudure.critere_generale(FvEd=50, FaxEd=100)
        assert taux[1] >= 0

    def test_soudure_discontinue(self, soudure):
        dims = soudure.soudure_discontinue(b=200, b1=100, corrosion=False)
        assert "Lwe" in dims
        assert "L1" in dims
        assert "L2" in dims
        assert dims["L1"].value > 0
        assert dims["L2"].value > 0

    def test_soudure_discontinue_corrosion(self, soudure):
        # En ambiance corrosive, la fonction retourne une exception
        with pytest.raises(ValueError):
            soudure.soudure_discontinue(b=200, b1=100, corrosion=True)


class Test_Tige_scellement_plaque_ancrage:
    def test_init(self, tige_scellement_plaque):
        assert tige_scellement_plaque.d == 16 * si.mm
        assert tige_scellement_plaque.d0 == 18 * si.mm
        assert tige_scellement_plaque.qualite == "4.6"
        assert tige_scellement_plaque.n == 4
        assert tige_scellement_plaque.n_traction == 2
        assert tige_scellement_plaque.tr == 9 * si.mm
        assert tige_scellement_plaque.dr == 50 * si.mm
        assert tige_scellement_plaque.l_ancr == 200 * si.mm
        assert tige_scellement_plaque.d1 == 100 * si.mm
        assert tige_scellement_plaque.e_tiges == 120 * si.mm
        assert tige_scellement_plaque.type_ancrage == "Plaque d'ancrage circulaire"
        # Vérifie que la classe de béton et fck sont bien définies
        assert tige_scellement_plaque.fck.value > 0

    def test_FtcRd(self, tige_scellement_plaque):
        FtcRd = tige_scellement_plaque.FtcRd
        assert FtcRd[1].value > 0

    def test_FvRd_specifique(self, tige_scellement_plaque):
        # FvRd est redéfini dans _Tige_scellement
        FvRd = tige_scellement_plaque.FvRd
        assert FvRd[1].value > 0

    def test_taux_ancrage_traction_seule(self, tige_scellement_plaque):
        latex, taux = tige_scellement_plaque.taux_ancrage(FtEd=100, FvEd=0)
        assert isinstance(latex, str)
        assert taux >= 0

    def test_taux_ancrage_combine(self, tige_scellement_plaque):
        latex, taux = tige_scellement_plaque.taux_ancrage(FtEd=100, FvEd=50)
        assert isinstance(latex, str)
        assert len(taux) == 2
        assert taux[0] >= 0
        assert taux[1] >= 0


class Test_Tige_scellement_crochet:
    def test_init(self, tige_scellement_crochet):
        assert tige_scellement_crochet.d == 16 * si.mm
        assert tige_scellement_crochet.d0 == 18 * si.mm
        assert tige_scellement_crochet.qualite == "4.6"
        assert tige_scellement_crochet.n == 4
        assert tige_scellement_crochet.n_traction == 2
        assert tige_scellement_crochet.l1 == 200 * si.mm
        assert tige_scellement_crochet.l2 == 30 * si.mm
        assert tige_scellement_crochet.rayon_crochet == 60 * si.mm
        assert tige_scellement_crochet.type_ancrage == "Crochet"
        assert tige_scellement_crochet.fck.value > 0

    def test_FtcRd(self, tige_scellement_crochet):
        FtcRd = tige_scellement_crochet.FtcRd
        assert FtcRd[1].value > 0

    def test_FvRd_specifique(self, tige_scellement_crochet):
        FvRd = tige_scellement_crochet.FvRd
        assert FvRd[1].value > 0

    def test_taux_ancrage_traction_seule(self, tige_scellement_crochet):
        latex, taux = tige_scellement_crochet.taux_ancrage(FtEd=80, FvEd=0)
        assert isinstance(latex, str)
        assert taux >= 0

    def test_taux_ancrage_combine(self, tige_scellement_crochet):
        latex, taux = tige_scellement_crochet.taux_ancrage(FtEd=80, FvEd=40)
        assert isinstance(latex, str)
        assert len(taux) == 2
        assert taux[0] >= 0
        assert taux[1] >= 0
