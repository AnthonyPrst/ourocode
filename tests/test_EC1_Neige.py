#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import pytest
import math

from ourocode.eurocode.EC1_Neige import Neige
import forallpeople as si


def _norm(v):
    """Normalise les valeurs pouvant être des tuples (latex, value) ou des scalaires.
    Retourne toujours la valeur numérique avec unités si présent.
    """
    if isinstance(v, tuple):
        # handcalc retourne (latex, value)
        return v[1]
    return v


@pytest.fixture
def neige_1v_normal_bas_alt():
    # Commune INSEE 1001 -> Zone_neige A2 (S_k_200 = 0.45 kN/m² ; S_Ad = 1.0 kN/m²)
    return Neige(
        exposition="Normal",
        type_toit="1 versant",
        blocage=False,
        # Paramètres Batiment
        h_bat=5, d_bat=10, b_bat=10,
        alpha_toit=10,  # deg
        code_INSEE=1001,
        alt=100,  # m
    )


@pytest.fixture
def neige_2v_protege_alt_moy():
    # Commune INSEE 1002 -> Zone_neige C2 (S_k_200 = 0.65 kN/m² ; S_Ad = 1.35 kN/m²)
    return Neige(
        exposition="Protégé",
        type_toit="2 versants",
        blocage=False,
        # Paramètres Batiment
        h_bat=6, d_bat=12, b_bat=18,
        alpha_toit=10,   # versant 1
        alpha_toit2=40,  # versant 2
        code_INSEE=1002,
        alt=400,  # m
    )


def test_region_Sk_SAd(neige_1v_normal_bas_alt):
    n = neige_1v_normal_bas_alt
    assert n.region_neige == "A2"
    assert pytest.approx(n.S_k_200.value, rel=1e-6) == 450
    assert pytest.approx(n.S_Ad.value, rel=1e-6) == 1000


def test_Ce_exposition(neige_1v_normal_bas_alt, neige_2v_protege_alt_moy):
    assert neige_1v_normal_bas_alt.C_e == 1
    assert neige_2v_protege_alt_moy.C_e == 1.25


def test_mu_1_versant(neige_1v_normal_bas_alt):
    # alpha_toit = 10° -> mu1 = 0.8
    mu1 = _norm(neige_1v_normal_bas_alt.mu["mu1"])
    assert pytest.approx(mu1, rel=1e-6) == 0.8


def test_mu_2_versants(neige_2v_protege_alt_moy):
    # alpha_toit1 = 10° -> 0.8 ; alpha_toit2 = 40° -> 0.8*(60-40)/30 = 0.533333...
    mu = neige_2v_protege_alt_moy.mu
    mu_v1 = _norm(mu["mu2 versant 1"])  # 0.8
    mu_v2 = _norm(mu["mu2 versant 2"])  # ~0.533333
    assert pytest.approx(mu_v1, rel=1e-6) == 0.8
    assert pytest.approx(mu_v2, rel=1e-6) == (0.8 * (60 - 40) / 30)


def test_Sk_alt_branches(neige_1v_normal_bas_alt, neige_2v_protege_alt_moy):
    # alt <= 200: Sk_alt = S_k_200
    sk1 = _norm(neige_1v_normal_bas_alt.Sk_alt)
    assert pytest.approx(sk1.value, rel=1e-6) == neige_1v_normal_bas_alt.S_k_200.value

    # 200 < alt <= 500 (zone C2) -> delta_S_2 (branche zone E) ou delta_S_1 sinon
    sk2 = _norm(neige_2v_protege_alt_moy.Sk_alt)
    # Vérifier cohérence: supérieur à S_k_200
    assert sk2.value > neige_2v_protege_alt_moy.S_k_200.value


def test_Sn_1_versant(neige_1v_normal_bas_alt):
    n = neige_1v_normal_bas_alt
    Ce = n.C_e
    Ct = n.C_t
    Sk = _norm(n.Sk_alt)
    mu1 = _norm(n.mu["mu1"])  # 0.8
    expected = Ce * Ct * Sk * mu1

    latex, val = n.Sn
    assert pytest.approx(val.value, rel=1e-6) == expected.value


def test_Sx_1_versant(neige_1v_normal_bas_alt):
    n = neige_1v_normal_bas_alt
    Ce = n.C_e
    Ct = n.C_t
    Sa = n.S_Ad
    mu1 = _norm(n.mu["mu1"])  # 0.8
    expected = Ce * Ct * Sa * mu1

    val = n.Sx
    # Sx retourne directement (latex, value) pour 1 versant
    assert pytest.approx(val[1].value, rel=1e-6) == expected.value


def test_Se_alt_threshold(neige_1v_normal_bas_alt):
    # alt <= 900 -> 0
    n = neige_1v_normal_bas_alt
    if n.type_toit == "1 versant":
        assert _norm(n.Se) == 0

    # alt > 900 -> > 0
    n2 = Neige(
        exposition="Normal",
        type_toit="1 versant",
        blocage=False,
        h_bat=5, d_bat=10, b_bat=10,
        alpha_toit=10,
        code_INSEE=1001,
        alt=950,
    )
    assert _norm(n2.Se).value > 0


def test_mu2_construction_attenante():
    n = Neige(
        exposition="Normal",
        type_toit="2 versants",
        blocage=False,
        h_bat=5, d_bat=10, b_bat=10,
        alpha_toit=15,
        alpha_toit2=15,
        code_INSEE=1001,
        alt=600,
    )

    res = n.mu2_construction_attenante(b1=10*si.m, b2=8*si.m, h=2*si.m, alpha_toit_attenant=10)
    assert set(res.keys()) == {"mu_W", "mu_S", "mu2", "ls"}
    assert res["mu_S"] == 0
    # ls borné par 5 m si 2h < 5 m
    assert res["ls"].value == pytest.approx(5.0)
    # mu2 = 0 + mu_W
    assert pytest.approx(res["mu2"].value, rel=1e-6) == res["mu_W"].value


def test_fs():
    n = Neige(
        exposition="Normal",
        type_toit="1 versant",
        blocage=False,
        h_bat=5, d_bat=10, b_bat=10,
        alpha_toit=10,
        code_INSEE=1001,
        alt=100,
    )
    out = n.fs(Sn=1.0, entraxe=1.5, alpha=30)
    # F_s = S_n * entraxe * sin(alpha)
    expected = 1.0 * si.kN/si.m**2 * 1.5 * si.m * math.sin(math.radians(30))
    assert pytest.approx(out[1].value, rel=1e-6) == expected.value
