#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import pytest

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
        d_bat=20,
        b_bat=30,
        alpha_toit=15
    )

class Test_Sismique:
    @pytest.fixture
    def sismique(self, batiment):
        return EC8.Sismique._from_parent_class(
            batiment,
            cat_importance="II",
        )

    def test_init(self, sismique):
        assert sismique.cat_importance == "II"
        assert sismique.region_sismique == "Zone 4"
        assert sismique.gamma_1 == 1
        assert sismique.h_bat == 10
        assert sismique.d_bat == 20
        assert sismique.b_bat == 30
        assert sismique.alpha_toit == 15
    
    @pytest.mark.parametrize(
        "cat_importance, expected_output",
        [
            ("I", False),
            ("II", True),
            ("III", True),
            ("IV", True),
        ],
    )
    def test_has_to_be_analyzed(self, sismique, cat_importance, expected_output):
        sismique.cat_importance = cat_importance
        assert sismique._has_to_be_analyzed() == expected_output