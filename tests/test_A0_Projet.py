#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
# pytest --cov=. --cov-report html
import sys
import os
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

# Ajout du chemin du projet au PYTHONPATH
sys.path.append(str(Path(__file__).parent.parent))

# Import des modules à tester
from ourocode.eurocode.A0_Projet import Wood_beam_model
import forallpeople as si


@pytest.fixture
def wood_beam_simple():
    wb = Wood_beam_model(
        longueur=6000,  # mm
        b=100,          # mm
        h=200,          # mm
        section="Rectangulaire",
        classe="C24",
        poids_propre=True,
        nbr_appuis=3,
        l_appuis=50,
        devers=5,
        inclinaison=0,
    )
    return wb


class Test_WoodBeamModel:
    def test_initialization_topology(self, wood_beam_simple):
        wb = wood_beam_simple

        # 3 appuis -> 3 noeuds, 2 membrures, 3 appuis définis
        assert len(wb.get_all_nodes()) == 3
        assert len(wb.get_all_members()) == 2
        assert len(wb.get_all_supports()) == 3

        # Section et matériau créés avec les bons identifiants
        assert "R100X200" in wb.get_all_sections()
        assert "C24" in wb.get_all_materials()

        # Rotation (dévers) transmise à la membrure
        assert wb.get_member("M1")["Rotation"] == 5

    def test_convert_beam_pos(self, wood_beam_simple):
        wb = wood_beam_simple

        pos_start, mem_start = wb._convert_beam_pos("start")
        pos_end, mem_end = wb._convert_beam_pos("end")

        # Début à 0 mm, fin à la longueur de la dernière membrure
        assert pos_start == 0
        last_mem = f"M{len(wb.get_all_members())}"
        assert mem_end == last_mem
        assert pos_end == wb.get_member_length(last_mem).value*10**3

    def test_create_dist_load(self, wood_beam_simple):
        wb = wood_beam_simple

        wb.create_dist_load(
            name="q_exploitation",
            start_load=-2.5,
            end_load=-2.5,
            start_pos="start",
            end_pos="end",
            action="Exploitation Q",
            direction="FY",
        )

        all_loads = wb.get_all_loads()
        load = all_loads["L4"]
        assert len(all_loads) == 4
        # Le wrapper applique la charge sur la membrure correspondant au dernier appel de position (ici "end")
        assert load["Type de charge"] == "Distribuée"
        assert load["Nom"] == "q_exploitation"
        assert load["Action"] == "Exploitation Q"
        assert load["Axe"] == "FY"
        assert "Position" in load and "start" in load["Position"] and "end" in load["Position"]
        assert hasattr(load["Position"]["start"], "value") and load["Position"]["start"].value == 0
        # Vérifie que la membrure ciblée est la dernière
        assert load["N° barre"] == f"M{len(wb.get_all_members())}"

    def test_create_point_load(self, wood_beam_simple):
        wb = wood_beam_simple

        p = wb.create_point_load(
            name="P_cent",
            load=-10.0,
            pos="58%",
            action="Permanente G",
            direction="FY",
        )

        # Type de charge (F -> Concentrée) et cohérence des champs
        assert p["Type de charge"] == "Concentrée"
        assert p["Nom"] == "P_cent"
        assert p["Action"] == "Permanente G"
        assert p["Axe"] == "FY"
        assert hasattr(p["Position"], "value")
        assert 0 < p["Position"].value < wb.get_member_length(p["N° barre"]).value

    def test_inclinaison_nodes(self):
        # Inclinaison 90° -> barre verticale: X ~ 0, Y croissant
        wb = Wood_beam_model(
            longueur=4000,
            b=80,
            h=160,
            section="Rectangulaire",
            classe="C24",
            nbr_appuis=3,
            inclinaison=90,
        )

        n1 = wb.get_node("N1")
        n2 = wb.get_node("N2")

        assert n1["X"].value == 0 and n1["Y"].value == 0
        assert n2["X"].value == 0
        assert n2["Y"].value > 0
