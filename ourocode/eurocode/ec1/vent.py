#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import logging
import os
import sys
from PIL import Image

import math as mt

import pandas as pd
import forallpeople as si
from ourocode.eurocode.core._renderer import handcalc
from ourocode.eurocode.mixins.math_utils import MathUtilsMixin

logger = logging.getLogger(__name__)

# Correspondance code département INSEE -> territoire d'outre-mer
DEPT_DOM = {"971": "Guadeloupe", "972": "Martinique", "973": "Guyane", "974": "Réunion", "976": "Mayotte"}

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# from eurocode.A0_Projet import Batiment
from ourocode.eurocode.core.batiment import Batiment

si.environment("structural")


class Vent(Batiment):
    RHO_AIR = 1.225 * si.kg / si.m**3
    CPI = [0.2, -0.3]
    VB0 = {
        "1": 22,
        "2": 24,
        "3": 26,
        "4": 28,
        "Guadeloupe": 36,
        "Guyane": 17,
        "Martinique": 32,
        "Mayotte": 30,
        "Réunion": 34,
    }

    CAT_TERRAIN = {
        "0": {"Z0": 0.005, "Zmin": 1},
        "II": {"Z0": 0.05, "Zmin": 2},
        "IIIa": {"Z0": 0.2, "Zmin": 5},
        "IIIb": {"Z0": 0.5, "Zmin": 9},
        "IV": {"Z0": 1, "Zmin": 15},
    }

    CAT_ORO = {"Aucun": 0, "Cas 1": 1, "Cas 2": 2}

    CFR = {"Lisse": 0.01, "Rugueuse": 0.02, "Très rugueuse": 0.04}

    def __init__(self, z: si.m, terrain: str = CAT_TERRAIN, oro: str = CAT_ORO, CsCd: float = 1, **kwargs):
        """Initialise le calcul de l'action du vent selon l'EN 1991-1-4 et son Annexe Nationale française.

        Hérite de Batiment. La zone de vent est déterminée automatiquement via le code_INSEE du projet.

        Args:
            z (float): Hauteur de référence z_e en m selon EN 1991-1-4 §7.2.2
                (hauteur à laquelle est étudiée la pression de vent).
            terrain (str): Catégorie de rugosité du terrain selon EN 1991-1-4 §4.3.2 :
                "0", "II", "IIIa", "IIIb" ou "IV". Defaults to "II".
            oro (str): Catégorie orographique selon AN français :
                - "Aucun" : pas d'effet orographique (C_0 = 1).
                - "Cas 1" : obstacles de hauteurs variées (calcul interactif des altitudes environnantes).
                - "Cas 2" : non implémenté (levée d'une ValueError).
                Defaults to "Aucun".
            CsCd (float): Coefficient structurel selon EN 1991-1-4 §6. Defaults to 1.
            **kwargs: Arguments transmis à la classe parent Batiment.
        """
        super().__init__(**kwargs)
        self.terrain = self._check_choice(terrain, self.CAT_TERRAIN.keys(), "terrain")
        self.oro = self._check_choice(oro, self.CAT_ORO.keys(), "oro")
        if self.oro == "Cas 1":
            self._calc_delta_AC()
        elif self.oro == "Cas 2":
            raise ValueError("Cas 2 non traité dans ce logiciel !")
        if CsCd <= 0:
            raise ValueError(f"CsCd doit être strictement positif, reçu : {CsCd}")
        self.z = z * si.m
        self.CsCd = CsCd

    @property
    def cat_terrain(self):
        return self.CAT_TERRAIN[self.terrain]

    @property
    def cat_oro(self):
        return self.CAT_ORO[self.oro]

    def _calc_delta_AC(self):
        """
        Calcule le delta AC en fonction des altitudes des obstacles environnants.
        """
        if not hasattr(self, "dict_alti_oro"):
            self.dict_alti_oro = {
                "An1": "",
                "An2": "",
                "Ae1": "",
                "Ae2": "",
                "As1": "",
                "As2": "",
                "Ao1": "",
                "Ao2": "",
            }
            # Demande des altitudes via des boîtes de dialogue QInputDialog (PySide6), autonome si aucune QApplication n'existe
            from PySide6.QtWidgets import QApplication, QInputDialog
            from PySide6.QtCore import Qt

            app = QApplication.instance()
            owns_app = False
            if app is None:
                app = QApplication(sys.argv)
                owns_app = True
            try:
                for cle, value in self.dict_alti_oro.items():
                    if cle.endswith("1"):
                        dist_msg = "à une distance de 500m"
                    else:
                        dist_msg = "à une distance de 1000m"
                    val, ok = QInputDialog.getInt(
                        None,
                        "Altitude orographique",
                        f"Altitude {cle} en m {dist_msg}:",
                        0,
                        -300,
                        3000,
                        1,
                        flags=Qt.WindowSystemMenuHint | Qt.WindowTitleHint,
                    )
                    if not ok:
                        raise RuntimeError("Saisie des altitudes annulée par l'utilisateur.")
                    self.dict_alti_oro[cle] = int(val) * si.m
            finally:
                if owns_app:
                    app.quit()

        Am = (
            2 * self.alt
            + self.dict_alti_oro["An1"]
            + self.dict_alti_oro["An2"]
            + self.dict_alti_oro["Ae1"]
            + self.dict_alti_oro["Ae2"]
            + self.dict_alti_oro["As1"]
            + self.dict_alti_oro["As2"]
            + self.dict_alti_oro["Ao1"]
            + self.dict_alti_oro["Ao2"]
        ) / 10

        self.delta_AC = self.alt - Am
        return self.delta_AC

    @property
    def zone_vent(self):
        """Zone de vent du projet : 1 à 4 en métropole, nom du territoire pour les DOM.

        Le fichier carte_action_region.csv code les communes d'outre-mer avec une zone
        "2.5" ou vide : le département INSEE (971-976) est alors utilisé pour déterminer
        la vitesse de référence cyclonique de l'Annexe Nationale.
        """
        file = "carte_action_region.csv"
        df = self._data_from_csv(file, index_col=1)
        row = df.loc[str(self.code_INSEE)]
        zone = row["Zone_vent"]
        try:
            zone_num = float(zone)
            if zone_num != int(zone_num):
                raise ValueError  # zone non entière ("2.5") : commune d'outre-mer
            return int(zone_num)
        except (TypeError, ValueError):
            dept = str(row["INSEE_Département"]).zfill(3)
            if dept in DEPT_DOM:
                return DEPT_DOM[dept]
            raise ValueError(
                f"Aucun zonage vent disponible pour la commune INSEE {self.code_INSEE} "
                f"(Zone_vent={zone!r})."
            )

    @property
    def Vb_0(self):
        """Retourne la vitesse de base de référence V_b,0 en m/s selon EN 1991-1-4 AN §4.2.

        Valeur déterminée par la zone de vent du projet (code_INSEE → zone 1 à 4, ou DOM).

        Returns:
            forallpeople.Physical: V_b,0 en m/s.
        """
        return self.VB0[str(self.zone_vent)] * si.m / si.s

    @property
    def Vb(self):
        """Retourne la vitesse de référence du vent V_b en m/s selon EN 1991-1-4 §4.2(2).

        Formule : V_b = C_dir × C_season × V_b,0.
        C_dir et C_season sont pris égaux à 1 (valeurs conservatrices de l'AN français).

        Returns:
            tuple: (latex_string, V_b) où V_b est la vitesse de référence en m/s (avec unité si.m/si.s).
        """
        V_b_0 = self.Vb_0
        Cdir = 1
        Cseason = 1

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            V_b = V_b_0 * Cdir * Cseason
            return V_b

        return val()

    @property
    def rayon_secteur_angu(self):
        """Retourne le rayon du secteur angulaire R à l'intérieur duquel la rugosité est à qualifier selon AN français.

        Formule : R = max(23 × z^1.2, 300) en m.

        Returns:
            tuple: (latex_string, R) où R est le rayon du secteur angulaire en m (avec unité si.m).
        """
        z = self.z.value

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            rayon_secteur_angu = max(23 * z**1.2, 300) * si.m
            return rayon_secteur_angu

        return val()

    @property
    def _Kr(self) -> float:
        z0II = self.CAT_TERRAIN["II"]["Z0"]
        return 0.19 * (self.cat_terrain["Z0"] / z0II) ** 0.07

    @property
    def _Cr_z(self):
        """Retourne le coefficient de rugosité, c r ( z ), tient compte de la variabilité de la vitesse moyenne du vent sur le site de la
        construction

        Returns:
                float: le coef de rugosité fonction de la hauteur au niveau du sol et de la rugosité de terrain
        """

        if self.cat_terrain["Zmin"] <= self.z.value <= 200:
            Cr_z = self._Kr * mt.log(self.z.value / self.cat_terrain["Z0"])
        elif self.z.value < self.cat_terrain["Zmin"]:
            Cr_z = self._Kr * mt.log(self.cat_terrain["Zmin"] / self.cat_terrain["Z0"])
        else:
            raise ValueError(
                "Erreur la hauteur du bâtiment ne peut dépasser 200m, impossible de calculer "
                "le coefficient de rugosité Cr"
            )

        # print("Kr :", self._Kr)
        # print("Cr_z :", Cr_z)
        return Cr_z

    @property
    def Co_z(self):
        """Retourne le coefficient orographique C_0(z) selon l'AN français de l'EN 1991-1-4 §4.3.3.

        - "Aucun" : C_0 = 1 (pas d'effet orographique).
        - "Cas 1" : C_0 calculé en fonction de delta_AC et de la hauteur z.

        Returns:
            tuple: (latex_string, C_0_z) où C_0_z est le coefficient orographique (sans unité, ≥ 1).
        """
        match self.cat_oro:
            case 1:
                Delta_A_C = self.delta_AC.value
                z_moins_10 = max(0.0, self.z.value - 10)

                @handcalc(
                    override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]"
                )
                def val():
                    C0_z = 1 + 0.004 * Delta_A_C * mt.exp(-0.014 * z_moins_10)  # 4.14-NA
                    res = max(C0_z, 1)
                    return res

            case 2:
                raise ValueError("Erreur la catégorie orographique 2 n'est pas encore implémentée")
            case _:

                @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
                def val():
                    C0_z = 1
                    return C0_z

        result = val()
        if result[1] > 1.15:
            logger.warning(
                "C0(z) = %.3f > 1.15 : la formule de l'AN (orographie complexe) n'est plus "
                "valable, une étude détaillée du relief est nécessaire (EN 1991-1-4 §A.3).",
                result[1],
            )
        return result

    @property
    # @handcalc(override="params", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
    def Vm_z(self):
        """Retourne la vitesse moyenne du vent V_m(z) en m/s selon EN 1991-1-4 §4.3.1.

        Formule : V_m(z) = C_r(z) × C_0(z) × V_b.

        Returns:
            tuple: (latex_string, V_m_z) où V_m_z est la vitesse moyenne en m/s (avec unité si.m/si.s).
        """
        C_r_z = self._Cr_z
        C_o_z = self.Co_z[1]
        V_b = self.Vb[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            V_m_z = C_r_z * C_o_z * V_b
            return V_m_z

        return val()

    @property
    def _Kl(self):
        """Retourne le coefficient de turbulence

        Returns:
                float: coef de turbulence
        """
        match self.cat_oro:
            case 1:
                return self.Co_z[1] * (1 - 2 * 10**-4 * (mt.log10(self.cat_terrain["Z0"]) + 3) ** 6)  # 4.19NA
            case _:
                return 1 - 2 * 10**-4 * (mt.log10(self.cat_terrain["Z0"]) + 3) ** 6  # 4.20NA

    @property
    def _sigma_v(self):
        """Retourne l'écart type de turbulence

        Returns:
                float: écart type
        """
        return self._Kr * self.Vb[1] * self._Kl  # 4.6

    @property
    def _Iv_z(self):
        """Retourne l'intensité de turbulence à la hauteur z

        Returns:
                float: intensité de la turbulence
        """
        if self.cat_terrain["Zmin"] <= self.z.value <= 200 or self.z.value < self.cat_terrain["Zmin"]:
            return self._sigma_v / self.Vm_z[1]  # 4.7
        else:
            raise ValueError(
                "Erreur la hauteur max du bâtiment ne peut dépasser 200m, impossible de calculer l'intensité de turbulence Iv"
            )

    @property
    def Qb(self):
        """Retourne la pression dynamique de référence Q_b en N/m² selon EN 1991-1-4 §4.5(1) (valeur informative).

        Formule : Q_b = 0.5 × ρ × V_b².

        Returns:
            tuple: (latex_string, Q_b) où Q_b est la pression dynamique de référence en N/m² (avec unité si.Pa).
        """
        V_b = self.Vb[1]
        rho_air = self.RHO_AIR

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            Q_b = 0.5 * rho_air * V_b**2  # 4.10
            return Q_b

        return val()

    @property
    def Qp_z(self):
        """Retourne la pression dynamique de pointe Q_p(z) en N/m² selon EN 1991-1-4 §4.5(1).

        Formule : Q_p(z) = (1 + 7 × I_v(z)) × 0.5 × ρ × V_m(z)².
        Tient compte de la turbulence via I_v(z).

        Returns:
            tuple: (latex_string, Q_p_z) où Q_p_z est la pression de pointe en N/m² (avec unité si.Pa).
        """
        I_v_z = self._Iv_z
        rho_air = self.RHO_AIR
        V_m_z = self.Vm_z[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            Q_p_z = (1 + 7 * I_v_z) * 0.5 * rho_air * V_m_z**2  # 4.8
            return Q_p_z

        return val()

    @property
    def Ce_z(self):
        """Retourne le coefficient d'exposition C_e(z) selon EN 1991-1-4 §4.5(2) (valeur informative).

        Formule : C_e(z) = Q_p(z) / Q_b.

        Returns:
            tuple: (latex_string, C_e_z) où C_e_z est le coefficient d'exposition (sans unité).
        """
        Q_p_z = self.Qp_z[1]
        Q_b = self.Qb[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            Ce_z = Q_p_z / Q_b  # 4.9
            return Ce_z

        return val()

    def We(self, Cpe: float):
        """Retourne la pression aérodynamique extérieure W_e en N/m² selon EN 1991-1-4 §5.2(1).

        Formule : W_e = Q_p(z) × C_pe.

        Args:
            Cpe (float): Coefficient de pression extérieure (C_pe,10 ou C_pe,1 selon l'aire chargée).

        Returns:
            tuple: (latex_string, W_e) où W_e est la pression extérieure en N/m² (avec unité si.Pa).
        """
        Q_p_z = self.Qp_z[1]
        C_pe = Cpe

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            W_e = Q_p_z * C_pe  # 5.1
            return W_e

        return val()

    def Wi(self, Cpi: float):
        """Retourne la pression aérodynamique intérieure W_i en N/m² selon EN 1991-1-4 §5.2(2).

        Formule : W_i = Q_p(z) × C_pi.
        Les valeurs standard de C_pi sont +0.2 (pression) et -0.3 (dépression) selon CPI.

        Args:
            Cpi (float): Coefficient de pression intérieure (C_pi).

        Returns:
            tuple: (latex_string, W_i) où W_i est la pression intérieure en N/m² (avec unité si.Pa).
        """
        Q_p_z = self.Qp_z[1]
        C_pi = Cpi

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            W_i = Q_p_z * C_pi  # 5.2
            return W_i

        return val()

    def Fw_e(self, Aref: float):
        pass  # 5.5

    def Fw_i(self, Aref: float):
        pass  # 5.6

    def Ffr(self, Afr: float, Cfr: str = CFR):
        """Retourne les forces de frottement sur le bâtiment en N selon l'EN 1991-1-4 §5.7/7.5.

        Args:
                Afr (float): aire de référence en m² correspondant à la surface d'application des forces de frottement.
                                        Il convient d'appliquer les forces de frottement sur la partie des surfaces extérieures
                                        parallèle au vent, située au-delà d'une certaine distance des bords au vent ou des angles
                                        au vent de la toiture, distance égale à la plus petite valeur de 2 · b ou 4 · h .

                Cfr (str): type de surface de frottement
                                                                        - Lisse / exemple: acier, béton lisse
                                                                        - Rugueuse / exemple: béton brut, bardeaux bitumés(shingles)
                                                                        - Très rugueuse / exemple: ondulation, nervures, pliures

        Returns:
                float: effort en N
        """
        A_fr = Afr * si.m**2
        C_fr = self.CFR[Cfr]
        Q_p_z = self.Qp_z[1]

        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY, left="\\[", right="\\]")
        def val():
            F_fr = C_fr * Q_p_z * A_fr
            return F_fr

        return val()

    def show_Ffr(self):
        """Affiche l'image du zonage pour un effet d'entrainement en toiture"""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Ffr.png")
        image = Image.open(file)
        image.show()

    def K_red_U(self):
        """Calcule et stocke le coefficient de défaut de corrélation K_red selon EN 1991-1-4 §7.2.2(3).

        Applicable uniquement aux faces D et E des murs verticaux.
        Formule dépendante de h_bat/d_bat :
        - h/d ≥ 5 : K_red = 1.0.
        - h/d ≤ 1 : K_red = 0.85.
        - 1 < h/d < 5 : interpolation linéaire.
        Le résultat est stocké dans self.Kred_U.
        """
        if self.h_bat / self.d_bat >= 5:
            self.Kred_U = 1
        elif self.h_bat / self.d_bat <= 1:
            self.Kred_U = 0.85
        else:
            self.Kred_U = 0.85 + 0.15 * (((self.h_bat / self.d_bat) - 1) / 4)


class Murs_verticaux(Vent):
    def __init__(self, load_area: float, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur les murs verticaux selon EN 1991-1-4 §7.2.2.

        Le calcul est effectué automatiquement pour les deux directions de vent (0° et 90°).
        Les zones sont A, B, C, D, E déterminées à partir de e = min(b, 2h) pour chaque direction.
        Les C_pe sont interpolés selon le rapport h/d et l'aire chargée.

        Args:
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
                Si 1 < load_area < 10, interpolation logarithmique entre C_pe,1 et C_pe,10.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        self.load_area = load_area
        self._wind_direction = {"0°": {}, "90°": {}}

        for cle in self._wind_direction.keys():
            match cle:
                case "0°":
                    self.e = min(self.b_bat, self.h_bat * 2)
                case "90°":
                    self.e = min(self.d_bat, self.h_bat * 2)
                    self.b_bat, self.d_bat = self.d_bat, self.b_bat

            self._wind_direction[cle]["geometrie"] = self._geo()
            self._wind_direction[cle]["Cpe"] = self._Cpe(cle)

            if cle == "90°":
                self.b_bat, self.d_bat = self.d_bat, self.b_bat

    def _geo(self):
        """Calcul les dimensions du zonage des parois verticale"""
        self._df = self._data_from_csv(os.path.join("vent", "vent_Cpe_mur_verticaux.csv"))
        self._df.reset_index(drop=False, inplace=True)

        if self.e < self.d_bat:
            geometrie = {"A": self.e / 5, "B": 0.8 * self.e, "C": self.d_bat - self.e}

        elif self.e >= 5 * self.d_bat:
            geometrie = {"A": self.d_bat}
        else:
            geometrie = {"A": self.e / 5, "B": self.d_bat - self.e / 5}
        geometrie["D"] = self.b_bat
        geometrie["E"] = self.b_bat

        return geometrie

    def _Cpe(self, cle: str):
        """Calcul les Cpe des parois verticales

        Args:
                cle (str): sens du vent sur le bâtiment 0° ou 90°
        """
        h_d = self.h_bat / self.d_bat
        list_CPE_line = []
        borne_inf = None
        if h_d > 0.25 and h_d < 1:
            df_index = [3, 4]
            borne_inf = 0.25
            borne_sup = 1

        elif h_d > 1 and h_d < 5:
            df_index = [2, 3]
            borne_inf = 1
            borne_sup = 5

        elif h_d <= 0.25:
            self._df = self._df[["Zone", "Cpe", "0.25"]]

        elif h_d >= 5:
            self._df = self._df[["Zone", "Cpe", "5"]]

        else:
            self._df = self._df[["Zone", "Cpe", "1"]]

        if borne_inf:
            for i in range(self._df.shape[0]):
                list_CPE_line.append(
                    MathUtilsMixin.interpolation_lineaire(
                        h_d,
                        borne_inf,
                        borne_sup,
                        self._df.iloc[i, df_index[1]],
                        self._df.iloc[i, df_index[0]],
                    )
                )
            self._df[round(h_d, 3)] = list_CPE_line
            self._df = self._df[["Zone", "Cpe", round(h_d, 3)]]

        if self.load_area > 1 and self.load_area < 10:
            for i in range(self._df.shape[0]):
                if i % 2 == 0:
                    self._df.loc[self._df.shape[0]] = [
                        self._df.iloc[i, 0],
                        "CPE " + str(round(self.load_area, 2)),
                        MathUtilsMixin.interpolation_logarithmique(
                            self.load_area,
                            1,
                            10,
                            self._df.iloc[i + 1, 2],
                            self._df.iloc[i, 2],
                        ),
                    ]

        self._df.set_index("Zone", inplace=True)
        self._df = self._df.loc[[zone for zone in self._wind_direction[cle]["geometrie"].keys()]]
        return self._df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self, dir=("0°", "90°")):
        """Retourne les dimensions des zones de pression pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            dict: {"A": ..., "B": ..., "C": ..., "D": ..., "E": ...} avec les dimensions en m.
        """
        return self._wind_direction[dir]["geometrie"]

    def get_Cpe(self, dir=("0°", "90°")):
        """Retourne les coefficients de pression extérieure C_pe pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            pandas.DataFrame: C_pe par zone (A, B, C, D, E), interpolés selon h/d et load_area.
        """
        return self._wind_direction[dir]["Cpe"]

    def show_zonage(self):
        """Affiche l'image du zonage pour les murs verticaux selon EN 1991-1-4 §7.2.2."""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cpe_mur_verticaux.png")
        image = Image.open(file)
        image.show()


class Toiture_terrasse_acrotere(Vent):
    def __init__(self, load_area: float, hp: si.m, h: si.m, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture terrasse avec acrotère selon EN 1991-1-4 §7.2.3.

        Le calcul est effectué pour la direction de vent 0° uniquement (symétrie).
        Les zones de pression sont F, G, H, I déterminées à partir de e = min(b, 2h).
        Les C_pe sont interpolés selon le rapport hp/h (plafonné à 0.1) et l'aire chargée.

        Args:
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
                Si 1 < load_area < 10, interpolation logarithmique entre C_pe,1 et C_pe,10.
            hp (float): Hauteur de l'acrotère en m.
            h (float): Hauteur du bâtiment sous l'acrotère en m.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        # hp = 0 correspond à une terrasse à rives vives (sans acrotère, tab. 7.2)
        self.type_terrasse = "Acrotères" if hp > 0 else "Rives vives"
        self.load_area = load_area
        self.hp = hp * si.m
        self.h = h * si.m
        self._hp_h = self.hp / self.h
        self._wind_direction = {"0°": {}}

        self.e = min(self.b_bat, self.h_bat * 2)

        self._wind_direction["0°"]["geometrie"] = self._geo()
        self._wind_direction["0°"]["Cpe"] = self._Cpe()

    def _geo(self):
        """Calcul les surfaces du zonage de la toiture"""
        geometrie = {
            "F": {
                "Longueur": self.e / 4,
                "Largeur": self.e / 10,
                "Surface": (self.e / 4) * (self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
            "G": {
                "Longueur": self.b_bat - (self.e / 4 * 2),
                "Largeur": self.e / 10,
                "Surface": (self.b_bat - (self.e / 4 * 2)) * (self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
            "H": {
                "Longueur": self.b_bat,
                "Largeur": self.e / 2 - self.e / 10,
                "Surface": self.b_bat * (self.e / 2 - self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
            "I": {
                "Longueur": self.b_bat,
                "Largeur": self.d_bat - self.e / 2,
                "Surface": self.b_bat * (self.d_bat - self.e / 2) / mt.cos(mt.radians(self.alpha_toit)),
            },
        }
        return geometrie

    def _Cpe(self):
        """Calcul les Cpe d'une toiture terrasse avec acrotère

        Args:
                direction (str): sens du vent sur le bâtiment 0° ou 90°
        """
        num_columns = 7
        df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_terrasse.csv"))
        name_csv_type = self.type_terrasse.replace("è", "e")
        df = df.loc[name_csv_type]
        df.reset_index(drop=False, inplace=True)
        df.drop("Type", axis=1, inplace=True)
        list_type_rapport = df["Rapport"].dropna().unique()

        # Pas de rapport hp/h pour les rives vives (bord sans acrotère)
        if list_type_rapport.size > 0:
            self._hp_h = min(self._hp_h, 0.1)
            if self._hp_h < list_type_rapport.min():
                # hp/h sous le plus petit rapport tabulé (0.025) : on retient la borne basse
                df = df[df["Rapport"] == list_type_rapport.min()]
                self._hp_h = list_type_rapport.min()
            elif self._hp_h not in list_type_rapport:
                minimum = list_type_rapport[list_type_rapport < self._hp_h].max()
                maximum = list_type_rapport[list_type_rapport > self._hp_h].min()
                df_min = df[df["Rapport"] == minimum]
                df_max = df[df["Rapport"] == maximum]

                df = df[df["Rapport"].isin([minimum, maximum])]
                df.reset_index(drop=True, inplace=True)

                for i, cpe in enumerate(["CPE 10", "CPE 1"]):
                    row = [round(self._hp_h, 3), cpe]
                    for j in range(2, num_columns):
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                self._hp_h,
                                minimum,
                                maximum,
                                df_min.iloc[i, j],
                                df_max.iloc[i, j],
                            )
                        )
                    df.loc[df.shape[0]] = row

        if self.load_area > 1 and self.load_area < 10:
            for i in range(df.shape[0]):
                if i % 2 == 0:
                    row = [df.iloc[i, 0], "CPE " + str(round(self.load_area, 2))]
                    for j in range(2, num_columns):
                        row.append(
                            MathUtilsMixin.interpolation_logarithmique(
                                self.load_area, 1, 10, df.iloc[i + 1, j], df.iloc[i, j]
                            )
                        )
                    df.loc[df.shape[0]] = row
            df.reset_index(drop=True, inplace=True)

        if list_type_rapport.size > 0:
            df = df[df["Rapport"] == round(self._hp_h, 3)]
            df.dropna(axis=1, how="all", inplace=True)
            df.set_index("Rapport", inplace=True)
        else:
            df.drop("Rapport", axis=1, inplace=True)
            df.dropna(axis=1, how="all", inplace=True)
            df.set_index("Cpe", inplace=True)
        return df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self):
        """Retourne les dimensions des zones de pression pour la direction de vent 0°.

        Returns:
            dict: {"F": ..., "G": ..., "H": ..., "I": ...} avec Longueur, Largeur et Surface en m / m².
        """
        return self._wind_direction["0°"]["geometrie"]

    def get_Cpe(self):
        """Retourne les coefficients de pression extérieure C_pe pour la direction de vent 0°.

        Returns:
            pandas.DataFrame: C_pe par zone (F, G, H, I), interpolés selon hp/h et load_area.
        """
        return self._wind_direction["0°"]["Cpe"]

    def show_zonage(self):
        """Affiche l'image du zonage pour une toiture terrasse avec acrotère selon EN 1991-1-4 §7.2.3."""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cpe_toiture_terrasse.png")
        image = Image.open(file)
        image.show()


class Toiture_1_pant(Vent):
    def __init__(self, load_area: float, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture à un versant selon EN 1991-1-4 §7.2.4.

        Le calcul est effectué pour les directions de vent 0°, 90° et 180°.
        Les zones de pression sont F, G, H (directions 0° et 180°) et F_up/F_low, G, H, I (direction 90°).
        Les C_pe sont interpolés en fonction de alpha_toit et de l'aire chargée.

        Note:
            Cette classe ne prend en compte qu'une seule valeur d'alpha_toit.

        Args:
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
                Si 1 < load_area < 10, interpolation logarithmique entre C_pe,1 et C_pe,10.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        self.load_area = load_area
        self._wind_direction = {"0°": {}, "180°": {}, "90°": {}}

        for cle in self._wind_direction.keys():
            match cle:
                case "0°" | "180°":
                    self.e = min(self.b_bat, self.h_bat * 2)

                case "90°":
                    self.e = min(self.d_bat, self.h_bat * 2)
                    self.b_bat, self.d_bat = self.d_bat, self.b_bat

            self._wind_direction[cle]["geometrie"] = self._geo(cle)
            self._wind_direction[cle]["Cpe"] = self._Cpe(cle)

            if cle == "90°":
                self.b_bat, self.d_bat = self.d_bat, self.b_bat

    def _geo(self, direction: str):
        """Calcul les surfaces du zonage de la toiture
        Attention on considère le faîtage centrée et les pentes égales entre les deux pants
        """
        match direction:
            case "0°" | "180°":
                geometrie = {
                    "F": {
                        "Longueur": self.e / 4,
                        "Largeur": self.e / 10,
                        "Surface": (self.e / 4) * (self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "G": {
                        "Longueur": self.b_bat - (self.e / 4 * 2),
                        "Largeur": self.e / 10,
                        "Surface": (self.b_bat - (self.e / 4 * 2))
                        * (self.e / 10)
                        / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "H": {
                        "Longueur": self.b_bat,
                        "Largeur": self.d_bat - self.e / 10,
                        "Surface": self.b_bat * (self.d_bat - self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                }
            case "90°":
                # Fig. 7.7(c) : F_up et F_low (e/4) aux deux extrémités du bandeau au vent (e/10),
                # G entre les deux, H jusqu'à e/2 puis I jusqu'au bord sous le vent
                geometrie = {
                    "F_up / F_low": {
                        "Longueur": self.e / 4,
                        "Largeur": self.e / 10,
                        "Surface": (self.e / 10) * (self.e / 4) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "G": {
                        "Longueur": self.b_bat - self.e / 2,
                        "Largeur": self.e / 10,
                        "Surface": (self.e / 10)
                        * (self.b_bat - self.e / 2)
                        / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "H": {
                        "Longueur": self.b_bat,
                        "Largeur": self.e / 2 - self.e / 10,
                        "Surface": (self.e / 2 - self.e / 10) * self.b_bat / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "I": {
                        "Longueur": self.b_bat,
                        "Largeur": self.d_bat - self.e / 2,
                        "Surface": (self.d_bat - self.e / 2) * self.b_bat / mt.cos(mt.radians(self.alpha_toit)),
                    },
                }
        return geometrie

    def _Cpe(self, direction: str):
        """Calcul les Cpe d'une toiture 2 versants

        Args:
                direction (str): sens du vent sur le bâtiment 0° ou 90°
        """
        match direction:
            case "0°":
                num_columns = 8
                df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_1_versant_0_degres.csv"))
            case "90°":
                num_columns = 7
                df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_1_versant_90_degres.csv"))
            case "180°":
                num_columns = 5
                df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_1_versant_180_degres.csv"))

        df.reset_index(drop=False, inplace=True)
        list_alpha_toit = df["alpha_toit"].unique()
        touch_min = self.alpha_toit < list_alpha_toit.min()
        touch_max = self.alpha_toit > list_alpha_toit.max()

        if self.alpha_toit not in list_alpha_toit:
            if touch_min:
                minimum = list_alpha_toit.min()
            else:
                minimum = list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
            if touch_max:
                maximum = list_alpha_toit.max()
            else:
                maximum = list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
            df_min = df[df["alpha_toit"] == minimum]
            df_max = df[df["alpha_toit"] == maximum]

            df = df[df["alpha_toit"].isin([minimum, maximum])]
            df.reset_index(drop=True, inplace=True)

            df.reset_index(drop=True, inplace=True)
            for i, cpe in enumerate(["CPE 10", "CPE 1"]):
                row = [round(self.alpha_toit, 2), cpe]
                for j in range(2, num_columns):
                    if touch_min:
                        row.append(df_min.iloc[i, j])
                    elif touch_max:
                        row.append(df_max.iloc[i, j])
                    else:
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                self.alpha_toit,
                                minimum,
                                maximum,
                                df_min.iloc[i, j],
                                df_max.iloc[i, j],
                            )
                        )
                df.loc[df.shape[0]] = row

        if self.load_area > 1 and self.load_area < 10:
            for i in range(df.shape[0]):
                if i % 2 == 0:
                    row = [df.iloc[i, 0], "CPE " + str(round(self.load_area, 2))]
                    for j in range(2, num_columns):
                        row.append(
                            MathUtilsMixin.interpolation_logarithmique(
                                self.load_area, 1, 10, df.iloc[i + 1, j], df.iloc[i, j]
                            )
                        )
                    df.loc[df.shape[0]] = row
            df.reset_index(drop=True, inplace=True)

        df = df[df["alpha_toit"] == round(self.alpha_toit, 2)]
        df.dropna(axis=1, how="all", inplace=True)
        df.set_index("alpha_toit", inplace=True)
        return df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self, direction=("0°", "90°", "180°")):
        """Retourne les dimensions des zones de pression pour la direction de vent donnée.

        Args:
            direction (str): Direction du vent ("0°", "90°" ou "180°"). Defaults to ("0°", "90°", "180°").

        Returns:
            dict: Zones avec Longueur, Largeur et Surface en m / m².
                - 0° / 180° : zones F, G, H.
                - 90° : zones F_up/F_low, G, H, I.
        """
        return self._wind_direction[direction]["geometrie"]

    def get_Cpe(self, direction=("0°", "90°", "180°")):
        """Retourne les coefficients de pression extérieure C_pe pour la direction de vent donnée.

        Args:
            direction (str): Direction du vent ("0°", "90°" ou "180°"). Defaults to ("0°", "90°", "180°").

        Returns:
            pandas.DataFrame: C_pe par zone, interpolés selon alpha_toit et load_area.
        """
        return self._wind_direction[direction]["Cpe"]

    def show_zonage(self):
        """Affiche l'image du zonage pour une toiture à 1 versant selon EN 1991-1-4 §7.2.4."""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cpe_toiture_1_versant.png")
        image = Image.open(file)
        image.show()


class Toiture_2_pants(Vent):
    def __init__(self, load_area: float, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture à deux versants selon EN 1991-1-4 §7.2.5.

        Le calcul est effectué pour les directions de vent 0° et 90°.
        Les zones de pression sont F, G, H, I, J (direction 0°) et F, G, H, I (direction 90°).
        Les C_pe sont interpolés en fonction de alpha_toit et de l'aire chargée.

        Note:
            Cette classe suppose des pentes identiques sur les deux versants et un faîtage centré.
            Un seul alpha_toit est pris en compte.

        Args:
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
                Si 1 < load_area < 10, interpolation logarithmique entre C_pe,1 et C_pe,10.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        self.load_area = load_area
        self._wind_direction = {"0°": {}, "90°": {}}

        for cle in self._wind_direction.keys():
            match cle:
                case "0°":
                    self.e = min(self.b_bat, self.h_bat * 2)

                case "90°":
                    self.e = min(self.d_bat, self.h_bat * 2)
                    self.b_bat, self.d_bat = self.d_bat, self.b_bat

            self._wind_direction[cle]["geometrie"] = self._geo(cle)
            self._wind_direction[cle]["Cpe"] = self._Cpe(cle)

            if cle == "90°":
                self.b_bat, self.d_bat = self.d_bat, self.b_bat

    def _geo(self, direction: str = "0°"):
        """Calcul les surfaces du zonage de la toiture
        Attention on considère le faîtage centrée et les pentes égales entre les deux pants
        """
        match direction:
            case "0°":
                geometrie = {
                    "F": {
                        "Longueur": self.e / 4,
                        "Largeur": self.e / 10,
                        "Surface": (self.e / 4) * (self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "G": {
                        "Longueur": self.b_bat - (self.e / 4 * 2),
                        "Largeur": self.e / 10,
                        "Surface": (self.b_bat - (self.e / 4 * 2))
                        * (self.e / 10)
                        / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "H": {
                        "Longueur": self.b_bat,
                        "Largeur": self.d_bat / 2 - self.e / 10,
                        "Surface": self.b_bat * (self.d_bat / 2 - self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "I": {
                        "Longueur": self.b_bat,
                        "Largeur": self.d_bat / 2 - self.e / 10,
                        "Surface": self.b_bat * (self.d_bat / 2 - self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "J": {
                        "Longueur": self.b_bat,
                        "Largeur": self.e / 10,
                        "Surface": self.b_bat * (self.e / 10) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                }
            case "90°":
                # Fig. 7.8(c) : sur chaque versant, F (e/4) au pignon et G jusqu'au faîtage
                geometrie = {
                    "F": {
                        "Longueur": self.e / 4,
                        "Largeur": self.e / 10,
                        "Surface": (self.e / 10) * (self.e / 4) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "G": {
                        "Longueur": self.b_bat / 2 - self.e / 4,
                        "Largeur": self.e / 10,
                        "Surface": (self.e / 10)
                        * (self.b_bat / 2 - self.e / 4)
                        / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "H": {
                        "Longueur": self.b_bat / 2,
                        "Largeur": self.e / 2 - self.e / 10,
                        "Surface": (self.e / 2 - self.e / 10) * (self.b_bat / 2) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                    "I": {
                        "Longueur": self.b_bat / 2,
                        "Largeur": self.d_bat - self.e / 2,
                        "Surface": (self.d_bat - self.e / 2) * (self.b_bat / 2) / mt.cos(mt.radians(self.alpha_toit)),
                    },
                }
        return geometrie

    def _Cpe(self, direction: str):
        """Calcul les Cpe d'une toiture 2 versants

        Args:
                direction (str): sens du vent sur le bâtiment 0° ou 90°
        """
        match direction:
            case "0°":
                num_columns = 12
                df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_2_versants_0_degres.csv"))
            case "90°":
                num_columns = 6
                df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_2_versants_90_degres.csv"))

        df.reset_index(drop=False, inplace=True)
        list_alpha_toit = df["alpha_toit"].unique()
        touch_min = self.alpha_toit < list_alpha_toit.min()
        touch_max = self.alpha_toit > list_alpha_toit.max()

        if self.alpha_toit not in list_alpha_toit:
            if touch_min:
                minimum = list_alpha_toit.min()
            else:
                minimum = list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
            if touch_max:
                maximum = list_alpha_toit.max()
            else:
                maximum = list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
            df_min = df[df["alpha_toit"] == minimum]
            df_max = df[df["alpha_toit"] == maximum]

            df = df[df["alpha_toit"].isin([minimum, maximum])]
            df.reset_index(drop=True, inplace=True)

            df.reset_index(drop=True, inplace=True)
            for i, cpe in enumerate(["CPE 10", "CPE 1"]):
                row = [round(self.alpha_toit, 2), cpe]
                for j in range(2, num_columns):
                    if touch_min:
                        row.append(df_min.iloc[i, j])
                    elif touch_max:
                        row.append(df_max.iloc[i, j])
                    else:
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                self.alpha_toit,
                                minimum,
                                maximum,
                                df_min.iloc[i, j],
                                df_max.iloc[i, j],
                            )
                        )
                df.loc[df.shape[0]] = row

        if self.load_area > 1 and self.load_area < 10:
            for i in range(df.shape[0]):
                if i % 2 == 0:
                    row = [df.iloc[i, 0], "CPE " + str(round(self.load_area, 2))]
                    for j in range(2, num_columns):
                        row.append(
                            MathUtilsMixin.interpolation_logarithmique(
                                self.load_area, 1, 10, df.iloc[i + 1, j], df.iloc[i, j]
                            )
                        )
                    df.loc[df.shape[0]] = row
            df.reset_index(drop=True, inplace=True)

        df = df[df["alpha_toit"] == round(self.alpha_toit, 2)]
        df.dropna(axis=1, how="all", inplace=True)
        df.set_index("alpha_toit", inplace=True)
        return df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self, dir=("0°", "90°")):
        """Retourne les dimensions des zones de pression pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            dict: Zones avec Longueur, Largeur et Surface en m / m².
                - 0° : zones F, G, H, I, J.
                - 90° : zones F, G, H, I.
        """
        return self._wind_direction[dir]["geometrie"]

    def get_Cpe(self, dir=("0°", "90°")):
        """Retourne les coefficients de pression extérieure C_pe pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            pandas.DataFrame: C_pe par zone, interpolés selon alpha_toit et load_area.
        """
        return self._wind_direction[dir]["Cpe"]

    def show_zonage(self):
        """Affiche l'image du zonage pour une toiture à 2 versants selon EN 1991-1-4 §7.2.5."""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cpe_toiture_2_versants.png")
        image = Image.open(file)
        image.show()


class Toiture_4_pants(Vent):
    def __init__(self, load_area: float, alpha_toit_90: float = None, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture à quatre pans (croupes) selon EN 1991-1-4 §7.2.6.

        Le calcul est effectué pour les directions de vent 0° et 90° (fig. 7.9).
        Zones : F, G, H (versant au vent), I, J, K, L, M (croupes et versant sous le vent) à 0°,
        et F, G, H, I, J, M, N à 90°. Les C_pe (tab. 7.5) sont les mêmes pour les deux directions
        mais dépendent de la pente du versant au vent : alpha_toit à 0°, alpha_toit_90 à 90°.

        Hypothèses : toiture symétrique, croupes à 45° en plan (retrait de croupe = d_bat/2,
        faitière de longueur b_bat − d_bat, donc b_bat > d_bat exigé).

        Args:
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
            alpha_toit_90 (float): Pente du versant de croupe en degrés. Si None, = alpha_toit.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        self.load_area = load_area
        self.alpha_toit_90 = alpha_toit_90 if alpha_toit_90 is not None else self.alpha_toit
        self._wind_direction = {"0°": {}, "90°": {}}

        if self.b_bat <= self.d_bat:
            raise ValueError(
                "Toiture à quatre pans : la dimension le long de la faitière (b_bat) doit être "
                f"strictement supérieure à la largeur de pignon (d_bat), reçu b={self.b_bat}, d={self.d_bat}."
            )

        for cle in self._wind_direction.keys():
            match cle:
                case "0°":
                    self.e = min(self.b_bat, self.h_bat * 2)
                case "90°":
                    self.e = min(self.d_bat, self.h_bat * 2)
                    self.b_bat, self.d_bat = self.d_bat, self.b_bat

            self._wind_direction[cle]["geometrie"] = self._geo(cle)
            self._wind_direction[cle]["Cpe"] = self._Cpe(cle)

            if cle == "90°":
                self.b_bat, self.d_bat = self.d_bat, self.b_bat

    def _cos_alpha(self, direction: str) -> float:
        alpha = self.alpha_toit if direction == "0°" else self.alpha_toit_90
        return mt.cos(mt.radians(alpha))

    def _geo(self, direction: str = "0°"):
        """Calcule les surfaces des zones de la toiture à quatre pans (fig. 7.9).

        Surfaces en plan divisées par cos(pente du versant au vent). Croupes à 45° en plan :
        les zones L et J des croupes sont les bandes de profondeur e/10 le long des arêtiers,
        K la bande e/10 du versant sous le vent le long de la faitière.
        """
        e = self.e
        b = self.b_bat
        d = self.d_bat
        cos_a = self._cos_alpha(direction)
        # Bandes e/10 le long des arêtiers, plafonnées pour que L+J n'excèdent pas la croupe
        s_hip = min(e / 10, min(b, d) / 4)
        match direction:
            case "0°":
                long_fait = b - d
                geometrie = {
                    "F": {"Longueur": e / 4, "Largeur": e / 10, "Surface": (e / 4) * (e / 10) / cos_a},
                    "G": {"Longueur": b - e / 2, "Largeur": e / 10, "Surface": (b - e / 2) * (e / 10) / cos_a},
                    "H": {
                        "Longueur": (b - e / 5 + long_fait) / 2,
                        "Largeur": d / 2 - e / 10,
                        "Surface": (b - e / 5 + long_fait) / 2 * (d / 2 - e / 10) / cos_a,
                    },
                    "K": {
                        "Longueur": long_fait + e / 10,
                        "Largeur": e / 10,
                        "Surface": (long_fait + e / 10) * (e / 10) / cos_a,
                    },
                    "I": {
                        "Longueur": (b - d / 2 + e / 10),
                        "Largeur": d / 2 - e / 10,
                        "Surface": (b - d / 2 + e / 10) * (d / 2 - e / 10) / cos_a,
                    },
                    "L": {"Longueur": d / 2, "Largeur": s_hip, "Surface": s_hip * (d / 2) / cos_a},
                    "M": {
                        "Longueur": d / 2,
                        "Largeur": d / 2 - 2 * s_hip,
                        "Surface": (d**2 / 4 - s_hip * d) / cos_a,
                    },
                    "J": {"Longueur": d / 2, "Largeur": s_hip, "Surface": s_hip * (d / 2) / cos_a},
                }
            case "90°":
                # Versant de croupe au vent : F, G, H ; versants latéraux : M, N ; croupe sous le vent : J, I
                triangle = b**2 / 4
                bande = (b + b - 2 * e / 10) / 2 * (e / 10)
                aire_versant = (2 * d - b) * b / 4
                m_plan = min(e / 2, b / 2) ** 2 / 2 + (e / 2 - min(e / 2, b / 2)) * (b / 2)
                geometrie = {
                    "F": {"Longueur": e / 4, "Largeur": e / 10, "Surface": (e / 4) * (e / 10) / cos_a},
                    "G": {"Longueur": b - e / 2, "Largeur": e / 10, "Surface": (b - e / 2) * (e / 10) / cos_a},
                    "H": {
                        "Longueur": b / 2,
                        "Largeur": b / 2 - e / 10,
                        "Surface": (triangle - bande) / cos_a,
                    },
                    "M": {"Longueur": e / 2, "Largeur": b / 2, "Surface": m_plan / cos_a},
                    "N": {"Longueur": d - e / 2, "Largeur": b / 2, "Surface": (aire_versant - m_plan) / cos_a},
                    "J": {"Longueur": b / 2, "Largeur": s_hip, "Surface": 2 * s_hip * (b / 2) / cos_a},
                    "I": {
                        "Longueur": b / 2,
                        "Largeur": b / 2 - s_hip,
                        "Surface": (triangle - 2 * s_hip * (b / 2)) / cos_a,
                    },
                }
        return geometrie

    def _Cpe(self, direction: str):
        """Calcule les Cpe d'une toiture à quatre pans (tab. 7.5), interpolés sur la pente
        du versant au vent (alpha_toit à 0°, alpha_toit_90 à 90°) et l'aire chargée."""
        alpha = self.alpha_toit if direction == "0°" else self.alpha_toit_90
        num_columns = 14
        df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_4_pans.csv"))
        df.reset_index(drop=False, inplace=True)
        list_alpha_toit = df["alpha_toit"].unique()
        touch_min = alpha < list_alpha_toit.min()
        touch_max = alpha > list_alpha_toit.max()

        if alpha not in list_alpha_toit:
            minimum = list_alpha_toit.min() if touch_min else list_alpha_toit[list_alpha_toit < alpha].max()
            maximum = list_alpha_toit.max() if touch_max else list_alpha_toit[list_alpha_toit > alpha].min()
            df_min = df[df["alpha_toit"] == minimum]
            df_max = df[df["alpha_toit"] == maximum]
            df = df[df["alpha_toit"].isin([minimum, maximum])]
            df.reset_index(drop=True, inplace=True)

            for i, cpe in enumerate(["CPE 10", "CPE 1"]):
                row = [round(alpha, 2), cpe]
                for j in range(2, num_columns):
                    if touch_min:
                        row.append(df_min.iloc[i, j])
                    elif touch_max:
                        row.append(df_max.iloc[i, j])
                    else:
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                alpha, minimum, maximum, df_min.iloc[i, j], df_max.iloc[i, j]
                            )
                        )
                df.loc[df.shape[0]] = row

        if self.load_area > 1 and self.load_area < 10:
            for i in range(df.shape[0]):
                if i % 2 == 0:
                    row = [df.iloc[i, 0], "CPE " + str(round(self.load_area, 2))]
                    for j in range(2, num_columns):
                        row.append(
                            MathUtilsMixin.interpolation_logarithmique(
                                self.load_area, 1, 10, df.iloc[i + 1, j], df.iloc[i, j]
                            )
                        )
                    df.loc[df.shape[0]] = row
            df.reset_index(drop=True, inplace=True)

        df = df[df["alpha_toit"] == round(alpha, 2)]
        df.dropna(axis=1, how="all", inplace=True)
        df.set_index("alpha_toit", inplace=True)
        return df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self, dir=("0°", "90°")):
        """Retourne les dimensions des zones de pression pour la direction de vent donnée (fig. 7.9).

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            dict: Zones avec Longueur, Largeur et Surface en m / m².
        """
        return self._wind_direction[dir]["geometrie"]

    def get_Cpe(self, dir=("0°", "90°")):
        """Retourne les C_pe pour la direction de vent donnée (tab. 7.5).

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            pandas.DataFrame: C_pe par zone (F à N), interpolés selon la pente et load_area.
        """
        return self._wind_direction[dir]["Cpe"]


class Toiture_multiple(Vent):
    def __init__(self, load_area: float, nb_travees: int, config: str = "shed", *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture à plusieurs versants (sheds) selon EN 1991-1-4 §7.2.7.

        Direction de vent 0° uniquement (vent normal aux travées).
        Facteurs de position par travée (fig. 7.10) appliqués aux C_pe de la travée au vent :
        - config="shed" (sheds monopentes, fig. 7.10 a/b) : travée 1 = C_pe monopente (tab. 7.3a),
          travée 2 = 0.8·C_pe, travées ≥3 = 0.6·C_pe ; si le cas en pression (C_pe > 0) est retenu,
          les travées suivantes prennent C_pe = −0.4.
        - config="duopitch" (versants symétriques multiples en sheds, fig. 7.10 c/d) :
          travée 1 = C_pe deux versants (tab. 7.4a), travées suivantes = 0.6·C_pe.

        Zones (§7.2.7(2)) : F, G, J uniquement sur la face au vent (travée 1), H et I sur chaque travée.

        Args:
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
            nb_travees (int): Nombre de travées (versants) de la toiture (≥ 2).
            config (str): "shed" pour sheds monopentes, "duopitch" pour versants symétriques multiples.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent (d_bat = profondeur totale).
        """
        super().__init__(*args, **kwargs)
        if nb_travees < 2:
            raise ValueError(f"Une toiture multiple exige au moins 2 travées, reçu : {nb_travees}")
        self.load_area = load_area
        self.nb_travees = nb_travees
        self.config = self._check_choice(config, ("shed", "duopitch"), "config")
        self.d_travee = self.d_bat / nb_travees
        self.e = min(self.b_bat, self.h_bat * 2)
        self._wind_direction = {"0°": {}}
        self._wind_direction["0°"]["geometrie"] = self._geo()
        self._wind_direction["0°"]["Cpe"] = self._Cpe()

    def _geo(self):
        """Zones par travée : F, G, J sur la face au vent de la travée 1 (bande e/10 en rive
        et en noue), H et I (demi-travées) sur chaque travée. Surfaces / cos(alpha_toit)."""
        e = self.e
        b = self.b_bat
        cos_a = mt.cos(mt.radians(self.alpha_toit))
        geometrie = {}
        for i in range(1, self.nb_travees + 1):
            zones = {}
            if i == 1:
                # Face au vent : F, G en rive au vent (e/10), J en noue (e/10), H au milieu
                zones["F"] = {"Longueur": e / 4, "Largeur": e / 10, "Surface": (e / 4) * (e / 10) / cos_a}
                zones["G"] = {
                    "Longueur": b - e / 2,
                    "Largeur": e / 10,
                    "Surface": (b - e / 2) * (e / 10) / cos_a,
                }
                zones["J"] = {"Longueur": b, "Largeur": e / 10, "Surface": b * (e / 10) / cos_a}
                # Travées très étroites : les bandes F/G/J couvrent toute la travée -> H et I à 0
                partie_centrale = max(self.d_travee - 2 * e / 10, 0 * si.m) / 2
                zones["H"] = {
                    "Longueur": b,
                    "Largeur": partie_centrale,
                    "Surface": b * partie_centrale / cos_a,
                }
                zones["I"] = {
                    "Longueur": b,
                    "Largeur": partie_centrale,
                    "Surface": b * partie_centrale / cos_a,
                }
            else:
                # Travées suivantes : H sur la moitié au vent, I sur la moitié sous le vent
                zones["H"] = {
                    "Longueur": b,
                    "Largeur": self.d_travee / 2,
                    "Surface": b * (self.d_travee / 2) / cos_a,
                }
                zones["I"] = {
                    "Longueur": b,
                    "Largeur": self.d_travee / 2,
                    "Surface": b * (self.d_travee / 2) / cos_a,
                }
            geometrie[f"Travée {i}"] = zones
        return geometrie

    def _cpe_base(self):
        """C_pe de la travée au vent : tab. 7.3a (shed) ou tab. 7.4a (duopitch), interpolés."""
        if self.config == "shed":
            df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_1_versant_0_degres.csv"))
        else:
            df = self._data_from_csv(os.path.join("vent", "vent_Cpe_toiture_2_versants_0_degres.csv"))
        df.reset_index(drop=False, inplace=True)
        list_alpha_toit = df["alpha_toit"].unique()
        touch_min = self.alpha_toit < list_alpha_toit.min()
        touch_max = self.alpha_toit > list_alpha_toit.max()
        num_columns = df.shape[1]

        if self.alpha_toit not in list_alpha_toit:
            minimum = list_alpha_toit.min() if touch_min else list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
            maximum = list_alpha_toit.max() if touch_max else list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
            df_min = df[df["alpha_toit"] == minimum]
            df_max = df[df["alpha_toit"] == maximum]
            df = df[df["alpha_toit"].isin([minimum, maximum])]
            df.reset_index(drop=True, inplace=True)
            for i, cpe in enumerate(["CPE 10", "CPE 1"]):
                row = [round(self.alpha_toit, 2), cpe]
                for j in range(2, num_columns):
                    if touch_min:
                        row.append(df_min.iloc[i, j])
                    elif touch_max:
                        row.append(df_max.iloc[i, j])
                    else:
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                self.alpha_toit, minimum, maximum, df_min.iloc[i, j], df_max.iloc[i, j]
                            )
                        )
                df.loc[df.shape[0]] = row

        if self.load_area > 1 and self.load_area < 10:
            for i in range(df.shape[0]):
                if i % 2 == 0:
                    row = [df.iloc[i, 0], "CPE " + str(round(self.load_area, 2))]
                    for j in range(2, num_columns):
                        row.append(
                            MathUtilsMixin.interpolation_logarithmique(
                                self.load_area, 1, 10, df.iloc[i + 1, j], df.iloc[i, j]
                            )
                        )
                    df.loc[df.shape[0]] = row
            df.reset_index(drop=True, inplace=True)

        df = df[df["alpha_toit"] == round(self.alpha_toit, 2)]
        df.dropna(axis=1, how="all", inplace=True)
        df.set_index("alpha_toit", inplace=True)
        return df

    def _Cpe(self):
        """C_pe par travée avec les facteurs de position de la fig. 7.10.

        Returns:
            dict: {"Travée 1": DataFrame, "Travée 2": DataFrame, ...} — pour "shed", les zones
            J et I reprennent le C_pe de la zone H (pas de valeurs dédiées au tab. 7.3a).
            "shed" : travée 2 → 0.8·C_pe, ≥3 → 0.6·C_pe, ou −0.4 dans le cas en pression.
            "duopitch" : travées ≥2 → 0.6·C_pe.
        """
        df = self._cpe_base().copy()
        if self.config == "shed":
            # Zones J et I déduites du C_pe de la zone H (pas de valeurs dédiées au tab. 7.3a)
            for sig in ("-", "+"):
                if f"H{sig}" in df.columns:
                    df[f"J{sig}"] = df[f"H{sig}"]
                    df[f"I{sig}"] = df[f"H{sig}"]

        result = {}
        zone_cols = [c for c in df.columns if c != "Cpe"]
        pos_cols = [c for c in zone_cols if c.endswith("+")]
        for i in range(1, self.nb_travees + 1):
            df_i = df.copy()
            if self.config == "shed":
                # Pression sur la travée au vent : cas "+" -> -0.4 partout ensuite (fig. 7.10 b)
                df_i[pos_cols] = -0.4 if i >= 2 else df_i[pos_cols]
                if i == 2:
                    df_i[[c for c in zone_cols if c.endswith("-")]] *= 0.8
                elif i >= 3:
                    df_i[[c for c in zone_cols if c.endswith("-")]] *= 0.6
            elif i >= 2:  # duopitch : 0.6·C_pe sur toutes les zones des travées suivantes
                df_i[zone_cols] *= 0.6
            df_i.insert(0, "Travée", i)
            df_i.set_index("Travée", append=True, inplace=True)
            result[f"Travée {i}"] = df_i
        return result

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self):
        """Retourne les zones de pression par travée pour la direction 0°.

        Returns:
            dict: {"Travée 1": {"F": ..., "G": ..., "J": ..., "H": ..., "I": ...}, "Travée 2": ...}.
        """
        return self._wind_direction["0°"]["geometrie"]

    def get_Cpe(self):
        """Retourne les C_pe par travée pour la direction 0° avec facteurs de position.

        Returns:
            dict: {"Travée 1": pandas.DataFrame, ...}.
        """
        return self._wind_direction["0°"]["Cpe"]


class Toiture_isolee_1_pant(Vent):
    def __init__(self, phi: float, load_area: float, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture isolée à un versant selon EN 1991-1-4 §7.3.

        Le calcul est effectué pour les directions de vent 0° et 90°.
        Les zones de pression nettes sont A, B, C, déterminées à partir des dimensions du bâtiment.
        Les C_p nets (pression et dépression simultanées) sont interpolés selon phi et alpha_toit.
        Utilise les coefficients C_p (pression nette) et non C_pe/C_pi séparément.

        Args:
            phi (float): Degré d'obstruction sous la toiture (0 = libre, 1 = totalement obstruée).
                Valeur comprise entre 0 et 1. Interpolation linéaire pour les valeurs intermédiaires.
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        self.phi = phi
        self.load_area = load_area
        self._wind_direction = {"0°": {}, "90°": {}}

        for cle in self._wind_direction.keys():
            match cle:
                case "90°":
                    self.b_bat, self.d_bat = self.d_bat, self.b_bat

            self._wind_direction[cle]["geometrie"] = self._geo()
            self._wind_direction[cle]["Cp"] = self._Cp()

            if cle == "90°":
                self.b_bat, self.d_bat = self.d_bat, self.b_bat

    def _geo(self):
        """Calcul les surfaces du zonage de la toiture"""
        self._df = self._data_from_csv(os.path.join("vent", "vent_Cp_toiture_isolee_1_versant.csv"))
        self._df.reset_index(drop=False, inplace=True)

        geometrie = {
            "A": {
                "Longueur": self.b_bat - self.b_bat / 10 * 2,
                "Largeur": self.d_bat - self.d_bat / 10 * 2,
                "Surface": (self.d_bat - self.d_bat / 10 * 2)
                * (self.b_bat - self.b_bat / 10 * 2)
                / mt.cos(mt.radians(self.alpha_toit)),
            },
            "B": {
                "Longueur": self.d_bat - self.d_bat / 10 * 2,
                "Largeur": self.b_bat / 10,
                "Surface": ((self.d_bat - self.d_bat / 10 * 2) * self.b_bat / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
            "C": {
                "Longueur": self.b_bat - self.b_bat / 10 * 2,
                "Largeur": self.d_bat / 10,
                "Surface": ((self.b_bat - self.b_bat / 10 * 2) * self.d_bat / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
        }
        return geometrie

    def _Cp(self):
        """Calcul les Cp d'une toiture isolée

        Args:
                cle (str): sens du vent sur le bâtiment 0° ou 90°
        """

        if self.phi > 0 and self.phi < 1:
            for i in range(0, self._df.shape[0], 3):
                row = [self._df.iloc[i, 0], self.phi]
                for j in range(2, 6):
                    row.append(
                        MathUtilsMixin.interpolation_lineaire(
                            self.phi,
                            0,
                            1,
                            self._df.iloc[i + 1, j],
                            self._df.iloc[i + 2, j],
                        )
                    )
                self._df.loc[self._df.shape[0]] = row
            self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

        elif self.phi >= 1:
            self._df = self._df[self._df["phi"].isin(["1", "max"])]
        else:
            self._df = self._df[self._df["phi"].isin(["0", "max"])]

        list_alpha_toit = self._df["alpha_toit"].unique()
        touch_min = self.alpha_toit < list_alpha_toit.min()
        touch_max = self.alpha_toit > list_alpha_toit.max()

        if self.alpha_toit not in list_alpha_toit:
            if touch_min:
                minimum = list_alpha_toit.min()
            else:
                minimum = list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
            if touch_max:
                maximum = list_alpha_toit.max()
            else:
                maximum = list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
            df_min = self._df[self._df["alpha_toit"] == minimum]
            df_max = self._df[self._df["alpha_toit"] == maximum]

            self._df.reset_index(drop=True, inplace=True)
            for i, phi in enumerate(["max", self.phi]):
                row = [round(self.alpha_toit, 2), phi]
                for j in range(2, 6):
                    if touch_min:
                        row.append(df_min.iloc[i, j])
                    elif touch_max:
                        row.append(df_max.iloc[i, j])
                    else:
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                self.alpha_toit,
                                minimum,
                                maximum,
                                df_min.iloc[i, j],
                                df_max.iloc[i, j],
                            )
                        )
                self._df.loc[self._df.shape[0]] = row
            self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

        self._df = self._df[self._df["alpha_toit"] == round(self.alpha_toit, 2)]
        self._df.set_index("alpha_toit", inplace=True)
        return self._df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self, dir=("0°", "90°")):
        """Retourne les dimensions des zones de pression pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            dict: {"A": ..., "B": ..., "C": ...} avec Longueur, Largeur et Surface en m / m².
        """
        return self._wind_direction[dir]["geometrie"]

    def get_Cp(self, dir=("0°", "90°")):
        """Retourne les coefficients de pression nette C_p pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            pandas.DataFrame: C_p nets par zone (A, B, C), interpolés selon phi et alpha_toit.
                Inclut les valeurs maximales ("max") et pour la valeur de phi donnée.
        """
        return self._wind_direction[dir]["Cp"]

    def show_zonage(self):
        """Affiche l'image du zonage pour une toiture isolée à un versant selon EN 1991-1-4 §7.3."""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cp_toiture_isolee_1_versant.png")
        image = Image.open(file)
        image.show()


class Toiture_isolee_2_pants(Vent):
    def __init__(self, phi: float, load_area: float, *args, **kwargs):
        """Crée un objet pour le calcul des pressions sur une toiture isolée à deux versants selon EN 1991-1-4 §7.3.

        Le calcul est effectué pour les directions de vent 0° et 90°.
        Les zones de pression nettes sont A, B, C, D, déterminées à partir des dimensions du bâtiment.
        Les C_p nets sont interpolés selon phi et alpha_toit.
        Utilise les coefficients C_p (pression nette) et non C_pe/C_pi séparément.

        Args:
            phi (float): Degré d'obstruction sous la toiture (0 = libre, 1 = totalement obstruée).
                Valeur comprise entre 0 et 1. Interpolation linéaire pour les valeurs intermédiaires.
            load_area (float): Aire chargée en m² pour le calcul des éléments ou des fixations.
            *args: Arguments positionnels transmis à la classe parent Vent.
            **kwargs: Arguments nommés transmis à la classe parent Vent.
        """
        super().__init__(*args, **kwargs)
        self.phi = phi
        self.load_area = load_area
        self._wind_direction = {"0°": {}, "90°": {}}

        for cle in self._wind_direction.keys():
            match cle:
                case "90°":
                    self.b_bat, self.d_bat = self.d_bat, self.b_bat

            self._wind_direction[cle]["geometrie"] = self._geo()
            self._wind_direction[cle]["Cp"] = self._Cp()

            if cle == "90°":
                self.b_bat, self.d_bat = self.d_bat, self.b_bat

    def _geo(self):
        """Calcul les surfaces du zonage de la toiture"""
        self._df = self._data_from_csv(os.path.join("vent", "vent_Cp_toiture_isolee_2_versants.csv"))
        self._df.reset_index(drop=False, inplace=True)

        geometrie = {
            "A": {
                "Longueur": self.d_bat - self.d_bat / 10 * 2,
                "Largeur": self.d_bat - (2 * self.d_bat / 10) - self.d_bat / 5,
                "Surface": ((self.d_bat - self.d_bat / 10 * 2) * (self.d_bat - (2 * self.d_bat / 10) - self.d_bat / 5))
                / mt.cos(mt.radians(self.alpha_toit)),
            },
            "B": {
                "Longueur": self.d_bat - self.d_bat / 10 * 2,
                "Largeur": self.b_bat / 10,
                "Surface": ((self.d_bat - self.d_bat / 10 * 2) * self.b_bat / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
            "C": {
                "Longueur": self.b_bat - self.b_bat / 10 * 2,
                "Largeur": self.d_bat / 10,
                "Surface": ((self.b_bat - self.b_bat / 10 * 2) * self.d_bat / 10) / mt.cos(mt.radians(self.alpha_toit)),
            },
            "D": {
                "Longueur": (self.b_bat - self.b_bat / 10 * 2) / 2,
                "Largeur": (self.d_bat / 5) / 2,
                "Surface": (
                    ((self.b_bat - self.b_bat / 10 * 2) * self.d_bat / 5) / mt.cos(mt.radians(self.alpha_toit)) / 2
                ),
            },
        }
        return geometrie

    def _Cp(self):
        """Calcul les Cp d'une toiture isolée

        Args:
                cle (str): sens du vent sur le bâtiment 0° ou 90°
        """
        if self.phi > 0 and self.phi < 1:
            for i in range(0, self._df.shape[0], 3):
                row = [self._df.iloc[i, 0], self.phi]
                for j in range(2, 7):
                    row.append(
                        MathUtilsMixin.interpolation_lineaire(
                            self.phi,
                            0,
                            1,
                            self._df.iloc[i + 1, j],
                            self._df.iloc[i + 2, j],
                        )
                    )
                self._df.loc[self._df.shape[0]] = row
            self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

        elif self.phi >= 1:
            self._df = self._df[self._df["phi"].isin(["1", "max"])]
        else:
            self._df = self._df[self._df["phi"].isin(["0", "max"])]

        list_alpha_toit = self._df["alpha_toit"].unique()
        touch_min = self.alpha_toit < list_alpha_toit.min()
        touch_max = self.alpha_toit > list_alpha_toit.max()

        if self.alpha_toit not in list_alpha_toit:
            if touch_min:
                minimum = list_alpha_toit.min()
            else:
                minimum = list_alpha_toit[list_alpha_toit < self.alpha_toit].max()
            if touch_max:
                maximum = list_alpha_toit.max()
            else:
                maximum = list_alpha_toit[list_alpha_toit > self.alpha_toit].min()
            df_min = self._df[self._df["alpha_toit"] == minimum]
            df_max = self._df[self._df["alpha_toit"] == maximum]

            self._df.reset_index(drop=True, inplace=True)
            for i, phi in enumerate(["max", self.phi]):
                row = [round(self.alpha_toit, 2), phi]
                for j in range(2, 7):
                    if touch_min:
                        row.append(df_min.iloc[i, j])
                    elif touch_max:
                        row.append(df_max.iloc[i, j])
                    else:
                        row.append(
                            MathUtilsMixin.interpolation_lineaire(
                                self.alpha_toit,
                                minimum,
                                maximum,
                                df_min.iloc[i, j],
                                df_max.iloc[i, j],
                            )
                        )
                self._df.loc[self._df.shape[0]] = row
            self._df = self._df[self._df["phi"].isin([self.phi, "max"])]

        self._df = self._df[self._df["alpha_toit"] == round(self.alpha_toit, 2)]
        self._df.set_index("alpha_toit", inplace=True)
        return self._df

    def get_wind_dict(self) -> dict:
        return self._wind_direction

    def get_geo(self, dir=("0°", "90°")):
        """Retourne les dimensions des zones de pression pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            dict: {"A": ..., "B": ..., "C": ..., "D": ...} avec Longueur, Largeur et Surface en m / m².
        """
        return self._wind_direction[dir]["geometrie"]

    def get_Cp(self, dir=("0°", "90°")):
        """Retourne les coefficients de pression nette C_p pour la direction de vent donnée.

        Args:
            dir (str): Direction du vent ("0°" ou "90°"). Defaults to ("0°", "90°").

        Returns:
            pandas.DataFrame: C_p nets par zone (A, B, C, D), interpolés selon phi et alpha_toit.
                Inclut les valeurs maximales ("max") et pour la valeur de phi donnée.
        """
        return self._wind_direction[dir]["Cp"]

    def show_zonage(self):
        """Affiche l'image du zonage pour une toiture isolée à deux versants selon EN 1991-1-4 §7.3."""
        file = os.path.join(Vent.PATH_CATALOG, "data", "vent", "vent_Cp_toiture_isolee_2_versants.png")
        image = Image.open(file)
        image.show()


if __name__ == "__main__":
    building = Batiment(h_bat=5, d_bat=15, b_bat=13.1, alpha_toit=15, alt=400, code_INSEE=73215)
    Action_wind = Vent._from_parent_class(building, terrain="IIIa", oro="Aucun", z=5)
    print(si.environment())
    print(Action_wind.Vb[1])
    qpz = Action_wind.Qp_z
    print(qpz)
    ffr = Action_wind.Ffr(15, "Lisse")
    print(ffr)
    print(Action_wind.rayon_secteur_angu)
    # Action_wind.show_Ffr()
    vertical = Toiture_terrasse_acrotere._from_parent_class(Action_wind, load_area=10, hp=1, h=10.5)
    # vertical.show_zonage()
    print(vertical.get_Cpe())
