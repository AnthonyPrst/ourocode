#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
import math as mt
import pandas as pd
from PySide6.QtWidgets import QFileDialog
from matplotlib import pyplot as plt

import forallpeople as si
si.environment("structural")
from handcalcs.decorator import handcalc

# sys.path.append(os.path.join(os.getcwd(), "ourocode"))
# from eurocode.EC3_Element_droit import Plat

from ourocode.eurocode.EC3_Element_droit import Plat


class _TemperatureGaz(object):
    """
    Courbes de température des gaz θ_g(t).
    Par défaut : feu nominal ISO 834 / EN 1991-1-2 §3.2.1(1):
        θ_g(t) = 20 + 345 * log10(8 * t_min + 1)  (°C)
    """

    @staticmethod
    def iso_834(t_seconds: float) -> float:
        """
        Température des gaz selon la courbe normale ISO 834 (t en secondes).
        """
        t_minutes = t_seconds / 60
        return 20.0 + 345.0 * mt.log10(8.0 * t_minutes + 1.0)

    def _show_element(self, element):
        """Affiche un élément"""
        plt.plot(element)
        plt.show()

    @staticmethod
    def courbe_iso(times_min: float) -> list[float]:
        """Calcule θ_g pour une liste de temps (min) sur ISO 834 avec intervalle de 5 secondes"""
        times_secondes = [dt for dt in range(0, times_min*60+5, 5)]
        return [_TemperatureGaz.iso_834(t) for t in times_secondes]


class _CoeffFeu(object):
    """
    Calcul des coefficients de réduction au feu selon l'EN 1993-1-2
    pour l'acier carbone :
        - k_y,θ  : réduction de la limite d'élasticité (Tableau 3.1)
        - k_b,θ  : réduction de la résistance des boulons (Tableau 3.4)

    Les valeurs sont interpolées linéairement entre les températures
    normalisées du code (en °C).
    """

    # Tableau 3.1 EN 1993-1-2 (acier carbone)
    KY_THETA = {
        20: 1,
        100: 1,
        200: 1,
        300: 1,
        400: 1,
        500: 0.78,
        600: 0.47,
        700: 0.23,
        800: 0.11,
        900: 0.06,
        1000: 0.04,
        1100: 0.02,
        1200: 0.0,
    }

    KP_THETA = {
        20: 1,
        100: 1,
        200: 0.807,
        300: 0.613,
        400: 0.420,
        500: 0.360,
        600: 0.180,
        700: 0.075,
        800: 0.050,
        900: 0.0375,
        1000: 0.0250,
        1100: 0.0125,
        1200: 0.0,
    }

    KE_THETA = {
        20: 1,
        100: 1,
        200: 0.9,
        300: 0.8,
        400: 0.7,
        500: 0.6,
        600: 0.31,
        700: 0.13,
        800: 0.09,
        900: 0.0675,
        1000: 0.045,
        1100: 0.0225,
        1200: 0.0,
    }

    # Tableau D.1 EN 1993-1-2 (boulons classes 4.6 à 10.9)
    KB_THETA = {
        20: 1,
        100: 0.968,
        150: 0.952,
        200: 0.935,
        300: 0.903,
        400: 0.775,
        500: 0.550,
        600: 0.220,
        700: 0.100,
        800: 0.067,
        900: 0.033,
        1000: 0.0
    }

    KW_THETA = {
        20: 1,
        100: 1,
        150: 1,
        200: 1,
        300: 1,
        400: 0.876,
        500: 0.627,
        600: 0.378,
        700: 0.130,
        800: 0.074,
        900: 0.018,
        1000: 0.0
    }


    def __init__(self, theta: float, **kwargs):
        """
        Args:
            theta (float): température de calcul en °C.
        """
        super().__init__(**kwargs)
        self.theta = theta

    @staticmethod
    def _interp(table: dict[int, float], theta: float) -> float:
        """Interpolation linéaire entre les points du tableau EN 1993-1-2."""
        points = sorted(table.items())
        if theta <= points[0][0]:
            return points[0][1]
        if theta >= points[-1][0]:
            return points[-1][1]

        for (t1, v1), (t2, v2) in zip(points[:-1], points[1:]):
            if t1 <= theta <= t2:
                # interpolation linéaire
                ratio = (theta - t1) / (t2 - t1)
                return v1 + ratio * (v2 - v1)
        # fallback de sécurité (ne devrait pas se produire)
        return points[-1][1]

    @property
    def ky_theta(self) -> float:
        """Coefficient k_y,θ (acier carbone) à la température theta."""
        return self._interp(self.KY_THETA, self.theta)

    @property
    def kp_theta(self) -> float:
        """Coefficient k_p,θ (acier carbone) à la température theta."""
        return self._interp(self.KP_THETA, self.theta)

    @property
    def kE_theta(self) -> float:
        """Coefficient k_E,θ (acier carbone) à la température theta."""
        return self._interp(self.KE_THETA, self.theta)

    @property
    def kb_theta(self) -> float:
        """Coefficient k_b,θ (boulons) à la température theta."""
        return self._interp(self.KB_THETA, self.theta)

    @property
    def kw_theta(self) -> float:
        """Coefficient k_w,θ (boulons) à la température theta."""
        return self._interp(self.KW_THETA, self.theta)
    

    def synthese(self) -> dict:
        """Renvoie un récapitulatif pratique sous forme de dict."""
        return {
            "theta_C": self.theta,
            "k_y_theta": self.ky_theta,
            "k_p_theta": self.kp_theta,
            "k_E_theta": self.kE_theta,
            "k_b_theta": self.kb_theta,
            "k_w_theta": self.kw_theta,
        }


class Feu_acier(Plat):
    """
    Classe permettant le calcul de la température de l'acier θ_a(t) au feu selon EN 1993-1-2, §4.2 et retourne les facteurs de réduction k.
    Cette classe est hérité de la classe Plat du module EC3_Element_droit.py.
    Hypothèses : acier non protégé, flux convectif et radiatif simplifiés.

    Formule incrémentale (Euler) :
        θ_a(t+Δt) = θ_a(t) + (k_sh * (A/V) / (ρ_a * c_a)) * (h_c*(θ_g-θ_a) + ε_m*Φ*σ*((θ_r+273)^4-(θ_a+273)^4)) * Δt

    Par défaut : ε_m = 0.7, Φ = 1.0, ρ_a = 7850 kg/m³, c_a = 600 J/kgK.
    Les températures sont en °C, le temps en minutes, A/V en m⁻¹.
    """

    SIGMA = 5.67e-8  # Constante de Stefan-Boltzmann (W/m²K⁴)
    ALPHA_C = 25 # W/m²K coefficient de transfert thermique par convection selon EN 1991-1-2 §3.2.1(2)

    def __init__(
        self,
        time: int,
        am_v: float,
        theta_a0: float = 20.0,
        epsilon_m: float = 0.7,
        rho_a: float = 7850.0,
        k_sh: float = 1.0,
        **kwargs,
    ):
        """
        Args:
            time (float): durée attendue totale en minutes du feu.
            am_v (float): facteur de section A/V en m⁻¹ (surface exposée / volume acier).
            theta_a0 (float): température initiale de l'acier en °C.
            epsilon_m (float): émissivité de l'acier (0.7 typique).
            rho_a (float): masse volumique de l'acier kg/m³ (7850 par défaut).
            k_sh (float): facteur de correction forme/ombrage (k_sh = 1 par défaut).
        """
        super().__init__(**kwargs)
        self.time = time
        self.am_v = max(am_v, 10)
        self.theta_a = float(theta_a0)
        self.epsilon_m = epsilon_m
        self.rho_a = rho_a
        self.k_sh = k_sh
        self._calculate()

    def _hdot_net_d(self, theta_g: float, theta_a: float) -> float:
        """Flux thermique net (conv + rad) en W/m²."""
        # Convection
        hdot_net_c = self.ALPHA_C * (theta_g - theta_a) #Flux thermique net par convection en W/m²

        #Flux thermique net par radiation en W/m²

        # Dans le cas d'éléments complètement immergés dans le feu, la température de rayonnement Θr(t) peut 
        # être représentée par la température des gaz Θg(t) les entourant.  
        theta_r = theta_g
        epsilon_f = 1 # Emissivité du feu
        phi = 1 # facteur de forme
        hdot_net_r = (
            phi
            * self.epsilon_m
            * epsilon_f
            * self.SIGMA
            * ((theta_r + 273) ** 4 - (theta_a + 273) ** 4)
        )
        return hdot_net_c + hdot_net_r

    @staticmethod
    def _c_a(theta_a: float) -> float:
        """
        Chaleur spécifique de l'acier c_a,θ (J/kgK) selon EN 1993-1-2 §3.4.1.

        Args:
            theta_a (float): température de l'acier en °C.

        Returns:
            float: chaleur spécifique de l'acier en J/kgK.
        """
        if theta_a < 20:
            raise ValueError("La température initiale de l'acier θ_a doit être supérieure ou égale à 20°C.")
        if theta_a < 600:
            return (
                425
                + 7.73e-1 * theta_a
                - 1.69e-3 * theta_a**2
                + 2.22e-6 * theta_a**3
            )
        elif theta_a < 735:
            return 666 + 13002 / (738 - theta_a)
        elif theta_a < 900:
            return 545 + 17820 / (theta_a - 731)
        elif theta_a <= 1200:
            return 650
        else:
            raise ValueError("La température initiale de l'acier θ_a doit être inférieure ou égale à 1200°C.")
            
    def _step(self, theta_g: float, dt: float) -> float:
        """
        Avance d'un pas de temps et met à jour θ_a.

        Args:
            theta_g (float): température des gaz à l'instant courant en °C.
            dt (float): incrément de temps en s maximum de 5s.

        Returns:
            float: nouvelle température acier θ_a en °C.
        """
        if dt > 5:
            raise ValueError("dt doit être inférieur ou égal à 5s.")
        c_a = self._c_a(self.theta_a)
        hdot_net = self._hdot_net_d(theta_g, self.theta_a)
        dtheta = (
            self.k_sh
            * self.am_v
            / (self.rho_a * c_a)
            * hdot_net
            * dt
        )
        self.theta_a += dtheta
        return dtheta

    def _calculate(self) -> list[float]:
        """
        Calcule l'évolution de θ_a(t) sur une courbe de feu donnée.

        Args:
            time (float): durée totale en minutes (pas de 5 s).

        Returns:
            list[float]: températures acier θ_a pour chaque pas de 5 s.
        """
        times_secondes = [dt for dt in range(0, self.time*60+5, 5)]
        theta_g = _TemperatureGaz.courbe_iso(self.time)
        self.fire_data = pd.DataFrame(index=times_secondes, columns=["t (s)", "θg (°C)","Δθa,t (C°)", "θa (°C)", "ky,θ", "kp,θ", "kE,θ", "kb,θ", "kw,θ"])
        t_prev = 0
        d_theta = 0
        for index, t_curr in enumerate(times_secondes):
            dt = t_curr - t_prev
            t_prev = t_curr
            if index:
                d_theta = self._step(theta_g=theta_g[index], dt=dt)
            coeff_feu = _CoeffFeu(self.theta_a)
            self.fire_data.loc[t_curr, "t (s)"] = t_curr
            self.fire_data.loc[t_curr, "θg (°C)"] = theta_g[index]
            self.fire_data.loc[t_curr, "Δθa,t (C°)"] = d_theta
            self.fire_data.loc[t_curr, "θa (°C)"] = self.theta_a
            self.fire_data.loc[t_curr, "ky,θ"] = coeff_feu.ky_theta
            self.fire_data.loc[t_curr, "kp,θ"] = coeff_feu.kp_theta
            self.fire_data.loc[t_curr, "kE,θ"] = coeff_feu.kE_theta
            self.fire_data.loc[t_curr, "kb,θ"] = coeff_feu.kb_theta
            self.fire_data.loc[t_curr, "kw,θ"] = coeff_feu.kw_theta
        return self.fire_data

    def get_fire_data(self):
        """Récupère les données du feu normalisées et de l'acier avec un intervalle de 5 secondes.
        """
        return self.fire_data

    def get_specific_time_data(self, time: float):
        """Récupère les données du feu normalisées et de l'acier au temps donné.

        Args:
            time (float): temps en minutes.
        """
        if time < 0 or time > self.time:
            raise ValueError(f"Le temps doit être compris entre 0 et la durée totale du feu, soit {self.time} minutes.")
        return self.fire_data.loc[time*60]

    def show_temperatures(self, screenshot: bool = ("False", "True"), filepath: str=None):
        """Affiche l'évolution de la température θg(t) et θa(t) en fonction du temps (minutes) avec pas 5 s.
        
        Args:
            screenshot (bool): si True, enregistre le graphique
            filepath (str): chemin d'enregistrement du graphique, si ce dernier est vide, 
                alors une boite de dialogue s'ouvre pour choisir le chemin.
        """
        plt.figure(figsize=(10, 5))
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["θg (°C)"], label="Température des gaz θg(t) (°C)")
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["θa (°C)"], label="Température de l'acier θa(t) (°C)")
        plt.xlabel("Temps (min)")
        plt.ylabel("Température (°C)")
        plt.grid(True, which="both", linestyle="--", alpha=0.4)
        plt.legend()
        plt.tight_layout()
        if screenshot:
            if not filepath:
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
            plt.savefig(filepath)
            return filepath
        else:
            plt.show()

    def show_reductions_factors(self, screenshot: bool = ("False", "True"), filepath: str=None):
        """Affiche le graphiques des facteurs de réduction ky,θ, kb,θ et kE,θ en fonction du temps (minutes) avec pas 5 s.
        
        Args:
            screenshot (bool): si True, enregistre le graphique
            filepath (str): chemin d'enregistrement du graphique, si ce dernier est vide, 
                alors une boite de dialogue s'ouvre pour choisir le chemin.
        """
        plt.figure(figsize=(10, 5))
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["ky,θ"], label="Facteur de réduction ky,θ")
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["kp,θ"], label="Facteur de réduction kp,θ")
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["kE,θ"], label="Facteur de réduction kE,θ")
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["kb,θ"], label="Facteur de réduction pour les boulons kb,θ")
        plt.plot(self.fire_data["t (s)"]/60, self.fire_data["kw,θ"], label="Facteur de réduction pour les soudures kw,θ")
        plt.xlabel("Temps (min)")
        plt.ylabel("Facteur de réduction")
        plt.grid(True, which="both", linestyle="--", alpha=0.4)
        plt.legend()
        plt.tight_layout()
        if screenshot:
            if not filepath:
                filepath = QFileDialog.getSaveFileName(
                    filter="PNG (*.png)",
                    selectedFilter=".png",
                )[0]
            plt.savefig(filepath)
            return filepath
        else:
            plt.show()


if __name__ == "__main__":
    plat = Plat(10, 200, 200, "S235", 3)
    modele = Feu_acier._from_parent_class(plat, am_v=100, theta_a0=20, epsilon_m=0.8, time=60)
    modele.show_temperatures(screenshot=False)
    modele.show_reductions_factors(screenshot=False)
    print(modele.get_specific_time_data(60))