#! env\Scripts\python.exe
# Encoding in UTF-8 by Anthony PARISOT
from math import sqrt, radians, cos, sin, floor
from copy import deepcopy
import pandas as pd
import numpy as np
import forallpeople as si

si.environment("structural")
from handcalcs.decorator import handcalc

from ourocode.eurocode.A0_Projet import Batiment
from ourocode.eurocode.EC5_Assemblage import Assemblage

class MOB(Batiment):
    CRITERE_ELS = 200
    POSITION_PANEL = ('Extérieur', 'Intérieur')
    COEFF_DOUBLE_PANEL = {"Non": 1, "Panneau et organe identiques": 1, "Panneau différent & glissement identique": 0.75, "Panneau différent & glissement différent": 0.5}
    def __init__(
        self,
        connecteurs_int: object=None, 
        connecteurs_ext: object=None, 
        **kwargs):
        """Classe permettant de définir des systemes de mur à ossature bois (MOB) d'un projet selon l'EN 1995-1-1 §9.2.4 méthode A.
        Cette classe est hérité de la classe Batiment, provenant du module A0_Projet.py.



        Un système de mur à ossature bois est composé de plusieurs murs servant aux contreventements, 
        ces murs sont intérompus par les portes et fenêtres. Ces murs sont eux même composés de plusieurs panneaux.

        Il convient donc d'ajouter une classe MOB puis d'ajouter via la méthode "add_sys_wall" chaque système de mur. 
        Ensuite d'ajouter les murs internes au système de mur via la méthode "add_wall_to_sys", puis d'ajouter ensuite les panneaux
        à chaque murs via la méthode "add_panel_to_wall".
        Une méthode plus efficace consiste à écrire dans un fichier JSON les informations pour chaque système de mur, 
        puis d'utiliser la méthode "load_walls_data" pour charger les informations.

        Ensuite on peu utiliser la méthode "calculate_loads_and_deformations" pour calculer les efforts de contreventements et les déformations,
        suivant le centre de gravité et la raideur des systèmes de mur. 
        Sinon on peux donner directement les efforts Fv_Ed_ELS et Fv_Ed_ELU par système de mur dans la méthode "add_sys_wall".

        Enfin on utilise les méthodes "taux_walls" et "taux_ancrage" pour récupérer les taux de contreventements et les taux d'ancrage.
        Attention cette classe ne calcule pas les montants bois notamment au flambement, ni la compression perpendiculaire des montants sur la lisse basse.
        Pour réaliser ces calculs, il faut utiliser la classe Barre et les classes de vérification associées du module EC5_Element_droit.py

        Args:
            connecteurs_int (Agrafe, Pointe, optional): Objet Agrafe ou Pointe, représentant les connecteurs intérieurs du mur à ossature bois provenant du module EC5_Assemblage.py.
                Attention: lors de l'instanciation de la classe Assemblage puis de la tige elle-même, il faut un n_file = 1 et un n = 1.

            connecteurs_ext (Agrafe, Pointe, optional): Objet Agrafe ou Pointe, représentant les connecteurs extérieurs du mur à ossature bois provenant du module EC5_Assemblage.py.
                Attention: lors de l'instanciation de la classe Assemblage puis de la tige elle-même, il faut un n_file = 1 et un n = 1.
        """
        super().__init__(**kwargs)
        self._connecteurs = {"Intérieur": connecteurs_int, "Extérieur": connecteurs_ext}
        self._data_MOB = {}
        self._data_connecteurs = {}
        self._list_on_top_of = []
        self._setup_MOB()
    
    def _setup_MOB(self):
        """Initialise les matériaux, sections et valeurs de calcul pour les murs à ossature bois."""
        for position, connecteur in self._connecteurs.items():
            if connecteur is None:
                continue
            effet_corde = False
            if connecteur.type_organe in ["Tirefond", "Pointe carrée lisse", "Pointe circulaire lisse", "Autres pointes"]:
                connecteur.Fax_Rk()
                effet_corde = True
            Fvrk = connecteur.FvRk(effet_corde=effet_corde)
            Kser = connecteur.Kser()
            self._data_connecteurs[position] = {
                "Ff,Rk": Fvrk[1],
                "Kser": Kser[1]
            }

    def _get_all_sys_walls_names(self):
        """Récupère tous les noms des systèmes de mur."""
        data = []
        for etage, sys in self.data.items():
            for sys_name in sys.keys():
                data.append(sys_name)
        return data
    
    def calculate_loads_and_deformations(self, Fw_Ed_ELS: si.kN, Fw_Ed_ELU: si.kN, XFw: si.m=0, YFw: si.m=0, alpha_Fw: float=0, etage: str=Batiment.ETAGE,):
        """Calcule automatiquement les efforts de contreventements et les déformations pour les murs à ossature bois,
        suivant le centre de gravité et la raideur des systèmes de mur. 
        Cette méthode est issue du guide AQCEN MOB et ne fonctionne que si le diaphragme de plancher ou de toiture est
        considéré parfaitement rigide, donc avec de faible ouverture et trémie. 

        Args:
            Fw_Ed_ELS (float): Force résultante de contreventement à la combinaison ELS en kN au niveau de l'étage.
            
            Fw_Ed_ELU (float): Force résultante de contreventement à la combinaison ELU en kN au niveau de l'étage.
            
            XFw (float): Position X en m depuis l'origine du repère cartésien choisi par le calculateur, 
                de la force résultante de contreventement sur un système de mur.
            
            YFw (float): Position Y en m depuis l'origine du repère cartésien choisi par le calculateur, 
                de la force résultante de contreventement sur un système de mur.
            
            alpha_Fw (float): Angle de la direction d'application de la force résultante de contreventement 
                sur un système de mur en degrés par rapport au repère cartésien choisi par le calculateur inscrit dans le cercle trigonométrique.
            
            etage (str): nom de l'étage où s'applique la force résultante de contreventement.
        """
        Fw_Ed_ELS = Fw_Ed_ELS * si.kN
        Fw_Ed_ELU = Fw_Ed_ELU * si.kN
        XFw = XFw * si.m
        YFw = YFw * si.m

        alpha_Fw_rad = radians(alpha_Fw)
        Ks = self.get_Kser_wall()
        sys_walls_names = self._get_all_sys_walls_names()
        df_sw = pd.DataFrame(index=sys_walls_names, columns=["dn_p_sw", "dn_sw", "a", "b", "c", "d", "e", "C1", "C2", "C3", "C6"])
        if not hasattr(self, "data_sys_walls_loads"):
            self.data_sys_walls_loads = pd.DataFrame(index=sys_walls_names, columns=["Étage", "delta,v,Ed,i ELS", "Fv,Ed,i ELU"])
            self.data_walls_loads = pd.DataFrame(index=self.data_walls.keys(), columns=["Étage", "Fv,Ed,i ELU", "Ft,Ed,i ELU (soulèvement)"])
            self.data_panels_loads = pd.DataFrame(index=self.data_panels.keys(), columns=["Étage", "Fv,Ed,i ELU"])

        for etage_data, sys in self.data.items():
            if etage_data != etage:
                continue
            for sys_name, sys_wall in sys.items():
                xg = sys_wall.get("Xg", 0)
                yg = sys_wall.get("Yg", 0)
                alpha_sw = radians(sys_wall.get("alpha_sw", 0))
                gamma_sw = alpha_sw + radians(90)
                Ksw = Ks.loc[sys_name]['Kser système de mur']
                Kp_ref = Ks.loc[sys_name]['Kser,p,ref']
                Ceq_sw = Ks.loc[sys_name]['Coeff. équivalence']

                dn_p_sw = xg * cos(gamma_sw) + yg * sin(gamma_sw)
                df_sw.loc[sys_name]["dn_p_sw"] = dn_p_sw
                df_sw.loc[sys_name]["a"] = Ksw * cos(gamma_sw)**2
                df_sw.loc[sys_name]["b"] = Ksw * cos(gamma_sw) * sin(gamma_sw)
                df_sw.loc[sys_name]["c"] = Ksw * sin(gamma_sw)**2
                df_sw.loc[sys_name]["d"] = Ksw * cos(gamma_sw) * dn_p_sw
                df_sw.loc[sys_name]["e"] = Ksw * sin(gamma_sw) * dn_p_sw
                df_sw.loc[sys_name]["C1"] = Ceq_sw * cos(alpha_sw)**2
                df_sw.loc[sys_name]["C2"] = Ceq_sw * cos(alpha_sw) * sin(alpha_sw)
                df_sw.loc[sys_name]["C3"] = Ceq_sw * sin(alpha_sw)**2


            a = df_sw["a"].sum()
            b = df_sw["b"].sum()
            c = df_sw["c"].sum()
            d = df_sw["d"].sum()
            e = df_sw["e"].sum()
            xO = (c*d - b*e) / (a*c - b**2) # Position X du centre de raideur
            yO = (a*e - b*d) / (a*c - b**2) # Position Y du centre de raideur
            dw_x = XFw - xO
            dw_y = YFw - yO
            dw_n = dw_x * sin(alpha_Fw_rad) - dw_y * cos(alpha_Fw_rad)

            for sys_name, sys_wall in sys.items():
                Ceq_sw = Ks.loc[sys_name]['Coeff. équivalence']
                dn_p_sw = df_sw.loc[sys_name]["dn_p_sw"]
                dn_sw = xO * cos(gamma_sw) + yO * sin(gamma_sw) - dn_p_sw
                df_sw.loc[sys_name]["dn_sw"] = dn_sw
                df_sw.loc[sys_name]["C6"] = Ceq_sw * dn_sw**2

            C1 = df_sw["C1"].sum()
            C2 = df_sw["C2"].sum()
            C3 = df_sw["C3"].sum()
            C4 = Fw_Ed_ELS * cos(alpha_Fw_rad) / Kp_ref
            C5 = Fw_Ed_ELS * sin(alpha_Fw_rad) / Kp_ref
            C6 = Kp_ref * df_sw["C6"].sum()

            u = (C4*C3 - C5*C2) / (C1*C3 - C2**2)
            v = (C1*C5 + C4*C2) / (C1*C3 - C2**2)
            teta = (Fw_Ed_ELS * dw_n) / C6

            for sys_name, sys_wall in sys.items():
                alpha_sw = radians(sys_wall.get("alpha_sw", 0))
                Ksw = Ks.loc[sys_name]['Kser système de mur']
                self.data_sys_walls_loads.loc[sys_name]["Étage"] = etage_data

                delta_v_i_ed = teta * df_sw.loc[sys_name]["dn_sw"] + u * cos(alpha_sw) + v * sin(alpha_sw)
                self.data_sys_walls_loads.loc[sys_name]["delta,v,Ed,i ELS"] = delta_v_i_ed
                
                Fv_i_Ed_sw = delta_v_i_ed * Ksw * (Fw_Ed_ELU / Fw_Ed_ELS)
                self.data_sys_walls_loads.loc[sys_name]["Fv,Ed,i ELU"] = Fv_i_Ed_sw

                for wall_name, wall in sys_wall['Murs'].items():
                    Kw = self.data_Kser_walls.loc[wall_name]['Kser mur']
                    Fv_i_Ed_wall = Fv_i_Ed_sw * Kw/Ksw
                    self.data_walls_loads.loc[wall_name]["Étage"] = etage_data
                    self.data_walls_loads.loc[wall_name]["Fv,Ed,i ELU"] = Fv_i_Ed_wall
                    # traction dans le montant d'extrémité
                    self.data_walls_loads.loc[wall_name]["Ft,Ed,i ELU (soulèvement)"] = self._FtEd(wall_name)

                    for panel_name, panel in wall['Panneaux'].items():
                        if not panel['Contreventement']:
                            continue
                        Kp = self.data_Kser_panels.loc[panel_name]['Kser,p unit.']
                        self.data_panels_loads.loc[panel_name]["Étage"] = etage_data
                        self.data_panels_loads.loc[panel_name]["Fv,Ed,i ELU"] = Fv_i_Ed_wall * Kp/Kw

        return self.data_sys_walls_loads
        
    @property
    def data(self):
        return self._data_MOB

    @data.setter
    def data(self, data: dict):
        self._data_MOB = data

    @property
    def data_sys_walls(self):
        sys_wall_data = {}
        _data = deepcopy(self._data_MOB)
        for etage, sys in _data.items():
            for sys_name, sys_wall in sys.items():
                sys_wall.pop("Murs", None)
                sys_wall_data[sys_name] = sys_wall
        return sys_wall_data
    
    @property
    def data_walls(self):
        wall_data = {}
        _data = deepcopy(self._data_MOB)
        for etage, sys in _data.items():
            for sys_name, sys_wall in sys.items():
                for wall_name, wall in sys_wall['Murs'].items():
                    wall.pop("Panneaux", None)
                    wall_data[wall_name] = wall
        return wall_data
    
    @property
    def data_panels(self):
        panel_data = {}
        for etage, sys in self._data_MOB.items():
            for sys_name, sys_wall in sys.items():
                for wall_name, wall in sys_wall['Murs'].items():
                    for panel_name, panel in wall['Panneaux'].items():
                        panel_data[panel_name] = panel
        return panel_data

    def add_sys_wall(self,
        h_etage: si.m, 
        h_sys_MOB: si.m, 
        l_sys_MOB: si.m,
        etage: str=Batiment.ETAGE,
        Fv_Ed_ELS: si.kN=None,
        Fv_Ed_ELU: si.kN=None,
        Xg: si.m=0,
        Yg: si.m=0,
        alpha_sw: float=0,
        ):
        """
        Ajoute un système de mur à ossature bois.

        Args:
            h_etage (float): hauteur de l'étage en m.
            h_sys_MOB (float): hauteur du système de mur à ossature bois en m.
            l_sys_MOB (float): longueur du système de mur à ossature bois en m.
            etage (str): nom de l'étage où se trouve le système de mur bois.

            Fv_Ed_ELS (float, optional): Effort de contreventement à l'ELS en kN.
                Cette valeur n'est à rentrer que si l'on souhaite écraser la valeur calculée automatiquement.

            Fv_Ed_ELU (float, optional): Effort de contreventement à l'ELU en kN.
                Cette valeur n'est à rentrer que si l'on souhaite écraser la valeur calculée automatiquement.

            Les arguments suivants sont optionnels, s'ils sont rentrés, ils permettent de calculer les déformations et les efforts de contreventement:

            Xg (float, optional): position X depuis l'origine du repère cartésien choisi par le calculateur, 
                du centre de gravité du système de mur à ossature bois en m.

            Yg (float, optional): position Y depuis l'origine du repère cartésien choisi par le calculateur, 
                du centre de gravité du système de mur à ossature bois en m.

            alpha_sw (float, optional): angle de rotation suivant le cercle trigonométrique du système de mur à ossature bois en degrés à l'origine O,i,g via le repère cartésien choisi par le calculateur.
                Cela veux donc dire que l'angle de rotation des systèmes de mur est mesuré dans le sens inverse des aiguilles d'une montre depuis l'origine gravitaire du système de mur.
        
        Returns:
            str: Nom du système de mur ajouté.
        """
        # Nommer les systèmes de mur de façon continue indépendamment des étages
        total_sys = sum(len(sys) for sys in self.data.values()) if self.data else 0
        name = f"S{total_sys + 1}"
        
        sys_wall = {
            'Etage': etage,
            'Hauteur étage': h_etage * si.m,
            'Hauteur du système de mur': h_sys_MOB * si.m,
            'Longueur du système de mur': l_sys_MOB * si.m,
            'Murs': {}
        }
        if Xg != 0 or Yg != 0 or alpha_sw != 0:
            sys_wall['Xg'] = Xg * si.m
            sys_wall['Yg'] = Yg * si.m
            sys_wall['alpha_sw'] = alpha_sw
        
        if Fv_Ed_ELU is not None:
            sys_wall['Fv_Ed_ELU'] = Fv_Ed_ELU * si.kN
        if Fv_Ed_ELS is not None:
            sys_wall['Fv_Ed_ELS'] = Fv_Ed_ELS * si.kN

        if self.data.get(etage) is None:
            self.data[etage] = {}
        self.data[etage][name] = sys_wall
        return name


    def add_wall_to_sys(self,
        sys_name: str,
        position: si.m,
        l_MOB: si.m,
        couturage_ext_rive: si.mm=150,
        couturage_ext_inter: si.mm=300,
        couturage_int_rive: si.mm=0,
        couturage_int_inter: si.mm=0,
        Ff_Rk_ext: si.N=None,
        Ff_Rk_int: si.N=None,
        Kser_ext: si.N/si.mm=None,
        Kser_int: si.N/si.mm=None,
        ):
        """
        Ajoute un mur interne au système de mur à ossature bois.

        Args:
            sys_name (str): Nom du système de mur auquel le mur est ajouté.

            position (float): Position X à partir du coin inférieur gauche du mur dans le système de mur en m.
                On regarde le système de mur comme un repère cartésien avec le coin inférieur gauche en (0,0), 
                avec la face externe du système de mur comme définition du sens des abscisses.

            l_MOB (float): Longueur du mur en m.

            couturage_ext_rive (float): Entraxe du couturage rive face extérieure en mm.
            couturage_ext_inter (float): Entraxe du couturage intermédiaire face extérieure en mm.
            couturage_int_rive (float): Entraxe du couturage rive face intérieure en mm.
            couturage_int_inter (float): Entraxe du couturage intermédiaire face intérieure en mm.:

            Ff_Rk_ext (float, optional): Capacité résistante caractéristique de cisaillement du connecteur en N sur la face extérieure.
                Attention, aucune majoration ne doit être faite sur cette valeur. Cette valeur n'est à rentrer que si l'on souhaite écraser la valeur calculée.
                
            Ff_Rk_int (float, optional): Capacité résistante caractéristique de cisaillement du connecteur en N sur la face intérieure.
                Attention, aucune majoration ne doit être faite sur cette valeur. Cette valeur n'est à rentrer que si l'on souhaite écraser la valeur calculée.

            Kser_ext (float, optional): Raideur du connecteur entre le panneau extérieur et les montants bois en N/mm.
                Cette valeur n'est à rentrer que si l'on souhaite écraser la valeur calculée.

            Kser_int (float, optional): Raideur du connecteur entre le panneau intérieur et les montants bois en N/mm.
                Cette valeur n'est à rentrer que si l'on souhaite écraser la valeur calculée.
        
        Returns:
            str: Nom du mur ajouté.
        """
        etage_found = None
        for etage, sys in self.data.items():
            if sys_name in sys.keys():
                etage_found = etage
                break
        if etage_found is None:
            raise ValueError(f"Le système de mur {sys_name} n'existe pas")

        name = f"{sys_name}_W{len(self.data[etage_found][sys_name]['Murs']) + 1}"
        wall = {
            'Position en X': position * si.m,
            'Longueur': l_MOB * si.m,
            'Panneaux': {}
        }

        if couturage_ext_rive:
            wall['Couturage ext. rive'] = couturage_ext_rive * si.mm
            wall['Couturage ext. intermédiaire'] = couturage_ext_inter * si.mm
            if Ff_Rk_ext is None:
                wall['Ff,Rk ext'] = self._data_connecteurs['Extérieur']['Ff,Rk']
            else:
                wall['Ff,Rk ext'] = Ff_Rk_ext * si.N
            if Kser_ext is None:
                wall['Kser ext'] = self._data_connecteurs['Extérieur']['Kser']
            else:
                wall['Kser ext'] = Kser_ext * si.N / si.mm

        if couturage_int_rive:
            wall['Couturage int. rive'] = couturage_int_rive * si.mm
            wall['Couturage int. intermédiaire'] = couturage_int_inter * si.mm
            if Ff_Rk_int is None:
                wall['Ff,Rk int'] = self._data_connecteurs['Intérieur']['Ff,Rk']
            else:
                wall['Ff,Rk int'] = Ff_Rk_int * si.N
            if Kser_int is None:
                wall['Kser int'] = self._data_connecteurs['Intérieur']['Kser']
            else:
                wall['Kser int'] = Kser_int * si.N / si.mm

        self.data[etage_found][sys_name]['Murs'][name] = wall
        return name
    
    def add_panel_to_wall(self,
        wall_name: str,
        h_panel: si.mm,
        b_panel: si.mm,
        number: int=1,
        position_panel: str=POSITION_PANEL,
        on_top_of: str=None,
        ):
        """
        Ajoute un panneau à un mur.

        Args:
            wall_name (str): Nom du mur auquel le panneau est ajouté.
            h_panel (float): Hauteur du panneau en mm (dimension la plus grande en théorie).
            b_panel (float): Largeur du panneau en mm (dimension la plus petite en théorie).
            number (int): Nombre de panneaux de même dimension.
            position_panel (str): Position du panneau par rapport au mur.
            on_top_of (str): Si le panneau est placé au-dessus d'un panneau, donner le nom du panneau au-dessus duquel le panneau est ajouté.
        
        Returns:
            str: Nom du panneau ajouté.
        """
        sys_name = wall_name.split("_")[0]
        etage_found = None
        for etage, sys in self.data.items():
            if sys_name in sys.keys():
                etage_found = etage
                break
        if etage_found is None:
            raise ValueError(f"Le système de mur {sys_name} n'existe pas")

        name = f"{wall_name}_P{len(self.data[etage_found][sys_name]['Murs'][wall_name]['Panneaux']) + 1}"
        wall_data = self.data[etage_found][sys_name]['Murs'][wall_name]
        
        panel = {
            'Type': self._connecteurs[position_panel].beam_1.classe,
            'Position': position_panel,
            'Panneau double': "Non",
            'Contreventement': True,
            'Nombre': number,
            'Epaisseur': self._connecteurs[position_panel].beam_1.b,
            'Hauteur': h_panel * si.mm,
            'Largeur': b_panel * si.mm,
        }

        # On test si le panneau est contreventant
        if h_panel > b_panel:
            if h_panel / 4 > b_panel:
                panel['Contreventement'] = False
        else:
            if b_panel / 4 > h_panel:
                panel['Contreventement'] = False
        
        if on_top_of:
            panel['On top of'] = on_top_of
            self._list_on_top_of.append(name)

        # Si le mur est composé potentiellement de deux panneaux, vérifier que les 2 côtés existent réellement
        if (
            self._connecteurs['Intérieur'] and self._connecteurs['Extérieur']
            and ('Kser ext' in wall_data) and ('Kser int' in wall_data)
        ):
            # Si les Kser sont différents
            if wall_data['Kser ext'] != wall_data['Kser int']:
                panel['Panneau double'] = "Panneau différent & glissement différent"
            # Si les classes sont identiques
            elif self._connecteurs['Intérieur'].beam_1.classe == self._connecteurs['Extérieur'].beam_1.classe:
                panel['Panneau double'] = "Panneau différent & glissement identique"
                # Si les largeurs sont identiques
                if self._connecteurs['Intérieur'].beam_1.b == self._connecteurs['Extérieur'].beam_1.b:
                    if ('Ff,Rk ext' in wall_data) and ('Ff,Rk int' in wall_data) and (wall_data['Ff,Rk ext'] != wall_data['Ff,Rk int']):
                        panel['Panneau double'] = "Panneau et organe identiques"
            else:
                panel['Panneau double'] = "Panneau différent & glissement identique"
        self.data[etage_found][sys_name]['Murs'][wall_name]['Panneaux'][name] = panel
        return name

    def save_walls_data(self, path: str=None):
        """Sauvegarde les données des murs dans un fichier JSON.

        Args:
            path (str, optional): Chemin du fichier à créer, s'il n'est pas fourni, une boite de dialogue s'ouvre 
                pour choisir le fichier. Defaults to None.
        """
        super().save_data(self.data, type_data="JSON", path=path)
    
    def load_walls_data(self, path: str=None):
        """Charge les données des murs depuis un fichier JSON.

        Args:
            path (str, optional): Chemin du fichier à charger, s'il n'est pas fourni, une boite de dialogue s'ouvre 
                pour choisir le fichier. Defaults to None.
        """
        self.data = super().load_data(type_data="JSON", path=path)
        return self.data

    def _get_Kser_panel(self, b: float, h: float, k_ser_r: float, k_ser_i: float):
        """
        Calcul la raideur du panneau.

        Args:
            b (float): largeur du panneau.
            h (float): hauteur du panneau.
            k_ser_r (float): raideur du panneau.
            k_ser_i (float): raideur du panneau.
        """ 
        K_p = (2 / (k_ser_r * (b + h/3) + k_ser_i * h/6) + 2*h**2 / (k_ser_r * b**2 * (h + b/3)))**-1
        return K_p

    def _get_Kser_panel_stack(self, top_panel_name: str, bottom_panel_name: str):
        """
        Calcul la raideur du panneau en pile.

        Args:
            top_panel_name (str): nom du panneau supérieur.
            bottom_panel_name (str): nom du panneau inférieur.
        """
        k2 = self.data_Kser_panels.loc[top_panel_name]['Kser,p unit.']
        k1 = self.data_Kser_panels.loc[bottom_panel_name]['Kser,p unit.']
        k_equi = (k1 * k2) / (k1 + k2)
        return k_equi

    def _k_panel(self):
        """
        Calcul la raideur des panneaux.
        """
        df_kp = pd.DataFrame(columns=['Étage', 'Type', 'Épaisseur', 'Largeur', 'Hauteur', 'Nombre', 'Kser,p unit.', 'Kser,p,ref', 'Coeff. équivalence'])
        df_kp.index.name = 'Panneau'
        K_p_ref = {}
        for etage, sys in self.data.items():
            for sys_name, sys_wall in sys.items():
                h_sys = sys_wall['Hauteur du système de mur']
                for wall_name, wall in sys_wall['Murs'].items():
                    if wall.get('Kser ext', None):
                        kser_ext = wall['Kser ext']
                        couturage_ext_rive = wall['Couturage ext. rive']
                        couturage_ext_inter = wall['Couturage ext. intermédiaire']
                    if wall.get('Kser int', None):
                        kser_int = wall['Kser int']
                        couturage_int_rive = wall['Couturage int. rive']
                        couturage_int_inter = wall['Couturage int. intermédiaire']
                            
                    for panel_name, panel in wall['Panneaux'].items():
                        # On limite la largeur efficace de contreventement des panneaux
                        if not panel['Contreventement']:
                            continue

                        if panel['Position'] == 'Extérieur':
                            k_ser_r = kser_ext / couturage_ext_rive
                            k_ser_i = kser_ext / couturage_ext_inter
                        else:
                            k_ser_r = kser_int / couturage_int_rive
                            k_ser_i = kser_int / couturage_int_inter
                        
                        h = panel['Hauteur']
                        b = panel['Largeur']
                        K_p = self._get_Kser_panel(b, h, k_ser_r, k_ser_i)
                        
                        df_kp.loc[panel_name] = [etage, panel['Type'], panel['Epaisseur'], b, h, panel['Nombre'], K_p, None, None]
                        if K_p > K_p_ref.get(etage, 0*si.N / si.mm):
                            K_p_ref[etage] = K_p
        
        # Trier par index (nommé 'Panneau') au lieu d'une colonne inexistante
        df_kp.sort_index(inplace=True)
        for panel_name, panel in df_kp.iterrows():
            etage = panel['Étage']
            df_kp.loc[panel_name,'Kser,p,ref'] = K_p_ref.get(etage, 0*si.N / si.mm)
            C_eq =  panel['Kser,p unit.'] / K_p_ref.get(etage, 0*si.N / si.mm)
            df_kp.loc[panel_name,'Coeff. équivalence'] = C_eq

        return df_kp

    def get_Kser_wall(self):
        """
        Calcul la raideur des murs du système de mur sélectionné.

        Args:
            sys_name (str): Nom du système de mur auquel le mur est ajouté.
        """
        self.data_Kser_panels = self._k_panel()
        # Colonnes alignées avec les valeurs calculées: [Étage, Kser, Coeff, Kser,p,ref]
        self.data_Kser_sys_walls = pd.DataFrame(columns=['Étage', 'Kser système de mur', 'Coeff. équivalence', 'Kser,p,ref'])
        self.data_Kser_sys_walls.index.name = 'Système de mur'
        self.data_Kser_walls = pd.DataFrame(columns=['Étage', 'Kser mur', 'Coeff. équivalence', 'Kser,p,ref'])
        self.data_Kser_walls.index.name = 'Mur'

        for etage, sys in self.data.items():
            for sys_name, sys_wall in sys.items():
                k_sys_wall = [etage, 0, 0, 0]

                for wall_name, wall in sys_wall['Murs'].items():
                    k_wall = [etage, 0, 0, 0]

                    # on gère d'abord les panneaux en pile
                    already_reported = []
                    for top_panel_name in self._list_on_top_of:
                        if top_panel_name not in wall['Panneaux'].keys():
                            continue

                        top_panel = wall['Panneaux'][top_panel_name]
                        bottom_panel_name = top_panel['On top of']
                        bottom_panel = wall['Panneaux'][bottom_panel_name]

                        # on vérifie que les deux panneaux sont des panneaux de contreventement, si ce n'est pas le cas alors, 
                        # les panneaux sont considérés comme des panneaux simples
                        if not top_panel['Contreventement'] or not bottom_panel['Contreventement']:
                            continue
                        kser_equi = self._get_Kser_panel_stack(top_panel_name, bottom_panel_name)
                        Kp_ref = self.data_Kser_panels.loc[bottom_panel_name]['Kser,p,ref']
                        coef_equi = kser_equi / Kp_ref
                        k_wall[1] = k_wall[1] + bottom_panel['Nombre'] * kser_equi
                        k_wall[2] = k_wall[2] + bottom_panel['Nombre'] * coef_equi
                        k_wall[3] = Kp_ref

                        already_reported.append(bottom_panel_name)
                        already_reported.append(top_panel_name)

                    # on gère ensuite les panneaux simples
                    for panel_name, panel in wall['Panneaux'].items():

                        if not panel['Contreventement'] or panel_name in already_reported:
                            continue
                        k_wall[1] = k_wall[1] + panel['Nombre'] * self.data_Kser_panels.loc[panel_name]['Kser,p unit.']
                        k_wall[2] = k_wall[2] + panel['Nombre'] * self.data_Kser_panels.loc[panel_name]['Coeff. équivalence']
                        k_wall[3] = self.data_Kser_panels.loc[panel_name]['Kser,p,ref']
                        k_sys_wall[3] = self.data_Kser_panels.loc[panel_name]['Kser,p,ref']
                    self.data_Kser_walls.loc[wall_name] = k_wall
                    k_sys_wall[1] = k_sys_wall[1] + k_wall[1]
                    k_sys_wall[2] = k_sys_wall[2] + k_wall[2]
                self.data_Kser_sys_walls.loc[sys_name] = k_sys_wall

        # Trier par index (nommés 'Système de mur' et 'Mur')
        self.data_Kser_sys_walls.sort_index(inplace=True)
        self.data_Kser_walls.sort_index(inplace=True)
        return self.data_Kser_sys_walls

    def _FivRd(self, b_i: si.mm, h: si.mm, s: si.mm, F_f_Rk: si.kN, coef_double_panel: float=1, K_mod: float=1, gamma_M: float=1.3):
        """
        Calcul la valeur de calcul (design) de la capacité résistante d'un panneau au contreventement.

        Args:
            b_i (si.mm): Largeur du panneau en mm.
            h (si.mm): Hauteur du panneau en mm.
            s (si.mm): Entraxe des organes en mm.
            F_f_Rk (si.kN): Capacité résistante caractéristique de cisaillement de l'organe en kN.
            K_mod (float, optional): Coefficient de modification Kmod. Defaults to 1.
            gamma_M (float, optional): Coefficient de matériaux. Par défaut 1.3 pour les assemblages.
        """
        if b_i >= h/2:
            c_i = 1
        else:
            c_i = b_i / (h/2)
        F_f_Rd = F_f_Rk * K_mod / gamma_M
        F_iv_Rd = (1.2 * coef_double_panel * F_f_Rd * b_i * c_i) / s
        return F_iv_Rd, c_i, coef_double_panel

    
    def _FtEd(self, wall_name: str):
        """
        Calcul l'effort de soulèvement d'un mur au contreventement.

        Args:
            wall_name (str): Nom du mur.
        """
        sys_wall = self.data_sys_walls[wall_name.split('_')[0]]
        wall = self.data_walls[wall_name]
        Fv_i_Ed_wall = self.data_walls_loads.loc[wall_name]['Fv,Ed,i ELU']
        return Fv_i_Ed_wall * sys_wall['Hauteur du système de mur'] / wall['Longueur']
        

    def taux_walls(self):
        """
        Calcul les taux de travail des MOB avec la vérification des critères ELS et ELU du contreventement.
        3 tables sont retournées:
            - taux_FvRd_panels: Taux de travail des panneaux avec:
                - Taux de vérification de la capacité de chaque panneau à reprendre l'effort de contreventement (%)
            - taux_FvRd_walls: Taux de travail des murs avec:
                - Taux de vérification de la capacité de chaque mur à reprendre l'effort de contreventement (%)
            - taux_FvRd_sys_walls: Taux de travail des systèmes de mur avec:
                - Taux de vérification de la capacité de chaque système de mur à reprendre l'effort de contreventement (%)
                - Taux du déplacement ELS de chaque système de mur (%)
        """
        self.taux_FvRd_panels = pd.DataFrame(columns=['Étage', 'Fiv,Ed', 'ci', 'Coeff. panneau double', 'Fiv,Rd', 'Taux cisaillement (%)'])
        self.taux_FvRd_panels.index.name = 'Panneau'
        self.taux_FvRd_walls = pd.DataFrame(columns=['Étage', 'Fw,v,Ed', 'Fw,v,Rd', 'Taux cisaillement (%)'])
        self.taux_FvRd_walls.index.name = 'Mur'
        self.taux_FvRd_sys_walls = pd.DataFrame(columns=['Étage', 'Fsw,v,Ed', 'Fsw,v,Rd', 'Taux cisaillement (%)', 'Déplacement', 'Critère ELS','Taux déplacement (%)'])
        self.taux_FvRd_sys_walls.index.name = 'Système de mur'

        self.get_Kser_wall()

        for etage, sys in self.data.items():
            for sys_name, sys_wall in sys.items():
                K_sw = self.data_Kser_sys_walls.loc[sys_name]['Kser système de mur']
                Fv_Rd_sys_walls = 0

                # Si des charges manuelles sont données
                direct_loads = sys_wall.get('Fv_Ed_ELU', None)
                if direct_loads:
                    F_sw_Ed = direct_loads
                    F_sw_Ed_ELS = sys_wall.get('Fv_Ed_ELS', 0*si.kN)
                    delta_v_sw_Ed = abs(K_sw * F_sw_Ed_ELS)
                else:
                    F_sw_Ed = self.data_sys_walls_loads.loc[sys_name]['Fv,Ed,i ELU']
                    delta_v_sw_Ed = abs(self.data_sys_walls_loads.loc[sys_name]['delta,v,Ed,i ELS'])

                for wall_name, wall in sys_wall['Murs'].items():
                    K_w = self.data_Kser_walls.loc[wall_name]['Kser mur']
                    Fv_Rd_walls = 0

                    if direct_loads:
                        F_w_Ed = direct_loads * K_w / K_sw
                        if not hasattr(self, "data_walls_loads"):
                            self.data_walls_loads = pd.DataFrame(index=self.data_walls.keys(), columns=["Étage", "Fv,Ed,i ELU", "Ft,Ed,i ELU (soulèvement)"])
                        self.data_walls_loads.loc[wall_name,'Fv,Ed,i ELU'] = F_w_Ed
                        self.data_walls_loads.loc[wall_name,'Ft,Ed,i ELU (soulèvement)'] = self._FtEd(wall_name)
                    else:
                        F_w_Ed = self.data_walls_loads.loc[wall_name]['Fv,Ed,i ELU']

                    for panel_name, panel in wall['Panneaux'].items():
                        if not panel['Contreventement']:
                            continue
                        
                        # Si des charges manuelles sont données
                        if direct_loads:
                            F_f_Ed = F_w_Ed * self.data_Kser_panels.loc[panel_name]['Kser,p unit.'] / self.data_Kser_walls.loc[wall_name]['Kser mur']
                        else:
                            F_f_Ed = self.data_panels_loads.loc[panel_name]['Fv,Ed,i ELU']

                        # Connecteur avec distance de couturage et capacité résistante caractéristique de cisaillement de l'organe
                        con = self._connecteurs[panel['Position']]
                        if panel['Position'] == 'Extérieur':
                            couturage_rive = wall['Couturage ext. rive']
                            F_f_Rk = wall['Ff,Rk ext']
                            K_mod = sqrt(con.beam_1._get_k_mod("Instantanee") * con.beam_2._get_k_mod("Instantanee"))
                        else:
                            couturage_rive = wall['Couturage int. rive']
                            F_f_Rk = wall['Ff,Rk int']
                            K_mod = sqrt(con.beam_1._get_k_mod("Instantanee") * con.beam_2._get_k_mod("Instantanee"))
                        
                        # Calcul de la valeur de calcul (design) de la capacité résistante d'un panneau au contreventement
                        F_f_Rd, c_i, coef_double_panel = self._FivRd(
                            panel['Largeur'], 
                            panel['Hauteur'], 
                            couturage_rive, 
                            F_f_Rk, 
                            self.COEFF_DOUBLE_PANEL[panel['Panneau double']], 
                            K_mod, 
                            gamma_M=1.3
                            )
                        
                        # Taux de travail panneau
                        self.taux_FvRd_panels.loc[panel_name] = [etage, F_f_Ed, c_i, coef_double_panel, F_f_Rd, round(abs(F_f_Ed * 100 / F_f_Rd), 2)]
                        Fv_Rd_walls =  Fv_Rd_walls + F_f_Rd * panel['Nombre']

                    # Taux de travail mur
                    self.taux_FvRd_walls.loc[wall_name] = [etage, F_w_Ed, Fv_Rd_walls, round(abs(F_w_Ed * 100 / Fv_Rd_walls), 2)]
                    Fv_Rd_sys_walls = Fv_Rd_sys_walls + Fv_Rd_walls

                # critère ELS
                crit_ELS = round(sys_wall['Hauteur étage'] / self.CRITERE_ELS, 2)
                # Taux de travail système de mur
                self.taux_FvRd_sys_walls.loc[sys_name] = [etage, F_sw_Ed, Fv_Rd_sys_walls, round(abs(F_sw_Ed * 100 / Fv_Rd_sys_walls), 2), delta_v_sw_Ed, crit_ELS, round(abs(delta_v_sw_Ed * 100 / crit_ELS), 2)]
                
        return self.taux_FvRd_walls

    def taux_ancrage(
        self,
        montant: object,
        lisse_basse: object,
        lisse_impl: object,
        Ft_Rk_wood: si.kN=0, 
        Ft_Rd_concrete: si.kN=0, 
        Fv_Rk_anc_lb: si.kN=0, 
        e_anc_lb: si.mm=0, 
        Fv_Rk_anc_li: si.kN=0, 
        e_anc_li: si.mm=0, 
        etage: str=None
    ):
        """ Calcul le taux de travail des ancrages, avec:
            - La vérification de l'equerre au soulèvement pour chaques murs
            - La vérification de l'ancrage de la lisse basse sur la lisse d'implantation
            - La vérification de l'ancrage de la lisse d'implantation au sol
        
        Args:
            montant (Barre): montant de l'ossature bois sur lequel l'équerre de traction est installée, 
                issue de la classe Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
            lisse_basse (Barre): lisse basse issue de la classe Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
            lisse_impl (Barre): lisse d'implantation issue de la classe Barre ou dérivé de cet objet provenant du module EC5_Element_droit.py
            Fv_Rk_anc_lb (float): Capacité résistante caractéristique de cisaillement de l'ancrage unitaire de la lisse basse en kN
                (attention à prendre le nef en compte).
            e_anc_lb (float): Entraxe de couturage de l'ancrage unitaire de la lisse basse en mm.
            Fv_Rk_anc_li (float): Capacité résistante caractéristique de cisaillement de l'ancrage unitaire de la lisse d'implantation en kN
                (attention à prendre le nef en compte).
            e_anc_li (float): Entraxe de couturage de l'ancrage unitaire de la lisse d'implantation en mm.
            etage (str): nom de l'étage à vérifier, si vide alors on vérifie tous les étages avec ces valeurs.
        """
        Ft_Rk_wood = Ft_Rk_wood * si.kN
        Ft_Rd_concrete = Ft_Rd_concrete * si.kN
        Fv_Rk_anc_lb = Fv_Rk_anc_lb * si.kN
        Fv_Rk_anc_li = Fv_Rk_anc_li * si.kN
        e_anc_lb = e_anc_lb * si.mm
        e_anc_li = e_anc_li * si.mm
        self.taux_ancr_walls = pd.DataFrame(columns=['Étage', 'Longueur mur', 'Ft,Ed', 'Ft,Rd', 'Taux ancr. soulèvement (%)', 'Fv,anc,Ed', 'Fv,anc,lb,Rd', 'Nbr ancr,lb', 'entraxe anc,lb',  'Taux ancr. lisse basse (%)', 'Fv,anc,li,Rd', 'Nbr ancr,li', 'entraxe anc,li', 'Taux ancr. lisse impl. (%)'])
        
        for wall_name, wall in self.data_walls_loads.iterrows():
            if etage and wall['Étage'] != etage:
                continue

            Ft_Ed = abs(wall['Ft,Ed,i ELU (soulèvement)'])
            Fv_Ed = abs(wall['Fv,Ed,i ELU'])
            K_mod_montant = montant._get_k_mod("Instantanee")
            K_mod_lb = lisse_basse._get_k_mod("Instantanee")
            K_mod_li = lisse_impl._get_k_mod("Instantanee")
            long_wall = self.data_walls[wall_name]['Longueur']
            nb_anc_lb = floor(long_wall / e_anc_lb)
            nb_anc_li = floor(long_wall / e_anc_li)
            Ft_Rd_wood = (Ft_Rk_wood * K_mod_montant / Assemblage.GAMMA_M_ASS)
            Fv_Rd_anc_lb = (Fv_Rk_anc_lb * K_mod_lb / Assemblage.GAMMA_M_ASS) * nb_anc_lb
            Fv_Rd_anc_li = (Fv_Rk_anc_li * K_mod_li / Assemblage.GAMMA_M_ASS) * nb_anc_li
            
            taux_soul = round(Ft_Ed * 100 / min(Ft_Rd_wood, Ft_Rd_concrete), 2)
            taux_lb = round(Fv_Ed * 100 / Fv_Rd_anc_lb, 2)
            taux_li = round(Fv_Ed * 100 / Fv_Rd_anc_li, 2)
            self.taux_ancr_walls.loc[wall_name] = [wall['Étage'], long_wall, Ft_Ed, min(Ft_Rd_wood, Ft_Rd_concrete), taux_soul, Fv_Ed, Fv_Rd_anc_lb, nb_anc_lb, e_anc_lb, taux_lb, Fv_Rd_anc_li, nb_anc_li, e_anc_li, taux_li]
        return self.taux_ancr_walls