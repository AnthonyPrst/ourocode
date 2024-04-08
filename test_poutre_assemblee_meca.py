
# %%
import sys
import os
import operator

from math import *
import inspect
import numpy as np
import pandas as pd
from copy import deepcopy
from matplotlib import pyplot as plt
import forallpeople as si
from PySide6.QtCore import Qt 
from PySide6.QtWidgets import QApplication, QInputDialog


sys.path.append(os.path.join(os.getcwd(), "catalog", "eurocode"))
import A0_Projet as A0
import A1_MEF as A1
import EC0_Combinaison as EC0
import EC5_Assemblage as EC5_Ass
import EC5_Element_droit as EC5_Ele

#############################   Vérification

ELU_STR = True
ELU_STR_ACC = False
ELS = True
#############################   Chargement
Charge1 = [1, '', 'Permanente G', 'Linéique', -10, '0/4000', 'Z']
Charge2 = [0, 'Poids propre', 'Permanente G', 'Linéique', -50, '0/4000', 'Z']
Charge3 = [2, '', 'Neige normale Sn', 'Linéique', -200, '0/4000', 'Z']
Charge4 = [3, '', 'Exploitation Q', 'Linéique', -150, '0/4000', 'Z']

_list_loads = [Charge2, Charge3]

#############################   Appuis
listdeplacement = [[1, "Rotule", 0, 0], [2, "Rotule", 4000, 0]]

#############################   Case à cocher
app = QApplication()
disposition, ok = QInputDialog.getItem(
    None, "Position des planches", "Position des planches :", ["Latérale", "Dessus / Dessous"], 0, False, flags=Qt.FramelessWindowHint)


type_batiment, ok = QInputDialog.getItem(
    None, "Type de bâtiment", "Type de bâtiment :", ["Bâtiments courants", "Bâtiments agricoles et similaires"], 0, False, flags=Qt.FramelessWindowHint)

#############################   Dalle
type_element = "Élément structuraux"
longueur = 4000 * si.mm
disposition
classe_service = 1
humidite_initiale = 12 #pourcent
humidite_finale = 12 #pourcent
recouvrement = -40 #mm
entraxe_connecteurs = 200 #mm

#############################   Maillage MEF
def calcs_number_node(long):
        """Retourne le nombre d'éléments pour le maillage MEF

        Returns:
            int: nombre d'éléments MEF
        """
        round_long = ceil(long/10)*10
        
        if round_long <= 1000:
            ele = int(round_long/10)
        elif round_long <= 10000:
            ele = int(round_long/100)
        else:
            ele = int(round_long/200)
        return ele
    
node = 600
#node = calcs_number_node(longueur.value*10**3)

#############################   Planche basse
b_planche_basse = 77 * si.mm
h_planche_basse = 360 * si.mm
classe_bois_planche_basse = "C24"

projet = A0.Projet(code_INSEE=73215, alt= 1200)
chargement = EC0.Chargement._from_parent_class(projet)
chargement.create_load_by_list(_list_loads)

barre2 = EC5_Ass.Barre._from_parent_class(projet, b=b_planche_basse.value*10**3, h=h_planche_basse.value*10**3, section="Rectangulaire", Hi=humidite_initiale, Hf=humidite_finale, classe=classe_bois_planche_basse, cs=classe_service, effet_systeme="False")

#############################   Planche haute
b_planche_intermediaire = 83 * si.mm
h_planche_intermediaire = 120 * si.mm
classe_bois_planche_intermediaire = "C24"

barre3 = EC5_Ass.Barre._from_parent_class(projet, b=b_planche_intermediaire.value*10**3, h=h_planche_intermediaire.value*10**3, section="Rectangulaire", Hi=humidite_initiale, Hf=humidite_finale, classe=classe_bois_planche_intermediaire, cs=classe_service, effet_systeme="False")
pd.DataFrame(barre3.caract_meca).T

#############################   Connecteur
d = 6 * si.mm
d1 = 3.9 * si.mm
ds = 4.3 * si.mm
dh = 12 * si.mm
l = 200 * si.mm
rho_a = 350 * si.kg / si.m**3
fhead = 10.5 * si.MPa
ftensk = 11300 * si.N
f_ax_k = 11.7 *si.MPa
MyRk = 9494 * si.N*si.mm
alpha1 = 90 #°
alpha2 = 90 #°


assemblage = EC5_Ass.Assemblage._from_parent_class(projet, beam_1=barre3, beam_2=barre2, nfile=1, nCis=2)

#############   Vérification vis-boulon ou vis-pointe selon d_ef
if d1.value*10**3*1.1 > 6:
    connecteur = EC5_Ass.Tirefond_sup_6._from_parent_class(assemblage, d=d.value*10**3, d1=d1.value*10**3, ds=ds.value*10**3, dh=dh.value*10**3, l=l.value*10**3, n=1, rho_a=rho_a.value, fhead=fhead.value*10**-6, ftensk=ftensk.value, MyRk=MyRk.value*10**3, alpha1=alpha1, alpha2=alpha2)
else:
    connecteur = EC5_Ass.Tirefond_inf_7._from_parent_class(assemblage, d=d.value*10**3, d1=d1.value*10**3, ds=ds.value*10**3, dh=dh.value*10**3, l=l.value*10**3, n=1, rho_a=rho_a.value, fhead=fhead.value*10**-6, ftensk=ftensk.value, MyRk=MyRk.value*10**3, alpha1=alpha1, alpha2=alpha2, percage = True)


#############   Kser par vis et par plan de cisaillement
Kser = connecteur.Kser

#############   Kser assemblage
Kser_ass = connecteur.Kser_ass
dalle = EC5_Ele.Poutre_assemblee_meca(beam_2=barre2, l=longueur.value*10**3, disposition=disposition, recouvrement=[0,recouvrement], Kser=[None,None,Kser_ass[1]], entraxe=[None,None,entraxe_connecteurs], psy_2=0, beam_3=barre3)
combinaison = A1.Combinaison._from_parent_class(chargement, ELU_STR=ELU_STR, ELU_STR_ACC=ELU_STR_ACC, ELS_C=ELS, ELS_QP=ELS, cat="Cat A : habitation", kdef=dalle.K_def[1])

#############   Kser final tenant compte du fluage
Kser_fin = dalle.Kser_fin

#############   Portance locale
fh1k = connecteur.fh1k
fh2k = connecteur.fh2k

#############   Nombre efficace
nef = connecteur.nef(a1_beam1=entraxe_connecteurs, a1_beam2=entraxe_connecteurs)
connecteur.nfile

#############   Résistance axiale
FaxRk = connecteur.FaxRk(faxk=f_ax_k, l_ef=connecteur.t1, alpha=alpha1, beam=barre3)

#############   Capacité résistante caractéristique par tige et par plan de cisaillement
FvRk = connecteur.FvRk(True)

#############   Gamma i
gamma_i = dalle.gamma_i
gamma_2 = gamma_i["gamma 2"]
gamma_3 = gamma_i["gamma 3"][0]

#############   EI efficace dalle
ei_eff = dalle.EI_eff


l = int(ceil(dalle.l.value * 10**3))
E_mean_fin = int(barre2.E_mean_fin.value*10**-6)
A = barre2.b * barre2.h + barre3.b * barre3.h
A = A.value*10**6
iy = ei_eff[1].value*10**6 / E_mean_fin
G = float(barre2.caract_meca.Gmoy)



dict_combi = {}
dict_taux_max = {type_verif: {"combinaison": "", "taux": 0} for type_verif in ("cisaillement", "Winst(Q)", "Wnet,fin")}
dict_taux_max_flexion = {dim_barre:{"combinaison": "", "taux": 0} for dim_barre in ("planche basse", "planche intermédiaire")}
        
for i, barre in enumerate([barre2,barre3]):
    index_barre = i+2
    if i==0:
        dim_barre = "planche basse"
    else:
        dim_barre = "planche intermédiaire"
        
        
    for combi in combinaison.get_list_combination():
        combinaison.get_combi_list_load(combi)
        dict_combi[combi]= {}
        dict_combi[combi]["EF"] = A1.MEF._from_parent_class(combinaison, long=l, E=E_mean_fin, A=A, G=G, J=1, Iy=iy, Iz=1, ele=node, alphaZ=0, alphaY=0, alphaX=0)
        dict_combi[combi]["EF"].create_supports_by_list(listdeplacement)
        dict_combi[combi]["EF"].calcul_1D()
        
        if combi[0:7] == 'ELU_STR':
            action = dict_combi[combi]["EF"].min_type_load(combi)
            barre2.K_mod = barre2.K_mod_table[action].iloc[0]
            barre3.K_mod = barre3.K_mod_table[action].iloc[0]

            dictEiMinMax = dict_combi[combi]["EF"].effort_interne_max()
            max_XY_values = (dictEiMinMax["My_max"]["Position"], dictEiMinMax["My_max"]["Effort"])
            min_XY_values = (dictEiMinMax["My_min"]["Position"], dictEiMinMax["My_min"]["Effort"])
                
            dict_effort_max = dict_combi[combi]["EF"].effort_interne_max()
            V_z = max(dict_effort_max['Vz_max']['Effort'],abs(dict_effort_max['Vz_min']['Effort']))
            tau_2 = dalle.tau_2_max(V_z)

            M_y = max(dict_effort_max['My_max']['Effort'],abs(dict_effort_max['My_min']['Effort']))
            sigma_i = dalle.sigma_i(M_y, beam=index_barre)
            sigma_mi = dalle.sigma_mi(M_y, beam=index_barre)

            F_i = dalle.F_i(V_z, connecteur=2)

            dict_combi[combi]["flexion"] = EC5_Ele.Flexion._from_parent_class(barre, lo=longueur.value*10**3, coeflef=0.9, pos="Charge sur fibre comprimée")
            dict_combi[combi]["flexion"].sigma_m_rd = {'y': sigma_mi[1], 'z': 0 * si.MPa}
            fmd = dict_combi[combi]["flexion"].f_m_d(action, "Fondamentales")

            if sigma_i[1] < 0:
                dict_combi[combi]["compression"] = EC5_Ele.Compression._from_parent_class(barre,  lo_y=longueur.value*10**3, lo_z=longueur.value*10**3, type_appuis= "Rotule - Rotule")
                fc0d = dict_combi[combi]["compression"].f_c_0_d(action, "Fondamentales")
                dict_combi[combi]["compression"].sigma_c_0_rd = abs(sigma_i[1])
                taux_compression = dict_combi[combi]["compression"].taux_c_0_d()
                compression=dict_combi[combi]["compression"]
                traction=None
            else:
                dict_combi[combi]["traction"] = EC5_Ele.Traction._from_parent_class(barre)
                ft0d = dict_combi[combi]["traction"].f_t_0_d(action, "Fondamentales")
                dict_combi[combi]["traction"].sigma_t_0_rd = sigma_i[1]
                taux_traction = dict_combi[combi]["traction"].taux_t_0_d()
                traction=dict_combi[combi]["traction"]
                compression=None
                
            taux_flexion = dict_combi[combi]["flexion"].taux_m_d(compression=compression, traction=traction)
            print(taux_flexion)
            dict_combi[combi]["cisaillement"] = EC5_Ele.Cisaillement._from_parent_class(barre)
            dict_combi[combi]["cisaillement"].tau_rd = tau_2[1]
            fmd = dict_combi[combi]["cisaillement"].f_v_d(action, "Fondamentales")
            taux_cisaillement = dict_combi[combi]["cisaillement"].taux_tau_d()
            
            
            
            dict_taux_flexion = dict_combi[combi]["flexion"].taux_m_rd
            max_taux_flexion = max(dict_taux_flexion.items(), key=operator.itemgetter(1))[1]
            
            if max_taux_flexion > dict_taux_max_flexion[dim_barre]["taux"]:
                dict_taux_max_flexion[dim_barre]["combinaison"] = combi
                dict_taux_max_flexion[dim_barre]["taux"] = max_taux_flexion
            
            dict_taux_cisaillement = dict_combi[combi]["cisaillement"].taux_tau_rd
            max_taux_cisaillement = max(dict_taux_cisaillement.items(), key=operator.itemgetter(1))[1]
            
            if max_taux_cisaillement > dict_taux_max["cisaillement"]["taux"]:
                dict_taux_max["cisaillement"]["combinaison"] = combi
                dict_taux_max["cisaillement"]["taux"] = max_taux_cisaillement
                
        elif combi[0:6] == 'W_inst':
            dictUMinMax = dict_combi[combi]["EF"].deplacement_max()
            max_XY_values = (dictUMinMax["Uz_max"][0], dictUMinMax["Uz_max"][1])
            min_XY_values = (dictUMinMax["Uz_min"][0], dictUMinMax["Uz_min"][1])
            
            fleche_max = max(dictUMinMax['Uz_max'][0],abs(dictUMinMax['Uz_min'][1]))
            barre.fleche(long=l, Ed_WinstQ=fleche_max, type_ele=type_element, type_bat=type_batiment)
            if barre.taux_ELS["Winst(Q)"] > dict_taux_max["Winst(Q)"]["taux"]:
                dict_taux_max["Winst(Q)"]["combinaison"] = combi
                dict_taux_max["Winst(Q)"]["taux"] = barre.taux_ELS["Winst(Q)"]
        
        elif combi[0:9] == 'W_net_fin':
            dictUMinMax = dict_combi[combi]["EF"].deplacement_max()
            max_XY_values = (dictUMinMax["Uz_max"][0], dictUMinMax["Uz_max"][1])
            min_XY_values = (dictUMinMax["Uz_min"][0], dictUMinMax["Uz_min"][1])
            
            fleche_max = max(dictUMinMax['Uz_max'][0],abs(dictUMinMax['Uz_min'][1]))
            barre.fleche(long=l, Ed_Wnetfin=fleche_max, type_ele=type_element, type_bat=type_batiment)
            if barre.taux_ELS["Wnet,fin"] > dict_taux_max["Wnet,fin"]["taux"]:
                dict_taux_max["Wnet,fin"]["combinaison"] = combi
                dict_taux_max["Wnet,fin"]["taux"] = barre.taux_ELS["Wnet,fin"]
        
    print("Combinaison avec le taux le plus élevé en flexion de la", dict_taux_max_flexion)
        
print("Combinaisons avec les taux les plus élevés",dict_taux_max)

dict_combi[dict_taux_max["cisaillement"]["combinaison"]]["EF"].show_graphique_Vz(dict_taux_max["cisaillement"]["combinaison"])
dict_combi[dict_taux_max_flexion[dim_barre]["combinaison"]]["EF"].show_graphique_My(dict_taux_max_flexion[dim_barre]["combinaison"])
dict_combi[dict_taux_max["Winst(Q)"]["combinaison"]]["EF"].show_graphique_fleche(dict_taux_max["Winst(Q)"]["combinaison"])
dict_combi[dict_taux_max["Wnet,fin"]["combinaison"]]["EF"].show_graphique_fleche(dict_taux_max["Wnet,fin"]["combinaison"])
    
Fv_Rd = connecteur.F_Rd(FvRk[1].value*10**-3)
Taux_connecteur = connecteur.taux_cisaillement(Fv_Ed = F_i[1].value*10**-3)
print("Taux connnecteur :", Taux_connecteur[1])
app.exit()