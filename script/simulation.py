import itertools
from itertools import accumulate
import numpy as np
from IPython.display import Math
import seaborn as sns
import matplotlib.pyplot as plt
import time
from matplotlib.colors import ListedColormap
from collections import Counter
import random


def H_delta(liste,spin,d):
    E = 0
    liste_2 = changement(liste,spin)
    voisin = voisin_delta(spin,d,liste_2)
       
    for k in range(len(voisin)): 
        E += 2*liste_2[spin]*voisin[k]
   
    return E,liste_2

def voisin_delta(i, d, liste):

    voisinnage = []

    if i >= d:
        voisinnage.append(liste[i-d])


    if i < d * (d - 1):
        voisinnage.append(liste[i+d])


    if i % d != 0:
        voisinnage.append(liste[i-1])


    if (i+1) % d != 0:
        voisinnage.append(liste[i+1])

    return voisinnage


def H(liste, d):
    somme = 0

    for i in range(len(liste)):
        
        produit = 0
        voisinnage = voisin_droite_bas(i, d, liste)
        
        for voisin in voisinnage:
            produit += liste[i] * voisin
            
        somme += produit
    return -somme

def H_exact(beta, liste, d):
    somme = 0

    for i in range(len(liste)):
        
        produit = 0
        voisinnage = voisin_droite_bas(i, d, liste)
        
        for voisin in voisinnage:
            produit += liste[i] * voisin
            
        somme += produit
    return -beta*somme

def voisin_droite_bas(i, d, liste):

    voisinnage = []


    if i < d * (d - 1):
        voisinnage.append(liste[i+d])


    if (i+1) % d != 0:
        voisinnage.append(liste[i+1])

    return voisinnage

def Z(d, beta):
    somme = 0
    combinaison = list(itertools.product([-1, 1], repeat=d*d))
    for i in range(0,2**(d*d)):
        somme += np.exp(-H_exact(beta,combinaison[i],d))
    return somme

def mu(grille, d, beta):
    return (1/Z(d,beta))*np.exp(-H_exact(beta,grille,d))

def affichage_ising(matrice):
    n = len(matrice)
    m = len(matrice[0])

    matrice_np = (np.array(matrice) + 1) / 2
    sns.set(style='white')
    
    custom_cmap = ListedColormap(['blue', 'red'])
    ax = sns.heatmap(matrice_np, cmap=custom_cmap, square=True, linewidths=0.5, linecolor='gray', cbar=False, xticklabels=False, yticklabels=False, vmin=0, vmax=1)
    
    legend_text = "-1 (bleu) +1 (rouge)"
    ax.annotate(legend_text, xy=(0.5, 1.1), xycoords='axes fraction', fontsize=12, ha='center', va='center', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    plt.show()

def conversion_matrice(liste_position,d):
    return [liste_position[i:i + d] for i in range(0, len(liste_position), d)]

def gray_flip(liste):
    k = liste[0]
    if k > len(liste)-1:
        return k
        
    liste[k - 1] = liste[k]
    liste[k] = k + 1
    
    if k != 1:
        liste[0] = 1
    return k

def changement(vrai_liste,entier):        
    vrai_liste[entier] = -vrai_liste[entier]           
    return vrai_liste

def apercu_configuration(liste_f,d,beta):
    liste_special = [i for i in range(1,d*d+2)]
    print(liste_f,mu(liste_f,d,beta),H_exact(beta,liste_f,d))
    s = mu(liste_f,d,beta)
    for i in range(0,2**(d*d)-1):
        entierA = gray_flip(liste_special) - 1
        final = changement(liste_f,entierA)
        print(final,mu(final,d,beta),H_exact(beta,liste_f,d))
        s += mu(final,d,beta)
        #print(final,'--',H(final),'--',mu(final,d))
    print(Z(d,beta),s)
    return 

def affichage_configuration_1(entier,d):
    
    configuration_depart = [1 for i in range(1,d*d+1)]
    initialisation_gray_flip = [i for i in range(1,d*d+2)]
        
    for k in range(0,entier):
        spin_a_changer = gray_flip(initialisation_gray_flip) - 1
        configuration_depart = changement(configuration_depart,spin_a_changer)
        
    return configuration_depart   

def affichage_configuration(d):
    
    configuration_depart = [1 for i in range(1,d*d+1)]
    initialisation_gray_flip = [i for i in range(1,d*d+2)]
    stock = []
    stock.append(tuple(configuration_depart))
    
    for k in range(1, 2**(d*d)):
        spin_a_changer = gray_flip(initialisation_gray_flip) - 1
        configuration_depart = changement(configuration_depart,spin_a_changer)
        stock.append(tuple(configuration_depart))
        
    return stock

def enumeration_configuration(d, vecteur_beta):


    configuration_depart = np.ones((d*d,))
    initialisation_gray_flip = np.arange(1, d*d + 2)
    resultat = np.zeros((2**(d*d), len(vecteur_beta), 3))

    energie = H(configuration_depart, d)
    vecteur_energie_beta = np.full((len(vecteur_beta),), energie) * vecteur_beta
    cdf = np.exp(-vecteur_energie_beta)

    resultat[0, :, 0] = 0
    resultat[0, :, 1] = np.exp(-vecteur_energie_beta)
    resultat[0, :, 2] = cdf


    for i in range(1, 2**(d*d)):
        

        spin_a_changer = gray_flip(initialisation_gray_flip) - 1
        energie_apres_changement = H_delta(configuration_depart, spin_a_changer, d)
        vecteur_energie_beta = vecteur_energie_beta - np.full((len(vecteur_beta),), energie_apres_changement[0]) * vecteur_beta
        cdf += np.exp(-vecteur_energie_beta)


        resultat[i, :, 0] = i
        resultat[i, :, 1] = np.exp(-vecteur_energie_beta)
        resultat[i, :, 2] = cdf

    return resultat


stockage = {}
cles_traitees = []

def stockage_fonction(d, vecteur_beta):
    global stockage, cles_traitees
    

    cles_non_traitees = [(d, beta) for beta in vecteur_beta if (d, beta) not in cles_traitees]
    

    if not cles_non_traitees:
        return np.concatenate([stockage[(d, beta)] for beta in vecteur_beta], axis=1)
    

    nouveaux_resultats = enumeration_configuration(d, [beta for _, beta in cles_non_traitees])
    

    for i, cle in enumerate(cles_non_traitees):
        stockage[cle] = nouveaux_resultats[:, i, :].reshape(-1, 1, 3)
        cles_traitees.append(cle)
    

    return np.concatenate([stockage[(d, beta)] for beta in vecteur_beta], axis=1)


def simulation(d, beta, affichage=False):
    

    donnees = stockage_fonction(d, [beta])
    cdf = donnees[:, 0, 2]
    cdf = np.concatenate(([0], cdf))

    u = np.random.uniform(0, cdf[-1])

    for i in range(1, len(cdf)):  
        
        if u < cdf[i] and u >= cdf[i-1]:
            
            break

    if affichage:
        
        return affichage_configuration(i-1, d)
    
    else:
        return i-1
    
def magnetisation(liste):
    somme = 0
    for i in range(len(liste)):
        somme += liste[i]
    return somme

def metropolis_hastings(configuration_initiale, energie, d, beta, nombre_simulations):
     
        
    tentatives = 0
    acceptations = 0
    proba_acceptations = [0]
    configurations = [tuple(configuration_initiale)]
    liste_magnetisations = [magnetisation(configuration_initiale)]
    
    mesure_points = [0]
    start_time = time.time_ns()
    temps = [start_time]
    
    for i in range(nombre_simulations):

        liste_magnetisations.append(magnetisation(configuration_initiale))
        k = np.random.randint(0, len(configuration_initiale))
        h = sum(voisin_delta(k, d, configuration_initiale))
        
        
        delta_energie = 2 * h * configuration_initiale[k] * beta

        r = np.exp(-delta_energie)
        u = np.random.uniform(0, 1)
        tentatives += 1
        
        if u < r:
            configuration_initiale[k] = -configuration_initiale[k]
            energie += delta_energie
            
            acceptations += 1
        
            
        configurations.append(tuple(configuration_initiale))
        proba_acceptations.append(acceptations/tentatives)
        
        if (i+1)%200 == 0:         
            temps.append(time.time_ns())
            mesure_points.append(i)
    
    temps2 = np.interp(np.arange(nombre_simulations+1), mesure_points, temps) - start_time

    return configurations, liste_magnetisations, proba_acceptations, temps2, np.array(temps), mesure_points


def heat_bath(configuration_initiale, energie, d, beta, nombre_simulations):
        
    tentatives = 0
    acceptations = 0
    proba_acceptations = [0]
    configurations = [tuple(configuration_initiale)]
    liste_magnetisations = [magnetisation(configuration_initiale)]
    
    mesure_points = [0]
    start_time = time.time_ns()
    temps = [start_time]

    for i in range(nombre_simulations):
        
        liste_magnetisations.append(magnetisation(configuration_initiale))
        k = np.random.randint(0, len(configuration_initiale))
        h = sum(voisin_delta(k, d, configuration_initiale))
        
        
        delta_energie = 2 * h * configuration_initiale[k] * beta
        r = 1/(1 + np.exp(delta_energie))
        u = np.random.uniform(0, 1)
        tentatives += 1
        
        if u < r:
            configuration_initiale[k] = -configuration_initiale[k]
            energie += delta_energie
            
            acceptations += 1
            
        configurations.append(tuple(configuration_initiale))
        proba_acceptations.append(acceptations/tentatives)
        
        
        if (i+1)%200 == 0:         
            temps.append(time.time_ns())
            mesure_points.append(i)
    
    temps2 = np.interp(np.arange(nombre_simulations+1), mesure_points, temps) - start_time


    return configurations, liste_magnetisations, proba_acceptations, temps2, np.array(temps), mesure_points

def metropolis_hastings_etoile(configuration_initiale, energie, d, beta, nombre_simulations):

        
    tentatives = 0
    acceptations = 0
    proba_acceptations = [0]
    configurations = [tuple(configuration_initiale)]
    liste_magnetisations = [magnetisation(configuration_initiale)]
    

    for _ in range(nombre_simulations):
        liste_magnetisations.append(magnetisation(configuration_initiale))
        k = np.random.randint(0, len(configuration_initiale))
        h = sum(voisin_delta(k, d, configuration_initiale))
        
        
        delta_energie = 2 * h * configuration_initiale[k] * beta
        r = np.exp(-delta_energie)
        if delta_energie == 0:
            r = 1/2
        u = np.random.uniform(0, 1)
        tentatives += 1
        
        if u < r:
            configuration_initiale[k] = -configuration_initiale[k]
            energie += delta_energie
            
            acceptations += 1
            
        configurations.append(tuple(configuration_initiale))
        proba_acceptations.append(acceptations/tentatives)


    return configurations, liste_magnetisations, proba_acceptations

def autocorrelation(magnetisation_liste, tau_max):
    
    magnetisation_moyenne = np.mean(magnetisation_liste)
    variance = np.var(magnetisation_liste)   
    autocorr = np.zeros(tau_max)

    for tau in range(tau_max):
        s = 0
        for t in range(len(magnetisation_liste) - tau):
            s += (magnetisation_liste[t] - magnetisation_moyenne) * (magnetisation_liste[t + tau] - magnetisation_moyenne)
        autocorr[tau] = s / ((len(magnetisation_liste) - tau) * variance)

    return autocorr

def generer_liens(configuration, d, p):
    lien_droite = configuration[:,:-1] == configuration[:,1:]
    lien_bas = configuration[:-1,:] == configuration[1:,:]
    
    lien_droite = (np.random.random(lien_droite.shape) < 1) *   lien_droite
    lien_bas    = (np.random.random(lien_bas.shape) < 1)    *   lien_bas
    
    return lien_droite, lien_bas


def parcours_sw(i,j, d,  matrice, cluster, lien_droite, lien_bas, visites, beta):
    if i < 0 or j < 0 or i >= d or j >= d or visites[i][j]:
        return
    
    visites[i][j] = True 
    cluster.append(i*d + j)

    if np.random.uniform(0,1) < 1 - np.exp(-2*beta) and j < d-1 and lien_droite[i][j]:
        parcours_sw(i,j+1, d, matrice, cluster, lien_droite, lien_bas, visites, beta)
        
    if np.random.uniform(0,1) < 1 - np.exp(-2*beta) and i < d-1 and lien_bas[i][j]:
        parcours_sw(i+1,j, d, matrice, cluster, lien_droite, lien_bas, visites, beta)
        
    if np.random.uniform(0,1) < 1 - np.exp(-2*beta) and j > 0 and lien_droite[i][j - 1]:
        parcours_sw(i, j - 1, d, matrice, cluster, lien_droite, lien_bas, visites, beta)
        
    if np.random.uniform(0,1) < 1 - np.exp(-2*beta) and i > 0 and lien_bas[i - 1][j]:
        parcours_sw(i - 1, j, d, matrice, cluster, lien_droite, lien_bas, visites, beta)
        
def generer_clusters(matrice, d, lien_droite, lien_bas, beta):
    clusters = []
    visites = np.full((d, d), False)
    for i in range(d):
        for j in range(d):
            if not visites[i][j]:
                cluster = []
                parcours_sw(i,j,d,  matrice, cluster, lien_droite, lien_bas, visites, beta)
                if cluster:  
                    clusters.append(cluster)   
                    
    return clusters
                    
def retourner_spins(cluster, configuration, d ):
    configuration = configuration.reshape(d*d)
    for indices in cluster:
        u = np.random.uniform(0, 1)
        if u < 0.5:
            for indice in indices:
                configuration[indice] = -configuration[indice]
    return configuration   

def swendsen_wang(configuration_initiale, d, beta, nombre_simulations):
    liste_magnetisations = [magnetisation(configuration_initiale)]
    matrice = np.array(configuration_initiale).reshape(d,d)
    configurations = [tuple(configuration_initiale)]

    mesure_points = [0]
    start_time = time.time_ns()
    temps = [start_time]
    
    for t in range(nombre_simulations):
        liens = generer_liens(matrice, d, 1-np.exp(-2*beta))
        clusters  =  generer_clusters(matrice, d, liens[0], liens[1], beta)
        nouvelle_configuration = retourner_spins(clusters , matrice, d)
        liste_magnetisations.append(magnetisation(nouvelle_configuration))
        configurations.append(tuple(nouvelle_configuration))
        
        if (t+1)%10 == 0:         
            temps.append(time.time_ns())
            mesure_points.append(t)
            
    temps2 = np.interp(np.arange(nombre_simulations+1), mesure_points, temps) - start_time
    
    return configurations, liste_magnetisations, temps2, np.array(temps), mesure_points

def parcours(i, j, d, configuration, cluster, beta, graphe_sommet_visite):
    
    if  i < 0 or j < 0 or j >= d or i >= d or graphe_sommet_visite[i][j]:
        return
       
    graphe_sommet_visite[i][j] = True 
    cluster.append(i*d+j)

    if  np.random.uniform(0,1) < 1 - np.exp(-2*beta) and j < d-1 and configuration[i][j+1] == configuration[i][j]:
        parcours(i, j+1, d, configuration, cluster, beta, graphe_sommet_visite)
        
    if  np.random.uniform(0,1) < 1 - np.exp(-2*beta) and i < d-1 and configuration[i+1][j] == configuration[i][j]:
        parcours(i+1, j, d, configuration, cluster, beta, graphe_sommet_visite)
        
    if  np.random.uniform(0,1) < 1 - np.exp(-2*beta) and i > 0 and configuration[i-1][j] == configuration[i][j]:
        parcours(i-1, j, d, configuration, cluster, beta, graphe_sommet_visite)
        
    if  np.random.uniform(0,1) < 1 - np.exp(-2*beta) and j > 0 and configuration[i][j-1] == configuration[i][j]:
        parcours(i, j-1, d, configuration, cluster, beta, graphe_sommet_visite)
        
    return cluster



def wolff(configuration_depart, d, beta, nombre_simulations):
    configurations = [tuple(configuration_depart)]
    liste_magnetisations = [magnetisation(configuration_depart)]
 
    
    mesure_points = [0]
    start_time = time.time_ns()
    temps = [start_time]
    
    for t in range(nombre_simulations):

        graphe_sommet_visite = np.full((d, d), False)
        cluster = []
        matrice = np.array(configuration_depart).reshape(d,d)
        i, j  = np.random.randint(0, d), np.random.randint(0, d)
        clst = parcours(i, j, d, matrice, cluster, beta, graphe_sommet_visite)
        if clst:
            configuration_depart = matrice.reshape(d*d)
            for sommet in clst:
                configuration_depart[sommet] = -configuration_depart[sommet]
                

                
        configurations.append(tuple(configuration_depart))
        liste_magnetisations.append(magnetisation(configuration_depart))
        
        if (t+1)%10 == 0:         
            temps.append(time.time_ns())
            mesure_points.append(t)
            
    temps2 = np.interp(np.arange(nombre_simulations+1), mesure_points, temps) - start_time
    
    return configurations, liste_magnetisations, temps2, np.array(temps), mesure_points

def autocorrelation_moyenne_mh(d, beta, iterations, nombre_iterations_corr, tau):
    autocorr_values = np.zeros((nombre_iterations_corr, tau))
    for i in range(nombre_iterations_corr):
        configuration_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        autocorr_values[i, :] = autocorrelation(metropolis_hastings(configuration_depart, H_exact(beta,configuration_depart,d),d,beta,iterations)[1], tau)
    mean_autocorr = autocorr_values.mean(axis=0)

    return mean_autocorr

def autocorrelation_moyenne_hb(d, beta, iterations, nombre_iterations_corr, tau):
    autocorr_values = np.zeros((nombre_iterations_corr, tau))
    for i in range(nombre_iterations_corr):
        configuration_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        autocorr_values[i, :] = autocorrelation(heat_bath(configuration_depart, H_exact(beta,configuration_depart,d),d,beta,iterations)[1], tau)
    mean_autocorr = autocorr_values.mean(axis=0)

    return mean_autocorr

def autocorrelation_moyenne_swendsen_wang(d, beta, iterations, nombre_iterations_corr, tau):
    autocorr_values = np.zeros((nombre_iterations_corr, tau))
    for i in range(nombre_iterations_corr):
        configuration_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        autocorr_values[i, :] = autocorrelation(swendsen_wang(configuration_depart, d, beta, iterations)[1],tau)
    mean_autocorr = autocorr_values.mean(axis=0)
    
    return mean_autocorr

def autocorrelation_moyenne_wolff(d, beta, iterations, nombre_iterations_corr, tau):
    autocorr_values = np.zeros((nombre_iterations_corr, tau))
    for i in range(nombre_iterations_corr):
        configuration_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        autocorr_values[i, :] = autocorrelation(wolff(configuration_depart, d, beta, iterations)[1],tau)
    mean_autocorr = autocorr_values.mean(axis=0)
    
    return mean_autocorr

def affichage_autocorrelation_moyenne_wf_sw(d,configuration_depart,nombre_iterations_corr,tau, iterations, beta):
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))   
    mh = autocorrelation_moyenne_wolff(d, beta, iterations, nombre_iterations_corr, tau)
    sw = autocorrelation_moyenne_swendsen_wang(d, beta, iterations, nombre_iterations_corr, tau)
    plt.plot(mh,label="Wolff", alpha=0.7, color ='orange')
    plt.plot(sw,label="Swendsen_wang", alpha=0.7, color='purple')
    plt.ylabel("Autocorrélation moyenne")
    plt.title(f"Autocorrelation moyenne, grille de taille {d}*{d}, beta = {beta}")
    plt.legend()

    plt.show()
    
def affichage_autocorrelation_moyenne_mh_hb(d,configuration_depart,nombre_iterations_corr,tau, iterations, beta):
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))    
    mh = autocorrelation_moyenne_mh(d, beta,iterations, nombre_iterations_corr, tau)
    sw = autocorrelation_moyenne_hb(d, beta, iterations, nombre_iterations_corr, tau)
    plt.plot(mh,label="Metropolis Hasting", alpha=0.7)
    plt.plot(sw,label="Heath Bath", alpha=0.7, color='green')
    plt.ylabel("Autocorrélation moyenne")
    plt.title(f"Autocorrelation moyenne, grille de taille {d}*{d}, beta = {beta}")
    plt.legend()
 
    plt.show()


def affichage_autocorrelation_moyenne_mh_sw(d,configuration_depart,nombre_iterations_corr,tau, iterations,beta):
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))    
    mh = autocorrelation_moyenne_mh(d, beta, iterations, nombre_iterations_corr, tau)
    sw = autocorrelation_moyenne_swendsen_wang(d, beta, iterations, nombre_iterations_corr, tau)
    plt.plot(mh,label="Metropolis-Hastings", alpha=0.7, color = 'orange')
    plt.plot(sw,label="Swendsen_wang", alpha=0.7, color='purple')
    plt.ylabel("Autocorrélation moyenne")
    plt.title(f"Autocorrelation moyenne, grille de taille {d}*{d}, beta = {beta}")
    plt.legend()

    plt.show()
    
def erreur_moyenne_cluster(d, beta, nombre_simulations, nombre_repetition, nom_fichier):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])

    somme_erreurs_wf = np.zeros(nombre_simulations)
    somme_erreurs_sw = np.zeros(nombre_simulations)
    #somme_erreurs_carre_mh = np.zeros(nombre_simulations)
    #somme_erreurs_carre_hb = np.zeros(nombre_simulations)
    #print(configurations_exactes, '\n', frequences_exactes)
    
    for k in range(nombre_repetition):
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_wf = wolff(config_depart, d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_sw = swendsen_wang(config_depart, d, beta, nombre_simulations)[0]

        dist_wf = np.zeros(len(configurations_exactes))
        dist_sw = np.zeros(len(configurations_exactes))

        for t in range(1, nombre_simulations + 1):
            config_wf = tuple(configurations_wf[t - 1])
            config_sw = tuple(configurations_sw[t - 1])

            idx_wf = configurations_exactes[config_wf]
            idx_sw = configurations_exactes[config_sw]

            dist_wf[idx_wf] += 1
            dist_sw[idx_sw] += 1

            norm_dist_wf = dist_wf / t
            norm_dist_sw = dist_sw / t

            erreur_wf = np.sqrt(np.sum((norm_dist_wf - frequences_exactes) ** 2))
            erreur_sw = np.sqrt(np.sum((norm_dist_sw - frequences_exactes) ** 2))

            somme_erreurs_wf[t-1] += erreur_wf
            somme_erreurs_sw[t-1] += erreur_sw
            
            #print('** WF ** \n norm_distance', norm_dist_wf, '\n erreur : ', erreur_wf)
            #print('** SW ** \n norm_distance', norm_dist_sw, '\n erreur : ', erreur_sw)

            #somme_erreurs_carre_mh[t-1] += erreur_mh ** 2
            #somme_erreurs_carre_hb[t-1] += erreur_hb ** 2

    moyenne_erreurs_wf = somme_erreurs_wf / nombre_repetition
    moyenne_erreurs_sw = somme_erreurs_sw / nombre_repetition
    #variance_erreurs_mh = somme_erreurs_carre_mh / nombre_repetition - np.square(moyenne_erreurs_mh)
    #variance_erreurs_hb = somme_erreurs_carre_hb / nombre_repetition - np.square(moyenne_erreurs_hb)


    #variance_erreurs_mh = np.maximum(variance_erreurs_mh, 0)
    #variance_erreurs_hb = np.maximum(variance_erreurs_hb, 0)

    #ecart_type_mh = np.sqrt(variance_erreurs_mh)
    #ecart_type_hb = np.sqrt(variance_erreurs_hb)

    #intervalle_confiance_mh = 1.96 * ecart_type_mh
    #intervalle_confiance_hb = 1.96 * ecart_type_hb

    np.savetxt(nom_fichier, (moyenne_erreurs_wf, moyenne_erreurs_sw))
    
    
def erreur_moyenne_tout(d, beta, nombre_simulations, nombre_repetition, nom_fichier):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])

    somme_erreurs_wf = np.zeros(nombre_simulations+1)
    somme_erreurs_sw = np.zeros(nombre_simulations+1)
    somme_erreurs_mh = np.zeros(nombre_simulations+1)
    somme_erreurs_hb = np.zeros(nombre_simulations+1)
    #somme_erreurs_carre_mh = np.zeros(nombre_simulations)
    #somme_erreurs_carre_hb = np.zeros(nombre_simulations)
    #print(configurations_exactes, '\n', frequences_exactes)
    
    for k in range(nombre_repetition):
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_wf = wolff(config_depart, d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_sw = swendsen_wang(config_depart, d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_mh = metropolis_hastings(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_hb = heat_bath(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]

        dist_mh = np.zeros(len(configurations_exactes))
        dist_hb = np.zeros(len(configurations_exactes))

        dist_wf = np.zeros(len(configurations_exactes))
        dist_sw = np.zeros(len(configurations_exactes))

        for i in range(0, nombre_simulations + 1):
            
            config_wf = tuple(configurations_wf[i])
            config_sw = tuple(configurations_sw[i])
            config_mh = tuple(configurations_mh[i])
            config_hb = tuple(configurations_hb[i])

            idx_wf = configurations_exactes[config_wf]
            idx_sw = configurations_exactes[config_sw]
            idx_mh = configurations_exactes[config_mh]
            idx_hb = configurations_exactes[config_hb]

            dist_wf[idx_wf] += 1
            dist_sw[idx_sw] += 1
            dist_mh[idx_mh] += 1
            dist_hb[idx_hb] += 1

            norm_dist_wf = dist_wf / (i+1)
            norm_dist_sw = dist_sw / (i+1)
            norm_dist_mh = dist_mh / (i+1)
            norm_dist_hb = dist_hb / (i+1)

            erreur_wf = np.sqrt(np.sum((norm_dist_wf - frequences_exactes) ** 2))
            erreur_sw = np.sqrt(np.sum((norm_dist_sw - frequences_exactes) ** 2))
            erreur_mh = np.sqrt(np.sum((norm_dist_mh - frequences_exactes) ** 2))
            erreur_hb = np.sqrt(np.sum((norm_dist_hb - frequences_exactes) ** 2))

            somme_erreurs_wf[i] += erreur_wf
            somme_erreurs_sw[i] += erreur_sw
            somme_erreurs_mh[i] += erreur_mh
            somme_erreurs_hb[i] += erreur_hb
            


    moyenne_erreurs_wf = somme_erreurs_wf / nombre_repetition
    moyenne_erreurs_sw = somme_erreurs_sw / nombre_repetition
    moyenne_erreurs_mh = somme_erreurs_mh / nombre_repetition
    moyenne_erreurs_hb = somme_erreurs_hb / nombre_repetition

    np.savetxt(nom_fichier, (moyenne_erreurs_mh, moyenne_erreurs_hb, moyenne_erreurs_sw, moyenne_erreurs_wf))
    
def affichage_erreur_moyenne_tout(moyenne_erreurs_mh, moyenne_erreurs_hb, moyenne_erreurs_sw, moyenne_erreurs_wf, d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))
    plt.plot(moyenne_erreurs_wf, label="Wolff", alpha=0.7, color = 'orange')
    plt.plot(moyenne_erreurs_mh, label="Metropolis-Hastings", alpha=0.7)
    plt.plot(moyenne_erreurs_hb, label="Heat Bath", alpha=0.7, color='green')
    plt.plot(moyenne_erreurs_sw, label="Swendsen Wang", alpha=0.7, color='purple')
    plt.xlabel("Itérations")
    plt.ylabel("Erreur")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.title(f"Erreur moyenne en fonction du nombre d'itérations (d={d}, beta={beta})")
    plt.show()
    

def affichage_erreur_moyenne_cluster(moyenne_erreurs_mh, moyenne_erreurs_hb, d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))
    plt.plot(moyenne_erreurs_mh, label="Wolff", alpha=0.7, color = 'orange')
    #plt.fill_between(range(nombre_simulations), moyenne_erreurs_mh - intervalle_confiance_mh, moyenne_erreurs_mh + intervalle_confiance_mh, alpha=0.2)
    plt.plot(moyenne_erreurs_hb, label="Swendsen Wang", alpha=0.7, color='purple')
    #plt.fill_between(range(nombre_simulations), moyenne_erreurs_hb - intervalle_confiance_hb, moyenne_erreurs_hb + intervalle_confiance_hb, alpha=0.2, color='green')
    plt.xlabel("Itérations")
    plt.ylabel("Erreur")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.title(f"Erreur moyenne en fonction du nombre d'itérations (d={d}, beta={beta})")
    plt.show()
    
def erreur_moyenne_single_flip(d, beta, nombre_simulations, nombre_repetition, nom_fichier):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])

    somme_erreurs_mh = np.zeros(nombre_simulations)
    somme_erreurs_hb = np.zeros(nombre_simulations)
    #somme_erreurs_carre_mh = np.zeros(nombre_simulations)
    #somme_erreurs_carre_hb = np.zeros(nombre_simulations)

    for k in range(nombre_repetition):
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]

        configurations_mh = metropolis_hastings(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        configurations_hb = heat_bath(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]

        dist_mh = np.zeros(len(configurations_exactes))
        dist_hb = np.zeros(len(configurations_exactes))

        for t in range(1, nombre_simulations + 1):
            config_mh = tuple(configurations_mh[t - 1])
            config_hb = tuple(configurations_hb[t - 1])

            idx_mh = configurations_exactes[config_mh]
            idx_hb = configurations_exactes[config_hb]

            dist_mh[idx_mh] += 1
            dist_hb[idx_hb] += 1

            norm_dist_mh = dist_mh / t
            norm_dist_hb = dist_hb / t

            erreur_mh = np.sqrt(np.sum((norm_dist_mh - frequences_exactes) ** 2))
            erreur_hb = np.sqrt(np.sum((norm_dist_hb - frequences_exactes) ** 2))

            somme_erreurs_mh[t-1] += erreur_mh
            somme_erreurs_hb[t-1] += erreur_hb
            #somme_erreurs_carre_mh[t-1] += erreur_mh ** 2
            #somme_erreurs_carre_hb[t-1] += erreur_hb ** 2

    moyenne_erreurs_mh = somme_erreurs_mh / nombre_repetition
    moyenne_erreurs_hb = somme_erreurs_hb / nombre_repetition
    #variance_erreurs_mh = somme_erreurs_carre_mh / nombre_repetition - np.square(moyenne_erreurs_mh)
    #variance_erreurs_hb = somme_erreurs_carre_hb / nombre_repetition - np.square(moyenne_erreurs_hb)


    #variance_erreurs_mh = np.maximum(variance_erreurs_mh, 0)
    #variance_erreurs_hb = np.maximum(variance_erreurs_hb, 0)

    #ecart_type_mh = np.sqrt(variance_erreurs_mh)
    #ecart_type_hb = np.sqrt(variance_erreurs_hb)

    #intervalle_confiance_mh = 1.96 * ecart_type_mh
    #intervalle_confiance_hb = 1.96 * ecart_type_hb

    np.savetxt(nom_fichier, (moyenne_erreurs_mh, moyenne_erreurs_hb))
    

def affichage_erreur_moyenne_single_flip(moyenne_erreurs_mh, moyenne_erreurs_hb, d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))
    plt.plot(moyenne_erreurs_mh, label="Metropolis-Hastings", alpha=0.7)
    #plt.fill_between(range(nombre_simulations), moyenne_erreurs_mh - intervalle_confiance_mh, moyenne_erreurs_mh + intervalle_confiance_mh, alpha=0.2)
    plt.plot(moyenne_erreurs_hb, label="Heat Bath", alpha=0.7, color='green')
    #plt.fill_between(range(nombre_simulations), moyenne_erreurs_hb - intervalle_confiance_hb, moyenne_erreurs_hb + intervalle_confiance_hb, alpha=0.2, color='green')
    plt.xlabel("Itérations")
    plt.ylabel("Erreur")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.title(f"Erreur moyenne en fonction du nombre d'itérations (d={d}, beta={beta})")
    plt.show()
    
def erreur_moyenne_single_flip_1(d, beta, nombre_simulations, nombre_repetition, nom_fichier):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])

    somme_erreurs_mh = np.zeros(nombre_simulations+1)
    somme_erreurs_hb = np.zeros(nombre_simulations+1)


    for k in range(nombre_repetition):
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]

        configurations_mh = metropolis_hastings(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_hb = heat_bath(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]

        dist_mh = np.zeros(len(configurations_exactes))
        dist_hb = np.zeros(len(configurations_exactes))

        for i in range(0, nombre_simulations + 1):
            config_mh = tuple(configurations_mh[i])
            config_hb = tuple(configurations_hb[i])

            idx_mh = configurations_exactes[config_mh]
            idx_hb = configurations_exactes[config_hb]

            dist_mh[idx_mh] += 1
            dist_hb[idx_hb] += 1

            norm_dist_mh = dist_mh / (i+1)
            norm_dist_hb = dist_hb / (i+1)

            erreur_mh = np.sqrt(np.sum((norm_dist_mh - frequences_exactes) ** 2))
            erreur_hb = np.sqrt(np.sum((norm_dist_hb - frequences_exactes) ** 2))

            somme_erreurs_mh[i] += erreur_mh
            somme_erreurs_hb[i] += erreur_hb


    moyenne_erreurs_mh = somme_erreurs_mh / nombre_repetition
    moyenne_erreurs_hb = somme_erreurs_hb / nombre_repetition
    

    np.savetxt(nom_fichier, (moyenne_erreurs_mh, moyenne_erreurs_hb))
    
def erreur_moyenne_wf_mh(d, beta, nombre_simulations, nombre_repetition, nom_fichier):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])

    somme_erreurs_mh = np.zeros(nombre_simulations+1)
    somme_erreurs_hb = np.zeros(nombre_simulations+1)


    for k in range(nombre_repetition):
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]

        configurations_mh = metropolis_hastings(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_hb = wolff(config_depart, d, beta, nombre_simulations)[0]

        dist_mh = np.zeros(len(configurations_exactes))
        dist_hb = np.zeros(len(configurations_exactes))

        for i in range(0, nombre_simulations + 1):
            config_mh = tuple(configurations_mh[i])
            config_hb = tuple(configurations_hb[i])

            idx_mh = configurations_exactes[config_mh]
            idx_hb = configurations_exactes[config_hb]

            dist_mh[idx_mh] += 1
            dist_hb[idx_hb] += 1

            norm_dist_mh = dist_mh / (i+1)
            norm_dist_hb = dist_hb / (i+1)

            erreur_mh = np.sqrt(np.sum((norm_dist_mh - frequences_exactes) ** 2))
            erreur_hb = np.sqrt(np.sum((norm_dist_hb - frequences_exactes) ** 2))

            somme_erreurs_mh[i] += erreur_mh
            somme_erreurs_hb[i] += erreur_hb


    moyenne_erreurs_mh = somme_erreurs_mh / nombre_repetition
    moyenne_erreurs_hb = somme_erreurs_hb / nombre_repetition
    

    np.savetxt(nom_fichier, (moyenne_erreurs_mh, moyenne_erreurs_hb))

def affichage_erreur_wf_mh(moyenne_erreurs_mh, moyenne_erreurs_hb, d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))
    plt.plot(moyenne_erreurs_mh, label="Metropolis-Hastings", alpha=0.7)
    plt.plot(moyenne_erreurs_hb, label="Wolff", alpha=0.7, color='orange')
    plt.xlabel("Itérations")
    plt.ylabel("Erreur")
    plt.legend()
    plt.title(f"Erreur moyenne en fonction du nombre d'itérations (d={d}, beta={beta})")
    plt.show()
    
def erreur_moyenne_single_flip_etoile(d, beta, nombre_simulations, nombre_repetition, nom_fichier):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])

    somme_erreurs_mh = np.zeros(nombre_simulations+1)
    somme_erreurs_hb = np.zeros(nombre_simulations+1)
    somme_erreurs_mh_ = np.zeros(nombre_simulations+1)


    for k in range(nombre_repetition):
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_mh = metropolis_hastings(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_hb = heat_bath(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        
        config_depart = [random.choice([-1, 1]) for _ in range(d*d)]
        configurations_mh_ = metropolis_hastings_etoile(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)[0]
        
        

        dist_mh = np.zeros(len(configurations_exactes))
        dist_hb = np.zeros(len(configurations_exactes))
        dist_mh_ = np.zeros(len(configurations_exactes))

        for i in range(0, nombre_simulations + 1):
            config_mh = tuple(configurations_mh[i])
            config_hb = tuple(configurations_hb[i])
            config_mh_ = tuple(configurations_mh_[i])

            idx_mh = configurations_exactes[config_mh]
            idx_hb = configurations_exactes[config_hb]
            idx_mh_ = configurations_exactes[config_mh_]

            dist_mh[idx_mh] += 1
            dist_mh_[idx_mh_] += 1
            dist_hb[idx_hb] += 1

            norm_dist_mh = dist_mh / (i+1)
            norm_dist_mh_ = dist_mh_ / (i+1)
            norm_dist_hb = dist_hb / (i+1)

            erreur_mh = np.sqrt(np.sum((norm_dist_mh - frequences_exactes) ** 2))
            erreur_mh_ = np.sqrt(np.sum((norm_dist_mh_ - frequences_exactes) ** 2))
            erreur_hb = np.sqrt(np.sum((norm_dist_hb - frequences_exactes) ** 2))

            somme_erreurs_mh[i] += erreur_mh
            somme_erreurs_mh_[i] += erreur_mh_
            somme_erreurs_hb[i] += erreur_hb


    moyenne_erreurs_mh = somme_erreurs_mh / nombre_repetition
    moyenne_erreurs_mh_ = somme_erreurs_mh_ / nombre_repetition
    moyenne_erreurs_hb = somme_erreurs_hb / nombre_repetition
    

    np.savetxt(nom_fichier, (moyenne_erreurs_mh, moyenne_erreurs_hb, moyenne_erreurs_mh_))
    
    
def affichage_erreur_moyenne_etoile(moyenne_erreurs_mh, moyenne_erreurs_hb, moyenne_erreurs_mh_, d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))
    plt.plot(moyenne_erreurs_mh, label="Metropolis-Hastings", alpha=0.7)
    plt.plot(moyenne_erreurs_hb, label="Heat Bath", alpha=0.7, color='green')
    plt.plot(moyenne_erreurs_mh_, label="Metropolis-Hastings *", alpha=0.7, color='red')
   
    plt.xlabel("Itérations")
    plt.ylabel("Erreur")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.title(f"Erreur moyenne en fonction du nombre d'itérations (d={d}, beta={beta})")
    plt.show()

def affichage_erreur_moyenne_single_flip_1(moyenne_erreurs_mh, moyenne_erreurs_hb,  d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(7, 5))
    plt.plot(moyenne_erreurs_mh, label="Metropolis-Hastings", alpha=0.7)
    plt.plot(moyenne_erreurs_hb, label="Heat Bath", alpha=0.7, color='green')
    plt.xlabel("Itérations")
    plt.ylabel("Erreur moyenne")
    #plt.yscale('log')
    #plt.xscale('log')
    plt.legend()
    plt.title(f"Convergence des trajectoires (d={d}, beta={beta})")
    plt.show()

def erreur_temps(d, beta, nombre_simulations, nombre_repetition, nom_fichier, T):
    
    enumerations = stockage_fonction(d, [beta])
    all_configs = affichage_configuration(d)
    enumeration_values = enumerations[:, :, 1][:, 0] / enumerations[2**(d*d)-1][0][2]

    stock = list(zip(all_configs, enumeration_values))
    configurations_exactes = {s[0]: idx for idx, s in enumerate(stock)}
    frequences_exactes = np.array([s[1] for s in stock])
    
    newarr = np.zeros(len(T))
    newarr_2 = np.zeros(len(T))


    for _ in range(nombre_repetition):
        
        somme_erreurs_wf  = np.zeros(nombre_simulations+1)
        somme_erreurs_mh  = np.zeros(nombre_simulations+1)
             
        config_depart     = [random.choice([-1, 1]) for _ in range(d*d)]
        dist_mh           = np.zeros(len(configurations_exactes))
        dist_wf           = np.zeros(len(configurations_exactes))
        
    
        simulation_mh     = metropolis_hastings(config_depart, H_exact(beta, config_depart, d), d, beta, nombre_simulations)
        configurations_mh = simulation_mh[0]
        somme_temps_mh    = simulation_mh[3]      
        
        
        
        simulation_wf     = wolff(config_depart, d, beta, nombre_simulations)
        configurations_wf = simulation_wf[0]
        somme_temps_wf    = simulation_wf[2]      
        

        
        
        for t in range(1, nombre_simulations + 1):
            
            config_mh = tuple(configurations_mh[t - 1])
            idx_mh = configurations_exactes[config_mh]
            dist_mh[idx_mh] += 1
            norm_dist_mh = dist_mh / t
            erreur_mh = np.sqrt(np.sum((norm_dist_mh - frequences_exactes) ** 2))
            somme_erreurs_mh[t-1] += erreur_mh
            
            
            config_wf = tuple(configurations_wf[t - 1])
            idx_wf = configurations_exactes[config_wf]
            dist_wf[idx_wf] += 1
            norm_dist_wf = dist_wf / t
            erreur_wf = np.sqrt(np.sum((norm_dist_wf - frequences_exactes) ** 2))
            somme_erreurs_wf[t-1] += erreur_wf


        newarr  += np.interp(T,somme_temps_mh, somme_erreurs_mh)
        newarr_2  += np.interp(T,somme_temps_wf, somme_erreurs_wf)
        
        
    newarr = newarr[newarr != 0]
    newarr_2 = newarr_2[newarr_2 != 0]
        
        
    

    

    np.savetxt(nom_fichier, (newarr[0:-1]/nombre_repetition, newarr_2[0:-1]/nombre_repetition, T[0:-1]))
    
def affichage_erreur_temps(newarr,newarr_2, T, d, beta):
    
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 5))
    plt.plot(T, newarr, label   =  "Metropolis Hasting", alpha=0.7)
    plt.plot(T, newarr_2, label =  "Wolff"             , alpha=0.7, color = 'orange')
    plt.xlabel("Temps")
    plt.ylabel("Erreur")
    plt.yscale('log')
    plt.legend()
    plt.title(f"Erreur moyenne en fonction du temps (d={d}, beta={beta})")

    
def temps_moyen_iteration_wf(d, beta,nombre_simulations,nombre_iteration_moyenne):
    
    delta_it = np.zeros((int(nombre_simulations/10) ))
    
    for _ in range(nombre_iteration_moyenne):   
        config   = [random.choice([-1, 1]) for _ in range(d*d)]
        temps_it = wolff(config, d, beta, nombre_simulations)[3]
        delta_it  += np.diff(temps_it) 
    
    return np.mean(delta_it/nombre_iteration_moyenne)/10


def temps_moyen_iteration_mh(d, beta,nombre_simulations,nombre_iteration_moyenne):
    

    delta_it = np.zeros((int(nombre_simulations/1000) ))
    for _ in range(nombre_iteration_moyenne):   
        
        config   = [random.choice([-1, 1]) for _ in range(d*d)]
        temps_it = metropolis_hastings(config, H_exact(beta, config, d), d, beta, nombre_simulations)[4]
        delta_it  += np.diff(temps_it)    
    
    return np.mean(delta_it/nombre_iteration_moyenne)/1000
