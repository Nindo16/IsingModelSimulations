# Simulation du modèle d'Ising avec des algorithmes MCMC

## Introduction

Le modèle d'Ising est essentiel en physique statistique. Il est couramment utilisé pour décrire les transitions de phase dans les systèmes magnétiques. Ce projet vise à simuler le modèle d'Ising en utilisant des algorithmes Monte Carlo par chaîne de Markov (MCMC). La comparaison est effectuée en se basant sur la loi exacte obtenue en la simulant sur des grilles de petites tailles.

## Objectif du projet

L'objectif est de comparer différents algorithmes MCMC pour la simulation du modèle d'Ising. Il est à noter que cette comparaison se limite aux grilles de petite taille, en raison des limitations inhérentes à la simulation de la loi exacte.

## Algorithmes comparés

- Metropolis-Hastings
- Heat Bath
- Swendsen-Wang (Cluster)
- Wolff (Cluster)

## Méthodologie

### Comparaison des algorithmes de type "single flip"

1. **Autocorrélation** : Évaluation basée sur la magnétisation.
2. **Convergence des trajectoires** : Mesure de la vitesse de convergence des différents algorithmes.
3. **Convergence en loi** : Analyse de la convergence de l'erreur relative à la taille de la chaîne de Markov.
4. **Influence de la probabilité d'acceptation** : Étude de son impact sur la convergence.

### Comparaison des algorithmes de cluster

1. **Autocorrélation** : Évaluation de la vitesse à laquelle les algorithmes de cluster oublient leur configuration initiale.
2. **Convergence des trajectoires en temps** : Mesure de la rapidité de convergence des algorithmes vers la loi cible.
3. **Convergence en loi** : Analyse de la convergence pour des valeurs élevées de \( \beta \).

## Résultats

### Comparaison des algorithmes

- **Autocorrélation**:
  - L'autocorrélation de Metropolis décroît plus rapidement vers 0 comparée à celle de Heat Bath.
  - L'algorithme de Metropolis présente une probabilité d'acceptation plus élevée.

![Autocorrélation](https://raw.githubusercontent.com/Nindo16/IsingModelSimulations/main/autocorr_s_b_1.png)

- **Convergence des trajectoires**:
  - Metropolis Hastings converge plus rapidement.

![Convergence des trajectoires](https://raw.githubusercontent.com/Nindo16/IsingModelSimulations/main/convergence_trajectoire_s_b_2.png)

### Comparaison des algorithmes de cluster

- **Autocorrélation**:
  - Les algorithmes de cluster oublient la configuration initiale plus rapidement.

- **Convergence des trajectoires en temps**:
  - L'algorithme de Wolff converge plus vite que celui de Swendsen-Wang.

![Convergence des trajectoires en temps de cluster](https://raw.githubusercontent.com/Nindo16/IsingModelSimulations/main/wf_sw_b_0_8_temps.png)

## Conclusion

Bien que les algorithmes MCMC autorisent la simulation sur de grandes grilles, leur comparaison a été principalement réalisée sur des grilles de petite taille. Cela est dû aux défis associés au calcul de la loi exacte pour des grilles plus vastes. L'algorithme de Wolff s'est distingué comme étant le plus performant parmi les méthodes évaluées.
