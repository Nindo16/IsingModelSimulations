# Simulation du modèle d'Ising avec des algorithmes MCMC

## Introduction

Ce dépôt concerne la simulation du modèle d'Ising, un modèle mathématique en physique statistique. En utilisant des techniques d'énumération et de calcul de l'énergie associées aux configurations binaires, nous explorons efficacement l'ensemble des configurations possibles. La simulation exacte, bien que précise, est limitée aux grilles de petite taille.

## Objectif du projet

L'objectif principal de ce travail était de comparer différents algorithmes de Monte Carlo par chaîne de Markov (MCMC) pour la simulation du modèle d'Ising.

## Algorithmes comparés

- Metropolis-Hastings
- Heat Bath
- Swendsen-Wang (Cluster)
- Wolff (Cluster)

## Méthodologie

### Comparaison des algorithmes de type "single flip"

1. **Autocorrélation** : Utilisé la magnétisation pour évaluer l'autocorrélation.
2. **Convergence des trajectoires** : Évalué la vitesse de convergence des algorithmes.
3. **Convergence en loi** : Observé la convergence de l'erreur en fonction de la taille de la chaîne de Markov.
4. **Influence de la probabilité d'acceptation** : Examined l'impact de la probabilité d'acceptation sur la convergence.

### Comparaison des algorithmes de cluster

1. **Autocorrélation** : Les algorithmes de cluster ont été évalués pour déterminer la rapidité avec laquelle ils oublient la configuration initiale.
2. **Convergence des trajectoires en temps** : La vitesse à laquelle les algorithmes convergent vers la loi cible a été mesurée.
3. **Convergence en loi** : Étudié la convergence en loi pour des valeurs de \( \beta \) élevées.

## Résultats clés

### Comparaison des algorithmes

- **Autocorrélation**:
  - Metropolis est généralement plus rapide que Heat Bath.
  - Probabilité d'acceptation plus élevée avec Metropolis comparé à Heat Bath.

![Insérez votre image d'autocorrélation ici](lien_vers_votre_image_autocorrélation)

- **Convergence des trajectoires**:
  - Metropolis Hastings converge généralement plus rapidement.

![Insérez votre image de convergence des trajectoires ici](lien_vers_votre_image_convergence_trajectoires)

- **Convergence en loi**:
  - Faible erreur pour des valeurs de \( \beta \) élevées.

![Insérez votre image de convergence en loi ici](lien_vers_votre_image_convergence_en_loi)

### Comparaison des algorithmes de cluster

- **Autocorrélation**:
  - Les algorithmes de cluster oublient rapidement la configuration initiale.

![Insérez votre image d'autocorrélation de cluster ici](lien_vers_votre_image_autocorrélation_cluster)

- **Convergence des trajectoires en temps**:
  - Wolff converge plus rapidement que Swendsen Wang.

![Insérez votre image de convergence des trajectoires en temps de cluster ici](lien_vers_votre_image_convergence_trajectoires_temps_cluster)

- **Convergence en loi**:
  - Les deux algorithmes montrent une faible erreur pour des valeurs élevées de \( \beta \).

![Insérez votre image de convergence en loi de cluster ici](lien_vers_votre_image_convergence_en_loi_cluster)

## Conclusion

Les algorithmes de Monte Carlo par chaîne de Markov permettent la simulation sur de grandes grilles. L'algorithme de Wolff semble être le plus performant parmi les méthodes testées.

