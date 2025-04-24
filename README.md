![Illustration du plan D-optimal](Graphiques/DONEcal.jpg)

# Plan d'expériences numériques pour la calibration de codes de calcul coûteux (à sortie scalaire)

Ce dépôt contient des stratégies de planification sequentielle d'experiences numériques pour la calibration d'un code de calculs coûteux à sortie scalaire. Le cadre statistique utilisé est celui de [Kennedy et O'Hagan (2001)](https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/kennedy01.pdf). 
------- en cours
## Description


## Installation 

Clonez ce dépôt pour télécharger les fichiers en local :

```bash
git clone https://github.com/TheseAdama/DONEcal.git
```
Vous pouvez également télécharger directement le fichier ZIP depuis GitHub.

## Package R : 
Exécutez le code suivant dans R pour installer les packages nécessaires : 

 ```r
install.packages(c("Dicekriging", "SimDesign", "lhs", "truncnorm"))
 ```

## Reference
**Adama Barry, François Bachoc, Sarah Bouquet, Miguel Munoz Zuniga, Clémentine Prieur. _Optimal Design of Physical and Numerical Experiments for Computer Code Calibration_. 2024. [hal-04615127v2](https://theses.hal.science/UNIV-UT3/hal-04615127v2)**

