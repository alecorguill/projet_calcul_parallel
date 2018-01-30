# Projet calcul parallèle
## Compilation
   
	make   

## Exécution
Le programme s'appuie sur un fichier de configuration nommé config.cfg. <br>
Il faut rentrer les valeurs désirées avant d'exécuter le programme.

     mpirun -np Np -- ./run config.cfg

## Visualisation et tests
**Tests unitaires** : <br>

     make test
     ./test_unitaire

**Courbe de speedup** : <br>
&nbsp;&nbsp; Le script  script.sh permet de générer un fichier csv que l'on peut visualiser avec le programme plot_csv.py.
Le makefile encapsule ces étapes comme ceci :

    make courbe


