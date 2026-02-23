# Projet_Informatique_Scientifique
L3 Info-maths, parcours de graphes, Julia
Dépendances utilisées : DataStructures
Depuis le répertoire racine du projet git :
>julia --project=.
>using Algo
Il est possible de lancer les commandes suivantes :
	algoBFS(fname, D, A)
	algoDijkstra(fname, D, A)
	algoGlouton(fname, D, A)
	algoAstar(fname, D, A)
avec:
-fname le nom du fichier de la map (certaines sont disponibles dans ./Algo/data)
-D le point de départ du parcours 
-A le point d'arrivée du parcours

Une execution d'un de ces quatres algorithmes pour une instance présentant D ou A sur un obstacle infranchissable retourna directement un message indiquant l'impossibilité de traverser cet obstacle

Une execution  d'un de ces quatres algorithmes pour une instance présentant D ou A sur une "île" se déroulera jusqu'au bout pouvant conduire à des performances dégradées, l'entièreté du graphe connexe des sommets accesibles depuis D sera visité, les algorithmes se terminent cependant.
	
	
Note: La version de Dijkstra et d'A* ont été implémentés selon une variante présente dans le livre Informatique MP2I de Jean-Christophe Filiâtre

