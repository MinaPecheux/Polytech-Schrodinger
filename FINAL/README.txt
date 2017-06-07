--------------
*** README ***
--------------
Projet :	Etude des équations de Schrödinger
Authors :	Anida Khizar, Ramatoulaye Ndiaye, Vincent Nicol, Mina Pêcheux
Date :		Janvier-Juin 2017


I. Contexte, motivation
-----------------------
Le but de ce projet était d'avoir une première approche de la physique quantique en étudiant les équations de Schrödinger et en les implémentant sous Python à l'aide d'outils d'analyse numérique. Les différentes problématiques ont été :
	- d'appréhender et de manipuler correctement les outils de base de la physique quantique, principalement la fonction d'onde et la probabilité de présence d'une particule
	- d'utiliser des méthodes numériques vues en cours pour passer du monde continu à un monde discret et ainsi trouver des solutions approchées des équations de Schrödinger
	- de découvrir les bibliothèques Python numpy (pour l'algèbre linéaire) et matplotlib (pour la représentation graphique)
Après avoir mis en place un cadre d'étude valide pour les équations de Schrödinger, nous nous en sommes servi pour illustrer des effets caractéristiques des fonctions d'onde. Nous avons également axé une partie de la réflexion sur les points communs et les différences entre mécanique quantique et mécanique classique.


II. Liste des fichiers
----------------------
README.txt							Ce fichier

animation-scripts/				Dossier contenant les fichiers étudiant les effets de réflexion et l'effet tunnel
	reflection_animation.py 		Implémentation de l'effet de réflexion et sortie d'un graphe animé (enregistré ou live)
	tunnel_animation.py 			Implémentation de l'effet tunnel et sortie d'un graphe animé (enregistré ou live)
basic-scripts/					Dossier contenant les fichiers essentiels pour la modélisation de la fonction d'onde en 1D et 2D
								indépendant du temps, et 1D dépendant du temps
	1D_indep.py 					Implémentation d'un paquet d'onde à partir de l'équation de Schrödinger indépendante du temps
									en 1D avec plusieurs potentiels
	2D_indep.py 					Implémentation d'un paquet d'onde à partir de l'équation de Schrödinger indépendante du temps
									en 2D avec le potentiel nul (cas simple)
	1D_dep.py 						Implémentation d'un paquet d'onde à partir de l'équation de Schrödinger dépendante du temps
									en 1D à l'aide de plusieurs schémas numériques (Euler explicite, Euler implicite, Cranck-Nicholson)
movies/							Dossier contenant des enregistrements du début de l'évolution temporelle d'un paquet d'onde gaussien
								(pour l'étude des effets de réflexion et de l'effet tunnel)
	reflection_animation.mp4 		Enregistrement des effets de réflexion (400 images)
	tunnel_animation.mp4 			Enregistrement de l'effet tunnel (400 images)


III. Installation & Lancement
-----------------------------
Le projet est codé en Python. Il nécessite les packages numpy et matplotlib.
Si nécessaire, ceux-ci peuvent être installés avec : 'pip install numpy matplotlib'.

Chaque script peut être lancé en utilisant la commande terminal : 'python script.py'
(Attention à vérifier que les packages ont été installés pour la bonne version de Python !)

IV. Pistes de réflexion
-----------------------
Ce projet nous a permis d'avoir une première idée des principes et des questions posés pour le monde quantique, cependant ce n'est qu'un premier pas ! De nombreuses améliorations et réflexions nouvelles sont possibles à partir de là, par exemple l'implémentation des équations de Schrödinger en 3D et des interactions entre plusieurs particules, l'optimisation de notre code et la mise en place de schémas numériques plus performant, l'étude de phénomènes quantiques plus complexes...
