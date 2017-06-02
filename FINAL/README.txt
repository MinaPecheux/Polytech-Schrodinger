--------------
*** README ***
--------------
Projet :	Etude des équations de Schrödinger
Authors :	Anida Khizar, Ramatoulaye Ndiyae, Vincent Nicol, Mina Pêcheux
Date :		Janvier-Juin 2017


0. Contexte, motivation
-----------------------
?


I. Liste des fichiers
---------------------
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


II. Installation & Lancement
----------------------------
Le projet est codé en Python. Il nécessite les packages numpy et matplotlib.
Si nécessaire, ceux-ci peuvent être installés avec : 'pip install numpy matplotlib'.

Chaque script peut être lancé en utilisant la commande terminal : 'python script.py'
(Attention à vérifier que les packages ont été installés pour la bonne version de Python !)