Introduction
============

Contexte
--------

SPHirit est un simulateur de mécanique des fluides 2D développé dans le cadre 
d'un projet académique au CESI, en collaboration avec des enseignants-chercheurs.

Il s'inscrit dans une démarche de développement et de validation de modèles 
numériques appliqués à des problématiques académiques et industrielles.

La méthode SPH
--------------

La méthode **SPH (Smoothed Particle Hydrodynamics)** est une méthode particulaire 
sans maillage. Au lieu de discrétiser l'espace en cellules fixes, le fluide est 
représenté par un ensemble de particules qui se déplacent librement.

Chaque particule porte des propriétés physiques (masse, vitesse, pression) et 
interagit avec ses voisines via une fonction de lissage appelée **kernel**.

La formulation **ALE (Arbitrary Lagrangian-Eulerian)** utilisée ici permet de 
combiner les avantages des approches Lagrangienne et Eulérienne pour améliorer 
la précision des simulations.

Références
----------

- Violeau, D. (2012). *Fluid Mechanics and the SPH Method*. Oxford University Press.
- Liu, M.B. & Liu, G.R. (2010). Smoothed Particle Hydrodynamics: an Overview and 
  Recent Developments. *Archives of Computational Methods in Engineering*, 17, 25-76.

Projets connexes
----------------

- `SPHinXsys <https://www.sphinxsys.org>`_ — bibliothèque SPH professionnelle en C++
- `OpenFOAM <https://www.openfoam.com>`_ — référence en simulation de fluides open-source