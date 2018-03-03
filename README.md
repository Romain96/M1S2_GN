# Projet de géométrie numérique

## Description

Projet réalisé dans le cadre de l'UE "Géométrie Numérique"
* Auteurs: 
    * Romain PERRIN (romain.perrin@etu.unistra.fr)
    * Maxime SEYER (maxime.seyer@etu.unistra.fr)
* Encadrement: Franck Hétroy-Wheeler
* Langage de programmation: C++ (compilé sour QtCreator version 5.9.x)
* Bibliothéques: 
    * glm (fonctions mathématiques)
    * Eigen (outils d'algèbre linéaire)
* Année: 2017 - 2018
* Etablissement: Université de Strasbourg
* Niveau: Master 1
* Semestre: S2
* UE: GN

## Description des classes

| Classe  | Description |
| -------- | -------- |
| DisjointSets | Représente des ensembles disjoints pour la création de MST (Kruskal) |
| Edge | Représente une arête dans le modèle HalfEdge |
| Face | Représente une face dans le modèle HalfEdge |
| Graph | Représente un graphe non orienté |
| HalfEdge | Représente une structure demi-arête dans le modèle HalfEdge |
| Mesh | Représente un ensemble de sommets + ensemble de faces dans le modèle HalfEdge |
| MeshReconstructor | Représente la reconstruction de surface maillée à partir d'un nuage de points |
| Node | Représente un noeud du graphe |
| Octree | Représente l'octree ( structure de données de type arbre dans laquelle chaque nœud compte exactement huit fils) |
| Plane | Représente les plans tangents (représentés par leur trois vecteurs propres) |
| Vertex | Représente les sommets dans le modèle HalfEdge |
