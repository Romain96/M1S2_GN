#ifndef MESH_H
#define MESH_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <iostream>
#include <vector>
#include <map>

#include "Vertex.h"
#include "Face.h"
#include "HalfEdge.h"

/**
 * @brief The Mesh class
 */
class Mesh
{
protected:
    // list of all Vertices
    std::vector<Vertex *> _vertices;
    // list of all Faces
    std::vector<Face *> _faces;
    // map of all Half Edges indexed by their origin Vertex
    std::map<int, std::vector<HalfEdge *>> _halfEdges;

    // Laplacian smoothing parameters
    float _laplacianLambda;
    unsigned int _laplacianIteration;

    // Taubin smoothing parameter
    float _taubinMu;

public:
    // constructor
    Mesh();

    // getters
    Vertex *getVertex(unsigned int i);
    Face *getFace(unsigned int i);

    float& getLaplacianLambda();
    unsigned int& getLaplacianIteration();
    float& getTaubinMu();

    // setters
    void setLaplacianLambda(float lambda);
    void setLaplacianIteration(unsigned int k);
    void setTaubinMu(float mu);

    // methods
    void importOFF(std::string filename);
    void exportOFF(std::string filename);
    void exportWithColors(std::string filename);

    void applyGaussianCurvature();
    void applyMeanCurvature();

    void laplacianSmoothing();
    void taubinSmoothing();
    void taubinSmoothingCotangent();

private:
    // internal methods
    void constructHalfEdgeMesh();
    void constructOppositeRelation();

};

#endif // MESH_H
