/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <iostream>
#include <vector>

#include "Vertex.h"
#include "Face.h"
#include "HalfEdge.h"
#include "Mesh.h"
#include "Octree.h"
#include "MeshReconstructor.h"

using namespace std;

int main()
{
    Mesh m;

    // loading the point cloud (loading only vertices of mesh)
    m.importOFF("../OFF/block.off");

    // MeshReconstructor test
    MeshReconstructor mr(m.getVertices());
    // parameters : k = 10, Octree built with 2 iterations (64 leaves)
    mr.setK(10);
    mr.setIterations(2);
    // parameters : 0-dense, 0.05-noisy

    mr.setDense(10000.f);
    mr.setNoisy(0.5f);
    mr.buildPointTreeWithIterations();
    mr.computeCentroidsAndTangentPlanes();
    mr.buildCentroidTreeWithIterations();
    /*
    mr.reorientateTangentPlanes();
    */
    // marching cube parameters : isolevel = 0, cell size 0.5
    mr.setIsolevel(0.f);
    mr.setSubdivisionFactor(1.f);
    mr.createIsosurface();

    // reading test
    //m.importOFF("../OFF/block.off");

    // Gaussian curvature + writting in colors
    //m.applyGaussianCurvature();
    //m.exportWithColors("../OFF/block_color.off");

    // Laplacian smoothing test (base parameters are k = 1, lambda = 0.5
    //m.setLaplacianIteration(25);
    //m.laplacianSmoothing();
    //m.exportWithColors("../OFF/block_laplacian.off");

    // Taubin smoothing test (base parameters are k = 1, lambda = 0.5, mu = 0.52)
    //m.taubinSmoothing();
    //m.exportWithColors("../OFF/block_taubin.off");

    // Taubin smoothing with cotangent weights test
    //m.taubinSmoothingCotangent();
    //m.exportWithColors("../OFF/block_taubin_cotangent.off");

    // NOTE : fast way to visualize results : http://masc.cs.gmu.edu/wiki/ObjViewer

    return 0;
}
