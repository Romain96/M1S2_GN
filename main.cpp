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

using namespace std;

int main()
{
    Mesh m;

    // testing on a hand made cube
    m.importOFF("../OFF/test.off");

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
