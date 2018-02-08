/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Vertex.h"
#include "Face.h"
#include "Mesh.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Mesh::Mesh constructs a new empty Mesh
 */
Mesh::Mesh() :
    _vertices(),
    _faces(),
    _halfEdges(),
    _laplacianLambda(0.5),
    _laplacianIteration(1),
    _taubinMu(0.52)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Mesh::getVertices
 * @return the vertices
 */
std::vector<Vertex *>& Mesh::getVertices()
{
    return _vertices;
}

/**
 * @brief Mesh::getVertex
 * @param i index of the Vertex to retrieve
 * @return a pointer to the desired Vertex (or null if outside of range)
 */
Vertex *Mesh::getVertex(unsigned int i)
{
    if (i > _vertices.size())
        return nullptr;
    return _vertices[i];
}

/**
 * @brief Mesh::getFace
 * @param i index of the Face to retrieve
 * @return a pointer to the desired Face (or null is outside of range)
 */
Face *Mesh::getFace(unsigned int i)
{
    if (i > _faces.size())
        return nullptr;
    return _faces[i];
}

/**
 * @brief Mesh::getLaplacianLambda
 * @return the lambda coefficient used in the laplacianSmoothing method
 */
float& Mesh::getLaplacianLambda()
{
    return _laplacianLambda;
}

/**
 * @brief Mesh::getLaplacianIteration
 * @return the number of iterations used in the laplacianSmoothing method
 */
unsigned int& Mesh::getLaplacianIteration()
{
    return _laplacianIteration;
}

/**
 * @brief Mesh::getTaubinMu
 * @return the mu coefficient used in the taubinSmoothing method
 */
float& Mesh::getTaubinMu()
{
    return _taubinMu;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Mesh::setLaplacianLambda
 * @param lambda new lambda coefficient used in the laplacianSmoothing method
 */
void Mesh::setLaplacianLambda(float lambda)
{
    _laplacianLambda = lambda;
}

/**
 * @brief Mesh::setLaplacianIteration
 * @param k new number of iterations used in the laplacianSmoothing method
 */
void Mesh::setLaplacianIteration(unsigned int k)
{
    _laplacianIteration = k;
}

/**
 * @brief Mesh::setTaubinMu
 * @param mu new mu coefficient used in the taubinSmoothing method
 */
void Mesh::setTaubinMu(float mu)
{
    _taubinMu = mu;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Mesh::importOFF loads a file in the OFF format
 * @param filename path where to load the file
 */
void Mesh::importOFF(std::string filename)
{
    int vertices = 0;
    int faces = 0;
    int vertices_processed = 0;
    //int faces_processed = 0;

    //int vertices_per_face = 0;
    float x, y, z;
    //int vertexIndex = 0;

    Vertex *vertex = nullptr;
    //Face *face = nullptr;

    std::string line;

    // opening file
    std::ifstream input(filename, std::ifstream::in);

    std::cout << "opening file named " << filename << std::endl;

    // line 1 : keyword "OFF" mandatory
    std::getline(input, line);

    if (line.compare("OFF") != 0)
    {
        std::cerr << "cannot open file" << std::endl;
        exit(1);
    }

    // line 2 : x y 0 (x number of vertices, y number of faces)
    std::getline(input, line);
    std::stringstream ss(line);
    ss >> vertices >> faces;
    std::cout << "Mesh contains " << vertices << " vertices and  " << faces << " faces." << std::endl;

    // following x lines : vertices in format x y z
    while (vertices_processed < vertices)
    {
        std::getline(input, line);
        std::stringstream ss(line);
        ss >> x >> y >> z;

        //std::cout << "vertex " << vertices_processed << " : x = " << x << ", y = " << y << ", z = " << z << std::endl;

        // storing vertices in std::vector (incident Half Edge is not yet known)
        vertex = new Vertex(vertices_processed, x, y, z);
        _vertices.push_back(vertex);

        vertices_processed++;
    }

    std::cout << vertices_processed << " vertices processed." << std::endl;

    /*// following y lines : faces index in format i1 i2 i3 (triangles)
    while (faces_processed < faces)
    {
        std::getline(input, line);
        std::stringstream ss(line);
        ss >> vertices_per_face;

        // creating the face
        face = new Face(faces_processed);

        // processing vertices index one by one as they depends on the number "vertices_per_face"
        for (int i = 0; i < vertices_per_face; i++)
        {
            ss >> vertexIndex;

            // retrieving vertex pointer corresponding to this index
            vertex = _vertices[vertexIndex];

            // adding the the vertices list of the face
            face->addVertex(vertex);
        }
        _faces.push_back(face);

        faces_processed++;
    }

    std::cout << faces_processed << " faces processed." << std::endl;

    // WARNING : experimental
    std::cout << "building Half Edge relations..." << std::endl;
    constructHalfEdgeMesh();
    std::cout << "done !" << std::endl;*/

    // closing file
    input.close();
}

/**
 * @brief Mesh::exportOFF saves the current Half Edge structure to a file using the OFF format
 * @param filename path where to write the file
 */
void Mesh::exportOFF(std::string filename)
{
    std::cout << "saving Mesh in format OFF in file " << filename << std::endl;

    // opening file
    std::ofstream output(filename);

    // writting the "OFF" header in first line
    output << "OFF\n";

    // writting the "x y 0" header in second line where x is the vertices number and y the faces number
    std::cout << "Mesh contains " << _vertices.size() << " vertices and " << _faces.size() << " faces." << std::endl;
    output << _vertices.size() << " " << _faces.size() << " 0\n";

    // writting x lines containing "x y z" coordinates of vertices
    std::vector<Vertex *>::iterator vertexIterator = _vertices.begin();

    while (vertexIterator != _vertices.end())
    {
        output << (*vertexIterator)->getX() << " " << (*vertexIterator)->getY() << " " << (*vertexIterator)->getZ() << "\n";
        vertexIterator++;
    }

    // writting y lines containing "i i1 i2 ... in" where i is the number of vertices per face
    // and i1 ... in the indices of vertices defining the polygonal faces
    std::vector<Face *>::iterator faceIterator = _faces.begin();

    while (faceIterator != _faces.end())
    {
        output << (*faceIterator)->getVerticesNumber();

        // writting vertices indices (number depends on the number of vertices per face)
        for (int i = 0; i < (*faceIterator)->getVerticesNumber(); i++)
        {
            output << " " << (*faceIterator)->getVertex(i)->getId();
        }
        output << "\n";

        faceIterator++;
    }

    // closing file
    output.close();
}

/**
 * @brief Mesh::exportWithColors
 * @param filename path where to write the file
 */
void Mesh::exportWithColors(std::string filename)
{
    std::cout << "saving Mesh with colors in format COFF in file " << filename << std::endl;

    // opening file
    std::ofstream output(filename);

    // writting the "COFF" header in first line
    output << "COFF\n";

    // writting the "x y 0" header in second line where x is the vertices number and y the faces number
    std::cout << "Mesh contains " << _vertices.size() << " vertices and " << _faces.size() << " faces." << std::endl;
    output << _vertices.size() << " " << _faces.size() << " 0\n";

    // writting x lines containing "x y z" coordinates of vertices
    std::vector<Vertex *>::iterator vertexIterator = _vertices.begin();

    while (vertexIterator != _vertices.end())
    {
        output << (*vertexIterator)->getX() << " " << (*vertexIterator)->getY() << " " << (*vertexIterator)->getZ();
        // adding colors (red [0-1] green [0-1] blue[0-1] transparency [0-1]
        output << " " << (*vertexIterator)->getRed() << " " << (*vertexIterator)->getGreen() << " " << (*vertexIterator)->getBlue()
               << " " << (*vertexIterator)->getTransparency() << "\n";
        vertexIterator++;
    }

    // writting y lines containing "i i1 i2 ... in" where i is the number of vertices per face
    // and i1 ... in the indices of vertices defining the polygonal faces
    std::vector<Face *>::iterator faceIterator = _faces.begin();

    while (faceIterator != _faces.end())
    {
        output << (*faceIterator)->getVerticesNumber();

        // writting vertices indices (number depends on the number of vertices per face)
        for (int i = 0; i < (*faceIterator)->getVerticesNumber(); i++)
        {
            output << " " << (*faceIterator)->getVertex(i)->getId();
        }
        output << "\n";

        faceIterator++;
    }

    // closing file
    output.close();
}

/**
 * @brief Mesh::applyGaussianCurvature applies a Gaussian curvature and gives a color to each Vertex
 */
void Mesh::applyGaussianCurvature()
{
    std::cout << "applying Gaussian curvature computation with colours" << std::endl;

    // applying the method "Vertex::gaussianCurvature" on each Vertex
    std::vector<Vertex *>::iterator vertexIterator = _vertices.begin();

    // keeping track of the minimal and maximal values of curvature to normalize the colors after procerssing
    float minCurvature = 0.f;
    float maxCurvature = 0.f;
    float curvature = 0.f;
    // simple vector to store temporary curvature values
    std::vector<float> curvatures;

    while (vertexIterator != _vertices.end())
    {
        curvature = (*vertexIterator)->gaussianCurvature();
        curvatures.push_back(curvature);

        if (curvature < minCurvature)
            minCurvature = curvature;
        if (curvature > maxCurvature)
            maxCurvature = curvature;

        vertexIterator++;
    }

    // now normalizing curvatures values and giving each Vertex a colour based on its curvature
    // normalizing curvatures to be between 0 and 1
    float factor = 1.f / (maxCurvature - minCurvature);
    std::vector<float>::iterator curvatureIterator = curvatures.begin();
    float red = 0.f;
    float green = 0.f;
    float blue = 0.f;
    float transparency = 0.f;
    float normalizedCurvature = 0.f;
    vertexIterator = _vertices.begin();

    while (curvatureIterator != curvatures.end())
    {
        // going from 0,0,0 to 1,1,1
        // or from black to white through blue, green and red
        normalizedCurvature = factor * ((*curvatureIterator) - minCurvature);   // [min...max] -> [0...1]

        if (normalizedCurvature <= 1.f/3.f)
        {
            red = 0.f;
            green = 0.f;
            blue = 3 * normalizedCurvature;
        }
        else if (normalizedCurvature > 1.f/3.f && normalizedCurvature <= 2.f/3.f)
        {
            red = 0.f;
            green = 3 * (normalizedCurvature - 1.f/3.f);
            blue = 1.f;
        }
        else if (normalizedCurvature > 2.f/3.f && normalizedCurvature <= 1.f)
        {
            red = 3 * (normalizedCurvature - 2.f/3.f);
            green = 1.f;
            blue = 1.f;
        }
        else
        {
            // errors : normalizeCurvature < 0 or > 1
            std::cerr << "normalized curvature outside range [0...1] : " << normalizedCurvature << std::endl;
        }

        (*vertexIterator)->setColor(red, green, blue, transparency);

        curvatureIterator++;
        vertexIterator++;
    }

    std::cout << "Gaussian curvature applied successfully" << std::endl;
}

/**
 * @brief Mesh::applyMeanCurvature
 */
void Mesh::applyMeanCurvature()
{
}

/**
 * @brief Mesh::laplacianSmoothing smooths k time the Mesh using the formula xi <- xi + lambda * L(xi)
 * where L is the method uniformLaplacian in class Vertex
 * lambda can be modified using the setLaplacianLambda() method in Mesh class
 * k can be modified using the setLaplacianIteration() method in Mesh class
 */
void Mesh::laplacianSmoothing()
{
    std::cout << "applying a Laplacian smoothing with parameters k = " << _laplacianIteration <<
                 " and lambda = " << _laplacianLambda << std::endl;

    // ATTENTION : in order to keep current values of Vertices until one
    // iteration is finished, we need to maintain a vector of the newest values
    // and instanciating the vertices after each iteration is complete

    std::vector<glm::vec3> temp;
    std::vector<Vertex *>::iterator vertexIterator;
    std::vector<glm::vec3>::iterator tempIterator;
    int vertexNum;

    // applying the smoothing algorithm k times
    for (unsigned int i = 0; i < _laplacianIteration; i++)
    {
        // iterating for each Vertex in _vertices
        vertexIterator = _vertices.begin();

        while (vertexIterator != _vertices.end())
        {
            // computing the Laplacian operator and storing the new position of the Vertex in temp
            glm::vec3 laplacian = (*vertexIterator)->uniformLaplacian();
            glm::vec3 newPosition = (*vertexIterator)->getPosition() + _laplacianLambda * laplacian;
            temp.push_back(newPosition);

            vertexIterator++;
        }


        // now replacing old positions with new positions stored in temp
        tempIterator = temp.begin();
        vertexNum = 0;

        while(tempIterator != temp.end())
        {
            _vertices[vertexNum]->setPosition((*tempIterator));

            vertexNum++;
            tempIterator++;
        }

        // don't forget to clear the temp values after updating
        temp.clear();
    }

    std::cout << "Laplacian smoothing applied successfully !" << std::endl;
}

/**
 * @brief Mesh::taubinSmoothing smooths k time the Mesh using alternatively
 * the formula xi <- xi + lambda * L(xi) and xi <- xi - mu * L(xi)
 * where L is the method uniformLaplacian in class Vertex
 * lambda can be modified using the setLaplacianLambda() method in Mesh class
 * k can be modified using the setLaplacianIteration() method in Mesh class
 * mu can be modified using the setTaubinMu() method in Mesh class
 */
void Mesh::taubinSmoothing()
{
    std::cout << "applying a Taubin smoothing with parameters k = " << _laplacianIteration <<
                 " , mu = " << _taubinMu << " and lambda = " << _laplacianLambda << std::endl;

    // ATTENTION : in order to keep current values of Vertices until one
    // iteration is finished, we need to maintain a vector of the newest values
    // and instanciating the vertices after each iteration is complete

    std::vector<glm::vec3> temp;
    std::vector<Vertex *>::iterator vertexIterator;
    std::vector<glm::vec3>::iterator tempIterator;
    int vertexNum;

    // same pricinciple as laplacianSmoothing but we use this to elternate the formula
    bool parity = true;

    // applying the smoothing algorithm k times
    for (unsigned int i = 0; i < _laplacianIteration; i++)
    {
        // iterating for each Vertex in _vertices
        vertexIterator = _vertices.begin();

        while (vertexIterator != _vertices.end())
        {
            // computing the Laplacian operator and storing the new position of the Vertex in temp
            glm::vec3 laplacian = (*vertexIterator)->uniformLaplacian();
            glm::vec3 newPosition;

            // the formula now depends on parity
            if (parity)
            {
                // Laplacian smoothing when true
                newPosition = (*vertexIterator)->getPosition() + _laplacianLambda * laplacian;
            }
            else
            {
                // and Taubin smoothing when false
                newPosition = (*vertexIterator)->getPosition() - _taubinMu * laplacian;
            }

            temp.push_back(newPosition);

            vertexIterator++;
        }


        // now replacing old positions with new positions stored in temp
        tempIterator = temp.begin();
        vertexNum = 0;

        while(tempIterator != temp.end())
        {
            _vertices[vertexNum]->setPosition((*tempIterator));

            vertexNum++;
            tempIterator++;
        }

        // don't forget to clear the temp values after updating
        temp.clear();
        // and to change the parity (otherwise it's just s regular Laplacian smoothing)
        parity = !parity;
    }

    std::cout << "Taubin smoothing applied sucessfully !" << std::endl;
}

/**
 * @brief Mesh::taubinSmoothingCotangent smooths k time the Mesh using alternatively
 * the formula xi <- xi + lambda * L(xi) and xi <- xi - mu * L(xi)
 * where L is the method uniformLaplacian in class Vertex
 * lambda can be modified using the setLaplacianLambda() method in Mesh class
 * k can be modified using the setLaplacianIteration() method in Mesh class
 * mu can be modified using the setTaubinMu() method in Mesh class
 * Uses cotangent weights when computing Laplacian operator
 * contrary to taubinSmoothing method which uses uniform weights
 */
void Mesh::taubinSmoothingCotangent()
{
    std::cout << "applying a Taubin smoothing (cotangent) with parameters k = " << _laplacianIteration <<
                 " , mu = " << _taubinMu << " and lambda = " << _laplacianLambda << std::endl;

    // ATTENTION : in order to keep current values of Vertices until one
    // iteration is finished, we need to maintain a vector of the newest values
    // and instanciating the vertices after each iteration is complete

    std::vector<glm::vec3> temp;
    std::vector<Vertex *>::iterator vertexIterator;
    std::vector<glm::vec3>::iterator tempIterator;
    int vertexNum;

    // same pricinciple as laplacianSmoothing but we use this to elternate the formula
    bool parity = true;

    // applying the smoothing algorithm k times
    for (unsigned int i = 0; i < _laplacianIteration; i++)
    {
        // iterating for each Vertex in _vertices
        vertexIterator = _vertices.begin();

        while (vertexIterator != _vertices.end())
        {
            // computing the Laplacian operator and storing the new position of the Vertex in temp
            glm::vec3 laplacian = (*vertexIterator)->cotangentLaplacian();
            glm::vec3 newPosition;

            // the formula now depends on parity
            if (parity)
            {
                // Laplacian smoothing when true
                newPosition = (*vertexIterator)->getPosition() + _laplacianLambda * laplacian;
            }
            else
            {
                // and Taubin smoothing when false
                newPosition = (*vertexIterator)->getPosition() - _taubinMu * laplacian;
            }

            temp.push_back(newPosition);

            vertexIterator++;
        }


        // now replacing old positions with new positions stored in temp
        tempIterator = temp.begin();
        vertexNum = 0;

        while(tempIterator != temp.end())
        {
            _vertices[vertexNum]->setPosition((*tempIterator));

            vertexNum++;
            tempIterator++;
        }

        // don't forget to clear the temp values after updating
        temp.clear();
        // and to change the parity (otherwise it's just s regular Laplacian smoothing)
        parity = !parity;
    }

    std::cout << "Taubin (cotangent) smoothing applied sucessfully !" << std::endl;
}

//-----------------------------------------------------------------------------
// Private method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Mesh::constructHalfEdgeMesh creates all Half Edges
 */
void Mesh::constructHalfEdgeMesh()
{
    HalfEdge *he = nullptr;
    std::vector<Face *>::iterator faceIterator = _faces.begin();

    // iterating on each Face
    while (faceIterator != _faces.end())
    {
        // creating Half Edges
        for (int i = 0; i < (*faceIterator)->getVerticesNumber() - 1; i++)
        {
            // creating Half Edge from current Vertex(i) to next Vertex(i+1)
            he = new HalfEdge((*faceIterator)->getVertex(i), (*faceIterator)->getVertex(i + 1));
            // he belongs to the current face
            he->setFace((*faceIterator));
            // he is added in class Vertex
            (*faceIterator)->getVertex(i)->setIncidentHE(he);
            // adding the Half Edge to the list
            _halfEdges[(*faceIterator)->getVertex(i)->getId()].push_back(he);
        }

        // adding the last Half Edge of the current Face (from Vertex(last) to Vertex(first))
        he = new HalfEdge((*faceIterator)->getVertex((*faceIterator)->getVerticesNumber() - 1), (*faceIterator)->getVertex(0));
        // he belongs to the current face
        he->setFace((*faceIterator));
        // let's use it to be the representative in class Face
        (*faceIterator)->setIncidentHE(he);
        // he is added in class Vertex
        (*faceIterator)->getVertex((*faceIterator)->getVerticesNumber() - 1)->setIncidentHE(he);
        // adding the Half Edge to the list
        _halfEdges[(*faceIterator)->getVertex((*faceIterator)->getVerticesNumber() - 1)->getId()].push_back(he);

        // now building the previous/next relation between created Half Edges
        for (int i = 1; i < (*faceIterator)->getVerticesNumber() - 1; i++)
        {
            he = (*faceIterator)->getVertex(i)->getIncidentHE();

            // previous is incidentHE of vertex(i - 1)
            he->setPrevious((*faceIterator)->getVertex(i - 1)->getIncidentHE());

            // next is incidentHE of vertex(i + 1)
            he->setNext((*faceIterator)->getVertex(i + 1)->getIncidentHE());
        }

        // linking first to last
        he = (*faceIterator)->getVertex(0)->getIncidentHE();
        // previous is incidentHE of vertex(last)
        he->setPrevious((*faceIterator)->getVertex((*faceIterator)->getVerticesNumber() - 1)->getIncidentHE());
        // next is incidentHE of vertex(1)
        he->setNext((*faceIterator)->getVertex(1)->getIncidentHE());

        // linking last to first
        he = (*faceIterator)->getVertex((*faceIterator)->getVerticesNumber() - 1)->getIncidentHE();
        // previous is incidentHE of vertex(last - 1)
        he->setPrevious((*faceIterator)->getVertex((*faceIterator)->getVerticesNumber() - 2)->getIncidentHE());
        // next is incidentHE of vertex(0)
        he->setNext((*faceIterator)->getVertex(0)->getIncidentHE());

        faceIterator++;
    }

    // building all opposite relations between Half Edges
    constructOppositeRelation();
}

/**
 * @brief Mesh::constructOppositeRelation constructs the opposite relations between the created Half Edges
 */
void Mesh::constructOppositeRelation()
{
    // algorithm principle :
    // FOR EACH Vertex v IN _halfEdges DO
    // ---- FOR EACH HalfEdge he IN _halfEdge[v] DO
    // ---- ---- FOR EACH Vertex v' IN _halfEdges AND HalfEdge he' WITH origin = he.destination AND destination = he.origin IN _halfEdges[v']
    // ---- ---- ---- SET _halfEdges[v'].he'.opposite = _halfEdges[v].he
    // ---- ---- ---- SET _halfEdges[v].he.opposite = _halfEdges[v'].he'
    // ---- ---- END FOR
    // ---- END FOR
    // END FOR

    std::map<int, std::vector<HalfEdge *>>::iterator vertexIterator = _halfEdges.begin();
    std::vector<HalfEdge *>::iterator halfEdgeIterator;
    std::vector<HalfEdge *>::iterator oppositeHalfEdgeIterator;
    std::vector<HalfEdge *> currentHalfEdges;
    std::vector<HalfEdge *> oppositeHalfEdges;
    Vertex *targetOrigin = nullptr;
    Vertex *targetDestination = nullptr;

    // iterating on each Vertex
    while (vertexIterator != _halfEdges.end())
    {
        currentHalfEdges = (*vertexIterator).second;
        halfEdgeIterator = currentHalfEdges.begin();

        // iterating on every Half Edge of whose origin is the current Vertex
        while (halfEdgeIterator != currentHalfEdges.end())
        {
            targetOrigin = (*halfEdgeIterator)->getDestination();
            targetDestination = (*halfEdgeIterator)->getOrigin();

            // finding opposites (aka every Vertex with origin = myDestination && destination = myOrigin)
            oppositeHalfEdges = _halfEdges[targetOrigin->getId()];
            oppositeHalfEdgeIterator = oppositeHalfEdges.begin();
            while (oppositeHalfEdgeIterator != oppositeHalfEdges.end())
            {
                if ((*oppositeHalfEdgeIterator)->getDestination()->getId() == targetDestination->getId())
                {
                    (*halfEdgeIterator)->setOpposite((*oppositeHalfEdgeIterator));
                    (*oppositeHalfEdgeIterator)->setOpposite((*halfEdgeIterator));
                }

                oppositeHalfEdgeIterator++;
            }

            halfEdgeIterator++;
        }

        vertexIterator++;
    }
}
