/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

#include <cmath>
#include <iostream>
#include <set>

// GLM
#include "glm/glm.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"

#include "HalfEdge.h"
#include "Vertex.h"

//-----------------------------------------------------------------------------
// Constant(s)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructor(s)
//-----------------------------------------------------------------------------

/**
 * @brief Vertex::Vertex constructs a new Vertex without an incident Half Edge
 * @param id unique id to identify the Vertex
 * @param x coordinate in the x axis
 * @param y coordinate in the y axis
 * @param z coordinate in the z axis
 */
Vertex::Vertex(int id, float x, float y, float z) :
    _id(id),
    _position(x, y, z),
    _incidentHE(nullptr),
    _colour(0.f, 0.f, 0.f, 0.f)
{
    // nothing
}

/**
 * @brief Vertex::Vertex constructs a new Vertex with an incident Half Edge
 * @param id unique id to identify the Vertex
 * @param x coordinate in the x axis
 * @param y coordinate in the y axis
 * @param z coordinate in the z axis
 * @param incidentHE one of the incident half Edges of the Vertex
 */
Vertex::Vertex(int id, float x, float y, float z, HalfEdge *incidentHE) :
    _id(id),
    _position(x, y, z),
    _incidentHE(incidentHE),
    _colour(0.f, 0.f, 0.f, 0.f)
{
    // nothing
}

//-----------------------------------------------------------------------------
// Getter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Vertex::getId
 * @return the unique id of the Vertex
 */
int& Vertex::getId()
{
    return _id;
}

/**
 * @brief Vertex::getPosition
 * @return the coordinates of the Vertex
 */
glm::vec3& Vertex::getPosition()
{
    return _position;
}

/**
 * @brief Vertex::getX
 * @return the coordinate in the x axis
 */
float& Vertex::getX()
{
    return _position.x;
}

/**
 * @brief Vertex::getY
 * @return the coordinate in the y axis
 */
float& Vertex::getY()
{
    return _position.y;
}

/**
 * @brief Vertex::getZ
 * @return the coordinate in the z axis
 */
float& Vertex::getZ()
{
    return _position.z;
}

/**
 * @brief Vertex::getIncidentHE
 * @return a pointer to one of the incident Half Edges in which the Vertex is the origin
 */
HalfEdge *Vertex::getIncidentHE()
{
    return _incidentHE;
}

/**
 * @brief Vertex::getRed
 * @return the red component of the color
 */
float& Vertex::getRed()
{
    return _colour.x;
}

/**
 * @brief Vertex::getGreen
 * @return the green component of the color
 */
float& Vertex::getGreen()
{
    return _colour.y;
}

/**
 * @brief Vertex::getBlue
 * @return the blue component of the color
 */
float& Vertex::getBlue()
{
    return _colour.z;
}

/**
 * @brief Vertex::getTransparency
 * @return the transparency component of the color
 */
float& Vertex::getTransparency()
{
    return _colour.t;
}

//-----------------------------------------------------------------------------
// Setter(s)
//-----------------------------------------------------------------------------

/**
 * @brief Vertex::setId
 * @param id new unique id of tht Vertex
 */
void Vertex::setId(int id)
{
    _id = id;
}

/**
 * @brief Vertex::setPosition
 * @param pos new coordinates of the Vertex
 */
void Vertex::setPosition(glm::vec3 &pos)
{
    _position = pos;
}

/**
 * @brief Vertex::setX
 * @param x new coordinate in the x axis
 */
void Vertex::setX(float &x)
{
    _position.x = x;
}

/**
 * @brief Vertex::setY
 * @param y new coordinate in the y axis
 */
void Vertex::setY(float &y)
{
    _position.y = y;
}

/**
 * @brief Vertex::setZ
 * @param z new coordinate in the z axis
 */
void Vertex::setZ(float &z)
{
    _position.z = z;
}

/**
 * @brief Vertex::setIncidentHE
 * @param incidentHE new pointer on an incident Half Edge in which the Vertex is the origin
 */
void Vertex::setIncidentHE(HalfEdge *incidentHE)
{
    _incidentHE = incidentHE;
}

/**
 * @brief Vertex::setRed
 * @param green new red component of the color (contrained to the interval [0-1])
 */
void Vertex::setRed(float &green)
{
    if (green < 0.f)
        _colour.x = 0.f;
    else if (green > 1.f)
        _colour.x = 1.f;
    else
        _colour.x = green;
}

/**
 * @brief Vertex::setGreen
 * @param green new green component of the color (contrained to the interval [0-1])
 */
void Vertex::setGreen(float &green)
{
    if (green < 0.f)
        _colour.y = 0.f;
    else if (green > 1.f)
        _colour.y = 1.f;
    else
        _colour.y = green;
}

/**
 * @brief Vertex::setBlue
 * @param blue new blue component of the color (contrained to the interval [0-1])
 */
void Vertex::setBlue(float &blue)
{
    if (blue < 0.f)
        _colour.z = 0.f;
    else if (blue > 1.f)
        _colour.z = 1.f;
    else
        _colour.z = blue;
}

/**
 * @brief Vertex::setTransparency
 * @param transparency new transparency component of the color (contrained to the interval [0-1])
 */
void Vertex::setTransparency(float &transparency)
{
    if (transparency < 0.f)
        _colour.t = 0.f;
    else if (transparency > 1.f)
        _colour.t = 1.f;
    else
        _colour.t = transparency;
}

//-----------------------------------------------------------------------------
// Method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Vertex::setColor an all-in-one method (alternative to using setRed, setGreen...)
 * @param red new red component of the color (contrained to the interval [0-1])
 * @param green new green component of the color (contrained to the interval [0-1])
 * @param blue new blue component of the color (contrained to the interval [0-1])
 * @param transparency new transparency component of the color (contrained to the interval [0-1])
 */
void Vertex::setColor(float &red, float &green, float &blue, float &transparency)
{
    setRed(red);
    setGreen(green);
    setBlue(blue);
    setTransparency(transparency);
}

/**
 * @brief Vertex::barycentricArea now using Heron's formula
 * @return the area of the barycentric cell associated with the Vertex
 */
float Vertex::barycentricArea()
{
    // retrieving the incident Half Edge of the Vertex
    HalfEdge *he = _incidentHE;

    // retrieving the previous and next vertices of the Vertex
    Vertex *previous = he->getPrevious()->getOrigin();
    Vertex *next = he->getNext()->getOrigin();

    // computing the center of the triangle
    glm::vec3 center = (_position + previous->getPosition() + next->getPosition()) / 3.f;

    // we now have two triangles
    // triangle 1 : center current_vertex previous_vertex
    // triangle 2 : center current_vertex next_vertex
    // area of the triangle (Heron's formula) : sqrt(s(s-a)(s-b)(s-c))
    // with s the semiperimeter of ABC and a,b,c length of side A,B,C

    // computing the length of each side of each triangle
    float previousSideA = Vertex::distance3(previous->getPosition(), _position);
    float previousSideB = Vertex::distance3(previous->getPosition(), center);
    float previousSideC = Vertex::distance3(center, _position);

    float nextSideA = Vertex::distance3(next->getPosition(), _position);
    float nextSideB = Vertex::distance3(next->getPosition(), center);
    float nextSideC = previousSideC;

    // computing the semi perimeter of each triangle
    float semiperimeterPrevious = (previousSideA + previousSideB + previousSideC) / 2.f;
    float semiperimeterNext = (nextSideA + nextSideB + nextSideC) / 2.f;

    // computing the area of each triangle (using Heron's formula)
    float areaPreviousTriangle = sqrt(semiperimeterPrevious *(semiperimeterPrevious - previousSideA) * (semiperimeterPrevious - previousSideB) * (semiperimeterPrevious - previousSideC));
    float areaNextTriangle = sqrt(semiperimeterNext *(semiperimeterNext - nextSideA) * (semiperimeterNext - nextSideB) * (semiperimeterNext - nextSideC));

    return areaPreviousTriangle + areaNextTriangle;
}

/**
 * @brief Vertex::voronoiArea
 * @return the area of the Voronoï cell associated to the Vertex
 */
float Vertex::voronoiArea()
{
    // retrieving every 1-neighbours of the current Vertex
    std::set<Vertex *> neighbours;
    std::set<Vertex *>::iterator it;
    std::pair<std::set<Vertex *>::iterator, bool> ret;

    // first neighbour
    HalfEdge *he = this->getIncidentHE()->getOpposite();
    Vertex *v = he->getOrigin();
    ret = neighbours.insert(v);

    while (ret.second)
    {
        // v is the origin of the opposite Half Edge
        he = he->getNext()->getOpposite();
        v = he->getOrigin();
        ret = neighbours.insert(v);
    }

    // computing each angle formed between two following neighbours and the current Vertex
    // in the mean time computing the area of each surface
    float angleSum = 0.f;
    float areaSum = 0.f;
    std::set<Vertex *>::iterator vertexIterator1 = neighbours.begin();
    std::set<Vertex *>::iterator vertexIterator2 = neighbours.begin();
    vertexIterator2++;

    float distanceA = 0.f;  // aka distance between current Vertex and first neighbouring Vertex
    float distanceB = 0.f;  // aka distance between current Vertex and second neighbouring Vertex
    float distanceC = 0.f;  // aka distance between the two neighbouring Vertices

    while (vertexIterator2 != neighbours.end())
    {
        float x1 = (*vertexIterator1)->getX();
        float y1 = (*vertexIterator1)->getY();
        float z1 = (*vertexIterator1)->getZ();

        float x2 = (*vertexIterator2)->getX();
        float y2 = (*vertexIterator2)->getY();
        float z2 = (*vertexIterator2)->getZ();

        glm::vec3 neighbour1(x1, y1, z1);
        glm::vec3 neighbour2(x2, y2, z2);

        distanceA = Vertex::distance3(neighbour1, _position);
        distanceB = Vertex::distance3(neighbour2, _position);
        distanceC = Vertex::distance3(neighbour1, neighbour2);

        // computing angle and adding it to the angle sum
        if (distanceA > distanceB)
        {
            // using distance B as adjacent side (A being the hypotenuse)
            angleSum += atan(distanceC / distanceB);
        }
        else
        {
            // using distance A as adjacent side (B being the hypotenuse)
            angleSum += atan(distanceC / distanceA);
        }

        areaSum += Vertex::barycentricArea(this, (*vertexIterator1), (*vertexIterator2));

        vertexIterator1++;
        vertexIterator2++;
    }

    return (2 * M_PI - angleSum) / areaSum;
}

/**
 * @brief Vertex::mixedArea
 * @return the area of the mixed Voronoï cell associated to the Vertex
 */
float Vertex::mixedArea()
{
    // TODO
    return 0.f;
}

/**
 * @brief Vertex::gaussianCurvature
 * @return the Gaussian curvature in the current Vertex
 */
float Vertex::gaussianCurvature()
{
    // retrieving every 1-neighbours of the current Vertex
    std::set<Vertex *> neighbours;
    std::pair<std::set<Vertex *>::iterator, bool> ret;

    // first neighbour
    HalfEdge *he = this->getIncidentHE()->getOpposite();
    Vertex *v = he->getOrigin();
    ret = neighbours.insert(v);

    while (ret.second)
    {
        // v is the origin of the opposite Half Edge
        he = he->getNext()->getOpposite();
        v = he->getOrigin();
        ret = neighbours.insert(v);
    }

    // computing each angle formed between two following neighbours and the current Vertex
    // in the mean time computing the area of each surface
    float angleSum = 0.f;
    float areaSum = 0.f;
    std::set<Vertex *>::iterator vertexIterator1 = neighbours.begin();
    std::set<Vertex *>::iterator vertexIterator2 = neighbours.begin();
    vertexIterator2++;

    float distanceA = 0.f;  // aka distance between current Vertex and first neighbouring Vertex
    float distanceB = 0.f;  // aka distance between current Vertex and second neighbouring Vertex
    float distanceC = 0.f;  // aka distance between the two neighbouring Vertices

    while (vertexIterator2 != neighbours.end())
    {
        float x1 = (*vertexIterator1)->getX();
        float y1 = (*vertexIterator1)->getY();
        float z1 = (*vertexIterator1)->getZ();

        float x2 = (*vertexIterator2)->getX();
        float y2 = (*vertexIterator2)->getY();
        float z2 = (*vertexIterator2)->getZ();

        glm::vec3 neighbour1(x1, y1, z1);
        glm::vec3 neighbour2(x2, y2, z2);

        distanceA = Vertex::distance3(neighbour1, _position);
        distanceB = Vertex::distance3(neighbour2, _position);
        distanceC = Vertex::distance3(neighbour1, neighbour2);

        // computing angle and adding it to the angle sum
        if (distanceA > distanceB)
        {
            // using distance B as adjacent side (A being the hypotenuse)
            angleSum += atan(distanceC / distanceB);
        }
        else
        {
            // using distance A as adjacent side (B being the hypotenuse)
            angleSum += atan(distanceC / distanceA);
        }

        areaSum += Vertex::barycentricArea(this, (*vertexIterator1), (*vertexIterator2));

        vertexIterator1++;
        vertexIterator2++;
    }

    return (2 * M_PI - angleSum) / areaSum;
}

/**
 * @brief Vertex::meanCurvature
 * @return the mean curvature in the current Vertex
 */
float Vertex::meanCurvature()
{
    // TODO
    return 0.f;
}

/**
 * @brief Vertex::uniformLaplacian
 * @return a vector being the Laplacian operator with uniform weights
 */
glm::vec3& Vertex::uniformLaplacian()
{
    // retrieving every 1-neighbours of the current Vertex
    std::set<Vertex *> neighbours;
    std::pair<std::set<Vertex *>::iterator, bool> ret;

    // first neighbour
    HalfEdge *he = this->getIncidentHE()->getOpposite();
    Vertex *v = he->getOrigin();
    ret = neighbours.insert(v);

    // until we try to insert an already processed Half Edge
    while (ret.second)
    {
        // v is the origin of the opposite Half Edge
        he = he->getNext()->getOpposite();
        v = he->getOrigin();
        ret = neighbours.insert(v);
    }

    // computing the Laplacian operator
    static glm::vec3 laplacian(0.f, 0.f, 0.f);

    std::set<Vertex *>::iterator neighbourIterator = neighbours.begin();

    while(neighbourIterator != neighbours.end())
    {
        laplacian = laplacian + ((*neighbourIterator)->getPosition() - _position);

        neighbourIterator++;
    }

    laplacian /= (float)(neighbours.size());

    return laplacian;
}

/**
 * @brief Vertex::cotangentLaplacian
 * @return a vector being the Laplacian operator with cotangent weights
 */
glm::vec3& Vertex::cotangentLaplacian()
{
    // retrieving every 1-neighbours of the current Vertex
    std::set<Vertex *> neighbours;
    std::pair<std::set<Vertex *>::iterator, bool> ret;
    std::set<Vertex *>::iterator tempIterator;
    float cotanPrevious;
    float cotanNext;

    // first neighbour
    HalfEdge *he = this->getIncidentHE()->getOpposite();
    Vertex *v = he->getOrigin();
    ret = neighbours.insert(v);

    // until we try to insert an already processed Half Edge
    while (ret.second)
    {
        // v is the origin of the opposite Half Edge
        he = he->getNext()->getOpposite();
        v = he->getOrigin();
        ret = neighbours.insert(v);
    }

    // computing the Laplacian operator
    static glm::vec3 laplacian(0.f, 0.f, 0.f);
    float laplacianFactor = 0.f;

    std::set<Vertex *>::iterator neighbourIterator = neighbours.begin();
    Vertex *previous = nullptr;
    Vertex *next = nullptr;

    while(neighbourIterator != neighbours.end())
    {
        // retrieving the previous and next Vertex
        if (neighbourIterator == neighbours.begin())
        {
            // special case : first boundary
            tempIterator = neighbours.end();
            tempIterator--;
            previous = (*tempIterator);
            tempIterator = neighbours.begin();
            tempIterator++;
            next = (*tempIterator);
        }
        else if (neighbourIterator == --neighbours.end())
        {
            // special case : last boundary
            tempIterator = neighbours.end();
            tempIterator--;
            previous = (*tempIterator);
            next = (*neighbours.begin());
        }
        else
        {
            tempIterator = neighbourIterator;
            tempIterator--;
            previous = (*tempIterator);
            tempIterator++;
            tempIterator++;
            next = (*tempIterator);
        }

        // using vectors, scalar and cross product
        glm::vec3 v1 = previous->getPosition() - _position;
        glm::vec3 v2 = (*neighbourIterator)->getPosition() - previous->getPosition();
        glm::vec3 v3 = next->getPosition() - _position;
        glm::vec3 v4 = (*neighbourIterator)->getPosition() - next->getPosition();

        // computing the two cotangents
        cotanPrevious = glm::dot(v1, v2) / Vertex::distance3(v1, v2);
        cotanNext = glm::dot(v3, v4) / Vertex::distance3(v3, v4);

        // adding to the sum
        laplacianFactor += (cotanPrevious + cotanNext);
        laplacian = laplacian + (cotanPrevious + cotanNext) * ((*neighbourIterator)->getPosition() - _position);

        neighbourIterator++;
    }

    laplacian = laplacian / laplacianFactor;

    return laplacian;
}

//-----------------------------------------------------------------------------
// Private method(s)
//-----------------------------------------------------------------------------

/**
 * @brief Vertex::barycentricArea private constrained version
 * @param v1 first Vertex of the triangle (used as reference point)
 * @param v2 second Vertex of the triangle
 * @param v3 third Vertex of the triangle
 * @return the area of the barycentric cell associated with the Vertex
 */
float Vertex::barycentricArea(Vertex *v1, Vertex *v2, Vertex *v3)
{
    // computing the center of the triangle
    float centerX = (v1->getX() + v2->getX() + v3->getX())/3.f;
    float centerY = (v1->getY() + v2->getY() + v3->getY())/3.f;
    float centerZ = (v1->getZ() + v2->getZ() + v3->getZ())/3.f;

    // we now have two triangles
    // triangle 1 : center current_vertex previous_vertex
    // triangle 2 : center current_vertex next_vertex
    // area of the triangle (Heron's formula) : sqrt(s(s-a)(s-b)(s-c))
    // with s the semiperimeter of ABC and a,b,c length of side A,B,C

    // computing the length of each side of each triangle
    float previousSideA = sqrt(pow((v2->getX() - v1->getX()), 2.f) + pow((v2->getY() - v1->getY()), 2.f) + pow((v2->getZ() - v1->getZ()), 2.f));
    float previousSideB = sqrt(pow((v2->getX() - centerX), 2.f) + pow((v2->getY() - centerY), 2.f) + pow((v2->getZ() - centerZ), 2.f));
    float previousSideC = sqrt(pow((centerX - v1->getX()), 2.f) + pow((centerY - v1->getY()), 2.f) + pow((centerZ - v1->getZ()), 2.f));

    float nextSideA = sqrt(pow((v3->getX() - v1->getX()), 2.f) + pow((v3->getY() - v1->getY()), 2.f) + pow((v3->getZ() - v1->getZ()), 2.f));
    float nextSideB = sqrt(pow((v3->getX() - centerX), 2.f) + pow((v3->getY() - centerY), 2.f) + pow((v3->getZ() - centerZ), 2.f));
    float nextSideC = previousSideC;

    // computing the semi perimeter of each triangle
    float semiperimeterPrevious = (previousSideA + previousSideB + previousSideC) / 2.f;
    float semiperimeterNext = (nextSideA + nextSideB + nextSideC) / 2.f;

    // computing the area of each triangle (using Heron's formula)
    float areaPreviousTriangle = sqrt(semiperimeterPrevious *(semiperimeterPrevious - previousSideA) * (semiperimeterPrevious - previousSideB) * (semiperimeterPrevious - previousSideC));
    float areaNextTriangle = sqrt(semiperimeterNext *(semiperimeterNext - nextSideA) * (semiperimeterNext - nextSideB) * (semiperimeterNext - nextSideC));

    return areaPreviousTriangle + areaNextTriangle;
}

/**
 * @brief Vertex::distance3
 * @param v1 first Vertex coordinates
 * @param v2 second Vertex coordinates
 * @return the distance (in dimension 3) between v1 and v2
 */
float Vertex::distance3(glm::vec3 &v1, glm::vec3 &v2)
{
    return sqrt(pow((v1.x - v2.x), 2.f) + pow((v1.y - v2.y), 2.f) + pow((v1.z - v2.z), 2.f));
}
