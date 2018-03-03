#ifndef VERTEX_H
#define VERTEX_H

/*
 * (C) Romain PERRIN
 * romain.perrin@etu.unistra.fr
 * UFR de Mathématiques-Informatique
 * 2018-2019
 */

// GLM
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"

#include "HalfEdge.h"

// forward declaration
class HalfEdge;

/**
 * @brief The Vertex class
 */
class Vertex
{
protected:
    // unique ID
    int _id;

    // coordinates XYZ
    glm::vec3 _position;

    // incident Half Edge
    HalfEdge *_incidentHE;  // in which it is the origin

    // colors RGBA
    glm::vec4 _colour;

public:
    // constructor
    Vertex();
    Vertex(int id, float x, float y, float z);
    Vertex(int id, float x, float y, float z, HalfEdge *incidentHE);

    // getters
    int& getId();

    glm::vec3& getPosition();
    float& getX();
    float& getY();
    float& getZ();

    HalfEdge *getIncidentHE();

    float& getRed();
    float& getGreen();
    float& getBlue();
    float& getTransparency();

    // setters
    void setId(int id);

    void setPosition(glm::vec3& pos);
    void setX(float& x);
    void setY(float& y);
    void setZ(float& z);

    void setIncidentHE(HalfEdge *incidentHE);

    void setRed(float& green);
    void setGreen(float& green);
    void setBlue(float& blue);
    void setTransparency(float& transparency);

    // methods
    void setColor(float& red, float& green, float& blue, float& transparency);

    float barycentricArea();
    float voronoiArea();
    float mixedArea();

    float gaussianCurvature();
    float meanCurvature();

    glm::vec3& uniformLaplacian();
    glm::vec3& cotangentLaplacian();

private:
    // methods
    static float barycentricArea(Vertex *v1, Vertex *v2, Vertex *v3);

public:
    static float distance3(glm::vec3& v1, glm::vec3& v2);
};

#endif // VERTEX_H
