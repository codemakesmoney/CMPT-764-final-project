//
//  Model.h
//  finalproject
//
//  Created by Ziyang Zhao on 14-8-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#ifndef __finalproject__Model__
#define __finalproject__Model__

#include <iostream>
#include <list>
#include <math.h>

#include <slapi/slapi.h>
#include <slapi/geometry.h>
#include <slapi/initialize.h>
#include <slapi/unicodestring.h>
#include <slapi/model/model.h>
#include <slapi/model/entities.h>
#include <slapi/model/face.h>
#include <slapi/model/edge.h>
#include <slapi/model/vertex.h>
#include <vector>

#if defined(GLUI_FREEGLUT)
#include <GL/freeglut.h>

#elif defined(GLUI_OPENGLUT)
#include <GL/openglut.h>

#else

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <GLUI/glui.h>// Header File For The GLUI Library
#include <OpenGL/gl.h> // Header File For The OpenGL Library
#include <OpenGL/glu.h> // Header File For The GLu Library
#else
#include <GL/glut.h>
#include <GL/glui.h>// Header File For The GLUI Library
#include <GL/gl.h> // Header File For The OpenGL Library
#include <GL/glu.h> // Header File For The GLu Library
#endif

#endif

#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
using namespace std;


class Model
{
public:
	Model(SUModelRef);
    ~Model();
    
    void SaveModel();//save the model to a smf file   
    void DrawModel();//draw the model
    void DecimateModel(double, double);//decimate the model
    
private:
	SUModelRef modelSU;
    
    double xMiddle, yMiddle, zMiddle;//middle point's coordinate
	double zMin, zMax;//minimum and maximum height of the model
    GLdouble range;//model maximum range

	list<Face*> faceList;
	list<Edge*> edgeList;
	list<Vertex*> vertexList;

	/*compute the minimun and maximum x,y,z of the model*/
	void ComputeModelPosition();
    /*User the max/min x/y/z coordinate to compute the model's middle point and range*/
    void ComputeMiddleAndRange(GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble);
	/*convert the model from skp to our data structure*/
	void ConvertSKPtoModel();

	/*Add a face to the face list*/
	void AddFace(SUFaceRef);
	/*Add a edge to the edge list*/
	Edge* AddEdge(SUEdgeRef);

	/*Get the profile of a height*/
	void GetSection(double, list<Edge*>*, list<Vertex*>*);
	/*Get the intersection of a face at a height*/
	void GetIntersectionOfFace(Face*, double, list<Edge*>*, list<Vertex*>*);
	/*Update the current profile according to the compare profile*/
	void UpdateNewModel(list<Vertex*>*, list<Vertex*>*, list<Vertex*>*, list<Edge*>*, list<Edge*>*, list<Edge*>*);
};

/*get a vertex according to the coordinate, if the vertex is not in the list, create it*/
Vertex* GetVertex(double x, double y, double z, list<Vertex*>*);
/*get a edge according to the start vertex and end vertex, if the edge is not in the list, create it*/
Edge* GetEdge(Vertex*, Vertex*, list<Edge*>*);
/*match the vertex from two profile*/
double GetCorrPoint(list<Vertex*>*, list<Vertex*>*, list<Edge*>*, list<Edge*>*);
/*given two coordinates, compute the 2D distance*/
double Get2DDist(double*, double*);

void AddNewEdge(list<Vertex*>*, list<Edge*>*);
void AddCylinder(list<Edge*>*, list<Vertex*>*, list<Edge*>*, list<Face*>*);
void UpdateCurrentSection(list<Edge*>*, list<Vertex*>*);
void SetSectionToDefault(list<Edge*>*, list<Vertex*>*);

#endif /* defined(__finalproject__Model__) */