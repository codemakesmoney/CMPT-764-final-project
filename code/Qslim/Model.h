//
//  Model.h
//  Qslim
//
//  Created by ziyang on 14-8-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#ifndef __finalproject__Model__
#define __finalproject__Model__

#include <iostream>
#include <list>
#include <math.h>
#include <time.h>

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
	Model(SUModelRef);//construct a model use SKP model
    ~Model();
    
    void SaveModel();//save the model to a smf file   
    void DrawModel();//draw the model
    void DecimateModel(int);//decimate the model
    
private:
	SUModelRef modelSU;
    
    double xMiddle, yMiddle, zMiddle;//middle point's coordinate
	double zMin, zMax;
    GLdouble range;//model maximum range

	list<Face*> faceList;
	list<Edge*> edgeList;
	list<Vertex*> vertexList;

	void ComputeModelPosition();
    /*User the max/min x/y/z coordinate to compute the model's middle point and range*/
    void ComputeMiddleAndRange(GLdouble, GLdouble, GLdouble, GLdouble, GLdouble, GLdouble);

	void ConvertSKPtoModel();//convert the SKP model to our data structure

	void AddFace(SUFaceRef);//add a face to the face list
	Edge* AddEdge(SUEdgeRef);//add a edge to the edge list

	void FindEdgeToDecimate();//decimate one edge
	double ComputeVertexAndError(Edge*, double*);//compute the error if decimate one edge, and the new vertex position
	Vertex* GetVertex(double x, double y, double z);//get a vertex from the vertex list, if this vertex does not exist, create it
	Edge* GetEdge(Vertex*, Vertex*);//get a edge from the edge list, if this edge does not exist, create it
};
bool ComputeCoordinate(double*, double*);
double ComputeError(double coor[4], double Q[][4]);


#endif /* defined(__finalproject__Model__) */