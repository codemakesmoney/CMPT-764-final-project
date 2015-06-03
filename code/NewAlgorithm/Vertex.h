//
// by ziyang @2014 Aug
//  Vertex.h
//
#ifndef __finalproject__Vertex__
#define __finalproject__Vertex__

#include<stdio.h>
#include<list>
using namespace std;

#define BIGNUM 2000000000
#define ERROR 0.00001

class Vertex
{
public:
	Vertex(double x, double y, double z);//construct a vertex by coordiates
	Vertex(double* coor);//construct a vertex by coordiates vector
	~Vertex();

	double *GetCoordinate();//get the coordiate vector of this vertex
	double GetDistance();//get the distance from this vertex to its associated vertex
	Vertex* GetCorrPoint();//get the associated vertex of this vertex
	list<Vertex*>* GetVertexList();//get the vertexs that are connected to this vertex

	void SetDistance(double);//set the distance from this vertex to its associated vertex
	void SetCorrPoint(Vertex*);//set the associated vertex of this vertex 
	void SetCoordinate(double x, double y, double z);//set the coordiates
	void AddVertex(Vertex *);//add a vertex to its connected vertex

private:
	double coordinate[3];//coordiate vector
	Vertex* corrVertex;//the associated vertex
	double distance;//the distance to its associated vertex
	list<Vertex*> vertexList;//the vertex list of its connected vertex
};

/***********bool CompareVertexCoordinate(Vertex*, Vertex*)*******/
/***********return 0:v1 == v2
								1: v1 > v2
								-1:v1
								v2**************************/
int CompareVertexCoordinate(Vertex*, Vertex*);

bool MergeVertexList(Vertex*, Vertex*);//the compare function to merge two vertex list together

#endif /* defined(__finalproject__Vertex__) */
