//
//  Vertex.h
//  by ziyang @2014 Aug
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
	Vertex();//construct vertex
	Vertex(double x, double y, double z);//construct vertex by coordinate
	~Vertex();

	double *GetCoordinate();//get the coordinate vector of the vertex
	double *GetQuadric();//get the quadric of the vertex
	list<Vertex*>* GetVertexList();//get the vertexs that are connected with this vertex
	void AddVertex(Vertex *);//add a vertex to the vertex list


private:
	double coordinate[4];//the coordinate vector of the vertex
	double quadric[10];//the quadric of the vertex
	list<Vertex*> vertexList;//all the vertex that are connected with this vertex
};

//compare two vertexs, the vertex with a smaller z, x, y is smaller
int CompareVertexCoordinate(Vertex* vertex1, Vertex* vertex2);


#endif /* defined(__finalproject__Vertex__) */
