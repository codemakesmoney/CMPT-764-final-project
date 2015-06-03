//
//  Edge.h
//  by ziyang @2014 Aug
//
#ifndef __finalproject__Edge__
#define __finalproject__Edge__

#include"Face.h"
#include"Vertex.h"
using namespace std;

class Face;

class Edge
{
public:
	Edge();//construct a edge
	Edge(Vertex *, Vertex *);//construct a edge by the start vertex and end vertex
	~Edge();

	Vertex* GetStartVertex();//get the start vertex
	Vertex* GetEndVertex();//get the end vertex
	void SetStartVertex(Vertex*);//set the start vertex
	void SetEndVertex(Vertex*);//set the end vertex
	void AddFace(Face*);//add a face to the face list
	list<Face*>* GetFaceList();//return the all the faces that is associated with this edge

private:
	Vertex *startVertex;
	Vertex *endVertex;
	list<Face*> faceList;//all the faces that is associated with this edge
};
#endif /* defined(__finalproject__Edge__) */