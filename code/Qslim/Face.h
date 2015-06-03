//
//  Face.h
//  by ziyang @2014 Aug
//
#ifndef __finalproject__Face__
#define __finalproject__Face__

#include<list>
#include "Edge.h"
using namespace std;

class Edge;
class Vertex;
class Face
{
public:
	Face();//construct a face
	~Face();
	void AddEdge(Edge*);//add an edge that form this face
	list<Edge*>* GetEdgeList();//get all the edges that form this face
	double* GetPanel();//return the plane of this face
	void ComputePanelAndQuadric();//compute the plane and quadric of this face

private:
	list<Edge*> edgeList;//all the edges that form this face
	double panel[4];//plane of this face
};
#endif