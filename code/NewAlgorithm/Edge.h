//
// by ziyang @2014 Aug
//  Edge.h
//
#ifndef __finalproject__Edge__
#define __finalproject__Edge__

#include"Vertex.h"
using namespace std;

class Edge
{
public:
	Edge();/*construct an edge*/
	Edge(Vertex *, Vertex *);/*construct an edge by two vertex*/
	~Edge();

	Vertex* GetStartVertex();/*get the start vertex*/
	Vertex* GetEndVertex();/*get the end vertex*/
	void SetStartVertex(Vertex*);/*set the start vertex*/
	void SetEndVertex(Vertex*);/*set the end vertex*/

private:
	Vertex *startVertex;
	Vertex *endVertex;
};

bool MergeEdgeList(Edge*, Edge*);/*The compare function to merge two edge list*/
#endif /* defined(__finalproject__Edge__) */