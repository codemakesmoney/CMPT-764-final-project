//
// by ziyang @2014 Aug
//  Face.h
//
#ifndef __finalproject__Face__
#define __finalproject__Face__

#include<list>
#include "Edge.h"
using namespace std;

class Face
{
public:
	Face();/*construct a face*/
	~Face();
	void AddEdge(Edge*);/*Add a edge that form this face*/
	list<Edge*>* GetEdgeList();/*Get all the edges that form this face*/
	double GetZMin();/*Get the lowest height of this face*/
	double GetZMax();/*Get the highest height of this face*/

private:
	list<Edge*> edgeList;/*all the edges that form this face*/
	double zMin, zMax;/* lowest height and highest height*/
};

bool MergeFaceList(Face*, Face*);/*The compare function to merge two face list*/
#endif