//
// by ziyang @2014 Aug
//  Edge.cpp
//
#include "Edge.h"
using namespace std;

Edge::Edge()
{
}

/***************Edge(Vertex* v1, Vertex* v2)*************/
/****v1 is the start vertex , v2 is the end vertex*********/
Edge::Edge(Vertex * v1, Vertex * v2)
{
	startVertex = v1;
	endVertex = v2;

	/*v1 is connect to v2, and v2 is connect to v1*/
	v1->AddVertex(v2);
	v2->AddVertex(v1);
}

Edge::~Edge()
{
	/*delete start vertex from end vertex's connect vertex list*/
	/*delete end vertex from start vertex's connect vertex list*/
	startVertex->GetVertexList()->remove(endVertex);
	endVertex->GetVertexList()->remove(startVertex);
}

/****************Vertex* GetStartVertex()****************/
Vertex* Edge::GetStartVertex()
{
	return startVertex;
}

/****************Vertex* GetEndVertex()*****************/
Vertex* Edge::GetEndVertex()
{
	return endVertex;
}

/***************void SetStartVertex(Vertex*)*************/
void Edge::SetStartVertex(Vertex* sv)
{
	startVertex = sv;
}

/***************void SetEndVertex(Vertex*)**************/
void Edge::SetEndVertex(Vertex* ev)
{
	endVertex = ev;
}


/****************bool MergeEdgeList(Edge*, Edge*)********/
bool MergeEdgeList(Edge* e1, Edge* e2)
{
	Vertex * startVertex1 = e1->GetStartVertex();
	Vertex *startVertex2 = e2->GetStartVertex();
	Vertex *endVertex1 = e1->GetEndVertex();
	Vertex *endVertex2 = e2->GetEndVertex();
	
	int result1 = CompareVertexCoordinate(startVertex1, startVertex2); 
	/*if the start vertex is smaller, the edge is smaller*/
	if(-1 == result1)
	{
		return true;
	}
	/*if the start vertex is bigger, the edge is bigger*/
	else if (0 == result1)
	{
		/*if the start vertex are the same, if the end vertex is smaller, the vertex is smaller*/
		if (-1 == CompareVertexCoordinate(endVertex1, endVertex2))
		{
			return true;
		}
		/*if the start vertex are the same, if the end vertex is not smaller, the vertex is bigger*/
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}