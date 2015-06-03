//
//  Edge.cpp
//  by ziyang @2014 Aug
//
#include "Edge.h"
using namespace std;

Edge::Edge()
{
}

/***************Edge(Vertex*, Vertex*)*************/
Edge::Edge(Vertex * v1, Vertex * v2)
{
	startVertex = v1;
	endVertex = v2;
	//v1 is connected with v2, and v2 is connected with v1 
	v1->AddVertex(v2);
	v2->AddVertex(v1);
}

Edge::~Edge()
{
	faceList.clear();
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

/************void Edge::SetStartVertex(Vertex* sv)*******/
void Edge::SetStartVertex(Vertex* sv)
{
	startVertex = sv;
}

/***********void Edge::SetEndVertex(Vertex* ev)*****************/
void Edge::SetEndVertex(Vertex* ev)
{
	endVertex = ev;
}

/************void Edge::AddFace(Face* face)*********************/
void Edge::AddFace(Face* face)
{
	for (list<Face*>::iterator f_iter = faceList.begin(); f_iter != faceList.end(); f_iter++)
	{
		if (face == *f_iter)
		{
			return;
		}
	}
	faceList.push_back(face);
}

/*************list<Face*>* Edge::GetFaceList()**********/
list<Face*>* Edge::GetFaceList()
{
	return &faceList;
}