//
//  Vertex.cpp
//  by ziyang @2014 Aug
//
#include "vertex.h"
using namespace std;

Vertex::Vertex()
{
}

/***********Vertex(double, double, double)***********/
Vertex::Vertex(double x, double y, double z)
{
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
	coordinate[3] = 1;

	for (int i = 0; i < 10; i++)
	{
		quadric[i] = 0;
	}
}

Vertex::~Vertex()
{
	vertexList.clear();
}

/***********double* GetCoordinate()****************/
double* Vertex::GetCoordinate()
{
	return coordinate;
}

/***********bool GetVertexList()***************/
list<Vertex*>* Vertex::GetVertexList()
{
	return &vertexList;
}

/***************AddEdge(Edge*)********************/
void Vertex::AddVertex(Vertex* vertex)
{
	int compareResult = 0;

	for (list<Vertex*>::iterator v_iter = vertexList.begin(); v_iter != vertexList.end(); v_iter++)
	{
		compareResult = CompareVertexCoordinate(vertex,*v_iter);
		if (1 == compareResult)
		{
			continue;
		}//end if
		else if (-1 == compareResult)
		{
			vertexList.insert(v_iter, vertex);
			return;
		}//end else if
		else
		{
			return;
		}//end else 
	}//end for
	vertexList.push_back(vertex);
	return;
}

int CompareVertexCoordinate(Vertex* vertex1, Vertex* vertex2)
{
	int result = 0;
	double *coordinate1 = vertex1->GetCoordinate();
	double *coordinate2 = vertex2->GetCoordinate();
	if (coordinate1[2] >coordinate2[2] + ERROR)
	{
		result = 1;
	}
	else if(coordinate1[2] <coordinate2[2] - ERROR)
	{
		result = -1;
	}
	else
	{
		if (coordinate1[0] > coordinate2[0] + ERROR)
		{
			result = 1;
		}
		else if(coordinate1[0] < coordinate2[0] - ERROR)
		{
			result = -1;
		}
		else
		{
			if (coordinate1[1] > coordinate2[1] + ERROR)
			{
				result = 1;
			}
			else if (coordinate1[1] < coordinate2[1] - ERROR)
			{
				result = -1;
			}
			else
			{
				result = 0;
			}
		}
	}
	return result;
}

double* Vertex::GetQuadric()
{
	return quadric;
}