//
// by ziyang @2014 Aug
//  Vertex.cpp
//
#include "vertex.h"
using namespace std;

/***********Vertex(double, double, double)***********/
Vertex::Vertex(double x, double y, double z)
{
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;

	corrVertex = NULL;
	distance = BIGNUM;
}

/*****************Vertex(double*)*********************/
Vertex::Vertex(double* coor)
{
	coordinate[0] = coor[0];
	coordinate[1] = coor[1];
	coordinate[2] = coor[2];

	corrVertex = NULL;
	distance = BIGNUM;
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

/***********double GetDistance()*******************/
double Vertex::GetDistance()
{
	return distance;
}

/***********void SetDistance(double)****************/
void Vertex::SetDistance(double d)
{
	distance = d;
}

/***********void SetCorrPoint(Vertex*)*************/
void Vertex::SetCorrPoint(Vertex* v)
{
	corrVertex = v;
}

/***********bool GetVertexList()***************/
list<Vertex*>* Vertex::GetVertexList()
{
	return &vertexList;
}

/***********bool GetHasCorrPoint()***************/
Vertex* Vertex::GetCorrPoint()
{
	return corrVertex;
}

/***************AddVertex(Vertex* vertex)********************/
void Vertex::AddVertex(Vertex* vertex)
{
	int compareResult = 0;

	for (list<Vertex*>::iterator v_iter = vertexList.begin(); v_iter != vertexList.end(); v_iter++)
	{
		//compare vertex to the vertexs in the list
		compareResult = CompareVertexCoordinate(vertex,*v_iter);
		if (1 == compareResult)//if the vertex to add is bigger
		{
			continue;
		}//end if
		else if (-1 == compareResult)//if the vertex to add is smaller
		{
			vertexList.insert(v_iter, vertex);//insert the vertex here
			return;
		}//end else if
		else//if this vertex is the vertex to insert
		{
			return;
		}//end else 
	}//end for
	vertexList.push_back(vertex);
	return;
}

/*******void SetCoordinate(double x, double y, double z)******/
void Vertex::SetCoordinate(double x, double y, double z)
{
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
}

/*CompareVertexCoordinate(Vertex* vertex1, Vertex* vertex2)*/
int CompareVertexCoordinate(Vertex* vertex1, Vertex* vertex2)
{
	int result = 0;
	double *coordinate1 = vertex1->GetCoordinate();
	double *coordinate2 = vertex2->GetCoordinate();

	/*if z is smaller, the vertex is smaller*/
	if (coordinate1[2] >coordinate2[2] + ERROR)
	{
		result = 1;
	}
	/*if z is bigger, the vertex is bigger*/
	else if(coordinate1[2] <coordinate2[2] - ERROR)
	{
		result = -1;
	}
	else
	{
		/*if z are the same, and x is smaller, the vertex is smaller*/
		if (coordinate1[0] > coordinate2[0] + ERROR)
		{
			result = 1;
		}
		/*if z are the same, and x is bigger, the vertex is bigger*/
		else if(coordinate1[0] < coordinate2[0] - ERROR)
		{
			result = -1;
		}
		else
		{
			/*if z, x are the same, and y is smaller, the vertex is smaller*/
			if (coordinate1[1] > coordinate2[1] + ERROR)
			{
				result = 1;
			}
			/*if z, x are the same, and y is bigger, the vertex is bigger*/
			else if (coordinate1[1] < coordinate2[1] - ERROR)
			{
				result = -1;
			}
			/*if z, x and y are all the same, then the two vertex are the same*/
			else
			{
				result = 0;
			}
		}
	}
	return result;
}

/****MergeVertexList(Vertex* v1, Vertex* v2)***/
bool MergeVertexList(Vertex* v1, Vertex* v2)
{
	/*if v1 < v2, return true, else return false*/
	return (-1 == CompareVertexCoordinate(v1, v2));
}