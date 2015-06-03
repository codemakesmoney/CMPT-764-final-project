//
// by ziyang @2014 Aug
//  Face.cpp
//
#include "Face.h"
using namespace std;

Face::Face()
{
	zMin = BIGNUM;
	zMax = 0 - BIGNUM;
}

/**********~Face()*******************************/
Face::~Face()
{
	edgeList.clear();
}

/***************AddEdge(Edge*)********************/
void Face::AddEdge(Edge * edge)
{
	int compareResult = 0;
	Vertex* startVertex = edge->GetStartVertex();
	Vertex* endVertex = edge->GetEndVertex();
	double zMinCurr = startVertex->GetCoordinate()[2];
	double zMaxCurr = endVertex->GetCoordinate()[2];

	/*update zMin and zMax*/
	if (zMinCurr < zMin)
	{
		zMin = zMinCurr;
	}
	if (zMaxCurr > zMax)
	{
		zMax = zMaxCurr;
	}

	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		compareResult = CompareVertexCoordinate(startVertex, (*e_iter)->GetStartVertex());
		/*if the new edge is bigger then the current edge*/
		if (1 == compareResult)
		{
			continue;
		}//end if
		/*if the new edge is smaller than the current edge*/
		else if (-1 == compareResult)
		{
			edgeList.insert(e_iter, edge);
			return;
		}//end else if
		else
		{
			compareResult = CompareVertexCoordinate(endVertex, (*e_iter)->GetEndVertex());
			/*if the new edge is smaller than the current edge*/
			if (1 == compareResult)
			{
				continue;
			}//end else-if
			/*if the new edge is bigger than the current edge*/
			else if (-1 == compareResult)
			{
				edgeList.insert(e_iter, edge);
				return;
			}//end else-else if
			else
			{
				return;
			}//end else-else
		}//end else 
	}//end for
	edgeList.push_back(edge);
	return;
}

/*****************list<Edge*>* GetEdgeList()************/
list<Edge*>* Face::GetEdgeList()
{
	return &edgeList;
}

/*****************double GetZMin()*******************/
double Face::GetZMin()
{
	return zMin;
}

/*****************double GetZMax()*******************/
double Face::GetZMax()
{
	return zMax;
}

/*******bool MergeFaceList(Face*, Face*)*****************/
bool MergeFaceList(Face* f1, Face*f2)
{
	double zMin1, zMin2, zMax1, zMax2;
	zMin1 = f1->GetZMin();
	zMin2 = f2->GetZMin();
	zMax1 = f1->GetZMax();
	zMax2 = f2->GetZMax();

	/*if zMin is smaller, this face is smaller*/
	if (zMin1 < zMin2 - ERROR)
	{
		return true;
	}
	/*if zMin is bigger, the face is bigger*/
	else if (zMin1 > zMin2 + ERROR)
	{
		return false;
	}
	else
	{
		/*if zMin are the same, zMax is smaller, the face is smaller*/
		if (zMax1 < zMax2 - ERROR)
		{
			return true;
		}
		/*if zMin are the same,zMax is not smaller, the face is bigger */
		else
		{
			return false;
		}
	}
}