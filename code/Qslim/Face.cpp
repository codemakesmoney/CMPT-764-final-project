//
//  Face.cpp
//  by ziyang @2014 Aug
//
#include "Face.h"
using namespace std;

/*************Face()***************/
Face::Face()
{
	for (int i = 0; i < 4; i++)
	{
		panel[i] = 0;
	}
}

/**********~Face()*******************************/
Face::~Face()
{
	/*remove this face from all the edge's face list that associated with this face*/
	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		(*e_iter)->GetFaceList()->remove(this);
	}
	edgeList.clear();
}

/***************AddEdge(Edge*)********************/
void Face::AddEdge(Edge * edge)
{
	int compareResult = 0;
	Vertex* startVertex = edge->GetStartVertex();
	Vertex* endVertex = edge->GetEndVertex();

	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		compareResult = CompareVertexCoordinate(startVertex, (*e_iter)->GetStartVertex());
		if (1 == compareResult)
		{
			continue;
		}//end if
		else if (-1 == compareResult)
		{
			edgeList.insert(e_iter, edge);
			return;
		}//end else if
		else
		{
			compareResult = CompareVertexCoordinate(endVertex, (*e_iter)->GetEndVertex());
			if (1 == compareResult)
			{
				continue;
			}//end else-if
			else if (-1 == compareResult)
			{
				edgeList.insert(e_iter, edge);
				return;
			}//end else-else if
			else
			{
				return;
			}
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

double* Face::GetPanel()
{
	return panel;
}

void Face::ComputePanelAndQuadric()
{
	Edge* edge = NULL;
	Vertex * v1 = NULL, *v2 = NULL, *v3 = NULL;
	double quadric[10];
	double *vQuadric1 = NULL, *vQuadric2 = NULL;

	/*get three vertex of this face*/
	edge = edgeList.front();
	v1 = edge->GetStartVertex();
	v2 = edge->GetEndVertex();

	edge = edgeList.back();
	v3 = edge->GetStartVertex();
	if (CompareVertexCoordinate(v1, v3) == 1 || CompareVertexCoordinate(v2, v3) == 1)
	{
		v3 = edge->GetEndVertex();
	}
	double *A = v1->GetCoordinate();
    double *B = v2->GetCoordinate();
    double *C = v3->GetCoordinate();

	/*compute the plane of this face*/
    panel[0] = (B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]);
    panel[1] = (B[2] - A[2]) * (C[0] - A[0]) - (C[2] - A[2]) * (B[0] - A[0]);
    panel[2] = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
    // normalize it with a^2 + b^2 + c^2 = 1
    double normalizer = sqrt( panel[0] * panel[0] + panel[1] * panel[1] + panel[2] * panel[2] );
    panel[0] /= normalizer;
    panel[1] /= normalizer;
    panel[2] /= normalizer;
    panel[3] = - (panel[0] * A[0] + panel[1] * A[1] + panel[2] * A[2]);

	/*compute the quadric of this face by plane*/
	quadric[0] = panel[0] * panel[0];
    quadric[1] = panel[0] * panel[1];
    quadric[2] = panel[0] * panel[2];
    quadric[3] = panel[0] * panel[3];
    quadric[4] = panel[1] * panel[1];
    quadric[5] = panel[1] * panel[2];
    quadric[6] = panel[1] * panel[3];
    quadric[7] = panel[2] * panel[2];
    quadric[8] = panel[2] * panel[3];
    quadric[9] = panel[3] * panel[3];

	/*add the value of quadric to each vertex that forms this face*/
	/*each vertex add two times, so we will divide it by two later*/
	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		vQuadric1 = (*e_iter)->GetStartVertex()->GetQuadric();
		vQuadric2 = (*e_iter)->GetEndVertex()->GetQuadric();
		for (int i = 0; i < 10; i++)
		{
			vQuadric1[i] += quadric[i] / 2;
			vQuadric2[i] += quadric[i] / 2;
		}
	}
    return;
}