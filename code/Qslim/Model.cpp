//
//  Model.cpp
//  qslim
//
//  Created by ziyang on 14-8-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#include "Model.h"

/****************Model()*************************/
/**********Construct a model from a SKP model********/
Model::Model(SUModelRef model)
{
	modelSU = model;
	ConvertSKPtoModel();
}

/****************~Model()************************/
Model::~Model()
{
	/*release the SKP mode*/
  	SUInitialize();
	SUModelRelease(&modelSU);
	SUTerminate();

	/*delete all the faces, edges and vertexs*/
	for (list<Face*>::iterator f_iter = faceList.begin(); f_iter != faceList.end(); f_iter++)
	{
		if (*f_iter)
		{
			delete *f_iter;
		}
	}

	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		if (*e_iter)
		{
			delete *e_iter;
		}
	}
	
	for (list<Vertex*>::iterator v_iter = vertexList.begin(); v_iter != vertexList.end(); v_iter++)
	{
		if (*v_iter)
		{
			delete *v_iter;
		}
	}

	/*clear the lists*/
	vertexList.clear();
	edgeList.clear();
	faceList.clear();
}

/***************GLfloat determinant(GLfloat[][], int)***********************************************/
double determinant(double b[][4],int m)
{
    int i,j;
    double sum = 0, c[4][4];
    if(m==2)
    {
        sum = b[0][0]*b[1][1] - b[0][1]*b[1][0];
        return sum;
    }
    for(int p=0;p<m;p++)
    {
        int h = 0,k = 0;
        for(i=1;i<m;i++)
        {
            for( j=0;j<m;j++)
            {
                if(j==p)
                    continue;
                c[h][k] = b[i][j];
                k++;
                if(k == m-1)
                {
                    h++;
                    k = 0;
                }
            }
        }
        sum = sum + b[0][p]*pow(-1,p)*determinant(c,m-1);
    }
    return sum;
}

/**************float round(float, float)**********************/
float round(float orig, float bit)
{
    for (int i = 0; i < bit; i++) {
        orig *= 10;
    }
    
    int temp1 = orig;
    float temp2 = temp1;
    for (int i = 0; i < bit; i++) {
        temp2 /= 10;
    }
    return temp2;
}

/****************void ConvertSKPtoModel******************/
/*******Convert the SKP model to our model struction**********/
void Model::ConvertSKPtoModel()
{
	SUInitialize();

	SUEntitiesRef entities = SU_INVALID;
	size_t faceCount = 0;
	Face *face = NULL;

	/*Get all the edges from the SKP model*/
	SUModelGetEntities(modelSU, &entities);
	SUEntitiesGetNumFaces(entities, &faceCount);

	vector<SUFaceRef> faces(faceCount);
	SUEntitiesGetFaces(entities, faceCount, &faces[0], &faceCount);

	for (int i = 0; i < faceCount; i++)
	{
		/*Add every face to the faceList*/
		AddFace(faces[i]);
	}
	
	SUTerminate();
}

/****************void AddFace()**************************/
/**Add each face to the faceList, and get the edges of this face***/
void Model::AddFace(SUFaceRef faceSU)
{
	SUInitialize();
	size_t edgeCount = 0;
	Face *face = new Face();
	Edge *edge = NULL;
	/*Add the face to the faceList*/
	faceList.push_back(face);
	
	/*Get all the edges of the face*/
	SUFaceGetNumEdges(faceSU, &edgeCount);
	vector<SUEdgeRef> edges(edgeCount);
	SUFaceGetEdges(faceSU, edgeCount, &edges[0], &edgeCount);
	for (int i = 0; i < edgeCount; i++)
	{
		/*Add each edge to the edgeList*/
		edge = AddEdge(edges[i]);
		if (!edge)
		{
			continue;
		}
		/*Add each edge to the edgeList of the face*/
		face->AddEdge(edge);
		edge->AddFace(face);
		edge = NULL;
	}
	face->ComputePanelAndQuadric();
	SUTerminate();
}

/*******Edge* AddEdge(SUEdgeRef)*****************/
/*******Add each edge to the edgeList****************/
Edge* Model::AddEdge(SUEdgeRef edgeSU)
{
	SUInitialize();
	SUVertexRef startSU, endSU;
	SUPoint3D startSU3D, endSU3D;
	Vertex *startVertex = NULL, *endVertex = NULL;
	Edge* edge = NULL;

	/*Get the start vertex and the end vertex of the edge*/
	SUEdgeGetStartVertex(edgeSU, &startSU);
	SUEdgeGetEndVertex(edgeSU, &endSU);
	SUVertexGetPosition(startSU, &startSU3D);
	SUVertexGetPosition(endSU, &endSU3D);

	/*Add this two edge into the vertexList*/
	startVertex = GetVertex(startSU3D.x, startSU3D.y, startSU3D.z);
	endVertex = GetVertex(endSU3D.x, endSU3D.y, endSU3D.z);

	/*Construct an edge between these two vertexs*/
	edge = GetEdge(startVertex, endVertex);
	SUTerminate();
	return edge;
}

/****************void AddPosition()************************/
void Model::ComputeMiddleAndRange(GLdouble xMin, GLdouble xMax, GLdouble yMin, GLdouble yMax, GLdouble zMin, GLdouble zMax)
{
    /*compute the middle point of the model*/
    xMiddle = (xMin + xMax) / 2;
    yMiddle = (yMin + yMax) / 2;
    zMiddle = (zMin + zMax) / 2;
    
    /*compute the range of the model: find the maximum among  x-range, y-range and z-range*/
    range = xMax - xMin;
    if ((xMax - xMin) < (yMax - yMin)){
        range = yMax - yMin;
    }
    if (range < (zMax - zMin)) {
        range = zMax - zMin;
    }
	if (0 == range)
	{
		range = 1;
	}
}

/****************void ComputeModelPosition()**************/
void Model::ComputeModelPosition()
{
	double xMin = BIGNUM, xMax = 0 - BIGNUM, yMin = BIGNUM, yMax = 0 - BIGNUM;
	zMin = BIGNUM, zMax = 0 - BIGNUM;
	double *coor = NULL;
	for (list<Vertex*>::iterator v_iter = vertexList.begin(); v_iter != vertexList.end(); v_iter++)
	{
		coor = (*v_iter)->GetCoordinate();
		if (coor[0] <xMin){
			xMin = coor[0];
		}else if(coor[0] > xMax){
			xMax = coor[0];
		}
		if (coor[1] < yMin){
			yMin = coor[1];
		}else if(coor[1] > yMax){
			yMax = coor[1];
		}
		if (coor[2] < zMin){
			zMin = coor[2];
		}else if(coor[2] > zMax){
			zMax = coor[2];
		}
	}
	ComputeMiddleAndRange(xMin, xMax, yMin, yMax, zMin, zMax);
}

/****************void DrawModel()***********************/
void Model::DrawModel()
{
	ComputeModelPosition();
	/*translate and scale the model*/
	GLdouble scale = 1/range;
	glScaled(scale, scale, scale);
	glTranslated(-xMiddle, -yMiddle, -zMiddle);

	/*draw the model*/
	glBegin(GL_LINES);
	for (list<Edge *>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		glVertex3dv((*e_iter)->GetStartVertex()->GetCoordinate());
		glVertex3dv((*e_iter)->GetEndVertex()->GetCoordinate());
	}
	glEnd();

	/*scale and translate the model back*/
	glTranslated(xMiddle, yMiddle, zMiddle);
	glScaled(range, range, range);
}

/****************void SaveModel()************************/
void Model::SaveModel()
{
} 

/****************void DecimateModel*********************/
void Model::DecimateModel(int number)//decimate the model
{
	if (number < 0 || number > edgeList.size())
	{
		cout<<"Wrong input!"<<endl;
	}
	/*each time decimate one edge*/
	for (int i = 0; i < number; i++)
	{
		FindEdgeToDecimate();
	}
	srand(time(0));
}

void Model::FindEdgeToDecimate()
{
	const int number = 8;
	int randNum[8] = {0}, count = 0, found = 0;
	double error = 0, minError = BIGNUM, coor[4], minCoor[4];
	Edge* minEdge = NULL, *edge = NULL, * newEdge = NULL;
	Vertex* newVertex = NULL, *startVertex = NULL, *endVertex = NULL;

	/*get the edge that will generate least error if decimate*/
	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++, count++)
	{
		error = ComputeVertexAndError(*e_iter, coor);
		if (error < minError)
		{
			minError = error;
			minEdge = *e_iter;
			for (int j = 0; j < 4; j++)
			{
				minCoor[j] = coor[j];
			}
		}
	}

	/*decimate the edge, and update all the associated values*/
	newVertex = GetVertex(minCoor[0], minCoor[1], minCoor[2]);//create the new vertex
	startVertex = minEdge->GetStartVertex();
	endVertex = minEdge->GetEndVertex();
	//set the quadric of the new vertex
	double* quadric = newVertex->GetQuadric();
	for (int i = 0; i < 10; i++)
	{
		quadric[i] = startVertex->GetQuadric()[i] + endVertex->GetQuadric()[i];
	}

	//update the values according to the start vertex
	if (newVertex != startVertex)
	{
		for (list<Vertex*>::iterator v_iter = startVertex->GetVertexList()->begin();
			v_iter != startVertex->GetVertexList()->end(); v_iter++)
		{
			if (*v_iter == endVertex)
			{
				continue;
			}
			if (*v_iter == newVertex)
			{
				edge = GetEdge(endVertex, *v_iter);
				edgeList.remove(edge);
				for (list<Face*>::iterator f_iter = edge->GetFaceList()->begin(); f_iter != edge->GetFaceList()->end(); f_iter++)
				{
					(*f_iter)->GetEdgeList()->remove(edge);
					if ((*f_iter)->GetEdgeList()->size() < 3)
					{
						delete *f_iter;
					}
				}
				(*v_iter)->GetVertexList()->remove(startVertex);
				delete edge;
				continue;
			}
			edge = GetEdge(startVertex, *v_iter);
			newEdge = GetEdge(newVertex, *v_iter);
			edgeList.remove(edge);
			for (list<Face*>::iterator f_iter = edge->GetFaceList()->begin(); f_iter != edge->GetFaceList()->end(); f_iter++)
			{
				(*f_iter)->GetEdgeList()->remove(edge);
				(*f_iter)->AddEdge(newEdge);
				newEdge->AddFace(*f_iter);
			}
			(*v_iter)->GetVertexList()->remove(startVertex);
			delete edge;
		}
	}

	//update the values according to the end vertex
	if (newVertex != endVertex)
	{
		for (list<Vertex*>::iterator v_iter = endVertex->GetVertexList()->begin();
			v_iter != endVertex->GetVertexList()->end(); v_iter++)
		{
			if (*v_iter == startVertex)
			{
				continue;
			}
			if (*v_iter == newVertex)
			{
				edge = GetEdge(endVertex, *v_iter);
				edgeList.remove(edge);
				for (list<Face*>::iterator f_iter = edge->GetFaceList()->begin(); f_iter != edge->GetFaceList()->end(); f_iter++)
				{
					(*f_iter)->GetEdgeList()->remove(edge);
					if ((*f_iter)->GetEdgeList()->size() < 3)
					{
						delete *f_iter;
					}
				}
				(*v_iter)->GetVertexList()->remove(endVertex);
				delete edge;
				continue;
			}
			edge = GetEdge(endVertex, *v_iter);
			newEdge = GetEdge(newVertex, *v_iter);
			edgeList.remove(edge);
			for (list<Face*>::iterator f_iter = edge->GetFaceList()->begin(); f_iter != edge->GetFaceList()->end(); f_iter++)
			{
				(*f_iter)->GetEdgeList()->remove(edge);
				(*f_iter)->AddEdge(newEdge);
				newEdge->AddFace(*f_iter);
			}
			(*v_iter)->GetVertexList()->remove(endVertex);
			delete edge;
		}
	}

	//update the faces that are associated with this edge
	for (list<Face*>::iterator f_iter = minEdge->GetFaceList()->begin(); f_iter != minEdge->GetFaceList()->end(); f_iter++)
	{
		(*f_iter)->GetEdgeList()->remove(minEdge);
		if ((*f_iter)->GetEdgeList()->size() < 4)
		{
			faceList.remove(*f_iter);
			delete *f_iter;
		}
	}

	//update this edge and the start vertex and the end vertex
	edgeList.remove(minEdge);
	startVertex->GetVertexList()->remove(endVertex);
	endVertex->GetVertexList()->remove(startVertex);
	delete minEdge;
	if (startVertex != newVertex)
	{
		vertexList.remove(startVertex);
		delete startVertex;
	}
	if (endVertex != newVertex)
	{
		vertexList.remove(endVertex);
		delete endVertex;
	}
}

/***double Model::ComputeVertexAndError(Edge* edge, double* coor)***/
/*return the error if decimate this edge*/
/*edge is the edge to decimate*/
/*coor is the coordinate of the new vetex*/
double Model::ComputeVertexAndError(Edge* edge, double* coor)
{
	double error = 0;
	double Q[10] = {0}, QTemp[4][4], matric[4][4];
	Vertex* startVertex = NULL, *endVertex = NULL;

	startVertex = edge->GetStartVertex();
	endVertex = edge->GetEndVertex();

	for (int i = 0; i < 10; i++)
	{
		Q[i] = startVertex->GetQuadric()[i] + endVertex->GetQuadric()[i];
	}
	QTemp[0][0] = Q[0];
    QTemp[0][1] = Q[1];
    QTemp[0][2] = Q[2];
    QTemp[0][3] = Q[3];
    QTemp[1][0] = Q[1];
    QTemp[1][1] = Q[4];
    QTemp[1][2] = Q[5];
    QTemp[1][3] = Q[6];
    QTemp[2][0] = Q[2];
    QTemp[2][1] = Q[5];
    QTemp[2][2] = Q[7];
    QTemp[2][3] = Q[8];
    QTemp[3][0] = Q[3];
    QTemp[3][1] = Q[6];
    QTemp[3][2] = Q[8];
    QTemp[3][3] = Q[9];
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            matric[i][j] = QTemp[i][j];
        }
    }
    matric[3][0] = matric[3][1] = matric[3][2] = 0;
    matric[3][3] = 1;

	double det = round(determinant(matric, 4),8);
	if ((0 != det) && ComputeCoordinate(Q, coor)) {//if the inverse of the matrix exist
		error = ComputeError(coor, QTemp);     
    }else{//if the inverse of the matrix does not exist
		double minError = 0;
		double newCoor[4] = {0}, *startCoor, *endCoor;

		startCoor = startVertex->GetCoordinate();
		endCoor = endVertex->GetCoordinate();
		/*compute the coordinate of the middle point of this edge*/
		for (int i = 0; i < 4; i++)
		{
			newCoor[i] = (startCoor[i] + endCoor[i])/2;
		}

		/*compute the error if we decimate the edge to the start vertex*/
		minError = error = ComputeError(startCoor, QTemp);
		for (int i = 0; i < 4; i++)
		{
			coor[i] = startCoor[i];
		}

		/*compute the error if we decimate the edge to the end vertex*/
		minError = ComputeError(endCoor, QTemp);
		if (minError < error)
		{
			for (int i = 0; i < 4; i++)
			{
				coor[i] = endCoor[i];
			}
			error = minError;
		}

		/*compute the error if we decimate the edge to the middle point*/
		minError = ComputeError(newCoor, QTemp);
		if(minError < error)
		{
			for (int i = 0; i < 4; i++)
			{
				coor[i] = newCoor[i];
			}
			error = minError;
		}
	}
	return error;
}

/****************Vertex GetVertex(double, double, double)***********/
Vertex* Model::GetVertex(double x, double y, double z)
{
	double *coor = NULL;
	Vertex* vertex = NULL;
	for (list<Vertex*>::iterator v_iter = vertexList.begin(); v_iter != vertexList.end(); v_iter++)
	{
		coor = (*v_iter)->GetCoordinate();
		if (z > coor[2] + ERROR)
		{
			continue;
		}//end if
		else if (z < coor[2] - ERROR)
		{
			vertex = new Vertex(x, y, z);
			vertexList.insert(v_iter, vertex);
			return vertex;
		}//end else if
		else
		{
			if (x > coor[0] + ERROR)
			{
				continue;
			}//end else-if
			else if (x < coor[0] - ERROR)
			{
				vertex = new Vertex(x, y, z);
				vertexList.insert(v_iter, vertex);
				return vertex;
			}//end else-else if
			else {
				if (y > coor[1] + ERROR)
				{
					continue;
				}//end else-else-if
				else if (y < coor[1] - ERROR)
				{
					vertex = new Vertex(x, y, z);
					vertexList.insert(v_iter, vertex);
					return vertex;
				}//end else-else-else if
				else
				{
					return *v_iter;
				}//end else-else-else
			}//end else-else
		}//end else
	}//end for
	vertex = new Vertex(x, y, z);
	vertexList.push_back(vertex);
	return vertex;
}

/****************Edge* GetEdge(Vertex*, Vertex*)***********/
Edge* Model::GetEdge(Vertex* v1, Vertex* v2)
{
	Vertex* startVertex = NULL, *endVertex = NULL;
	Edge* edge = NULL;
	int compareResult = 0;
	compareResult = CompareVertexCoordinate(v1, v2) ;
	if (compareResult >0)
	{
		startVertex = v2;
		endVertex = v1;
	}
	else if(compareResult < 0)
	{
		startVertex = v1;
		endVertex = v2;
	}
	else
	{
		return NULL;
	}

	for (list<Edge*>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		compareResult = CompareVertexCoordinate(startVertex, (*e_iter)->GetStartVertex());
		if (1 == compareResult)
		{
			continue;
		}//end if
		else if (-1 == compareResult)
		{
			edge = new Edge(startVertex, endVertex);
			edgeList.insert(e_iter, edge);
			return edge;
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
				edge = new Edge(startVertex, endVertex);
				edgeList.insert(e_iter, edge);
				return edge;
			}//end else-else if
			else 
			{
				return *e_iter;
			}//end else-else
		}//end else 
	}//end for
	edge = new Edge(startVertex, endVertex);
	edgeList.push_back(edge);
	return edge;
}

/**************ComputeCoordinate(float*, float *)**************************/
/***Compute the optimal Coordinate of the new vertex**************************/
bool ComputeCoordinate(double* Q, double* coordinate)
{
     double d, x, y, z;
     d = Q[0]*Q[4]*Q[7] + Q[1]*Q[5]*Q[2] + Q[2] * Q[1] * Q[5] - Q[2]*Q[4]*Q[2] - Q[0]*Q[5]*Q[5] - Q[1]*Q[1]*Q[7];
     if (d == 0) {
         return false;
     }
    x = 0-Q[3] * Q[4] * Q[7] - Q[1]*Q[5]*Q[8] - Q[2] * Q[6] * Q[5] + Q[2]*Q[4]*Q[8] + Q[3]*Q[5]*Q[5] + Q[1]*Q[6]*Q[7];
    y = 0- Q[0]*Q[6]*Q[7] - Q[3]*Q[5]*Q[2] - Q[2]*Q[1]*Q[8] + Q[2]*Q[6]*Q[2] + Q[3]*Q[1]*Q[7] + Q[0]*Q[5]*Q[8];
    z = 0-Q[0]*Q[4]*Q[8] - Q[1]*Q[6]*Q[2] -Q[3]*Q[5]*Q[1] +Q[3]*Q[4]*Q[2] + Q[1]*Q[1]*Q[8] + Q[0]*Q[6]*Q[5];
    
    coordinate[0] = x/d;
    coordinate[1] = y/d;
    coordinate[2] = z/d;
    coordinate[3] = 1;
    return true;
}

double ComputeError(double coordinate[4], double Q[][4])
{
    double error = 0, errorTemp[4] = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            errorTemp[i] += coordinate[j] * Q[j][i];
        }
        error += errorTemp[i] * coordinate[i];
    }
    return error;
}