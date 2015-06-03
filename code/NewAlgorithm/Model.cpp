//
//  Model.cpp
//  finalproject
//
//  Created by Ziyang on 14-8-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#include "Model.h"
double H = 500;
double DISTANCE = 8000;

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
	/*Sort the faceList*/
	faceList.sort(MergeFaceList);
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
		/*Add each edge to the edgeList of the face*/
		face->AddEdge(edge);
		edge = NULL;
	}
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
	startVertex = GetVertex(startSU3D.x, startSU3D.y, startSU3D.z, &vertexList);
	endVertex = GetVertex(endSU3D.x, endSU3D.y, endSU3D.z, &vertexList);

	/*Construct an edge between these two vertexs*/
	edge = GetEdge(startVertex, endVertex, &edgeList);
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

	/*draw the mesh*/
	glBegin(GL_LINES);
	for (list<Edge *>::iterator e_iter = edgeList.begin(); e_iter != edgeList.end(); e_iter++)
	{
		glVertex3dv((*e_iter)->GetStartVertex()->GetCoordinate());
		glVertex3dv((*e_iter)->GetEndVertex()->GetCoordinate());
	}
	glEnd();

	/*translate and scale the mesh back*/
	glTranslated(xMiddle, yMiddle, zMiddle);
	glScaled(range, range, range);
}

/****************void SaveModel()************************/
void Model::SaveModel()
{
} 

/****************void DecimateModel*********************/
void Model::DecimateModel(double h, double e)//decimate the model
{
	H = h;
	DISTANCE = e;
	list<Face*> resultFaceList;
	list<Edge*> resultEdgeList;
	list<Vertex*> resultVertexList;

	list<Vertex*> compareSectionVertex, currentSectionVertex;
	list<Edge*> compareSectionEdge, currentSectionEdge;

	double distance = 0;

	/*Get the bottom of this model*/
	GetSection(zMin, &currentSectionEdge, &currentSectionVertex);
	compareSectionVertex = currentSectionVertex;
	compareSectionEdge = currentSectionEdge;

	for (double currHeight = zMin + H; currHeight < zMax + H; currHeight += H)
	{
		currentSectionVertex.clear();
		currentSectionEdge.clear();
		/*For each height get a section*/
		if (currHeight >= zMax)
		{
			currHeight = zMax;
		}
		GetSection(currHeight, &currentSectionEdge, &currentSectionVertex);
		/*Get the corresponding point from the compare section*/
		distance = GetCorrPoint(&compareSectionVertex, &currentSectionVertex, &compareSectionEdge, &currentSectionEdge);
		/*if the error is bigger than the threshold*/
		if (distance > DISTANCE || currHeight == zMax)
		{
			/*build the outline surface*/
			AddCylinder(&compareSectionEdge, &resultVertexList, &resultEdgeList, &resultFaceList);
			UpdateCurrentSection(&currentSectionEdge, &currentSectionVertex);

			/*delete the compare profile*/
			for (list<Edge*>::iterator e_iter = compareSectionEdge.begin();
				e_iter != compareSectionEdge.end();
				e_iter++)
			{
				if (*e_iter)
				{
					delete *e_iter;
				}
			}

			for (list<Vertex*>::iterator v_iter = compareSectionVertex.begin();
				v_iter != compareSectionVertex.end();
				v_iter++)
			{
				if (*v_iter)
				{
					delete *v_iter;
				}
			}
			compareSectionVertex.clear();
			compareSectionEdge.clear();
			/*set the current profile as the compare profile*/
			compareSectionVertex = currentSectionVertex;
			compareSectionEdge = currentSectionEdge;
		}
		SetSectionToDefault(&compareSectionEdge, &compareSectionVertex);
	} 

	/*replace the face list, edge list and vertex list*/
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

	faceList.clear();
	edgeList.clear();
	vertexList.clear();

	resultFaceList.sort(MergeFaceList);
	faceList = resultFaceList;
	edgeList = resultEdgeList;
	vertexList = resultVertexList;
}

/***************list<Face*> GetSection(double)*****************/
void Model::GetSection(double height, list<Edge*> *resultEdgeList, list<Vertex*> *resultVertexList)
{
	/*Get all the intersection points*/
	for (list<Face*>::iterator f_iter = faceList.begin(); f_iter != faceList.end(); f_iter++)
	{
		/*Get all the faces that intersect with current height*/
		if (((*f_iter)->GetZMin() < height + ERROR) && ((*f_iter)->GetZMax() > height - ERROR))
		{
			GetIntersectionOfFace((*f_iter), height, resultEdgeList, resultVertexList);
		}
		else if((*f_iter)->GetZMin() > height + ERROR)
		{
			break;
		}
	}
}

/**********list<Edge*>* GetIntersectionOfFace(double)***********/
void  Model::GetIntersectionOfFace(Face* face, double height, list<Edge*>* resultEdgeList, list<Vertex*>* resultVertexList)
{
	list<Edge*> *edgeListOfFace = NULL;
	list<Vertex*> currVertexList;
	Vertex* startVertex = NULL, *endVertex = NULL, *currVertex = NULL;
	Edge* edge = NULL;
	double* coorStart = NULL, *coorEnd = NULL, x, y;

	/*Get all the intersection points*/
	edgeListOfFace = face->GetEdgeList();
	for (list<Edge*>::iterator e_iter = (*edgeListOfFace).begin(); 
		e_iter != (*edgeListOfFace).end();
		e_iter++)
	{
		startVertex = (*e_iter)->GetStartVertex();
		endVertex = (*e_iter)->GetEndVertex();
		coorStart = startVertex->GetCoordinate();
		coorEnd = endVertex->GetCoordinate();

		if (((height - ERROR) < coorStart[2]) && (coorStart[2] < (height + ERROR)) 
			&& ((height - ERROR) < coorEnd[2]) && (coorEnd[2] < (height + ERROR)))
		{
			startVertex = GetVertex(coorStart[0], coorStart[1], height, resultVertexList);
			endVertex = GetVertex(coorEnd[0], coorEnd[1], height, resultVertexList);
			GetEdge(startVertex, endVertex, resultEdgeList);
		}
		else if ((coorStart[2] < (height + ERROR)) && (coorEnd[2] > (height - ERROR)))
		{
			x = (height - coorStart[2]) / (coorEnd[2] - coorStart[2]) * (coorEnd[0] - coorStart[0]) + coorStart[0];
			y = (height - coorStart[2]) / (coorEnd[2] - coorStart[2]) * (coorEnd[1] - coorStart[1]) + coorStart[1];
			currVertex = GetVertex(x, y, height, resultVertexList);
			currVertexList.push_back(currVertex);
		}
		else if (startVertex->GetCoordinate()[2] > (height + ERROR))
		{
			break;
		}
		startVertex = NULL;
		endVertex = NULL;
		coorStart = NULL;
		coorEnd = NULL;
	}

	if (currVertexList.size() > 0)
	{
		currVertexList.sort(MergeVertexList);
		/*Get the intersection edges*/
		for (list <Vertex*>::iterator v_iter = currVertexList.begin(); v_iter != currVertexList.end(); v_iter++)
		{
			startVertex = GetVertex((*v_iter)->GetCoordinate()[0], (*v_iter)->GetCoordinate()[1], height, resultVertexList);
			v_iter++;
			endVertex = GetVertex((*v_iter)->GetCoordinate()[0], (*v_iter)->GetCoordinate()[1], height, resultVertexList);
			if (CompareVertexCoordinate(startVertex, endVertex) != 0)
			{
				GetEdge(startVertex, endVertex, resultEdgeList);
			}
			startVertex = NULL;
			endVertex = NULL;
		}
		currVertexList.clear();
	}
	return;
}

/*void UpdateNewModel(list<Vertex*>* vList1, list<Vertex*> vList2, list<Edge*>* eList1, list<Edge*>* eList2)*/
/*vList1: original vertexList, vList2:current vertexList
	eList1:original edgeList, eList2:current edgeList*/
void Model::UpdateNewModel(list<Vertex*>* vList1, list<Vertex*>* vList2, list<Vertex*>* vList3, 
						   list<Edge*>* eList1, list<Edge*>* eList2, list<Edge*>* eList3)
{
	Vertex* vertex = NULL, *corrVertex = NULL;
	list<Vertex*> newVertexList;
	list<Edge*> newEdgeList;
	for (list<Vertex*>::iterator v_iter = vList2->begin(); v_iter != vList2->end(); v_iter ++)
	{
		corrVertex = (*v_iter)->GetCorrPoint();
		if ((*v_iter)->GetDistance() > DISTANCE/10)//the error is big enough,
		{//currently, we set the threshold as DISTANCE/10
			/*put this vertex into the result vertex list*/
			vertex =GetVertex(corrVertex->GetCoordinate()[0], 
				corrVertex->GetCoordinate()[1],
				(*v_iter)->GetCoordinate()[2],
				vList3);
			/*put the edge into the result edge list*/
			GetEdge((*v_iter)->GetCorrPoint(), vertex, eList3);

			vertex =GetVertex(corrVertex->GetCoordinate()[0], 
				corrVertex->GetCoordinate()[1],
				(*v_iter)->GetCoordinate()[2],
				&newVertexList);
			vertex->SetCorrPoint(*v_iter);
			(*v_iter)->SetCorrPoint(vertex);
		}
		vertex = NULL;
	}
	AddNewEdge(&newVertexList, &newEdgeList);

}

/****************Vertex GetVertex(double, double, double)***********/
Vertex* GetVertex(double x, double y, double z, list<Vertex*>* vertexList)
{
	double *coor = NULL;
	Vertex* vertex = NULL;
	for (list<Vertex*>::iterator v_iter = (*vertexList).begin(); v_iter != (*vertexList).end(); v_iter++)
	{
		coor = (*v_iter)->GetCoordinate();
		if (z > coor[2] + ERROR)
		{
			continue;
		}//end if
		else if (z < coor[2] - ERROR)
		{
			vertex = new Vertex(x, y, z);
			(*vertexList).insert(v_iter, vertex);
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
				(*vertexList).insert(v_iter, vertex);
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
					(*vertexList).insert(v_iter, vertex);
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
	(*vertexList).push_back(vertex);
	return vertex;
}

/****************Edge* GetEdge(Vertex*, Vertex*)***********/
Edge* GetEdge(Vertex* v1, Vertex* v2, list<Edge*>* edgeList)
{
	Vertex* startVertex = NULL, *endVertex = NULL;
	Edge* edge = NULL;
	int compareResult = 0;
	if (CompareVertexCoordinate(v1, v2) >0)
	{
		startVertex = v2;
		endVertex = v1;
	}
	else
	{
		startVertex = v1;
		endVertex = v2;
	}

	for (list<Edge*>::iterator e_iter = (*edgeList).begin(); e_iter != (*edgeList).end(); e_iter++)
	{
		compareResult = CompareVertexCoordinate(startVertex, (*e_iter)->GetStartVertex());
		if (1 == compareResult)
		{
			continue;
		}//end if
		else if (-1 == compareResult)
		{
			edge = new Edge(startVertex, endVertex);
			(*edgeList).insert(e_iter, edge);
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
				(*edgeList).insert(e_iter, edge);
				return edge;
			}//end else-else if
			else 
			{
				return *e_iter;
			}//end else-else
		}//end else 
	}//end for
	edge = new Edge(startVertex, endVertex);
	(*edgeList).push_back(edge);
	return edge;
}

/************double GetCorrPoing(list<Vertex*>* vList1, list<Vertex*>*vList2)********/
/************vList1 is the original list, vList2 is the current list********************/
double GetCorrPoint(list<Vertex*>* vList1, list<Vertex*>* vList2, list<Edge*>* eList1, list<Edge*>* eList2)
{
	double distance = 0, minDistance = 0, totalDistance = 0;
	Vertex* tempVertex, *closestVertex = NULL;
	list<Vertex*> tempVertexListCurr, tempVertexListOrig;
	/* Find the corresponding point of the original list*/
	for (list<Vertex*>::iterator v_iter1 = vList1->begin(); v_iter1 != vList1->end(); v_iter1++)
	{
		minDistance =  (*v_iter1)->GetDistance();
		closestVertex = (*v_iter1)->GetCorrPoint();
		for (list<Vertex*>::iterator v_iter2 = vList2->begin(); v_iter2 != vList2->end(); v_iter2++)
		{
			distance = Get2DDist((*v_iter1)->GetCoordinate(), (*v_iter2)->GetCoordinate());
			if (distance < minDistance)
			{
				minDistance = distance;
				closestVertex = *v_iter2;
			}//end for-for-if
		}//end for-for

		if (closestVertex->GetCorrPoint())
		{
			tempVertex = new Vertex((closestVertex)->GetCoordinate());
			tempVertexListCurr.push_back(tempVertex);
		}//end for-if
		else
		{
			tempVertex = closestVertex;
		}//end for-else
		(*v_iter1)->SetDistance(minDistance);
		(*v_iter1)->SetCorrPoint(tempVertex);
		tempVertex->SetDistance(minDistance);
		tempVertex->SetCorrPoint(*v_iter1);
		totalDistance += minDistance;
		tempVertex = NULL;
		closestVertex = NULL;
		distance = minDistance = 0;
	}//end for

	/* Find the corresponding point of the current list*/
	for (list<Vertex*>::iterator v_iter2 = vList2->begin(); v_iter2 != vList2->end(); v_iter2++)
	{
		if (!(*v_iter2)->GetCorrPoint())
		{
			minDistance = (*v_iter2)->GetDistance();
			closestVertex = *v_iter2;
			for (list<Vertex*>::iterator v_iter1 = vList1->begin(); v_iter1 != vList1->end(); v_iter1++)
			{
				distance = Get2DDist((*v_iter2)->GetCoordinate(), (*v_iter1)->GetCoordinate());
				if (distance < minDistance)
				{
					minDistance = distance;
					closestVertex = *v_iter1;
				}
			}
			tempVertex = new Vertex(closestVertex->GetCoordinate());
			tempVertexListOrig.push_back(tempVertex);
			(*v_iter2)->SetDistance(minDistance);
			(*v_iter2)->SetCorrPoint(tempVertex);
			tempVertex->SetDistance(minDistance);
			tempVertex->SetCorrPoint(*v_iter2);
			totalDistance += minDistance;
			tempVertex = NULL;
			closestVertex = NULL;
			distance = minDistance = 0;
		}//end for-if
	}//end for

	if (tempVertexListOrig.size() > 0)
	{
		tempVertexListOrig.sort(MergeVertexList);
		vList1->sort(MergeVertexList);
		vList1->merge(tempVertexListOrig, MergeVertexList);
		AddNewEdge(&tempVertexListOrig, eList1);
	}
	if (tempVertexListCurr.size() > 0)
	{
		tempVertexListCurr.sort(MergeVertexList);
		vList2->sort(MergeVertexList);
		vList2->merge(tempVertexListCurr, MergeVertexList);
		AddNewEdge(&tempVertexListCurr, eList2);
	}
	return totalDistance;
}

/**********double GetDist(double* coor, double* coor)********************/
double Get2DDist(double *coor1, double* coor2)
{
	return sqrt((coor1[0] - coor2[0]) * (coor1[0] - coor2[0]) +(coor1[1] - coor2[1]) * (coor1[1] - coor2[1]));
}

/**********double AddNewEdge(list<Vertex*>*, list<Edge*>*)*******************/
void AddNewEdge(list<Vertex*>* vList, list<Edge*>* eList)
{
	list<Vertex*>::iterator v_iterTemp;
	Vertex* tempVertex = NULL;
	for (list<Vertex*>::iterator v_iter1 = vList->begin(); v_iter1 != vList->end(); v_iter1++)
	{
		tempVertex = (*v_iter1)->GetCorrPoint();
		v_iterTemp = tempVertex->GetVertexList()->end();
		v_iterTemp--;
		GetEdge(*v_iter1, (*v_iterTemp)->GetCorrPoint(), eList);
		v_iterTemp--;
		GetEdge(*v_iter1, (*v_iterTemp)->GetCorrPoint(), eList);
	}
}

/*******void AddSide(list<Edge*>*, list<Vertex*>*, list<Edge*>*, list<Face*>*)******/
void AddCylinder(list<Edge*>* edgeList, list<Vertex*>* resultVertexList, list<Edge*>* resultEdgeList, list<Face*>* resultFaceList)
{
	Face *face = NULL, *faceTop, *faceBottom = NULL;
	Edge *edge = NULL;
	Vertex* vertexStart = NULL, *vertexEnd = NULL, *vertexStart2 = NULL, * vertexEnd2 = NULL;

	faceTop = new Face();
	faceBottom = new Face();
	resultFaceList->push_back(faceBottom);

	for (list<Edge*>::iterator e_iter = edgeList->begin();
		e_iter != edgeList->end();
		e_iter++)
	{
		face = new Face();
		vertexStart = GetVertex((*e_iter)->GetStartVertex()->GetCoordinate()[0], 
			(*e_iter)->GetStartVertex()->GetCoordinate()[1], 
			(*e_iter)->GetStartVertex()->GetCoordinate()[2],
			resultVertexList);

		vertexEnd = GetVertex((*e_iter)->GetEndVertex()->GetCoordinate()[0],
			(*e_iter)->GetEndVertex()->GetCoordinate()[1],
			(*e_iter)->GetEndVertex()->GetCoordinate()[2],
			resultVertexList);

		edge = GetEdge(vertexStart, vertexEnd, resultEdgeList);
		face->AddEdge(edge);//Add the bottom edge to the side face
		faceBottom->AddEdge(edge);//Add the bottom edge into the bottom face

		//build two side edges, add to the side face
		vertexStart2 = GetVertex((*e_iter)->GetStartVertex()->GetCoordinate()[0],
			(*e_iter)->GetStartVertex()->GetCoordinate()[1],
			(*e_iter)->GetStartVertex()->GetCorrPoint()->GetCoordinate()[2],
			resultVertexList);

		edge = GetEdge(vertexStart, vertexStart2, resultEdgeList);
		face->AddEdge(edge);

		vertexEnd2 = GetVertex((*e_iter)->GetEndVertex()->GetCoordinate()[0],
			(*e_iter)->GetEndVertex()->GetCoordinate()[1],
			(*e_iter)->GetEndVertex()->GetCorrPoint()->GetCoordinate()[2],
			resultVertexList);

		edge = GetEdge(vertexEnd, vertexEnd2, resultEdgeList);
		face->AddEdge(edge);

		edge = GetEdge(vertexStart2, vertexEnd2, resultEdgeList);
		face->AddEdge(edge);//add the top edge to the side face
		faceTop->AddEdge(edge);//add the top edge to the top face
		resultFaceList->push_back(face);
	}
	resultFaceList->push_back(faceTop);
}

/******void UpdateSection(list<Edge*>*, list<Vertex*>*)********/
void UpdateCurrentSection(list<Edge*>* edgeList, list<Vertex*>* vertexList)
{
	Vertex *startVertex = NULL, *endVertex = NULL;
	for (list<Vertex*>::iterator v_iter = vertexList->begin();
		v_iter != vertexList->end();
		v_iter++)
	{
		if ((*v_iter)->GetDistance() < (DISTANCE/16))
		{
			(*v_iter)->SetCoordinate((*v_iter)->GetCorrPoint()->GetCoordinate()[0],
				(*v_iter)->GetCorrPoint()->GetCoordinate()[1],
				(*v_iter)->GetCoordinate()[2]);
		}
	}
}

/***********void SetSectionToDefault(list<Edge*>* edgeList, list<Vertex*>* vertexList)***********/
void SetSectionToDefault(list<Edge*>* edgeList, list<Vertex*>* vertexList)
{
	list<Vertex*> newVertexList;
	list<Edge*> newEdgeList;
	Vertex* startVertexOld = NULL, *endVertexOld = NULL, *startVertexNew = NULL, *endVertexNew = NULL;
	double* coordinate = NULL;
	for (list<Edge*>::iterator e_iter = edgeList->begin(); e_iter != edgeList->end(); e_iter++)
	{
		startVertexOld = (*e_iter)->GetStartVertex();
		endVertexOld = (*e_iter)->GetEndVertex();
		if (0 != CompareVertexCoordinate(startVertexOld, endVertexOld))
		{
			coordinate = startVertexOld->GetCoordinate();
			startVertexNew = GetVertex(coordinate[0], coordinate[1], coordinate[2], &newVertexList);
			coordinate = endVertexOld->GetCoordinate();
			endVertexNew = GetVertex(coordinate[0], coordinate[1], coordinate[2], &newVertexList);
			GetEdge(startVertexNew, endVertexNew, &newEdgeList);
		}
	}

	for (list<Edge*>::iterator e_iter = edgeList->begin(); e_iter != edgeList->end(); e_iter++)
	{
		if (*e_iter)
		{
			delete *e_iter;
		}
	}

	for (list<Vertex*>::iterator v_iter = vertexList->begin(); v_iter != vertexList->end(); v_iter++)
	{
		if (*v_iter)
		{
			delete *v_iter;
		}
	}

	vertexList->clear();
	edgeList->clear();
	*vertexList = newVertexList;
	*edgeList = newEdgeList;
}