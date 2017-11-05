/**
 * clusters.c
 * Author: Ting Li <tingli3@illinois.edu>
 * Date: 08/07/2017
 */ 


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "clusters.h"

/**
 * NAME:	PossionTest
 * DESCRIPTION:	calculate the probability to get a value equal or larger than nP under a Poisson (lambda) distribution
 * PARAMETERS:
 * 	int nP:	the value from Poisson distribution
 * 	double lambda: the mean of Poisson distribution
 * RETURN:
 * 	TYPE:	double
 * 	VALUE:	the probability to get a value equal or larger than nP
 */
double PossionTest(int nP, double lambda)
{
	double sum = 1.0;
	double element = 1;
	for(int i = 1; i < nP; i++)
	{
		element = element * lambda / i;
		sum += element;
	}

	sum *= exp(-lambda);
	return 1 - sum;
}

/**
 * NAME:	BinomialTest
 * DESCRIPTION:	calculate the probability to get equal or more cases than nCas under a Binomial (nCas, (nCas+nCon), p) distribution
 * PARAMETERS:
 * 	int nCas: the number of cases
 * 	int nCon: the number of controls
 * 	double p: the p of Binomial distribution (e.g., the probability of any point to be a case)
 * RETURN:
 * 	TYPE:	double
 * 	VALUE:	the probability to get a value equal or larger than nCas
 */
double BinomialTest(int nCas, int nCon, double p)
{
	double q = 1 - p;
	int n = nCas + nCon;
	double logElement = n * log(q);
	 
	double sum = exp(logElement);

	for(int i = 1; i < nCas; i++)
	{
		logElement = logElement + log(n+1-i) + log(p) - log(i) - log(q);
		sum += exp(logElement);
	}
	return 1 - sum;
}

/**
 * NAME:	doClusterPoi
 * DESCRIPTION:	cluster all event points based on a Possion Test
 * PARAMETERS:
 * 	double * x: 		the array of points' X values
 * 	double * y: 		the array of points' Y values
 * 	int * ind:			the array of points' type indicator (1: events, 0: background)
 * 	int * index:		the index of all event points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int countB:			the number of background points
 *	int countE:			the number of event points
 *	int * eC:			the number of event (1) points (within radius) near each point
 *	double * lambda:	the local lambda of Possion distribution of each event points
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 *	struct clusterInfo ** pCInfo: the resulting output clusterInfo
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	the cluster ID of each point
 */
int * doClusterPoi(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countB, int countE, int * eC, double * lambda, double significance, int minCore, bool nonCorePoints, struct clusterInfo ** pCInfo)
{
	int count = index[nBlockX * nBlockY];

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		if(PossionTest(eC[i], lambda[i]) < significance)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	int * inCluster;

	if(NULL == (inCluster = (int *)malloc(sizeof(int) * count))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	
	for(int i = 0; i < count; i++) {
		inCluster[i] = -1;
	}

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;
	int iNbEnd;

	int coreCount;

	int nEInCluster;
	int nBInCluster;
	struct clusterInfo * curInfo;
//	printf("ClusterID,Events,expEvents,LL\n");

	for(int i = 0; i < count; i++)
	{
		if(clusterID[i] != 0 || ind[i] == 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;

		inCluster[i] = cID;
		nEInCluster = 1;
		nBInCluster = 0;

		while(nPToDo > 0) {
			nPToDo --;
			cX = x[pointsToDo[nPToDo]];		
			cY = y[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = index[row * nBlockX + colMin]; iNb < index[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(inCluster[iNb] != cID) {
						if(dist2 >= ((x[iNb] - cX) * (x[iNb] - cX) + (y[iNb] - cY) * (y[iNb] - cY))) {
							if(clusterID[iNb] == 0) {
								clusterID[iNb] = cID;
								if(ind[iNb] == 0) {
									pointsToDo[nPToDo] = iNb;
									nPToDo ++;
									nBInCluster ++;
								}
								else {
									nEInCluster ++;
									coreCount ++;
								}
							}
							else if(clusterID[iNb] == -1 && nonCorePoints) {
								clusterID[iNb] = cID;
								if(ind[iNb] == 0) {
									nBInCluster ++;
								}
								else {
									nEInCluster ++;
								}
							}

							inCluster[iNb] = cID;
						}
					}

				}
			}
		
		}

		if(coreCount <= minCore)
		{
			for(int j = 0; j < count; j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
		else {
			double expEventInCluster = (double)(nBInCluster) / countB * countE;
			double LL = nEInCluster * log(nEInCluster/expEventInCluster);
			if(nEInCluster < countE) {
				LL += (countE - nEInCluster) * log((countE - nEInCluster) / (countE - expEventInCluster));
			}

			if(cID == 1) {
				*pCInfo = (struct clusterInfo *) malloc (sizeof (struct clusterInfo));
				curInfo = *pCInfo;
			}
			else {
				curInfo->next = (struct clusterInfo *) malloc (sizeof (struct clusterInfo));
				curInfo = curInfo->next; 
			}
			curInfo->clusterID = cID;
			curInfo->count1 = nEInCluster;
			curInfo->expCount1 = expEventInCluster;
			curInfo->ll = LL;
			curInfo->next = NULL; 

//			printf("%d,%d,%lf,%lf\n", cID, nEInCluster, expEventInCluster, LL);
		}
	}

	free(inCluster);

	free(pointsToDo);
	return clusterID; 
}


/**
 * NAME:	doClusterBer
 * DESCRIPTION:	cluster all event points based on a Binomial Test
 * PARAMETERS:
 * 	double * x: 		the array of points' X values
 * 	double * y: 		the array of points' Y values
 * 	int * ind:			the array of points' type indicator (1: case, 0: control)
 * 	int * index:		the index of all event points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int countCas:		the number of case points
 *	int countCon:		the number of control points
 *	int * casC:			the number of case points (within radius) near each case points
 *	int * conC:			the number of control points (within radius) near each case points
 *	double p:			the p of Possion distribution
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 *	struct clusterInfo ** pCInfo: the resulting output clusterInfo
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	the cluster ID of each point
 */
int * doClusterBer(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countCas, int countCon, int * casC, int * conC, double p, double significance, int minCore, bool nonCorePoints, struct clusterInfo ** pCInfo)
{
	int count = index[nBlockX * nBlockY];

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		if(BinomialTest(casC[i], conC[i], p) < significance)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	int * inCluster;

	if(NULL == (inCluster = (int *)malloc(sizeof(int) * count))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++) {
		inCluster[i] = -1;
	}

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;

	int coreCount;

	int nCasInCluster;
	int nConInCluster;
	struct clusterInfo * curInfo;
//	printf("ClusterID,nCas,nCon,LL\n");

	for(int i = 0; i < count; i++)
	{
		if(clusterID[i] != 0 || ind[i] == 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;

		inCluster[i] = cID;
		nCasInCluster = 1;
		nConInCluster = 0;

		while(nPToDo > 0) {
			nPToDo --;
			cX = x[pointsToDo[nPToDo]];		
			cY = y[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = index[row * nBlockX + colMin]; iNb < index[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(inCluster[iNb] != cID) {
						if(dist2 >= ((x[iNb] - cX) * (x[iNb] - cX) + (y[iNb] - cY) * (y[iNb] - cY))) {
							if(clusterID[iNb] == 0) {
								clusterID[iNb] = cID;

								if(ind[iNb] == 0) {
									nConInCluster ++;
								}
								else {
									pointsToDo[nPToDo] = iNb;
									nPToDo ++;
									nCasInCluster ++;
									coreCount ++;
								}
							}
							else if(clusterID[iNb] == -1 && nonCorePoints) {
								clusterID[iNb] = cID;
								if(ind[iNb] == 0) {
									nConInCluster ++;
								}
								else {
									nCasInCluster ++;
								}
							}

							inCluster[iNb] = cID;
							
						}
					}
				}
			}
		}

		if(coreCount <= minCore)
		{
			for(int j = 0; j < (count); j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
		else {
			double countInCl = nCasInCluster + nConInCluster;
			double LL = 0;
			if(nCasInCluster > 0) {
				LL += nCasInCluster * log(nCasInCluster/countInCl);
			}
			if(nConInCluster > 0) {
				LL += nConInCluster * log(nConInCluster/countInCl);
			}
			if(countCas > nCasInCluster) {
				LL += (countCas - nCasInCluster) * log((countCas - nCasInCluster)/(count-countInCl));
			}
			if(countCon > nConInCluster) {
				LL += (countCon - nConInCluster) * log((countCon - nConInCluster)/(count-countInCl));
			}
			if(cID == 1) {
				*pCInfo = (struct clusterInfo *) malloc (sizeof (struct clusterInfo));
				curInfo = *pCInfo;
			}
			else {
				curInfo->next = (struct clusterInfo *) malloc (sizeof (struct clusterInfo));
				curInfo = curInfo->next; 
			}
			curInfo->clusterID = cID;
			curInfo->count1 = nCasInCluster;
			curInfo->count0 = nConInCluster;
			curInfo->ll = LL;
			curInfo->next = NULL; 
//			printf("%d,%d,%d,%lf\n", cID, nCasInCluster, nConInCluster, LL);
		}
	}

	free(pointsToDo);
	free(inCluster);

	return clusterID; 
}


/**
 * NAME:	doClusterDBSCAN
 * DESCRIPTION:	cluster all event points using DBSCAN algorithm
 * PARAMETERS:
 * 	double * x: 		the array of events' X values
 * 	double * y: 		the array of events' Y values
 * 	int * index:		the index of all event points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	int minPts:		the minimum points to form a core points
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int * eC:		the number of event points (within radius) near each event points
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	an array of length count: the cluster ID of each case and control point
 */
int * doClusterDBSCAN(double * x, double * y, int * index, int nBlockX, int nBlockY, double radius, int minPts, double xMin, double yMin, int * eC, int minCore, bool nonCorePoints) {

	int count = index[nBlockX * nBlockY];

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		if(eC[i] >= minPts)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;

	int coreCount;

	for(int i = 0; i < count; i++)
	{
		if(clusterID[i] != 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;	

		while(nPToDo > 0) {
			nPToDo --;
			cX = x[pointsToDo[nPToDo]];		
			cY = y[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = index[row * nBlockX + colMin]; iNb < index[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(clusterID[iNb] < 1)
					{
						if(dist2 >= ((x[iNb] - cX) * (x[iNb] - cX) + (y[iNb] - cY) * (y[iNb] - cY)))
						{
							if(clusterID[iNb] != -1)
							{
								pointsToDo[nPToDo] = iNb;
								nPToDo ++;
								coreCount ++;
								clusterID[iNb] = cID;
							}
							else if(nonCorePoints)
								clusterID[iNb] = cID;
						}
					}
				}
			}
		
		}

		if(coreCount <= minCore)
		{
			for(int j = 0; j < count; j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
	}


	free(pointsToDo);
	return clusterID;
}

/**
 * NAME:	berMaximumLL
 * DESCRIPTION:	find the maximum log likelihood of any cluster in a Bernoulli model
 * PARAMETERS:
 * 	double * x: 		the array of points' X values
 * 	double * y: 		the array of points' Y values
 * 	int * ind:			the array of points' type indicator (1: case, 0: control)
 * 	int * index:		the index of all event points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 *	double radius:		the search radius, which is also the block size
 *	double xMin:		the minimum X of all points
 *	double yMin:		the minimum Y of all points
 *	int countCas:		the number of case points
 *	int countCon:		the number of control points
 *	int * casC:			the number of case points (within radius) near each case points
 *	int * conC:			the number of control points (within radius) near each case points
 *	double p:			the p of Possion distribution
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:		the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:	whether a cluster include non-core points
 * RETURN:
 * 	TYPE:	double 
 * 	VALUE:	the maximum log likelihood of any clusters
 */
double berMaximumLL(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countCas, int countCon, int * casC, int * conC, double p, double significance, int minCore, bool nonCorePoints)
{
	int count = index[nBlockX * nBlockY];

	double resultLL = 1;

	int * clusterID;
	if(NULL == (clusterID = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		if(BinomialTest(casC[i], conC[i], p) < significance)
			clusterID[i] = 0;
		else
			clusterID[i] = -1;
//		printf("%d,%d,%d\n", casC[i], conC[i], clusterID[i]);
	}

	int * pointsToDo;
	if(NULL == (pointsToDo = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int nPToDo = 0;
	int cID = 0;

	int * inCluster;

	if(NULL == (inCluster = (int *)malloc(sizeof(int) * count))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++) {
		inCluster[i] = -1;
	}

	double dist2 = radius * radius;

	double cX, cY;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	int iNb;

	int coreCount;

	int nCasInCluster;
	int nConInCluster;

	for(int i = 0; i < count; i++)
	{
		if(clusterID[i] != 0 || ind[i] == 0)
			continue;
		pointsToDo[0] = i;
		nPToDo = 1;
		cID ++;
		clusterID[i] = cID;
		
		coreCount = 1;

		inCluster[i] = cID;
		nCasInCluster = 1;
		nConInCluster = 0;

		while(nPToDo > 0) {
			nPToDo --;
			cX = x[pointsToDo[nPToDo]];		
			cY = y[pointsToDo[nPToDo]];

			colID = (int)((cX - xMin) / radius);
			rowID = (int)((cY - yMin) / radius);

			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);

			for(int row = rowMin; row <= rowMax; row ++)
			{
				for(iNb = index[row * nBlockX + colMin]; iNb < index[row * nBlockX + colMax + 1]; iNb ++)
				{
					if(inCluster[iNb] != cID) {
						if(dist2 >= ((x[iNb] - cX) * (x[iNb] - cX) + (y[iNb] - cY) * (y[iNb] - cY))) {
							if(clusterID[iNb] == 0) {
								clusterID[iNb] = cID;

								if(ind[iNb] == 0) {
									nConInCluster ++;
								}
								else {
									pointsToDo[nPToDo] = iNb;
									nPToDo ++;
									nCasInCluster ++;
									coreCount ++;
								}
							}
							else if(clusterID[iNb] == -1 && nonCorePoints) {
								clusterID[iNb] = cID;
								if(ind[iNb] == 0) {
									nConInCluster ++;
								}
								else {
									nCasInCluster ++;
								}
							}

							inCluster[iNb] = cID;
							
						}
					}
				}
			}
		}

		if(coreCount <= minCore)
		{
//			printf("Here");
			for(int j = 0; j < (count); j++)
			{
				if(clusterID[j] == cID)
					clusterID[j] = -1;
			}
			cID --;
		}
		else {
			double countInCl = nCasInCluster + nConInCluster;
			double LL = 0;
			if(nCasInCluster > 0) {
				LL += nCasInCluster * log(nCasInCluster/countInCl);
			}
			if(nConInCluster > 0) {
				LL += nConInCluster * log(nConInCluster/countInCl);
			}
			if(countCas > nCasInCluster) {
				LL += (countCas - nCasInCluster) * log((countCas - nCasInCluster)/(count-countInCl));
			}
			if(countCon > nConInCluster) {
				LL += (countCon - nConInCluster) * log((countCon - nConInCluster)/(count-countInCl));
			}
			if(resultLL > 0 || resultLL < LL) {
				resultLL = LL;
			}
		}
	}

	free(pointsToDo);
	free(inCluster);
	free(clusterID);

	return resultLL; 
}


