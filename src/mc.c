#include <stdlib.h>
#include <stdio.h>
#include <random>
#include "clusters.h"
#include "countPoints.h"

using namespace std;
/**
 * NAME:	simBerCase
 * DESCRIPTION:	simulate cases for Monte Carlo Simulation based on Bernoulli Model
 * PARAMETERS:
 * 	int * ind:			the array of points' type indicator (1: case, 0: control), will be randomly shuffled in the simulation
 *	int countCas:		the number of case points
 *	int count:			the number of all points
 */
void simBerCase(int * ind, int countCas, int count) {

	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::uniform_int_distribution<int> uni(0, count - 1);	

	for(int i = 0; i < count; i++) {
		ind[i] = 0;
	}

	int casID;

	for(int i = 0; i < countCas; i++) {
		casID = uni(rng);
		while(ind[casID] == 1)
			casID = uni(rng);
		ind[casID] = 1;
	}

	return;
}

/**
 * NAME:	monteCarloBer
 * DESCRIPTION:	calculate the P-Value of each cluster in a Bernoulli model
 * PARAMETERS:
 * 	double * x: 			the array of points' X values
 * 	double * y: 			the array of points' Y values
 * 	int * ind:				the array of points' type indicator (1: case, 0: control), will be randomly shuffled in the simulation
 * 	int * index:			the index of all event points
 * 	int nBlockX:			the number of index blocks along X dimension
 * 	int nBlockY:			the number of index blocks along Y dimension
 *	double radius:			the search radius, which is also the block size
 *	double xMin:			the minimum X of all points
 *	double yMin:			the minimum Y of all points
 *	int countCas:			the number of case points
 *	int countCon:			the number of control points
 *	double p:				the p of Possion distribution
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:			the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:		whether a cluster include non-core points
 *	int nSim:				the number of simulation to be conducted
 *	struct clusterInfo * cInfo:		the info of detected clusters, resulting p-values will be written to it
 */

void monteCarloBer(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countCas, int countCon, double p, double significance, int minCore, bool nonCorePoints, int nSim, struct clusterInfo * cInfo) {

	int nClusters = 0;
	int count = countCas + countCon;
	struct clusterInfo * curInfo = cInfo;
	while (curInfo!=NULL) {
		nClusters ++;
		curInfo = curInfo->next;
	}

	int llAbove[nClusters];
	double cLL[nClusters];

	curInfo = cInfo;
	for(int i = 0; i < nClusters; i++) {
		llAbove[i] = 0;
		cLL[i] = curInfo->ll;
		curInfo = curInfo->next;
	}
		
	int * countPointsCas;
	int * countPointsCon;

	if(NULL == (countPointsCas = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (countPointsCon = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	double simMaxLL;	

	for(int i = 0; i < nSim; i++) {
		//SimulateCases
		simBerCase(ind, countCas, count);
	
		//CalcCount
		countInDistance(x, y, ind, index, nBlockX, nBlockY, radius, countPointsCon, countPointsCas);

		//GetMaxLL
		simMaxLL = berMaximumLL(x, y, ind, index, nBlockX, nBlockY, radius, xMin, yMin, countCas, countCon, countPointsCas, countPointsCon, p, significance, minCore, nonCorePoints);

		printf("Simulation: %d\tLL: %lf\n", i, simMaxLL);
		//CompareLL
		if(simMaxLL<0) {
			for(int j = 0; j < nClusters; j++) {
				if(cLL[j] <= simMaxLL) {
					llAbove[j] ++;
				}
			}
		}
	}

	free(countPointsCas);
	free(countPointsCon);

	curInfo = cInfo;
	for(int i = 0; i < nClusters; i++) {
		curInfo->pValue = (double)(1+llAbove[i]) / (1+nSim);
		curInfo = curInfo->next;
	}

}

/**
 * NAME:	monteCarloPoi
 * DESCRIPTION:	calculate the P-Value of each cluster in a Poisson model
 * PARAMETERS:
 * 	double * xB: 			the array of background points' X values
 * 	double * yB: 			the array of background points' Y values
 * 	int * indexB:			the index of all background points
 * 	int nBlockX:			the number of index blocks along X dimension
 * 	int nBlockY:			the number of index blocks along Y dimension
 *	double radius:			the search radius, which is also the block size
 *	double xMin:			the minimum X of all points
 *	double yMin:			the minimum Y of all points
 *	int countE:				the number of event points
 *	int countB:				the number of background points
 *	double baseLineRatio:	the ratio null hypothesis to complete randomness baseline 1 means the same as baseline, 2 means twice the baseline
 *	double significance: 	the significane level to tell a cluste core point
 *	int minCore:			the minimum number of core points in each cluster (each cluste should have more core points than minCore)
 *	bool nonCorePoints:		whether a cluster include non-core points
 *	int nSim:				the number of simulation to be conducted
 *	struct clusterInfo * cInfo:		the info of detected clusters, resulting p-values will be written to it
 */

void monteCarloPoi(double * xB, double * yB, int * indexB, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countE, int countB, double baseLineRatio, double significance, int minCore, bool nonCorePoints, int nSim, struct clusterInfo * cInfo) {

	int nClusters = 0;
	struct clusterInfo * curInfo = cInfo;
	while (curInfo!=NULL) {
		nClusters ++;
		curInfo = curInfo->next;
	}

	int llAbove[nClusters];
	double cLL[nClusters];

	curInfo = cInfo;
	for(int i = 0; i < nClusters; i++) {
		llAbove[i] = 0;
		cLL[i] = curInfo->ll;
		curInfo = curInfo->next;
	}

	int * ind;
	if(NULL == (ind = (int *)malloc(sizeof(int) * countB))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int * countPointsB = countInDistance_Single(xB, yB, indexB, nBlockX, nBlockY, radius);
	int * countPointsE;
	if(NULL == (countPointsE = (int *)malloc(sizeof(int) * countB))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	double * lambda;
	if(NULL == (lambda = (double *)malloc(sizeof(double) * countB))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	for(int i = 0; i < countB; i++) {
		lambda[i] = (double)(countPointsB[i]) * countE * baseLineRatio / countB;
	}

	double simMaxLL;

	for(int i = 0; i < nSim; i++) {
		//Simulate case
		simBerCase(ind, countE, countB);

		//CountEvent
		countInDistance_EventsInPop(xB, yB, ind, indexB, nBlockX, nBlockY, radius, countPointsE);	
//		for(int j = 0; j < countB; j++) {
//			printf("%d,%lf\n", countPointsE[j], lambda[j]);
//		}

		//GetTopLikelihood
		simMaxLL = poiMaximumLL(xB, yB, ind, indexB, nBlockX, nBlockY, radius, xMin, yMin, countB, countE, countPointsE, lambda, significance, minCore, nonCorePoints);

		//Compare and update
		printf("Simulation: %d\tLL: %lf\n", i, simMaxLL);
		for(int j = 0; j < nClusters; j++) {
			if(cLL[j] <= simMaxLL) {
				llAbove[j] ++;
			}
		}
	}

	free(ind);
	free(countPointsB);
	free(countPointsE);
	free(lambda);

	curInfo = cInfo;
	for(int i = 0; i < nClusters; i++) {
		curInfo->pValue = (double)(1+llAbove[i]) / (1+nSim);
		curInfo = curInfo->next;
	}
}
