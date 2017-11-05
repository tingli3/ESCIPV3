#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io.h"
#include "countPoints.h"
#include "clusters.h"
#include "mc.h"

int main(int argc, char ** argv) {

	if(argc != 10) {
		printf("ERROR! Incorrect number of input arguments\n");
		printf("ESCIB_Poisson inputBackground inputEvents output searchRadius significance(alpha) baselineRatio minCorPointsInEachCluster nonCorePoints nSim\n");
		return 1;
	}

	double xMin = 999999999, yMin = 999999999, xMax = -999999999, yMax = -999999999;

	FILE * inputB;
	FILE * inputE;
	FILE * output;

	double radius = atof(argv[4]);
	double significance = atof(argv[5]);

	double baseLineRatio = atof(argv[6]);
	double minCore = atof(argv[7]);
	bool nonCorePoints = true;
	if(atoi(argv[8]) == 0)
		nonCorePoints = false;
	int nSim = atoi(argv[9]);

	if(NULL == (inputB = fopen(argv[1], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}
	if(NULL == (inputE = fopen(argv[2], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}

		
	int countB = getCount(inputB, xMin, xMax, yMin, yMax);
	int countE = getCount(inputE, xMin, xMax, yMin, yMax);
	int count = countE + countB;

	double * x;
	double * y;
	int * ind;

	if(NULL == (x = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (y = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (ind = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	printf("Number of background points: %d\n", countB);
	printf("Number of event points: %d\n", countE);
	printf("X Range: %lf - %lf\n", xMin, xMax);
	printf("Y Range: %lf - %lf\n", yMin, yMax);
	printf("Search radius %lf\n", radius);
	
	int nBlockX = ceil((xMax - xMin) / radius);
	int nBlockY = ceil((yMax - yMin) / radius);

	readPoints(inputB, x, y);
	readPoints(inputE, x + countB, y + countB);

	for(int i = 0; i < countB; i++) {
		ind[i] = 0;
	}
	for(int i = countB; i < count; i++) {
		ind[i] = 1;
	}

	double * xB;
	double * yB;
	int * indexB;
	
	if(nSim > 0) {
		//Point index for MC

		if(NULL == (xB = (double *)malloc(sizeof(double) * countB))) {
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}
		if(NULL == (yB = (double *)malloc(sizeof(double) * countB))) {
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}

		for(int i = 0; i < countB; i++) {
			xB[i] = x[i];
			yB[i] = y[i];
		}

		indexB = indexPoints(xB, yB, countB, xMin, yMin, nBlockX, nBlockY, radius);
	}



//	printf("Index blocks: %d * %d\n", nBlockX, nBlockY);

	int * index;


	index = indexPoints(x, y, ind, count, xMin, yMin, nBlockX, nBlockY, radius);

	fclose(inputB);
	fclose(inputE);

	int * countPointsE;
	int * countPointsB;

	if(NULL == (countPointsE = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (countPointsB = (int *)malloc(sizeof(int) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	countInDistance(x, y, ind, index, nBlockX, nBlockY, radius, countPointsB, countPointsE);

	double * lambda;
	if(NULL == (lambda = (double *)malloc(sizeof(double) * count)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < count; i++)
	{
		//lambda[i] = (double)(countPointsB[i]) * countE / countB;
		lambda[i] = (double)(countPointsB[i]) * countE * baseLineRatio / countB;
	}
	
	struct clusterInfo * cInfo;

	int * clusters = doClusterPoi(x, y, ind, index, nBlockX, nBlockY, radius, xMin, yMin, countB, countE, countPointsE, lambda, significance, minCore, nonCorePoints, &cInfo);

	//Output 
	if(NULL == (output = fopen(argv[3], "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}


	for(int i = 0; i < count; i++) {
		if(ind[i] == 1) {
			fprintf(output, "%lf,%lf,%d\n", x[i], y[i], clusters[i]);
		}
	}

	fclose(output);
	free(lambda);
	free(x);
	free(y);
	free(ind);
	free(index);
	free(countPointsE);
	free(countPointsB);
	free(clusters);

	if(nSim > 0) {
		//MC
		monteCarloPoi(xB, yB, indexB, nBlockX, nBlockY, radius, xMin, yMin, countE, countB, baseLineRatio, significance, minCore, nonCorePoints, nSim, cInfo);


		free(xB);
		free(yB);
		free(indexB);
	}


	char * outputCInfo = (char *) malloc((strlen(argv[3]) + 10) * sizeof(char));
	outputCInfo[0] = '\0';
	strcat(outputCInfo, argv[3]);
	strcat(outputCInfo, "_Info");

	if(NULL == (output = fopen(outputCInfo, "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}

	if(nSim > 0) {
		fprintf(output, "ClusterID,Events,expEvents,LL,PValue\n");
	}
	else {
		fprintf(output, "ClusterID,Events,expEvents,LL\n");
	}
	struct clusterInfo * curInfo = cInfo;
	while(curInfo != NULL) {
		cInfo = curInfo->next;
		if(nSim > 0) {
			fprintf(output, "%d,%d,%lf,%lf,%lf\n", curInfo->clusterID, curInfo->count1, curInfo->expCount1, curInfo->ll, curInfo->pValue);
		}
		else {
			fprintf(output, "%d,%d,%lf,%lf\n", curInfo->clusterID, curInfo->count1, curInfo->expCount1, curInfo->ll);
		}
		free(curInfo);
		curInfo = cInfo;		
	}
	

	fclose(output);
	free(outputCInfo);	
	
	return 0;
}
