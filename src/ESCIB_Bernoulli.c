#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io.h"
#include "countPoints.h"
#include "clusters.h"

int main(int argc, char ** argv) {

	if(argc != 9) {
		printf("ERROR! Incorrect number of input arguments\n");
		printf("ESCIB_Bernoulli inputCase inputControl output searchRadius significance(alpha) baselineRatio minCorPointsInEachCluster nonCorePoints\n");
		return 1;
	}

	double xMin = 999999999, yMin = 999999999, xMax = -999999999, yMax = -999999999;

	FILE * inputCas;
	FILE * inputCon;
	FILE * output;

	double radius = atof(argv[4]);
	double significance = atof(argv[5]);

	double baseLineRatio = atof(argv[6]);
	double minCore = atof(argv[7]);
	bool nonCorePoints = true;
	if(atoi(argv[8]) == 0)
		nonCorePoints = false;

	if(NULL == (inputCas = fopen(argv[1], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}
	if(NULL == (inputCon = fopen(argv[2], "r")))
	{
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}

		
	int countCas = getCount(inputCas, xMin, xMax, yMin, yMax);
	int countCon = getCount(inputCon, xMin, xMax, yMin, yMax);
	int count = countCas + countCon;

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

	printf("Number of cases: %d\n", countCas);
	printf("Number of controls: %d\n", countCon);
	printf("X Range: %lf - %lf\n", xMin, xMax);
	printf("Y Range: %lf - %lf\n", yMin, yMax);

	readPoints(inputCas, x, y);
	readPoints(inputCon, x + countCas, y + countCas);

	for(int i = 0; i < countCas; i++) {
		ind[i] = 1;
	}
	for(int i = countCas; i < count; i++) {
		ind[i] = 0;
	}

	int nBlockX = ceil((xMax - xMin) / radius);
	int nBlockY = ceil((yMax - yMin) / radius);

//	printf("Index blocks: %d * %d\n", nBlockX, nBlockY);

	int * index;

	index = indexPoints(x, y, ind, count, xMin, yMin, nBlockX, nBlockY, radius);

	fclose(inputCas);
	fclose(inputCon);

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

	countInDistance(x, y, ind, index, nBlockX, nBlockY, radius, countPointsCon, countPointsCas);

	double p = baseLineRatio * countCas / (countCas + countCon); 

	struct clusterInfo * cInfo;

	int * clusters = doClusterBer(x, y, ind, index, nBlockX, nBlockY, radius, xMin, yMin, countCas, countCon, countPointsCas, countPointsCon, p, significance, minCore, nonCorePoints, &cInfo);
		//Output 
	if(NULL == (output = fopen(argv[3], "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}

	fprintf(output, "X,Y,CaseOrCon,ClusterID\n");
	for(int i = 0; i < count; i++) {
		fprintf(output, "%lf,%lf,%d,%d\n", x[i], y[i], ind[i], clusters[i]);
	}

	fclose(output);

	char * outputCInfo = (char *) malloc((strlen(argv[3]) + 10) * sizeof(char));
	outputCInfo[0] = '\0';
	strcat(outputCInfo, argv[3]);
	strcat(outputCInfo, "_Info");

	if(NULL == (output = fopen(outputCInfo, "w"))) {
		printf("ERROR: Can't open the output file.\n");
		exit(1);
	}

	fprintf(output, "ClusterID,nCas,nCon,LL\n");
	struct clusterInfo * curInfo = cInfo;
	while(curInfo != NULL) {
		cInfo = curInfo->next;
		fprintf(output, "%d,%d,%d,%lf\n", curInfo->clusterID, curInfo->count1, curInfo->count0, curInfo->ll);
		
		free(curInfo);
		curInfo = cInfo;		
	}
	

	fclose(output);

	free(outputCInfo);	



	free(x);
	free(y);
	free(ind);
	free(index);
	free(countPointsCas);
	free(countPointsCon);


	free(clusters);

	return 0;
}
