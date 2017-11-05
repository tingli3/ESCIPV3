#ifndef CH
#define CH

struct clusterInfo {
	int clusterID;
	int count0;
	int count1;
	double expCount1;
	double ll;
	double pValue;
	struct clusterInfo * next;
};

//Poisson
int * doClusterPoi(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countB, int countE, int * eC, double * lambda, double significance, int minCore, bool nonCorePoints, struct clusterInfo ** pCInfo);
double poiMaximumLL(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countB, int countE, int * eC, double * lambda, double significance, int minCore, bool nonCorePoints);
//Bernoulli
int * doClusterBer(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countCas, int countCon, int * casC, int * conC, double p, double significance, int minCore, bool nonCorePoints, struct clusterInfo ** pCInfo);
double berMaximumLL(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countCas, int countCon, int * casC, int * conC, double p, double significance, int minCore, bool nonCorePoints);
//DBSCAN
int * doClusterDBSCAN(double * x, double * y, int * index, int nBlockX, int nBlockY, double radius, int minPts, double xMin, double yMin, int * eC, int minCore, bool nonCorePoints);

#endif
