#ifndef MCH
#define MCH
void monteCarloBer(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countCas, int countCon, double p, double significance, int minCore, bool nonCorePoints, int nSim, struct clusterInfo * cInfo);
void monteCarloPoi(double * xB, double * yB, int * indexB, int nBlockX, int nBlockY, double radius, double xMin, double yMin, int countE, int countB, double baseLineRatio, double significance, int minCore, bool nonCorePoints, int nSim, struct clusterInfo * cInfo);

#endif
