#ifndef CPH
#define CPH

void countInDistance(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double distance, int * count0, int * count1);
int * countInDistance_Single(double * xE, double * yE, int * indexE, int nBlockX, int nBlockY, double distance);
int * countInDistance_Double(double * xE, double * yE, double * xB, double * yB, int * indexE, int * indexB, int nBlockX, int nBlockY, double distance);
void countInDistance_EventsInPop(double * xB, double * yB, int * ind, int * indexB, int nBlockX, int nBlockY, double distance, int * countPointsE);

#endif
