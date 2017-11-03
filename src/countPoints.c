/**
 * countPoints.c
 * Author: Ting Li <tingli3@illinois.edu>
 * Date: 08/07/2017
 */ 


#include <stdio.h>
#include <stdlib.h>

/**
 * NAME:	countInDistance
 * DESCRIPTION:	get the number of type 2 type of points within a distance of each point
 * PARAMETERS:
 * 	double * x:			points' X values 
 * 	double * y:			points' Y values 
 * 	int * ind:			points' type indicator
 * 	int * index:		the index of type A points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 * 	double distance:	the distance, which is also the size (side length) of each index block
 * 	int * count0:		the output array of the numbers of first type of points within the distance , ordered the same as x and y
 * 	int * count1:		the output array of the numbers of first type of points within the distance , ordered the same as x and y
 */
void countInDistance(double * x, double * y, int * ind, int * index, int nBlockX, int nBlockY, double distance, int * count0, int * count1)
{
	int count = index[nBlockX * nBlockY];
	double xi, yi;
	double dist2 = distance * distance;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;

	for(rowID = 0; rowID < nBlockY; rowID ++)
	{
		for(colID = 0; colID < nBlockX; colID ++)
		{
			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);
			for(int i = index[rowID * nBlockX + colID]; i < index[rowID * nBlockX + colID + 1]; i++) {	
				xi = x[i];
				yi = y[i];
				count0[i] = 0;
				count1[i] = 0;
				for(int row = rowMin; row <= rowMax; row ++)
				{
					for(int j = index[row * nBlockX + colMin]; j < index[row * nBlockX + colMax + 1]; j ++)
					{
						if(dist2 >= ((x[j] - xi) * (x[j] - xi) + (y[j] - yi) * (y[j] - yi)))
						{
							if(ind[j] == 0)
							{
								count0[i] ++;
							}
							else
							{
								count1[i] ++;
							}
						}
					}	
				}
			}	
		}
	}
}


/**
 * NAME:	countInDistance_Single
 * DESCRIPTION:	get the number of type A points within a distance of each type A point
 * PARAMETERS:
 * 	double * xE:		type A points' X values 
 * 	double * yE:		type A points' Y values 
 * 	int * indexE:		the index of type A points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 * 	double distance:	the distance, which is also the size (side length) of each index block
 * RETURN:
 * 	TYPE:	int * 
 * 	VALUE:	an array of the numbers of points within the distance, ordered the same as xE and yE
 */
int * countInDistance_Single(double * xE, double * yE, int * indexE, int nBlockX, int nBlockY, double distance)
{
	int countE = indexE[nBlockX * nBlockY];
	int * count;
	
	if(NULL == (count = (int *)malloc(sizeof(int) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	double x, y;
	double dis2 = distance * distance;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;
	int iC, iP;
	int pCol, pRow;
	for(rowID = 0; rowID < nBlockY; rowID ++)
	{
		for(colID = 0; colID < nBlockX; colID ++)
		{
			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);
			for(iC = indexE[rowID * nBlockX + colID]; iC < indexE[rowID * nBlockX + colID + 1]; iC++)
			{
				x = xE[iC];
				y = yE[iC];
				count[iC] = 0;
				for(int row = rowMin; row <= rowMax; row ++)
				{
					for(iP = indexE[row * nBlockX + colMin]; iP < indexE[row * nBlockX + colMax + 1]; iP ++)
					{
						if(dis2 >= ((xE[iP] - x) * (xE[iP] - x) + (yE[iP] - y) * (yE[iP] - y)))
							count[iC] ++;
					}

				}

			}
		}
	}
	return count;
}

/**
 * NAME:	countInDistance_Double
 * DESCRIPTION:	get the number of type B points within a distance of each type A point
 * PARAMETERS:
 * 	double * xE:		type A points' X values 
 * 	double * yE:		type A points' Y values 
 * 	double * xB:		type B points' X values 
 * 	double * yB:		type B points' Y values 
 * 	int * indexE:		the index of type A points
 * 	int * indexB:		the index of type B points
 * 	int nBlockX:		the number of index blocks along X dimension
 * 	int nBlockY:		the number of index blocks along Y dimension
 * 	double distance:	the distance, which is also the size (side length) of each index block
 * RETURN:
 * 	TYPE:	int * 
 * 	VALUE:	an array of the numbers of points within the distance
 */

int * countInDistance_Double(double * xE, double * yE, double * xB, double * yB, int * indexE, int * indexB, int nBlockX, int nBlockY, double distance)
{
	int countE = indexE[nBlockX * nBlockY];
	int countB = indexB[nBlockX * nBlockY];

	int * count;
	
	if(NULL == (count = (int *)malloc(sizeof(int) * countE)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	double x, y;
	double dis2 = distance * distance;
	int colID, rowID;
	int colMin, colMax, rowMin, rowMax;
	int iC, iP;
	int pCol, pRow;
	for(rowID = 0; rowID < nBlockY; rowID ++)
	{
		for(colID = 0; colID < nBlockX; colID ++)
		{
			colMin = (colID == 0) ? 0 : (colID - 1);
			colMax = (colID == nBlockX - 1) ? (nBlockX - 1) : (colID + 1);
			rowMin = (rowID == 0) ? 0 : (rowID - 1);
			rowMax = (rowID == nBlockY - 1) ? (nBlockY - 1) : (rowID + 1);
			for(iC = indexE[rowID * nBlockX + colID]; iC < indexE[rowID * nBlockX + colID + 1]; iC++)
			{
				x = xE[iC];
				y = yE[iC];
				count[iC] = 0;
				for(int row = rowMin; row <= rowMax; row ++)
				{
					for(iP = indexB[row * nBlockX + colMin]; iP < indexB[row * nBlockX + colMax + 1]; iP ++)
					{
						if(dis2 >= ((xB[iP] - x) * (xB[iP] - x) + (yB[iP] - y) * (yB[iP] - y)))
							count[iC] ++;
					}

				}

			}
		}
	}
	return count;
}

