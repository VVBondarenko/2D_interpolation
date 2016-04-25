#include <stdio.h>
#include <stdlib.h>

void init1DArr	(double **array, int size)
{
	*array = malloc(size * sizeof(double));
}

void free1DArr	(double **array)
{
	free((void *)(*array));
}


void init2DArr	(double ***array, int sizeX, int sizeY)
{
	int i,j;
	*array = malloc(sizeX * sizeof(double*));
	for(i=0; i<sizeX; i++)
		(*array)[i] = malloc(sizeY * sizeof(double));
		
	for(i = 0; i < sizeX; i++)
		for(j = 0; j < sizeY; j++)
			(*array)[i][j] = 0;
}

void free2DArr	(double ***array, int sizeX)
{
	int i;
	for(i=0; i<sizeX; i++)
		free((void *)(*array)[i]);
	free((void *)(*array));
}

