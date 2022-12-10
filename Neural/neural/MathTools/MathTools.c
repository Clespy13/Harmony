#include "MathTools.h"
#include <stdio.h>

//return the image of x by the sigmoid function
double sigmoid(double x)
{
	return 1.0f / (1.0f + exp(-x));
}

double dSigmoid(double x)
{
	return x * (1.0f - x);
}

double ReLu(double x)
{
	if (x < 0)
		return 0;
	else
		return x;
}

//returns a random double between 0.0 and 1.0
double init_weight()
{
	double tmp = ((double) rand() / ((double) RAND_MAX));
	printf("random = %f ", tmp);
	return tmp;
}

//initialize random seed with current time as seed
void init_random()
{
	time_t t;
	srand((unsigned)time(&t));
}

double get_random(double x)
{
	x = (double)rand() / (double)RAND_MAX;
	return x;
}