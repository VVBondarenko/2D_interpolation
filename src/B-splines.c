#include <math.h>

double f_B_3(double X)
{
	double absV=fabs(X);
	if(absV<2)
	{
		if(absV>=1)
			return 0.25*(2.0-absV)*(2.0-absV)*(2.0-absV); 
		else
			return 1.0 - 1.5*absV*absV*(1.0 - 0.5*absV);
	}
	return 0.0;
}

double f_d_B_3(double X)
{
	double absV=fabs(X);
	if(X>=0)
	{
		if(absV<2)
		{
			if(absV>=1)
				return -0.75*(2.0-absV)*(2.0-absV); 
			else
				return -1.5*(2.0*absV - 1.5*absV*absV);
		}
	}
	else
	{
		if(absV<2)
		{
			if(absV>=1)
				return 0.75*(2.0-absV)*(2.0-absV); 
			else
				return 1.5*(2.0*absV - 1.5*absV*absV);
		}
	}
	return 0.0;
}

double f_dd_B_3(double X)
{
	double absV=fabs(X);
	if(absV<2)
	{
		if(absV>=1)
			return 1.5*(2.0-absV); 
		else
			return -1.5*(2.0 - 3.0*absV);
	}
	return 0.0;
}
