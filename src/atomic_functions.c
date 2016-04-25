/*
 *	Atomic Functions set
 *	version 1.0 
 *	
 *	Functions: up(x), cup(x)
 *
 *	April, 2016
 */

#include <math.h>
#include <stdint.h>
#define STEPS 12


double  x	[4096]; //values of x where up(x) is avaliable
double upx	[4096]; //function  up (x)
double d_up	[4096]; //function up' (x) (is it req.?)
double ddup	[4096]; //function up''(x) (is it req.?)

double 	cup [8192];	//function cup(x) values
double dcup [8192];


void init_up(){	//calculates up(x), up`(x), up``(x)
	int  m=0,n=0, k,line, size_l=pow(2,STEPS);
	long int  data[size_l][2];
	for(m=0;m<size_l;m++)
	{
		data[m][0] = 0.0;
		data[m][1] = 0.0;
	}
	data [0][0] = 1.0;
	
	data [0][1] = 1.0;
	data [1][0] = 3.0;
	data [2][0] = 5.0;
	data [3][0] = 7.0;
	data [4][0] = 8.0;
	data [5][0] = 8.0;
	data [6][0] = 8.0;
	data [7][0] = 8.0;
	data [8][0] = 7.0;
	data [9][0] = 5.0;
	data [10][0] = 3.0;
	data [11][0] = 1.0;

	for(m=0;m<16;m++){
		data[m+16][0] = -data[m][0];
	}

	for(m=5;m<STEPS;m++){
		 line = m%2;
		 k = pow(2,m);
		data[1][line] = (float)m;

		for(n=1;n<k;n++)
			data[n][line] = data[n-1][line] + data[n][(m-1)%2];
		for(n=0;n<k;n++)
			data[k+n][line] = /*data[k][line]*/ - data[n][line];	
	}

	for(n=1;n<size_l;n++)
		data[n][line] = data[n-1][line] + data[n][(m-1)%2];
		
	for(n=0;n<size_l;n++)
	{
		upx[n] = ((double)(data[n][line]))/((double)(data[size_l/2-1][line]));
		ddup[n]= 0.0;
		x[n] = -1.0+2.0/4096.0*(n+1);
	}
	for(n=0;n<size_l/2;n++)
	{
			d_up[n] = 2.0*upx[2*n];
			d_up[n+size_l/2] = -2.0*upx[2*n];
			if(n<size_l/4)
			{
				ddup[n] = 			 4.0*upx[4*n];
				ddup[n+size_l/4] = 	-4.0*upx[4*n];	
				ddup[n+size_l/2] = 	-4.0*upx[4*n];
				ddup[n+3*size_l/4] = 4.0*upx[4*n];
			}
	}
}
//returns value of function up(x) in particular point, if its possible
double f_up(double X){ 
	int i;
	if(X>-1.0 && X<1.0)
	{
		//i=((int)((X+1.0)/2.0*4084.0))-1;
		i=((int)((X+1.0)*2048.0))-1;
		return upx[i];
	}
	else
	{
		return 0.0;
	}
	
}

double interpolate(double X, double Values[][2], double step, int N)
{
	int i;
	if(X>Values[0][0] && X<Values[N-1][0])
		for(i=0;i<N-1;i++)
		{
			if(Values[i][0]<=X && Values[i+1][0]>X)
				return Values[i][1]*f_up((X-Values[i][0])/step) + Values[i+1][1]*f_up((X-Values[i+1][0])/step);
		}
	return 0.0;
}

double interpolateCIP(double X, double Values[][5], double step, int N)
{
	int i;
	double toRet=0.0;
	if(X>Values[0][0] && X<Values[N-1][0])
		for(i=0;i<N-1;i++){
			if(Values[i][0]<=X && Values[i+1][0]>X)
				toRet =  f_up((X-Values[i][0])/step)*(Values[i][1] + (X-Values[i][0])*(Values[i][2]+0.5*(X-Values[i][0])*(Values[i][3]+0.333333333333333333*(X-Values[i][0])*Values[i][4])));
				toRet += f_up((X-Values[i+1][0])/step)*(Values[i+1][1] + (X-Values[i+1][0])*(Values[i+1][2]+0.5*(X-Values[i+1][0])*(Values[i+1][3]+0.333333333333333333*(X-Values[i+1][0])*Values[i+1][4])));
		}
	return toRet;
}

double interpolateDD(double X, double Values[][4], double step, int N)
{
	int i;
	double toRet=0.0;
	//toRet  = 0.5 * (/*Values[0][1]*f_up(0.5*(X-Values[0][0]-2.0)/step) + */Values[0][1]*f_up(0.5*(X-Values[0][0]-1.0)/step) +  Values[0][1]*f_up(0.5*(X-Values[0][0])/step)+Values[1][1]*f_up(0.5*(X-Values[1][0])/step)+Values[2][1]*f_up(0.5*(X-Values[2][0])/step));
	//toRet = Values[0][1]*f_up((X-Values[0][0]-2.0)/step);
	if(X>=Values[0][0] && X<Values[1][0])
		toRet = Values[0][1]+(X-Values[0][0])*Values[0][2]+(X-Values[0][0])*(X-Values[0][0])*0.5*Values[0][3];

	if(X>=Values[N-2][0] && X<=Values[N-1][0])
		toRet = Values[N-2][1]+(X-Values[N-2][0])*Values[N-2][2]+(X-Values[N-2][0])*(X-Values[N-2][0])*0.5*Values[N-2][3];

	if(X>=Values[1][0] && X<Values[N-2][0])
		for(i=1;i<N-1;i++)
		{
				if(Values[i][0]<=X && Values[i+1][0]>X && i!=0)
					toRet  = 0.5 * (Values[i-2][1]*f_up(0.5*(X-Values[i-2][0])/step) + Values[i-1][1]*f_up(0.5*(X-Values[i-1][0])/step) +  Values[i][1]*f_up(0.5*(X-Values[i][0])/step)+Values[i+1][1]*f_up(0.5*(X-Values[i+1][0])/step)+Values[i+2][1]*f_up(0.5*(X-Values[i+2][0])/step));
				//toRet += 0.5 * ((X-Values[i-2][0])*Values[i-2][2]*f_up(0.5*(X-Values[i-2][0])/step) + (X-Values[i-1][0])*Values[i-1][2]*f_up(0.5*(X-Values[i-1][0])/step) + (X-Values[i][0])*Values[i][2]*f_up(0.5*(X-Values[i][0])/step)+(X-Values[i+1][0])*Values[i+1][2]*f_up(0.5*(X-Values[i+1][0])/step) +(X-Values[i+2][0])*Values[i+2][2]*f_up(0.5*(X-Values[i+2][0])/step));
		}
	//if(X> 0 && X<Values[1][0])
	//	toRet  = 0.5 * (/*Values[0][1]*f_up(0.5*(X-Values[0][0]-2.0)/step) + */Values[0][1]*f_up(0.5*(X-Values[0][0]-1.0)/step) +  Values[0][1]*f_up(0.5*(X-Values[0][0])/step)+Values[1][1]*f_up(0.5*(X-Values[1][0])/step)+Values[2][1]*f_up(0.5*(X-Values[2][0])/step));

	return toRet;
}

double interpolateDDd(double X, double Values[][4], double step, int N)
{
	int i;
	double toRet=0.0;
	//toRet  = 0.5 * (/*Values[0][1]*f_up(0.5*(X-Values[0][0]-2.0)/step) + */Values[0][1]*f_up(0.5*(X-Values[0][0]-1.0)/step) +  Values[0][1]*f_up(0.5*(X-Values[0][0])/step)+Values[1][1]*f_up(0.5*(X-Values[1][0])/step)+Values[2][1]*f_up(0.5*(X-Values[2][0])/step));
	//toRet = Values[0][1]*f_up((X-Values[0][0]-2.0)/step);
	if(X>=Values[0][0] && X<Values[1][0])
		toRet = Values[0][1]+(X-Values[0][0])*Values[0][2]+(X-Values[0][0])*(X-Values[0][0])*0.5*Values[0][3];

	if(X>=Values[N-2][0] && X<=Values[N-1][0])
		toRet = Values[N-2][1]+(X-Values[N-2][0])*Values[N-2][2]+(X-Values[N-2][0])*(X-Values[N-2][0])*0.5*Values[N-2][3];

	if(X>=Values[1][0] && X<Values[N-2][0])
		for(i=1;i<N-1;i++)
		{
				if(Values[i][0]<=X && Values[i+1][0]>X && i!=0)
					toRet  = 0.5 * (Values[i-2][1]*f_up(0.5*(X-Values[i-2][0])/step) + Values[i-1][1]*f_up(0.5*(X-Values[i-1][0])/step) +  Values[i][1]*f_up(0.5*(X-Values[i][0])/step)+Values[i+1][1]*f_up(0.5*(X-Values[i+1][0])/step)+Values[i+2][1]*f_up(0.5*(X-Values[i+2][0])/step));
					//toRet -= 0.25 * (X-Values[i][0])*(Values[i-1][2]*f_up((X-Values[i-1][0])/step) + Values[i][2]*f_up((X-Values[i][0])/step) + Values[i+1][2]*f_up((X-Values[i+1][0])/step));
					toRet -= 0.25 * (Values[i-1][2]*(X-Values[i-1][0])*f_up((X-Values[i-1][0])/step) + Values[i][2]*(X-Values[i][0])*f_up((X-Values[i][0])/step) + Values[i+1][2]*(X-Values[i+1][0])*f_up((X-Values[i+1][0])/step));
		}
	//if(X> 0 && X<Values[1][0])
	//	toRet  = 0.5 * (/*Values[0][1]*f_up(0.5*(X-Values[0][0]-2.0)/step) + */Values[0][1]*f_up(0.5*(X-Values[0][0]-1.0)/step) +  Values[0][1]*f_up(0.5*(X-Values[0][0])/step)+Values[1][1]*f_up(0.5*(X-Values[1][0])/step)+Values[2][1]*f_up(0.5*(X-Values[2][0])/step));

	return toRet;
}

//
//
//


void init_cup()
{
    //double * result = new double[N + M - 1];
    //memset(result, 0, sizeof(double) * (N + M - 1));
	const int N=4096;
	int i, j;
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
			 cup[i + j] += upx[i] * upx[N-1-j];
            dcup[i + j] += d_up[i] * upx[N-1-j];
        }
    }

}

double f_cup(double X)
{
	int i;
	if(X>-2.0 && X<2.0)
	{
		//i=((int)((X+1.0)/2.0*4084.0))-1;
		i=((int)((X+2.0)*2048.0))-1;
		//i=((int)((X+1.0)*4085.0))-1;
		//return cup[i]/cup[4085];
		return cup[i]/cup[4096];
	}
	else
	{
		return 0.0;
	}
}

double f_d_cup(double X)
{
	//return (f_cup(X+0.02)-f_cup(X-0.02))/0.04;
	int i;
	if(X>-2.0 && X<2.0)
	{
		//i=((int)((X+1.0)/2.0*4084.0))-1;
		i=((int)((X+2.0)*2048.0))-1;
		return 2.0*dcup[i]/cup[4096];
	}
	else
	{
		return 0.0;
	}
}

double f_dd_cup1(double X)
{
	return (f_d_cup(X+0.03)-f_d_cup(X-0.03))/0.06;
}
double f_dd_cup2(double X)
{
	return (f_cup(X+0.03)-2.0*f_cup(X)+f_cup(X-0.03))/0.009;
}

