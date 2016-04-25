#include <stdio.h>
#include <stdlib.h>
#include "B-splines.c"
#include <gsl/gsl_linalg.h>
#include "atomic_functions.c"


double phi(double x, double y)
{
	return f_cup(x)*f_cup(y);
}

void insertA(gsl_matrix * sys, int i0, int j0, int n)
{
	int i;
	const double a=f_B_3(0.0),b=f_B_3(1.0);
	for(i=0;i<n;i++)
	{
		if(i==0)
		{
			gsl_matrix_set(sys, i0+0,j0+0,a);
			gsl_matrix_set(sys, i0+0,j0+1,b);
		}
		
		if(i!=0 && i!=n-1)
		{
			gsl_matrix_set(sys,i0+i,j0+i-1,b);
			gsl_matrix_set(sys,i0+i,j0+i,a);
			gsl_matrix_set(sys,i0+i,j0+i+1,b);
		}
		
		
		if(i==n-1)
		{
			gsl_matrix_set(sys,i0+n+1,j0+n-2,b);
			gsl_matrix_set(sys,i0+n+1,j0+n-1,a);
		}	
	}
}

void insertB(gsl_matrix * sys, int i0, int j0, int n)
{
	int i;
	const double a=f_B_3(1.0), b=f_B_3(sqrt(2.0));
	for(i=0;i<n;i++)
	{
		if(i==0)
		{
			gsl_matrix_set(sys, i0+0,j0+0,a);
			gsl_matrix_set(sys, i0+0,j0+1,b);
		}
		
		if(i!=0 && i!=n-1)
		{
			gsl_matrix_set(sys,i0+i,j0+i-1,b);
			gsl_matrix_set(sys,i0+i,j0+i,a);
			gsl_matrix_set(sys,i0+i,j0+i+1,b);
		}
		
		
		if(i==n-1)
		{
			gsl_matrix_set(sys,i0+n+1,j0+n-2,b);
			gsl_matrix_set(sys,i0+n+1,j0+n-1,a);
		}	
	}
}

void interp2D(double x1, double x2, double y1, double y2, int xs, int ys)
{
	int i,j,k;
	double outp, arg;
	const double 	xh=(x2-x1)/(xs-1.0), yh=(y2-y1)/(xs-1.0);
	//const double a=f_B_3(0.0),b=f_B_3(1.0),c=f_B_3(sqrt(2.0)) ;
	
	gsl_matrix * sys = gsl_matrix_alloc (xs*ys, xs*ys);
	gsl_vector * x = gsl_vector_alloc (xs*ys);
	gsl_vector * b = gsl_vector_alloc (xs*ys);
	
	for( i = 0; i < ys; i++ )
	{
		if(i==0)
		{
			insertA(sys,0,0,xs);
			insertB(sys,0,xs,xs);
		}		
		if(i!=0 && i!=ys-1)
		{
			insertB(sys,i*xs,(i-1)*xs,xs);
			insertA(sys,i*xs,i*xs,xs);
			insertB(sys,i*xs,(i+1)*xs,xs);
		}				
		if(i==ys-1)
		{
			insertB(sys,i*xs,(i-1)*xs,xs);
			insertA(sys,i*xs,i*xs,xs);
		}
	}
	
	
}

