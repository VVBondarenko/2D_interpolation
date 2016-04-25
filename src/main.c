#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic_functions.c"
#include "B-splines.c"
#include "smplDynArray.c"

double f(double x, double y)
{
	return exp(x*(x-1.))*exp(y*(y-1.));
}

double phi(double x, double y)
{
	return f_B_3(x)*f_B_3(y);
}

void interp(double x0, double x1, int nx, double y0, double y1, int ny)
{
	const double hx = (x1-x0)/(nx-1.), hy = (y1-y0)/(ny-1.),
			coef = 4.*(phi(1.0,0.0)+phi(1.0,1.0)),
			c1 = phi(0.,1.), c2 = phi(1.,1.);
	int i,j,k;
	double **a,**an;
	
	init2DArr(&a,nx,ny);
	init2DArr(&an,nx,ny);

	for(i = 0; i < nx; i++)
		for(j = 0; j < ny; j++)
			a[i][j] = f(x0+hx*(double)i, y0+hy*(double)j)/coef;
	
	for(k = 0; k <= 1000; k++)
	{

/*			for(j=1;j<=ny-1;j++)	*/
/*				a[i][j] = 	 a[i-1][j+1]*phi(-1.,1.) + a[i][j+1]*phi(0.,1.) + a[i+1][j+1]*phi(1.,1.)	*/
/*						+a[i-1][j]*phi(-1.,0.) 				+ a[i+1][j]*phi(-1.,0.)		*/
/*						+a[i-1][j+1]*phi(-1.,1.) + a[i][j-1]*phi(0.,1.) + a[i+1][j-1]*phi(1.,-1.)	*/
/*						-f(hx*(double)i, hy*(double)j);	*/
		//внутри области
		for(i=1;i<nx-1;i++)
			for(j=1;j<=ny-1;j++)
				an[i][j] = 	-a[i-1][j+1]*c2 - a[i][j+1]*c1  - a[i+1][j+1]*c2
						-a[i-1][j]*c1 			- a[i+1][j]*c1
						-a[i-1][j-1]*c2 - a[i][j-1]*c1  - a[i+1][j-1]*c2
						+f(hx*(double)i, hy*(double)j);
		//на границах
		i=0;
		for(j=1;j<ny-1;j++)
			an[i][j] = 	 - a[i][j+1]*c1  - a[i+1][j+1]*c2
					 		 - a[i+1][j]*c1
					 - a[i][j-1]*c1  - a[i+1][j-1]*c2
					 +f(hx*(double)i, hy*(double)j);
		i=nx-1;
		for(j=1;j<ny-1;j++)
			an[i][j] = 	-a[i-1][j+1]*c2 - a[i][j+1]*c1  
					-a[i-1][j]*c1
					-a[i-1][j-1]*c2 - a[i][j-1]*c1
					+f(hx*(double)i, hy*(double)j);
		
		j=0;
		for(i=1;i<nx-1;i++)
			an[i][j] = 	-a[i-1][j+1]*c2 - a[i][j+1]*c1  - a[i+1][j+1]*c2
					-a[i-1][j]*c1 			- a[i+1][j]*c1
					
					+f(hx*(double)i, hy*(double)j);
		
		j=ny-1;
		for(i=1;i<nx-1;i++)
			an[i][j] = 	
					-a[i-1][j]*c1 			- a[i+1][j]*c1
					-a[i-1][j-1]*c2 - a[i][j-1]*c1  - a[i+1][j-1]*c2
					+f(hx*(double)i, hy*(double)j);
					
		
		i=0;j=0;
			an[i][j] = 	- a[i][j+1]*c1  - a[i+1][j+1]*c2
							- a[i+1][j]*c1
					
					+f(hx*(double)i, hy*(double)j);
					
		i=0;j=ny;
			an[i][j] = 	
							 - a[i+1][j]*c1
					 - a[i][j-1]*c1  - a[i+1][j-1]*c2
					+f(hx*(double)i, hy*(double)j);
		
		i=nx;j=0;
			an[i][j] = 	-a[i-1][j+1]*c2 - a[i][j+1]*c1  
					-a[i-1][j]*c1 			
					
					+f(hx*(double)i, hy*(double)j);
					
		i=nx;j=ny;
			an[i][j] = 	
					-a[i-1][j]*c1 			
					-a[i-1][j-1]*c2 - a[i][j-1]*c1  
					+f(hx*(double)i, hy*(double)j);
		
		//копирование
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				a[i][j] = an[i][j];		
	}
	double x, y;
	for(x=x0+hx; x<x1-hx; x+=0.01)
		for(y=y0+hy; y<y1-hy; y+=0.01)
		{
			i = (int)(x-x0)/hx;
			j = (int)(y-y0)/hy;

printf("%f %f %f\n",x,y,a[i-1][j+1]*phi((x-x0)/hx-1.,(y-y0)/hy+1.) + a[i][j+1]*phi((x-x0)/hx,(y-y0)/hy+1.) + a[i+1][j+1]*phi((x-x0)/hx+1.,(y-y0)/hy+1.)	
			+a[i-1][j]*phi((x-x0)/hx-1.,(y-y0)/hy) 	+ a[i][j]*phi((x-x0)/hx,(y-y0)/hy) + a[i+1][j]*phi((x-x0)/hx+1.,(y-y0)/hy)		
			+a[i-1][j+1]*phi((x-x0)/hx-1.,(y-y0)/hy-1.) + a[i][j-1]*phi((x-x0)/hx,(y-y0)/hy-1.) + a[i+1][j-1]*phi((x-x0)/hx+1.,(y-y0)/hy-1.)	
					-f(hx*(double)i, hy*(double)j));				
			
		}	

}


int main()
{
	
	interp(0.,1.,8,0.,1.,8);
	
	return 0;
}
