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
			a[i][j] = f(hx*(double)i, hy*(double)j)/coef;
	/*
	//0.807221
	for(i = 0; i < nx; i++)
	{
		a[i][0] /= 0.807221;
		a[i][ny-1] /= 0.807221;
	}	
	
	for(j = 0; j < nx; j++)
	{
		a[0][j] /= 0.807221;
		a[nx-1][j] /= 0.807221;
	}
	
	a[0][0] 	/= 0.652648;
	a[nx-1][0] 	/= 0.652648;
	a[0][ny-1] 	/= 0.652648;
	a[nx-1][ny-1] 	/= 0.652648;*/
	
	double x, y, z;//z[nx*8][ny*8];
	FILE *outp;
	outp = fopen("./output/plot", "w");
	int m,n,p;
	for(m=0;m<(nx-1)*8;m++)
		for(n=0;n<(ny-1)*8;n++)
		{
			i=m/8;
			j=n/8;

			z=0.;
			
			if(i>0 && i<nx-2 && j>0 && j<ny-2)
				for(k=-1;k<=2;k++)
					for(p=-1;p<=2;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);
			
			if(i==0 && j>0 && j<ny-2)
			{
				for(k=0;k<=2;k++)
					for(p=-1;p<=2;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);
				z+=a[i][j+p]*phi((double)(m%8)/8.-2.0,(double)(n%8)/8.-(double)p);
			}

			if(i==nx-2 && j>0 && j<ny-2)
				for(k=-1;k<=1;k++)
					for(p=-1;p<=2;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);
			
			if(i>0 && i<nx-2 && j==0)
				for(k=-1;k<=2;k++)
					for(p=0;p<=2;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);

			if(i>0 && i<nx-2 && j==ny-2)
				for(k=-1;k<=2;k++)
					for(p=-1;p<=1;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);
					
			if(i==0 && j==0)
				for(k=0;k<=2;k++)
					for(p=0;p<=2;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);
						
			if(i==0 && j==ny-2)
				for(k=0;k<=2;k++)
					for(p=-1;p<=1;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);

			if(i==nx-2 && j==0)
				for(k=-1;k<=1;k++)
					for(p=0;p<=2;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);

			if(i==nx-2 && j==ny-2)
				for(k=-1;k<=1;k++)
					for(p=-1;p<=1;p++)
						z+=a[i+k][j+p]*phi((double)(m%8)/8.-(double)k,(double)(n%8)/8.-(double)p);


			fprintf(outp,"%f %f %f\n", (hx+1.)/64.*(double)m, (hy+1.)/64.*(double)n, z/1.828);
			//fprintf(outp,"%f %f %f\n",(hx+1.)/64.*(double)m,(hy+1.)/64.*(double)n,z/f((hx+1.)/64.*(double)m,(hy+1.)/64.*(double)n));
		}
}


int main()
{
	
	interp(0.,1.,8,0.,1.,8);
	
	return 0;
}
