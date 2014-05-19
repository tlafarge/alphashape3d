#include <stdio.h>
#include <stdlib.h>
#include<R.h>
#include<Rmath.h>


//Return 1 if the ray intersect the triangle -1 if it does backward 0 if it doesn t and 2 if the point lie on the triangle
int rayTriangleIntersection(double * point,double * direction,double * triangles)
{
	double e[3]={0};
	double e2[3]={0};
	double h[3];
	double f[3];
	double d[3];
	int i=0,j=2,ii;
	int inside=0;
	while(i<3)
	{

		for( ii=0; ii<3;++ii)
		{
			e2[ii]=e[ii];
			e[ii]=triangles[3*i+ii]-triangles[3*j+ii];
		}
		for( ii=0; ii<3;++ii)
			h[ii]=point[ii]-triangles[3*j+ii];
		if(h[0]==0 && h[1]==0 && h[2]==0) //not on a vertices
			return(2);

		f[0] = e[1]*direction[2]-(e[2]*direction[1]);
	    f[1] = e[2]*direction[0]-(e[0]*direction[2]);
        f[2] = e[0]*direction[1]-(e[1]*direction[0]);
		d[i] = h[0]*f[0] + h[1]*f[1] + h[2]*f[2];
		j=i;
		++i;
	}
    if((d[0]>=0&&d[1]>=0&&d[2]>=0)||(d[0]<=0&&d[1]<=0&&d[2]<=0))
    {
    	inside=1;
		f[0] = e2[1]*e[2]-(e2[2]*e[1]);
	    f[1] = e2[2]*e[0]-(e2[0]*e[2]);
        f[2] = e2[0]*e[1]-(e2[1]*e[0]);
        if(direction[0]*f[0] + direction[1]*f[1] + direction[2]*f[2] < 0)
        {
        	for(ii=0; ii<3;++ii)
        		f[ii]=-f[ii];
        }
        double cote=h[0]*f[0] + h[1]*f[1] + h[2]*f[2];

    	 if(cote>0)
    		 inside=-1;
    	 if(cote==0)
    		 inside=2;
    }
	return(inside);
}



void pointinashape(int* triangles, int* nbTrin,double * coords,int* nbPin,  double * points,int* nin,int*inside,double * direction)
{
	int nbTr      = nbTrin[0];
	int n         = nin[0];
	int nbP       = nbPin[0];

	double x[3];

	double triangle[9];
	int i,j;
	int in;
	int res;
	for(i = 0; i < n; i++)
	{
		in=0;
		x[0]=points[i];
		x[1]=points[i+n];
		x[2]=points[i+2*n];
		j=0;
		while(j<nbTr && in!=-1)
		{
			triangle[0]=coords[triangles[j]-1];
			triangle[1]=coords[triangles[j]+nbP-1];
			triangle[2]=coords[triangles[j]+2*nbP-1];
			triangle[3]=coords[triangles[j+nbTr]-1];
			triangle[4]=coords[triangles[j+nbTr]+nbP-1];
			triangle[5]=coords[triangles[j+nbTr]+2*nbP-1];
			triangle[6]=coords[triangles[j+2*nbTr]-1];
			triangle[7]=coords[triangles[j+2*nbTr]+nbP-1];
			triangle[8]=coords[triangles[j+2*nbTr]+2*nbP-1];
			res=rayTriangleIntersection(x,direction,triangle);
			if(res==1)
				++in;
			if(res==2)
				in=-1;
			j++;
		}
		inside[i]=in;
	}

}
