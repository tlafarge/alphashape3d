#include <stdio.h>
#include <stdlib.h>
#include<R.h>
#include<Rmath.h>


int includedin(int* value,int* tab)
{

	if (tab[0]==value[0]&&tab[1]==value[1]&&tab[2]==value[2])
		return(3);
	if (tab[0]==value[0]&&tab[1]==value[1]&&tab[3]==value[2])
		return(2);
	if (tab[0]==value[0]&&tab[2]==value[1]&&tab[3]==value[2])
		return(1);
	if (tab[1]==value[0]&&tab[2]==value[1]&&tab[3]==value[2])
		return(0);

	return(-1);
}

void triangleNormals(int* triangles,int* nbTrin, int*tetra, int*nbTetrain, double * points,int* nbPin, double* normalMat, double* middlePoint)
{
	int nbTr      = nbTrin[0];
	int nbTetra = nbTetrain[0];
	int nbP       = nbPin[0];

	int tr[3];
	int ta[4];
	double v1[3];
	double v2[3];
	int i;
	int signe;
	for(i = 0; i < nbTr; i++)
	{
		 tr[0]=triangles[i];
		 tr[1]=triangles[i+nbTr];
		 tr[2]=triangles[i+2*nbTr];
		 int j = 0;
		 ta[0]=tetra[j];
		 ta[1]=tetra[j+nbTetra];
		 ta[2]=tetra[j+2*nbTetra];
		 ta[3]=tetra[j+3*nbTetra];


		 int position= includedin(tr,ta);

		 while ( j< nbTetra && position==-1)
		 {
		 	 j++;
		 	 ta[0]=tetra[j];
		 	 ta[1]=tetra[j+nbTetra];
		 	 ta[2]=tetra[j+2*nbTetra];
		 	 ta[3]=tetra[j+3*nbTetra];
		 	 position= includedin(tr,ta);
		 }
		 int point4 = tetra[j+position*nbTetra];

		 v1[0]=points[tr[1] -1]-points[tr[0]-1];
		 v1[1]=points[tr[1] +nbP   -1]-points[tr[0]+nbP-1];
		 v1[2]=points[tr[1]+2*nbP  -1]-points[tr[0]+2*nbP-1];
		 v2[0]=points[tr[2] -1]-points[tr[0]-1];
		 v2[1]=points[tr[2] +nbP   -1]-points[tr[0]+nbP-1];
		 v2[2]=points[tr[2]+2*nbP  -1]-points[tr[0]+2*nbP-1];

		 normalMat[i]       = v1[1]*v2[2]-(v1[2]*v2[1]);
	     normalMat[i+nbTr]  = v1[2]*v2[0]-(v1[0]*v2[2]);
         normalMat[i+2*nbTr]= v1[0]*v2[1]-(v1[1]*v2[0]);

		 v1[0]=(points[triangles[i]-1]       + points[triangles[i+nbTr]-1]       + points[triangles[i+2*nbTr]-1]       )/3;
		 v1[1]=(points[triangles[i]+nbP-1]   + points[triangles[i+nbTr]+nbP-1]   + points[triangles[i+2*nbTr]+nbP-1]   )/3;
		 v1[2]=(points[triangles[i]+2*nbP-1] + points[triangles[i+nbTr]+2*nbP-1] + points[triangles[i+2*nbTr]+2*nbP-1] )/3;
		 v2[0]=v1[0]-points[point4-1];
		 v2[1]=v1[1]-points[point4+nbP-1];
		 v2[2]=v1[2]-points[point4+2*nbP-1];

		 middlePoint[i]       = v1[0];
		 middlePoint[i+nbTr]  = v1[1];
		 middlePoint[i+2*nbTr]= v1[2];

		 signe = (normalMat[i]*v2[0] +  normalMat[i+nbTr]*v2[1] + normalMat[i+2*nbTr]*v2[2]) > 0;
         if ( !signe)
         {
        	 normalMat[i]       = -normalMat[i]  ;
        	 normalMat[i+nbTr]  = -normalMat[i+nbTr];
        	 normalMat[i+2*nbTr]= -normalMat[i+2*nbTr];
         }
         triangles[i]=signe;



	}

}

