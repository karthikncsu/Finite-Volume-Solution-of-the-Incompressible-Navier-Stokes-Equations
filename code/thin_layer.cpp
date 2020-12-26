// Function to discretize the poissons equation using Thin Layer Gradient model

#include "functions.h"

void thin_layer()
{

double *hx,*hy;
double volave,area2;

int indx =max(imax,jmax);
hx=new double[indx];
hy=new double[indx];

//"i" direction flux 
for(int j=1; j<jmax; j++)
{
for(int i=0; i<imax;i++)
{
area2=pow(del[i][j][0][0],2)+pow(del[i][j][0][1],2);	
volave=0.5*(vol[i+1][j]+vol[i][j]);
hx[i] = visc[parameter_num]*(usol[i+1][j]-usol[i][j])*area2/volave;
hy[i] = visc[parameter_num]*(vsol[i+1][j]-vsol[i][j])*area2/volave;
}
for(int i=1;i<imax;i++)
{
res[i][j][1]=res[i][j][1]-(hx[i]-hx[i-1]);
res[i][j][2]=res[i][j][2]-(hy[i]-hy[i-1]);
}
}

//"j" direction flux 
for(int i=1; i<imax; i++)
{
for(int j=0; j<jmax;j++)
{
area2=pow(del[i][j][1][0],2)+pow(del[i][j][1][1],2);	
volave=0.5*(vol[i][j+1]+vol[i][j]);
hx[j] = visc[parameter_num]*(usol[i][j+1]-usol[i][j])*area2/volave;
hy[j] = visc[parameter_num]*(vsol[i][j+1]-vsol[i][j])*area2/volave;
}
for(int j=1;j<jmax;j++)
{
res[i][j][1]=res[i][j][1]-(hx[j]-hx[j-1]);
res[i][j][2]=res[i][j][2]-(hy[j]-hy[j-1]);
}
}

delete [] hx;
delete [] hy;

}
