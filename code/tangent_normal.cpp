// Function to discretize the poissons equation using Tangent/normal
// Decomposition method

#include "functions.h"

void tangent_normal()
{

double ***grad;
grad=new double** [imax+1];
for(int i=0;i<imax+1;i++) grad[i]=new double* [jmax+1];
for(int i=0;i<imax+1;i++) for(int j=0;j<jmax+1;j++) grad[i][j]=new double[2];

int indx=max(imax,jmax);
double *gflux;
gflux=new double[indx];

//X Momentum Equation
//---------------------------------------------------------
//Calculating Cell centered gradients using Green's theorem
//---------------------------------------------------------

for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
double f1,f2,f3,f4;
f1 = 0.5*(usol[i][j]+usol[i+1][j]);
f2 = 0.5*(usol[i][j]+usol[i][j+1]);
f3 = 0.5*(usol[i][j]+usol[i-1][j]);
f4 = 0.5*(usol[i][j]+usol[i][j-1]);
for(int k=0;k<2;k++)
{
grad[i][j][k] = f1*del[i][j][0][k] - f3*del[i-1][j][0][k]
			  + f2*del[i][j][1][k] - f4*del[i][j-1][1][k];
grad[i][j][k] = grad[i][j][k]/vol[i][j];
}
}

//----------------------------------
//Calculating I fluxes
//-----------------------------------
for(int j=1;j<jmax; j++)
{
for(int i=0; i<imax; i++)
{
double f1,f3;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double areaface,taudotn,graddottau;
double gradave1,gradave2;
double gradf1,gradf2;
f1 = usol[i+1][j];
f3 = usol[i][j];
delt1=xc[i+1][j]-xc[i][j];
delt2=yc[i+1][j]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
areaface=sqrt(pow(del[i][j][0][0],2)+pow(del[i][j][0][1],2));
xn1=del[i][j][0][0]/areaface;
xn2=del[i][j][0][1]/areaface;
taudotn=tau1*xn1+tau2*xn2;
gradave1= 0.5*(grad[i][j][0]+grad[i+1][j][0]);
gradave2= 0.5*(grad[i][j][1]+grad[i+1][j][1]);
graddottau = (gradave1*tau1+ gradave2*tau2)*taudotn;
gradf1=(f1-f3)*taudotn*xn1/dels+gradave1-graddottau*xn1;
gradf2=(f1-f3)*taudotn*xn2/dels+gradave2-graddottau*xn2;
gflux[i]=gradf1*del[i][j][0][0]+gradf2*del[i][j][0][1];
}
for(int i=1; i<imax; i++) res[i][j][1]=res[i][j][1]-visc[parameter_num]*(gflux[i]-gflux[i-1]);
}

//----------------------------------
//Calculating J fluxes
//-----------------------------------
for(int i=1;i<imax; i++)
{
for(int j=0;j<jmax; j++)
{
double f2,f4;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double areaface,taudotn,graddottau;
double gradave1,gradave2;
double gradf1,gradf2;
f2 = usol[i][j+1];
f4 = usol[i][j];
delt1=xc[i][j+1]-xc[i][j];
delt2=yc[i][j+1]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
areaface=sqrt(pow(del[i][j][1][0],2)+pow(del[i][j][1][1],2));
xn1=del[i][j][1][0]/areaface;
xn2=del[i][j][1][1]/areaface;
taudotn=tau1*xn1+tau2*xn2;
gradave1= 0.5*(grad[i][j][0]+grad[i][j+1][0]);
gradave2= 0.5*(grad[i][j][1]+grad[i][j+1][1]);
graddottau = (gradave1*tau1+ gradave2*tau2)*taudotn;
gradf1=(f2-f4)*taudotn*xn1/dels+gradave1-graddottau*xn1;
gradf2=(f2-f4)*taudotn*xn2/dels+gradave2-graddottau*xn2;
gflux[j]=gradf1*del[i][j][1][0]+gradf2*del[i][j][1][1];
}
for(int j=1; j<jmax; j++) res[i][j][1]=res[i][j][1]-visc[parameter_num]*(gflux[j]-gflux[j-1]);
}

//Y Momentum Equation
//---------------------------------------------------------
//Calculating Cell centered gradients using Green's theorem
//---------------------------------------------------------

for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
double f1,f2,f3,f4;
f1 = 0.5*(vsol[i][j]+vsol[i+1][j]);
f2 = 0.5*(vsol[i][j]+vsol[i][j+1]);
f3 = 0.5*(vsol[i][j]+vsol[i-1][j]);
f4 = 0.5*(vsol[i][j]+vsol[i][j-1]);
for(int k=0;k<2;k++)
{
grad[i][j][k] = f1*del[i][j][0][k] - f3*del[i-1][j][0][k]
			  + f2*del[i][j][1][k] - f4*del[i][j-1][1][k];
grad[i][j][k] = grad[i][j][k]/vol[i][j];
}
}

//----------------------------------
//Calculating I fluxes
//-----------------------------------
for(int j=1;j<jmax; j++)
{
for(int i=0; i<imax; i++)
{
double f1,f3;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double areaface,taudotn,graddottau;
double gradave1,gradave2;
double gradf1,gradf2;
f1 = vsol[i+1][j];
f3 = vsol[i][j];
delt1=xc[i+1][j]-xc[i][j];
delt2=yc[i+1][j]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
areaface=sqrt(pow(del[i][j][0][0],2)+pow(del[i][j][0][1],2));
xn1=del[i][j][0][0]/areaface;
xn2=del[i][j][0][1]/areaface;
taudotn=tau1*xn1+tau2*xn2;
gradave1= 0.5*(grad[i][j][0]+grad[i+1][j][0]);
gradave2= 0.5*(grad[i][j][1]+grad[i+1][j][1]);
graddottau = (gradave1*tau1+ gradave2*tau2)*taudotn;
gradf1=(f1-f3)*taudotn*xn1/dels+gradave1-graddottau*xn1;
gradf2=(f1-f3)*taudotn*xn2/dels+gradave2-graddottau*xn2;
gflux[i]=gradf1*del[i][j][0][0]+gradf2*del[i][j][0][1];
}
for(int i=1; i<imax; i++) res[i][j][2]=res[i][j][2]-visc[parameter_num]*(gflux[i]-gflux[i-1]);
}

//----------------------------------
//Calculating J fluxes
//-----------------------------------
for(int i=1;i<imax; i++)
{
for(int j=0;j<jmax; j++)
{
double f2,f4;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double areaface,taudotn,graddottau;
double gradave1,gradave2;
double gradf1,gradf2;
f2 = vsol[i][j+1];
f4 = vsol[i][j];
delt1=xc[i][j+1]-xc[i][j];
delt2=yc[i][j+1]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
areaface=sqrt(pow(del[i][j][1][0],2)+pow(del[i][j][1][1],2));
xn1=del[i][j][1][0]/areaface;
xn2=del[i][j][1][1]/areaface;
taudotn=tau1*xn1+tau2*xn2;
gradave1= 0.5*(grad[i][j][0]+grad[i][j+1][0]);
gradave2= 0.5*(grad[i][j][1]+grad[i][j+1][1]);
graddottau = (gradave1*tau1+ gradave2*tau2)*taudotn;
gradf1=(f2-f4)*taudotn*xn1/dels+gradave1-graddottau*xn1;
gradf2=(f2-f4)*taudotn*xn2/dels+gradave2-graddottau*xn2;
gflux[j]=gradf1*del[i][j][1][0]+gradf2*del[i][j][1][1];
}
for(int j=1; j<jmax; j++) res[i][j][2]=res[i][j][2]-visc[parameter_num]*(gflux[j]-gflux[j-1]);
}

for(int i=0;i<imax+1;i++) for(int j=0;j<jmax+1;j++) delete [] grad[i][j];
for(int i=0;i<imax+1;i++) delete [] grad[i];
delete [] grad;

delete [] gflux;

}
