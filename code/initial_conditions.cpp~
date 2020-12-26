//Program for applying Initial Conditions

#include "functions.h"

void initial_conditions()
{
//----------------------------------------------------------
// Assigning initial conditions for Velocities and pressure
//----------------------------------------------------------

usol=new double* [imax+1];
vsol=new double* [imax+1];
psol=new double* [imax+1];

for (int i=0; i < imax+1; i++)
usol[i]= new double[jmax+1];

for (int i=0; i < imax+1; i++)
vsol[i]= new double[jmax+1];

for (int i=0; i < imax+1; i++)
psol[i]= new double[jmax+1];

if(restart==1)
{
for(int i=0;i<imax+1;i++)
for(int j=0;j<jmax+1;j++)
usol[i][j]=0;

for(int i=0;i<imax+1;i++)
for(int j=0;j<jmax+1;j++)
vsol[i][j]=0;

for(int i=0;i<imax+1;i++)
for(int j=0;j<jmax+1;j++)
psol[i][j]=0;

}


if(restart==2)
{
string line;
ifstream fin;
fin.open("restart.dat",std::ios::in);

// Check for successful opening of file
if(fin.fail())
{
    cout<<"Error Opening File: restart.dat"<<endl;
    exit(1);
}

getline(fin,line);
getline(fin,line);

double xc_d,yc_d;
for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
fin>>usol[i][j]>>vsol[i][j]>>psol[i][j];

fin.close();

}

flag=new double* [imax+1];
for (int i=0; i < imax+1; i++)
flag[i]= new double[jmax+1];
for(int i=0;i<imax+1;i++)
for(int j=0;j<jmax+1;j++)
flag[i][j]=1;

// Finding the cells of points used to monitor the solution

if(npoints>0)
{
ipoint =new int[npoints];
jpoint =new int[npoints];

for(int n=0;n<npoints;n++)
{
int pointflag=0;

for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
{
/*double m1,c1;
m1=(y[i+1][j+1]-y[i+1][j])/(x[i+1][j+1]-x[i+1][j]+1e);
c1=y[i+1][j]-m1*x[i+1][j];
double m2,c2;
m2=(y[i][j+1]-y[i+1][j+1])/(x[i][j+1]-x[i+1][j+1]);
c2=y[i+1][j+1]-m2*x[i+1][j+1];
double m3,c3;
m3=(y[i][j]-y[i][j+1])/(x[i][j]-x[i][j+1]);
c3=y[i][j+1]-m3*x[i][j+1];
double m4,c4;
m4=(y[i+1][j]-y[i][j])/(x[i+1][j]-x[i][j]);
c4=y[i][j]-m3*x[i][j];
double valr,valt,vall,valb;
valr=ypoint[n]-m1*xpoint[n]-c1;
valt=ypoint[n]-m2*xpoint[n]-c2;
vall=ypoint[n]-m3*xpoint[n]-c3;
valb=ypoint[n]-m4*xpoint[n]-c4;
*/
double xmin,xmax,ymin,ymax;
xmin=min(x[i][j],x[i+1][j]);
xmin=min(xmin,x[i+1][j+1]);
xmin=min(xmin,x[i][j+1]);
ymin=min(y[i][j],y[i+1][j]);
ymin=min(ymin,y[i+1][j+1]);
ymin=min(ymin,y[i][j+1]);
xmax=max(x[i][j],x[i+1][j]);
xmax=max(xmax,x[i+1][j+1]);
xmax=max(xmax,x[i][j+1]);
ymax=max(y[i][j],y[i+1][j]);
ymax=max(ymax,y[i+1][j+1]);
ymax=max(ymax,y[i][j+1]);

if(xpoint[n]<=xmax&&xpoint[n]>=xmin&&
   ypoint[n]<=ymax&&ypoint[n]>=ymin)
{
ipoint[n]=i;
jpoint[n]=j;
pointflag=1;
}
}

if(pointflag==0)
{
ipoint[n]=imax/2;
jpoint[n]=jmax/2;
}
}
}
for(int n=0;n<npoints;n++)
cout<<ipoint[n]<<" "<<jpoint[n]<<endl;

}
