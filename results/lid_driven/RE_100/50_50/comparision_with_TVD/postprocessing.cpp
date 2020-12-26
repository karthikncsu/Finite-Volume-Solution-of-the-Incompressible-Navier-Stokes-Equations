//Function to write data

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cmath>
using namespace std;

int main()
{
int imax,jmax;

string line,filename,filename1,filename2;

cout<<"Enter the file name:"<<endl;
cin>>filename;

ifstream fin;
fstream fout1,fout2;

fin.open(filename.c_str(),std::ios::in);

getline(fin,line);
fin>>line;
fin>>line;
fin>>imax;
fin>>line;
fin>>jmax;
getline(fin,line);

double xc[imax+2][jmax+2],yc[imax+2][jmax+2], usol[imax+2][jmax+2], vsol[imax+2][jmax+2];
double umag[imax+1][jmax+1],psol[imax+1][jmax+1];

cout<<imax<<" "<<jmax<<endl;

//double xmin=1e10,ymin=1e10,xmax=1e-10,ymax=1e-10;
for(int j=1;j<jmax+1;j++)
for(int i=1;i<imax+1;i++)
{
fin>>xc[i][j]>>yc[i][j]>>usol[i][j]>>vsol[i][j]>>umag[i][j]>>psol[i][j];
cout<<xc[i][j]<<" "<<yc[i][j]<<" "<<usol[i][j]<<" "<<vsol[i][j]<<" "<<umag[i][j]<<" "<<psol[i][j]<<endl;
}

//Assigning Ghost Cells
for(int j=1; j<jmax+1;j++)
{
xc[0][j]=2*xc[1][j]-xc[2][j];
yc[0][j]=2*yc[1][j]-yc[2][j];
xc[imax+1][j]=2*xc[imax][j]-xc[imax-1][j];
yc[imax+1][j]=2*yc[imax][j]-yc[imax-1][j];
}

for(int i=1; i<imax+1;i++)
{
xc[i][0]=2*xc[i][1]-xc[i][2];
yc[i][0]=2*yc[i][1]-yc[i][2];
xc[i][jmax+1]=2*xc[i][jmax]-xc[i][jmax-1];
yc[i][jmax+1]=2*yc[i][jmax]-yc[i][jmax-1];
}
xc[0][0]=2*xc[1][0]-xc[2][0];
yc[0][0]=2*yc[0][1]-yc[0][2];
xc[imax+1][0]=2*xc[imax][0]-xc[imax-1][0];
yc[imax+1][0]=2*yc[imax+1][1]-yc[imax+1][2];
xc[0][jmax+1]=2*xc[1][jmax+1]-xc[2][jmax+1];
yc[0][jmax+1]=2*yc[0][jmax]-yc[0][jmax-1];
xc[imax+1][jmax+1]=2*xc[imax][jmax+1]-xc[imax-1][jmax+1];
yc[imax+1][jmax+1]=2*yc[imax+1][jmax]-yc[imax+1][jmax-1];

fin.close();

double m,c;
int num_points;
double xstart,ystart,xend,yend;
//double length;

cout<<"Enter the start and end point of the line: "<<endl;
cout<<"xstart: ";
cin>>xstart;
cout<<"ystart: ";
cin>>ystart;
cout<<"xend: ";
cin>>xend;
cout<<"yend: ";
cin>>yend;
m=(yend-ystart)/(xend-xstart+1e-10);
c=ystart-m*xstart;
//length=sqrt(pow(yend-ystart,2)+pow(xend-xstart,2));
cout<<"Enter the number of points along the line for plotting: ";
cin>>num_points;

double xintop[num_points],yintop[num_points];
double u[num_points],v[num_points],p[num_points],uabs[num_points];
for(int intop=0;intop<num_points;intop++)
{
if(xend!=xstart)
{
xintop[intop]=xstart+(xend-xstart)/(num_points-1)*intop;
yintop[intop]=m*xintop[intop]+c;
}
else if(xend==xstart)
{
xintop[intop]=xstart;
yintop[intop]=ystart+(yend-ystart)/(num_points-1)*intop;
}

}

double xmin,xmax,ymin,ymax;

for(int intop=0;intop<num_points;intop++)
{

for(int j=1;j<jmax+1;j++)
{
for(int i=1;i<imax+1;i++)
{
xmin=min(xc[i][j],xc[i+1][j]);
xmin=min(xmin,xc[i+1][j+1]);
xmin=min(xmin,xc[i][j+1]);
ymin=min(yc[i][j],yc[i+1][j]);
ymin=min(ymin,yc[i+1][j+1]);
ymin=min(ymin,yc[i][j+1]);
xmax=max(xc[i][j],xc[i+1][j]);
xmax=max(xmax,xc[i+1][j+1]);
xmax=max(xmax,xc[i][j+1]);
ymax=max(yc[i][j],yc[i+1][j]);
ymax=max(ymax,yc[i+1][j+1]);
ymax=max(ymax,yc[i][j+1]);
if(xintop[intop]<=xmax&&xintop[intop]>=xmin&&yintop[intop]<=ymax&&yintop[intop]>=ymin)
{
//cout<<intop<<" "<<(xintop[intop]<=xmax)<<" "<<(xintop[intop]>=xmin)<<" "<<(yintop[intop]<=ymax)<<" "<<(yintop[intop]>=ymin)<<endl;	
//cout<<"entering the loop "<<xintop[intop]<<" "<<xmin<<" "<<xmax<<" "<<yintop[intop]<<" "<<ymin<<" "<<ymax<<" "<<intop<<endl;
u[intop]=usol[i][j];
v[intop]=vsol[i][j];
uabs[intop]=umag[i][j];
p[intop]=psol[i][j];
}
}
}

}


fstream fout;
fout.open("results.dat", ios::out | ios::trunc);

for(int intop=0;intop<num_points;intop++)
fout<<xintop[intop]<<" "<<yintop[intop]<<" "<<u[intop]<<" "<<v[intop]<<" "<<uabs[intop]<<" "<<p[intop]<<endl;;

}
