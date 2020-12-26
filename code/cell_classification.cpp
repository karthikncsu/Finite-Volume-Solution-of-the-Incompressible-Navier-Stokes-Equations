//Functions for cell classification in immerssed boundary method

#include "functions.h"

void cell_classification()
{
double dd;
double dotp;
double xnxd,xnyd,xcd,ycd,xcp,ycp,dx,dy;
double deldsum1,deldsum1d;
double distance,deld,dcross,dratio,diste,distd,term;
int iq,jq;
int dummy;

//--------------------------------------------
// Step 1a:  determine local signed distance
//---------------------------------------------
for(int n=0;n<nobjs;n++)
{
for(int j=0;j<jmax+1;j++)
for(int i=0;i<imax+1;i++)
{
dist[i][j][n]=1000;
for(int k=0; k<nlist[n]; k++)	
{
dd=pow(xc[i][j]-xs[k][n],2)+pow(yc[i][j]-ys[k][n],2);
	if(dd<dist[i][j][n])
	{
	itag[i][j][n]=k;
	dist[i][j][n]=dd;
	}
}
dist[i][j][n]=sqrt(dist[i][j][n]);
dotp=(xc[i][j]-xs[itag[i][j][n]][n])*xnxs[itag[i][j][n]][n];
	+(yc[i][j]-ys[itag[i][j][n]][n])*xnys[itag[i][j][n]][n];
dist[i][j][n]=dist[i][j][n]*sign(1,dotp);
}
}


//--------------------------------------------------
// Step 1b:  determine global signed distance
//--------------------------------------------------

for(int j=1;j<jmax;j++)
{
for(int i=1;i<imax;i++)
{
distg[i][j]=1e12;
lpri[i][j]=1;
for(int n=0;n<nobjs;n++)
{
distg[i][j]=min(distg[i][j],dist[i][j][n]);
if(distg[i][j]==dist[i][j][n]) lpri[i][j]=n;
}
}
}

//----------------------------------------------------
//  Step 1c:  determine Heaviside function
//----------------------------------------------------

for(int j=1;j<jmax;j++)
{
for(int i=1;i<imax;i++)
{
int icflag=0;
hh[i][j]=0;
for(int k=0;k<nintpts;k++)
{
iq=i+idif[k];
jq=j+jdif[k];
if(distg[i][j]>0.0 && distg[iq][jq] <0.0) icflag=1;
}
if(icflag==1) hh[i][j]=1.0;
if(distg[i][j]<0.0) hh[i][j]=1.0;
}
}

for(int j=1;j<jmax;j++)
{
hh[0][j]=hh[1][j];
hh[imax][j]=hh[imax-1][j];
}
for(int i=1;i<imax;i++)
{
hh[i][0]=hh[i][0];
hh[i][jmax]=hh[i][jmax-1];
}

fstream fout;
fout.open("domain.dat", ios::out | ios::trunc);
fout<<"Variables = x, y, hh"<<endl;
fout<<"ZONE I= "<<imax-1<<" J= "<<jmax-1<<" F=point"<<endl;

for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
{
fout<<xc[i][j]<<" "<<yc[i][j]<<" "<<hh[i][j]<<endl;
}
fout.close();


/*
/*c -------------------------------------------------------------
c ---- now you have four arrays that define cell classification
c
c      distg(i,j) = global signed distance function
c      hh(i,j) = Heaviside function (1.0 if interior or band; 0.0 otherwise)
c      lpri(i,j) = tag that tells which object is closest
c      itag(i,j,l) = index for nearest point on object 'l'
c -------------------------------------------------------------

c ---------------------------------------------------------
c ---- Step 2:  Determine interpolation data
c ---------------------------------------------------------
*/



for(int j=1;j<jmax;j++)
{
for(int i=1;i<imax;i++)
{
	
if(hh[i][j]==1.0 && distg[i][j]>0.0)
{

for(int k=0;k<nintpts;k++) weight[i][j][k]=0.0;       

xnxd = xnxs[itag[i][j][lpri[i][j]]][lpri[i][j]];
xnyd = xnys[itag[i][j][lpri[i][j]]][lpri[i][j]];
xcd = xc[i][j];
ycd = yc[i][j]; 

deldsum1 = 0.0;
deldsum1d = 0.0;

for(int k=0;k<nintpts ;k++)
{
iq = i + idif[k];
jq = j + jdif[k];
xcp = xc[iq][jq];
ycp = yc[iq][jq];
dx = xcp - xcd;
dy = ycp - ycd;
distance = sqrt(pow(dx,2) + pow(dy,2));
deld = dx*xnxd + dy*xnyd;
dcross = sqrt(abs(pow(distance,2)-pow(deld,2)));

if(deld>0.0&&hh[iq][jq]==0.0) //Consider only field cells
{
dcross = 1.0/(dcross + 1e-12);
weight[i][j][k] = dcross;
deldsum1 = deldsum1 + dcross;
deldsum1d = deldsum1d + dcross*deld;
}
}

if(deldsum1==0.0)                    // consider band and field cells
{
deldsum1 = 0.0;
deldsum1d = 0.0;

for(int k=0;k<nintpts;k++)
{
iq = i + idif[k];
jq = j + jdif[k];
xcp = xc[iq][jq];
ycp = yc[iq][jq];
dx = xcp - xcd;
dy = ycp - ycd;
distance = sqrt(pow(dx,2) + pow(dy,2));
deld = dx*xnxd + dy*xnyd;
dcross = sqrt(abs(pow(distance,2)-pow(deld,2)));

if(deld>0.0 && hh[iq][jq]>=0.0)
{
dcross = 1.0/(dcross + 1e-12);
weight[i][j][k] = dcross;
deldsum1 = deldsum1 + dcross;
deldsum1d = deldsum1d + dcross*deld;
}
}
}

if(deldsum1!=0.0)
{
for(int k=0;k<nintpts;k++)
weight[i][j][k] = weight[i][j][k]/deldsum1;
deld = deldsum1d/deldsum1;
}

else
{
cout<<"No interpolation point found"<<endl;
for(int k=0;k<nintpts;k++)
weight[i][j][k] = 0.0;
deld = 0.0;
}

dratio = deld/distg[i][j];
acoef[i][j] = pow(1.0/(1.0+dratio),power);
diste = distg[i][j] + 0.5*deld;
distd = 0.5*distg[i][j];
term = pow(distd/diste,power)/dratio;
bcoef[i][j] = term/(1.0 + term);           //Attempts to satisfy continuity
// bcoef[i][j] = distg[i][j]/(deld+distg[i][j]) // ! linear behavior
}

}
}


for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
if(hh[i][j]==1) cout<<i<<" "<<j<<" "<<hh[i][j]<<" "<<acoef[i][j]<<" "<<bcoef[i][j]<<endl;


}
