// Calculating the residue for immerssed boundary method

#include "functions.h"

void immersed_boundary()
{

int iq,jq;
double ubc,vbc,pbc,uave,vave,pave,unor,vnor,utan,vtan;
double udotn,xnxd,xnyd,dtv;
	
for(int j=1;j<jmax;j++)
{
for(int i=1;i<imax;i++)
{

if(hh[i][j]==1)
{
if(distg[i][j]<=0.0) // interior cell
{
ubc = 0.0;
vbc = 0.0;
pbc = 0.0;
}
else                          //band cell
{
uave = 0.0;
vave = 0.0;
pave = 0.0;
for(int k=1;k<nintpts;k++)
{
iq = i+idif[k];
jq = j+jdif[k];
uave = uave + weight[i][j][k]*usol[iq][jq];
vave = vave + weight[i][j][k]*vsol[iq][jq];
pave = pave + weight[i][j][k]*psol[iq][jq];
}
xnxd = xnxs[itag[i][j][lpri[i][j]]][lpri[i][j]];
xnyd = xnys[itag[i][j][lpri[i][j]]][lpri[i][j]];
udotn = uave*xnxd + vave*xnyd;
unor = udotn*xnxd;
vnor = udotn*xnyd;
utan = uave - unor;
vtan = vave - vnor;
ubc = utan*acoef[i][j] + unor*bcoef[i][j];
vbc = vtan*acoef[i][j] + vnor*bcoef[i][j];
pbc = pave;
}



if(time_accurate==1)
{
if(dt_crit==3) dtv = vol[i][j]/dt[i][j];
else dtv = vol[i][j]/dt_global;	
res[i][j][0] = dtv*(psol[i][j]-pbc)/beta2[i][j]; // beta is your AC parameter
res[i][j][1] = dtv*rho*(usol[i][j]-ubc); 
res[i][j][2] = dtv*rho*(vsol[i][j]-vbc);
}
else if(time_accurate==2)
{
dtv = vol[i][j]/tstep;
res[i][j][0] = dtv*(psol[i][j]-pbc)/beta2[i][j];
res[i][j][1] = dtv*rho*(usol[i][j]-ubc);
res[i][j][2] = dtv*rho*(vsol[i][j]-vbc);
	
}

}

}
} 


}
