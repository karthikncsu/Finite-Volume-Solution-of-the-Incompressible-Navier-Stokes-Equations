//Function to impelment time accurate method

#include "functions.h"

void timeaccurate_ib_method()
{

double area=0, xnx=0, xny=0;
double uave=0, vave=0, pave=0;
double delp=0,pleft=0, pright=0;
int ind_max, Nconv;
double a11ave,a12ave;
double b11ave,b12ave;
double areaaave,areabave;
double vdotna,vdotnb,eiga,eigb,voldt;
double dt_global=1e10;
double dtv,deltap,deltau,deltav;
double xmassflux,udotn,btave2=0,dissip;
double a,b,signa,signb,minab;
double dtmax=0,dtmin=1e10;

ind_max=max(imax,jmax);
double **h;	
h=new double* [ind_max];
for (int i=0; i < ind_max; i++)
h[i]= new double[3];

double *du;
du=new double[ind_max];

double **beta2;
beta2=new double* [imax+1];
for (int i=0; i < imax+1; i++)
beta2[i]=new double[jmax+1];

res=new double** [imax];
for (int i=0; i < imax; i++) res[i]= new double* [jmax];
for(int i=0; i <imax; i++) for(int j=0; j<jmax; j++)
res[i][j]=new double[3];

double **dt;
dt=new double* [imax+1];
for (int i=0; i < imax+1; i++) dt[i]= new double[jmax+1];

double mindelta_global=1e10,maxeiga=-1e10;

//Calcualitng Beta values
for(int i=1;i<imax;i++)
for(int j=1;j<jmax;j++)
{
double deltax1,deltax2,deltax3,deltax4;
deltax1=sqrt(pow(x[i+1][j+1]-x[i+1][j],2)
			+pow(y[i+1][j+1]-y[i+1][j],2));
deltax2=sqrt(pow(x[i+1][j+1]-x[i][j+1],2)
			+pow(y[i+1][j+1]-y[i][j+1],2));
deltax3=sqrt(pow(x[i][j+1]-x[i][j],2)
			+pow(y[i][j+1]-y[i][j],2));
deltax4=sqrt(pow(x[i][j]-x[i+1][j],2)
			+pow(y[i][j]-y[i+1][j],2));
mindelta_global=min(mindelta_global,deltax1);
mindelta_global=min(mindelta_global,deltax2);
mindelta_global=min(mindelta_global,deltax3);
mindelta_global=min(mindelta_global,deltax4);
}

string filename;
fstream fres;
filename=mesh_name+"_residue_RE_";
ostringstream convert;
convert << RE;
filename=filename+convert.str();
filename=filename+".dat";
fres.open(filename.c_str(), ios::out | ios::trunc);
fres<<"Variables = Iteration, log_Relative_Residue_momentum, Residue_momentum, log_RR_cont, Residue_continuity "<<endl;

fstream fpoint;
filename=mesh_name+"_point_values_RE_";
filename=filename+convert.str();
filename=filename+".dat";
fpoint.open(filename.c_str(), ios::out | ios::trunc);

double **uold;
uold=new double* [imax+1];
for (int i=0; i < imax+1; i++)
uold[i]=new double[jmax+1];

double **vold;
vold=new double* [imax+1];
for (int i=0; i < imax+1; i++)
vold[i]=new double[jmax+1];

int converged=1;
int num_nonconv=0;

int solsave_flag=-1;
double deltaT_store=tstep;

tstep=1e-5;

for(phy_time=tstart; phy_time<tend; phy_time=phy_time+tstep)
{
solsave_flag++;
if(converged==1)
{	
for(int j=0; j<jmax;j++) for(int i=0;i<imax;i++) uold[i][j]=usol[i][j];
for(int j=0; j<jmax;j++) for(int i=0;i<imax;i++) vold[i][j]=vsol[i][j];
num_nonconv=0;
converged=0;
}
else 
{
for(int j=0; j<jmax;j++) for(int i=0;i<imax;i++) usol[i][j]=uold[i][j];
for(int j=0; j<jmax;j++) for(int i=0;i<imax;i++) vsol[i][j]=vold[i][j];
num_nonconv++;
if(num_nonconv<10)
{
cout<<"Solution did not converge for mimimum time step"<<endl;
}

}	
//Iterations Started
for(int iter=0;iter<Niter;iter++)
{
	
boundary_conditions();

//Calcualitng Beta values
for(int i=0;i<imax+1;i++)
for(int j=0;j<jmax+1;j++)
{
double deltax1,deltax2,deltax3,deltax4,mindelta;
double Umag=sqrt(pow(usol[i][j],2)+pow(vsol[i][j],2));
beta2[i][j]=max(Uref,Umag);
deltax1=sqrt(pow(x[i+1][j+1]-x[i+1][j],2)
			+pow(y[i+1][j+1]-y[i+1][j],2));
deltax2=sqrt(pow(x[i+1][j+1]-x[i][j+1],2)
			+pow(y[i+1][j+1]-y[i][j+1],2));
deltax3=sqrt(pow(x[i][j+1]-x[i][j],2)
			+pow(y[i][j+1]-y[i][j],2));
deltax4=sqrt(pow(x[i][j]-x[i+1][j],2)
			+pow(y[i][j]-y[i+1][j],2));

mindelta=min(deltax1,deltax2);
mindelta=min(mindelta,deltax3);
mindelta=min(mindelta,deltax4);
beta2[i][j]=max(beta2[i][j],visc[parameter_num]/rho/mindelta);
beta2[i][j]=pow(beta2[i][j],2);
}

//End of Calculaitng hte beta values	

for(int i=0; i <imax; i++)
for(int j=0; j<jmax; j++)
for(int k=0; k<3; k++)
res[i][j][k]=0;

//--------------------------------------------------
// Calculating the Inviscid fluxes
//--------------------------------------------------
//"i" direction flux 
for(int j=1; j<jmax; j++)
{
//--------------------------------------------------
// Rhie-chow calculation for pressure velocity coupling
//--------------------------------------------------
if(TVD_method==1)
{ 
for(int i=1;i<imax;i++)
du[i] = 0.5*(psol[i+1][j]-psol[i][j])+0.5*(psol[i][j]-psol[i-1][j]);
du[0] = du[1];
du[imax] = du[imax-1];
}

if(TVD_method==2)
{ 
for(int i=1;i<imax;i++)
{
a=psol[i+1][j]-psol[i][j];
b=psol[i][j]-psol[i-1][j];

if(a!=0) signa=a/abs(a);
else signa=0;
if(b!=0) signb=b/abs(b);
else signb=0;

minab=min(abs(a),abs(b));

du[i]=0.5*(signa+signb)*minab;
}
du[0] = du[1];
du[imax] = du[imax-1];
}
//--------------------------------------------------
// End of Rhie-chow calculation for pressure velocity coupling
//--------------------------------------------------

for(int i=0; i<imax;i++)
{
area=sqrt(pow(del[i][j][0][0],2) + pow(del[i][j][0][1],2));
xnx = del[i][j][0][0]/area;
xny = del[i][j][0][1]/area;
uave = 0.5*(usol[i][j]+usol[i+1][j]);
vave = 0.5*(vsol[i][j]+vsol[i+1][j]);
pave = 0.5*(psol[i][j]+psol[i+1][j]);
udotn = uave*xnx + vave*xny;
btave2 = 0.5*(beta2[i][j]+beta2[i+1][j]);
xmassflux = area*rho*udotn*flag[i][j];
eiga = 0.5*(abs(udotn) + sqrt(pow(udotn,2) + 4.0*btave2));

//--------------------------------------------------
//	Implementation of Rhie-chow for pressure velocity coupling
//--------------------------------------------------
if(TVD_method==0) 
{
delp = psol[i][j]-psol[i+1][j];
}
else if(TVD_method==1 || TVD_method==2)
{
pleft = psol[i][j] + 0.5*du[i];
pright = psol[i+1][j] - 0.5*du[i+1];
delp = pleft-pright;
}
//--------------------------------------------------
//End of implementation
//--------------------------------------------------
dissip  = 0.5*area*eiga*delp/btave2;
h[i][0] = xmassflux + dissip;
h[i][1] = xmassflux*uave + xnx*area*pave;
h[i][2] = xmassflux*vave + xny*area*pave;
}
for(int i=1;i<imax;i++)
for(int k=0; k<3;k++)
res[i][j][k]=res[i][j][k]+h[i][k]-h[i-1][k];

}

//"j" direction flux 
for(int i=1;i<imax;i++)
{
//--------------------------------------------------
// Rhie-chow calculation for pressure velocity coupling
//--------------------------------------------------
if(TVD_method==1)
{ 
for(int j=1;j<jmax;j++)
du[j] = 0.5*(psol[i][j+1]-psol[i][j])+0.5*(psol[i][j]-psol[i][j-1]);
du[0] = du[1];
du[jmax] = du[jmax-1];
}
if(TVD_method==2)
{ 
for(int j=1;j<jmax;j++)
{
a=psol[i][j+1]-psol[i][j];
b=psol[i][j]-psol[i][j-1];

if(a!=0) signa=a/abs(a);
else signa=0;
if(b!=0) signb=b/abs(b);
else signb=0;

minab=min(abs(a),abs(b));

du[j]=0.5*(signa+signb)*minab;
}
du[0] = du[1];
du[jmax] = du[jmax-1];
}
//--------------------------------------------------
// End of Rhie-chow calculation for pressure velocity coupling
//--------------------------------------------------

for(int j=0;j<jmax;j++)
{
area=sqrt(pow(del[i][j][1][0],2) + pow(del[i][j][1][1],2));
xnx = del[i][j][1][0]/area;
xny = del[i][j][1][1]/area;
uave = 0.5*(usol[i][j]+usol[i][j+1]);
vave = 0.5*(vsol[i][j]+vsol[i][j+1]);
pave = 0.5*(psol[i][j]+psol[i][j+1]);
udotn = uave*xnx + vave*xny;
btave2 = 0.5*(beta2[i][j]+beta2[i][j+1]);
xmassflux = area*rho*udotn*flag[i][j];
eiga = 0.5*(abs(udotn) + sqrt(pow(udotn,2) + 4.0*btave2));

//--------------------------------------------------
//	Implementation of Rhie-chow for pressure velocity coupling
//--------------------------------------------------
if(TVD_method==0) 
{
delp = psol[i][j]-psol[i][j+1];
}
else if(TVD_method==1 || TVD_method==2)
{
pleft = psol[i][j] + 0.5*du[j];
pright = psol[i][j+1] - 0.5*du[j+1];
delp = pleft-pright;
}
//--------------------------------------------------
//End of implementation
//--------------------------------------------------
dissip  = 0.5*area*eiga*delp/btave2;

h[j][0] = xmassflux+dissip;
h[j][1] = xmassflux*uave + xnx*area*pave;
h[j][2] = xmassflux*vave + xny*area*pave;
}

for(int j=1;j<jmax;j++)
for(int k=0; k<3;k++)
res[i][j][k]=res[i][j][k]+h[j][k]-h[j-1][k];
}

//--------------------------------------------------
// End of calculating the Inviscid fluxes
//--------------------------------------------------

//--------------------------------------------------
//      Addition of viscid fluxes
//-------------------------------------------------- 


if(fv_method==1) thin_layer();
else if(fv_method==2) tangent_normal();
//--------------------------------------------------
//      End of addition of viscid fluxes
//--------------------------------------------------

//--------------------------------------------------
// Calculation of local time step (average areas to cell centers)
//--------------------------------------------------
dt_global=1e10;
maxeiga=-1e10;
for(int j=1;j<jmax; j++)
{
for(int i=1;i<imax;i++)
{
a11ave = 0.5*(del[i][j][0][0] + del[i-1][j][0][0]);
a12ave = 0.5*(del[i][j][0][1] + del[i-1][j][0][1]);
b11ave = 0.5*(del[i][j][1][0] + del[i][j-1][1][0]);
b12ave = 0.5*(del[i][j][1][1] + del[i][j-1][1][1]);
areaaave = sqrt(pow(a11ave,2) + pow(a12ave,2));
areabave = sqrt(pow(b11ave,2) + pow(b12ave,2));
vdotna = (usol[i][j]*a11ave + vsol[i][j]*a12ave)/areaaave;
vdotnb = (usol[i][j]*b11ave + vsol[i][j]*b12ave)/areabave;
eiga = 0.5*areaaave*(abs(vdotna) + sqrt(pow(vdotna,2) + 4.0*beta2[i][j]));
eigb = 0.5*areabave*(abs(vdotnb) + sqrt(pow(vdotnb,2) + 4.0*beta2[i][j]));
voldt = eiga + eigb;
dt[i][j] = CFL*vol[i][j]/voldt;
dt_global=min(dt_global,dt[i][j]);
maxeiga=max(maxeiga,eiga);
maxeiga=max(maxeiga,eiga);
}
}

if(dt_crit==2) dt_global=CFL*mindelta_global/maxeiga;

//Ib addition
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

dtv = /*vol[i][j]/tstep;*/ vol[i][j]/dt[i][j];
res[i][j][0] = dtv*(psol[i][j]-pbc)/beta2[i][j];
dtv = vol[i][j]/tstep; // vol[i][j]/dt[i][j];
res[i][j][1] = dtv*rho*(usol[i][j]-ubc);
res[i][j][2] = dtv*rho*(vsol[i][j]-vbc);
}
}
} 

//end of ib addition


//--------------------------------------------------
//Checking the residue for termination
//--------------------------------------------------
double resnormx=0,resnormy=0;
double ratioc=0,ratiom=0,resnorm=0,resnormc=0;

for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
res[i][j][1]=res[i][j][1]+rho*vol[i][j]*(usol[i][j]-uold[i][j])/tstep;
res[i][j][2]=res[i][j][2]+rho*vol[i][j]*(vsol[i][j]-vold[i][j])/tstep;
resnormc=resnormc+res[i][j][0];
resnormx=resnormx+res[i][j][1]*res[i][j][1];
resnormy=resnormy+res[i][j][2]*res[i][j][2];
}

resnormc=resnormc/imax/jmax;
resnormx=resnormx/imax/jmax;
resnormy=resnormy/imax/jmax;
resnorm=sqrt(resnormx+resnormy);

if(iter==0)
{
if(parameter_num==0&&phy_time==0)
{
resnorm_ref=resnorm+1e-30;
resnormc_ref=resnormc+1e-30;
}
}

ratioc=abs(resnormc/(resnormc_ref));
ratiom=abs(resnorm/(resnorm_ref*visc[parameter_num]/visc[0]));

cout<<"For Reynolds number:"<<RE<<endl;
cout<<"Time:"<<phy_time<<endl;
cout<<"Time Step:"<<tstep<<endl;
cout<<"Iteration Number: "<<iter<<endl;
cout<<"Absolute Error in continuity:"<<resnormc<<endl;
cout<<"Reference Residue for continuity:"<<resnormc_ref<<endl;
cout<<"Relative Error in continuity:"<<ratioc<<endl;
cout<<"Absolute Error in X-Momentum:"<<sqrt(resnormx)<<endl;
cout<<"Absolute Error in Y-Momentum:"<<sqrt(resnormy)<<endl;
cout<<"Reference Residue for Momentum:"<<resnorm_ref*visc[parameter_num]/visc[0]<<endl;
cout<<"Relative Error in Momentum:"<<ratiom<<endl;
cout<<"-----------------------------------------------------"<<endl;
if(ratioc<tol_con && ratiom<tol_mom)
{
fres<<phy_time<<" "<<iter<<" "<<ratiom<<" "<<resnorm<<" "<<ratioc<<" "<<resnormc<<endl;
converged=1;
break;
}

//---------------------------------------------------
//Solution update
//---------------------------------------------------

for(int j=1;j<jmax;j++)
{
for(int i=1;i<imax;i++)
{
	double vdtl,vdtg,t11,t22;
        vdtl = vol[i][j]/dt[i][j];  //local 
        vdtg = vol[i][j]/tstep;       //physical
        t11 = ((1.0-0*hh[i][j])*vdtl + 0*hh[i][j]*vdtg)/beta2[i][j];
        t22 = ((1.0-0*hh[i][j])*vdtl + vdtg)*rho;

        deltap = -res[i][j][0]/t11;
        deltau = -res[i][j][1]/t22;
        deltav = -res[i][j][2]/t22;

        psol[i][j] = psol[i][j] + deltap;
        usol[i][j] = usol[i][j] + deltau;
        vsol[i][j] = vsol[i][j] + deltav;
}
}

//End of Solution update

fpoint<<phy_time<<" ";
for(int n=0;n<npoints;n++)
{
fpoint<<sqrt(pow(usol[ipoint[n]][jpoint[n]],2)+pow(vsol[ipoint[n]][jpoint[n]],2))<<" ";
}
fpoint<<endl;



if(solsave_flag==nsave)
{
	cout<<"Writing solution at t="<<phy_time<<endl;
	results();
	solsave_flag=0;
}

}

if(converged==0)
{
results();
cout<<"Solution did not converge with the time step:"<<tstep<<endl;
exit(1);
}

maxeiga=-1e10;
//tstep=1e10;
for(int j=1;j<jmax; j++)
{
for(int i=1;i<imax;i++)
{
a11ave = 0.5*(del[i][j][0][0] + del[i-1][j][0][0]);
a12ave = 0.5*(del[i][j][0][1] + del[i-1][j][0][1]);
b11ave = 0.5*(del[i][j][1][0] + del[i][j-1][1][0]);
b12ave = 0.5*(del[i][j][1][1] + del[i][j-1][1][1]);
areaaave = sqrt(pow(a11ave,2) + pow(a12ave,2));
areabave = sqrt(pow(b11ave,2) + pow(b12ave,2));
vdotna = (usol[i][j]*a11ave + vsol[i][j]*a12ave)/areaaave;
vdotnb = (usol[i][j]*b11ave + vsol[i][j]*b12ave)/areabave;
eiga = 0.5*areaaave*(abs(vdotna) + sqrt(pow(vdotna,2) + 4.0*beta2[i][j]));
eigb = 0.5*areabave*(abs(vdotnb) + sqrt(pow(vdotnb,2) + 4.0*beta2[i][j]));
voldt = eiga + eigb;
dt[i][j] = CFL*vol[i][j]/voldt;
//tstep=min(tstep,dt[i][j]);
maxeiga=max(maxeiga,eiga);
maxeiga=max(maxeiga,eiga);
}
}

//tstep=CFL*mindelta_global/maxeiga;
tstep=deltaT_store;
}
fres.close();
fpoint.close();

}
