//Program for solving Navier-Stokes equation

#include "functions.h"

void incomp_navier_stokes()
{

double area=0, xnx=0, xny=0;
double uave=0, vave=0, pave=0;
double delp=0,pleft=0, pright=0;
int ind_max, Nconv;
double a11ave,a12ave;
double b11ave,b12ave;
double areaaave,areabave;
double vdotna,vdotnb,eiga,eigb,voldt;
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

beta2=new double* [imax+1];
for (int i=0; i < imax+1; i++)
beta2[i]=new double[jmax+1];


res=new double** [imax];
for (int i=0; i < imax; i++) res[i]= new double* [jmax];
for(int i=0; i <imax; i++) for(int j=0; j<jmax; j++)
res[i][j]=new double[3];

dt=new double* [imax];
for (int i=1; i < imax; i++) dt[i]= new double[jmax];

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
if(restart==2)
{
fres.open(filename.c_str(), ios::out | ios::app);
}
else 
{
fres.open(filename.c_str(), ios::out | ios::trunc);
fres<<"Variables = Iteration, log_Relative_Residue_momentum, Residue_momentum, log_RR_cont, Residue_continuity "<<endl;
}

fstream frestart;
filename=mesh_name+"_RE_";
ostringstream convert3;
convert3 << RE;
filename=filename+convert3.str();
filename=filename+"_Restart";
filename=filename+".dat";
frestart.open(filename.c_str(), ios::out | ios::trunc);

//Iterations Started
for(int iter=0;iter<Niter;iter++)
{
boundary_conditions();

//Calculating Beta values
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

//End of Calculating the beta values	

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
maxeiga=max(maxeiga,abs(eiga));
maxeiga=max(maxeiga,abs(eigb));
}
}
if(dt_crit==2) dt_global=CFL*mindelta_global/maxeiga;

if(ibmethod==2)
{
cout<<"Immerssed Boundary method used"<<endl;
immersed_boundary();
}
//--------------------------------------------------
//Checking the residue for termination
//--------------------------------------------------
double resnormx=0,resnormy=0;
double ratioc=0,ratiom=0,resnorm=0,resnormc=0;

for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
resnormc=resnormc+res[i][j][0];
resnormx=resnormx+res[i][j][1]*res[i][j][1];
resnormy=resnormy+res[i][j][2]*res[i][j][2];
}

resnormc=resnormc/imax/jmax+1e-30;
resnormx=resnormx/imax/jmax;
resnormy=resnormy/imax/jmax;
resnorm=sqrt(resnormx+resnormy)+1e-30;

if(iter==0)
{
//Restarting
if(restart==2)
{
ifstream fin;
string line;
double log_RR_dummy;
fin.open("restart.dat",std::ios::in);

// Check for successful opening of file
if(fin.fail())
{
    cout<<"Error Opening File: restart.dat"<<endl;
    exit(1);
}

fin>>line>>resnormc_ref;
fin>>line>>resnorm_ref;
fin.close();

resnormc_ref=resnormc_ref+1e-30;
resnorm_ref=resnorm_ref+1e-30;
}
//End of restarting
//Starting from scratch
else
{
if(parameter_num==0)
{
resnorm_ref=resnorm;
resnormc_ref=resnormc;
}	
}
}

ratioc=abs(resnormc/(resnormc_ref));
//ratiom=abs(resnorm/(resnorm_ref*visc[parameter_num]/visc[0]));
ratiom=abs(resnorm/(resnorm_ref));

fres<<iter<<" "<<log10(ratiom)<<" "<<resnorm<<" "<<log10(ratioc)<<" "<<resnormc<<endl;

cout<<"For Reynolds number: "<<RE<<endl;
cout<<"Iteration Number: "<<iter<<endl;
cout<<"Absolute Error in continuity: "<<resnormc<<endl;
cout<<"Reference Residue for continuity: "<<resnormc_ref<<endl;
cout<<"Relative Error in continuity: "<<ratioc<<endl;
cout<<"Absolute Error in X-Momentum: "<<sqrt(resnormx)<<endl;
cout<<"Absolute Error in Y-Momentum: "<<sqrt(resnormy)<<endl;
//cout<<"Reference Residue for Momentum: "<<resnorm_ref*visc[parameter_num]/visc[0]<<endl;
cout<<"Reference Residue for Momentum: "<<resnorm_ref<<endl;
cout<<"Relative Error in Momentum: "<<ratiom<<endl;
cout<<"-----------------------------------------------------"<<endl;
if(ratioc<tol_con && ratiom<tol_mom)
{
frestart<<"Reference_residue_cont "<<resnormc_ref<<endl;
frestart<<"Reference_residue_mom "<<resnorm_ref<<endl;
break;
}
if(iter==Niter-1)
{
cout<<"Solution did not converge"<<endl;
cout<<"Maximum number of iterations reached, N="<<Niter<<endl;
frestart<<"Reference_residue_cont: "<<resnormc_ref<<endl;
frestart<<"Reference_residue_mom: "<<resnorm_ref<<endl;
frestart.close();
results();
exit(2);
}
//---------------------------------------------------
//Solution update
//---------------------------------------------------


if(dt_crit==3)
{
for(int j=1;j<jmax; j++)
for(int i=1;i<imax; i++)
{
double dtv,deltap,deltau,deltav;
dtv = dt[i][j]/vol[i][j];
deltap = -beta2[i][j]*dtv*res[i][j][0];
deltau = -dtv*res[i][j][1]/rho;
deltav = -dtv*res[i][j][2]/rho;
psol[i][j] = psol[i][j] + deltap;
usol[i][j] = usol[i][j] + deltau;
vsol[i][j] = vsol[i][j] + deltav;
}
}
else
{
for(int j=1;j<jmax; j++)
for(int i=1;i<imax; i++)
{
dtv = dt_global/vol[i][j];
deltap = -beta2[i][j]*dtv*res[i][j][0];
deltau = -dtv*res[i][j][1]/rho;
deltav = -dtv*res[i][j][2]/rho;
psol[i][j] = psol[i][j] + deltap;
usol[i][j] = usol[i][j] + deltau;
vsol[i][j] = vsol[i][j] + deltav;
}
}

//-------------------------------------------
// End of solution update
//------------------------------------------

}

results();

fres.close();
frestart.close();


}
