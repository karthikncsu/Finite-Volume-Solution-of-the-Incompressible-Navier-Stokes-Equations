//Program for solving Poission equation using different Finite Volume Methods

#include "functions.h"
#include <time.h>
#include <sys/time.h>

string mesh_name; //Mesh file name
int imax; // Number of mesh nodes along x direction
int jmax; // Number of mesh nodes along y direction
double **x,**y; //X and Y coordinates
double **xc,**yc; //X and Y coordinates of cente of the cell
double ****del; // Node numbers of each element
double **vol; // Volumes of Cells
double **usol, **vsol, **psol; //Velocity and pressure
double rho; //Material Property: Density
int viscparas; // Number of viscous parameters to solve for
double *visc; //Material Propert: Viscosity
double ***res; //Residue
int TVD_method; //Variable indicating the accuracy using TVD
double **flag; //flag for boundary condition
int fv_method; //Variable indicating finite volume method used
int rb_type; //Boundary condtion type of right boundary
int tb_type; //Boundary condtion type of top boundary
int lb_type; //Boundary condtion type of left boundary
int bb_type; //Boundary condtion type of bottom boundary
int rb_flow_cond; //Flow condtion type of right boundary
int tb_flow_cond; //Flow condtion type of top boundary
int lb_flow_cond; //Flow condtion type of left boundary
int bb_flow_cond; //Flow condtion type of bottom boundary
double rb_value; //Value of boundary condition for bottom boundary
double tb_value; //Value of boundary condition for bottom boundary
double lb_value; //Value of boundary condition for bottom boundary
double bb_value; //Value of boundary condition for bottom boundary
int Niter;
double tol_con; //Tolerance for continuity equation
double tol_mom; //Tolerance for momentum equation
double *CFLnums;
double CFL; //CFL Number
double Uref; //Reference velocity
double Lref; //Reference Length
double dt_crit; //Time Step Criteria
int restart;
int ibmethod; //2: Immeressed boundary method else not
int time_accurate; //2: time accurate else not
double tstart; //start time
double tstep; //Time step
double tend; //Final time
int nsave; //Solution output_time_step
int npoints; //Number of points to monitor the results
double *xpoint,*ypoint; //Indices of the npoints
int *ipoint,*jpoint;
double resnorm_ref,resnormc_ref;
int parameter_num; 
double phy_time; //Physical Time
double RE; //Reynolds Number
double **dt; //local time step;
double **beta2;
double dt_global;



//Variables used in immerssed boundary method
int nobjs;
int *nlist;
int *idif;
int *jdif;
double **xs;
double **ys;
double **xnxs;
double **xnys;
double ***dist;
double **distg;
double **hh;
double ***weight;
double **acoef;
double **bcoef;
int ***itag;
int **lpri;
int nintpts;
double power;
//Variables used in immerssed boundary method



int main()
{
double startcputime, endcputime,cpu_time;
string filename;
fstream fout;

cout<<"chek1"<<endl;
read_input();
cout<<"check2"<<endl;
cout<<"check2"<<endl;
reading_grid(mesh_name);
cout<<"check4";
cout<<"check3";
initial_conditions();
cout<<"chek4";

if(ibmethod==2)
{
cout<<"entered ibmethod"<<endl;
reading_objects();
cell_classification();
}

if(time_accurate==2&&ibmethod==1)
{
for(parameter_num=0;parameter_num<viscparas;parameter_num++)
{
startcputime = (double)clock();
RE=rho*Uref*Lref/(visc[parameter_num]+1e-16);
CFL=CFLnums[parameter_num];
timeaccurate_function();
endcputime=(double)clock();
cpu_time=(endcputime-startcputime)/CLOCKS_PER_SEC;
filename=mesh_name+"_cputime_RE_";
ostringstream convert;
convert << RE;
filename=filename+convert.str();

fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"simulation using Reynolds number: "<<RE<<"using time accurate method"<<endl;
fout<<"Solution time: "<<cpu_time; 
fout.close();
}
}
if(time_accurate==1)
{
phy_time=100;
for(parameter_num=0;parameter_num<viscparas;parameter_num++)
{
startcputime = (double)clock();
RE=rho*Uref*Lref/(visc[parameter_num]+1e-16);
CFL=CFLnums[parameter_num];
incomp_navier_stokes();
endcputime=(double)clock();
cpu_time=(endcputime-startcputime)/CLOCKS_PER_SEC;
filename=mesh_name+"_cputime_RE_";
ostringstream convert;
convert << RE;
filename=filename+convert.str();
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"simulation using Reynolds number: "<<RE<<" using steady state method and with CFL number: "<<CFL<<endl;

fout<<"Solution time: "<<cpu_time; 
fout.close();
}
}
if(time_accurate==2&&ibmethod==2)
{
for(parameter_num=0;parameter_num<viscparas;parameter_num++)
{
startcputime = (double)clock();
RE=rho*Uref*Lref/(visc[parameter_num]+1e-16);
CFL=CFLnums[parameter_num];
timeaccurate_ib_method();
endcputime=(double)clock();
cpu_time=(endcputime-startcputime)/CLOCKS_PER_SEC;
filename=mesh_name+"_cputime_RE_";
ostringstream convert;
convert << RE;
filename=filename+convert.str();

fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"simulation using Reynolds number: "<<RE<<"using time accurate method"<<endl;
fout<<"Solution time: "<<cpu_time; 
fout.close();
}
}

return(0);

}



