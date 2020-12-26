#ifndef HEADER_H
#define HEADER_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cmath>

using namespace std;

extern string mesh_name; //Mesh file name
extern int imax; // Number of mesh nodes along x direction
extern int jmax; // Number of mesh nodes along y direction
extern double **x,**y; //X and Y coordinates
extern double **xc,**yc; //X and Y coordinates of cente of the cell
extern double ****del; // Node numbers of each element
extern double **vol; // Volumes of Cells
extern double **usol, **vsol, **psol; //Velocity and pressure.
extern double rho; //Material Property: Density
extern int viscparas; // Number of viscous parameters to solve for
extern double *visc; //Material Propert: Viscosity
extern double ***res; //Residue
extern int TVD_method; //Variable indicating the accuracy using TVD
extern double **flag; //flag for boundary condition
extern int fv_method; //Variable indicating finite volume method used
extern int rb_type; //Boundary condtion type of right boundary
extern int tb_type; //Boundary condtion type of top boundary
extern int lb_type; //Boundary condtion type of left boundary
extern int bb_type; //Boundary condtion type of bottom boundary
extern int rb_flow_cond; //Flow condtion type of right boundary
extern int tb_flow_cond; //Flow condtion type of top boundary
extern int lb_flow_cond; //Flow condtion type of left boundary
extern int bb_flow_cond; //Flow condtion type of bottom boundary
extern double rb_value; //Value of condition for right boundary
extern double tb_value; //Value of condition for top boundary
extern double lb_value; //Value of condition for left boundary
extern double bb_value; //Value of condition for bottom boundary
extern int Niter;
extern double tol_con; //Tolerance for continuity equation
extern double tol_mom; //Tolerance for momentum equation
extern double *CFLnums;
extern double CFL; //CFL Number
extern double Uref; //Reference velocity
extern double Lref; //Reference Length
extern double dt_crit; //Time Step Criteria 
extern int restart;
extern int ibmethod; //2: Immeressed boundary method else not
extern int time_accurate; //2: time accurate else not
extern double tstart; //start time
extern double tstep; //Time step
extern double tend; //Final time
extern int nsave; //Solution output_time_step
extern int npoints; //Number of points to monitor the results
extern double *xpoint,*ypoint;
extern int *ipoint,*jpoint;
extern double resnorm_ref,resnormc_ref;
extern int parameter_num;
extern double phy_time; //Physical Time
extern double RE; //Reynolds Number
extern double **dt; //local time step;
extern double **beta2;
extern double dt_global;



//Variables used in immerssed boundary method
extern int nobjs;
extern int *nlist;
extern int *idif;
extern int *jdif;
extern double **xs;
extern double **ys;
extern double **xnxs;
extern double **xnys;
extern double ***dist;
extern double **distg;
extern double **hh;
extern double ***weight;
extern double **acoef;
extern double **bcoef;
extern int ***itag;
extern int **lpri;
extern int nintpts;
extern double power;
//Variables used in immerssed boundary method



void read_input();
void reading_grid(string);
double triangle_area(double,double,double,double,double,double);
void initial_conditions();
void timeaccurate_function();
void incomp_navier_stokes();
void boundary_conditions();
void thin_layer();
void tangent_normal();
double min(double,double);
double max(double,double);
void results();
double sign(double,double);
void reading_objects();
void cell_classification();
void immersed_boundary();
void timeaccurate_ib_method();


#endif







