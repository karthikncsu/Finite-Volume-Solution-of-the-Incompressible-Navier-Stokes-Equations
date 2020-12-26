//Function for assigning Boundary conditions and source terms

#include "functions.h"

void boundary_conditions()
{

//---------------------------------------------
//Right and left boundary condition;
//---------------------------------------------
for(int j=1; j<jmax; j++) 
{
//---------------------------------------------
// Left boundary condition implementation
//---------------------------------------------
if(lb_type==0)
{
double ubx,uby,area;
flag[0][j]=0;
area=sqrt(pow(x[1][j+1]-x[1][j],2)+pow(y[1][j+1]-y[1][j],2));
ubx=lb_value*(x[1][j+1]-x[1][j])/area;
uby=lb_value*(y[1][j+1]-y[1][j])/area;
usol[0][j]=2*ubx-usol[1][j];
vsol[0][j]=2*uby-vsol[1][j];
psol[0][j]=psol[1][j];
}
else if(lb_type==1)
{
double nx,ny,area,vdotn;
flag[0][j]=0;
area=sqrt(pow(del[0][j][0][0],2)+pow(del[0][j][0][1],2));
nx=del[0][j][0][0]/area;
ny=del[0][j][0][1]/area;
vdotn=usol[1][j]*nx+vsol[1][j]*ny;
usol[0][j]=usol[1][j]-2*vdotn*nx;
vsol[0][j]=vsol[1][j]-2*vdotn*ny;
psol[0][j]=psol[1][j];
}

else if(lb_type==2)
{
if(lb_flow_cond==0)
{
double nx,ny,area;
area=sqrt(pow(del[0][j][0][0],2)+pow(del[0][j][0][1],2));
nx=del[0][j][0][0]/area;
ny=del[0][j][0][1]/area;
double tsize=4*tstep;
//double f=1-pow((phy_time-tsize)/tsize,2);
usol[0][j]=4*lb_value*((yc[1][j]-y[1][1])/(y[1][jmax]-y[1][1])
				  -pow((yc[1][j]-y[1][1])/(y[1][jmax]-y[1][1]),2))*nx;
vsol[0][j]=4*lb_value*((yc[1][j]-y[1][1])/(y[1][jmax]-y[1][1])
				  -pow((yc[1][j]-y[1][1])/(y[1][jmax]-y[1][1]),2))*ny;
psol[0][j]=psol[1][j];
}
else if(lb_flow_cond==1)
{
usol[0][j]=usol[1][j];
vsol[0][j]=vsol[1][j];
psol[0][j]=lb_value-0.5*rho*(pow(usol[0][j],2)+pow(vsol[0][j],2));
}
}

else if(lb_type==3)
{
usol[0][j]=usol[1][j];
vsol[0][j]=vsol[1][j];
psol[0][j]=lb_value;
}
//---------------------------------------------
// End of left boundary condition implementation
//---------------------------------------------

//---------------------------------------------
// Right boundary condition implementation
//---------------------------------------------
if(rb_type==0)
{
double ubx,uby,area;
flag[imax-1][j]=0;
area=sqrt(pow(x[imax][j+1]-x[imax][j],2)
		+pow(y[imax][j+1]-y[imax][j],2));
ubx=rb_value*(x[imax][j+1]-x[imax][j])/area;
uby=rb_value*(y[imax][j+1]-y[imax][j])/area;
usol[imax][j]=2*ubx-usol[imax-1][j];
vsol[imax][j]=2*uby-vsol[imax-1][j];
psol[imax][j]=psol[imax-1][j];
}
else if(rb_type==1)
{
double nx,ny,area,vdotn;
flag[imax-1][j]=0;
area=sqrt(pow(del[imax-1][j][0][0],2)+pow(del[imax-1][j][0][1],2));
nx=del[imax-1][j][0][0]/area;
ny=del[imax-1][j][0][1]/area;
vdotn=usol[imax-1][j]*nx+vsol[imax-1][j]*ny;
usol[imax][j]=usol[imax-1][j]-2*vdotn*nx;
vsol[imax][j]=vsol[imax-1][j]-2*vdotn*ny;
psol[imax][j]=psol[imax-1][j];
}

else if(rb_type==2)
{
if(rb_flow_cond==0)
{
double nx,ny,area;
area=sqrt(pow(del[imax-1][j][0][0],2)+pow(del[imax-1][j][0][1],2));
nx=del[imax-1][j][0][0]/area;
ny=del[imax-1][j][0][1]/area;
usol[imax][j]=4*rb_value*((yc[imax][j]-y[imax][1])/(y[imax][jmax]-y[imax][1])
       				 -pow((yc[imax][j]-y[imax][1])/(y[imax][jmax]-y[imax][1]),2))*nx;
vsol[imax][j]=4*rb_value*((yc[imax][j]-y[imax][1])/(y[imax][jmax]-y[imax][1])
				     -pow((yc[imax][j]-y[imax][1])/(y[imax][jmax]-y[imax][1]),2))*ny;
psol[imax][j]=psol[imax-1][j];
}
else if(rb_flow_cond==1)
{
usol[imax][j]=usol[imax-1][j];
vsol[imax][j]=vsol[imax-1][j];
psol[imax][j]=rb_value-0.5*rho*(pow(usol[imax-1][j],2)
							   +pow(vsol[imax-1][j],2));
}
}

else if(rb_type==3)
{
usol[imax][j]=usol[imax-1][j];
vsol[imax][j]=vsol[imax-1][j];
psol[imax][j]=rb_value;
}
//---------------------------------------------
// End of left boundary condition implementation
//---------------------------------------------

}
//----------------------------------------------
// End of Right and left boundary condition
//----------------------------------------------

//---------------------------------------------
//Top and bottom boundary condition;
//---------------------------------------------
for(int i=1; i<imax; i++) 
{
//---------------------------------------------
// Bottom boundary condition implementation
//---------------------------------------------
if(bb_type==0)
{
double ubx,uby,area;
flag[i][0]=0;
area=sqrt(pow(x[i+1][1]-x[i][1],2)+pow(y[i+1][1]-y[i][1],2));
ubx=bb_value*(x[i+1][1]-x[i][1])/area;
uby=bb_value*(y[i+1][1]-y[i][1])/area;
usol[i][0]=2*ubx-usol[i][1];
vsol[i][0]=2*uby-vsol[i][1];
psol[i][0]=psol[i][1];
}
else if(bb_type==1)
{
double nx,ny,area,vdotn;
flag[i][0]=0;
area=sqrt(pow(del[i][0][1][0],2)+pow(del[i][0][1][1],2));
nx=del[i][0][1][0]/area;
ny=del[i][0][1][1]/area;
vdotn=usol[i][1]*nx+vsol[i][1]*ny;
usol[i][0]=usol[i][1]-2*vdotn*nx;
vsol[i][0]=vsol[i][1]-2*vdotn*ny;
psol[i][0]=psol[i][1];
}
else if(bb_type==2)
{
if(bb_flow_cond==0)
{
double nx,ny,area;
area=sqrt(pow(del[i][0][1][0],2)+pow(del[i][0][1][1],2));
nx=del[i][0][1][0]/area;
ny=del[i][0][1][1]/area;
usol[i][0]=4*bb_value*((xc[i][1]-x[1][1])/(x[imax][1]-x[1][1])
				  -pow((xc[i][1]-x[1][1])/(x[imax][1]-x[1][1]),2))*nx;
vsol[i][0]=4*bb_value*((xc[i][1]-x[1][1])/(x[imax][1]-x[1][1])
				  -pow((xc[i][1]-x[1][1])/(x[imax][1]-x[1][1]),2))*ny;
psol[i][0]=psol[i][1];
}
else if(bb_flow_cond==1)
{
usol[i][0]=usol[i][1];
vsol[i][0]=vsol[i][1];
psol[i][0]=bb_value-0.5*rho*(pow(usol[i][0],2)+pow(vsol[i][0],2));
}
}

else if(bb_type==3)
{
usol[i][0]=usol[i][1];
vsol[i][0]=vsol[i][1];
psol[i][0]=bb_value;
}
//---------------------------------------------
// End of bottom boundary condition implementation
//---------------------------------------------

//---------------------------------------------
// Top boundary condition implementation
//---------------------------------------------
if(tb_type==0)
{
double ubx,uby,area;
flag[i][jmax-1]=0;
area=sqrt(pow(x[i+1][jmax]-x[i][jmax],2)
	 +pow(y[i+1][jmax]-y[i][jmax],2));
ubx=tb_value*(x[i+1][jmax]-x[i][jmax])/area;
uby=tb_value*(y[i+1][jmax]-y[i][jmax])/area;
usol[i][jmax]=2*ubx-usol[i][jmax-1];
vsol[i][jmax]=2*uby-vsol[i][jmax-1];
psol[i][jmax]=psol[i][jmax-1];
}
else if(tb_type==1)
{
double nx,ny,area,vdotn;
flag[i][jmax-1]=0;
area=sqrt(pow(del[i][jmax-1][1][0],2)+pow(del[i][jmax-1][1][1],2));
nx=del[i][jmax-1][1][0]/area;
ny=del[i][jmax-1][1][1]/area;
vdotn=usol[i][jmax-1]*nx+vsol[i][jmax-1]*ny;
usol[i][jmax]=usol[i][jmax-1]-2*vdotn*nx;
vsol[i][jmax]=vsol[i][jmax-1]-2*vdotn*ny;
psol[i][jmax]=psol[i][jmax-1];
}

else if(tb_type==2)
{
if(tb_flow_cond==0)
{
double nx,ny,area;
area=sqrt(pow(del[i][jmax-1][1][0],2)+pow(del[i][jmax-1][1][1],2));
nx=del[i][jmax-1][1][0]/area;
ny=del[i][jmax-1][1][1]/area;
usol[i][jmax]=4*tb_value*((xc[i][jmax]-x[1][jmax])/(x[imax][jmax]-x[1][jmax])
				     -pow((xc[i][jmax]-x[1][jmax])/(x[imax][jmax]-x[1][1]),2))*nx;
vsol[i][jmax]=4*tb_value*((xc[i][jmax]-x[1][jmax])/(x[imax][jmax]-x[1][jmax])
				     -pow((xc[i][jmax]-x[1][jmax])/(x[imax][jmax]-x[1][1]),2))*ny;
psol[i][jmax]=psol[i][jmax-1];
}
else if(tb_flow_cond==1)
{
usol[i][jmax]=usol[i][jmax-1];
vsol[i][jmax]=vsol[i][jmax-1];
psol[i][jmax]=tb_value-0.5*rho*(pow(usol[i][jmax],2)
							   +pow(vsol[i][jmax],2));
}
}

else if(tb_type==3)
{
usol[i][jmax]=usol[i][jmax-1];
vsol[i][jmax]=vsol[i][jmax-1];
psol[i][jmax]=tb_value;
}
//---------------------------------------------
// End of bottom boundary condition implementation
//---------------------------------------------

}
//----------------------------------------------
// End of top and bottom boundary condition
//----------------------------------------------

}
