//Function which plots the results

#include "functions.h"

void results()
{
fstream fout;
string phifile;
ostringstream convert1;
convert1 <<phy_time;
phifile=mesh_name+"_"+convert1.str();
phifile=phifile+"_results_flow_field_RE_";
ostringstream convert2;
convert2 << RE;
phifile=phifile+convert2.str();
phifile=phifile+".dat";

fout.open(phifile.c_str(), ios::out | ios::trunc);

fout<<"Variables = x, y, u, v, umag, p"<<endl;
fout<<"ZONE I= "<<imax<<" J= "<<jmax<<" F=point"<<endl;

double unode,vnode,pnode;
for(int j=1;j<jmax+1;j++)
for(int i=1;i<imax+1;i++)
{
unode=0.25*(usol[i][j]+usol[i-1][j]+usol[i-1][j-1]+usol[i][j-1]);
vnode=0.25*(vsol[i][j]+vsol[i-1][j]+vsol[i-1][j-1]+vsol[i][j-1]);
pnode=0.25*(psol[i][j]+psol[i-1][j]+psol[i-1][j-1]+psol[i][j-1]);
fout<<x[i][j]<<" "<<y[i][j]<<" "<<unode<<" "<<vnode<<" "
	<<sqrt(pow(unode,2)+pow(vnode,2))<<" "<<pnode<<endl;
}
fout.close();

if(time_accurate==1)
{
phifile=mesh_name+"_RE_";
ostringstream convert;
convert << RE;
phifile=phifile+convert.str();
phifile=phifile+"_Restart";
phifile=phifile+".dat";

fout.open(phifile.c_str(), ios::out | ios::app);

for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
fout<<usol[i][j]<<" "<<vsol[i][j]<<" "<<psol[i][j]<<endl;
}

fout.close();
}
