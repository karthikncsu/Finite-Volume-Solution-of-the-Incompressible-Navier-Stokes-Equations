// Function to read the input file

#include "functions.h"

void read_input()
{
string line;
ifstream fin;
fstream fout;
fout.open("input replicate", ios::out | ios::trunc);

fin.open("input.dat",std::ios::in);

if(fin.fail())
{
cout<<"Error Opening File: "<<"input.dat"<<endl;
exit(1);
}
getline(fin,line);
fin>>mesh_name;
fout<<mesh_name<<endl;
getline(fin,line);
fin>>fv_method;
fout<<fv_method<<endl;
getline(fin,line);
fin>>rb_type;
fout<<rb_type<<endl;
getline(fin,line);
fin>>rb_flow_cond;
fout<<rb_flow_cond<<endl;
getline(fin,line);
fin>>rb_value;
fout<<rb_value<<endl;
getline(fin,line);
fin>>tb_type;
fout<<tb_type<<endl;
getline(fin,line);
fin>>tb_flow_cond;
fout<<tb_flow_cond<<endl;
getline(fin,line);
fin>>tb_value;
fout<<tb_value<<endl;
getline(fin,line);
fin>>lb_type;
fout<<lb_value<<endl;
getline(fin,line);
fin>>lb_flow_cond;
fout<<lb_flow_cond<<endl;
getline(fin,line);
fin>>lb_value;
fout<<lb_value<<endl;
getline(fin,line);
fin>>bb_type;
fout<<bb_type<<endl;
getline(fin,line);
fin>>bb_flow_cond;
fout<<bb_flow_cond<<endl;
getline(fin,line);
fin>>bb_value;
fout<<bb_value<<endl;
getline(fin,line);
fin>>Niter;
fout<<Niter<<endl;
getline(fin,line);
fin>>tol_con;
fout<<tol_con<<endl;
getline(fin,line);
fin>>tol_mom;
fout<<tol_mom<<endl;
getline(fin,line);
fin>>rho;
fout<<rho<<endl;
getline(fin,line);
fin>>viscparas;
fout<<viscparas<<endl;
getline(fin,line);

visc=new double[viscparas];

for( int i=0;i<viscparas; i++)
{
fin>>visc[i];
fout<<visc[i]<<" ";
}
fout<<endl;
getline(fin,line);
fin>>TVD_method;
fout<<TVD_method<<endl;
getline(fin,line);

CFLnums=new double[viscparas];

for( int i=0;i<viscparas; i++)
{
fin>>CFLnums[i];
fout<<CFLnums[i]<<" ";
}
fout<<endl;

getline(fin,line);
fin>>Uref;
fout<<Uref<<endl;
getline(fin,line);
fin>>Lref;
fout<<Lref<<endl;
getline(fin,line);
fin>>dt_crit;
fout<<dt_crit<<endl;
getline(fin,line);
fin>>restart;
fout<<restart<<endl;
getline(fin,line);
fin>>ibmethod;
fout<<ibmethod<<endl;
getline(fin,line);
fin>>time_accurate;
fout<<time_accurate<<endl;
getline(fin,line);
fin>>tstart;
fout<<tstart<<endl;
getline(fin,line);
fin>>tstep;
fout<<tstep<<endl;
getline(fin,line);
fin>>tend;
fout<<tend<<endl;
getline(fin,line);
fin>>nsave;
fout<<nsave<<endl;
getline(fin,line);
fin>>npoints;
fout<<npoints<<endl;
getline(fin,line);

xpoint=new double[npoints];
ypoint=new double[npoints];

for( int i=0;i<npoints;i++)
{
fin>>xpoint[i]>>ypoint[i];
fout<<xpoint[i]<<" "<<ypoint[i]<<endl;
}


fin.close();
/*
//Printing the Input File
cout<<"mesh_name :"<<mesh_name<<endl;
cout<<"fv_method :"<<fv_method<<endl;
cout<<"rb_type :"<<rb_type<<endl;
cout<<"rb_flow_cond :"<<rb_flow_cond<<endl;
cout<<"rb_value :"<<rb_value<<endl;
cout<<"tb_type :"<<tb_type<<endl;
cout<<"tb_flow_cond :"<<tb_flow_cond<<endl;
cout<<"tb_value :"<<tb_value<<endl;
cout<<"lb_type :"<<lb_type<<endl;
cout<<"lb_flow_cond :"<<lb_flow_cond<<endl;
cout<<"lb_value :"<<lb_value<<endl;
cout<<"bb_type :"<<bb_type<<endl;
cout<<"bb_flow_cond :"<<bb_flow_cond<<endl;
cout<<"bb_value :"<<bb_value<<endl;
cout<<"Niter :"<<Niter<<endl;
cout<<"tol_con :"<<tol_con<<endl;
cout<<"tol_mom :"<<tol_mom<<endl;
cout<<"rho :"<<rho<<endl;
cout<<"visc :"<<visc<<endl;
cout<<"TVD_method :"<<TVD_method<<endl;
cout<<"CFL :"<<CFL<<endl;
cout<<"Uref :"<<Uref<<endl;
cout<<"dt_crit :"<<dt_crit<<endl;
*/


}
