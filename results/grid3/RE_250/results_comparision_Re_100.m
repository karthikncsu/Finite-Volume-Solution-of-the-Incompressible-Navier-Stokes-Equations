clc
close all
clear all

fileIDlinear = fopen('results_x_3.5','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A50 = fscanf(fileIDlinear,formatSpec,sizeA);
fclose(fileIDlinear);
A50=A50';
x35=A50(:,1);
y35=A50(:,2);
u35=A50(:,3);
v35=A50(:,4);
umag35=A50(:,5);
p35=A50(:,6);

fileIDlinear = fopen('results_x_4.5','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A50 = fscanf(fileIDlinear,formatSpec,sizeA);
fclose(fileIDlinear);
A50=A50';
x45=A50(:,1);
y45=A50(:,2);
u45=A50(:,3);
v45=A50(:,4);
umag45=A50(:,5);
p45=A50(:,6);

fileIDlinear = fopen('results_x_5','r');




fileIDlinear = fopen('results_x_6.0','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A60 = fscanf(fileIDlinear,formatSpec,sizeA);
fclose(fileIDlinear);
A60=A60';
x60=A60(:,1);
y60=A60(:,2);
u60=A60(:,3);
v60=A60(:,4);
umag60=A60(:,5);
p60=A60(:,6);

fig=figure('Name',['Comparision of u-values along verticle line']);
title(['Comparision of u-values at various vertical lines for Re=250 with grid3']);
xlabel('u velocity')
ylabel('y')
hold on
plot(u35,y35,'color','k','LineWidth',2);
hold on
%plot(u40,y40,'LineWidth',2);
%hold on
plot(u45,y45,'--','color','k','LineWidth',2);
hold on
%plot(u50,y50,'LineWidth',2);
%hold on
%plot(u55,y55,'LineWidth',2);
%hold on
plot(u60,y60,'-o','color','k','LineWidth',2);
hold on

grid on
hold on
hleg1 = legend('x=3.5','x=4.5','x=6.0');
set(hleg1,'Location','NorthEast')

fig=figure('Name',['Comparision of v-values along verticle line']);
title(['Comparision of v-values at various vertical lines for Re=250 with grid3']);
xlabel('v velocity')
ylabel('y')
hold on
plot(v35,y35,'color','k','LineWidth',2);
hold on
plot(v45,y45,'--','color','k','LineWidth',2);
hold on
plot(v60,y60,'-o','color','k','LineWidth',2);
hold on

grid on
hold on
hleg1 = legend('x=3.5','x=4.5','x=6.0');
set(hleg1,'Location','NorthEast')

clear all

fileID = fopen('grid3_re_250_residue.dat','r');
formatSpec = '%lf %lf %lf %lf %lf';
sizeA = [5 Inf];
A50 = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A50=A50';
iter=A50(:,1);
mom_rr=A50(:,2);
mom_ar=A50(:,3);
mom_ar=log10(abs(mom_ar));
mass_rr=A50(:,4);
mass_ar=log10(abs(A50(:,5)));


fig=figure('Name',['Logorithm of momemtum and mass absolute residue']);
title(['Logorithm of momentum and mass absolute residue for RE=250']);
xlabel('Iteration No')
ylabel('log10(absolute residue)')
hold on
plot(iter,mom_ar,'color','k','LineWidth',2);
hold on
plot(iter,mass_ar,'--','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Momentum Residue','Mass Residue');
set(hleg1,'Location','NorthEast')

fig=figure('Name',['Logorithm of momemtum and mass relative residue']);
title(['Logorithm of momentum and mass relative residue for RE=250']);
xlabel('Iteration No')
ylabel('log10(relative residue)')
hold on
plot(iter,mom_rr,'color','k','LineWidth',2);
hold on
plot(iter,mass_rr,'--','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Momentum Residue','Mass Residue');
set(hleg1,'Location','NorthEast')


