clc
close all
clear all

fileIDlinear = fopen('results_linear_y.dat','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A50 = fscanf(fileIDlinear,formatSpec,sizeA);
fclose(fileIDlinear);
A50=A50';
xlinear=A50(:,1);
ylinear=A50(:,2);
ulinear=A50(:,3);
vlinear=A50(:,4);
umaglinear=A50(:,5);
plinear=A50(:,6);

fileIDtvd = fopen('results_TVD_y.dat','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A50 = fscanf(fileIDtvd,formatSpec,sizeA);
fclose(fileIDtvd);
A50=A50';
xtvd=A50(:,1);
ytvd=A50(:,2);
utvd=A50(:,3);
vtvd=A50(:,4);
umagtvd=A50(:,5);
ptvdr=A50(:,6);

fig=figure('Name',['Comparision of v-values along horizontal Line through Geometric Center of Cavity']);
title(['Comparision of v-values along horizontal line through geometric center of cavity for Re=250 and with square grid 50x50']);
ylabel('v velocity')
xlabel('x')
hold on
plot(xlinear,vlinear,'color','k','LineWidth',2);
hold on
plot(xtvd,vtvd,'--','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Linear ','TVD');
set(hleg1,'Location','NorthEast')

clear all
fileIDlinear = fopen('results_linear_x.dat','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A50 = fscanf(fileIDlinear,formatSpec,sizeA);
fclose(fileIDlinear);
A50=A50';
xlinear=A50(:,1);
ylinear=A50(:,2);
ulinear=A50(:,3);
vlinear=A50(:,4);
umaglinear=A50(:,5);
plinear=A50(:,6);

fileIDtvd = fopen('results_TVD_x.dat','r');
formatSpec = '%lf %lf %lf %lf %lf %lf';
sizeA = [6 Inf];
A50 = fscanf(fileIDtvd,formatSpec,sizeA);
fclose(fileIDtvd);
A50=A50';
xtvd=A50(:,1);
ytvd=A50(:,2);
utvd=A50(:,3);
vtvd=A50(:,4);
umagtvd=A50(:,5);
ptvdr=A50(:,6);

fig=figure('Name',['Comparision of u-values along vertical Line through Geometric Center of Cavity']);
title(['Comparision of u-values along vertical line through geometric center of cavity for Re=250 and with square grid 50x50']);
xlabel('u velocity')
ylabel('y')
hold on
plot(ulinear,ylinear,'color','k','LineWidth',2);
hold on
plot(utvd,ytvd,'--','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Linear ','TVD');
set(hleg1,'Location','NorthWest')

clear all

fileID = fopen('residue_linear.dat','r');
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
mass_ar=A50(:,5);


fileID = fopen('residue_tvd.dat','r');
formatSpec = '%lf %lf %lf %lf %lf';
sizeA = [5 Inf];
A50 = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A50=A50';
iter_tvd=A50(:,1);
mom_rr_tvd=A50(:,2);
mom_ar_tvd=A50(:,3);
mom_ar_tvd=log10(abs(mom_ar_tvd));
mass_rr_tvd=A50(:,4);
mass_ar_tvd=A50(:,5);

fig=figure('Name',['Logorithm of momemtum reative residue']);
title(['Logorithm of momentum relative residue for RE=250 and square grid 50x50']);
xlabel('Iteration No')
ylabel('log10(momentum relative residue)')
hold on
plot(iter,mom_rr,'color','k','LineWidth',2);
hold on
plot(iter_tvd,mom_rr_tvd,'-.','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Linear ','TVD');
set(hleg1,'Location','NorthEast')


fig=figure('Name',['Logorithm of momemtum absolute residue']);
title(['Logorithm of momentum absolute residue for RE=250 and square grid 50x50']);
xlabel('Iteration No')
ylabel('log10(momentum absolute residue)')
hold on
plot(iter,mom_ar,'color','k','LineWidth',2);
hold on
plot(iter_tvd,mom_ar_tvd,'-.','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Linear ','TVD');
set(hleg1,'Location','NorthEast')





