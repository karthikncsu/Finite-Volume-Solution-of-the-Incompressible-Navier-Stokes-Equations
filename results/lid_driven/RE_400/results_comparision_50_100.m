clc
close all
clear all

fileIDlinear = fopen('results_y.dat','r');
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

fileID = fopen('reference_vvel','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];
Auvel = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
Auvel=Auvel';
xvvel=Auvel(:,1);
vsolvvel=Auvel(:,2);


fig=figure('Name',['Comparision of v-values along horizontal Line through Geometric Center of Cavity']);
title(['Comparision of v-values along horizontal line through geometric center of cavity for Re=400']);
ylabel('v velocity')
xlabel('x')
hold on
plot(xlinear,vlinear,'color','k','LineWidth',2);
hold on
plot(xvvel,vsolvvel,'o','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('computation','Reference');
set(hleg1,'Location','NorthEast')

clear all
fileIDlinear = fopen('results_x.dat','r');
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

fileID = fopen('reference_uvel','r');
formatSpec = '%lf %lf';
sizeA = [2 Inf];
Auvel = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
Auvel=Auvel';
yuvel=Auvel(:,1);
usoluvel=Auvel(:,2);

fig=figure('Name',['Comparision of u-values along vertical Line through Geometric Center of Cavity']);
title(['Comparision of u-values along vertical line through geometric center of cavity for Re=400']);
xlabel('u velocity')
ylabel('y')
hold on
plot(ulinear,ylinear,'color','k','LineWidth',2);
hold on
plot(usoluvel,yuvel,'o','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('computation','Reference');
set(hleg1,'Location','NorthWest')