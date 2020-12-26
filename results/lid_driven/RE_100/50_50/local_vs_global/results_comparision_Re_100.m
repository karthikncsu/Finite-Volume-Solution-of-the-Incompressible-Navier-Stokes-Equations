clc
close all
clear all

fileID = fopen('re_100_50_50_local.dat','r');
formatSpec = '%lf %lf %lf %lf %lf';
sizeA = [5 Inf];
A50 = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A50=A50'
iter_tvd=A50(:,1);
mom_rr_tvd=A50(:,2);
mom_ar_tvd=A50(:,3);
mom_ar_tvd=log10(abs(mom_ar_tvd));
mass_rr_tvd=A50(:,4);
mass_ar_tvd=A50(:,5);

fileID = fopen('global_file.txt','r');
formatSpec = '%lf %lf %lf %lf %lf';
sizeA = [5 Inf];
A60 = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A60=A60';
iter_tvdtg=A60(:,1);
mom_rr=A60(:,2);
mom_ar=A60(:,3);
mom_ar=log10(abs(mom_ar));
mass_rr=A60(:,4);
mass_ar=A60(:,5);

fig=figure('Name',['Logorithm of momemtum reative residue']);
title(['Logorithm of momentum relative residue for RE=100 and square grid 50x50']);
xlabel('Iteration No')
ylabel('log10(momentum relative residue)')
hold on
plot(iter_tvdtg,mom_rr,'color','k','LineWidth',2);
hold on
plot(iter_tvd,mom_rr_tvd,'--','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('gloabl ','local');
set(hleg1,'Location','NorthEast')


fig=figure('Name',['Logorithm of momemtum absolute residue']);
title(['Logorithm of momentum absolute residue for RE=100 and square grid 50x50']);
xlabel('Iteration No')
ylabel('log10(momentum absolute residue)')
hold on
plot(iter_tvdtg,mom_ar,'color','k','LineWidth',2);
hold on
plot(iter_tvd,mom_ar_tvd,'--','color','k','LineWidth',2);
hold on
grid on
hold on
hleg1 = legend('Global ','Local');
set(hleg1,'Location','NorthEast')





