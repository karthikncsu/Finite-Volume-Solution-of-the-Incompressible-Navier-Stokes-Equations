fileID = fopen('grid3.dat_residue_RE_100.dat','r');
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
title(['Logorithm of momentum and mass absolute residue for RE=100']);
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
title(['Logorithm of momentum and mass relative residue for RE=100']);
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
