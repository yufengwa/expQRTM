clc,clear;
nx=269;
nz=412;
dx=2.5;
dz=2.5;
dt=0.025;
nt=3201;


dk=2*pi*1/(dz*nz);

it=1:nt;
%plot snapshots trace 
fid1=fopen('../data/acofullns/Final_image_cor_type0.dat','rb');
fid2=fopen('../data/attfullns/Final_image_cor_type1.dat','rb');
fid3=fopen('../data/compfullns/Final_image_cor_type2.dat','rb');
fid4=fopen('../data/stafullns/Final_image_cor_type5.dat','rb');

ima1=fread(fid1,[nz,nx],'float');
ima2=fread(fid2,[nz,nx],'float');
ima3=fread(fid3,[nz,nx],'float');
ima4=fread(fid4,[nz,nx],'float');



%% reftrace = 60

figure;

Reftrace=60;

trace1=ima1(:,Reftrace);
trace2=ima2(:,Reftrace);
trace3=ima3(:,Reftrace);
trace4=ima4(:,Reftrace);


xi=1:0.1:length(trace1);

trace11=interp1(trace1,xi,'spline');
trace22=interp1(trace2,xi,'spline');
trace33=interp1(trace3,xi,'spline')/1e6;
trace44=interp1(trace4,xi,'spline');


plot(trace11, 'k', 'Linewidth', 2); hold on;
plot(trace22, 'r','Linewidth', 1); hold on;
plot(trace33, 'b', 'Linewidth', 0.5); hold on;
plot(trace44, 'g', 'Linewidth', 0.5); hold on;

h = legend('Reference', 'Attenuated', ...
    'Compensated','Stabilized');

set(h,'fontsize',10);
set(h,'box','off');


axis([0 10*nz -2.5e10 2.5e10]);
xlabel('Depth (ft)', 'fontsize',12);
ylabel('Amplitude', 'fontsize',12);
set(gca,'Xtick',[0:500:4200],'Ytick',[-2.5e10:1e10:2.5e10],'fontsize',10);
set(gca,'Xticklabel',{[0:125:1050]},'fontsize',10);
set(gca,'YAxisLocation','right');
set(gcf,'Position',[100 100 300 800]); 
set(gca,'Position',[.2 .05 .78 .88]); 

view([90 90]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 8.0]);
print ./Fig/cwtrace1.eps -depsc2 -r600;


klength=200;
figure;
spec1=abs((fft(trace1)));
spec2=abs((fft(trace2)));
spec3=abs((fft(trace3)));
spec4=abs((fft(trace4)));

specmax=max([max(spec1),max(spec4)]);

plot(spec1/specmax, 'k', 'Linewidth', 1); hold on;
plot(spec2/specmax, 'r','Linewidth', 1); hold on;
plot(spec3/max(spec3), 'b', 'Linewidth', 0.5); hold on;
plot(spec4/specmax, 'g', 'Linewidth', 0.5); hold on;

axis([0 200 0 1.0]);
xlabel('Wavenumber ({ft}^{-1})', 'fontsize',12);
ylabel('Normalized amplitude', 'fontsize',12);
set(gca,'Xtick',[0:50:200],'Ytick',[0:0.2:1.0],'fontsize',10);
set(gca,'Xticklabel',{[0:0.3:1.2]},'fontsize',10);
%grid on;


h = legend('Reference', 'Attenuated', ...
    'Compensated','Stabilized');

set(h,'fontsize',7);
set(h,'box','off');


set(gcf,'Position',[100 100 400 280]); 
set(gca,'Position',[.15 .16 .80 .80]); 

set(gcf,'paperpositionmode','auto');
print ./Fig/cwspec1.eps -depsc2 -r600;



%% reftrace=150

figure;

Reftrace=150;

trace1=ima1(:,Reftrace);
trace2=ima2(:,Reftrace);
trace3=ima3(:,Reftrace);
trace4=ima4(:,Reftrace);


xi=1:0.1:length(trace1);

trace11=interp1(trace1,xi,'spline');
trace22=interp1(trace2,xi,'spline');
trace33=interp1(trace3,xi,'spline')/1e6;
trace44=interp1(trace4,xi,'spline');


plot(trace11, 'k', 'Linewidth', 2); hold on;
plot(trace22, 'r','Linewidth', 1); hold on;
plot(trace33, 'b', 'Linewidth', 0.5); hold on;
plot(trace44, 'g', 'Linewidth', 0.5); hold on;


axis([0 10*nz -2.5e10 2.5e10]);
xlabel('Depth (ft)', 'fontsize',12);
ylabel('Amplitude', 'fontsize',12);
set(gca,'Xtick',[0:500:4200],'Ytick',[-2.5e10:1e10:2.5e10],'fontsize',10);
set(gca,'Xticklabel',{[0:125:1050]},'fontsize',10);
set(gca,'YAxisLocation','right');
set(gcf,'Position',[100 100 300 800]); 
set(gca,'Position',[.2 .05 .78 .88]); 

view([90 90]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 8.0]);
print ./Fig/cwtrace2.eps -depsc2 -r600;

klength=200;
figure;
spec1=abs((fft(trace1)));
spec2=abs((fft(trace2)));
spec3=abs((fft(trace3)));
spec4=abs((fft(trace4)));

specmax=max([max(spec1),max(spec4)]);

plot(spec1/specmax, 'k', 'Linewidth', 1); hold on;
plot(spec2/specmax, 'r','Linewidth', 1); hold on;
plot(spec3/max(spec3), 'b', 'Linewidth', 0.5); hold on;
plot(spec4/specmax, 'g', 'Linewidth', 0.5); hold on;

axis([0 200 0 1.0]);
xlabel('Wavenumber ({ft}^{-1})', 'fontsize',12);
ylabel('Normalized amplitude', 'fontsize',12);
set(gca,'Xtick',[0:50:200],'Ytick',[0:0.2:1.0],'fontsize',10);
set(gca,'Xticklabel',{[0:0.3:1.2]},'fontsize',10);
%grid on;


set(gcf,'Position',[100 100 400 280]); 
set(gca,'Position',[.15 .16 .80 .80]); 

set(gcf,'paperpositionmode','auto');
print ./Fig/cwspec2.eps -depsc2 -r600;

%% reftrace=220
figure;

Reftrace=220;

trace1=ima1(:,Reftrace);
trace2=ima2(:,Reftrace);
trace3=ima3(:,Reftrace);
trace4=ima4(:,Reftrace);


xi=1:0.1:length(trace1);

trace11=interp1(trace1,xi,'spline');
trace22=interp1(trace2,xi,'spline');
trace33=interp1(trace3,xi,'spline')/1e7;
trace44=interp1(trace4,xi,'spline');


plot(trace11, 'k', 'Linewidth', 2); hold on;
plot(trace22, 'r','Linewidth', 1); hold on;
plot(trace33, 'b', 'Linewidth', 0.5); hold on;
plot(trace44, 'g', 'Linewidth', 0.5); hold on;


axis([0 10*nz -2.5e10 2.5e10]);
xlabel('Depth (ft)', 'fontsize',12);
ylabel('Amplitude', 'fontsize',12);
set(gca,'Xtick',[0:500:4200],'Ytick',[-2.5e10:1e10:2.5e10],'fontsize',10);
set(gca,'Xticklabel',{[0:125:1050]},'fontsize',10);
set(gca,'YAxisLocation','right');
set(gcf,'Position',[100 100 300 800]); 
set(gca,'Position',[.2 .05 .78 .88]); 

view([90 90]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 8.0]);
print ./Fig/cwtrace3.eps -depsc2 -r600;

klength=200;
figure;
spec1=abs((fft(trace1)));
spec2=abs((fft(trace2)));
spec3=abs((fft(trace3)));
spec4=abs((fft(trace4)));

specmax=max([max(spec1),max(spec4)]);

plot(spec1/specmax, 'k', 'Linewidth', 1); hold on;
plot(spec2/specmax, 'r','Linewidth', 1); hold on;
plot(spec3/max(spec3), 'b', 'Linewidth', 0.5); hold on;
plot(spec4/specmax, 'g', 'Linewidth', 0.5); hold on;

axis([0 200 0 1.0]);
xlabel('Wavenumber ({ft}^{-1})', 'fontsize',12);
ylabel('Normalized amplitude', 'fontsize',12);
set(gca,'Xtick',[0:50:200],'Ytick',[0:0.2:1.0],'fontsize',10);
set(gca,'Xticklabel',{[0:0.3:1.2]},'fontsize',10);
%grid on;


set(gcf,'Position',[100 100 400 280]); 
set(gca,'Position',[.15 .16 .80 .80]); 

set(gcf,'paperpositionmode','auto');
print ./Fig/cwspec3.eps -depsc2 -r600;

