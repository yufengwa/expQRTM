clc,clear,close all;
nx=663;
nz=234;
dx=10;
dz=10;
dt=0.001;
nt=2001;

dk=2*pi*1/(dz*nz);

it=1:nt;
%plot snapshots trace 
fid1=fopen('../data/Final_image_cor_type0.dat','rb');
fid2=fopen('../data/Final_image_cor_type3.dat','rb');
fid3=fopen('../data/Final_image_cor_type4.dat','rb');
fid4=fopen('../data/Final_image_cor_type5.dat','rb');


ima1=fread(fid1,[nz,nx],'float');
ima2=fread(fid2,[nz,nx],'float');
ima3=fread(fid3,[nz,nx],'float');
ima4=fread(fid4,[nz,nx],'float');

%% reftrace=100
figure;

Reftrace=100;

trace1=ima1(:,Reftrace);
trace2=ima2(:,Reftrace);
trace3=ima3(:,Reftrace);
trace4=ima4(:,Reftrace);

xi=1:0.1:length(trace1);

trace11=interp1(trace1,xi,'spline');
trace22=interp1(trace2,xi,'spline');
trace33=interp1(trace3,xi,'spline');
trace44=interp1(trace4,xi,'spline');


plot(trace11, 'k', 'Linewidth', 1); hold on;
plot(trace22, 'r','Linewidth', 0.5); hold on;
plot(trace33, 'b', 'Linewidth', 0.5); hold on;
plot(trace44, 'g', 'Linewidth', 0.5); hold on;


h = legend('Reference', 'Stabilized (\alpha=1)', ...
    'Stabilized (\alpha=2)', ...
    'Stabilized (\alpha=8)');

set(h,'fontsize',8);
set(h,'box','off');

axis([0 10*nz -8e5 8e5]);
xlabel('Depth (m)', 'fontsize',12);
ylabel('Amplitude', 'fontsize',12);
set(gca,'Xtick',[0:500:2500],'Ytick',[-6e5:3e5:6e6],'fontsize',10);
set(gca,'Xticklabel',{[0:500:2000]},'fontsize',10);
set(gca,'YAxisLocation','right');
set(gcf,'Position',[100 100 300 600]); 
set(gca,'Position',[.20 .08 .78 .80]); 

view([90 90]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 6.0]);
print ./Fig/martrace1.eps -depsc2 -r600;


klength=100;
figure;
spec1=abs((fft(trace1)));
spec2=abs((fft(trace2)));
spec3=abs((fft(trace3)));
spec4=abs((fft(trace4)));

specmax=max([max(spec1),max(spec4)]);

plot(spec1/specmax, 'k', 'Linewidth', 1); hold on;
plot(spec2/specmax, 'r','Linewidth', 1); hold on;
plot(spec3/specmax, 'b', 'Linewidth', 0.5); hold on;
plot(spec4/specmax, 'g', 'Linewidth', 0.5); hold on;

axis([0 100 0 1.0]);
xlabel('Wavenumber ({m}^{-1})', 'fontsize',12);
ylabel('Normalized amplitude', 'fontsize',12);
set(gca,'Xtick',[0:20:100],'Ytick',[0:0.2:1.0],'fontsize',10);
set(gca,'Xticklabel',{[0:0.05:0.25]},'fontsize',10);
%grid on;

h = legend('Reference', 'Stabilized (\alpha=1)', ...
    'Stabilized (\alpha=2)', ...
    'Stabilized (\alpha=8)');

set(h,'fontsize',11);
set(h,'box','off');


set(gcf,'Position',[100 100 400 280]); 
set(gca,'Position',[.15 .16 .80 .80]); 

set(gcf,'paperpositionmode','auto');
print ./Fig/marspec1.eps -depsc2 -r600;


%% reftrace=320
figure;

Reftrace=320;

trace1=ima1(:,Reftrace);
trace2=ima2(:,Reftrace);
trace3=ima3(:,Reftrace);
trace4=ima4(:,Reftrace);

xi=1:0.1:length(trace1);

trace11=interp1(trace1,xi,'spline');
trace22=interp1(trace2,xi,'spline');
trace33=interp1(trace3,xi,'spline');
trace44=interp1(trace4,xi,'spline');


plot(trace11, 'k', 'Linewidth', 1); hold on;
plot(trace22, 'r','Linewidth', 0.5); hold on;
plot(trace33, 'b', 'Linewidth', 0.5); hold on;
plot(trace44, 'g', 'Linewidth', 0.5); hold on;


axis([0 10*nz -6e5 6e5]);
xlabel('Depth (m)', 'fontsize',12);
ylabel('Amplitude', 'fontsize',12);
set(gca,'Xtick',[0:500:2500],'Ytick',[-6e5:3e5:6e6],'fontsize',10);
set(gca,'Xticklabel',{[0:500:2000]},'fontsize',10);
set(gca,'YAxisLocation','right');
set(gcf,'Position',[100 100 300 600]); 
set(gca,'Position',[.20 .08 .78 .80]); 

view([90 90]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 6.0]);
print ./Fig/martrace2.eps -depsc2 -r600;


klength=100;
figure;
spec1=abs((fft(trace1)));
spec2=abs((fft(trace2)));
spec3=abs((fft(trace3)));
spec4=abs((fft(trace4)));

specmax=max([max(spec1),max(spec4)]);

plot(spec1/specmax, 'k', 'Linewidth', 1); hold on;
plot(spec2/specmax, 'r','Linewidth', 1); hold on;
plot(spec3/specmax, 'b', 'Linewidth', 0.5); hold on;
plot(spec4/specmax, 'g', 'Linewidth', 0.5); hold on;

axis([0 100 0 1.0]);
xlabel('Wavenumber ({m}^{-1})', 'fontsize',12);
ylabel('Normalized amplitude', 'fontsize',12);
set(gca,'Xtick',[0:20:100],'Ytick',[0:0.2:1.0],'fontsize',10);
set(gca,'Xticklabel',{[0:0.05:0.25]},'fontsize',10);
%grid on;


set(gcf,'Position',[100 100 400 280]); 
set(gca,'Position',[.15 .16 .80 .80]); 

set(gcf,'paperpositionmode','auto');
print ./Fig/marspec2.eps -depsc2 -r600;


%% reftrace=580
figure;

Reftrace=580;

trace1=ima1(:,Reftrace);
trace2=ima2(:,Reftrace);
trace3=ima3(:,Reftrace);
trace4=ima4(:,Reftrace);

xi=1:0.1:length(trace1);

trace11=interp1(trace1,xi,'spline');
trace22=interp1(trace2,xi,'spline');
trace33=interp1(trace3,xi,'spline');
trace44=interp1(trace4,xi,'spline');


plot(trace11, 'k', 'Linewidth', 1); hold on;
plot(trace22, 'r','Linewidth', 0.5); hold on;
plot(trace33, 'b', 'Linewidth', 0.5); hold on;
plot(trace44, 'g', 'Linewidth', 0.5); hold on;


axis([0 10*nz -8e5 8e5]);
xlabel('Depth (m)', 'fontsize',12);
ylabel('Amplitude', 'fontsize',12);
set(gca,'Xtick',[0:500:2500],'Ytick',[-1e6:3e5:1e6],'fontsize',10);
set(gca,'Xticklabel',{[0:500:2000]},'fontsize',10);
set(gca,'YAxisLocation','right');
set(gcf,'Position',[100 100 300 600]); 
set(gca,'Position',[.20 .08 .78 .80]); 

view([90 90]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.0 6.0]);
print ./Fig/martrace3.eps -depsc2 -r600;


klength=100;
figure;
spec1=abs((fft(trace1)));
spec2=abs((fft(trace2)));
spec3=abs((fft(trace3)));
spec4=abs((fft(trace4)));

specmax=max([max(spec1),max(spec4)]);

plot(spec1/specmax, 'k', 'Linewidth', 1); hold on;
plot(spec2/specmax, 'r','Linewidth', 1); hold on;
plot(spec3/specmax, 'b', 'Linewidth', 0.5); hold on;
plot(spec4/specmax, 'g', 'Linewidth', 0.5); hold on;

axis([0 100 0 1.0]);
xlabel('Wavenumber ({m}^{-1})', 'fontsize',12);
ylabel('Normalized amplitude', 'fontsize',12);
set(gca,'Xtick',[0:20:100],'Ytick',[0:0.2:1.0],'fontsize',10);
set(gca,'Xticklabel',{[0:0.05:0.25]},'fontsize',10);
%grid on;


set(gcf,'Position',[100 100 400 280]); 
set(gca,'Position',[.15 .16 .80 .80]); 

set(gcf,'paperpositionmode','auto');
print ./Fig/marspec3.eps -depsc2 -r600;


