%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D thermalization test using myPIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%addpath('~angus1/Programs/COGENT_matlabTools/');
%addpath('../matlabScripts/');

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

rootPath = '../myPIC/thermalization_test0/';

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxeslinewidth',1.5);
set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesfontweight','bold');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh.h5'];
%fileinfo = hdf5info(meshFile);
group = 2; ghosts = 0;
data = import2Ddata_singleFile(meshFile,group,ghosts);
Xcc = squeeze(data.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);
Xce = squeeze(data.block(1).Fce(:,:,1));
Zce = squeeze(data.block(1).Fce(:,:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles and density
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%  loop over files and create movie
%

close(figure(1));
f1=figure(1); set(f1,'position',[570 450 900 420]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('./pistonCar.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);


fileList = dir([rootPath,'particle_data/part*']);
ListLength = length(fileList)

step = zeros(size(fileList));
index = zeros(size(fileList));
for n=1:ListLength
    n
    thisFile = fileList(n).name
    step(n) = str2num(thisFile(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
totalParts = zeros(1,iLmax);
numberDen = zeros(nX,nZ,iLmax);
momentumDenX = zeros(nX,nZ,iLmax);
momentumDenY = zeros(nX,nZ,iLmax);
momentumDenZ = zeros(nX,nZ,iLmax);
energyDenX = zeros(nX,nZ,iLmax);
energyDenY = zeros(nX,nZ,iLmax);
energyDenZ = zeros(nX,nZ,iLmax);

for iL=1:iLmax


    partsFile = [rootPath,'particle_data/',fileList(index(iL)).name];
    %fileinfo = hdf5info(partsFile);
    partData = hdf5read(partsFile,'/level_0/particles:data');
    SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
    numParts = h5readatt(partsFile,'/level_0','num_particles');
    time(iL) = h5readatt(partsFile,'/level_0','time');
    if(iL==1)
        Mass = h5readatt(partsFile,'/level_0','mass');
    end
    partData = reshape(partData,2+3*SpaceDim+3,numParts);
    partData = partData';
    totalParts(iL) = numParts;
    if(SpaceDim==2)
       particle.x    = partData(:,1);
       particle.y    = partData(:,2);
       particle.z    = partData(:,3); % virtual dimesion
       particle.vx   = partData(:,4);
       particle.vy   = partData(:,5);
       particle.vz   = partData(:,6); % virtual dimension
       particle.ax   = partData(:,7);
       particle.ay   = partData(:,8);
       particle.az   = partData(:,9); % virtual dimension
       particle.weight = partData(:,10);
       particle.ID   = partData(:,11);
    end

    %%%   reading density from part file
    %
    group = 2; ghosts = 0;
    data = import2Ddata_singleFile(partsFile,group,ghosts);
    numberDen(:,:,iL) = squeeze(data.Fcc(:,:,1));
    momentumDenX(:,:,iL) = squeeze(data.Fcc(:,:,2));
    momentumDenY(:,:,iL) = squeeze(data.Fcc(:,:,3));
    momentumDenZ(:,:,iL) = squeeze(data.Fcc(:,:,4));
    energyDenX(:,:,iL) = squeeze(data.Fcc(:,:,5));
    energyDenY(:,:,iL) = squeeze(data.Fcc(:,:,6));
    energyDenZ(:,:,iL) = squeeze(data.Fcc(:,:,7));
    %
    rhoAvg = sum(numberDen(:,:,iL),2)/length(Zcc(1,:));
    momAvgX = sum(momentumDenX(:,:,iL),2)/length(Zcc(1,:));
    momAvgY = sum(momentumDenY(:,:,iL),2)/length(Zcc(1,:));
    momAvgZ = sum(momentumDenZ(:,:,iL),2)/length(Zcc(1,:));
    eneAvgX = sum(energyDenX(:,:,iL),2)/length(Zcc(1,:));
    eneAvgY = sum(energyDenY(:,:,iL),2)/length(Zcc(1,:));
    eneAvgZ = sum(energyDenZ(:,:,iL),2)/length(Zcc(1,:));
    %
    subplot(1,2,1); %hold on;
    p1=plot(particle.x,particle.vx,'*'); box on;
    title('particle x-vx phase space'); %axis([0 1 -10 10]);
    %set(gca,'xtick',0:0.2:1); set(gca,'ytick',-10:4:10);
    xlabel('x/R'); ylabel('vx/vp'); axis('square');
    %
    if(iL==1) 
        rho0 = sum(rhoAvg)/length(rhoAvg);
    end
    subplot(1,2,2); %hold on; 
    p2=plot(Xcc(:,1),rhoAvg/rho0); box on;
    title('plasma density'); axis([0 1 0 2]);
    %set(gca,'xtick',0:0.2:1); set(gca,'ytick',0:2:10);
    xlabel('x/R'); ylabel('\rho/\rho_0'); axis('square');


    NameString = 'pistonCar';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p1);
        delete(p2);
    end

end
close(v);


%%%   compute the temperature
%
velX = momentumDenX./numberDen/Mass;
velY = momentumDenY./numberDen/Mass;
velZ = momentumDenZ./numberDen/Mass;
tempX = 2*(energyDenX - momentumDenX.^2./numberDen/Mass/2.0)./numberDen*me/qe;
tempY = 2*(energyDenY - momentumDenY.^2./numberDen/Mass/2.0)./numberDen*me/qe;
tempZ = 2*(energyDenZ - momentumDenZ.^2./numberDen/Mass/2.0)./numberDen*me/qe;

for i=1:nX
    for j=1:nZ
        for n=1:iLmax
            thisDen = numberDen(i,j,n);
            if(thisDen==0)
                velX(i,j,n) = 0;
                velY(i,j,n) = 0;
                velZ(i,j,n) = 0;
                tempX(i,j,n) = 0;
                tempY(i,j,n) = 0;
                tempZ(i,j,n) = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%
%%%
%%%   bin up the particles in velocity space and look at distribution
%%%
%%%%%%%%%%%%%%%%%%

Vmax = max(particle.vy);
Vmin = min(particle.vy);
Vgrid = linspace(Vmin,Vmax,100);
dVy = Vgrid(2)-Vgrid(1);

dfn = zeros(size(Vgrid));
for i=1:numParts
    thisv = particle.vy(i);
    [~,index]=min(abs(thisv-Vgrid));
    dfn(index) = dfn(index) + 1;
end
normC = sum(dfn)*dVy;
dfn = dfn/normC;

close(figure(44)); f4=figure(44); 

T1 = 25.9e-3;
VT1 = 4.19e5*sqrt(T1/Mass);
plot(Vgrid,exp(-(Vgrid/sqrt(2)/VT1).^2)/sqrt(2*pi)/VT1); hold on;
plot(Vgrid,dfn); 
xlabel('y-velocity'); ylabel('dfn-vy');
title('vy-distribution function');
legend('maxwellian','data'); axis('tight');


close(figure(7)); f7=figure(7); 
plot(particle.x,particle.vz);
xlabel('x'); ylabel('vz velocity');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   compute global integrals to check conservations
%%%

numberTot = zeros(1,iLmax);
momentumTotX = zeros(1,iLmax);
momentumTotY = zeros(1,iLmax);
momentumTotZ = zeros(1,iLmax);
energyTotX = zeros(1,iLmax);
energyTotY = zeros(1,iLmax);
energyTotZ = zeros(1,iLmax);
for n=1:iLmax
    numberTot(n) = sum(sum(numberDen(:,:,n)))*dX*dZ;
    momentumTotX(n) = sum(sum(momentumDenX(:,:,n)))*dX*dZ;
    momentumTotY(n) = sum(sum(momentumDenY(:,:,n)))*dX*dZ;    
    momentumTotZ(n) = sum(sum(momentumDenZ(:,:,n)))*dX*dZ;
    energyTotX(n) = sum(sum(energyDenX(:,:,n)))*dX*dZ;
    energyTotY(n) = sum(sum(energyDenY(:,:,n)))*dX*dZ;
    energyTotZ(n) = sum(sum(energyDenZ(:,:,n)))*dX*dZ;
end
tempTotX = energyTotX - momentumTotX.^2./numberTot/Mass/2.0;
tempTotY = energyTotY - momentumTotY.^2./numberTot/Mass/2.0;
tempTotZ = energyTotZ - momentumTotZ.^2./numberTot/Mass/2.0;

close(figure(8)); 
f8=figure(8); set(f8,'position',[316 424 1450 405]);

subplot(1,3,1);
plot(time,numberTot/numberTot(1),'displayName','Mass/Mass(t=0)');
xlabel('time [s]'); title('total mass'); ylim([0.999 1.001]);
legend('show','location','northwest');

subplot(1,3,2);
plot(time,me*momentumTotX,'displayName','x-momentum'); hold on;
plot(time,me*momentumTotY,'displayName','y-momentum'); hold on;
plot(time,me*momentumTotZ,'displayName','z-momentum'); hold off;
xlabel('time [s]'); title('total momentum');
legend('show','location','best');

subplot(1,3,3);
% plot(time,tempTotX,'displayName','Temperature-X'); hold on;
% plot(time,tempTotY,'displayName','Temperature-Y'); hold on;
% plot(time,tempTotZ,'displayName','Temperature-Z'); hold off;
plot(time,me*energyTotX,'displayName','x-energy'); hold on;
plot(time,me*energyTotY,'displayName','y-energy'); hold on;
plot(time,me*energyTotZ,'displayName','z-energy'); hold off;
xlabel('time [s]'); title('total energy [Joules]');
legend('show','location','best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot global average velocity and temperatures
%%%
%%%

velXavg = squeeze(mean(velX,1)); velXavg = squeeze(mean(velXavg,1));
velYavg = squeeze(mean(velY,1)); velYavg = squeeze(mean(velYavg,1));
velZavg = squeeze(mean(velZ,1)); velZavg = squeeze(mean(velZavg,1));
%
tempXavg = squeeze(mean(tempX,1)); tempXavg = squeeze(mean(tempXavg,1));
tempYavg = squeeze(mean(tempY,1)); tempYavg = squeeze(mean(tempYavg,1));
tempZavg = squeeze(mean(tempZ,1)); tempZavg = squeeze(mean(tempZavg,1));


close(figure(9)); 
f9=figure(9); set(f9,'position',[316 424 1050 405]);

subplot(1,2,1);
plot(time,velXavg,'displayName','x-velocity'); hold on;
plot(time,velYavg,'displayName','y-velocity'); hold on;
plot(time,velZavg,'displayName','z-velocity'); hold off;
xlabel('time [s]'); title('global average velocity');
ylabel('velocity [m/s]'); axis('tight');
legend('show','location','best');

subplot(1,2,2);
plot(time,tempXavg,'displayName','x-temperature'); hold on;
plot(time,tempYavg,'displayName','y-temperature'); hold on;
plot(time,tempZavg,'displayName','z-temperature'); hold off;
xlabel('time [s]'); title('global average temperature');
ylabel('temperature [eV]');
legend('show','location','best'); axis('tight');
tempAvg0 = (tempXavg(1)+tempYavg(1)+tempZavg(1))/3;
hold on; l1=line([time(1) time(end)],[tempAvg0 tempAvg0],'linestyle','--');
set(l1,'color','black');


