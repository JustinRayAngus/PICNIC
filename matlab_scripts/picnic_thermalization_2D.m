%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D thermalization test using PICNIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]

species = 0;
rootPath = '../fromQuartz/thermalizationTests/test0/';
rootPath = '../fromQuartz/thermalizationTests/test1/';
rootPath = '../fromQuartz/thermalizationTests/test2_oldNorm/';
rootPath = '../fromQuartz/thermalizationTests/test2_newNorm/';
rootPath = '../fromQuartz/thermalizationTests/test2_newAdvance/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh_data/mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = squeeze(data.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);
%
fileinfo = hdf5info(meshFile);
length_scale = 1;
try
    length_scale = h5readatt(meshFile,'/','length_scale_SI');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles and density
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partList = dir([rootPath,'particle_data/species',num2str(species), ...
                         '_data/part*']);
momentList = dir([rootPath,'mesh_data/species',num2str(species), ...
                           '_data/moments*']);
ListLength = length(partList);
assert(ListLength==length(momentList));

step = zeros(size(partList));
index = zeros(size(partList));
for n=1:ListLength
    n
    thisFile = partList(n).name
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


close(figure(11));
f1=figure(11); set(f1,'position',[570 450 900 450]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('./figs/thermalization.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);


for iL=1:iLmax

    partsFile = [rootPath,'particle_data/species',num2str(species), ...
                          '_data/',partList(index(iL)).name];
    fileinfo = hdf5info(partsFile);
    
    SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
    partData = hdf5read(partsFile,'/species_data/particles:data');
    numParts = h5readatt(partsFile,'/species_data','num_particles');
    time(iL) = h5readatt(partsFile,'/species_data','time');
    if(iL==1)
        time_scale = 1;
        try
            time_scale = h5readatt(partsFile,'/species_data','time_scale_SI');
        end
        Mass = h5readatt(partsFile,'/species_data','mass');
        Charge = double(h5readatt(partsFile,'/species_data','charge'));
        Uint = h5readatt(partsFile,'/species_data','Uint'); % [eV]
        numPartComps = h5readatt(partsFile,'/species_data','numPartComps');
    end
    partData = reshape(partData,numPartComps,numParts);
    partData = partData';
    totalParts(iL) = numParts;
    if(SpaceDim==2)
       particle.weight = partData(:,1);
       particle.x    = partData(:,2);
       particle.y    = partData(:,3);
       particle.z    = partData(:,4); 
       particle.vx   = partData(:,5)*cvac;
       particle.vy   = partData(:,6)*cvac;
       particle.vz   = partData(:,7)*cvac;
     %  particle.ax   = partData(:,8);
     %  particle.ay   = partData(:,9);
     %  particle.az   = partData(:,10);
       particle.ID   = partData(:,numPartComps);
    end

    momentFile = [rootPath,'mesh_data/species',num2str(species), ...
                           '_data/',momentList(index(iL)).name];
    
    %%%   reading moments from part file
    %
    groupName = '/species_data'; ghosts = 0;
    data = import2Ddata_singleFile(momentFile,groupName,ghosts);
    numberDen(:,:,iL) = squeeze(data.Fcc(:,:,1));
    momentumDenX(:,:,iL) = squeeze(data.Fcc(:,:,2))*cvac;
    momentumDenY(:,:,iL) = squeeze(data.Fcc(:,:,3))*cvac;
    momentumDenZ(:,:,iL) = squeeze(data.Fcc(:,:,4))*cvac;
    energyDenX(:,:,iL) = squeeze(data.Fcc(:,:,5))*cvac^2;
    energyDenY(:,:,iL) = squeeze(data.Fcc(:,:,6))*cvac^2;
    energyDenZ(:,:,iL) = squeeze(data.Fcc(:,:,7))*cvac^2;
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
        display(rho0);
    end
    subplot(1,2,2); %hold on; 
    p2=plot(Xcc(:,1),rhoAvg/rho0); box on;
    title('plasma density'); axis('tight');
    ylim([0 2]);
    %set(gca,'xtick',0:0.2:1); set(gca,'ytick',0:2:10);
    xlabel('x/R'); ylabel('\rho/\rho_0'); axis('square');


    NameString = './figs/therm2D';
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

close(figure(4)); f4=figure(4); 

TY0 = mean(mean(nonzeros(tempY(:,:,end))));
UY0 = mean(mean(nonzeros(velY(:,:,end))));
VTY = 4.19e5*sqrt(TY0/Mass);
plot(Vgrid,exp(-((Vgrid-UY0)/sqrt(2)/VTY).^2)/sqrt(2*pi)/VTY); hold on;
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
dV = dX*dZ*length_scale^2;
for n=1:iLmax
    numberTot(n) = sum(sum(numberDen(:,:,n)))*dV;
    momentumTotX(n) = sum(sum(momentumDenX(:,:,n)))*dV;
    momentumTotY(n) = sum(sum(momentumDenY(:,:,n)))*dV;    
    momentumTotZ(n) = sum(sum(momentumDenZ(:,:,n)))*dV;
    energyTotX(n) = sum(sum(energyDenX(:,:,n)))*dV;
    energyTotY(n) = sum(sum(energyDenY(:,:,n)))*dV;
    energyTotZ(n) = sum(sum(energyDenZ(:,:,n)))*dV;
end
tempTotX = energyTotX - momentumTotX.^2./numberTot/Mass/2.0;
tempTotY = energyTotY - momentumTotY.^2./numberTot/Mass/2.0;
tempTotZ = energyTotZ - momentumTotZ.^2./numberTot/Mass/2.0;

close(figure(88)); 
f8=figure(88); set(f8,'position',[316 424 1450 405]);

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
energyTot = me*(energyTotX+energyTotY+energyTotZ);
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


close(figure(99)); 
f9=figure(99); set(f9,'position',[316 424 1050 405]);

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


