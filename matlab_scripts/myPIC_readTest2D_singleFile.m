%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D reading test for myPIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

rootPath = '../myPIC/thermalization_test0/';

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

partsFile = [rootPath,'particle_data/parts0000.h5'];
%fileinfo = hdf5info(partsFile);
partData = hdf5read(partsFile,'/level_0/particles:data');
SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
numParts = h5readatt(partsFile,'/level_0','num_particles');
time = h5readatt(partsFile,'/level_0','time');
Mass = h5readatt(partsFile,'/level_0','mass');

partData = reshape(partData,2+3*SpaceDim+3,numParts);
partData = partData';
if(SpaceDim==2)
   particle.x    = partData(:,1);
   particle.y    = partData(:,2);
   particle.z    = partData(:,3);
   particle.vx   = partData(:,4);
   particle.vy   = partData(:,5);
   particle.vz   = partData(:,6);
   particle.ax   = partData(:,7);
   particle.ay   = partData(:,8);
   particle.az   = partData(:,9);
   particle.weight = partData(:,10);
   particle.ID   = partData(:,11);
end

%%%   reading moments from part file
%
group = 2; ghosts = 0;
data = import2Ddata_singleFile(partsFile,group,ghosts);
numberDen = squeeze(data.Fcc(:,:,1));
momentumDenX = squeeze(data.Fcc(:,:,2));
momentumDenY = squeeze(data.Fcc(:,:,3));
momentumDenZ = squeeze(data.Fcc(:,:,4));
energyDenX = squeeze(data.Fcc(:,:,5));
energyDenY = squeeze(data.Fcc(:,:,6));
energyDenZ = squeeze(data.Fcc(:,:,7));
%
rhoAvg = mean(numberDen,2);
momAvgX = mean(momentumDenX,2);
momAvgY = mean(momentumDenY,2);
momAvgZ = mean(momentumDenZ,2);
eneAvgX = mean(energyDenX,2);
eneAvgY = mean(energyDenY,2);
eneAvgZ = mean(energyDenZ,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    plot some results
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(figure(1));
f1=figure(1); set(f1,'position',[570 450 900 420]);
set(gcf,'color','white');


subplot(1,2,1); %hold on;
p1=plot(particle.x,particle.vx,'*'); box on;
title('particle x-vx phase space');
xlabel('x/R'); ylabel('vx/vp'); axis('square');
%
subplot(1,2,2);
rho0 = sum(rhoAvg)/length(rhoAvg);
p2=plot(Xcc(:,1),rhoAvg/rho0); box on;
title('plasma density'); axis([0 1 0 2]);
xlabel('x/R'); ylabel('\rho/\rho_0'); axis('square');




