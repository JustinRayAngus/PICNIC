%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D thermalization test using myPIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

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
ListLength = length(fileList);

step = zeros(size(fileList));
for n=1:ListLength
    thisFile = fileList(n).name;
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


