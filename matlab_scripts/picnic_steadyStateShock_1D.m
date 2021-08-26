%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D steady state shock using picnic
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

species = 1; gamma = 5/3;
rootPath = '../fromQuartz/1D/steadyStateShock/neutrals/test0_noBCs/';
rootPath = '../fromQuartz/1D/steadyStateShock/neutrals/test0_withBCs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh_data/mesh.h5'];
% vecData  = hdf5read(thisFile,[groupName,'/data:datatype=0']);
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = data.Fcc; nX = length(Xcc); dX = Xcc(2)-Xcc(1);

groupName = '/face_centered_grid'; ghosts = 0;
data = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xfc = data.Ffc0;
X0 = Xfc(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles and density
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partList   = dir([rootPath,'particle_data/species',num2str(species), ...
                           '_data/parts*']);
momentList = dir([rootPath,'mesh_data/species',num2str(species), ...
                           '_data/moments*']);
ListLength = length(partList);
assert(ListLength==length(momentList));
% 
step = zeros(size(partList));
index = zeros(size(partList));
for n=1:ListLength
    thisFile = partList(n).name;
    step(n) = str2num(thisFile(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
totalParts = zeros(1,iLmax);
numberDen = zeros(nX,iLmax);
momentumDenX = zeros(nX,iLmax); velX = zeros(nX,iLmax);
momentumDenY = zeros(nX,iLmax); velY = zeros(nX,iLmax);
momentumDenZ = zeros(nX,iLmax); velZ = zeros(nX,iLmax);
energyDenX = zeros(nX,iLmax); tempX = zeros(nX,iLmax);
energyDenY = zeros(nX,iLmax); tempY = zeros(nX,iLmax);
energyDenZ = zeros(nX,iLmax); tempZ = zeros(nX,iLmax);
%
tempAvg = zeros(nX,iLmax);
pressAvg = zeros(nX,iLmax);
eintAvg = zeros(nX,iLmax);

%%%  loop over files and create movie
%
close(figure(11));
f1=figure(11); set(f1,'position',[860 240 900 760]);
set(gcf,'color','white');
%
images = cell(1,1);
v=VideoWriter('./figs/sodShock.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

%iLmax = 1;
for iL=1:iLmax

    partsFile = [rootPath,'particle_data/species',num2str(species), ...
                          '_data/',partList(index(iL)).name];
    fileinfo = hdf5info(partsFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;
    
    partData = hdf5read(partsFile,'/species_data/particles:data');
    SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
    numParts = h5readatt(partsFile,'/species_data','num_particles');
    time(iL) = h5readatt(partsFile,'/species_data','time');
    if(iL==1)
        Mass = h5readatt(partsFile,'/species_data','mass');
        Charge = double(h5readatt(partsFile,'/species_data','charge'));
        Uint = h5readatt(partsFile,'/species_data','Uint'); % [eV]
        numPartComps = h5readatt(partsFile,'/species_data','numPartComps');
    end
    partData = reshape(partData,numPartComps,numParts);
    partData = partData';
    totalParts(iL) = numParts;

    particle.weight = partData(:,1);
    particle.x    = partData(:,2);
    particle.y    = partData(:,3);
    particle.z    = partData(:,4); 
    particle.vx   = partData(:,5)*cvac;
    particle.vy   = partData(:,6)*cvac;
    particle.vz   = partData(:,7)*cvac;
    particle.ID   = partData(:,numPartComps);


    momentFile = [rootPath,'mesh_data/species',num2str(species), ...
                           '_data/',momentList(index(iL)).name];
    
    %%%   reading density from species moment file
    %
    groupName = '/species_data'; ghosts = 0;
    data = import1Ddata_singleFile(momentFile,groupName,ghosts);
    numberDen(:,iL) = data.Fcc(:,1);
    momentumDenX(:,iL) = data.Fcc(:,2)*cvac;
    momentumDenY(:,iL) = data.Fcc(:,3)*cvac;
    momentumDenZ(:,iL) = data.Fcc(:,4)*cvac;
    energyDenX(:,iL) = squeeze(data.Fcc(:,5))*cvac^2; tempX(:,iL) = 2*(energyDenX(:,iL)*me - me*Mass*numberDen(:,iL).*velX(:,iL).^2/2.0)./numberDen(:,iL)/qe;
    energyDenY(:,iL) = squeeze(data.Fcc(:,6))*cvac^2;
    energyDenZ(:,iL) = squeeze(data.Fcc(:,7))*cvac^2;
    
    %%%   compute physical variables
    
    velX(:,iL) = momentumDenX(:,iL)./numberDen(:,iL)/Mass; % [m/s]
    velY(:,iL) = momentumDenY(:,iL)./numberDen(:,iL)/Mass; % [m/s]
    velZ(:,iL) = momentumDenZ(:,iL)./numberDen(:,iL)/Mass; % [m/s]
    
    tempX(:,iL) = 2*me/qe*(energyDenX(:,iL) - Mass*numberDen(:,iL).*velX(:,iL).^2/2.0)./numberDen(:,iL);
    tempY(:,iL) = 2*me/qe*(energyDenY(:,iL) - Mass*numberDen(:,iL).*velY(:,iL).^2/2.0)./numberDen(:,iL);
    tempZ(:,iL) = 2*me/qe*(energyDenZ(:,iL) - Mass*numberDen(:,iL).*velZ(:,iL).^2/2.0)./numberDen(:,iL);
   
    tempAvg(:,iL) = (tempX(:,iL) + tempY(:,iL) + tempZ(:,iL))/3;
    pressAvg(:,iL) = tempAvg(:,iL).*numberDen(:,iL);
    eintAvg(:,iL) = tempAvg(:,iL)/(gamma-1);
    
    if(iL==1) 
        n0 = min(numberDen(:,iL));
        T0 = mean(tempAvg(1:nX/2,iL));
        P0 = n0*T0;
        V0 = sqrt(qe*T0/(me*Mass)*gamma);
    end
    
    subplot(2,2,1); %hold on; 
    p1=plot(Xcc/X0,numberDen(:,iL)/n0); box on;
    title('density'); grid on; 
    set(gca,'xtick',0:0.25:1); %ylim([0 1.2]);
    xlabel('x/x_0'); ylabel('n/n_0'); axis('square');
    %
    subplot(2,2,2); %hold on;
    p2=plot(Xcc/X0,velX(:,iL)/V0); box on;
    title('x-velocity'); grid on; 
    set(gca,'xtick',0:0.25:1); %ylim([-0.2 1.2]);
    xlabel('x/x_0'); ylabel('v/v_0'); axis('square');
    %
    subplot(2,2,3); %hold on;
    p3=plot(Xcc/X0,pressAvg(:,iL)/P0); box on;
    title('pressure'); grid on;
    set(gca,'xtick',0:0.25:1); %ylim([0 1.2]);
    xlabel('x/x_0'); ylabel('p/p_0'); axis('square');
    %
    subplot(2,2,4); %hold on;
    p4=plot(Xcc/X0,eintAvg(:,iL)/T0); box on;
    title('internal energy'); grid on;
    set(gca,'xtick',0:0.25:1); %ylim([0.8 2.2]);
    xlabel('x/x_0'); ylabel('T/(\gamma-1)/T_0'); axis('square');

    NameString = './figs/sodShock1D';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    
    frame = getframe(gcf);
    %frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p1); delete(p2);
        delete(p3); delete(p4);
    end

end
close(v);


for i=1:nX
    for n=1:iLmax
        thisDen = numberDen(i,n);
        if(thisDen==0)
            velX(i,n) = 0;
            velY(i,n) = 0;
            velZ(i,n) = 0;
            tempX(i,n) = 0;
            tempY(i,n) = 0;
            tempZ(i,n) = 0;
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

TY0 = mean(nonzeros(tempY(:,iLmax)));
UY0 = mean(nonzeros(velY(:,iLmax)));
VTY = 4.19e5*sqrt(TY0/Mass);
plot(Vgrid,exp(-((Vgrid-UY0)/sqrt(2)/VTY).^2)/sqrt(2*pi)/VTY); hold on;
plot(Vgrid,dfn); 
xlabel('y-velocity'); ylabel('dfn-vy');
title('vy-distribution function');
legend('maxwellian','data'); axis('tight');




