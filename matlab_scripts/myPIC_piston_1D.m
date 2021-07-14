%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D piston module using PICNIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

species = 1;
rootPath = '../fromQuartz/1D/pistonTests/piston_collisionless/'; vpiston = 100;
rootPath = '../fromQuartz/1D/pistonTests/piston_collisionless_fixedSeed/'; vpiston = 100;
rootPath = '../fromQuartz/1D/pistonTests/piston_vp1e4/'; vpiston = 1e4;

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
momentumDenX = zeros(nX,iLmax);
momentumDenY = zeros(nX,iLmax);
momentumDenZ = zeros(nX,iLmax);
energyDenX = zeros(nX,iLmax);
energyDenY = zeros(nX,iLmax);
energyDenZ = zeros(nX,iLmax);


%%%  loop over files and create movie
%
close(figure(1));
f1=figure(1); set(f1,'position',[888 570 900 450]);
set(gcf,'color','white');
%
images = cell(1,1);
v=VideoWriter('./figs/pistonCar.mp4', 'MPEG-4');
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
    energyDenX(:,iL) = squeeze(data.Fcc(:,5))*cvac^2;
    energyDenY(:,iL) = squeeze(data.Fcc(:,6))*cvac^2;
    energyDenZ(:,iL) = squeeze(data.Fcc(:,7))*cvac^2;
    %
    subplot(1,2,1); %hold on;
    p1=plot(particle.x,particle.vx/vpiston,'*'); box on;
    title('particle x-vx phase space'); axis([0 1 -10 10]); %axis([0 1 -10 10]);
    set(gca,'xtick',0:0.2:1); %set(gca,'ytick',-10:4:10);
    xpiston = 1-vpiston*time(iL);
    hold on; l1=line([xpiston xpiston],[-10 10],'linestyle','--','color','black');
    hold off;
    xlabel('x/R'); ylabel('vx/vp'); axis('square');
    %
    if(iL==1) 
        n0 = mean(numberDen(:,iL));
    end
    subplot(1,2,2); %hold on; 
    p2=plot(Xcc,numberDen(:,iL)/n0); box on;
    hold on; l2=line([xpiston xpiston],[0 10],'linestyle','--','color','black');
    hold off;
    title('density'); axis([0 1 0 10]);
    set(gca,'xtick',0:0.2:1); set(gca,'ytick',0:2:10);
    xlabel('x/R'); ylabel('n/n_0'); axis('square');


    NameString = './figs/pistonCar';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p1); delete(l1);
        delete(p2); delete(l2);
    end

end
close(v);


%%%   compute the temperature
%
velX = momentumDenX./numberDen/Mass;
velY = momentumDenY./numberDen/Mass;
velZ = momentumDenZ./numberDen/Mass;
tempX = 2*(energyDenX*me - me*Mass*numberDen.*velX.^2/2.0)./numberDen/qe;
tempY = 2*(energyDenY*me - me*Mass*numberDen.*velY.^2/2.0)./numberDen/qe;
tempZ = 2*(energyDenZ*me - me*Mass*numberDen.*velZ.^2/2.0)./numberDen/qe;

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


close(figure(7)); f7=figure(7); 
plot(particle.x,particle.vz);
xlabel('x'); ylabel('vz velocity');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%         compute analytic shock jump conditions
%%%         and compare with the simulation
%%%

g = 5/3;
Tavg = (tempX + tempY + tempZ)/3;
Cs = sqrt(g*qe*Tavg/(me*Mass)); % sound speed [m/s]
Cs0 = mean(Cs(:,1));
VpovCs0 = abs(vpiston)/Cs0; 

%%%   express shock jump conditions on 1D pressure grid
%%%   and find the values corresponding to vp/Cs0
%
P2ovP1 = 1:0.01:1e3; % pressure jump
N2ovN1 = ((g+1)*P2ovP1 + (g-1))./((g-1)*P2ovP1 + (g+1)); % density jump
USovC1 = sqrt(1+(g+1)/2/g*(P2ovP1-1)); % shock speed / upstream sound speed
U2ovC1 = (1-1./N2ovN1).*USovC1;

%%%   find index corresonding to Vpiston/Cs0 and compute the
%%%   expected shock-jump ratios for pressure and density
%
[~,i0] = min(abs(U2ovC1-VpovCs0));
P2ovP10 = P2ovP1(i0);
N2ovN10 = N2ovN1(i0);
U2ovC10 = U2ovC1(i0);
display(VpovCs0);
display(P2ovP10);
display(N2ovN10);

Navg = squeeze(mean(numberDen,2));
Pavg = Tavg.*numberDen;
P0 = mean(Tavg(:,1))*n0;

close(figure(10));
f10=figure(10); set(f10,'position',[800 30 1000 400]);
for m=2:2:14
    subplot(1,2,1); hold on;
    plot(Xcc,numberDen(:,m)/n0);
    %
    subplot(1,2,2); hold on;
    plot(Xcc,Pavg(:,m)/P0);    
end
subplot(1,2,1);
line([0 1],[N2ovN10 N2ovN10],'color','black','linestyle','--'); hold off; 
box on; grid on;
xlabel('x'); ylabel('\rho/\rho_0');
title('density profiles');
%
subplot(1,2,2);
line([0 1],[P2ovP10 P2ovP10],'color','black','linestyle','--'); hold off; 
box on; grid on;
xlabel('x'); ylabel('P/P_0');
title('pressure profiles');

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

