%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D test of CIC method for particle deposition
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

species = 1;
rootPath = '../myPIC/depositTests/test0_2D/';
%rootPath = '../myPIC/depositTests/test1_2D/';
ghosts = 1; % must have ghosts in output for this test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh.h5'];

groupName = '/cell_centered_grid';
data1 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = squeeze(data1.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data1.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);

groupName = '/face_centered_grid';
data2 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xfc0 = squeeze(data2.Ffc0(:,:,1)); 
Zfc0 = squeeze(data2.Ffc0(:,:,2)); 
Xfc1 = squeeze(data2.Ffc1(:,:,1)); 
Zfc1 = squeeze(data2.Ffc1(:,:,2)); 

groupName = '/edge_centered_grid';
data2 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xec0 = squeeze(data2.Fec0(:,:,1)); 
Zec0 = squeeze(data2.Fec0(:,:,2)); 
Xec1 = squeeze(data2.Fec1(:,:,1)); 
Zec1 = squeeze(data2.Fec1(:,:,2)); 

groupName = '/node_centered_grid';
data2 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xnc = squeeze(data2.Fnc(:,:,1)); 
Znc = squeeze(data2.Fnc(:,:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles and moments
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%  loop over files and create movie
%
close(figure(1));
f1=figure(1); set(f1,'position',[740 540 1000 500]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('./deposit2D.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fileList = dir([rootPath,'species',num2str(species),'_data/part*']);
ListLength = length(fileList);

step = zeros(size(fileList));
index = zeros(size(fileList));
for n=1:ListLength
    thisFile = fileList(n).name;
    step(n) = str2num(thisFile(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
totalWeight = zeros(2,iLmax);
totalParts = zeros(1,iLmax);
numberDen = zeros(nX,nZ,iLmax);
chargeDen = zeros(nX,nZ,iLmax);
Jnc = zeros(nX+1,nZ+1,iLmax);

%%%   compute the weights manually to verify
%
W00 = zeros(1,iLmax);
W01 = zeros(1,iLmax);
W10 = zeros(1,iLmax);
W11 = zeros(1,iLmax);
error = zeros(4,iLmax);

%iLmax = 1;
for iL=1:iLmax

    partsFile = [rootPath,'species',num2str(species),'_data/',fileList(index(iL)).name];
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
    if(SpaceDim==2)
       particle.weight = partData(:,1);
       particle.x    = partData(:,2);
       particle.y    = partData(:,3);
       particle.z    = partData(:,4); 
       particle.vx   = partData(:,5);
       particle.vy   = partData(:,6);
       particle.vz   = partData(:,7);
     %  particle.ax   = partData(:,8);
     %  particle.ay   = partData(:,9);
     %  particle.az   = partData(:,10);
       particle.ID   = partData(:,numPartComps);
    end

    %%%   reading current density from part file
    %
    groupName = '/species_data';
    data3 = import2Ddata_singleFile(partsFile,groupName,ghosts);
    numberDen(:,:,iL) = squeeze(data3.Fcc(:,:,1));
    
    groupName = '/cell_centered_charge_density';
    data3p1 = import2Ddata_singleFile(partsFile,groupName,ghosts);
    chargeDen(:,:,iL) = squeeze(data3p1.Fcc(:,:,end));
    if(iL==1)
        rho0 = sum(sum(chargeDen(2:end-1,2:end-1,1)));
    end
    
    %%% make contour plot movie of deposited particle weight
    %
    i0 = 1+data3p1.numGhostX; i1 = nX-data3p1.numGhostX;
    j0 = 1+data3p1.numGhostY; j1 = nZ-data3p1.numGhostY;    
    chargeDenC = squeeze(chargeDen(i0:i1,j0:j1,iL));
    chargeDenC(end+1,:) = chargeDenC(end,:);
    chargeDenC(:,end+1) = chargeDenC(:,end);  
    
    subplot(1,2,1);
    p1=pcolor(Zfc1(1,j0:j1+1), ...
              Xfc0(i0:i1+1,1),chargeDenC/rho0); colorbar; box on; 
    hold on; p2=plot(particle.y,particle.x,'r*');
    title('CIC deposit on cells'); axis('equal'); axis('tight');
    xlabel('Z'); ylabel('X');
    %
    %
    %
    groupName = '/virtual_current_density';
    data3p2 = import2Ddata_singleFile(partsFile,groupName,ghosts);
    Jnc(:,:,iL) = squeeze(data3p2.Fnc(:,:,end));
    if(iL==1)
        Jnc0 = sum(sum(Jnc(2:end-2,2:end-2,1)));
    end
    
    subplot(1,2,2);
    p3=pcolor(Zcc(1,1:end),Xcc(1:end,1),Jnc(2:end,2:end,iL)/Jnc0); 
    colorbar; box on; 
    hold on; p4=plot(particle.y,particle.x,'r*');
    hold on; l1=line([2 3],[1 1],'linestyle','--','color','r'); % Xhi face
    hold on; l2=line([2 3],[0 0],'linestyle','--','color','r'); % Xlo face
    hold on; l3=line([3 3],[0 1],'linestyle','--','color','r'); % Zhi face
    hold on; l4=line([2 2],[0 1],'linestyle','--','color','r'); % Zlo face
    title('CIC deposit on nodes'); axis('equal'); axis('tight');
    xlabel('Z'); ylabel('X');   
    
    
    NameString = 'deposit2D';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p1); delete(p2);
        delete(p3); delete(p4);
        delete(l1); delete(l2); delete(l3); delete(l4);
    end

     
    %%%   compute the weights manually using CIC to verify
    %
    
    totalWeightCells(iL) = sum(sum(chargeDen(i0:i1,j0:j1,iL)))/rho0;
    totalWeightNodes(iL) = sum(sum(Jnc(2:end-2,2:end-2,iL)))/Jnc0;
    %
    if(length(particle.x)==1)
    [~,i1] = min(abs(particle.x-Xfc0(:,1)));
    [~,j1] = min(abs(particle.y-Zfc1(1,:)));
    Wx1 = 1-abs(particle.x-Xcc(i1,1))/dX;
    Wz1 = 1-abs(particle.y-Zcc(1,j1))/dZ;
        
    i0 = i1-1;
    j0 = j1-1;
     
    Wx0 = 1-abs(particle.x-Xcc(i0,1))/dX;
    Wz0 = 1-abs(particle.y-Zcc(1,j0))/dZ;

    W00(iL) = Wx0*Wz0;
    W01(iL) = Wx0*Wz1;
    W10(iL) = Wx1*Wz0;
    W11(iL) = Wx1*Wz1;
    %Sum = W00 + W01 + W10 + W11; % should be one

    error(1,iL) = abs(chargeDen(i0,j0,iL)/rho0-W00(iL));
    error(2,iL) = abs(chargeDen(i0,j1,iL)/rho0-W01(iL));
    error(3,iL) = abs(chargeDen(i1,j0,iL)/rho0-W10(iL));
    error(4,iL) = abs(chargeDen(i1,j1,iL)/rho0-W11(iL));
    end
    
end
close(v);

close(figure(3));
f3=figure(3);
plot(time(1:iLmax),totalWeightCells(1:iLmax),'displayName','cells'); hold on; 
plot(time(1:iLmax),totalWeightNodes(1:iLmax),'displayName',...
                                       'nodes','linestyle','--'); hold off;
xlabel('time'); title('total Weight change');
legend('show');

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


display(['max error = ',num2str(max(max(error)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%             plot the deposit on cell faces
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupName = '/face_centered_charge_density'; patchData=1;
data4 = import2Ddata_singleFile(partsFile,groupName,ghosts,patchData);
Jfc0 = data4.Ffc0;
Jfc1 = data4.Ffc1;

close(figure(4));
f4=figure(4); set(f4,'position',[1160 50 560 800]);
%
subplot(2,1,1);
pcolor(Zfc1(1,1:end-1),Xfc1(1:end,1),Jfc0(2:end,:)/rho0); colorbar;
hold on; plot(particle.y,particle.x,'r*');
hold on; line([2 3],[1 1],'linestyle','--','color','r'); % Xhi face
hold on; line([2 3],[0 0],'linestyle','--','color','r'); % Xlo face
xlabel('Z'); ylabel('X');
title('CIC deposit on X-faces');
axis([2 3 0-dX/2 1+dX/2]);
%
subplot(2,1,2);
pcolor(Zfc0(1,1:end),Xfc0(1:end-1,1),Jfc1(:,2:end)/rho0); colorbar;
hold on; plot(particle.y,particle.x,'r*');
hold on; line([3 3],[0 1],'linestyle','--','color','r'); % Zhi face
hold on; line([2 2],[0 1],'linestyle','--','color','r'); % Zlo face
xlabel('Z'); ylabel('X'); 
title('CIC deposit on Z-faces');
axis([2-dZ/2 3+dZ/2 0 1]);


%%%   compute particle weights for X-faces grid
%
if(length(particle.x)==1)
Wz0_stagX = Wz0;
Wz1_stagX = Wz1;
[~,i0] = min(abs(particle.x-Xcc(:,1)));
i1 = i0+1;
Wx0_stagX = 1-abs(particle.x-Xfc0(i0,1))/dX;
Wx1_stagX = 1-abs(particle.x-Xfc0(i1,1))/dX;


Wx0_stagZ = Wx0;
Wx1_stagZ = Wx1;
[~,j0] = min(abs(particle.y-Zcc(1,:)));
j1 = j0+1;
Wz0_stagZ = 1-abs(particle.y-Zfc1(1,j0))/dZ;
Wz1_stagZ = 1-abs(particle.y-Zfc1(1,j1))/dZ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%             plot the deposit on cell faces
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupName = '/current_density'; patchData=0;
data5 = import2Ddata_singleFile(partsFile,groupName,ghosts,patchData);
Jec0 = data5.Fec0;
Jec1 = data5.Fec1;

close(figure(5));
f5=figure(5); set(f5,'position',[680 50 560 800]);
%
subplot(2,1,1);
pcolor(Zec0(1,1:end-1),Xec0(1:end,1),abs(Jec1(2:end,:))/rho0); colorbar;
hold on; plot(particle.y,particle.x,'r*');
hold on; line([2 3],[1 1],'linestyle','--','color','r'); % Xhi face
hold on; line([2 3],[0 0],'linestyle','--','color','r'); % Xlo face
xlabel('Z'); ylabel('X');
title('CIC deposit on X-edge');
axis([2 3 0-dX/2 1+dX/2]);
%
subplot(2,1,2);
pcolor(Zec1(1,1:end),Xec1(1:end-1,1),Jec0(:,2:end)/rho0); colorbar;
hold on; plot(particle.y,particle.x,'r*');
hold on; line([3 3],[0 1],'linestyle','--','color','r'); % Zhi face
hold on; line([2 2],[0 1],'linestyle','--','color','r'); % Zlo face
xlabel('Z'); ylabel('X'); 
title('CIC deposit on Z-edges');
axis([2-dZ/2 3+dZ/2 0 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


