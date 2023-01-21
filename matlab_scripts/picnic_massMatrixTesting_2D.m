%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D mass matrix implementation testing
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]
Mi   = 1836.15*me;
mu0 = 4*pi*1e-7;
ep0 = 1/mu0/cvac^2;

sp = 1; ghosts = 0;

% rootPath = '../fromQuartz/2D/numericalEnergyTests/theta_implicit/theta0p5/massMatrix/';
% dataPath = [rootPath,'test0_MM_TEST/']; thisFig = 5;
% %dataPath = [testPath,'test0_MM_TEST_noSigmaxy/']; thisFig = 6;
% dataPath = [rootPath,'testing_CIC/']; thisFig = 5;
% dataPath = [rootPath,'testing_TSC/']; thisFig = 6;
% dataPath = [rootPath,'testing_CC1_2Xings/']; thisFig = 7;
% dataPath = [rootPath,'testing_CC1_1Xing/']; thisFig = 8; 

rootPath = '../fromQuartz/2D/numericalEnergyTests/theta_implicit/theta0p5/massMatrix/';
dataPath = [rootPath,'CIC/']; thisFig = 3;
dataPath = [rootPath,'TSC/']; thisFig = 4;
%dataPath = [rootPath,'CC1_1Xing/']; thisFig = 5;
%dataPath = [rootPath,'CC1_2Xing/']; thisFig = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [dataPath,'mesh_data/mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = squeeze(data.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);

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
%
fileinfo = hdf5info(meshFile);
length_scale = 1;
try
    length_scale = h5readatt(meshFile,'/','length_scale_SI');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   load the particle and field data
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loadData = 1;
if(loadData)

    species_folders = dir([dataPath,'mesh_data/species*']);
    numSpecies = length(species_folders);

    momentList = dir([dataPath,'mesh_data/',species_folders(sp).name,'/moment*']);

    ListLength = length(momentList)

    fieldList = dir([dataPath,'mesh_data/field_data/field*']);
    ListLength_fields = length(fieldList);
    assert(ListLength_fields==ListLength)

    step = zeros(size(fieldList));
    index = zeros(size(fieldList));
    for n=1:ListLength
        thisFile_fields = fieldList(n).name;
        step(n) = str2num(thisFile_fields(7:end-3));
    end
    [step,index] = sort(step);

    iLmax = ListLength;
    time = zeros(1,iLmax);
    totalParts = zeros(1,iLmax);
    %
    JX = zeros(nX,nZ+1,iLmax);
    JY = zeros(nX+1,nZ,iLmax);
    JZ = zeros(nX+1,nZ+1,iLmax);
    %
    JX_test = zeros(nX,nZ+1,iLmax);
    JY_test = zeros(nX+1,nZ,iLmax);
    JZ_test = zeros(nX+1,nZ+1,iLmax);

    for iL=1:iLmax

        %%%   reading moments from part file for species 1
        %
        momentFile_0 = [dataPath,'mesh_data/',species_folders(sp).name, ...
                 '/',momentList(index(iL)).name];
        fileinfo = hdf5info(momentFile_0);

        SpaceDim = h5readatt(momentFile_0,'/Chombo_global','SpaceDim');
        time(iL) = h5readatt(momentFile_0,'/species_data','time');
        if(iL==1)
            time_scale = h5readatt(momentFile_0,'/species_data','time_scale_SI');
            Mass_0 = h5readatt(momentFile_0,'/species_data','mass');
            Charge_0 = double(h5readatt(momentFile_0,'/species_data','charge'));
        end

        fieldFile = [dataPath,'mesh_data/field_data/',fieldList(index(iL)).name];
        fileinfo = hdf5info(fieldFile);
        fileinfo.GroupHierarchy.Groups(2).Attributes.Name;

        if(iL==1)
            Escale = h5readatt(fieldFile,'/field_data','electric_field_scale_SI');
            Bscale = h5readatt(fieldFile,'/field_data','magnetic_field_scale_SI');
        end

        %%%   read current density from field file
        %
        groupName = '/current_density';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        JX(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac; % [Amps/m^2]
        JY(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac; % [Amps/m^2]

        %%%   read virtual current density from field file
        %
        groupName = '/virtual_current_density';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        JZ(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac;  % [Amps/m^2]
        
        %%%   read current density test from field file
        %
        groupName = '/current_density_test';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        JX_test(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac; % [Amps/m^2]
        JY_test(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac; % [Amps/m^2]

        %%%   read virtual current density test from field file
        %
        groupName = '/current_density_virtual_test';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        JZ_test(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac;  % [Amps/m^2]
        
        display(iL);
    end

          
end

diff_JX = JX_test-JX;
diff_JY = JY_test-JY;
diff_JZ = JZ_test-JZ;

it = length(time);

close(figure(thisFig));
f1 = figure(thisFig); set(f1,'position',[200 100 1130 900]);
set(gcf,'color','white');
%
%
%
subplot(3,3,7);
pcolor(Xnc,Znc,JZ(:,:,it)); shading flat; colorbar; 
axis('equal'); axis('tight');
title('J_z'); xlabel('X'); ylabel('Y');
%
subplot(3,3,8);
pcolor(Xnc,Znc,JZ_test(:,:,it)); shading flat; colorbar;
axis('equal'); axis('tight');
title('J_z from MM'); xlabel('X'); ylabel('Y');
%
subplot(3,3,9);
pcolor(Xnc,Znc,diff_JZ(:,:,it)); shading flat; colorbar;
axis('equal'); axis('tight');
title('diff J_z'); xlabel('X'); ylabel('Y');
%
%
%
subplot(3,3,4);
pcolor(Xec1,Zec1,JY(:,:,it)); shading flat; colorbar; 
axis('equal'); axis('tight');
title('J_y'); xlabel('X'); ylabel('Y');
%
subplot(3,3,5);
pcolor(Xec1,Zec1,JY_test(:,:,it)); shading flat; colorbar;
axis('equal'); axis('tight');
title('J_y from MM'); xlabel('X'); ylabel('Y');
%
subplot(3,3,6);
pcolor(Xec1,Zec1,diff_JY(:,:,it)); shading flat; colorbar;
axis('equal'); axis('tight');
title('diff J_y'); xlabel('X'); ylabel('Y');
%
%
%
subplot(3,3,1);
pcolor(Xec0,Zec0,JX(:,:,it)); shading flat; colorbar; 
axis('equal'); axis('tight');
title('J_x'); xlabel('X'); ylabel('Y');
%
subplot(3,3,2);
pcolor(Xec0,Zec0,JX_test(:,:,it)); shading flat; colorbar;
axis('equal'); axis('tight');
title('J_x from MM'); xlabel('X'); ylabel('Y');
%
subplot(3,3,3);
pcolor(Xec0,Zec0,diff_JX(:,:,it)); shading flat; colorbar;
axis('equal'); axis('tight');
title('diff J_x'); xlabel('X'); ylabel('Y');
