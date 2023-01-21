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
rootPath = '../fromQuartz/1D/numericalEnergyTests/theta_implicit/theta0p5/massMatrix/';
dataPath = [rootPath,'CIC/']; thisFig = 3;
dataPath = [rootPath,'TSC/']; thisFig = 4;
dataPath = [rootPath,'CC1_1Xing/']; thisFig = 5;
dataPath = [rootPath,'CC1_2Xing/']; thisFig = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [dataPath,'mesh_data/mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid';
data = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = data.Fcc(:,1); nX = length(Xcc(:,1)); dX = Xcc(2)-Xcc(1);

groupName = '/face_centered_grid';
data2 = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xfc0 = data2.Ffc0(:,1); 

groupName = '/edge_centered_grid';
data2 = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xec0 = data2.Fec0(:,1);

groupName = '/node_centered_grid';
data2 = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xnc = data2.Fnc(:,1); 
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
    BX = zeros(nX+1,iLmax);  curlE_X = zeros(size(BX));
    BY = zeros(nX,iLmax);    curlE_Y = zeros(size(BY));
    BZ = zeros(nX,iLmax);    curlE_Z = zeros(size(BZ));
    EX = zeros(nX,iLmax);    curlB_X = zeros(size(EX));
    EY = zeros(nX+1,iLmax);  curlB_Y = zeros(size(EY));
    EZ = zeros(nX+1,iLmax);  curlB_Z = zeros(size(EZ));
    JX = zeros(nX,iLmax);
    JY = zeros(nX+1,iLmax);
    JZ = zeros(nX+1,iLmax);
    %
    JX_test = zeros(nX,iLmax);
    JY_test = zeros(nX+1,iLmax);
    JZ_test = zeros(nX+1,iLmax);

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

        %%%   reading moments from part file for species 2
        %
        momentFile_1 = [dataPath,'mesh_data/',species_folders(sp+1).name, ...
                 '/',momentList(index(iL)).name];  
        if(iL==1)
            Mass_1 = h5readatt(momentFile_1,'/species_data','mass');
            Charge_1 = double(h5readatt(momentFile_1,'/species_data','charge'));
        end

        fieldFile = [dataPath,'mesh_data/field_data/',fieldList(index(iL)).name];
        fileinfo = hdf5info(fieldFile);
        fileinfo.GroupHierarchy.Groups(2).Attributes.Name;

        if(iL==1)
            Escale = h5readatt(fieldFile,'/field_data','electric_field_scale_SI');
            Bscale = h5readatt(fieldFile,'/field_data','magnetic_field_scale_SI');
        end

        %%%   read magnetic field from field file
        %
        groupName = '/magnetic_field';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        BX(:,iL) = data.Ffc0*Bscale; % [Tesla]

        %%%   read virtual magnetic field from field file
        %
        groupName = '/virtual_magnetic_field';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        BY(:,iL) = squeeze(data.Fcc(:,1))*Bscale; % [Tesla]
        BZ(:,iL) = squeeze(data.Fcc(:,2))*Bscale; % [Tesla]

        %%%   read electric field from field file
        %
        groupName = '/electric_field';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        EX(:,iL) = data.Fec0*Escale; % [V/m]

        %%%   read virtual electric field from field file
        %
        groupName = '/virtual_electric_field';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        EY(:,iL) = squeeze(data.Fnc(:,1))*Escale; % [V/m]
        EZ(:,iL) = squeeze(data.Fnc(:,2))*Escale; % [V/m]

        %%%   read current density from field file
        %
        groupName = '/current_density';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        JX(:,iL) = data.Fec0*qe*cvac; % [Amps/m^2]

        %%%   read virtual current density from field file
        %
        groupName = '/virtual_current_density';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        JY(:,iL) = squeeze(data.Fnc(:,1))*qe*cvac; % [Amps/m^2]
        JZ(:,iL) = squeeze(data.Fnc(:,2))*qe*cvac;  % [Amps/m^2]
        
        %%%   read current density test from field file
        %
        groupName = '/current_density_test';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        JX_test(:,iL) = data.Fec0*qe*cvac; % [Amps/m^2]

        %%%   read virtual current density test from field file
        %
        groupName = '/current_density_virtual_test';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        JY_test(:,iL) = squeeze(data.Fnc(:,1))*qe*cvac; % [Amps/m^2]
        JZ_test(:,iL) = squeeze(data.Fnc(:,2))*qe*cvac;  % [Amps/m^2]
        
        display(iL);
    end

          
end

diff_JX = JX_test-JX;
diff_JY = JY_test-JY;
diff_JZ = JZ_test-JZ;

it = length(time);

close(figure(thisFig));
f1 = figure(thisFig); set(f1,'position',[950 25 800 980]);
%
subplot(3,2,1);
plot(Xec0,JX(:,it)); box on; grid on; hold on;
plot(Xec0,JX_test(:,it),'r--'); title('J_x');
%
subplot(3,2,3);
plot(Xnc,JY(:,it)); box on; grid on; hold on;
plot(Xnc,JY_test(:,it),'r--'); title('J_y');
%
subplot(3,2,5);
plot(Xnc,JZ(:,it)); box on; grid on; hold on;
plot(Xnc,JZ_test(:,it),'r--'); title('J_z');
%
subplot(3,2,2);
plot(Xec0,diff_JX(:,it)); box on; grid on;
title('diff J_x');
%
subplot(3,2,4);
plot(Xnc,diff_JY(:,it)); box on; grid on;
title('diff J_y');
%
subplot(3,2,6);
plot(Xnc,diff_JZ(:,it)); box on; grid on;
title('diff J_z');



