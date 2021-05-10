%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D numerical energy conservation tests using PICNIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]
mu0 = 4*pi*1e-7;
ep0 = 1/mu0/cvac^2;

species = 1; thisFig = 1;

rootPath = '../fromQuartz/numericalEnergyTests/explicit/test2p0_fieldsOn_noCollisions/'; thisFig = 1;
rootPath = '../fromQuartz/numericalEnergyTests/explicit/test2p0_fieldsOff_withCollisions/'; thisFig = 2;
rootPath = '../fromQuartz/numericalEnergyTests/explicit/test2p0_fieldsOn_withCollisions/'; thisFig = 3;
%
%rootPath = '../fromQuartz/numericalEnergyTests/semi_implicit/noCollisions/test2p0_iterMax2/'; thisFig = 4;
%rootPath = '../fromQuartz/numericalEnergyTests/semi_implicit/noCollisions/test2p0_iterMax5/'; thisFig = 5;
%rootPath = '../fromQuartz/numericalEnergyTests/semi_implicit/noCollisions/test2p0_iterMax10/'; thisFig = 6;
%
%rootPath = '../fromQuartz/numericalEnergyTests/semi_implicit/withCollisions/test2p0_iterMax2/'; thisFig = 7;
%rootPath = '../fromQuartz/numericalEnergyTests/semi_implicit/withCollisions/test2p0_iterMax5/'; thisFig = 8;
%rootPath = '../fromQuartz/numericalEnergyTests/semi_implicit/withCollisions/test2p0_iterMax10/'; thisFig = 9;
%
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/noCollisions/test2p0_iterMax2/'; thisFig = 10;
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/noCollisions/test2p0_iterMax6/'; thisFig = 11;
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/noCollisions/test2p0_iterMax12/'; thisFig = 12;
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/noCollisions/test2p0_iterMax18/'; thisFig = 13;
%
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/withCollisions/test2p0_iterMax2/'; thisFig = 14;
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/withCollisions/test2p0_iterMax6/'; thisFig = 15;
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/withCollisions/test2p0_iterMax12/'; thisFig = 16;
%rootPath = '../fromQuartz/numericalEnergyTests/fully_implicit/withCollisions/test2p0_iterMax18/'; thisFig = 17;

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


partList = dir([rootPath,'particle_data/species',num2str(species), ...
                               '_data/part*']);
momentList = dir([rootPath,'mesh_data/species',num2str(species), ...
                           '_data/moments*']);
ListLength = length(partList);
assert(ListLength==length(momentList));

fieldList = dir([rootPath,'mesh_data/field_data/field*']);
ListLength_fields = length(fieldList);
assert(ListLength_fields==ListLength)

step = zeros(size(partList));
index = zeros(size(partList));
for n=1:ListLength
    thisFile_parts = partList(n).name;
    thisFile_fields = fieldList(n).name;
    step(n) = str2num(thisFile_parts(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
totalParts = zeros(1,iLmax);
%
numberDen_1 = zeros(nX,nZ,iLmax);
momentumDenX_1 = zeros(nX,nZ,iLmax);
momentumDenY_1 = zeros(nX,nZ,iLmax);
momentumDenZ_1 = zeros(nX,nZ,iLmax);
energyDenX_1 = zeros(nX,nZ,iLmax);
energyDenY_1 = zeros(nX,nZ,iLmax);
energyDenZ_1 = zeros(nX,nZ,iLmax);
JX_1 = zeros(nX,nZ+1,iLmax);
JY_1 = zeros(nX+1,nZ,iLmax);
JZ_1 = zeros(nX+1,nZ+1,iLmax);
%
numberDen_2 = zeros(nX,nZ,iLmax);
momentumDenX_2 = zeros(nX,nZ,iLmax);
momentumDenY_2 = zeros(nX,nZ,iLmax);
momentumDenZ_2 = zeros(nX,nZ,iLmax);
energyDenX_2 = zeros(nX,nZ,iLmax);
energyDenY_2 = zeros(nX,nZ,iLmax);
energyDenZ_2 = zeros(nX,nZ,iLmax);
JX_2 = zeros(nX,nZ+1,iLmax);
JY_2 = zeros(nX+1,nZ,iLmax);
JZ_2 = zeros(nX+1,nZ+1,iLmax);
%
BX = zeros(nX+1,nZ,iLmax);
BY = zeros(nX,nZ+1,iLmax);
BZ = zeros(nX,nZ,iLmax);
EX = zeros(nX,nZ+1,iLmax);
EY = zeros(nX+1,nZ,iLmax);
EZ = zeros(nX+1,nZ+1,iLmax);
JX = zeros(nX,nZ+1,iLmax);
JY = zeros(nX+1,nZ,iLmax);
JZ = zeros(nX+1,nZ+1,iLmax);

for iL=1:iLmax

    %%%   reading moments from part file for species 1
    %
    partsFile_1 = [rootPath,'particle_data/species',num2str(species), ...
                            '_data/',partList(index(iL)).name];
    fileinfo = hdf5info(partsFile_1);
    
    SpaceDim = h5readatt(partsFile_1,'/Chombo_global','SpaceDim');
    time(iL) = h5readatt(partsFile_1,'/species_data','time');
    if(iL==1)
        time_scale = h5readatt(partsFile_1,'/species_data','time_scale_SI');
        Mass_1 = h5readatt(partsFile_1,'/species_data','mass');
        Charge_1 = double(h5readatt(partsFile_1,'/species_data','charge'));
    end

    momentFile_1 = [rootPath,'mesh_data/species',num2str(species), ...
                           '_data/',momentList(index(iL)).name];
    
    groupName = '/species_data'; ghosts = 0;
    data = import2Ddata_singleFile(momentFile_1,groupName,ghosts);
    numberDen_1(:,:,iL) = squeeze(data.Fcc(:,:,1));            % [1/m^3]
    momentumDenX_1(:,:,iL) = squeeze(data.Fcc(:,:,2))*me*cvac; % [m/s/m^3]
    momentumDenY_1(:,:,iL) = squeeze(data.Fcc(:,:,3))*me*cvac; % [m/s/m^3]
    momentumDenZ_1(:,:,iL) = squeeze(data.Fcc(:,:,4))*me*cvac; % [m/s/m^3]
    energyDenX_1(:,:,iL) = squeeze(data.Fcc(:,:,5))*me*cvac^2; % [J/m^3]
    energyDenY_1(:,:,iL) = squeeze(data.Fcc(:,:,6))*me*cvac^2; % [J/m^3]
    energyDenZ_1(:,:,iL) = squeeze(data.Fcc(:,:,7))*me*cvac^2; % [J/m^3]

    try
        groupName = '/current_density';
        data = import2Ddata_singleFile(momentFile_1,groupName,ghosts);
        JX_1(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac;  
        JY_1(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac;     

        groupName = '/virtual_current_density';
        data = import2Ddata_singleFile(momentFile_1,groupName,ghosts);
        JZ_1(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac;  
    end
    
    %%%   reading moments from part file for species 2
    %
    momentFile_2 = [rootPath,'mesh_data/species2_data/',momentList(index(iL)).name];
    if(iL==1)
        Mass_2 = h5readatt(momentFile_2,'/species_data','mass');
        Charge_2 = double(h5readatt(momentFile_2,'/species_data','charge'));
    end

    groupName = '/species_data';
    data = import2Ddata_singleFile(momentFile_2,groupName,ghosts);
    numberDen_2(:,:,iL) = squeeze(data.Fcc(:,:,1));            % [1/m^3]
    momentumDenX_2(:,:,iL) = squeeze(data.Fcc(:,:,2))*me*cvac; % [m/s/m^3]
    momentumDenY_2(:,:,iL) = squeeze(data.Fcc(:,:,3))*me*cvac; % [m/s/m^3]
    momentumDenZ_2(:,:,iL) = squeeze(data.Fcc(:,:,4))*me*cvac; % [m/s/m^3]
    energyDenX_2(:,:,iL) = squeeze(data.Fcc(:,:,5))*me*cvac^2; % [J/m^3]
    energyDenY_2(:,:,iL) = squeeze(data.Fcc(:,:,6))*me*cvac^2; % [J/m^3]
    energyDenZ_2(:,:,iL) = squeeze(data.Fcc(:,:,7))*me*cvac^2; % [J/m^3]
    
    try
        groupName = '/current_density';
        data = import2Ddata_singleFile(momentFile_2,groupName,ghosts);
        JX_2(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac;  
        JY_2(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac;     

        groupName = '/virtual_current_density'; 
        data = import2Ddata_singleFile(momentFile_2,groupName,ghosts);
        JZ_2(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac; 
    end
    
    %
    %
    %
    
    fieldFile = [rootPath,'mesh_data/field_data/',fieldList(index(iL)).name];
    fileinfo = hdf5info(fieldFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;

    if(iL==1)
        Escale = h5readatt(fieldFile,'/field_data','electric_field_scale_SI');
        Bscale = h5readatt(fieldFile,'/field_data','magnetic_field_scale_SI');
    end
    
    %%%   read magnetic field from field file
    %
    groupName = '/magnetic_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    BX(:,:,iL) = squeeze(data.Ffc0(:,:))*Bscale; % [Tesla]
    BY(:,:,iL) = squeeze(data.Ffc1(:,:))*Bscale; % [Tesla]
    
    %%%   read virtual magnetic field from field file
    %
    groupName = '/virtual_magnetic_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    BZ(:,:,iL) = squeeze(data.Fcc(:,:))*Bscale; % [Tesla]
    
    %%%   read electric field from field file
    %
    groupName = '/electric_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    EX(:,:,iL) = squeeze(data.Fec0(:,:))*Escale; % [V/m]
    EY(:,:,iL) = squeeze(data.Fec1(:,:))*Escale; % [V/m]
    
    %%%   read virtual electric field from field file
    %
    groupName = '/virtual_electric_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    EZ(:,:,iL) = squeeze(data.Fnc(:,:))*Escale; % [V/m]

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
    
    
    display(iL);
end

%%%   compute the total energy in the fields
%
Energy_Efield = zeros(size(time));
Energy_Bfield = zeros(size(time));
dV2D_SI = dX*dZ*length_scale^2;
for iL=1:length(time)
    Energy_Efield(iL) = (sum(sum(EX(:,1:end-1,iL).^2)) ...
                      +  sum(sum(EY(1:end-1,:,iL).^2)) ...
                      +  sum(sum(EZ(1:end-1,1:end-1,iL).^2)))*ep0/2*dV2D_SI; % [Joules]
    Energy_Bfield(iL) = (sum(sum(BX(1:end-1,:,iL).^2)) ...
                      +  sum(sum(BY(:,1:end-1,iL).^2)) ...
                      +  sum(sum(BZ(:,:,iL).^2)))/mu0/2*dV2D_SI; % [Joules]
end

%%%   compute the temperature
%
velX_1 = momentumDenX_1./numberDen_1/(me*Mass_1);
velY_1 = momentumDenY_1./numberDen_1/(me*Mass_1);
velZ_1 = momentumDenZ_1./numberDen_1/(me*Mass_1);
tempX_1 = 2*(energyDenX_1 - momentumDenX_1.^2./numberDen_1/(me*Mass_1)/2.0)./numberDen_1/qe; % [eV]
tempY_1 = 2*(energyDenY_1 - momentumDenY_1.^2./numberDen_1/(me*Mass_1)/2.0)./numberDen_1/qe; % [eV]
tempZ_1 = 2*(energyDenZ_1 - momentumDenZ_1.^2./numberDen_1/(me*Mass_1)/2.0)./numberDen_1/qe; % [eV]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   compute global integrals to check conservations
%%%

numberTot_1 = zeros(1,iLmax);
momentumTotX_1 = zeros(1,iLmax);
momentumTotY_1 = zeros(1,iLmax);
momentumTotZ_1 = zeros(1,iLmax);
energyTotX_1 = zeros(1,iLmax);
energyTotY_1 = zeros(1,iLmax);
energyTotZ_1 = zeros(1,iLmax);
%
numberTot_2 = zeros(1,iLmax);
momentumTotX_2 = zeros(1,iLmax);
momentumTotY_2 = zeros(1,iLmax);
momentumTotZ_2 = zeros(1,iLmax);
energyTotX_2 = zeros(1,iLmax);
energyTotY_2 = zeros(1,iLmax);
energyTotZ_2 = zeros(1,iLmax);
%
dV = dX*dZ*length_scale^2;
for n=1:iLmax
    numberTot_1(n) = sum(sum(numberDen_1(:,:,n)))*dV2D_SI;
    momentumTotX_1(n) = sum(sum(momentumDenX_1(:,:,n)))*dV2D_SI;
    momentumTotY_1(n) = sum(sum(momentumDenY_1(:,:,n)))*dV2D_SI;    
    momentumTotZ_1(n) = sum(sum(momentumDenZ_1(:,:,n)))*dV2D_SI;
    energyTotX_1(n) = sum(sum(energyDenX_1(:,:,n)))*dV2D_SI;
    energyTotY_1(n) = sum(sum(energyDenY_1(:,:,n)))*dV2D_SI;
    energyTotZ_1(n) = sum(sum(energyDenZ_1(:,:,n)))*dV2D_SI;
    %
    numberTot_2(n) = sum(sum(numberDen_2(:,:,n)))*dV2D_SI;
    momentumTotX_2(n) = sum(sum(momentumDenX_2(:,:,n)))*dV2D_SI;
    momentumTotY_2(n) = sum(sum(momentumDenY_2(:,:,n)))*dV2D_SI;    
    momentumTotZ_2(n) = sum(sum(momentumDenZ_2(:,:,n)))*dV2D_SI;
    energyTotX_2(n) = sum(sum(energyDenX_2(:,:,n)))*dV2D_SI;
    energyTotY_2(n) = sum(sum(energyDenY_2(:,:,n)))*dV2D_SI;
    energyTotZ_2(n) = sum(sum(energyDenZ_2(:,:,n)))*dV2D_SI;
end
Energy_species1 = energyTotX_1 + energyTotY_1 + energyTotZ_1; % [Joules]
Energy_species2 = energyTotX_2 + energyTotY_2 + energyTotZ_2; % [Joules]


tempTotX = energyTotX_1 - momentumTotX_1.^2./numberTot_1/(me*Mass_1)/2.0;
tempTotY = energyTotY_1 - momentumTotY_1.^2./numberTot_1/(me*Mass_1)/2.0;
tempTotZ = energyTotZ_1 - momentumTotZ_1.^2./numberTot_1/(me*Mass_1)/2.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Energy_parts = Energy_species1 + Energy_species2;
Energy_fields = Energy_Efield + Energy_Bfield;
Energy_total = Energy_parts + Energy_fields;

close(figure(thisFig)); f1=figure(thisFig);
plot(time,(Energy_total-Energy_total(1))/Energy_total(1),'displayName', ...
           '\Delta\epsilon_t_o_t/\epsilon_0'); hold on;
plot(time,(Energy_species1-Energy_species1(1))/Energy_total(1),'displayName', ...
           '\Delta\epsilon_e/\epsilon_0','color',[0.47 0.67 0.19]); hold on;
plot(time,(Energy_species2-Energy_species2(1))/Energy_total(1),'displayName', ...
           '\Delta\epsilon_i/\epsilon_0','color',[0.85 0.33 0.10]); hold on;
plot(time,(Energy_fields-Energy_fields(1))/Energy_total(1),'displayName', ...
           '\Delta\epsilon_E_M/\epsilon_0','color',[0.93 0.59 0.13]); hold off;
xlabel('\omega_p_et'); title('relative energy change');
legend('show','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(figure(88)); 
f8=figure(88); set(f8,'position',[316 424 1450 405]);

subplot(1,3,1);
plot(time,numberTot_1/numberTot_1(1),'displayName','Mass/Mass(t=0)');
xlabel('time [s]'); title('total mass'); ylim([0.999 1.001]);
legend('show','location','northwest');

subplot(1,3,2);
plot(time,me*momentumTotX_1,'displayName','x-momentum'); hold on;
plot(time,me*momentumTotY_1,'displayName','y-momentum'); hold on;
plot(time,me*momentumTotZ_1,'displayName','z-momentum'); hold off;
xlabel('time [s]'); title('total momentum');
legend('show','location','best');

subplot(1,3,3);
energyTot = me*(energyTotX_1+energyTotY_1+energyTotZ_1);
plot(time,me*energyTotX_1,'displayName','x-energy'); hold on;
plot(time,me*energyTotY_1,'displayName','y-energy'); hold on;
plot(time,me*energyTotZ_1,'displayName','z-energy'); hold off;
xlabel('time [s]'); title('total energy [Joules]');
legend('show','location','best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot global average velocity and temperatures
%%%
%%%

velXavg = squeeze(mean(velX_1,1)); velXavg = squeeze(mean(velXavg,1));
velYavg = squeeze(mean(velY_1,1)); velYavg = squeeze(mean(velYavg,1));
velZavg = squeeze(mean(velZ_1,1)); velZavg = squeeze(mean(velZavg,1));
%
tempXavg = squeeze(mean(tempX_1,1)); tempXavg = squeeze(mean(tempXavg,1));
tempYavg = squeeze(mean(tempY_1,1)); tempYavg = squeeze(mean(tempYavg,1));
tempZavg = squeeze(mean(tempZ_1,1)); tempZavg = squeeze(mean(tempZavg,1));


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it = 4;
% 
% close(figure(12)); f12 = figure(12);
% set(f12,'position',[420 312 1310 700]);
% 
% subplot(2,3,1); pcolor(JX_1(:,:,it)+JX_2(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,2); pcolor(JY_1(:,:,it)+JY_2(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,3); pcolor(JZ_1(:,:,it)+JZ_2(:,:,it)); shading flat; colorbar; axis('square');
% %
% subplot(2,3,4); pcolor(JX(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,5); pcolor(JY(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,6); pcolor(JZ(:,:,it)); shading flat; colorbar; axis('square');
% 
% 

