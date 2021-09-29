%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D steady state shock in plasma using picnic
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;  % electron mass [kg]
mp  = 1.67262192e-27;     % proton mass [kg]
qe   = 1.602176634e-19;   % electron charge [C]
ep0  = 8.85418782e-12;
cvac = 2.99792458e8;      % speed of light [m/s]


basePath = '../fromQuartz/1D/steadyStateShock/plasma/hydrogen1/mach6/uniformFlow/'; Z = 1; Mach = 6;
%basePath = '../fromQuartz/1D/steadyStateShock/plasma/hydrogen1/mach24/uniformFlow/'; Z = 1; Mach = 24;
%rootPath = [basePath,'L200/inflow_outflow/']; thisFig = 8;
% rootPath = [basePath,'L200/inflow_outflow_old/']; thisFig = 9;
% rootPath = [basePath,'L200/locally_periodic/']; thisFig = 10;
%rootPath = [basePath,'L200/periodic/']; thisFig = 11;
%rootPath = [basePath,'L200/inflow_outflow_noForces/']; thisFig = 14;
%rootPath = [basePath,'L200/inflow_outflow_smallDt/']; thisFig = 15;
% rootPath = [basePath,'L200/inflow_outflow_implicit/']; thisFig = 11;
% rootPath = [basePath,'L200/inflow_outflow_implicit2/']; thisFig = 12;

rootPath = [basePath,'L200/testing/']; thisFig = 13;
rootPath = [basePath,'L200/testing2/']; thisFig = 14; Mach = 12;
rootPath = [basePath,'L200/testing2p1/']; thisFig = 15; Mach = 12;
rootPath = [basePath,'L200/testing2p2/']; thisFig = 16; Mach = 12;
rootPath = [basePath,'L200/testing2p2_refined/']; thisFig = 17; Mach = 12;
rootPath = [basePath,'L200/testing2p2_halfDt/']; thisFig = 18; Mach = 12;
rootPath = [basePath,'L200/testing2p2_twiceDt/']; thisFig = 19; Mach = 12;
rootPath = [basePath,'L200/testing2p2_moreParts/']; thisFig = 20; Mach = 12;
rootPath = [basePath,'L200/testing2p2_periodic/']; thisFig = 21; Mach = 12;
rootPath = [basePath,'L200/testing2p2_symmetryBClo/']; thisFig = 22; Mach = 12;
rootPath = [basePath,'L200/testing2p2_wtf/']; thisFig = 23; Mach = 12;


% basePath = '../fromQuartz/1D/steadyStateShock/plasma/hydrogen1/mach6/'; Z = 1; Mach = 6;
%basePath = '../fromQuartz/1D/steadyStateShock/plasma/hydrogen1/mach24/'; Z = 1; Mach = 24;
% 
%rootPath = [basePath,'L200/Np200_dt0p20_testing/']; thisFig = 1;
% rootPath = [basePath,'L200/Np400_dt0p20/']; thisFig = 2;
% %rootPath = [basePath,'L200/Np200_dt1p00_implicit/']; thisFig = 3;


% basePath = '../fromQuartz/1D/steadyStateShock/plasma/helium4/mach6/'; Z = 2;
% % rootPath = [basePath,'L200/Np200_dt0p20/']; thisFig = 6;
% % rootPath = [basePath,'L400/Np200_dt0p20/']; thisFig = 7;
% rootPath = [basePath,'L800/Np200_dt0p20/']; thisFig = 8;
% %rootPath = [basePath,'L800/Np200_dt0p20_noInterCollisions/']; thisFig = 9;
% % rootPath = [basePath,'L800/Np200_dt0p50_implicit/']; thisFig = 10;
% % rootPath = [basePath,'L800/Np200_dt1p00_implicit/']; thisFig = 11;
% rootPath = [basePath,'L800/Np400_dt1p00_implicit/']; thisFig = 12;


% basePath = '../fromQuartz/1D/steadyStateShock/plasma/beryllium9/mach6/';  Z = 4;
% rootPath = [basePath,'L800/Np100_dt0p20/']; thisFig = 1;
% rootPath = [basePath,'L1600/Np100_dt0p20/']; thisFig = 3;

basePath = '../fromQuartz/1D/steadyStateShock/plasma/hydrogen1/mach6/';
rootPath = [basePath,'fullyEM/test0_periodic/']; thisFig = 4; Mach = 6;
rootPath = [basePath,'fullyEM/test0_inflow/']; thisFig = 5;
rootPath = [basePath,'fullyEM/test0_symmetry/']; thisFig = 6;
rootPath = [basePath,'fullyEM/test0_implicit/']; thisFig = 7;
%rootPath = [basePath,'fullyEM/test0/']; thisFig = 7;


%%% see planaShockJump.m
%
gamma = 5/3; gp1 = gamma+1; gm1 = gamma-1;
N2oN1 = gp1*Mach^2/(gm1*Mach^2+2);
P2oP1 = (2*gamma*Mach^2 - gm1)/gp1;
T2oT1 = P2oP1/N2oN1;
pistonMach = (1-1./N2oN1)*Mach;


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

species_folders = dir([rootPath,'mesh_data/species*']);
numSpecies = length(species_folders);

%partList_1 = dir([rootPath,'particle_data/',species_folders(1).name,'/part*']);
%partList_2 = dir([rootPath,'particle_data/',species_folders(1).name,'/part*']);
momentList = dir([rootPath,'mesh_data/',species_folders(1).name,'/moment*']);
momentList_2 = dir([rootPath,'mesh_data/',species_folders(2).name,'/moment*']);

ListLength = length(momentList);
assert(ListLength==length(momentList_2));

fieldList = dir([rootPath,'mesh_data/field_data/field*']);
ListLength_fields = length(fieldList);
assert(ListLength_fields==ListLength)

step = zeros(size(momentList));
index = zeros(size(momentList));
for n=1:ListLength
    thisFile = fieldList(n).name;
    step(n) = str2num(thisFile(7:end-3));
end
[step,index] = sort(step);

iLmax = min(ListLength,301)
dt_out = 2; iLvect = 1:dt_out:iLmax;
iLmax = length(iLvect);

time = zeros(1,iLmax);
%
EX = zeros(nX+2,iLmax); % +2 is for ghosts
JX = zeros(nX+2,iLmax); % +2 is for ghosts
EY = zeros(nX+1,iLmax);
EZ = zeros(nX+1,iLmax);
%
%totalParts_1 = zeros(1,iLmax);
numberDen_1 = zeros(nX,iLmax);
momDenX_1 = zeros(nX,iLmax); velX_1 = zeros(nX,iLmax);
momDenY_1 = zeros(nX,iLmax); velY_1 = zeros(nX,iLmax);
momDenZ_1 = zeros(nX,iLmax); velZ_1 = zeros(nX,iLmax);
eneDenX_1 = zeros(nX,iLmax); tempX_1 = zeros(nX,iLmax); EmeanX_1 = zeros(nX,iLmax);
eneDenY_1 = zeros(nX,iLmax); tempY_1 = zeros(nX,iLmax); EmeanY_1 = zeros(nX,iLmax);
eneDenZ_1 = zeros(nX,iLmax); tempZ_1 = zeros(nX,iLmax); EmeanZ_1 = zeros(nX,iLmax);
%
tempAvg_1 = zeros(nX,iLmax); presAvg_1 = zeros(nX,iLmax); 
eintAvg_1 = zeros(nX,iLmax);
%
%totalParts_2 = zeros(1,iLmax);
numberDen_2 = zeros(nX,iLmax);
momDenX_2 = zeros(nX,iLmax); velX_2 = zeros(nX,iLmax);
momDenY_2 = zeros(nX,iLmax); velY_2 = zeros(nX,iLmax);
momDenZ_2 = zeros(nX,iLmax); velZ_2 = zeros(nX,iLmax);
eneDenX_2 = zeros(nX,iLmax); tempX_2 = zeros(nX,iLmax); EmeanX_2 = zeros(nX,iLmax);
eneDenY_2 = zeros(nX,iLmax); tempY_2 = zeros(nX,iLmax); EmeanY_2 = zeros(nX,iLmax);
eneDenZ_2 = zeros(nX,iLmax); tempZ_2 = zeros(nX,iLmax); EmeanZ_2 = zeros(nX,iLmax);
%
tempAvg_2 = zeros(nX,iLmax); presAvg_2 = zeros(nX,iLmax); 
eintAvg_2 = zeros(nX,iLmax);
%

%%%  loop over files and create movie
%
close(figure(thisFig));
f1=figure(thisFig); set(f1,'position',[860 240 900 760]);
set(gcf,'color','white');
%
images = cell(1,1);
v=VideoWriter('./figs/plasmaShock.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

for iL = 1:length(iLvect)
    it = iLvect(iL);
%for iL = [1,iLmax]   

    %
    %   load field data
    %
    fieldFile = [rootPath,'mesh_data/field_data/',fieldList(index(it)).name];
    
    if(iL==1)
        time_scale = h5readatt(fieldFile,'/field_data','time_scale_SI');
        length_scale = h5readatt(meshFile,'/','length_scale_SI');
        Escale = h5readatt(fieldFile,'/field_data','electric_field_scale_SI');
        Bscale = h5readatt(fieldFile,'/field_data','magnetic_field_scale_SI');
    end

    %%%   read magnetic field from field file
    %
    groupName = '/electric_field';
    data = import1Ddata_singleFile(fieldFile,groupName,1);
    EX(:,iL) = data.Fec0*Escale; % [V/m]

    %%%   read virtual magnetic field from field file
    %
    groupName = '/virtual_electric_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    EY(:,iL) = squeeze(data.Fnc(:,1))*Escale; % [V/m]
    EZ(:,iL) = squeeze(data.Fnc(:,2))*Escale; % [V/m]

    try
        groupName = '/div_electric_field';
        data = import1Ddata_singleFile(fieldFile,groupName,0);
        divE_code(:,iL) = ep0*data.Fnc*Escale/length_scale;  % SI
    end
    
    groupName = '/current_density';
    data = import1Ddata_singleFile(fieldFile,groupName,1);
    JX(:,iL) = data.Fec0*qe*cvac; % [A/m^2]
    
    %
    %   load data for species1
    %
    momentFile_1 = [rootPath,'mesh_data/',species_folders(1).name, ...
                 '/',momentList(index(iL)).name]; 
    if(iL==1)
        time_scale = h5readatt(momentFile_1,'/species_data','time_scale_SI');
        Mass_1 = h5readatt(momentFile_1,'/species_data','mass');
        Charge_1 = double(h5readatt(momentFile_1,'/species_data','charge'));
    end
    time(iL) = h5readatt(momentFile_1,'/species_data','time');
    
    %%%   reading density from species moment file
    %
    groupName = '/species_data'; ghosts = 0;
    data = import1Ddata_singleFile(momentFile_1,groupName,ghosts);
    numberDen_1(:,iL) = data.Fcc(:,1);
    momDenX_1(:,iL) = data.Fcc(:,2)*cvac;
    momDenY_1(:,iL) = data.Fcc(:,3)*cvac;
    momDenZ_1(:,iL) = data.Fcc(:,4)*cvac;
    eneDenX_1(:,iL) = squeeze(data.Fcc(:,5))*me*cvac^2/qe; % [eV/m^3]
    eneDenY_1(:,iL) = squeeze(data.Fcc(:,6))*me*cvac^2/qe; % [eV/m^3]
    eneDenZ_1(:,iL) = squeeze(data.Fcc(:,7))*me*cvac^2/qe; % [eV/m^3]
    
    %%%   compute physical variables
    
    velX_1(:,iL) = momDenX_1(:,iL)./numberDen_1(:,iL)/Mass_1; % [m/s]
    velY_1(:,iL) = momDenY_1(:,iL)./numberDen_1(:,iL)/Mass_1; % [m/s]
    velZ_1(:,iL) = momDenZ_1(:,iL)./numberDen_1(:,iL)/Mass_1; % [m/s]
    
    EmeanX_1(:,iL) = me*Mass_1*numberDen_1(:,iL).*velX_1(:,iL).^2/2.0/qe; % [eV/m^3]
    EmeanY_1(:,iL) = me*Mass_1*numberDen_1(:,iL).*velY_1(:,iL).^2/2.0/qe; % [eV/m^3]
    EmeanZ_1(:,iL) = me*Mass_1*numberDen_1(:,iL).*velZ_1(:,iL).^2/2.0/qe; % [eV/m^3]
    
    tempX_1(:,iL) = 2*(eneDenX_1(:,iL) - EmeanX_1(:,iL))./numberDen_1(:,iL);
    tempY_1(:,iL) = 2*(eneDenY_1(:,iL) - EmeanY_1(:,iL))./numberDen_1(:,iL);
    tempZ_1(:,iL) = 2*(eneDenZ_1(:,iL) - EmeanZ_1(:,iL))./numberDen_1(:,iL);
   
    tempAvg_1(:,iL) = (tempX_1(:,iL) + tempY_1(:,iL) + tempZ_1(:,iL))/3;
    presAvg_1(:,iL) = tempAvg_1(:,iL).*numberDen_1(:,iL);
    eintAvg_1(:,iL) = tempAvg_1(:,iL)/(gamma-1);
    
    %
    %   load data for species2
    %

    momentFile_2 = [rootPath,'mesh_data/',species_folders(2).name, ...
                 '/',momentList(index(iL)).name]; 
    if(iL==1)
        time_scale = h5readatt(momentFile_2,'/species_data','time_scale_SI');
        Mass_2 = h5readatt(momentFile_2,'/species_data','mass');
        Charge_2 = double(h5readatt(momentFile_2,'/species_data','charge'));
    end
    
    %%%   reading density from species moment file
    %
    groupName = '/species_data'; ghosts = 0;
    data = import1Ddata_singleFile(momentFile_2,groupName,ghosts);
    numberDen_2(:,iL) = data.Fcc(:,1);
    momDenX_2(:,iL) = data.Fcc(:,2)*cvac;
    momDenY_2(:,iL) = data.Fcc(:,3)*cvac;
    momDenZ_2(:,iL) = data.Fcc(:,4)*cvac;
    eneDenX_2(:,iL) = squeeze(data.Fcc(:,5))*me*cvac^2/qe; % [eV/m^3]
    eneDenY_2(:,iL) = squeeze(data.Fcc(:,6))*me*cvac^2/qe; % [eV/m^3]
    eneDenZ_2(:,iL) = squeeze(data.Fcc(:,7))*me*cvac^2/qe; % [eV/m^3]
    
    %%%   compute physical variables
    
    velX_2(:,iL) = momDenX_2(:,iL)./numberDen_2(:,iL)/Mass_2; % [m/s]
    velY_2(:,iL) = momDenY_2(:,iL)./numberDen_2(:,iL)/Mass_2; % [m/s]
    velZ_2(:,iL) = momDenZ_2(:,iL)./numberDen_2(:,iL)/Mass_2; % [m/s]
    
    EmeanX_2(:,iL) = me*Mass_2*numberDen_2(:,iL).*velX_2(:,iL).^2/2.0/qe; % [eV/m^3]
    EmeanY_2(:,iL) = me*Mass_2*numberDen_2(:,iL).*velY_2(:,iL).^2/2.0/qe; % [eV/m^3]
    EmeanZ_2(:,iL) = me*Mass_2*numberDen_2(:,iL).*velZ_2(:,iL).^2/2.0/qe; % [eV/m^3]
    
    tempX_2(:,iL) = 2*(eneDenX_2(:,iL) - EmeanX_2(:,iL))./numberDen_2(:,iL); % [eV]
    tempY_2(:,iL) = 2*(eneDenY_2(:,iL) - EmeanY_2(:,iL))./numberDen_2(:,iL); % [eV]
    tempZ_2(:,iL) = 2*(eneDenZ_2(:,iL) - EmeanZ_2(:,iL))./numberDen_2(:,iL); % [eV]
   
    tempAvg_2(:,iL) = (tempX_2(:,iL) + tempY_2(:,iL) + tempZ_2(:,iL))/3;
    presAvg_2(:,iL) = tempAvg_2(:,iL).*numberDen_2(:,iL);
    eintAvg_2(:,iL) = tempAvg_2(:,iL)/(gamma-1);   
    
    try
    momentFile_3 = [rootPath,'mesh_data/',species_folders(3).name, ...
                 '/',momentList(index(iL)).name]; 
    if(iL==1)
        time_scale = h5readatt(momentFile_3,'/species_data','time_scale_SI');
        Mass_3 = h5readatt(momentFile_3,'/species_data','mass');
        Charge_3 = double(h5readatt(momentFile_3,'/species_data','charge'));
    end
    
    %%%   reading density from species moment file
    groupName = '/species_data'; ghosts = 0;
    data = import1Ddata_singleFile(momentFile_3,groupName,ghosts);
    numberDen_3(:,iL) = data.Fcc(:,1);
    momDenX_3(:,iL) = data.Fcc(:,2)*cvac;
    momDenY_3(:,iL) = data.Fcc(:,3)*cvac;
    momDenZ_3(:,iL) = data.Fcc(:,4)*cvac;
    eneDenX_3(:,iL) = squeeze(data.Fcc(:,5))*cvac^2;
    eneDenY_3(:,iL) = squeeze(data.Fcc(:,6))*cvac^2;
    eneDenZ_3(:,iL) = squeeze(data.Fcc(:,7))*cvac^2;
    
    %%%   compute physical variables
    velX_3(:,iL) = momDenX_3(:,iL)./numberDen_3(:,iL)/Mass_3; % [m/s]
    velY_3(:,iL) = momDenY_3(:,iL)./numberDen_3(:,iL)/Mass_3; % [m/s]
    velZ_3(:,iL) = momDenZ_3(:,iL)./numberDen_3(:,iL)/Mass_3; % [m/s]
    
    tempX_3(:,iL) = 2*me/qe*(eneDenX_3(:,iL) - Mass_3*numberDen_3(:,iL).*velX_3(:,iL).^2/2.0)./numberDen_3(:,iL);
    tempY_3(:,iL) = 2*me/qe*(eneDenY_3(:,iL) - Mass_3*numberDen_3(:,iL).*velY_3(:,iL).^2/2.0)./numberDen_3(:,iL);
    tempZ_3(:,iL) = 2*me/qe*(eneDenZ_3(:,iL) - Mass_3*numberDen_3(:,iL).*velZ_3(:,iL).^2/2.0)./numberDen_3(:,iL);
   
    tempAvg_3(:,iL) = (tempX_3(:,iL) + tempY_3(:,iL) + tempZ_3(:,iL))/3;
    presAvg_3(:,iL) = tempAvg_3(:,iL).*numberDen_3(:,iL);
    eintAvg_3(:,iL) = tempAvg_3(:,iL)/(gamma-1);  
    end
    
    if(iL==1) 
        n0_1 = min(numberDen_1(:,iL)); 
        n0_2 = min(numberDen_2(:,iL)); %n0_3 = min(numberDen_3(:,iL));
        %n0 = 9.082e24; n0_2 = n0;
        %Mass_avg = (n0*Mass_1 + n0_2*Mass_2 + n0_3*Mass_3)/(n0 + n0_2 + n0_3);
        Mass_avg = (n0_1*Mass_1 + n0_2*Mass_2)/(n0_1 + n0_2);
        T0 = mean(tempAvg_2(1:nX/2,iL));
        T2_soln = T0*T2oT1;
        P0_1 = n0_1*T0;
        P0_2 = n0_2*T0;
        V0 = sqrt(qe*T0/(me*(Mass_avg))*gamma);
        %
        Clog = 10;
        mu = sqrt(Mass_2*me/mp);
        nui = 4.80e-8*Z^4*(n0_2*N2oN1/1e6)*Clog/(T2_soln)^1.5*sqrt(mu); % [Hz]
    end
    
    Xshift = Xfc(end); Xmin = -1*200; dXplot = abs(Xmin)/5;
    Xplot = Xcc - Xshift; Xfcplot = Xfc - Xshift; L = Xfc(end)-Xfc(1);
    subplot(2,2,1); hold off; axis('tight');
    p11=plot(Xplot,numberDen_2(:,iL)/n0_2,'displayName','ion'); box on; hold on;
    %p13=plot(Xplot,numberDen_3(:,iL)/n0,'displayName','ion 2'); box on; hold on;
    p12=plot(Xplot,numberDen_1(:,iL)/n0_1,'displayName','ele','linestyle','--');
    title('density'); grid on; legend('show','location','northwest');
    set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xmin Xfcplot(end)]); ylim([0 5]);
    xlabel('x/\lambda_i_i^D^S'); ylabel('n/n_0'); axis('square');
    hold on; l1=line([Xplot(1) Xplot(end)],[N2oN1 N2oN1],'linestyle','-.','color','black');

    %
    subplot(2,2,2); hold off; axis('tight');
    p21=plot(Xplot,velX_2(:,iL)/V0,'displayName','ion'); box on; hold on;
    %p22=plot(Xplot,velX_1(:,iL)/V0,'displayName','ele','linestyle','--');
    hold on; l2=line([Xplot(1) Xplot(end)],[pistonMach pistonMach],'linestyle','-.','color','black');
    title('x-velocity'); grid on; 
    set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xmin Xfcplot(end)]); ylim([0 Mach]);
    xlabel('x/\lambda_i_i^D^S'); ylabel('v/c_1'); axis('square');
    %
    subplot(2,2,3); hold off; axis('tight');
    p31=plot(Xplot,tempAvg_2(:,iL),'displayName','ion'); box on; hold on;
   %p33=plot(Xplot,tempAvg_3(:,iL),'displayName','ion 2'); box on; hold on;
    p32=plot(Xplot,tempAvg_1(:,iL),'displayName','ele','linestyle','--');
    title('temperature'); grid on; legend('show','location','best');
    set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xmin Xfcplot(end)]); ylim([0 14]);
    xlabel('x/\lambda_i_i^D^S'); ylabel('temperature [eV]'); axis('square'); Ymax = max(15,max(tempAvg_2(:,iL)));
    hold on; plot(Xplot,tempAvg_2(:,1),'black--'); ylim([0 Ymax]); set(gca,'ytick',0:3:Ymax)
    hold on; line([Xfcplot(1) Xfcplot(end)],[T2_soln T2_soln],'linestyle','-.','color','black');
    %
    subplot(2,2,4); hold off; axis('tight');
    p41=plot(Xplot,presAvg_2(:,iL)/P0_2,'displayName','ion'); box on; hold on;
    p42=plot(Xplot,presAvg_1(:,iL)/P0_1,'displayName','ele','linestyle','--');
    title('pressure'); grid on; legend('show','location','best'); Ymax = max(60,max(presAvg_2(:,iL)/P0_2));
    set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xmin Xfcplot(end)]); ylim([0 Ymax]);
    xlabel('x/\lambda_i_i^D^S'); ylabel('P/P_0'); axis('square');
    hold on; line([Xfcplot(1) Xfcplot(end)],[P2oP1 P2oP1],'linestyle','-.','color','black');
    %
%     subplot(2,2,4); hold off; axis('tight');
%     p41=plot(Xplot,EX(2:end-1,iL)); box on;
%     title('electric field [V/m]'); grid on; legend('show','location','best');
%     set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xmin Xfcplot(end)]); %ylim([-6e7 6e7]);
%     xlabel('x/\lambda_i_i^D^S'); ylabel('E [V/m]'); axis('square');

    % Create textbox
    t0 = length_scale/(6*V0); % time for piston to go 1 lambda_ii
    a1=annotation(f1,'textbox',...
    [0.444 0.935 0.13 0.043],...
    'String',{['t/t_0 = ',num2str(time(iL)*time_scale/t0,3)]},...
    'FitBoxToText','off','BackgroundColor',[1 1 0]);
    
    NameString = './figs/plasmShock1D';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    
    frame = getframe(gcf);
    %frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p11); delete(p12);
        delete(p21); %delete(p22);
        delete(p31); delete(p32);
        delete(p41); %delete(p42);
        delete(a1);
    end
    display(it);

end
close(v);

f7=figure(771); hold on;
plot(Xplot,(tempAvg_1(:,end)-tempAvg_2(:,end))/T2_soln); grid on;
xlim([-200 -70]); title('temp diff'); xlabel('x/\lambda_i_i^D^S');
ylabel('(T_e - T_i)/T_s_o_l_n^D^S'); box on;

f8=figure(881); hold on;
plot(time*time_scale/t0,tempAvg_2(end,:)/T2_soln,'displayName','ion'); hold on;
plot(time*time_scale/t0,tempAvg_1(end,:)/T2_soln,'displayName','ele');
grid on; box on; title('temperature at Xmax'); ylabel('T/T_s_o_l_n^D^S');
xlabel('t/t_0'); legend('show');


Pii_xx = numberDen_2.*(tempX_2-tempAvg_2); % ion viscosity tensor xx
Pii_yy = numberDen_2.*(tempY_2-tempAvg_2);
Pii_zz = numberDen_2.*(tempZ_2-tempAvg_2);
Pie_xx = numberDen_1.*(tempX_1-tempAvg_1); % ele viscosity tensor xx
Pie_yy = numberDen_1.*(tempY_1-tempAvg_1);
Pie_zz = numberDen_1.*(tempZ_1-tempAvg_1);
%
f9=figure(993); set(gcf,'position',[860 640 900 380]);
%
subplot(1,2,1);
plot(Xplot,presAvg_1(:,end)/P0_2,'displayName','scalar P'); hold on;
plot(Xplot,Pie_xx(:,end)/P0_2,'displayName','\Pi_x_x'); hold on;
grid on; box on; title('electron pressure'); ylabel('P/P_i_0');
xlabel('x/\lambda_i_i^D^S'); legend('show','location','best');
set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xfcplot(1) Xfcplot(end)]);
plot_ylim = get(gca,'ylim');
%
subplot(1,2,2);
plot(Xplot,presAvg_2(:,end)/P0_2,'displayName','scalar P'); hold on;
plot(Xplot,Pii_xx(:,end)/P0_2,'displayName','\Pi_x_x'); hold on;
grid on; box on; title('ion pressure'); ylabel('P/P_i_0');
xlabel('x/\lambda_i_i^D^S'); legend('show','location','best');
set(gca,'xtick',Xfcplot(1):dXplot:Xfcplot(end)); xlim([Xfcplot(1) Xfcplot(end)]);
set(gca,'ylim',plot_ylim);

%%%   compute divE and charge density

chargeDen = qe*(Charge_1*numberDen_1 + Charge_2*numberDen_2);
divE = zeros(length(Xfc),length(time));
for i=1:length(Xfc)
    divE(i,:) = ep0*(EX(i+1,:) - EX(i,:))/dX/length_scale; % [SI]
end

divE_cc = zeros(length(Xcc),length(time));
for i=1:length(Xcc)
    divE_cc(i,:) = ep0*(EX(i+2,:) - EX(i,:))/dX/2/length_scale; % [SI]
end

% [a,i0]=max(tempAvg_2(:,end));
% hold on; line([Xplot(i0)-0*sqrt(Mass_2/Mass_1) Xplot(i0)-0*sqrt(Mass_2/Mass_1)],[0 1e23]);
% hold on; line([Xplot(i0)-sqrt(Mass_2/Mass_1) Xplot(i0)-sqrt(Mass_2/Mass_1)],[0 1e23]);

