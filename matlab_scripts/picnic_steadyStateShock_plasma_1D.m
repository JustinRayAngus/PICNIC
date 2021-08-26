%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D steady state shock in plasma using picnic
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
ep0  = 8.85418782e-12;
cvac = 2.99792458e8;       % speed of light [m/s]

gamma = 5/3;

rootPath = '../fromQuartz/1D/steadyStateShock/plasma/test0/';

basePath = '../fromQuartz/1D/steadyStateShock/plasma/uniformMovingFrame/explicit/';
%basePath = '../fromQuartz/1D/steadyStateShock/plasma/uniformMovingFrame/implicit/';
rootPath = [basePath,'test1_periodic/'];
%rootPath = [basePath,'test0_labFrame/'];
rootPath = [basePath,'test1/'];

basePath = '../fromQuartz/1D/steadyStateShock/plasma/counterBeams/';
rootPath = [basePath,'L400/Np200_dt0p025/'];
%rootPath = [basePath,'L400/Np200_dt0p050/'];
rootPath = [basePath,'L400/Np200_dt0p100/'];
rootPath = [basePath,'L400/test_symmetry/'];
%rootPath = [basePath,'L400/Np400_dt0p050/'];
%rootPath = [basePath,'L400/Np800_dt0p050/'];
%
%rootPath = [basePath,'L800/Np200_dt0p050/'];
%rootPath = [basePath,'L800/Np400_dt0p050/'];
%
%rootPath = [basePath,'L1600/Np200_dt0p050/'];
%rootPath = [basePath,'L1600/Np400_dt0p050/'];

%%% see planaShockJump.m
%
T2_soln = 5.3716;
P2oP1_soln = 82.2069; 
N2oN1_soln = 3.826;


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

%partList_1   = dir([rootPath,'particle_data/species1_data/parts*']);
%partList_2   = dir([rootPath,'particle_data/species2_data/parts*']);
momentList_1 = dir([rootPath,'mesh_data/species1_data/moments*']);
momentList_2 = dir([rootPath,'mesh_data/species2_data/moments*']);
ListLength = length(momentList_1);
assert(ListLength==length(momentList_2));

fieldList = dir([rootPath,'mesh_data/field_data/field*']);
ListLength_fields = length(fieldList);
assert(ListLength_fields==ListLength)

step = zeros(size(momentList_1));
index = zeros(size(momentList_1));
for n=1:ListLength
    thisFile = fieldList(n).name;
    step(n) = str2num(thisFile(7:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
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
eneDenX_1 = zeros(nX,iLmax); tempX_1 = zeros(nX,iLmax);
eneDenY_1 = zeros(nX,iLmax); tempY_1 = zeros(nX,iLmax);
eneDenZ_1 = zeros(nX,iLmax); tempZ_1 = zeros(nX,iLmax);
%
tempAvg_1 = zeros(nX,iLmax); presAvg_1 = zeros(nX,iLmax); eintAvg_1 = zeros(nX,iLmax);
%
%totalParts_2 = zeros(1,iLmax);
numberDen_2 = zeros(nX,iLmax);
momDenX_2 = zeros(nX,iLmax); velX_2 = zeros(nX,iLmax);
momDenY_2 = zeros(nX,iLmax); velY_2 = zeros(nX,iLmax);
momDenZ_2 = zeros(nX,iLmax); velZ_2 = zeros(nX,iLmax);
eneDenX_2 = zeros(nX,iLmax); tempX_2 = zeros(nX,iLmax);
eneDenY_2 = zeros(nX,iLmax); tempY_2 = zeros(nX,iLmax);
eneDenZ_2 = zeros(nX,iLmax); tempZ_2 = zeros(nX,iLmax);
%
tempAvg_2 = zeros(nX,iLmax); presAvg_2 = zeros(nX,iLmax); eintAvg_2 = zeros(nX,iLmax);
%

%%%  loop over files and create movie
%
close(figure(1));
f1=figure(1); set(f1,'position',[860 240 900 760]);
set(gcf,'color','white');
%
images = cell(1,1);
v=VideoWriter('./figs/plasmaShock.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

dt_out = 1;
if ( (dt_out==2) && (floor(iLmax/dt_out)==iLmax/dt_out) )
   iLmax = iLmax-1; 
end
%iLmax = 16; %iLmax = 141;
display(iLmax);
%for iL = 1:dt_out:iLmax
for iL = [1,iLmax]   

    %
    %   load field data
    %
    fieldFile = [rootPath,'mesh_data/field_data/',fieldList(index(iL)).name];
    
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
    
%     partsFile_1 = [rootPath,'particle_data/species1_data/', ...
%                  partList_1(index(iL)).name];
%     fileinfo = hdf5info(partsFile_1);
%     fileinfo.GroupHierarchy.Groups(2).Attributes.Name;
%     
%     partData = hdf5read(partsFile_1,'/species_data/particles:data');
%     SpaceDim = h5readatt(partsFile_1,'/Chombo_global','SpaceDim');
%     numParts = h5readatt(partsFile_1,'/species_data','num_particles');
%     time(iL) = h5readatt(partsFile_1,'/species_data','time');
%     if(iL==1)
%         numPartComps_1 = h5readatt(partsFile_1,'/species_data','numPartComps');
%     end
%     partData = reshape(partData,numPartComps_1,numParts);
%     partData = partData';
%     totalParts_1(iL) = numParts;
% 
%     particle.weight = partData(:,1);
%     particle.x    = partData(:,2);
%     particle.y    = partData(:,3);
%     particle.z    = partData(:,4); 
%     particle.vx   = partData(:,5)*cvac;
%     particle.vy   = partData(:,6)*cvac;
%     particle.vz   = partData(:,7)*cvac;
%     particle.ID   = partData(:,numPartComps_1);

    momentFile_1 = [rootPath,'mesh_data/species1_data/',momentList_1(index(iL)).name];
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
    eneDenX_1(:,iL) = squeeze(data.Fcc(:,5))*cvac^2; tempX_1(:,iL) = 2*(eneDenX_1(:,iL)*me - me*Mass_1*numberDen_1(:,iL).*velX_1(:,iL).^2/2.0)./numberDen_1(:,iL)/qe;
    eneDenY_1(:,iL) = squeeze(data.Fcc(:,6))*cvac^2;
    eneDenZ_1(:,iL) = squeeze(data.Fcc(:,7))*cvac^2;
    
    %%%   compute physical variables
    
    velX_1(:,iL) = momDenX_1(:,iL)./numberDen_1(:,iL)/Mass_1; % [m/s]
    velY_1(:,iL) = momDenY_1(:,iL)./numberDen_1(:,iL)/Mass_1; % [m/s]
    velZ_1(:,iL) = momDenZ_1(:,iL)./numberDen_1(:,iL)/Mass_1; % [m/s]
    
    tempX_1(:,iL) = 2*me/qe*(eneDenX_1(:,iL) - Mass_1*numberDen_1(:,iL).*velX_1(:,iL).^2/2.0)./numberDen_1(:,iL);
    tempY_1(:,iL) = 2*me/qe*(eneDenY_1(:,iL) - Mass_1*numberDen_1(:,iL).*velY_1(:,iL).^2/2.0)./numberDen_1(:,iL);
    tempZ_1(:,iL) = 2*me/qe*(eneDenZ_1(:,iL) - Mass_1*numberDen_1(:,iL).*velZ_1(:,iL).^2/2.0)./numberDen_1(:,iL);
   
    tempAvg_1(:,iL) = (tempX_1(:,iL) + tempY_1(:,iL) + tempZ_1(:,iL))/3;
    presAvg_1(:,iL) = tempAvg_1(:,iL).*numberDen_1(:,iL);
    eintAvg_1(:,iL) = tempAvg_1(:,iL)/(gamma-1);
    
    %
    %   load data for species2
    %
    
%     partsFile_2 = [rootPath,'particle_data/species2_data/', ...
%                  partList_2(index(iL)).name];
%     fileinfo = hdf5info(partsFile_2);
%     fileinfo.GroupHierarchy.Groups(2).Attributes.Name;
    
%     partData = hdf5read(partsFile_2,'/species_data/particles:data');
%     SpaceDim = h5readatt(partsFile_2,'/Chombo_global','SpaceDim');
%     numParts = h5readatt(partsFile_2,'/species_data','num_particles');
%     time(iL) = h5readatt(partsFile_2,'/species_data','time');
%     if(iL==1)
%         Mass_2 = h5readatt(partsFile_2,'/species_data','mass');
%         Charge_2 = double(h5readatt(partsFile_2,'/species_data','charge'));
%         Uint_2 = h5readatt(partsFile_2,'/species_data','Uint'); % [eV]
%         numPartComps_2 = h5readatt(partsFile_2,'/species_data','numPartComps');
%     end
%     partData = reshape(partData,numPartComps_2,numParts);
%     partData = partData';
%     totalParts_2(iL) = numParts;
% 
%     particle.weight = partData(:,1);
%     particle.x    = partData(:,2);
%     particle.y    = partData(:,3);
%     particle.z    = partData(:,4); 
%     particle.vx   = partData(:,5)*cvac;
%     particle.vy   = partData(:,6)*cvac;
%     particle.vz   = partData(:,7)*cvac;
%     particle.ID   = partData(:,numPartComps_1);

    momentFile_2 = [rootPath,'mesh_data/species2_data/',momentList_2(index(iL)).name];
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
    eneDenX_2(:,iL) = squeeze(data.Fcc(:,5))*cvac^2; tempX_2(:,iL) = 2*(eneDenX_2(:,iL)*me - me*Mass_2*numberDen_2(:,iL).*velX_2(:,iL).^2/2.0)./numberDen_2(:,iL)/qe;
    eneDenY_2(:,iL) = squeeze(data.Fcc(:,6))*cvac^2;
    eneDenZ_2(:,iL) = squeeze(data.Fcc(:,7))*cvac^2;
    
    %%%   compute physical variables
    
    velX_2(:,iL) = momDenX_2(:,iL)./numberDen_2(:,iL)/Mass_2; % [m/s]
    velY_2(:,iL) = momDenY_2(:,iL)./numberDen_2(:,iL)/Mass_2; % [m/s]
    velZ_2(:,iL) = momDenZ_2(:,iL)./numberDen_2(:,iL)/Mass_2; % [m/s]
    
    tempX_2(:,iL) = 2*me/qe*(eneDenX_2(:,iL) - Mass_2*numberDen_2(:,iL).*velX_2(:,iL).^2/2.0)./numberDen_2(:,iL);
    tempY_2(:,iL) = 2*me/qe*(eneDenY_2(:,iL) - Mass_2*numberDen_2(:,iL).*velY_2(:,iL).^2/2.0)./numberDen_2(:,iL);
    tempZ_2(:,iL) = 2*me/qe*(eneDenZ_2(:,iL) - Mass_2*numberDen_2(:,iL).*velZ_2(:,iL).^2/2.0)./numberDen_2(:,iL);
   
    tempAvg_2(:,iL) = (tempX_2(:,iL) + tempY_2(:,iL) + tempZ_2(:,iL))/3;
    presAvg_2(:,iL) = tempAvg_2(:,iL).*numberDen_2(:,iL);
    eintAvg_2(:,iL) = tempAvg_2(:,iL)/(gamma-1);   
    
    
    if(iL==1) 
        n0 = min(numberDen_2(:,iL));
        T0 = mean(tempAvg_2(1:nX/2,iL));
        P0 = n0*T0;
        V0 = sqrt(qe*2*T0/(me*(Mass_1+Mass_2))*gamma);
        %
        Z = 1;
        Clog = 10;
        T2_eV = 5.4;
        nui = 4.80e-8*Z^4*(n0/1e6)*Clog/(T2_eV)^1.5; % [Hz]
    end
    
    subplot(2,2,1); hold off; axis('tight');
    p11=plot(Xcc-X0/2,numberDen_2(:,iL)/n0,'displayName','ion'); box on; hold on;
    p12=plot(Xcc-X0/2,numberDen_1(:,iL)/n0,'displayName','ele','linestyle','--');
    title('density'); grid on; legend('show','location','northwest');
    set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([0 5]);
    xlabel('x/\lambda_i_i'); ylabel('n/n_0'); axis('square');
    hold on; line([Xcc(1)-X0/2 Xcc(end)-X0/2],[N2oN1_soln N2oN1_soln],'linestyle','-.','color','black');

    %
    subplot(2,2,2); hold off; axis('tight');
    p21=plot(Xcc-X0/2,velX_2(:,iL)/V0,'displayName','ion'); box on; hold on;
    %p22=plot(Xcc,velX_1(:,iL)/V0,'displayName','ele','linestyle','--');
    title('x-velocity'); grid on; 
    set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([-8 8]);
    xlabel('x/\lambda_i_i'); ylabel('v/c_1'); axis('square');
    %
    subplot(2,2,3); hold off; axis('tight');
    p31=plot(Xcc-X0/2,tempAvg_2(:,iL),'displayName','ion'); box on; hold on;
    p32=plot(Xcc-X0/2,tempAvg_1(:,iL),'displayName','ele','linestyle','--');
    title('temperature'); grid on; legend('show','location','best');
    set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([0 7]);
    xlabel('x/\lambda_i_i'); ylabel('temperature [eV]'); axis('square');
    hold on; plot(Xcc-X0/2,tempAvg_2(:,1),'black--');
    hold on; line([Xcc(1)-X0/2 Xcc(end)-X0/2],[T2_soln T2_soln],'linestyle','-.','color','black');
    %
%     subplot(2,2,4); hold off; axis('tight');
%     p41=plot(Xcc-X0/2,presAvg_2(:,iL)/P0,'displayName','ion'); box on; hold on;
%     p42=plot(Xcc-X0/2,presAvg_1(:,iL)/P0,'displayName','ele','linestyle','--');
%     title('pressure'); grid on; legend('show','location','best');
%     set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([0 100]);
%     xlabel('x/\lambda_i_i'); ylabel('P/P_0'); axis('square');
%     hold on; line([Xcc(1)-X0/2 Xcc(end)-X0/2],[P2oP1_soln P2oP1_soln],'linestyle','-.','color','black');
    %
    subplot(2,2,4); hold off; axis('tight');
    p41=plot(Xcc-X0/2,EX(2:end-1,iL)); box on;
    title('electric field [V/m]'); grid on; legend('show','location','best');
    set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([-6e7 6e7]);
    xlabel('x/\lambda_i_i'); ylabel('E [V/m]'); axis('square');
    %
%     subplot(2,2,4); hold off;
%     p41=plot(Xcc,JX(2:end-1,iL)); box on;
%     title('current density [A/m^2]'); grid on; legend('show','location','northwest');
%     set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([-1e10 1e10]);
%     xlabel('x/\lambda_i_i'); ylabel('J [A/m^2]'); axis('square');

%     subplot(2,2,4); hold off;
%     p41=plot(Xcc,JX(2:end-1,iL).*EX(2:end-1,iL)); box on;
%     title('J*E'); grid on; legend('show','location','northwest');
%     set(gca,'xtick',-X0/2:80:X0/2); xlim([-240 240]); ylim([-1e17 1e17]);
%     xlabel('x/\lambda_i_i'); ylabel('J*E'); axis('square');

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
    display(iL);

end
close(v);

% close(figure(7));
% f7=figure(7); 
% 
% plot(Xcc,tempAvg_2(:,end),'displayName','ion'); hold on;
% plot(Xcc,tempAvg_1(:,end),'displayName','ele'); hold off;
% xlabel('x/\lambda_i_i'); ylabel('T [eV]'); title('temperatures');
% hold on; plot(Xcc,tempAvg_2(:,1),'black--');
% hold on; line([Xcc(1) Xcc(end)],[T2_soln T2_soln],'linestyle','-.','color','black');


for i=1:nX
    for n=1:iLmax
        thisDen = numberDen_1(i,n);
        if(thisDen==0)
            velX_1(i,n) = 0;
            velY_1(i,n) = 0;
            velZ_1(i,n) = 0;
            tempX_1(i,n) = 0;
            tempY_1(i,n) = 0;
            tempZ_1(i,n) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%
%%%
%%%   bin up the particles in velocity space and look at distribution
%%%
%%%%%%%%%%%%%%%%%%

% Vmax = max(particle.vy);
% Vmin = min(particle.vy);
% Vgrid = linspace(Vmin,Vmax,100);
% dVy = Vgrid(2)-Vgrid(1);
% 
% dfn = zeros(size(Vgrid));
% for i=1:numParts
%     thisv = particle.vy(i);
%     [~,index]=min(abs(thisv-Vgrid));
%     dfn(index) = dfn(index) + 1;
% end
% normC = sum(dfn)*dVy;
% dfn = dfn/normC;

% close(figure(44)); f4=figure(44); 
% 
% TY0 = mean(nonzeros(tempY_1(:,iLmax)));
% UY0 = mean(nonzeros(velY_1(:,iLmax)));
% VTY = 4.19e5*sqrt(TY0/Mass_1);
% plot(Vgrid,exp(-((Vgrid-UY0)/sqrt(2)/VTY).^2)/sqrt(2*pi)/VTY); hold on;
% plot(Vgrid,dfn); 
% xlabel('y-velocity'); ylabel('dfn-vy');
% title('vy-distribution function');
% legend('maxwellian','data'); axis('tight');


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

