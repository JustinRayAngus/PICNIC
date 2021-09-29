%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D test of exchange for NodeFArrayBox
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

rootPath = '../fromQuartz/fieldTests/2D/testingNode3x3/'; thisFig = 3;
ghosts = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh_data/mesh.h5'];

groupName = '/cell_centered_grid';
data1 = import2Ddata_singleFile(meshFile,groupName,1);
nGX = data1.numGhostX; nGY = data1.numGhostY;
Xcc = squeeze(data1.Fcc(:,:,1)); nX = length(Xcc(:,1))-2*nGX; dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data1.Fcc(:,:,2)); nZ = length(Zcc(1,:))-2*nGY; dZ = Zcc(1,2)-Zcc(1,1);

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



groupName = '/node_centered_test'; patchData = 1;
data5 = import2Ddata_singleFile(meshFile,groupName,ghosts,patchData);
Fnc = data5.Fnc(:,:); 
patchData = data5.patchData;
numProcs = length(patchData.proc);

close(figure(thisFig)); 
f3=figure(thisFig); set(f3,'position',[830 300 960 700]);

if(numProcs==4)
    procPlotMapping = [3 4 1 2];
    procTitleMapping = [0 2 1 3];
end
if(numProcs==9)
   procPlotMapping = [7 8 9 4 5 6 1 2 3];  
   procTitleMapping = [0 3 6 1 4 7 2 5 8];
end

for n=1:numProcs
    proc_data = patchData.proc(n).data;
    proc_data(end+1,:) = proc_data(end,:);
    proc_data(:,end+1) = proc_data(:,end);
    
    subplot(sqrt(numProcs),sqrt(numProcs),procPlotMapping(n)); 
    pcolor(proc_data); colorbar; caxis([0 numProcs-1]);
    title(['proc ',num2str(procTitleMapping(n))]);
end
    
    
% data_proc0 = patchData.proc(1).data; 
% data_proc0(end+1,:) = data_proc0(end,:);
% data_proc0(:,end+1) = data_proc0(:,end);
% data_proc1 = patchData.proc(3).data;
% data_proc1(end+1,:) = data_proc1(end,:);
% data_proc1(:,end+1) = data_proc1(:,end);
% data_proc2 = patchData.proc(2).data;
% data_proc2(end+1,:) = data_proc2(end,:);
% data_proc2(:,end+1) = data_proc2(:,end);
% data_proc3 = patchData.proc(4).data;
% data_proc3(end+1,:) = data_proc3(end,:);
% data_proc3(:,end+1) = data_proc3(:,end);
% 
% close(figure(333)); 
% f3=figure(333); set(f3,'position',[830 300 960 700]);
% subplot(2,2,4); pcolor(data_proc2); colorbar; title('proc 2');
% subplot(2,2,2); pcolor(data_proc3); colorbar; title('proc 3');
% subplot(2,2,3); pcolor(data_proc0); colorbar; title('proc 0');
% subplot(2,2,1); pcolor(data_proc1); colorbar; title('proc 1');
% 



%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


