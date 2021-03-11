function [DataStr] = import2Ddata_singleFile(fileName,group,ghostCells)

if nargin<2
  group = 2;      % group 2 is /level_0
  ghostCells = 0; % do not keep ghost cells in output by default
end

if nargin<3 
  ghostCells = 0; % do not keep ghost cells in output by default
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   import 2D data from myPIC output files
%%%
%%%   input parameters:
%%%   fileName = 'string/to/plt_something_plots/thisFile' 
%%%   ghostCells = 0 (do not retain ghost cells in output)
%%%              = 1 (retain ghost cells in output)
%%%
%%%   output parameter is a structure with:
%%%   block structure with Xce, Yce, and Fce for each block
%%%   time found in output file
%%%   cell-edge grid (Xce and Yce) for all blocks
%%%   cell-center data (Fcc) for all blocks
%%%   number of components
%%%   number of processors
%%%   number of blocks
%%%   number of ghost cells at each end in X
%%%   number of ghost cells at each end in Y
%%%
%%%   block structure has local Xce, Yce, and Fce for each block, where 
%%%   Fce is cell centered data extended by one in each direction
%%%   to be compatible with Xce and Yce for pcolor and contour plots
%%%
%%%   Note that the function value is defined at cell-center
%%%   while the grid is defined at cell-edge
%%%
%%%   Last updated July 11, 2018
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisFile = fileName;
fileinfo = hdf5info(thisFile);
groupName = fileinfo.GroupHierarchy.Groups(group).Name;

% time = h5readatt(thisFile,'/level_0','time');
% dx = h5readatt(thisFile,'/level_0','dx');
prob_domain = h5readatt(thisFile,groupName,'prob_domain');
ghost = h5readatt(thisFile,[groupName,'/data_attributes'],'ghost');
outputGhost = h5readatt(thisFile,[groupName,'/data_attributes'],'outputGhost');
%nDims = length(fileinfo.GroupHierarchy.Groups(2).Attributes(5).Value.Data)/2;
nComps = h5readatt(thisFile,'/','num_components');
nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
ny = double(prob_domain.hi_j-prob_domain.lo_j+1);
procs = hdf5read(thisFile,[groupName,'/Processors']);
numProcs = length(procs);


%%%   loop over time and load the data
%
%vecData  = hdf5read(thisFile,[groupName,'/density:datatype=0']);
%offsets  = hdf5read(thisFile,[groupName,'/density:offsets=0']);
vecData  = hdf5read(thisFile,[groupName,'/data:datatype=0']);
offsets  = hdf5read(thisFile,[groupName,'/data:offsets=0']);
boxes    = hdf5read(thisFile,[groupName,'/boxes']);
for iP = 1:numProcs
    lo_i(iP) = boxes(iP).Data{1} + 1 - ghost.intvecti;
    lo_j(iP) = boxes(iP).Data{2} + 1 - ghost.intvectj;       
    hi_i(iP) = boxes(iP).Data{3} + 1 + ghost.intvecti;
    hi_j(iP) = boxes(iP).Data{4} + 1 + ghost.intvectj;
end
%display(lo_j(2));

%%%   some indices are negative, so shift for MATLAB
%
min_lo_i = min(lo_i);
lo_i = lo_i-double(min_lo_i)+1;
hi_i = hi_i-double(min_lo_i)+1;
min_lo_j = min(lo_j);
lo_j = lo_j-double(min_lo_j)+1;
hi_j = hi_j-double(min_lo_j)+1;
    
    
% thisFileMap = [fileName(1:end-5),'.map.hdf5'];
% fileinfoMap = hdf5info(thisFileMap);
% vecMap      = hdf5read(thisFileMap,'/level_0/data:datatype=0');
% offsetsMap  = hdf5read(thisFileMap,'/level_0/data:offsets=0');
% boxesMap    = hdf5read(thisFileMap,'/level_0/boxes');
% lo_i_Map = zeros(1,numProcs);
% lo_j_Map = zeros(1,numProcs);
% hi_i_Map = zeros(1,numProcs);
% hi_j_Map = zeros(1,numProcs);
% for iP = 1:numProcs
%     lo_i_Map(iP) = boxesMap(iP).Data{1} + 1 - ghost.intvecti;
%     lo_j_Map(iP) = boxesMap(iP).Data{2} + 1 - ghost.intvectj;       
%     hi_i_Map(iP) = boxesMap(iP).Data{3} + 1 + ghost.intvecti;
%     hi_j_Map(iP) = boxesMap(iP).Data{4} + 1 + ghost.intvectj;
% end
% min_lo_i_Map = min(lo_i_Map);
% lo_i_Map = lo_i_Map-double(min_lo_i_Map)+1;
% hi_i_Map = hi_i_Map-double(min_lo_i_Map)+1;
% min_lo_j_Map = min(lo_j_Map);
% lo_j_Map = lo_j_Map-double(min_lo_j_Map)+1;
% hi_j_Map = hi_j_Map-double(min_lo_j_Map)+1;
    

%%%   Note that boxes and boxesMap are always the same, while   %%%
%%%   data is at cell center and grid at cell-edges             %%%

if(ghostCells)
    ig = 0;
    jg = 0;
else
    ig = ghost.intvecti;
    jg = ghost.intvectj;
end


%%%   map to reshaped grid
%   
map0cc = zeros(1,1);
for m=1:numProcs
    
%     gridMap = vecMap(offsetsMap(m)+1:offsetsMap(m+1));
%     thisMapX = gridMap(1:end/2);
%     thisMapY = gridMap(end/2+1:end);

%     %%%   formulate grid at cell edges
%     %
%     i0 = lo_i_Map(m);
%     i1 = hi_i_Map(m)+1;
%     nXsub = i1-i0+1;
%     j0 = lo_j_Map(m);
%     j1 = hi_j_Map(m)+1;
%     nYsub = j1-j0+1;
%     thisX(i0:i1,j0:j1) = reshape(thisMapX,nXsub,nYsub);
%     thisY(i0:i1,j0:j1) = reshape(thisMapY,nXsub,nYsub);
%     %
%     X(i0+ig:i1-ig,j0+jg:j1-jg) = thisX(i0+ig:i1-ig,j0+jg:j1-jg);
%     Y(i0+ig:i1-ig,j0+jg:j1-jg) = thisY(i0+ig:i1-ig,j0+jg:j1-jg);
    
    
    %%%   formulate function matrix on cell-center grid
    %
    i0data = offsets(m)+1;
    i1data = offsets(m+1);
    i0 = lo_i(m);
    i1 = hi_i(m);
    nXsub = i1-i0+1;
    j0 = lo_j(m);
    j1 = hi_j(m);
    nYsub = j1-j0+1;
   % size(vecData(i0data:i1data))

    thisvecData = vecData(i0data:i1data);
    subSize = length(thisvecData)/nComps;
    for n=1:nComps
        nlow = 1+(n-1)*subSize;
        nup = n*subSize;
        thisvecData0 = thisvecData(nlow:nup);
%         size(thisvecData0)
%         nXsub
%         nYsub
        thisdataPatch = reshape(thisvecData0,nXsub,nYsub);
        data0cc(i0+ig:i1-ig,j0+jg:j1-jg,n) = ...
                thisdataPatch(1+ig:end-ig,1+jg:end-jg); 
    end
    map0cc(i0+ig:i1-ig,j0+jg:j1-jg) = ...
           ones(size(squeeze(data0cc(i0+ig:i1-ig,j0+jg:j1-jg,1))));

end
%figure(100); contourf(data0cc);
%display(size(data0cc));

%%%   use binary map0cc to determine indices for different blocks
%
nx = length(data0cc(:,1,1));
ny = length(data0cc(1,:,1));
thisBox = 1;
for i=1:nx
    for j=1:ny
        if(map0cc(i,j)==1) % get lower indices for thisBox
            i0(thisBox) = i;
            j0(thisBox) = j;
            %
            for j2=j:ny    % get upper y index for thisBox
                if(map0cc(i,j2)==0)
                    j1(thisBox) = j2;
                    break;
                end
                if(j2==ny)  
                    j1(thisBox) = j2+1;
                    break;
                end
            end
            for i2=i:nx    % get upper x index for thisBox
                if(map0cc(i2,j)==0)
                    i1(thisBox) = i2;
                    break;
                end
                if(i2==nx)
                    i1(thisBox) = i2+1;
                    break;
                end
            end    
            
            %%%   zero out thisBox so don't find it twice
            %
            i00 = i0(thisBox); i11 = i1(thisBox);
            j00 = j0(thisBox); j11 = j1(thisBox);
            map0cc(i00:i11,j00:j11) = zeros(i11-i00+1,j11-j00+1);
            thisBox = thisBox+1;
        end
    end
end
%display(i0); display(i1);
%display(j0); display(j1);


%%%  write grid and data in output structure for each block
%
numBlocks  = length(i0); % number of blocks
for b=1:numBlocks
    nx = i1(b)-i0(b);
    ny = j1(b)-j0(b);
     
    %%% extend cc data by one in each dim
    %%% for compatability with ce grid
    %
    data0ce = zeros(nx+1,ny+1,nComps);
    data0ce(1:end-1,1:end-1,:) = data0cc(i0(b):i1(b)-1,j0(b):j1(b)-1,:);
    data0ce(end,:,:) = data0ce(end-1,:,:);
    data0ce(:,end,:) = data0ce(:,end-1,:);

    
 %   DataStr.block(b).Xce = X(i0(b):i1(b),j0(b):j1(b));
 %   DataStr.block(b).Yce = Y(i0(b):i1(b),j0(b):j1(b));
    DataStr.block(b).Fce = data0ce;
    DataStr.block(b).i0 = i0(b);
    DataStr.block(b).i1 = i1(b);
    DataStr.block(b).j0 = j0(b);
    DataStr.block(b).j1 = j1(b);
end


%%%   write output information to a data structure
%
%DataStr.time = time;     % time vector
%DataStr.Xce = X;         % X-grid at cell edges
%DataStr.Yce = Y;         % Y-grid at cell edges
DataStr.Fcc = squeeze(data0cc(i0(b):i1(b)-1,j0(b):j1(b)-1,:));   % data values at cell-center
DataStr.nComps = nComps;
DataStr.numProcs = numProcs;   % number of processors
DataStr.numBlocks = numBlocks;
%DataStr.numGhostX = ig; %ghost.intvecti;
%DataStr.numGhostY = jg; %ghost.intvectj;
if(ghostCells)
    DataStr.numGhostX = ghost.intvecti;
    DataStr.numGhostY = ghost.intvectj;
else
    DataStr.numGhostX = 0;
    DataStr.numGhostY = 0;
end

end


