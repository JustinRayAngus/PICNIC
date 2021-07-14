function [DataStr] = import1Ddata_singleFile( fileName, ...
                                              groupName, ...
                                              ghostCells, ...
                                              writePatchData )

if nargin<4
  writePatchData = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   import 1D data from PICNIC output files
%%%
%%%   input parameters:
%%%   fileName = 'string/to/plt_something_plots/thisFile'
%%%   groupName = '/goupName'
%%%   ghostCells = 0 (do not retain ghost cells in output)
%%%              = 1 (retain ghost cells in output)
%%%   writePatchData = 0 (do not retain structure of all data on each
%%%                       proccessor in output)
%%%   writePatchData = 1 (do retain structure of all data on each
%%%                       processor in output)
%%%
%%%   Last updated March 2021
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisFile = fileName;
fileinfo = hdf5info(thisFile);
GH = fileinfo.GroupHierarchy;
groupName0 = GH.Groups(1).Name; % /Chombo_global
SpaceDim = h5readatt(thisFile,groupName0,'SpaceDim');

numGroups = length(GH.Groups);
group = 0;
for n=1:numGroups
    thisGroupName = GH.Groups(n).Name;
    if(strcmp(thisGroupName,groupName))
        group = n;
        break;
    end
end
if(group==0)
    for n=1:numGroups
        thisGroupName = GH.Groups(n).Name;
        %display(thisGroupName);
    end
    error(['groupName = ',groupName,' not found']);
end
%group = 2;
%groupName = fileinfo.GroupHierarchy.Groups(group).Name

% time = h5readatt(thisFile,'/level_0','time');
% dx = h5readatt(thisFile,'/level_0','dx');
prob_domain = h5readatt(thisFile,groupName,'prob_domain');
ghost = h5readatt(thisFile,[groupName,'/data_attributes'],'ghost');
outputGhost = h5readatt(thisFile,[groupName,'/data_attributes'],'outputGhost');

nComps = h5readatt(thisFile,groupName,'num_components');
procs = hdf5read(thisFile,[groupName,'/Processors']);
numProcs = length(procs);

ig = ghost.intvecti;

vecData  = hdf5read(thisFile,[groupName,'/data:datatype=0']);
offsets  = hdf5read(thisFile,[groupName,'/data:offsets=0']);
boxes    = hdf5read(thisFile,[groupName,'/boxes']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   check if this is a fluxBox
%%%

isCellBox = 1;
isFluxBox = 0;
isEdgeBox = 0;
isNodeBox = 0;

try
    isFluxBox = h5readatt(thisFile,groupName,'is_fluxbox');
catch
    %display('not a fluxbox');
end
isCellBox = isCellBox - isFluxBox;

try
    isEdgeBox = h5readatt(thisFile,groupName,'is_edgebox');
catch
    %display('not a edgebox');
end
isCellBox = isCellBox - isEdgeBox;

try
    isNodeBox = h5readatt(thisFile,groupName,'is_nodebox');
catch
    %display('not a edgebox');
end
isCellBox = isCellBox - isNodeBox;

assert(isCellBox+isEdgeBox+isFluxBox+isNodeBox==1);

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isCellBox || isEdgeBox)
    nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
else
    prob_domain.hi_i = prob_domain.hi_i + 1;
    nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
end


for iP = 1:numProcs
    lo_i(iP) = boxes(iP).Data{1} + 1 - ig;
    hi_i(iP) = boxes(iP).Data{2} + 1 + ig ...
             + isFluxBox ...
             + isNodeBox;
end

%%%   some indices are negative, so shift for MATLAB
%
min_lo_i = double(min(lo_i));
lo_i = lo_i - min_lo_i + 1;
hi_i = hi_i - min_lo_i + 1;


%%%   map to reshaped grid
%
isStag = isFluxBox + isEdgeBox;
totalComps = nComps + isStag*nComps*(SpaceDim-1);
data0cc = zeros(max(nx)+2*ig,totalComps);
size(data0cc);
for d=1
    for m=1:numProcs

        i0data = offsets(m)+1;
        i1data = offsets(m+1);
                
        i0 = lo_i(m);
        i1 = hi_i(m);

        thisvecData = vecData(i0data:i1data);
        subSize = length(thisvecData)/nComps;
        if(isStag)
            subSize = subSize/SpaceDim;
        end
        for n=1:nComps
            N = n + nComps*(d-1);
            nlow = 1+(N-1)*subSize;
            nup = N*subSize;
            thisvecData0 = thisvecData(nlow:nup);
            thisdataPatch = thisvecData0;
            if(isCellBox || isNodeBox)
                patchData.proc(m).data(:,N) = thisdataPatch;
            else
                patchData.proc(m).dir(d).data(:,n) = thisdataPatch;
            end
            data0cc(i0:i1,N) = thisdataPatch(1:end);
        end
        if(isCellBox || isNodeBox)
            patchData.proc(m).i0 = i0;
            patchData.proc(m).i1 = i1;
        else
            patchData.proc(m).dir(d).i0 = i0;
            patchData.proc(m).dir(d).i1 = i1;        
        end

    end
end
size(data0cc);

% The following is needed to replace some internal cells with valid
% data from ghost cell data that was used above. This is only needed
% for the case where exchange() was not called on the dataset such
% that the data in the interior ghost cells does not match the corresponding
% valid data of the neighboring processor

if(isCellBox || isNodeBox)
    for m=1:numProcs
        thisdata = patchData.proc(m).data;
        i0 = patchData.proc(m).i0 + ig;
        i1 = patchData.proc(m).i1 - ig;
        data0cc(i0:i1,n) = thisdata(1+ig:end-ig,n);
    end
else
    for d=1:SpaceDim
        for m=1:numProcs
            N = n + nComps*(d-1);
            thisdata = patchData.proc(m).dir(d).data;
            i0 = patchData.proc(m).dir(d).i0 + ig;
            i1 = patchData.proc(m).dir(d).i1 - ig;
            data0cc(i0:i1,N) = thisdata(1+ig:end-ig,n);
        end
    end
end

%%%   write output information to a data structure
%
if(writePatchData)
    patchData.nComps = nComps;
    DataStr.patchData = patchData; % contains all data on each processor
end

if(ghostCells)
    ig = 0;
end

if(isFluxBox)
    DataStr.Ffc0 = squeeze(data0cc(1+ig:end-ig,:));
elseif(isEdgeBox)
    DataStr.Fec0 = squeeze(data0cc(1+ig:end-ig,:)); 
elseif(isNodeBox)
    DataStr.Fnc = squeeze(data0cc(1+ig:end-ig,:)); 
else
    DataStr.Fcc = squeeze(data0cc(1+ig:end-ig,:));
end       
     
DataStr.isCellBox = isCellBox;
DataStr.isFluxBox = isFluxBox;
DataStr.isEdgeBox = isEdgeBox;
DataStr.isNodeBox = isNodeBox;

DataStr.numComps = nComps;
DataStr.numProcs = numProcs;

if(ghostCells)
    DataStr.numGhostX = ghost.intvecti;
else
    DataStr.numGhostX = 0;
end

fclose('all');

end


