function BA_Rep_Thresh()
% function to report multiple imgs brain area
% function based on xjview
% written by hongshengcheng.math@gmail.com
% 
% v7:20141109,fix when img not supass thresh and not report the res
% v6:20141101,fix error when report the mean cluster value and peak coordinate
% v5:20141030,when threshold is positive then get the binearal gap
% v4: 20140918,support 4D img
% v3: add mean cluster value and peak cluster value&coord
% created date: 20140828


% set up paras
para_dlg = inputdlg({'Reslut Dir Name','cluster size',...
    'threshold (e.g.0)'},...
    'Please input', ...
    [1 50;1 50;1 50],{'BA_Report_Res','10','0'}) ;

ResDir = para_dlg{1};
mkdir(ResDir);
clustersize = str2num(para_dlg{2});
threshold = str2num(para_dlg{3});

% select multiple 3D-imgs or just 1 4D-img
ImgList = ScanImg2Cell('Choose Img Files');
save(fullfile(pwd,ResDir,'ImgFileList.mat'),'ImgList');

ImgNum = size(ImgList,1);
ListMeanValue = cell(ImgNum,1);

% add a wait bar
hwait = waitbar(0,'Please Wait ...');
step = ImgNum/100;

for ii = 1:ImgNum
    V = spm_vol(ImgList{ii});
    [pathstr,filename,~] = fileparts(ImgList{ii});
    [~,subfolder] = fileparts(pathstr);
    outputname = ['BA_Rep_',subfolder,'_',filename];
    
    if size(V,1) == 1 %3D Img
        SubMeanValue = BA_Report(V,threshold,clustersize,[outputname,'.txt']);
        ListMeanValue{ii,1} = SubMeanValue;
        
    else % 4D Img
        SetNum = size(V,1);
        SetMeanValue = cell(SetNum,1);
        for jj = 1:SetNum
            SigleImgValue = BA_Report(V(jj,1),threshold,clustersize,[outputname,int2str(jj),'.txt']);
            SetMeanValue{jj,1} = SigleImgValue;
        end
        ListMeanValue{ii,1} = SetMeanValue;
    end
    
    % wait bar
    PerStr = fix(ii/step);
    str = ['Running ...',num2str(PerStr),'%'];
    waitbar(ii/ImgNum,hwait,str);
    pause(0.05);
end

close(hwait);% close the wait bar

save AllMeanValue ListMeanValue

res_list = dir('BA_Rep_*.txt');
if ~isempty(res_list)
    movefile('BA_Rep_*.txt',fullfile(pwd,ResDir))% the res report txt file may like fine in matlab
end
movefile('AllMeanValue.mat',fullfile(pwd,ResDir))

clear;clc;
disp('All Work Done');

end

function Cluster_Mean_All = BA_Report(V,Threshold,ClusterSize,output_txt)
%% function to report an img brain area info
% ImgFilePath. e.g. 'd:/temp/test.nii'
% Thereshold,e.g.2.58
% ClusterSize,e.g.10
% output_txt,e.g. 'BA_Report.txt'

%     V = spm_vol(ImgFilePath);
    [Y,XYZ] = spm_read_vols(V);
    TransMat = V.mat;
    if Threshold == 0
        Index = find(Y~=0);
    elseif Threshold > 0
        Index = sort([find(Y >= Threshold);find(Y <= -Threshold)]);
    elseif Threshold < 0
        Index = find(Y <= Threshold);
    end
        
if ~isempty(Index)%when the voxel value supass the thresh
    
    MNI = XYZ(:,Index)'; % XYZ is a 3*n matrix,so MNI is n*3
    % how many clusters in one ROI
    Coord = mni2cor(MNI,TransMat);
    ClusterInfo = spm_clusters(Coord');
    
% ********************** delete small clusters ****************************
    % get the frequency table
    ClusterTable = tabulate(ClusterInfo);
    SmallClusterIdx = find(ClusterTable(:,2)<ClusterSize);
    if ~isempty(SmallClusterIdx)
        for kk = 1:length(SmallClusterIdx)
            if kk == 1
            Idx2Del = find(ClusterInfo == SmallClusterIdx(kk));
            else
                Idx2Del_temp = find(ClusterInfo == SmallClusterIdx(kk));
                Idx2Del = [Idx2Del Idx2Del_temp];
            end
        end
        % delete small cluster mni idx
        MNI(Idx2Del,:) = [];
        % delete small cluster id
        ClusterInfo(Idx2Del) = [];
        % modified 20141101
        Index(Idx2Del) = [];
    end
% *************************************************************************

    ClusterID = unique(ClusterInfo);
    Cluster_Mean_All = zeros(1,length(ClusterID));
    idx = 1;
    % check filterd cluster info again
    if ~isempty(ClusterInfo)
        for clusteridx = 1:length(ClusterID)
            pos = find(ClusterInfo == ClusterID(clusteridx));
            numVoxels = length(pos);
            MNI_cluster = MNI(pos,:);

            % comput the mean value in cluster
            cluster_idx = Index(pos);
            Y_cluster = Y(cluster_idx);
            Cluster_Mean = mean(Y_cluster);
            Cluster_Mean_All(idx) = Cluster_Mean;idx = idx +1;

            % get the peak coord in the cluster
            % 20141101
            % fit both for positive and negative
            if max(Y_cluster) > 0
                Cluster_Peak = max(Y_cluster);
            else
                Cluster_Peak = min(Y_cluster);
            end
            PeakIndex = Y_cluster == Cluster_Peak;
            Coord_Peak = MNI_cluster(PeakIndex,:);

            % 20141101 Coord_Peak may be n*3 format matrix
            % while length(Coord_Peak) == 3
            % so change to use function "size" to fix the bug
            if size(Coord_Peak,1)>1
                PeakCoord = round(mean(Coord_Peak));
            else
                PeakCoord =  Coord_Peak;
            end
            [Report1,Report2] = CoordReport(MNI_cluster);    % Report2 is 1*n mat

            Report{1,1} = [{'----------------------'};...
                {['Cluster ' num2str(clusteridx)]};...
                {['Number of voxels: ' num2str(numVoxels)]};...
                {['Cluster Mean Value: ' num2str(Cluster_Mean)]};...
                {['Cluster Peak Value: ' num2str(Cluster_Peak)]};...
                {['Cluster Peak Coord: ' num2str(PeakCoord)]};Report2'];

            cell2txt(Report,output_txt)

        end
    else % if the cluster is empty,
        msg = fprintf('%s is empty after threshold\n',V.fname);
        Cluster_Mean_All = [];
        disp(msg);
    end
else
    msg = sprintf('%s is empty after threshold\n',V.fname);
    Cluster_Mean_All = [];
    disp(msg);
end

end


%% cuixuFindStructure
function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
%% function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
%
% this function converts MNI coordinate to a description of brain structure
% in aal
%
%   mni: the coordinates (MNI) of some points, in mm.  It is Nx3 matrix
%   where each row is the coordinate for one point
%   DB (optional): the database.  This variable is available if you load
%   TDdatabase.mat (available to download at www.alivelearn.net/xjview. Or
%   you omit this parameter, then you need to have TDdatabase.mat in your
%   matlab path.
%
%   onelinestructure: description of the position, one line for each point
%   cellarraystructure: description of the position, a cell array for each point
%
%   Example:
%   cuixuFindStructure([72 -34 -2; 50 22 0], DB)
%
% Xu Cui
% 2007-11-20
%
if nargin==1
    load('TDdatabase.mat');
    TALMNI = 1;
elseif nargin == 2
    TALMNI = 1;
end


N = size(mni, 1);

% round the coordinates
mni = round(mni/2) * 2;

T = [...
     2     0     0   -92
     0     2     0  -128
     0     0     2   -74
     0     0     0     1];

index = mni2cor(mni, T); 

cellarraystructure = cell(N, length(DB));
onelinestructure = cell(N, 1);

for ii=1:N
    for jj=1:length(DB)
        graylevel = DB{jj}.mnilist(index(ii, 1), index(ii, 2),index(ii, 3));
        if graylevel == 0
            thelabel = 'undefined';
        else
            if jj==length(DB); tmp = ' (aal)'; else tmp = ''; end
            thelabel = [DB{jj}.anatomy{graylevel} tmp];
        end
        cellarraystructure{ii, jj} = thelabel;
        onelinestructure{ii} = [ onelinestructure{ii} ' // ' thelabel ];
    end
end
end

%% mni2cor
function coordinate = mni2cor(mni, T)
% function coordinate = mni2cor(mni, T)
% convert mni coordinate to matrix coordinate
%
% mni: a Nx3 matrix of mni coordinate
% T: (optional) transform matrix
% coordinate is the returned coordinate in matrix
%
% caution: if T is not specified, we use:
% T = ...
%     [-4     0     0    84;...
%      0     4     0  -116;...
%      0     0     4   -56;...
%      0     0     0     1];
%
% xu cui
% 2004-8-18
%

if isempty(mni)
    coordinate = [];
    return;
end

if nargin == 1
	T = ...
        [-4     0     0    84;...
         0     4     0  -116;...
         0     0     4   -56;...
         0     0     0     1];
end

coordinate = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(T))';
coordinate(:,4) = [];
coordinate = round(coordinate);
return;
end

%% cell2txt
function cell2txt(Cell2Trans,TxtName)
% Cell2Trans is n*1 format cell,Cell2Trans{1,1} is a cell too
% TxtName e.g. 'output.txt'

% codes below from
% http://stackoverflow.com/questions/8565617/print-a-cell-array-as-txt-in-matlab
% modified by hscheng
    fid = fopen(TxtName, 'a');
    for ii = 1:size(Cell2Trans,1)
        Temp = Cell2Trans{ii,1};
        fstr = '';
        
        % generate a print format
        for jj=1:size(Temp,2)
           
           switch class(Temp{jj})
               case 'char'
                   fstr = [fstr '%s'];
               otherwise
                   % Assume numeric
                   fstr = [fstr '%g'];
           end
           
           if jj < size(Temp,2)
               fstr = [fstr '\t'];
           else
               fstr = [fstr '\n'];
           end
        
        end
        
        C = Temp.';

        fprintf(fid, fstr,C{:});
    end

    fclose(fid);
end

%% CoordReport
function [Report,Report2] = CoordReport(XYZ)


    %*************************************************************************%
    %    codes below from xjview.m(offered by Xu Cui-www.alivelearn.net)
    %*************************************************************************%

    % list structure of voxels in this cluster
    [a, b] = cuixuFindStructure(XYZ);
    names = unique(b(:));
    index = NaN*zeros(length(b(:)),1);
    for ii=1:length(names)
        pos = strcmp(b(:),names{ii});
        index(pos) = ii;
    end

    Report = {};

    for ii=1:max(index)
        Report{ii,1} = names{ii};
        Report{ii,2} = length(find(index==ii));
    end
    for ii=1:size(Report,1)
        for jj=ii+1:size(Report,1)
            if Report{ii,2} < Report{jj,2}
                tmp = Report(ii,:);
                Report(ii,:) = Report(jj,:);
                Report(jj,:) = tmp;
            end
        end
    end
    Report = [{'structure','# voxels'}; {'--TOTAL # VOXELS--', length(a)}; Report];

    Report2 = {sprintf('%s\t%s',Report{1,2}, Report{1,1}),''};
    for ii=2:size(Report,1)
        if strcmp('undefined', Report{ii,1}); continue; end
        Report2 = [Report2, {sprintf('%5d\t%s',Report{ii,2}, Report{ii,1})}];
    end
end