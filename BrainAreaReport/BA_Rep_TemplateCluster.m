function BA_Rep_TemplateCluster()
% written by hscheng [hongshengcheng.math@gmail.com]
% 
% 20140820 add the cluster num
% 20140627 add time bar and modified the cls xlsx info check
% modified date 20140624
% create date 20140421
%
% function to merge a cluster from template
% and report the brain area in each cluster
% OR just report the roi regins in the template


    %choose the template
    [TemplateName,TemplatePath] = uigetfile({'*.img;*.nii'},'Please Select your Template');
    TemplateLoc = fullfile(TemplatePath,TemplateName);

    %read template
    Vtem = spm_vol(TemplateLoc);
    % choose
    choice = questdlg('Which work ?','Quest','ROI','Cluster','Cluster');

    switch choice
        case 'ROI'
            % make a folder to store results
            if exist('ROI_BA_Reports','dir')~=7
                mkdir ROI_BA_Reports
            end
            ROIndex = input('ROI Index (e.g.1:1024) : ');
            ROI_Report(Vtem,ROIndex);
        case 'Cluster'
            % make a folder to store results
            if exist('Cluster_BA_Reports','dir')~=7
                mkdir Cluster_BA_Reports
            end
            CluterReport(Vtem)
    end
end

function ROI_Report(Vtem,ROIndex)
    % read the img info
    [Ytem, XYZtem] = spm_read_vols(Vtem);
    Ytem(isnan(Ytem)) = 0;
    Ytem = round(Ytem);
    %set the roi index
    
    ROINum=length(ROIndex);

    ROI_Reports = cell(ROINum,2);
    ROI_Center = zeros(ROINum,3);
    ROI_ClusterNum = zeros(ROINum,1);
    TransMat = Vtem.mat;


    for nn = 1:ROINum
        Region = ROIndex(nn);
        Index = Ytem==Region;

        XYZ = XYZtem(:,Index)'; % XYZ is a n*3 matrix
        % how many clusters in one ROI
        Coord = mni2cor(XYZ,TransMat);
        ClusterInfo = spm_clusters(Coord');

        ClusterID = unique(ClusterInfo);
        ROI_ClusterNum(nn) = length(ClusterID);
        % report the brain area info
        [roi_report,roi_report2] = CoordReport(XYZ);

        ROI_Reports{nn,1} = roi_report;
        ROI_Reports{nn,2} = roi_report2;
        
        ROI_Center(nn,:) = round(mean(XYZ));

    end
    
    save ROI_Reports ROI_Reports ROI_Center ROI_ClusterNum
    
    %save the ROI information to txt format
    txtTrans = ROI_Reports(:,2);
    txtname = 'ROI_Reports';
    cell2txt(txtTrans,txtname)
    
    movefile('ROI_Reports*.mat','ROI_BA_Reports\')
    movefile('ROI_Reports*.txt','ROI_BA_Reports\')
    
    disp('All work done,Please check your result files')
%*************************************************************************%
end

function CluterReport(Vtem)

    [Ytem, XYZtem] = spm_read_vols(Vtem);
    Ytem(isnan(Ytem)) = 0;

    Vmask = Vtem; % copy input info for output image
    %select the cluster xls file
    [ClusterName,ClusterPath] = uigetfile({'*.xls;*.xlsx;*.txt'},'Please choose your Cluster File');
    [~,~,ext] = fileparts(ClusterName);
    if strcmp(ext,'.txt')
        % the cluster.txt file only contain numeric values
        ClusterMat = importdata(fullfile(ClusterPath,ClusterName));
    else
        [ClusterMat,ClusterNames,raw] = xlsread(fullfile(ClusterPath,ClusterName),1);
    end
    cls_name_mark = isempty(ClusterNames);
    
    ClusterNum = size(ClusterMat,1);

    ClusterReports = cell(ClusterNum,2);
    % add a wait bar
    % http://blog.sina.com.cn/s/blog_4d633dc70100nwcw.html
    hwait = waitbar(0,'Please Wait ...');
    step = ClusterNum/100;
    for jj = 1:ClusterNum

        ROIIndex = ClusterMat(jj,:);
        ROIIndex = ROIIndex(isfinite(ROIIndex));%filter the NAs
        ROINum = length(ROIIndex); % how many ROIs in this cluster

        for kk = 1:ROINum
            TempIndex = Ytem==ROIIndex(kk);
            if kk==1
                Ymask = TempIndex;
            else
                Ymask = Ymask+TempIndex;
            end

        end
        if cls_name_mark == 1
            Vmask.fname = fullfile(pwd, ['Cluster_' num2str(jj) '.nii']);
        else
%             Vmask.fname = fullfile(pwd, ['Cluster_',ClusterNames{jj},'.nii']);
            Vmask.fname = fullfile(pwd, [ClusterNames{jj},'.nii']);
        end
        spm_write_vol(Vmask,Ymask);

        Index = Ymask==1;
        %Index = find(Ymask==1);
        XYZ = XYZtem(:,Index)';

        %get the Region Report
        [cluster_report,cluster_report2] = CoordReport(XYZ);
        ClusterReports{jj,1} = cluster_report;
        ClusterReports{jj,2} = cluster_report2;
    %     save(['Cluster_' num2str(jj) '_Report'],'ROIReport','ROIReport2');
    
    % wait bar
    PerStr = fix(jj/step);
    str = ['Running ...',num2str(PerStr),'%'];
    waitbar(jj/ClusterNum,hwait,str);
    pause(0.05);

    end
    close(hwait);
    
    save ClusterReports ClusterReports
    
    %save the Clusters information to txt format
%     txtTrans=ClusterReports(:,2);
    txtTrans = strcat(ClusterNames,ClusterReports(:,2));
    txtname = 'ClusterReports';
    cell2txt(txtTrans,txtname)
    
    movefile('Cluster*.nii','Cluster_BA_Reports\')
    movefile('Cluster*.mat','Cluster_BA_Reports\')
    movefile('Cluster*.txt','Cluster_BA_Reports\')
%     movefile('Cluster*.*','Results\')
    disp('All work done,Please check your result files')
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
%% cell2txt
function cell2txt(Cell2Trans,TxtName)
%     TxtName=input('Please Enter the Result TXT File Name (e.g. res) : ','s');

% codes below from
% http://stackoverflow.com/questions/8565617/print-a-cell-array-as-txt-in-matlab
% modified by hscheng
    fid = fopen([TxtName,'.txt'], 'a');
    for ii = 1:size(Cell2Trans,1)
        Temp=Cell2Trans{ii,1};
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
%% cuixuFindStructure
function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
% function [onelinestructure, cellarraystructure] = cuixuFindStructure(mni, DB)
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

