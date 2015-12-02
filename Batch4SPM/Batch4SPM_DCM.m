function Batch4SPM_DCM()
%% function to run dynamic causal modeling
% written by hongshengcheng.math@gmail.com
% created date:20141201
%
% SPECIFICATION DCM 
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

fprintf('Welcome to our lab:  <a href="http://www.qiujlab.com/">Qiu Jiang Lab</a>\nhttp://www.qiujlab.com/\n');
fprintf('<a href="http://www.creativitybrain.com/">Center for Creativity and Brain</a>\nhttp://www.creativitybrain.com/\n');
disp('Any question, Please contact me at hongshengcheng.math@gmail.com');
disp('****************************************************************');

    spm_mat_list = ScanImg2Cell('Select SPM Mat Files','mat');
    
    para = inputdlg('voi names','input',[10 30],{''});
    voi_names_str = para{1};
    voi_num = size(voi_names_str,1);
    voi_names = cell(voi_num,1);
    for kk = 1:voi_num
        str = voi_names_str(kk,:);
        str(isspace(str)) = [];
        voi_names{kk} = str;
    end
    
    dcm_temp_list = ScanImg2Cell('select DCM model template files','mat');
    model_num = length(dcm_temp_list);
    
    for ii = 1:length(spm_mat_list)
        
        sub_spm_mat = spm_mat_list{ii};
        sub_dir = fileparts(sub_spm_mat);
        
        if ispc
            voi_mat_list = strcat(sub_dir,'\',voi_names);
        else %suit for linux
            voi_mat_list = strcat(sub_dir,'/',voi_names);
        end
        
        for jj = 1:model_num
            DCM_template = importdata(dcm_temp_list{jj});
            [sub_dir,model_name] = fileparts(dcm_temp_list{jj});
            DCM_mod_est = Sub_DCM(sub_spm_mat,voi_mat_list,DCM_template,model_name);
            save([sub_dir,DCM_mod_est],'DCM_mod_est');
        end
    end
    
end

function DCM_mod_est = Sub_DCM(spm_mat_path,voi_mat_list,DCM_template,model_name)

clear DCM

% SPECIFICATION DCM 
%--------------------------------------------------------------------------
% To specify a DCM, you might want to create a template one using the GUI
% then use spm_dcm_U.m and spm_dcm_voi.m to insert new inputs and new
% regions. The following code creates a DCM file from scratch, which
% involves some technical subtleties and a deeper knowledge of the DCM
% structure.

SPM = importdata(spm_mat_path);
sub_dir = fileparts(spm_mat_path);

voi_num = length(voi_mat_list);

for ii = 1:voi_num
    load(voi_mat_list{ii},'xY');
    DCM.xY(ii) = xY;
    DCM.Y.y(:,ii)  = xY.u;
    DCM.Y.name{ii} = xY.name;
end
% load(fullfile(data_path,'GLM','VOI_V1_1.mat'),'xY');
% DCM.xY(1) = xY;
% load(fullfile(data_path,'GLM','VOI_V5_1.mat'),'xY');
% DCM.xY(2) = xY;
% load(fullfile(data_path,'GLM','VOI_SPC_1.mat'),'xY');
% DCM.xY(3) = xY;

DCM.n = voi_num; % length(DCM.xY),number of regions
DCM.v = length(DCM.xY(1).u); % number of time points

DCM.Y.dt  = SPM.xY.RT;
DCM.Y.X0  = DCM.xY(1).X0;

% for i = 1:voi_num % DCM.n
%     DCM.Y.y(:,i)  = DCM.xY(i).u;
%     DCM.Y.name{i} = DCM.xY(i).name;
% end

DCM.Y.Q    = spm_Ce(ones(1,voi_num)*DCM.v);

DCM.U.dt   =  SPM.Sess.U(1).dt;
DCM.U.name = [SPM.Sess.U.name];
DCM.U.u    = [SPM.Sess.U(1).u(33:end,1) ...
              SPM.Sess.U(2).u(33:end,1) ...
              SPM.Sess.U(3).u(33:end,1)];

DCM.delays = repmat(SPM.xY.RT,voi_num,1);% voi_num*1 format mat
DCM.TE     = 0.04; % echo time,TE

DCM.options.nonlinear  = 0; % bilinear or nonlinear
DCM.options.two_state  = 0; % state per region:one/two
DCM.options.stochastic = 0; % stochastic effects : no/yes
DCM.options.nograph    = 1;

DCM.a = DCM_template.a;
DCM.b = DCM_template.b;
DCM.c = DCM_template.c;
DCM.d = DCM_template.d;

dcm_mod_mat = fullfile(sub_dir,[model_name,'.mat']);
save(dcm_mod_mat,'DCM'); % 

DCM_mod_est = spm_dcm_estimate(dcm_mod_mat);

% specific intrinsic connections,mat size (voi_num,voi_num)
% from col_idx to row_idx
% from voi_2 to voi_1 [0 1 0;0 0 0;0 0 0]
DCM.a = [1 1 0;
         1 1 1;
         0 1 1];

% SPECIFICATION DCM "attentional modulation of backward connection"
% effects of conditions on regions and connections
% DCM.b = zeros(voi_num,voi_num,cond_num)
DCM.b = zeros(3,3,3);  DCM.b(2,1,2) = 1;  DCM.b(2,3,3) = 1;

% voi idx mark
DCM.c = [1 0 0; 0 0 0; 0 0 0];
DCM.d = zeros(3,3,0);% only for nonlinear,otherwize should be empty

save(fullfile(sub_dir,'DCM_mod_bwd.mat'),'DCM'); % 

% SPECIFICATION DCM "attentional modulation of forward connection"
%--------------------------------------------------------------------------
DCM.b = zeros(3,3,3);  DCM.b(2,1,2) = 1;  DCM.b(2,1,3) = 1;

save(fullfile(sub_dir,'DCM_mod_fwd.mat'),'DCM');

% ESTIMATION
%--------------------------------------------------------------------------
DCM_bwd = spm_dcm_estimate(fullfile(sub_dir,'DCM_mod_bwd.mat'));
DCM_fwd = spm_dcm_estimate(fullfile(sub_dir,'DCM_mod_fwd.mat'));

% BAYESIAN MODEL COMPARISON
%--------------------------------------------------------------------------
fprintf('Model evidence: %f (bwd) vs %f (fwd)\n',DCM_bwd.F,DCM_fwd.F);
end