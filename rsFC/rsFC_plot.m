roi_labels = {'ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','ROI7','ROI8','ROI9','ROI10','ROI11','ROI12'};

%% threshold p value
p_thresh = 0.2;

raw_r_3d = rsFC_Res.beh_fc_r_3d;
raw_p_3d = rsFC_Res.beh_fc_p_3d;

p_thresh_idx = 0<raw_p_3d & raw_p_3d<p_thresh;
r_threshed = raw_r_3d.*p_thresh_idx;
p_threshed = raw_p_3d.*p_thresh_idx;

%% plot thresholded r p mat
mkdir('beh_fc_img');

for ii = 1:length(beh_label)
    corr_mat = r_threshed(:,:,ii);
    matrixplot(corr_mat,'XVarNames',roi_labels,'YVarNames',roi_labels,'TextColor',[0.6,0.6,0.6],'ColorBar','on');
    title(beh_label{ii});
    % xlabel(beh_label{ii});
    saveas(gcf,['beh_fc_img\',beh_label{ii},'.jpg'])
    close
end

% r_beh = p_threshed(:,:,16);%r_threshed(:,:,4).*r_threshed(:,:,16);
% %     p_threshed(:,:,7);
% matrixplot(r_beh,'XVarNames',roi_labels,'YVarNames',roi_labels,'TextColor',[0.6,0.6,0.6],'ColorBar','on');

%% 2 behavor fc correlation overlap
idx1 = 4; idx2 = [7 8 16];

for jj = 1:length(idx2)
    r_beh = r_threshed(:,:,idx1).*r_threshed(:,:,idx2(jj));
    matrixplot(r_beh,'XVarNames',roi_labels,'YVarNames',roi_labels,'TextColor',[0.6,0.6,0.6],'ColorBar','on');
    res_mark = ['roi_',int2str(idx1),'_',int2str(idx2(jj))];
    title(res_mark);
    saveas(gcf,[res_mark,'.jpg'])
    close
end