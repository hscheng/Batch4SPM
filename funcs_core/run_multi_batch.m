function run_multi_batch()

batch_mat_list = ScanImg2Cell('Select batch mat files','mat');

for ii = 1:length(batch_mat_list)
    matlabbatch = importdata(batch_mat_list{ii});
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end