% Last modified: 

clear all;
close all;

%% Changeable parameters
% base folder directory. Contains folder Data, Code
% Base_dir = 'D:\ResGenQSM';
Base_dir = 'D:\LIST\Research\2_Resolution_Generalized_NN_QSM\Github_upload';

% Test name that contains the data folders. within the Testname folder, individual subject data must 
% be named test_subject1, test_subject2, ...
% All data in same test name should have same resolution
Testname = 'Dummy_Test';

% current test data directory
data_dir = [Base_dir filesep 'Data\1_InputData\' Testname filesep];

% save directory for the re-sampled data for network infrerence
save_dir = [Base_dir filesep 'Data\2_DataForInference\' ];

% matfilename: name of ,mat file containing 
% voxel_size_input: input data resolution [mm^3]
% local_f_ppm: local field map [ppm]
% sus_label: label susceptibility map
% Mask: brain mask
% csf_mask: csf mask
matfilename = 'data'; 

% voxel size used for training the network
voxel_size_train = [1 1 1]*1.0; 

% number of orientations of data.
numorient = 1; 

% Number of re-sampling per dimension. Current version use the same
% Nsamp for all three dimensions.
Nsamp = 3;

%%
testfolderlist = dir([data_dir 'test_subject*']);
if(~exist('test_sub_list','var'))
    test_sub_list = 1:length(testfolderlist); % if not defined, generate for all subjects
end


Nsamp = [Nsamp Nsamp Nsamp];
mkdir(save_dir);
mkdir([save_dir filesep Testname]);
mkdir([save_dir filesep Testname filesep 'Proposed']);
mkdir([save_dir filesep Testname filesep 'Proposed' filesep 'Input']);
mkdir([save_dir filesep Testname filesep 'Proposed' filesep 'Label']);
mkdir([save_dir filesep Testname filesep 'Interpolation']);
mkdir([save_dir filesep Testname filesep 'Interpolation' filesep 'Input']);
mkdir([save_dir filesep Testname filesep 'Interpolation' filesep 'Label']);
mkdir([save_dir filesep Testname filesep 'Naive']);
mkdir([save_dir filesep Testname filesep 'Naive' filesep 'Input']);
mkdir([save_dir filesep Testname filesep 'Naive' filesep 'Label']);


%%
for sub=test_sub_list
    for orient = 1:numorient
        clear chi_cosmos        
        clear phs_tissue chi_cosmos;
        
        load([testfolderlist(sub).folder filesep testfolderlist(sub).name filesep matfilename '.mat']);
        eval(['local_f_infer = local_f_ppm(:,:,:,' num2str(orient) ');']);
        eval(['sus_infer = sus_label(:,:,:,' num2str(orient) ');']);
        eval(['mask_infer = Mask(:,:,:,' num2str(orient) ');']);
        eval(['csf_mask_infer = csf_mask(:,:,:,' num2str(orient) ');']);
        sus_gt = sus_label(:,:,:,1:numorient);
        mask = Mask(:,:,:,1:numorient);
        matrix_size_infer = size(local_f_infer);
        Mask_CSF = csf_mask(:,:,:,1:numorient);
        
        matrix_size_train = round(size(local_f_infer).*voxel_size_input./voxel_size_train/2)*2;
        voxel_size_train = size(local_f_infer)./matrix_size_train.*voxel_size_input;
        disp(['actual voxel size is ' num2str(voxel_size_train)]);
        
        % Naive
        
        phs_tissue(:,:,:,orient) = single(local_f_infer);
        chi_cosmos(:,:,:,orient) = single(sus_infer);
        
        if(orient == numorient)
            save([save_dir filesep Testname filesep  'Naive/Input/test_input_sub' num2str(sub) '_1_1_1.mat'], 'phs_tissue');
            save([save_dir filesep Testname filesep 'Naive/Label/test_label_sub' num2str(sub) '_1_1_1.mat'], 'chi_cosmos');
        end
        
        clear phs_tissue chi_cosmos;
    end
    for orient = 1:numorient
        % Interpolation
        clear chi_cosmos
        load([testfolderlist(sub).folder filesep testfolderlist(sub).name filesep matfilename '.mat']);
        eval(['local_f_infer = local_f_ppm(:,:,:,' num2str(orient) ');']);
        eval(['sus_infer = sus_label(:,:,:,' num2str(orient) ');']);
        eval(['mask_infer = Mask(:,:,:,' num2str(orient) ');']);
 
        
        undersampdim = single(matrix_size_infer>matrix_size_train);
        upsampledim = 1-undersampdim;
        matrix_size_tmp = matrix_size_train.*undersampdim+matrix_size_infer.*upsampledim;
        
        k_chi = fft3c(sus_infer);
        k_chi_to = k_chi(end/2-matrix_size_tmp(1)/2+1:end/2+matrix_size_tmp(1)/2,end/2-matrix_size_tmp(2)/2+1:end/2+matrix_size_tmp(2)/2,end/2-matrix_size_tmp(3)/2+1:end/2+matrix_size_tmp(3)/2)/(matrix_size_infer(1)/matrix_size_tmp(1)*matrix_size_infer(2)/matrix_size_tmp(2)*matrix_size_infer(3)/matrix_size_tmp(3));
        sus_tmp = real(ifft3c(k_chi_to));
        
        k_lf = fft3c(local_f_infer);
        k_localf_to = k_lf(end/2-matrix_size_tmp(1)/2+1:end/2+matrix_size_tmp(1)/2,end/2-matrix_size_tmp(2)/2+1:end/2+matrix_size_tmp(2)/2,end/2-matrix_size_tmp(3)/2+1:end/2+matrix_size_tmp(3)/2)/(matrix_size_infer(1)/matrix_size_tmp(1)*matrix_size_infer(2)/matrix_size_tmp(2)*matrix_size_infer(3)/matrix_size_tmp(3));
        local_f_tmp = real(ifft3c(k_localf_to));
        
        matrix_size_diff = (matrix_size_train-matrix_size_tmp)/2;
        scale = matrix_size_train./matrix_size_tmp;
        scalefac = scale(1)*scale(2)*scale(3);
        k_tmp = zeros(matrix_size_train);
        k_tmp(matrix_size_diff(1)+1:end-matrix_size_diff(1),matrix_size_diff(2)+1:end-matrix_size_diff(2),matrix_size_diff(3)+1:end-matrix_size_diff(3)) = fft3c(sus_tmp)*scalefac;
        sus_undersamp = real(ifft3c(k_tmp));
        
        
        k_tmp = zeros(matrix_size_train);
        k_tmp(matrix_size_diff(1)+1:end-matrix_size_diff(1),matrix_size_diff(2)+1:end-matrix_size_diff(2),matrix_size_diff(3)+1:end-matrix_size_diff(3)) = fft3c(local_f_tmp)*scalefac;
        local_f_undersamp = real(ifft3c(k_tmp));
        
        chi_cosmos_(:,:,:,orient) = single(sus_undersamp);
        phs_tissue_(:,:,:,orient) = single(local_f_undersamp);
        
        if(orient==numorient)
            clear chi_cosmos phs_tissue;
            chi_cosmos = chi_cosmos_;
            phs_tissue = phs_tissue_;
            save([save_dir filesep Testname filesep  'Interpolation/Input/test_input_sub' num2str(sub) '_1_1_1.mat'], 'phs_tissue');
            save([save_dir filesep Testname filesep 'Interpolation/Label/test_label_sub' num2str(sub) '_1_1_1.mat'], 'chi_cosmos');
        end
        clear phs_tissue chi_cosmos;
    end
    
    %% Proposed
    shift_tmpx = [0:1/Nsamp(1):1];
    shift_tmpy = [0:1/Nsamp(2):1];
    shift_tmpz = [0:1/Nsamp(3):1];
    shiftvalx = shift_tmpx(1:end-1).*voxel_size_train(1)./voxel_size_input(1);
    shiftvaly = shift_tmpy(1:end-1).*voxel_size_train(2)./voxel_size_input(2);
    shiftvalz = shift_tmpz(1:end-1).*voxel_size_train(3)./voxel_size_input(3);
    
    
    for shiftx=1:length(shiftvalx)
        for shifty=1:length(shiftvaly)
            for shiftz=1:length(shiftvalz)
                for orient = 1:numorient
                    clear chi_cosmos
                    load([testfolderlist(sub).folder filesep testfolderlist(sub).name filesep matfilename '.mat']);
                    eval(['local_f_infer = local_f_ppm(:,:,:,' num2str(orient) ');']);
                    eval(['sus_infer = sus_label(:,:,:,' num2str(orient) ');']);
                    eval(['mask_infer = Mask(:,:,:,' num2str(orient) ');']);
                    
                    
                    gaussfac=1.5;
                    sus_tmp = shift_usingkspace(sus_infer,[shiftvalx(shiftx),shiftvaly(shifty),shiftvalz(shiftz)]);
                    local_f_tmp = shift_usingkspace(local_f_infer,[shiftvalx(shiftx),shiftvaly(shifty),shiftvalz(shiftz)]);
                    
                    chi_cosmos_(:,:,:,orient) = single(interp3_gausssinc(sus_tmp,matrix_size_infer,matrix_size_train,voxel_size_input,voxel_size_train,gaussfac));
                    phs_tissue_(:,:,:,orient) = single(interp3_gausssinc(local_f_tmp,matrix_size_infer,matrix_size_train,voxel_size_input,voxel_size_train,gaussfac));
                    

                    if(orient == numorient)
                        clear chi_cosmos phs_tissue;
                        chi_cosmos = chi_cosmos_;
                        phs_tissue = phs_tissue_;
                        save([save_dir filesep Testname filesep  'Proposed/Input/test_input_sub' num2str(sub) '_' num2str(shiftx) '_' num2str(shifty) '_'  num2str(shiftz) '.mat'], 'phs_tissue');
                        save([save_dir filesep Testname filesep 'Proposed/Label/test_label_sub' num2str(sub) '_' num2str(shiftx) '_' num2str(shifty) '_'  num2str(shiftz) '.mat'], 'chi_cosmos');
                    end

                end
            end
        end
        save([save_dir filesep Testname filesep '/metadata.mat'], 'voxel_size_train',  'voxel_size_input','shiftvalx', 'shiftvaly', 'shiftvalz');

        save([save_dir filesep Testname filesep '/sub' num2str(sub) '.mat'], 'sus_gt', 'mask', 'Mask_CSF');
        % figure;imshow3Dfull(local_f_cat,[-.1,.1]);
    end
end