%% Changeable parameters
Base_dir = 'D:\LIST\Research\2_Resolution_Generalized_NN_QSM\Github_upload';

voxel_size_train = [1 1 1]*1.0;
numtestsub =1;
numorient = 1;
testname = 'Dummy_Test';
test_type = 'Proposed';

NetworkOutputFolder = [Base_dir filesep 'Data\\3_NetworkOutput\' testname filesep test_type filesep];
FinalOutputFolder = [Base_dir filesep 'Data\\4_FinalOutput\' testname filesep test_type filesep];
sus_var_name = 'sus';
%%
if(~strcmp(test_type,'Proposed'))
    shiftnum = 1;
end

testdatafolder = [testname '\' test_type];
load([Base_dir filesep 'Data\2_DataForInference\' testname  '\metadata.mat']);
outputmats = dir([NetworkOutputFolder '*.mat']);
load([outputmats(1).folder filesep outputmats(1).name]);

eval(['sus=' sus_var_name ';']);

voxel_size_from = voxel_size_train;
voxel_size_to = voxel_size_input;
matrix_size_from = size(sus,[1,2,3]);
matrix_size_to = round(matrix_size_from.*voxel_size_from./voxel_size_to/2)*2;
clear tmp;

%%
if(strcmp(test_type,'Proposed'))
    lam = 0.2;
    B0_dir = [0 0 1];
    D2  = fftshift(dipole_kernel(matrix_size_from, voxel_size_from, B0_dir));
    D1 =  fftshift(dipole_kernel(matrix_size_to, voxel_size_to, B0_dir));
    
    D1(0<(D1)&D1<lam) = lam;
    D1(-lam<(D1)&D1<=0) = -lam;

    Dtmp = repmat(D2,[3,3,3]);
    D2new = (Dtmp(end/2-matrix_size_to(1)/2+1:end/2+matrix_size_to(1)/2,end/2-matrix_size_to(2)/2+1:end/2+matrix_size_to(2)/2,end/2-matrix_size_to(3)/2+1:end/2+matrix_size_to(3)/2));

    scale = ((D2new./D1));
    scale(isnan(scale)) = 1;
    scale(isinf(scale)) = 1;
    
    matrix_size_ = min(matrix_size_from,size(scale));
    scale(end/2-matrix_size_(1)/2+1:end/2+matrix_size_(1)/2,end/2-matrix_size_(2)/2+1:end/2+matrix_size_(2)/2,end/2-matrix_size_(3)/2+1:end/2+matrix_size_(3)/2)=1;
%     scale = ones(size(scale));
    for sub=1:numtestsub
        load([Base_dir filesep 'Data\2_DataForInference\' testname  '\metadata.mat']);
        load([Base_dir filesep 'Data\2_DataForInference\' testname  '\sub' num2str(sub) '.mat']);

        sus_orig = sus_gt;
        sus_recon = zeros([matrix_size_to, numorient]);
        for orient = 1:numorient
            for shiftx=1:length(shiftvalx)
                for shifty=1:length(shiftvaly)
                    for shiftz=1:length(shiftvalz)
                        filename = dir([NetworkOutputFolder 'sub' num2str(sub) '_shift' num2str(shiftx) '_' num2str(shifty) '_' num2str(shiftz) '_*.mat']);
                        load([filename.folder filesep filename.name]);
%                         sus_recon(:,:,:,orient) = sus_recon(:,:,:,orient) + real(shift_usingkspace(ifft3c(scale.*fft3c(recover_resolution(sus(:,:,:,orient),matrix_size_to))),-1*[shiftvalx(shiftx),shiftvaly(shifty),shiftvalz(shiftz)]))/(length(shiftvalx)*length(shiftvaly)*length(shiftvalz));
                          sus_recon(:,:,:,orient) = sus_recon(:,:,:,orient) + real(shift_usingkspace(ifft3c(scale.*fft3c(interp3_gausssinc(sus,matrix_size_from,matrix_size_to,voxel_size_from,voxel_size_to,5))),-1*[shiftvalx(shiftx),shiftvaly(shifty),shiftvalz(shiftz)]))/(length(shiftvalx)*length(shiftvaly)*length(shiftvalz));
                    end
                end
            end
        end
        save([FinalOutputFolder 'recon_sub' num2str(sub) '_dipolecompensated_lambda' strrep(sprintf('%1.1f',lam),'.','_') '.mat'], 'sus_recon', 'sus_orig','mask');
%         savenii(sus_recon,voxel_size_infer,[DLoutputfolder 'recon_sub' num2str(sub) '_dipolecompensated_lambda' strrep(sprintf('%1.1f',lam),'.','_')]);
    end
else
    for sub=1:numtestsub
        load([NetworkOutputFolder testname '\metadata.mat']);
        load([NetworkOutputFolder testname '\sub' num2str(sub) '.mat']);
        filename = dir([NetworkOutputFolder 'sub' num2str(sub) '_shift' num2str(1) '_' num2str(1) '_' num2str(1) '_*.mat']);
        load([filename.folder filesep filename.name]);
        sus_orig = sus_gt;
        sus_recon = sus_gt;
        if(size(sus_gt,1)~=size(sus,1) || size(sus_gt,2)~=size(sus,2) || size(sus_gt,3)~=size(sus,3))
            matrix_size_train = size(sus);
            matrix_size_infer = size(sus_gt);
            [Xq,Yq,Zq] = meshgrid([-matrix_size_infer(2)/2:matrix_size_infer(2)/2-1]*voxel_size_infer(2),[-matrix_size_infer(1)/2:matrix_size_infer(1)/2-1]*voxel_size_infer(1),[-matrix_size_infer(3)/2:matrix_size_infer(3)/2-1]*voxel_size_infer(3));
            [X,Y,Z] = meshgrid([-matrix_size_train(2)/2:matrix_size_train(2)/2-1]*voxel_size_train(2),[-matrix_size_train(1)/2:matrix_size_train(1)/2-1]*voxel_size_train(1),[-matrix_size_train(3)/2:matrix_size_train(3)/2-1]*voxel_size_train(3));
            for orient=1:numorient
                sus_recon(:,:,:,orient) = interp3(X,Y,Z,sus(:,:,:,orient),Xq,Yq,Zq,'spline');
            end
        else
            sus_recon = sus;
        end
        save([FinalOutputFolder 'recon_sub' num2str(sub) '.mat'], 'sus_recon', 'sus_orig','mask');
%         savenii(sus_recon,voxel_size_infer,[DLoutputfolder 'recon_sub' num2str(sub) '']);
    end
end