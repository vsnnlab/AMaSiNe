clear all; close all;
warning('off')
STEP_0_Parameters;
img_name=Img_filename_list;
load('Step1_Outline_result.mat')

no_images=size(anc_img_IDs, 2);
starting_ang_resolution=8;
denominator_step=2;
yaw_stage1=0; pitch_stage1=0;

BlockSz_Stg1=[13 13];
BlockSz_Stg2=[13 13];
BlockSz_Stg3=[13 13];
BlockSz_Stg4=[13 13];
BlockSz_Stg5=[13 13];
BlockSz_Stg6=[13 13];
BlockSz_Stg7=[13 13];

%% LOAD ALLEN ATLAS 

if strcmp(Structure_stain,'DAPI') || strcmp(Structure_stain,'Nissl')
    [VOL, metaVOL] = nrrdread('ara_nissl_25_2017.nrrd');
    VOL=uint8(uint16(rot90(permute(VOL,[3 1 2]),3))/(2^8));
elseif strcmp(Structure_stain,'AutoF')
    load('AutofluoresenceAtlas.mat');
    VOL=rot90(permute(VOL,[3 1 2]),3);
else
    error('Staining method not recognized');
end

ref_resc=1/0.5;
VOL=imresize3(VOL,1/ref_resc);
downscaled_xy_pix=ref_atlas_vox_res*ref_resc;


%% STAGE 1 : Contrast Adjustment and AP Positioning 
computational_time_checker_stage1=tic;
disp('STAGE 1 : Contrast Adjustment and AP Positioning')

tform_yaw_1=transform_matrix_0822( yaw_stage1,[0 1 0]);
tform_pitch_1=transform_matrix_0822(pitch_stage1,[1 0 0] );
tform_combined_1=mtimes(tform_yaw_1, tform_pitch_1);
tf = affine3d(tform_combined_1);
VOL_rot=imwarp(VOL,tf,'cubic');

img_library(1:no_images)=struct('img',[]);
for img_no_ii=1:no_images
    
    %%%% load raw image %%%%
    img_act=(imread(img_name{anc_img_IDs(img_no_ii),Color_Channel_Structure}));
    try
        img_act=rgb2gray(img_act);
    end
    
    %%%% Leave the slice part only %%%%
    img_act_pad=zeros(size(img_act));
    img_act=img_act(img_info(anc_img_IDs(img_no_ii)).slice_window(1):img_info(anc_img_IDs(img_no_ii)).slice_window(2),...
        img_info(anc_img_IDs(img_no_ii)).slice_window(3):img_info(anc_img_IDs(img_no_ii)).slice_window(4));
    
    img_act_pad=img_act_pad(img_info(anc_img_IDs(img_no_ii)).slice_window(1):img_info(anc_img_IDs(img_no_ii)).slice_window(2),...
        img_info(anc_img_IDs(img_no_ii)).slice_window(3):img_info(anc_img_IDs(img_no_ii)).slice_window(4));
    img_act_pad(img_info(anc_img_IDs(img_no_ii)).bnd_pix_ind)=1;

    img_act_pad=(imfill(img_act_pad));
    img_act_pad=uint8(logical(img_act_pad));
    img_act=img_act.*(img_act_pad);
    
    %     I_adj=imadjust(img_act,stretchlim(img_act,0.01),[0 1]);
    I_adj=imadjust(img_act,stretchlim(img_act,0.0),[0 1]);
    
    I_adj= I_adj+( I_adj-imgaussfilt( I_adj,1*100,'FilterSize',[3 3]*301,...
        'FilterDomain','frequency'))*5;
    
    I_adj=imresize(I_adj, (xy_pix)/(downscaled_xy_pix) );
    
    for current_pos_idx=1:52
        current_ap= ((current_pos_idx)*250/(downscaled_xy_pix))-size(VOL,3)/2;
        current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage1)*cosd(pitch_stage1));
        img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
        img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
        
        rough_initial_search(current_pos_idx) =...
            try_matching0805(I_adj, img_ref, downscaled_xy_pix,BlockSz_Stg1);
    end
    
    AP_arr_at_max_angle_stage1=rough_initial_search;
    [AP_rough_best]=find(rough_initial_search==max(max(rough_initial_search)));
    idx=find(AP_arr_at_max_angle_stage1==max(AP_arr_at_max_angle_stage1,[],2));
    max_APpos_stage1(img_no_ii)=AP_rough_best;
    %     figure; plot(AP_arr_at_max_angle_stage1)
    
    %%%% Get the slice img from the atlas that is least dissimilar to the experiment-obtained img %%%%
    current_ap= ((AP_rough_best)*250/(downscaled_xy_pix))-size(VOL,3)/2;
    current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage1)*cosd(pitch_stage1));
    img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
    img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
    
    I_adjusted=I_adj;
    figure; imshow(I_adj)
    img_library(img_no_ii).img=I_adjusted;
    
end

computational_time_stage1=toc(computational_time_checker_stage1);
disp(strcat({'   Stage 1   :   '}, {num2str(computational_time_stage1)}, {'  seconds   '}))

%% STAGE 2 : Angle Search at Search Resolution = 8 deg
computational_time_checker_stage2=tic;
disp(strcat('Searching Stage 2 - Find Angles: Angle Resolution = ', num2str(starting_ang_resolution)))

AP_array_2=repmat(-10:10,[no_images, 1])*2/ref_resc;   %% +/- 500um search
AP_array_2=max_APpos_stage1'*250/(downscaled_xy_pix)+AP_array_2; 

yaw_stage2=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]+yaw_stage1; 
pitch_stage2=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]+pitch_stage1;
[yaw_stage2, pitch_stage2]=meshgrid(yaw_stage2,pitch_stage2);
yaw_stage2= yaw_stage2(:); pitch_stage2=pitch_stage2(:);

parfor ii=1:9
    
    tform_yaw_2=transform_matrix_0822( yaw_stage2(ii),[0 1 0]);
    tform_pitch_2=transform_matrix_0822(pitch_stage2(ii),[1 0 0] );
    
    tform_combined_2=mtimes(tform_yaw_2, tform_pitch_2);
    tf = affine3d(tform_combined_2);
    VOL_rot=imwarp(VOL,tf,'cubic');
    
    for img_no_ii=1:no_images
        img_act=img_library(img_no_ii).img;
        for current_pos_idx=1:21
            current_ap= AP_array_2(img_no_ii,current_pos_idx)-size(VOL,3)/2;
            current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage2(ii))*cosd(pitch_stage2(ii)));
            img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
       
            img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
            [ img_similarity]= try_matching0805(img_act, img_ref, downscaled_xy_pix,BlockSz_Stg2);
            angle_2{ii}.img_similarity_arr_AP(img_no_ii,current_pos_idx)=img_similarity;
        end
    end
end

for ii=1:9
    tempo_img_sim_arr= angle_2{ii}.img_similarity_arr_AP;
    for jj=1:size(tempo_img_sim_arr,1)
        tempo_img_sim_arr_single_img=tempo_img_sim_arr(jj,:);
        angle_2_filtered{ii}.img_similarity_arr_AP(jj,:)= tempo_img_sim_arr_single_img;
    end
end

for ii=1:9
    max_val_stage2(ii)=sum(max(angle_2_filtered{ii}.img_similarity_arr_AP,[],2));
end

max_val_stage2=reshape(max_val_stage2,[3 3]);
max_val_stage2_smooth=max_val_stage2(:);

max_ang_idx_2=find(max_val_stage2_smooth==max(max_val_stage2_smooth));
yaw_stage2_max=yaw_stage2(max_ang_idx_2)
pitch_stage2_max=pitch_stage2(max_ang_idx_2)

%%% find Max Sim AP IDX
AP_arr_at_max_angle_stage2= angle_2_filtered{max_ang_idx_2}.img_similarity_arr_AP;

for ii=1:no_images
    tempo_AP_arr_2=AP_arr_at_max_angle_stage2(ii,:);
    tempo_max_sim_2=find(tempo_AP_arr_2==max(tempo_AP_arr_2));
    
    if length(tempo_max_sim_2)>1   %% in case with more than 1 max values
        tempo_padded_AP_arr_2=padarray(tempo_AP_arr_2,[0 1],'replicate','both');
        for jj=1:length(tempo_max_sim_2)
            tempo_whos_higher_2(jj)=...
                nanmean(tempo_padded_AP_arr_2(tempo_max_sim_2(jj):tempo_max_sim_2(jj)+2));
        end
        whos_higher_idx_2=find(tempo_whos_higher_2==max(tempo_whos_higher_2));
        max_APpos_stage2_idx(ii)=tempo_max_sim_2(whos_higher_idx_2);
    else
        max_APpos_stage2_idx(ii)=tempo_max_sim_2;
    end
end

for ii=1:no_images
    max_APpos_stage2(ii)=AP_array_2(ii,max_APpos_stage2_idx(ii));
end

computational_time_stage2=toc(computational_time_checker_stage2);
disp(strcat({'   Stage 2   :   '}, {num2str(computational_time_stage2)}, {'  seconds   '}))

%% STAGE 3 : Angle Search at Search Resolution = 4 deg

computational_time_checker_stage3=tic;
disp(strcat('Searching Stage 3 - Find Angles: Angle Resolution = ', num2str(starting_ang_resolution/denominator_step)))

AP_array_3=repmat(-10:10,[no_images, 1])*2/ref_resc; 
AP_array_3=max_APpos_stage2'+AP_array_3; 

yaw_stage3=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]/denominator_step+yaw_stage2_max;
pitch_stage3=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]/denominator_step+pitch_stage2_max;
[yaw_stage3, pitch_stage3]=meshgrid(yaw_stage3,pitch_stage3);
yaw_stage3= yaw_stage3(:); pitch_stage3=pitch_stage3(:);

parfor ii=1:9
    
    tform_yaw_3=transform_matrix_0822( yaw_stage3(ii),[0 1 0]);
    tform_pitch_3=transform_matrix_0822(pitch_stage3(ii),[1 0 0] );
    
    tform_combined_3=mtimes(tform_yaw_3, tform_pitch_3);
    tf = affine3d(tform_combined_3);
    VOL_rot=imwarp(VOL,tf,'cubic');
    
    for img_no_ii=1:no_images
        img_act=img_library(img_no_ii).img;
        for current_pos_idx=1:21
            current_ap= AP_array_3(img_no_ii,current_pos_idx)-size(VOL,3)/2;
            current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage3(ii))*cosd(pitch_stage3(ii)));
            img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
    
            img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
            [ img_similarity]= try_matching0805(img_act, img_ref, downscaled_xy_pix,BlockSz_Stg3); 
            angle_3{ii}.img_similarity_arr_AP(img_no_ii,current_pos_idx)=img_similarity;
        end
    end
end

for ii=1:9
    tempo_img_sim_arr= angle_3{ii}.img_similarity_arr_AP;
    for jj=1:size(tempo_img_sim_arr,1)
        tempo_img_sim_arr_single_img=tempo_img_sim_arr(jj,:);
        angle_3_filtered{ii}.img_similarity_arr_AP(jj,:)= tempo_img_sim_arr_single_img;
    end
end

for ii=1:9
    max_val_stage3(ii)=sum(max(angle_3_filtered{ii}.img_similarity_arr_AP,[],2));
end

max_val_stage3=reshape(max_val_stage3,[3 3]);
max_val_stage3_smooth=max_val_stage3;
max_val_stage3_smooth=max_val_stage3_smooth(:);

max_ang_idx_3=find(max_val_stage3_smooth==max(max_val_stage3_smooth));

yaw_stage3_max=yaw_stage3(max_ang_idx_3)
pitch_stage3_max=pitch_stage3(max_ang_idx_3)

%%% find Max Sim AP IDX
AP_arr_at_max_angle_stage3=angle_3_filtered{max_ang_idx_3}.img_similarity_arr_AP;

for ii=1:no_images
    tempo_AP_arr_3=AP_arr_at_max_angle_stage3(ii,:);
    tempo_max_sim_3=find(tempo_AP_arr_3==max(tempo_AP_arr_3));
    
    if length(tempo_max_sim_3)>1
        tempo_padded_AP_arr_3=padarray(tempo_AP_arr_3,[0 1],'replicate','both');
        for jj=1:length(tempo_max_sim_3)
            tempo_whos_higher_3(jj)=...
                nanmean(tempo_padded_AP_arr_3(tempo_max_sim_3(jj):tempo_max_sim_3(jj)+2));
        end
        
        whos_higher_idx_3=find(tempo_whos_higher_3==max(tempo_whos_higher_3));
        max_APpos_stage3_idx(ii)=tempo_max_sim_3(whos_higher_idx_3);     
    else
        max_APpos_stage3_idx(ii)=tempo_max_sim_3;
    end 
end

for ii=1:no_images
    max_APpos_stage3(ii)=AP_array_3(ii,max_APpos_stage3_idx(ii));
end

computational_time_stage3=toc(computational_time_checker_stage3);
disp(strcat({'   Stage 3   :   '}, {num2str(computational_time_stage3)}, {'  seconds   '}))


%% STAGE 4: Angle Search at Search Resolution = 2 deg
computational_time_checker_stage4=tic;
disp(strcat('Searching Stage 4 - Find Angles: Angle Resolution = ', num2str(starting_ang_resolution/denominator_step^2)))

AP_array_4=repmat(-10:10,[no_images, 1])*2/ref_resc;  
AP_array_4=max_APpos_stage3'+AP_array_4;  

yaw_stage4=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]/denominator_step^2+yaw_stage3_max;
pitch_stage4=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]/denominator_step^2+pitch_stage3_max;
[yaw_stage4, pitch_stage4]=meshgrid(yaw_stage4,pitch_stage4);
yaw_stage4= yaw_stage4(:); pitch_stage4=pitch_stage4(:);

parfor ii=1:9
    tform_yaw_4=transform_matrix_0822( yaw_stage4(ii),[0 1 0]);
    tform_pitch_4=transform_matrix_0822(pitch_stage4(ii),[1 0 0] );
    
    tform_combined_4=mtimes(tform_yaw_4, tform_pitch_4);
    tf = affine3d(tform_combined_4);
    VOL_rot=imwarp(VOL,tf,'cubic');
    
    for img_no_ii=1:no_images
        img_act=img_library(img_no_ii).img;
        for current_pos_idx=1:21
            current_ap= AP_array_4(img_no_ii,current_pos_idx)-size(VOL,3)/2;
            current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage4(ii))*cosd(pitch_stage4(ii)));
            img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
  
            img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
            [ img_similarity]= try_matching0805(img_act, img_ref, downscaled_xy_pix,BlockSz_Stg4);
            angle_4{ii}.img_similarity_arr_AP(img_no_ii,current_pos_idx)=img_similarity;
        end
    end
end

for ii=1:9
    tempo_img_sim_arr= angle_4{ii}.img_similarity_arr_AP;
    for jj=1:size(tempo_img_sim_arr,1)
        tempo_img_sim_arr_single_img=tempo_img_sim_arr(jj,:);
        angle_4_filtered{ii}.img_similarity_arr_AP(jj,:)= tempo_img_sim_arr_single_img;
    end
end

for ii=1:9
    max_val_stage4(ii)=sum(max(angle_4{ii}.img_similarity_arr_AP,[],2));
end
max_val_stage4=reshape(max_val_stage4,[3 3]);
max_val_stage4_smooth=max_val_stage4; 
max_val_stage4_smooth=max_val_stage4_smooth(:);

max_ang_idx_4=find(max_val_stage4_smooth==max(max_val_stage4_smooth));
yaw_stage4_max=yaw_stage4(max_ang_idx_4)
pitch_stage4_max=pitch_stage4(max_ang_idx_4)

%%% find Max Sim AP IDX
AP_arr_at_max_angle_stage4=angle_4_filtered{max_ang_idx_4}.img_similarity_arr_AP;
for ii=1:no_images
    tempo_AP_arr_4=AP_arr_at_max_angle_stage4(ii,:);
    tempo_max_sim_4=find(tempo_AP_arr_4==max(tempo_AP_arr_4));
    
    if length(tempo_max_sim_4)>1
        tempo_padded_AP_arr_4=padarray(tempo_AP_arr_4,[0 1],'replicate','both');
        for jj=1:length(tempo_max_sim_4)
            tempo_whos_higher_4(jj)=...
                nanmean(tempo_padded_AP_arr_4(tempo_max_sim_4(jj):tempo_max_sim_4(jj)+2));
        end
        whos_higher_idx_4=find(tempo_whos_higher_4==max(tempo_whos_higher_4));
        max_APpos_stage4_idx(ii)=tempo_max_sim_4(whos_higher_idx_4);
    else
        max_APpos_stage4_idx(ii)=tempo_max_sim_4;
    end
end

for ii=1:no_images
    max_APpos_stage4(ii)=AP_array_4(ii,max_APpos_stage4_idx(ii));
end

computational_time_stage4=toc(computational_time_checker_stage4);
disp(strcat({'   Stage 4   :   '}, {num2str(computational_time_stage4)}, {'  seconds   '}))

%% STAGE 5: Angle Search at Search Resolution = 1 deg
computational_time_checker_stage5=tic;
disp(strcat('Searching Stage 5 - Find Angles: Angle Resolution = ', num2str(starting_ang_resolution/denominator_step^3)))

AP_array_5=repmat(-10:10,[no_images, 1])*2/ref_resc;  
AP_array_5=max_APpos_stage4'+AP_array_5; 

yaw_stage5=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]/denominator_step^3+yaw_stage4_max;
pitch_stage5=[-starting_ang_resolution:starting_ang_resolution:starting_ang_resolution]/denominator_step^3+pitch_stage4_max;
[yaw_stage5, pitch_stage5]=meshgrid(yaw_stage5,pitch_stage5);
yaw_stage5= yaw_stage5(:); pitch_stage5=pitch_stage5(:);

parfor ii=1:9
    
    tform_yaw_5=transform_matrix_0822( yaw_stage5(ii),[0 1 0]);
    tform_pitch_5=transform_matrix_0822(pitch_stage5(ii),[1 0 0] );
    
    tform_combined_5=mtimes(tform_yaw_5, tform_pitch_5);
    tf = affine3d(tform_combined_5);
    VOL_rot=imwarp(VOL,tf,'cubic');
    
    for img_no_ii=1:no_images
        img_act=img_library(img_no_ii).img;
        for current_pos_idx=1:21
            current_ap= AP_array_5(img_no_ii,current_pos_idx)-size(VOL,3)/2;
            current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage5(ii))*cosd(pitch_stage5(ii)));
            img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
            img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);    
            [ img_similarity]= try_matching0805(img_act, img_ref, downscaled_xy_pix,BlockSz_Stg5);
            angle_5{ii}.img_similarity_arr_AP(img_no_ii,current_pos_idx)=img_similarity;
        end
    end
end

for ii=1:9
    tempo_img_sim_arr= angle_5{ii}.img_similarity_arr_AP;
    for jj=1:size(tempo_img_sim_arr,1)
        tempo_img_sim_arr_single_img=tempo_img_sim_arr(jj,:);
        angle_5_filtered{ii}.img_similarity_arr_AP(jj,:)= tempo_img_sim_arr_single_img;
    end
end

for ii=1:9
    max_val_stage5(ii)=sum(max(angle_5{ii}.img_similarity_arr_AP,[],2));
end

max_val_stage5=reshape(max_val_stage5,[3 3]);
max_val_stage5_smooth=max_val_stage5;
max_val_stage5_smooth=max_val_stage5_smooth(:);

max_ang_idx_5=find(max_val_stage5_smooth==max(max_val_stage5_smooth));
yaw_stage5_max=yaw_stage5(max_ang_idx_5)
pitch_stage5_max=pitch_stage5(max_ang_idx_5)

%%% find Max Sim AP IDX
AP_arr_at_max_angle_stage5=angle_5_filtered{max_ang_idx_5}.img_similarity_arr_AP;

for ii=1:no_images
    tempo_AP_arr_5=AP_arr_at_max_angle_stage5(ii,:);
    tempo_max_sim_5=find(tempo_AP_arr_5==max(tempo_AP_arr_5));
    
    if length(tempo_max_sim_5)>1
        tempo_padded_AP_arr_5=padarray(tempo_AP_arr_5,[0 1],'replicate','both');
        for jj=1:length(tempo_max_sim_5)
            tempo_whos_higher_5(jj)=...
                nanmean(tempo_padded_AP_arr_5(tempo_max_sim_5(jj):tempo_max_sim_5(jj)+2));
        end
        whos_higher_idx_5=find(tempo_whos_higher_5==max(tempo_whos_higher_5));
        max_APpos_stage5_idx(ii)=tempo_max_sim_5(whos_higher_idx_5);
    else
        max_APpos_stage5_idx(ii)=tempo_max_sim_5;
    end
    
end
for ii=1:no_images
    max_APpos_stage5(ii)=AP_array_5(ii,max_APpos_stage5_idx(ii));
end

computational_time_stage5=toc(computational_time_checker_stage5);
disp(strcat({'   Stage 5   :   '}, {num2str(computational_time_stage5)}, {'  seconds   '}))

%% STAGE 6 : Rough AP searching with angles found

computational_time_checker_stage6=tic;
disp(strcat('Searching Stage 6 - Rough AP searching with angles found'))

tform_yaw_6=transform_matrix_0822( yaw_stage5_max,[0 1 0]);
tform_pitch_6=transform_matrix_0822(pitch_stage5_max,[1 0 0] );
tform_combined_6=mtimes(tform_yaw_6, tform_pitch_6);
tf = affine3d(tform_combined_6);
VOL_rot=imwarp(VOL,tf,'cubic');

parfor img_no_ii=1:no_images
    disp(num2str( img_no_ii))
    img_act=img_library(img_no_ii).img;
    
    rough_search=zeros(1,52);
    
    for current_pos_idx=1:52

        current_ap= ((current_pos_idx)*250/(downscaled_xy_pix))-size(VOL,3)/2;
        current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage5_max)*cosd(pitch_stage5_max));
        img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
        img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
        
        rough_search(current_pos_idx) =...
            try_matching0805(img_act, img_ref, downscaled_xy_pix,BlockSz_Stg6);
    end

    AP_arr_at_max_angle_stage6=rough_search;
    [AP_rough_best]=find(rough_search==max(max(rough_search)));
    idx=find(AP_arr_at_max_angle_stage6==max(AP_arr_at_max_angle_stage6,[],2));
    max_APpos_stage6(img_no_ii)=AP_rough_best;
    
end

computational_time_stage6=toc(computational_time_checker_stage6);
disp(strcat({'   Stage 6   :   '}, {num2str(computational_time_stage6)}, {'  seconds   '}))

%% STAGE 7 : Rough AP searching with angles found
computational_time_checker_stage7=tic;
disp('Searching Stage 7 - Final AP searching')

if strcmp(Structure_stain,'DAPI') || strcmp(Structure_stain,'Nissl')
     [VOL, metaVOL] = nrrdread('ara_nissl_25_2017.nrrd');
     VOL_ori=uint8(uint16(rot90(permute(VOL,[3 1 2]),3))/(2^8));
elseif strcmp(Structure_stain,'AutoF') || strcmp(Structure_stain,'BrightField')
    load('AutofluoresenceAtlas.mat');
    VOL_ori=rot90(permute(VOL,[3 1 2]),3);
end

VOL_rot=imwarp(VOL_ori,tf,'cubic');
AP_array_7=repmat(-20:20,[no_images, 1]);  
AP_array_7=max_APpos_stage6'*250/(ref_atlas_vox_res)+AP_array_7; 

parfor ii=1:no_images
    img_act=img_library(ii).img;
    for jj=1:41
        current_ap= AP_array_7(ii,jj)-size(VOL_ori,3)/2;
        current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(pitch_stage5_max)*cosd(yaw_stage5_max));
        img_ref=imresize(uint8(squeeze(VOL_rot(:,:,current_ap))), 1/ref_resc);
        img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
          
        [ img_similarity]= try_matching0805(img_act, img_ref,downscaled_xy_pix,BlockSz_Stg7);
        img_similarity_AP_final(ii,jj)=img_similarity;
    end
end

for ii=1:no_images
    tempo_AP_arr_final=img_similarity_AP_final(ii,:);
    set(gca,'FontSize',16)
    
    tempo_max_sim_final=find(tempo_AP_arr_final==max(tempo_AP_arr_final));
    max(tempo_AP_arr_final)
    if length(tempo_max_sim_final)>1
        tempo_padded_AP_arr_final=padarray(tempo_AP_arr_final,[0 1],'replicate','both');
        for jj=1:length(tempo_max_sim_final)
            tempo_whos_higher_final(jj)=...
                nanmean(tempo_padded_AP_arr_final(tempo_max_sim_final(jj):tempo_max_sim_final(jj)+2));
        end
        
        whos_higher_idx_final=find(tempo_whos_higher_final==max(tempo_whos_higher_final));
        max_APpos_stage_final_idx(ii)=tempo_max_sim_final(whos_higher_idx_final);
        
    else
        max_APpos_stage_final_idx(ii)=tempo_max_sim_final;
    end
    
end

for ii=1:no_images
    max_APpos_stage_final(ii)=AP_array_7(ii, max_APpos_stage_final_idx(ii));
end

computational_time_stage7=toc(computational_time_checker_stage7);
disp(strcat({'   Stage 7   :   '}, {num2str(computational_time_stage7)}, {'  seconds   '}))

%%%%%%%%%%%%%%%%%%% SHOW RESULT %%%%%%%%%%%%%%%%%%%%%%
for ii=1:no_images
    img_act=img_library(ii).img;
    
    current_ap= max_APpos_stage_final(ii)-size(VOL_ori,3)/2;
    current_ap= round(size(VOL_rot,3)/2+current_ap*cosd(yaw_stage5_max)*cosd(pitch_stage5_max));
    img_ref=imresize(uint8(squeeze(VOL_rot(:,:,current_ap))),1/ref_resc);
    
    img_similarity= try_matching0805(img_act, img_ref, downscaled_xy_pix,BlockSz_Stg7);
    figure; imshowpair(img_act,img_ref,'montage')
    title(strcat('ARA Slice Corresponding to the Anchor Img No.',num2str(ii)));
end

clearvars VOL VOL_rot VOL_ori
save('Step_4_Angle_Search_Result','-v7.3');
