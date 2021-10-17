profile on
close all
clear all

warning('off')
STEP_0_Parameters;
load('ANO_roi_edge');
load('Step_4_Angle_Search_Result');
img_name=Img_filename_list;
mkdir Image_Analysed
mkdir Image_Analysed_ROI_absent
STEP_0_Parameters;
%% Prepartion : indices for both the anchor and non-anchor imgs
img_idx=anc_img_IDs;
ap_found=max_APpos_stage_final;
img_AP=[];
if length(img_idx) == 1
    img_AP = ap_found;
else
    for img_ID=1:length(img_idx)-1
        img_AP=[img_AP, linspace(ap_found(img_ID),ap_found(img_ID+1),img_idx(img_ID+1)-img_idx(img_ID)+1)];
    end
end

img_AP=round(unique(img_AP));
if Slice_AP_orPA==1
    img_idx=min(img_idx):max(img_idx);
else
    img_idx=max(img_idx):-1:min(img_idx);
end

%% LOAD ATLAS
if strcmp(Structure_stain,'DAPI') || strcmp(Structure_stain,'Nissl')
    [VOL, metaVOL] = nrrdread('ara_nissl_25_2017.nrrd');
    VOL=uint8(uint16(rot90(permute(VOL,[3 1 2]),3))/(2^8));
elseif strcmp(Structure_stain,'AutoF')
    load('AutofluoresenceAtlas.mat');
    VOL=rot90(permute(VOL,[3 1 2]),3);
else
    error('Staining method not recognized');
end

yaw_found=yaw_stage5_max;
pitch_found=pitch_stage5_max;

tform_yaw=transform_matrix_0822( yaw_found,[0 1 0]);
tform_pitch=transform_matrix_0822(pitch_found,[1 0 0] );
tform_combined=mtimes(tform_yaw, tform_pitch);
tf_atlas= affine3d(tform_combined);
VOL_rot=imwarp(VOL,tf_atlas,'cubic');
ANO_rot=imwarp(ANO_roi_edge,tf_atlas,'nearest');

%% Load Cell detection From STEP4to5
load('Step_4to5B_Axon_Detection_Result.mat')
STEP_0_Parameters;
%% MAIN PART
parfor_progress(length(img_AP));
errorneous=false(1,length(img_AP));

for img_ID=1:length(img_AP)
    img_ID
    if ~isempty(img_info(img_idx(img_ID)).slice_window)
        
        img_essence(img_ID).img_AP_pos= img_AP(img_ID);
        
        current_ap=img_AP(img_ID)-size(VOL,3)/2;
        current_ap=round(size(VOL_rot,3)/2+current_ap*cosd(pitch_found)*cosd(yaw_found));
        
        %%%%%%%%%%%%%%%%%%%%%%%% Atlas Slice Prepartion for Transformation %%%%%%%%%%%%%%%%%%%%
        tform_general_resc_factor=1;
        img_ref=uint8(squeeze(VOL_rot(:,:,current_ap)));
        img_ref=imadjust(img_ref,stretchlim(img_ref,0),[0 1]);
        img_ref=imresize(img_ref,tform_general_resc_factor);
        img_ref=padarray(img_ref,round([3000 3000]/(ref_atlas_vox_res/tform_general_resc_factor)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Annotated Slice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        img_ANO=uint8(squeeze(ANO_rot(:,:,current_ap)));
        img_ANO=padarray(img_ANO,round([3000 3000]/(ref_atlas_vox_res)));
        img_ANO=imresize(img_ANO,xy_pix_resc_factor*(ref_atlas_vox_res/xy_pix),'method','nearest','Antialiasing',false);
        
        %%%%%%%%%%%%%%%%%%%%%%%% Slice Preparation for Transformation %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        img_act=imread(img_name{img_idx(img_ID), Color_Channel_Structure});
        try
            img_act=rgb2gray(img_act);
        end
        
        img_act_pad=zeros(size(img_act));
        img_act=img_act(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
            img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
        
        img_act_pad=img_act_pad(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
            img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
        img_act_pad(img_info(img_idx(img_ID)).bnd_pix_ind)=1;
        img_act_pad=(imfill(img_act_pad));
        img_act_pad=uint8(logical(img_act_pad));
        
        img_act=img_act.*(img_act_pad);
        img_act=padarray(img_act,round([3000 3000]/(xy_pix)));
        img_act_reserve=img_act;
        
        img_act=imadjust(img_act,stretchlim(img_act,0.00),[0 1]);
        
        img_act=img_act+( img_act-imgaussfilt( img_act,0.5*201,'FilterSize',[3 3]*603,...
            'FilterDomain','frequency'))*5;
        
        img_act=imresize(img_act,xy_pix/ref_atlas_vox_res*tform_general_resc_factor);
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Find Transformation Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        try
            [matched_img, tform_1, tform_2, tform_3, ...
                matchedPoints_act_400, matchedPoints_ref_400, ...
                matchedPoints_act_200, matchedPoints_ref_200, matchedPoints_act_100, matchedPoints_ref_100, ...
                matchedPoints_act_50, matchedPoints_ref_50] = ...
                tform_finder1015( img_act, img_ref, ref_atlas_vox_res/tform_general_resc_factor, [13 13] );
        catch
            errorneous(img_ID)=true;
            continue;
        end
        
        if isnan(matchedPoints_act_50)   %% if image matching is unsuccessful, skip current iteration
            disp(strcat({'Not enough matching points for img #'},num2str(img_ID)))
            continue;
        end
        
        T1=tform_1.T; T2=tform_2.T; T3=tform_3.T;
        tform_combined=mtimes(T1, T2);
        tform_combined=mtimes(tform_combined, T3);
        tf_general = affine2d(tform_combined);
        
        %%%%%%%%%%%%%%%%%%%%%%%% Find translation %%%%%%%%%%%%%%%%%%%%%%%%%
        img_act_trans=imwarp(img_act,tf_general);
        trans_xcorr = conv2(matched_img,rot90(conj(img_act_trans),2));
        [max_cc, imax] = max(abs(trans_xcorr(:)));
        [ypeak, xpeak] = ind2sub(size(trans_xcorr),imax(1));
        translation_vec=([ypeak, xpeak]-size(img_act_trans));
        translation_vec=[translation_vec(2), translation_vec(1)];
        
        translation_vec=translation_vec*(ref_atlas_vox_res/tform_general_resc_factor)/xy_pix;
        
        img_essence(img_ID).tform.affine=tf_general;
        img_essence(img_ID).tform.trans=translation_vec;
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Get DownScaled Img for 3D Rec. %%%%%%%%%%%%%%%%%%%%%%%
        img_downscaled=(imread(img_name{img_idx(img_ID), Color_Channel_Structure}));
        try
            img_downscaled=rgb2gray(img_downscaled);
        end
                
        img_downscaled=img_downscaled(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
            img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
        img_downscaled=img_downscaled.*img_act_pad;
        img_downscaled=padarray(img_downscaled,round([3000 3000]/(xy_pix)));
        
        img_act_downscale = imwarp(img_downscaled,tf_general);
        img_act_downscale = imtranslate(img_act_downscale,translation_vec);
        
        %%% Match size of images
        whos_bigger=size(img_ANO)-size(img_act_downscale);
        if whos_bigger(1)>=0; x_dim=size(img_ANO,1); else; x_dim=size(img_act_downscale,1); end
        if whos_bigger(2)>=0; y_dim=size(img_ANO,2); else; y_dim=size(img_act_downscale,2); end
        
        img_act_sizeMatch=[x_dim,y_dim]-size(img_act_downscale); img_act_sizeMatch(img_act_sizeMatch<0)=0;
        img_ANO_sizeMatch=[x_dim,y_dim]-size(img_ANO); img_ANO_sizeMatch(img_ANO_sizeMatch<0)=0;
        
        img_act_downscale=padarray(img_act_downscale,img_act_sizeMatch,'post');
        
        img_essence(img_ID).image_downscaled=imresize(img_act_downscale,xy_pix/ref_atlas_vox_res);
                
        %%%%%%%%%%%%%%%%%%%%%%%%% Find Cells in ROI %%%%%%%%%%%%%%%%%%%%%%%
        if max(max(img_ANO))~=0
            %%% Find Cells in Each ROI %%%
            for color_ch_ID=1:length(Color_Channel_Interest)                
                cells_in_this_image = cell_detection_rs(img_ID).Color_Cells(color_ch_ID).cell_locations;
                img_Color=(imread(img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)}));
                try
                    img_Color=rgb2gray(img_Color);
                end
                
                img_Color=img_Color(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
                    img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
                img_Color=img_Color.*img_act_pad;
                img_Color=padarray(img_Color,round([3000 3000]/(xy_pix)));
                
                [img_act_no_scale, tf_general_RB] = imwarp(img_Color,tf_general);
                [xdataT_cell, ydataT_cell] = transformPointsForward(tf_general, cells_in_this_image(:,1), cells_in_this_image(:,2));
                [xdataI_cell, ydataI_cell] = worldToIntrinsic(tf_general_RB, xdataT_cell, ydataT_cell);
                
                [img_act_no_scale, tf_translation_RB] = imtranslate(img_act_no_scale,translation_vec);                
                cell_warped_x = xdataI_cell+translation_vec(1);
                cell_warped_y = ydataI_cell+translation_vec(2);
                
                resc_warp = size(img_act_no_scale);
                
                atlas_exp_scale_fac=(25/xy_pix)*xy_pix_resc_factor;
                
                img_act_no_scale=imresize(img_act_no_scale,xy_pix_resc_factor);
                cell_warped_x = cell_warped_x * xy_pix_resc_factor;
                cell_warped_y = cell_warped_y * xy_pix_resc_factor;           
                
                img_ref_no_scale=imresize(img_ref,atlas_exp_scale_fac);
                
                outputView_no_scale = imref2d(size(img_ref_no_scale));

                matchedPoints_act_400_tuned=matchedPoints_act_400*atlas_exp_scale_fac;
                matchedPoints_ref_400_tuned=matchedPoints_ref_400*atlas_exp_scale_fac;
                
                matchedPoints_act_200_tuned=matchedPoints_act_200*atlas_exp_scale_fac;
                matchedPoints_ref_200_tuned=matchedPoints_ref_200*atlas_exp_scale_fac;
                
                matchedPoints_act_100_tuned=matchedPoints_act_100*atlas_exp_scale_fac;
                matchedPoints_ref_100_tuned=matchedPoints_ref_100*atlas_exp_scale_fac;
                
                matchedPoints_act_50_tuned=matchedPoints_act_50*atlas_exp_scale_fac;
                matchedPoints_ref_50_tuned=matchedPoints_ref_50*atlas_exp_scale_fac;
                                
                try
                    tform_noScale_400 = fitgeotrans(matchedPoints_act_400_tuned,matchedPoints_ref_400_tuned,'pwl');
                    [img_warped_no_scale,~]  = imwarp_custom((img_act_no_scale), tform_noScale_400,...
                        'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);
                    if ~isempty(cell_warped_x)
                        distfun = @(p) norm(transformPointsInverse(tform_noScale_400, p) - double([cell_warped_x, cell_warped_y]));
                        options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
                        [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
                        cell_warped_x = cell_warpedXY(:,1);
                        cell_warped_y = cell_warpedXY(:,2);
                    end
                catch
%                     disp([img_name{img_idx(img_ID), Color_Channel_Structure}, ': PWL fails in 400um scale - Trying LWM']);
%                     try
%                         tform_noScale_400 = fitgeotrans(matchedPoints_act_400_tuned,matchedPoints_ref_400_tuned,'lwm', 12);
%                         [ img_warped_no_scale,~]  = imwarp_custom((img_act_no_scale), tform_noScale_400,...
%                             'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);                        
%                         if ~isempty(cell_warped_x)
%                             distfun = @(p) norm(transformPointsInverse(tform_noScale_400, p) - double([cell_warped_x, cell_warped_y]));
%                             options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
%                             [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
%                             cell_warped_x = cell_warpedXY(:,1);
%                             cell_warped_y = cell_warpedXY(:,2);
%                         end
%                     catch
                        disp(['Error during analyzing : ', img_name{img_idx(img_ID), Color_Channel_Structure}]);
                        disp('Error in tform_400, File will be saved, but the most crude matching was not made properly, double check the original image and process image');
                        img_warped_no_scale = img_act_no_scale;
%                     end
                end
                
                try
                    tform_noScale_200 = fitgeotrans(matchedPoints_act_200_tuned,matchedPoints_ref_200_tuned,'pwl');
                    [ img_warped_no_scale,~]  = imwarp_custom((img_warped_no_scale), tform_noScale_200,...
                        'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);                    
                    if ~isempty(cell_warped_x)
                        distfun = @(p) norm(transformPointsInverse(tform_noScale_200, p) - double([cell_warped_x, cell_warped_y]));
                        options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
                        [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
                        cell_warped_x = cell_warpedXY(:,1);
                        cell_warped_y = cell_warpedXY(:,2);
                    end
                catch
%                     disp([img_name{img_idx(img_ID), Color_Channel_Structure}, ': PWL fails in 200um scale - Trying LWM']);
%                     try
%                         tform_noScale_200 = fitgeotrans(matchedPoints_act_200_tuned,matchedPoints_ref_200_tuned,'lwm', 12);
%                         [ img_warped_no_scale,~]  = imwarp_custom((img_warped_no_scale), tform_noScale_200,...
%                             'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);
%                         if ~isempty(cell_warped_x)
%                             distfun = @(p) norm(transformPointsInverse(tform_noScale_200, p) - double([cell_warped_x, cell_warped_y]));
%                             options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
%                             [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
%                             cell_warped_x = cell_warpedXY(:,1);
%                             cell_warped_y = cell_warpedXY(:,2);
%                         end
%                     catch
                        disp(['PWL fails during analyzing : ', img_name{img_idx(img_ID), Color_Channel_Structure}]);
                        disp('File will be saved, but the 200um matching was not made properly, double check the original image and process image');
%                     end
                end
                
                try
                    tform_noScale_100 = fitgeotrans(matchedPoints_act_100_tuned,matchedPoints_ref_100_tuned,'pwl');
                    [ img_warped_no_scale,~]  = imwarp_custom(img_warped_no_scale, tform_noScale_100,...
                        'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);
                    
                    if ~isempty(cell_warped_x)
                        distfun = @(p) norm(transformPointsInverse(tform_noScale_100, p) - double([cell_warped_x, cell_warped_y]));
                        options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
                        [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
                        cell_warped_x = cell_warpedXY(:,1);
                        cell_warped_y = cell_warpedXY(:,2);
                    end
                catch
%                     disp([img_name{img_idx(img_ID), Color_Channel_Structure}, ': PWL fails in 100um scale - Trying LWM']);
%                     try
%                         tform_noScale_100 = fitgeotrans(matchedPoints_act_100_tuned,matchedPoints_ref_100_tuned,'lwm', 12);
%                         [ img_warped_no_scale,~]  = imwarp_custom(img_warped_no_scale, tform_noScale_100,...
%                             'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);
%                         
%                         if ~isempty(cell_warped_x)
%                             distfun = @(p) norm(transformPointsInverse(tform_noScale_100, p) - double([cell_warped_x, cell_warped_y]));
%                             options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
%                             [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
%                             cell_warped_x = cell_warpedXY(:,1);
%                             cell_warped_y = cell_warpedXY(:,2);
%                         end
%                     catch
                        disp(['PWL fails during analyzing : ', img_name{img_idx(img_ID), Color_Channel_Structure}]);
                        disp('File will be saved, but please check the image afterwards (100um warping was not made)');
%                     end
                end
                
                try
                    tform_noScale_50 = fitgeotrans(matchedPoints_act_50_tuned,matchedPoints_ref_50_tuned,'pwl');
                    [ img_warped_no_scale,~]  = imwarp_custom(img_warped_no_scale, tform_noScale_50,...
                        'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);
                    
                    if ~isempty(cell_warped_x)
                        distfun = @(p) norm(transformPointsInverse(tform_noScale_50, p) - double([cell_warped_x, cell_warped_y]));
                        options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
                        [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
                        cell_warped_x = cell_warpedXY(:,1);
                        cell_warped_y = cell_warpedXY(:,2);
                    end
                catch
%                     disp([img_name{img_idx(img_ID), Color_Channel_Structure}, ': PWL fails in 50um scale - Trying LWM']);
%                     try
%                         tform_noScale_50 = fitgeotrans(matchedPoints_act_50_tuned,matchedPoints_ref_50_tuned,'lwm', 12);
%                         [ img_warped_no_scale,~]  = imwarp_custom(img_warped_no_scale, tform_noScale_50,...
%                             'cubic','OutputView',outputView_no_scale,'WarpRescale',resc_warp);
%                         
%                         if ~isempty(cell_warped_x)
%                             distfun = @(p) norm(transformPointsInverse(tform_noScale_50, p) - double([cell_warped_x, cell_warped_y]));
%                             options = optimset('Display','none','MaxFunEvals', numel(cell_warped_x)*500, 'MaxIter', numel(cell_warped_x)*500);
%                             [cell_warpedXY, err] = fminsearch(distfun, double([cell_warped_x, cell_warped_y]), options);
%                             cell_warped_x = cell_warpedXY(:,1);
%                             cell_warped_y = cell_warpedXY(:,2);
%                         end
%                     catch
                        disp(['PWL fails during analyzing : ', img_name{img_idx(img_ID), Color_Channel_Structure}]);
                        disp('File will be saved, but please check the image afterwards (50um warping was not made)');
%                     end
                end

                %%% Match size of images
                whos_bigger=size(img_ANO)-size(img_warped_no_scale);
                if whos_bigger(1)>=0; x_dim=size(img_ANO,1); else; x_dim=size(img_warped_no_scale,1); end
                if whos_bigger(2)>=0; y_dim=size(img_ANO,2); else; y_dim=size(img_warped_no_scale,2); end
                
                img_warped_sizeMatch=[x_dim,y_dim]-size(img_warped_no_scale); img_warped_sizeMatch(img_warped_sizeMatch<0)=0;
                img_ANO_sizeMatch=[x_dim,y_dim]-size(img_ANO); img_ANO_sizeMatch(img_ANO_sizeMatch<0)=0;                
                
                img_ANO=padarray(img_ANO,img_ANO_sizeMatch,'post');
                img_warped_no_scale=padarray(img_warped_no_scale,img_warped_sizeMatch,'post');
                
                img_ch_thumbnail=imresize(img_warped_no_scale,xy_pix/ref_atlas_vox_res);
                [injection_site_pts] = injection_volume_pix( img_ch_thumbnail );
                img_essence(img_ID).injection.Color(color_ch_ID).valid_pts=injection_site_pts;
                img_essence(img_ID).injection.Color(color_ch_ID).image_thumbnail=img_ch_thumbnail;
                img_essence(img_ID).Color_image_downscaled(color_ch_ID).img=img_warped_no_scale;
                
                %%%% for process monitoring %%%%
                
                monitoring_img=img_warped_no_scale;
                monitoring_img_raw=monitoring_img;
                
                ROI_map=cat(3,false(size(img_ANO)),logical(img_ANO),logical(img_ANO));
                
                monitoring_img=imadjust(monitoring_img,stretchlim(monitoring_img,0),[0 1]);
                monitoring_img=cat(3,monitoring_img,monitoring_img,monitoring_img);
                monitoring_img(ROI_map)=255;
                               
                %%% detect cells across the whole slice image
                cell_detected_all = round([cell_warped_x, cell_warped_y]);
                
                if ~isempty(cell_detected_all)
                    out_bnd_alpha = ref_boundarypad_0809_step5( img_warped_no_scale, xy_pix/xy_pix_resc_factor );
                    out_bnd=inShape(out_bnd_alpha,double(cell_detected_all(:,2)),double(cell_detected_all(:,1)));
                    cell_detected_all=cell_detected_all(~out_bnd,:);
                end                
                
                if ~isempty(cell_detected_all)
                    cell_detected_all_pos_ind=sub2ind(size(img_warped_no_scale),...
                        cell_detected_all(:,2),cell_detected_all(:,1));
                else
                    cell_detected_all_pos_ind=[];
                end
                
                img_essence(img_ID).Color_Cells(color_ch_ID).cell_locations=cell_detected_all;                
                
                image_analyzed_ROI=figure;
                imshow(monitoring_img,[]); hold on
                title(strcat({'Image Name : '}, img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},...
                    {'   - Ch No. '},num2str(color_ch_ID)),'Interpreter', 'none');
                if ~isempty(cell_detected_all)
                    scatter(cell_detected_all(:,1),cell_detected_all(:,2),9,'r','filled')
                end
                
                image_analyzed=figure;
                imshow(monitoring_img_raw,[]); hold on
                title(strcat({'Image Name : '}, img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},...
                    {'   - Ch No. '},num2str(color_ch_ID)),'Interpreter', 'none');
                if ~isempty(cell_detected_all)
                    scatter(cell_detected_all(:,1),cell_detected_all(:,2),9,'r','filled')
                end
                                                
                try
                    save_name=strcat('/Image_Analysed/',img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},'.fig');
                    saveas( image_analyzed_ROI ,[pwd  save_name]);
                    close(image_analyzed_ROI)
                    
                    save_name=strcat('/Image_Analysed_ROI_absent/',img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},'.fig');
                    saveas( image_analyzed,[pwd  save_name]);
                    close(image_analyzed)
                catch
                    disp('If this error message occurs without any error message above, something is wrong during saving stage, not from image processing stage');
                end                
            end
        end
    end
    parfor_progress;
end

delete parfor_progress.txt
save('Step_5_Axon_Detection_Result','img_essence','-v7.3');

profile viewer