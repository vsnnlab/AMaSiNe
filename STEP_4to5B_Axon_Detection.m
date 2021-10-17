%% Code written by Woochul Choi at 2021.Oct
%%% This script determines threshold, to detect axon-like fibers in the images
%%% Threshold will be detected in two alternative ways:
%%% 1) From raw histogram of luminance
%%% 2) From a) fibermetric function then b) thresholding fibermetric func. output
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
xy_pix_resc_factor = 1;
%% MAIN PART
for img_ID=1:length(img_AP)
    disp(['Cell Detection: ', num2str(img_ID), '/', num2str(length(img_AP))]);
    img_ID
    for color_ch_ID=1:length(Color_Channel_Interest)
        if ~isempty(img_info(img_idx(img_ID)).slice_window)
            cell_detection_rs(img_ID).img_AP_pos= img_AP(img_ID);
            
            img_Color=(imread(img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)}));
            try
                img_Color=rgb2gray(img_Color);
            end                        
            img_Color_pad=zeros(size(img_Color));      
            img_Color=img_Color(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
                img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
            
            img_Color_pad=img_Color_pad(img_info(img_idx(img_ID)).slice_window(1):img_info(img_idx(img_ID)).slice_window(2),...
                img_info(img_idx(img_ID)).slice_window(3):img_info(img_idx(img_ID)).slice_window(4));
            img_Color_pad(img_info(img_idx(img_ID)).bnd_pix_ind)=1;
            img_Color_pad=(imfill(img_Color_pad));
            img_Color_pad=uint8(logical(img_Color_pad));
        
            img_Color=img_Color.*img_Color_pad;
            img_Color=padarray(img_Color,round([3000 3000]/(xy_pix)));
            
            atlas_exp_scale_fac=(25/xy_pix)*xy_pix_resc_factor;
            
            img_Color=imresize(img_Color,xy_pix_resc_factor);
            
            close all;
            %% Method 1. thresholding from raw luminance
            disp('You chose method 1: Thresholding from raw luminance value');
            disp('Please choose the right upper_bnd value to detect fibers');
            luminance_value = img_Color(:);
            luminance_value = luminance_value(luminance_value~=0);
            upper_bnd1 = 99; %% Upper (100-X)% value will be set as threshold
            thres_1 = prctile(luminance_value, upper_bnd1); 
            
            figure; imagesc(img_Color); colormap gray; axis tight image equal; hold on;  title('Detected Axons from Method 1');
            [row,col] = find(img_Color > thres_1);
            plot(col,row, 'r.');
            
            cell_detected_all = [col, row]; %% treat detected pixel as cells!
            cell_detection_rs(img_ID).Color_Cells(color_ch_ID).cell_locations= cell_detected_all;
                        
            %% Method 2. Apply Fibermetric function, then threshold
%             disp('You chose method 2: Apply fibermetric function then thresholding ');
%             axon_width = 5;
%             disp(['A) Check the fiber width in pixel value: current value = ', num2str(axon_width), 'px']);
%             disp('It might take a while (~min)');
%             tempImg = fibermetric(img_Color, 'ObjectPolarity', 'bright', 'StructureSensitivity', axon_width);
%             figure; imagesc(tempImg); colormap gray; axis tight image equal; hold on; colorbar;
%             title('Detected Fiber from fibermetric func'); colorbar;
%             
%             upper_bnd2 = 99; %% Lower X % value is set as a threshold, increase this value to be more 'accurate' with 'less' number
%             disp(['B) Check the right threshold value: ', num2str(upper_bnd2), '% Upper Bound ']);
%             tempVal = tempImg(tempImg ~= 0);
%             thres_2 = prctile(tempVal, upper_bnd2); 
%             
%             figure; imagesc(img_Color); colormap gray; axis tight image equal; hold on; colorbar;
%             [row,col] = find(tempImg > thres_2);
%             plot(col, row, 'r.');
%             title('Detected Fiber from fibermetric func');            
%             
%             cell_detected_all = [col, row]; %% treat detected pixel as cells!
%             cell_detection_rs(img_ID).Color_Cells(color_ch_ID).cell_locations= cell_detected_all;
                        
            %% Below is same as before
            figure; imshow(img_Color);
            hold on;
            plot(cell_detected_all(:,1), cell_detected_all(:,2), 'r.');
            
            image_analyzed=figure;
            imshow(img_Color,[]); hold on
            title(strcat({'Detected axon from Image: '}, img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},...
                {'   - Ch No. '},num2str(color_ch_ID)),'Interpreter', 'none');
            if ~isempty(cell_detected_all)
                scatter(cell_detected_all(:,1),cell_detected_all(:,2),9,'r','filled')
            end
                        
            save_name=strcat('/Image_Analysed_ROI_absent/Axon_',img_name{img_idx(img_ID), Color_Channel_Interest(color_ch_ID)},'.fig');
            saveas( image_analyzed,[pwd  save_name]);
            close(image_analyzed)            
        end
    end
end
save('Step_4to5B_Axon_Detection_Result','cell_detection_rs','-v7.3');