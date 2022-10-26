%% STEP_4B_Determine_Parameter_cell_det_thres
profile on
close all; clear all
addpath(genpath(pwd));

warning('off')
STEP_0_Parameters;
img_name=Img_filename_list;
load('Step1_Outline_result.mat')
%% Section 0: Image Selection and parameter range
disp('Image List on the folder')
disp(img_name(:, Color_Channel_Interest))
img_idx = input('Choose the sample image to determine the cell_det_thres:', 's');
img_ID = str2double(img_idx)
disp(['You chose the image named : ', img_name{img_ID, Color_Channel_Interest}]);

cell_det_thresh_candidate = [0.01 0.1 0.15 0.2 0.25];
disp(['This script will show the cell detection result with various cell_det_thres value : ', num2str(cell_det_thresh_candidate)])
disp('If you want to test another value, please edit the parameter "cell_det_thres_candidate", from the Section 0')
%% Section 1: detection from test images
for color_ch_ID=1:length(Color_Channel_Interest)
    if ~isempty(img_info(img_ID).slice_window)
        
        img_Color=(imread(img_name{img_ID, Color_Channel_Interest(color_ch_ID)}));
        try
            img_Color=rgb2gray(img_Color);
        end
        img_Color_pad=zeros(size(img_Color));
        img_Color=img_Color(img_info((img_ID)).slice_window(1):...
            img_info((img_ID)).slice_window(2),...
            img_info((img_ID)).slice_window(3):...
            img_info((img_ID)).slice_window(4));
        img_Color_pad=img_Color_pad(img_info((img_ID)).slice_window(1):...
            img_info((img_ID)).slice_window(2),...
            img_info((img_ID)).slice_window(3):...
            img_info((img_ID)).slice_window(4));
        img_Color_pad(img_info((img_ID)).bnd_pix_ind)=1;
        img_Color_pad=(imfill(img_Color_pad));
        img_Color_pad=uint8(logical(img_Color_pad));
        
        img_Color=img_Color.*img_Color_pad;
        img_Color=padarray(img_Color,round([3000 3000]/(xy_pix)));
        
        for param_iter = 1:numel(cell_det_thresh_candidate)
            cell_det_thresh = cell_det_thresh_candidate(param_iter)
            %%% detect cells across the whole slice image
            cell_detected_all= SomaDetection0827(...
                img_Color, xy_pix, soma_radius,cell_det_thresh);
            
            if ~isempty(cell_detected_all)
                out_bnd_alpha = ref_boundarypad_0809_step5(...
                    img_Color, xy_pix);
                out_bnd=inShape(out_bnd_alpha,cell_detected_all(:,2),...
                    cell_detected_all(:,1));
                cell_detected_all=cell_detected_all(~out_bnd,:);
            end
            
            cell_detected_all= round(cell_detected_all);
            
            if ~isempty(cell_detected_all)
                cell_detected_all_pos_ind=sub2ind(size(img_Color),...
                    cell_detected_all(:,2),cell_detected_all(:,1));
            else
                cell_detected_all_pos_ind=[];
            end            
            
            image_analyzed=figure;
            imshow(img_Color,[]); hold on
            title(strcat({'Image Name : '}, img_name{(img_ID),...
                Color_Channel_Interest(color_ch_ID)},...
                {'   - Ch No. '},num2str(color_ch_ID),...
                {'   cell_det_thres='}, num2str(cell_det_thresh)),...
                'Interpreter', 'none');
            if ~isempty(cell_detected_all)
                scatter(cell_detected_all(:,1),...
                    cell_detected_all(:,2),9,'r','filled')
            end
        end
    end
end