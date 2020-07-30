%% 1. Image Directory
main_folder_dir='C:\Users\junhs\Downloads\AMaSiNe-master\AMaSiNe-with_Elastix_07162020';   
cd(main_folder_dir);
addpath(genpath(main_folder_dir));

%% 2. Image names and slice order
img_format ='jpg';  % if tif = 'tif' ; jpg='jpg'
Slice_AP_orPA= 1;               % If the brain is sliced from anterior to posterior, set this value = 1
                                %                         posterior to anterior, set this value = -1
slide_digit=3;
scene_digit=4;
channel_digit= 5;

%% 3. Anchor Image IDs for Angle Finding (for STEP_2 and 3)
anc_img_IDs= sort([1 4]);
img_IDs_reBoundary=[2];
threshold_scale = 1.5;          %for step 3 only - increase this number if your slice boundary is smaller than what you expect;
                                %decrease this number if slice
                                %boundary is not clean enough.
                                
%% 4. Image Parameters
xy_pix=0.653 * 2;                 % Pixel size = um/pixel   
Name_Channels={'eGFP','DAPI'}; % In the right order
Color_Channel_Structure=2;      % ID of Color Channel to be used for angle finding process
                                % DAPI or NISSL is very strongly recommended for
                                % the angle finding process (Step_4)
Structure_stain={'DAPI'};       % Choose one of the three : 'DAPI','AutoF','Nissl'
Color_Channel_Interest=[1];     % Color channel in which labelled cells are imaged (e.g. eGFP, tdTomato)

%% 5. Detection Parameters
soma_radius=[10 16];            % range of "radius" of labelled soma (in um) to be searched for

cell_det_thresh = 0.25;         % Intensity difference between a cell and its background for a cell to be detected as a cell
                                % Lower this value, you'd get a better chance
                                % of detecting cells dim, but you also risk detecting noise as a cell

%% 6. Allen Atlas Info
size_vol = [528 320 456];       % Matrix size of the Allen 3D atlas
ref_atlas_vox_res=25;           % Allen 3D Ref atlas : 1 voxel= 25um x 25um x 25um
