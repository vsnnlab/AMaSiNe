# AMaSiNe
This is the script used in "Precise Mapping of Single Neurons by Calibrated 3D Reconstruction of Brain Slices Reveals Topographic Projection in Mouse Visual Cortex" by JH Song, W Choi et al., Cell Reports 2020 (https://doi.org/10.1016/j.celrep.2020.107682)

Version 1.1.1 release! 
You may find the old version in the following:
https://github.com/vsnnlab/AMaSiNe/releases

Before use, please download "Core Functions.zip" and "Images.zip" from following

Core Functions.zip : https://drive.google.com/file/d/1eeXTxUyN9UmdxMT9KAq8HjVq7HNY0jkM/view?usp=sharing

Images.zip : https://drive.google.com/file/d/1UvWJSARRrennZVRuoGgbRWQwFjrAKvfM/view?usp=sharing

and unzip the files in the same directory.

Please read the manual carefully if error occurs

Possible questions #1: If you encounter toolbox version error, please download the essential toolbox for running;
If you encounter toolbox version error even if you have the toolbox, comment the line 8 on STEP_1_Slice_Outline (%toolbox_chk)

Possible questions #2: If you are using Mac environment, you may encounter some error messages in STEP5,
such as elastix error or something like "Unrecognized function or variable 'matched_img'".
In this case, we highly recommend you to use Windows or Linux environment,
or use previous version of our tool (https://github.com/vsnnlab/AMaSiNe/releases/tag/v1.0)


Updates in 17 Oct. 2021
Several pilot functions are uploaded:
1) "STEP_4to5B_Axon_Detection.m" detects axon-fibers using luminance thresholding: If you imaged axon fibers, you may try this function.
The following scripts are "STEP_5B_ForAxon_Transform_and_ROI_drawing.m" and "STEP_7B_Axon_Reconstruction.m"

2) "STEP_7X_Reconstruction_with_ROI_area_and_density.m" calculates number of detected cells, ROI area and density of detected cells in specific ROIs
in each slice image.
If you do not need visualization but only excel results were expected, you may skip the STEP8 and use this script.

Updates in 27 Oct. 2022
Function "STEP_4B_Determine_Parameter_cell_det_thres" is uploaded.
This script let you test some cell_det_thres parameters, before you proceed to step4to5_Cell_Detection. 
You choose a sample image, and let you test several cell_det_thres parameters, and show the results.
