% E158rep2_Tissue_extract
% E158rep2_Tissue_raw was run before this script
%% Inputs
addpath('/Users/iandardani/Code/dentist2') % dentist2
addpath('/Users/iandardani/Code/rajlabimagetools') % if don't have this, wil get this error: Undefined function 'scale' for input arguments of type 'uint16'. Note that it includes a function 'scale.m' which is also the name of functions in various MATLAB toolboxes

% with sigma=0.4
sigma=0.4;
aTrousMinThreshFactor=1.5; % only output spots into spots.csv that are 1.5-fold lower than the threshold (if one is provided), otherwise 1.5-fold lower than the autothreshold for a given block). Eg. if threshold provided is 45, then every spots 30 or greater will be in spots.csv, although all spots <45 will have valid=false
launchGUI=true; % false to run a continuous loop to process spots for all subregions in a batch. Then, afterwards, turn to true to QC check these. The second time you run launchD2ThresholdGUI there will be a spots.csv table in the folder, and it will take these instead of finding spots again

preStitchedScanFilelist={...
'R1_DAPI_50ms.tif',...
'R1_YFP_500ms_UBC.tif',...
'R1_CY3_250ms_NGFR.tif',...
'R1_A594_500ms_MMP1.tif',...
'R1_CY5_250ms_AXL.tif',...
'R2_DAPI_50ms.tif',...
'R2_YFP_500ms_UBC.tif',...
'R2_CY3_250ms_ITGA3.tif',...
'R2_A594_1000ms_FN1.tif',...
'R2_CY5_500ms_EGFR.tif',...
'R3_DAPI_50ms.tif',...
'R3_YFP_500ms_UBC.tif',...
'R3_CY3_1000ms_WNT5A.tif',...
'R3_A594_250ms_DDX58.tif',...
'R3_CY5_100ms_MITF.tif',...
}
channelTypes={'other','FISH','FISH','FISH','FISH','other','FISH','FISH','FISH','FISH','dapi','FISH','FISH','FISH','FISH'};
% dentist2 expects 1 channelType to be 'dapi', which it will use to find nuclei. 'other' DAPI channels are display-only

% spot intensity thresholds are only for the FISH channels.
%              R1_YFP_500ms_UBC          R1_CY3_250ms_NGFR  R1_A594_500ms_MMP1     R1_CY5_250ms_AXL    R2_YFP_500ms_UBC    R2_CY3_250ms_ITGA3    R2_A594_1000ms_FN1  R2_CY5_500ms_EGFR   R3_YFP_500ms_UBC        R3_CY3_1000ms_WNT5A     R3_A594_250ms_DDX58     R3_CY5_100ms_MITF                  
thresholds=   [150                       50                 60                     60                  150                 50                    100                 80                  150                     70                      35                      35                     ];
%% call dentist2 to find spots (if not already found) and allow user to mask regions through GUI
% navigate to a subregion folder (to make it the current folder) and run this section to find spots and open GUI for QC steps
%
%
% subregions of interest:
%
%   scan1_TumorFreshFrozen_Drug
%       Subregion_11_r3_c3
%
%   scan2_TumorFreshFrozen_NoDrug
%       Subregion_4_r1_c4

h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);    
pause(1)