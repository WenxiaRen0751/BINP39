# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia 
Date: 2025-10-01
Description: 
This script is to extract the test data. The classification of PBMC an non-PBMC generated from the finalized gating strategy 
is served as the ground truth.

"""
# Import packages
import aligater as ag
import yaml,os
from math import inf
import pandas as pd
import numpy as np
import torch

####################################### this is from the finalized gating strategy for PBMC gating  #####################################################
# gating strategy for AliGater
SSC="SSC 488/10-A"
FSC="FSC 488/10-A"
IgA="FITC IgA-A"
IgD="PerCP-Cy5.5 IgD-A"
CD141="PE-Cy7 CD141-A"
CD123="PE-Cy5 CD123-A"
CD27="PE (R-phycoerythrin) CD27-A"
CD34="PE-Dazzle 594 CD34-A"
HLADR="BV650 HLA-DR-A"
CD56="BV711 CD56-A"
CD14="BV786 CD14-A"
CD38="BV421 CD38-A"
CD19="BV605 CD19-A"
CD24="BV510 CD24-A"
CD3="Alexa 700 CD3-A"
CD45="APC-H7 CD45-A"
CD16="APC CD16-A"

# Read the corrections yaml file
corrections_yaml_file ="/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_B_Panel_1006/Scripts/correction.yaml"

with open(corrections_yaml_file, 'r') as stream:
    corrections_dict_raw = yaml.safe_load(stream)
    correction_1 =  corrections_dict_raw['PBMC Gating']['CUSTOM_100000']
    correction_2 = corrections_dict_raw['PBMC Gating']['CUSTOM_95000']
    correction_3 = corrections_dict_raw['PBMC Gating']['CUSTOM_90000']
    correction_4 = corrections_dict_raw['PBMC Gating']['Ylim_100K']
    correction_5 = corrections_dict_raw['PBMC Gating']['Ylim_110K']
    correction_6 = corrections_dict_raw['PBMC Gating']['Ylim_130K']
    correction_7 = corrections_dict_raw['PBMC Gating']['Ylim_90K']
    correction_8 = corrections_dict_raw['PBMC Gating']['xlim_65K']
    correction_9 = corrections_dict_raw['PBMC Gating']['xlim_70K']
    correction_10 = corrections_dict_raw['PBMC Gating']['xlim_80K']
    correction_11 = corrections_dict_raw['PBMC Gating']['No_clutter_240K']
    correction_12 = corrections_dict_raw['PBMC Gating']['Ylim_75K']
    correction_13 = corrections_dict_raw['PBMC Gating']['Dijkstra failure']
    correction_14 = corrections_dict_raw['Singlet Gating']['adjustAngle_10']
    correction_15 = corrections_dict_raw['Singlet Gating']['adjustAngle_5']
    correction_16 = corrections_dict_raw['Singlet Gating']['heightScale_3']
    correction_17 = corrections_dict_raw['Singlet Gating']['CUSTOM_2']
    correction_18 = corrections_dict_raw['Singlet Gating']['widthScale_4']
    correction_19 = corrections_dict_raw['Singlet Gating']['adjustAngle_14']
    correction_20 = corrections_dict_raw['Singlet Gating']['heightScale_1']
    correction_21 = corrections_dict_raw['Singlet Gating']['widthScale_6']
    correction_22 = corrections_dict_raw['CD45_CD19 Gating']['CUSTOM_x_1900']
    correction_23 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+100']
    correction_24 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+200']
    correction_25 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+500']
    correction_26 = corrections_dict_raw['CD45_CD19 Gating']['xlim_+1000']
    correction_27 = corrections_dict_raw['CD45_CD19 Gating']['xlim_-100']
    correction_28 = corrections_dict_raw['CD45_CD19 Gating']['CUSTOM_y_to_2500']
    correction_29 = corrections_dict_raw['CD45_CD19 Gating']['ylim_-100']
    correction_30 = corrections_dict_raw['CD45_CD19 Gating']['ylim_-200']

def gateGeneralDataSet(fcs, *args, **kwargs):
    ag.agconf.ag_verbose=False
    date_plate = ag.getFileName(ag.getParent(fcs.filePath))
    date = date_plate.split("B panel")[0]
    plate = date_plate.split("B panel")[1]
    sampleName=ag.getFileName(fcs.filePath)
    sample_ids = date+"-"+plate+"-"+sampleName
    # Next: PBMC gating
    if  sample_ids in correction_11:
        no_clutter_thresh = 240000
    else:
        no_clutter_thresh = 214000

    no_clutter1=ag.gateThreshold(fcs,"no_clutter", "FSC 488/10-A", "FSC 488/10-H",thresh= no_clutter_thresh, orientation='vertical',population="lower")
    no_clutter=ag.gateThreshold(fcs,"no_clutter","FSC 488/10-A", "FSC 488/10-H", parentGate=no_clutter1,thresh=214000, orientation='horizontal',population="lower")

    halfcut_middle = ag.gateThreshold(fcs, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=110000,
                                    parentGate=no_clutter, orientation='vertical', population='upper')
    ylim_bot = ag.densityDelimitation(fcs, xCol='SSC 488/10-A', parentGate=halfcut_middle, interval=[10000,20000], limit_threshold=0.05, direction='left',scale='linear')
    if ylim_bot == inf:
        ylim_bot = 10000

    halfcut_tail = ag.gateThreshold(fcs, name="right_tail",  xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=180000, parentGate=no_clutter, orientation='vertical', population='upper')
    ylim_top = ag.densityDelimitation(fcs, xCol='SSC 488/10-A', parentGate=halfcut_tail, interval=[50000, 125000], limit_threshold=0.2, direction='right',scale='linear')+25000
    if ylim_top == inf:
        ylim_top = 150000

    PBMC_step1 = ag.gateCorner(fcs, name="cut_corner", xCol='FSC 488/10-A', yCol='SSC 488/10-A', xThresh=65000, yThresh=ylim_bot+20000, 
                            xOrientation='lower', yOrientation='lower', Outer=True, parentGate=no_clutter)
    
    if  sample_ids in correction_13:
        startY_1= 30000
        xboundaries_1 = [70000,214000]
        direction_1 = 'up'
    else:
        startY_1= ylim_bot+10000
        xboundaries_1 = None
        direction_1 = 'both'

    PBMC_step2 = ag.horizontalPath(fcs, name="PBMC",
                    xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                    startY=startY_1, endY=ylim_top, xboundaries=xboundaries_1,
                    yboundaries=[ylim_bot-5000,ylim_top+5000],
                    leftRight=True , direction=direction_1,
                    maxStep=2, phi=0.1, bins=100, sigma=1,
                    scale='linear', parentGate=PBMC_step1)
    if  sample_ids in correction_8:
        cut_off_xlim = 65000
    elif sample_ids in correction_9:
        cut_off_xlim = 70000
    elif sample_ids in correction_10:
        cut_off_xlim = 80000
    elif sample_ids in correction_1:
        cut_off_xlim = 100000
    elif sample_ids in correction_2: 
        cut_off_xlim = 95000
    elif sample_ids in correction_3:
        cut_off_xlim = 90000
    else:
        cut_off_xlim = 75000
    outline = ag.gateThreshold(fcs, name="outline", xCol=FSC, yCol=SSC, thresh=cut_off_xlim,parentGate=PBMC_step2,
                                    orientation='vertical', population="upper",filePlot=None)

    if  sample_ids in correction_4:
        cut_off_ylim = 100000
    elif sample_ids in correction_5:
        cut_off_ylim = 110000
    elif sample_ids in correction_6:
        cut_off_ylim = 130000
    elif sample_ids in correction_7:
        cut_off_ylim = 90000
    elif sample_ids in correction_12:  
        cut_off_ylim = 75000
    else:
        cut_off_ylim = 120000
    outline_1 = ag.gateThreshold(fcs, name="outline_1", xCol=FSC, yCol=SSC, thresh=cut_off_ylim,parentGate=outline,
                                    orientation='horizontal', population="lower",filePlot=None)
    PBMC=ag.gatePC(fcs,name="PBMC",xCol=FSC,yCol=SSC,center='centroid', adjustAngle=2, widthScale=8, heightScale = 6, parentGate=outline_1,filePlot=None)
    ag.backGate(fcs,population=PBMC,background_population=None,xCol=FSC, yCol=SSC,markersize=0.1, filePlot=None)
    fcs.update(ag.AGgate(PBMC,None,FSC,SSC,"PBMC"), QC=True, xlim=[0,215000], ylim=[0,215000],MFI=True, MFI_type='all', extra_MFI=None)

    ################################################ extracting the ground truth ##########################################################
    # extract data
    df_all =fcs()
    all_index = df_all.index.tolist() 
    # the index number 
    PBMC_index = PBMC() # label as 1
    non_PBMC_index = list(set(all_index) - set(PBMC_index)) # label as 0
    features = ["SSC 488/10-A","FSC 488/10-A","FSC 488/10-H"]
    df_0 = df_all.loc[:, features].copy()
    # add ground truth
    df_0.loc[non_PBMC_index, "label"] = 0
    df_0.loc[PBMC_index, "label"] = 1
    # prepare data for pytorch
    x = df_0[features].values
    y = df_0["label"].values
    X = np.vstack(x)
    Y = np.hstack(y)
    # convert to tensor
    X_tensor = torch.tensor(X, dtype=torch.float32)
    Y_tensor = torch.tensor(Y, dtype=torch.long)
    # define the path for the output
    out_dir = "/home/wenxiaren/Part_3_ML/Dataset/PBMC_dataset/test_data"
    os.makedirs(out_dir, exist_ok=True)
    out_pt = os.path.join(out_dir, f"{sample_ids}.pt")
    torch.save({"X": X_tensor, "y": Y_tensor, "ID": sample_ids }, out_pt)

    return {
    "ID": sample_ids,
    "path": out_pt
    }

# Main function to process test data
if __name__ == '__main__':
    results = []
    with open("/home/wenxiaren/Part_3_ML/Data_Processing/test_path.csv") as f:   
        test_paths = [line.strip() for line in f if line.strip()]
        for i in test_paths:
            fcs = ag.loadFCS(path=i, return_type="agsample",
                            compensate=True, flourochrome_area_filter=True)
            if fcs is None:
                print(f"Skipped {i} (failed to load)")
                continue
            info = gateGeneralDataSet(fcs)
            results.append(info)
        print("All samples processed:")
    
