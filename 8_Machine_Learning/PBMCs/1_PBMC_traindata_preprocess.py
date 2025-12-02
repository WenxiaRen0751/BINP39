# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia
Date: 2025-10-01
Description:
This script is to preprocess PBMC data to extract the traindataset for PBMC and non-PBMC classification model training.
It includes gating strategy implementation, feature extraction, normalization, and downsampling.

"""
# Import packages
import yaml
import aligater as ag
from math import inf
import pandas as pd
import numpy as np
import torch
from sklearn.preprocessing import StandardScaler

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
    correction_4 = corrections_dict_raw['PBMC Gating']['Ylim_100K'
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

    ########################################### extract the data for ML model training ##########################################################
    # Next: downsampling the data
    # returns the index number for the gated population
    # extract data
    df_all =fcs() 
    all_index = df_all.index.tolist() 
    # the index number 
    PBMC_index = PBMC() # label as 1
    non_PBMC_index = list(set(all_index) - set(PBMC_index )) # label as 0
    # normalize the features
    scaler = StandardScaler()
    scaler.fit(df_all[["FSC 488/10-A", "SSC 488/10-A", "FSC 488/10-H"]])
    df_all[["FSC_A_Normal", "SSC_A_Normal", "FSC_H_Normal"]] = scaler.transform(df_all[["FSC 488/10-A", "SSC 488/10-A", "FSC 488/10-H"]])
    features =["SSC 488/10-A","FSC 488/10-A","FSC 488/10-H","FSC_A_Normal", "SSC_A_Normal", "FSC_H_Normal"]
    # get the feature dataframe
    df_0 = df_all.loc[:, features].copy()
    # determine the downsample size
    if len(PBMC_index) < 30000:
        PBMC_index_size = len(PBMC_index)
    else:
        PBMC_index_size = 30000
    # downsample to balanced dataset
    pbmc_idxs =  np.random.default_rng(seed=30).choice(PBMC_index, size=PBMC_index_size, replace=False)
    df_PBMC = df_0.loc[pbmc_idxs].sample(frac=1)
    df_PBMC["label"] = 1
    # there are more non-PBMC events, downsample to the same size as PBMC
    non_pbmc_idxs =  np.random.default_rng(seed=30).choice(non_PBMC_index, size=PBMC_index_size, replace=False)
    df_non_PBMC = df_0.loc[non_pbmc_idxs].sample(frac=1)
    df_non_PBMC["label"] = 0
    # combining the downsampled
    df = pd.concat([df_PBMC,df_non_PBMC])
    df = df.reset_index(drop=True)
    x = df[features].values
    y = df["label"].values
    return x,y,sampleName

def traindata_extract(paths, out_pt, out_id_csv):

    """
    Process a list of FCS file paths, gate PBMC populations, extract features and labels,
    and save the combined dataset as a PyTorch tensor file along with an index CSV.
    1. For each FCS file:
       - Load the FCS file using AliGater.
       - Apply PBMC gating strategy to identify PBMC populations.
       - Extract relevant features and labels (PBMC vs non-PBMC).       
    2. Combine data from all files into a single dataset.
    3. Save the combined dataset 

    """
    sample_X, sample_y = [], []
    meta = []  
    offset = 0
    n_th = 0
    for fpath in paths:
        fcs = ag.loadFCS(path=fpath, return_type="agsample",compensate=True, flourochrome_area_filter=True)
        if fcs is None:
            print(f"Skipped {fpath} (failed to load)")
            continue
        x, y, sample_id = gateGeneralDataSet(fcs)  
        n = x.shape[0]
        # append to the list
        sample_X.append(x)
        sample_y.append(y)
        meta.append({"sample_id": sample_id, "file": fpath,
                     "start": offset, "end": offset + n})
        offset += n
        n_th +=1
        print("sample: "+f"{n_th}")
    # combine all samples
    X = np.vstack(sample_X)
    Y = np.hstack(sample_y)
    # convert to tensor and save
    X_tensor = torch.tensor(X, dtype=torch.float32)
    Y_tensor = torch.tensor(Y, dtype=torch.long)
    torch.save({"X": X_tensor, "y": Y_tensor}, out_pt)
    pd.DataFrame(meta).to_csv(out_id_csv, index=False)
    print(f"Saved {out_pt} and {out_id_csv} | total rows={len(X)}; samples={len(meta)}")

if __name__ == '__main__':
    train_data_path = "/home/wenxiaren/Part_3_ML/Data_Processing/train_path.csv"
    with open(train_data_path) as f:
        train_paths = [line.strip() for line in f if line.strip()]
    traindata_extract(train_paths,
        "/home/wenxiaren/Part_3_ML/Dataset/PBMC_dataset/train_data/train_dataset.pt",
        "/home/wenxiaren/Part_3_ML/Dataset/PBMC_dataset/train_data/train_index.csv"
    )

    