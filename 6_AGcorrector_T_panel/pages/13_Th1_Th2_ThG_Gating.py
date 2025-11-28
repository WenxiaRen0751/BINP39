# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Title:13_Th1_Th2_ThG_Gating.py
Date: 2025-05-27
Author: Wenxia Ren
"""
from matplotlib import pyplot as plt
import streamlit as st
import aligater as ag
from math import inf
import os
from PIL import Image
import re
import yaml
import json
import scipy.stats as stats
import matplotlib.pyplot as plt

corrections_yaml_file =  "/home/wenxiaren/AGcorrector_T_panel/tmp/T_PBMC_correction.yaml"
with open(corrections_yaml_file, 'r') as stream:
    corrections_dict_raw = yaml.safe_load(stream)

def gateGeneralDataSet(fcs, CCR10_sigma, Th1_2_G_sigma_1, Th1_2_G_sigma_2):
    SSC = "SSC 488/10-A"
    FSC = "FSC 488/10-A"
    FSC_H = "FSC 488/10-H"
    CD39 = "BB515 CD39-A"
    CCR10_A = "PerCP-Cy5.5 CCR10-A"
    CD25 = "PE-Cy7 CD25-A"
    CD127 = "PE CD127-A"
    CCR6 = "PE-Dazzle 594 CCR6-A"
    HLADR = "BV650 HLA-DR-A"
    CCR7 = "BV711 CCR7-A"
    CXCR5 = "BV786 CXCR5-A"
    CXCR3 = "BV421 CXCR3-A"
    CD194 = "BV605 CD194-A"
    CD4 = "BV510 CD4-A"
    CD3 = "Alexa 700 CD3-A"
    CD8 = "APC-H7 CD8-A"
    CD45RA = "APC CD45RA-A"
    extra_MFI_list = [SSC, FSC, FSC_H, CD39, CD25, CD127, CCR6, HLADR, CCR7, CXCR5, CXCR3, CD194, CD3, CD4, CD8, CD45RA]
    
    ag.agconf.ag_verbose=False
    base_path = "/home/wenxiaren/cbio4/data/Aitzkoa/ImmSight/"
    root_path = os.path.relpath(fcs.filePath,base_path)
    Image_ID_prefix_1 = root_path.replace("/","+")
    Image_ID_prefix =Image_ID_prefix_1.replace(".fcs", "")
    ID =  Image_ID_prefix

    def para_extracter(prefix, pattern):
        mathched_pattern = 0
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith(prefix):
                pattern_match = re.search(pattern, keys)
                if pattern_match and ID in ids:
                    mathched_pattern = float(pattern_match.group())
                    break
        return mathched_pattern
    # para_extracter_1 returns a range of value, if not matched, then return None
    def para_extracter_1(prefix, pattern):
        mathched_pattern = None
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith(prefix):
                pattern_match = re.search(pattern, keys)
                if pattern_match and ID in ids:
                    mathched_pattern = list(map(int,pattern_match.groups()))
                    break
        return mathched_pattern
        # PBMC gating

    no_clutter_thresh = 214000
    no_clutter1=ag.gateThreshold(fcs,"no_clutter", "FSC 488/10-A", "FSC 488/10-H",thresh= no_clutter_thresh, orientation='vertical',population="lower")
    no_clutter=ag.gateThreshold(fcs,"no_clutter","FSC 488/10-A", "FSC 488/10-H", parentGate=no_clutter1,thresh=214000, orientation='horizontal',population="lower")
    halfcut_middle = ag.gateThreshold(fcs, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=110000,
                                    parentGate=no_clutter, orientation='vertical', population='upper')
    ylim_bot = ag.densityDelimitation(fcs, xCol='SSC 488/10-A', parentGate=halfcut_middle, interval=[10000,20000], limit_threshold=0.05, direction='left',scale='linear')
    if ylim_bot == inf:
        ylim_bot = 20000
    halfcut_tail = ag.gateThreshold(fcs, name="right_tail",  xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=180000, parentGate=no_clutter, orientation='vertical', population='upper')
    ylim_top = ag.densityDelimitation(fcs, xCol='SSC 488/10-A', parentGate=halfcut_tail, interval=[50000, 125000], limit_threshold=0.2, direction='right',scale='linear')+25000
    if ylim_top == inf:
        ylim_top = 150000
    PBMC_step1 = ag.gateCorner(fcs, name="cut_corner", xCol='FSC 488/10-A', yCol='SSC 488/10-A', xThresh=65000, yThresh=ylim_bot+20000, 
                            xOrientation='lower', yOrientation='lower', Outer=True, parentGate=no_clutter)
    
    # // define the following three parameters: xboundaries, direction and startY
    # // there are four combination for the xboundaries and direction
    # // extract the parameters from YAML file
    xboundaries_1 = None
    direction_1 = 'up'
    # // different startY value
    if ID in corrections_dict_raw["PBMC Gating"]["startY_1 35000"]:
        startY_1 = 35000
    elif ID in corrections_dict_raw["PBMC Gating"]["startY_1 50000"]:
        startY_1 = 50000
    else:
        startY_1= ylim_bot+10000
    PBMC_step2 = ag.horizontalPath(fcs, name="PBMC",
                    xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                    startY=startY_1, endY=ylim_top, xboundaries=xboundaries_1,
                    yboundaries=[ylim_bot-5000,ylim_top+5000],
                    leftRight=True , direction=direction_1,
                    maxStep=2, phi=0.1, bins=100, sigma=1,
                    scale='linear', parentGate=PBMC_step1)
    # // different cut_off_xlim value
    if ID in corrections_dict_raw["PBMC Gating"]["X_Lim 75000"]:
        cut_off_xlim = 75000
    elif ID in corrections_dict_raw["PBMC Gating"]["X_Lim 70000"]:
        cut_off_xlim = 70000
    else:
        cut_off_xlim = 65000
    outline = ag.gateThreshold(fcs, name="outline", xCol=FSC, yCol=SSC, thresh=cut_off_xlim,parentGate=PBMC_step2,orientation='vertical', population="upper",filePlot=None)
    # // different cut_off_ylim value
    if ID in corrections_dict_raw["PBMC Gating"]["Y_Lim 110000"]:
        cut_off_ylim = 110000
    elif ID in corrections_dict_raw["PBMC Gating"]["Y_Lim 120000"]:
        cut_off_ylim = 120000
    elif ID in corrections_dict_raw["PBMC Gating"]["Y_Lim 130000"]:
        cut_off_ylim = 130000
    elif ID in corrections_dict_raw["PBMC Gating"]["Y_Lim 135000"]:
        cut_off_ylim = 135000
    elif ID in corrections_dict_raw["PBMC Gating"]["Y_Lim 140000"]:
        cut_off_ylim = 140000
    else:
        cut_off_ylim = 115000
    outline_1 = ag.gateThreshold(fcs, name="outline_1", xCol=FSC, yCol=SSC, thresh=cut_off_ylim,parentGate=outline,orientation='horizontal', population="lower",filePlot=None)
    # // save the image,PBMCs.png
    PBMC=ag.gatePC(fcs,name="PBMC",xCol=FSC,yCol=SSC,center='centroid', adjustAngle=2, widthScale=8, heightScale = 6, parentGate=outline_1,filePlot=None)
    # // save the image,PBMCs_backgate.png
    ag.backGate(fcs,population=PBMC,background_population=None,xCol=FSC, yCol=SSC,markersize=0.1, filePlot=None)
    fcs.update(ag.AGgate(PBMC,None,FSC,SSC,"PBMC"), QC=True, xlim=[0,215000], ylim=[0,215000],MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    

    # Next: 2)  Singlet gating
    # // setting the adjustAngle
    if ID in corrections_dict_raw["Singlets Gating"]["heightScale 2"]:
         heightScale_1 = 2
    elif ID in corrections_dict_raw["Singlets Gating"]["heightScale 3"]:
        heightScale_1 = 3
    else:
        heightScale_1 = 4
    
    if ID in corrections_dict_raw["Singlets Gating"]["adjustAngle 3"]:
        adjustAngle_1 = 3
    elif ID in corrections_dict_raw["Singlets Gating"]["adjustAngle 4"]:
        adjustAngle_1 = 4
    elif ID in corrections_dict_raw["Singlets Gating"]["adjustAngle 6"]:
        adjustAngle_1 = 6
    else:
        adjustAngle_1 = 1
        
    widthScale_1 = 8

    singlets=ag.gatePC(fcs,xCol=FSC, yCol="FSC 488/10-H",name="singlets",center='density', adjustAngle= adjustAngle_1, widthScale=widthScale_1, heightScale=heightScale_1,parentGate=PBMC,filePlot=None)
    fcs.update(ag.AGgate(singlets,PBMC,FSC,"FSC 488/10-H","singlets"), QC=True,MFI=True, MFI_type='current', extra_MFI=extra_MFI_list) 


    # Next: 3)  CD3pos/neg gating
    # // setting the interval for CD3 
    interval_1 = [300, 2000]
    xlim = ag.valleySeek(fcs=fcs, xCol=CD3, parentGate=singlets, interval=interval_1, require_local_min=True,
                        scale='bilog', T=200)
    if xlim == ag.np.inf:
        xlim = 400

    CD3pos = ag.gateThreshold(fcs=fcs, name='CD3pos', xCol=CD3, yCol="FSC 488/10-H", parentGate=singlets, thresh=xlim, population='upper', scale='bilog', T=200, filePlot=None)
    CD3neg = ag.gateThreshold(fcs=fcs, name='CD3neg', xCol=CD3, yCol="FSC 488/10-H", parentGate=singlets, thresh=xlim,population='lower', scale='bilog', T=200)
    fcs.update(ag.AGgate(CD3pos, singlets, xCol=CD3, yCol="FSC 488/10-H", name="CD3pos"), QC=True, MFI=True,MFI_type='current', extra_MFI=extra_MFI_list,scale='bilog', T=200)
    fcs.update(ag.AGgate(CD3neg, singlets, xCol=CD3, yCol="FSC 488/10-H", name="CD3neg"), QC=True, MFI=True,MFI_type='current', extra_MFI=extra_MFI_list,scale='bilog', T=200)

    # Next: 4)  CD4/CD8 gating
    # // Setting the interval of the CD4 
    if ID in corrections_dict_raw['CD4_CD8 Gating']["CD4_interval [300, 2000]"]:
        interval_2 = [300, 2000]
    elif ID in corrections_dict_raw['CD4_CD8 Gating']["CD4_interval [5000, 10000]"]:
        interval_2 = [5000, 10000]
    elif ID in corrections_dict_raw['CD4_CD8 Gating']["CD4_interval [4000, 10000]"]:
        interval_2 = [4000, 10000]
    elif ID in corrections_dict_raw['CD4_CD8 Gating']["CD4_interval [300, 2000]"]:
        interval_2 = [300, 2000]
    elif ID in corrections_dict_raw['CD4_CD8 Gating']["CD4_interval [1000, 10000]"]:
        interval_2 =  [1000, 10000]
    else:
        interval_2 = [0, 10000]
    xlim = ag.valleySeek(fcs, xCol=CD8, parentGate=CD3pos, interval=interval_2, require_local_min=True, scale='bilog',T=200)
    if xlim == ag.np.inf:
        xlim = 1500
    CD8neg_tmp = ag.gateThreshold(fcs=fcs, name='CD8neg_temp', xCol=CD8, parentGate=CD3pos, thresh=xlim,population='lower', scale='bilog', T=200)
    # // Setting the interval of the CD8

    if ID in corrections_dict_raw['CD4_CD8 Gating']["CD8_interval [200, 1000]"]:
        interval_3 = [200, 1000]
    else:
        interval_3 = [50, 1000]
    ylim = ag.valleySeek(fcs, xCol=CD4, parentGate=CD8neg_tmp, interval=interval_3, require_local_min=True, scale='bilog', T=200)
    if ylim == ag.np.inf:
        ylim = 200
    if ylim < 100: 
        ylim = 200  

    CD4pos, CD4CD8_doublepos, CD8pos, CD4CD8_doubleneg = ag.quadGate(fcs, names=["CD4pos", "CD4posCD8pos", "CD8pos","CD4negCD8neg"], xCol=CD8, yCol=CD4,
                                                                    parentGate=CD3pos, xThresh=xlim, yThresh=ylim,scale='bilog', T=200, filePlot=None)
    fcs.update(ag.AGgate(CD4pos, CD3pos, xCol=CD8, yCol=CD4, name="CD4pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list, xlim=[-5000, 25000], ylim=[0, 1000], scale='bilog', T=200)
    fcs.update(ag.AGgate(CD8pos, CD3pos, xCol=CD8, yCol=CD4, name="CD8pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doublepos, CD3pos, xCol=CD8, yCol=CD4, name="CD4posCD8pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doubleneg, CD3pos, xCol=CD8, yCol=CD4, name="CD4negCD8neg"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    
    # Next: 6)  CD4pos-N_CM_E_EM
    # the naive, the central memory, the memoryeffector, the effector
    CD45RA_interval = [50,400]
    
    output = para_extracter_1("CD45RA_interval", r'\[(\d+),\s*(\d+)\]')
    if output:
        CD45RA_interval = output
    
    ylim = ag.valleySeek(fcs, xCol=CD45RA, interval=CD45RA_interval, require_local_min=True, parentGate=CD4pos, scale='bilog', T=200)
    if ylim == ag.np.inf:
        ylim=200
    if ylim < 200:
        ylim=200
    # // now, back to gating the naive, the central memory, the memoryeffector, the effector from CD4pos
    tempCD4quad_top = ag.gateThreshold(fcs,xCol=CCR7, yCol=CD45RA,name="topTmp",parentGate=CD4pos, 
                                    thresh = ylim, population='upper',orientation='horizontal',scale='bilog', T=200,update=False)
    mean, median, sigma, maxVal = ag.axisStats(fcsDF = fcs(), xCol = CCR7,vI = tempCD4quad_top(), scale = 'bilog', T = 100)
    if maxVal > median:
        
        CCR7_lim_down_sigma = 3
    
        # // CCR7_lim_down_sigma 
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("CCR7_lim_down_sigma"):
                down_sigma_match = re.search(r'-?\d+\.?\d*$', keys)
                if down_sigma_match and ID in ids:
                    CCR7_lim_down_sigma = float(down_sigma_match.group())
                    break
        
        lim_down = ag.inverseTransformWrapper(maxVal-CCR7_lim_down_sigma*abs(sigma),scale = 'bilog', T = 100)
        # // CCR7_lim_upper_sigma
        CCR7_lim_upper_sigma = 0
        
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("CCR7_lim_upper_sigma"):
                upper_sigma_match = re.search(r'-?\d+\.?\d*$', keys)
                if upper_sigma_match and ID in ids:
                    CCR7_lim_upper_sigma = float(upper_sigma_match.group())
                    break
        
        lim_upper = ag.inverseTransformWrapper(maxVal+CCR7_lim_upper_sigma*abs(sigma),scale = 'bilog', T = 100)
    else:
        CCR7_lim_down_sigma = 1
        # // CCR7_lim_down_sigma
    
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("CCR7_lim_down_sigma"):
                down_sigma_match = re.search(r'-?\d+\.?\d*$', keys)
                if down_sigma_match and ID in ids:
                    CCR7_lim_down_sigma = float(down_sigma_match.group())
                    break
        
        lim_down = ag.inverseTransformWrapper(maxVal-CCR7_lim_down_sigma*abs(sigma),scale = 'bilog', T = 100)
        # // CCR7_lim_upper_sigma
        CCR7_lim_upper_sigma = 0
        
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("CCR7_lim_upper_sigma"):
                upper_sigma_match = re.search(r'-?\d+\.?\d*$', keys)
                if upper_sigma_match and ID in ids:
                    CCR7_lim_upper_sigma = float(upper_sigma_match.group())
                    break
    
        lim_upper = ag.inverseTransformWrapper(maxVal+CCR7_lim_upper_sigma*abs(sigma),scale = 'bilog', T = 100)
    xlim = ag.valleySeek(fcs,xCol=CCR7,interval=[lim_down,lim_upper],require_local_min=True, parentGate=tempCD4quad_top,scale='bilog', T=1000)
    if xlim == ag.np.inf:
        xlim=600
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_N_CM_E_EM/"+Image_ID_prefix+"-CD4pos-N_CM_E_EM.png"
    effector_CD4pos, naive_CD4pos, central_memory_CD4pos, effector_memory_CD4pos = ag.quadGate(fcs=fcs, names=['effector_CD4pos', 'naive_CD4pos', 'central_memory_CD4pos', 'effector_memory_CD4pos'],xCol=CCR7, yCol=CD45RA, xThresh=xlim, yThresh=ylim, parentGate=CD4pos, scale='bilog', T=200, filePlot=None)
    fcs.update(ag.AGgate(effector_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="effector_CD4pos"), QC=True, MFI=True,MFI_type='current', extra_MFI=extra_MFI_list, 
                        xlim=[0, 3000], ylim=[-200, 1000],scale='bilog', T=200)
    fcs.update(ag.AGgate(naive_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="naive_CD4pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(central_memory_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="central_memory_CD4pos"), QC=True,MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(effector_memory_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="effector_memory_CD4pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    
    # Next: 8)  Internediate_Th gating
    # gating the internediate Th1, Th2, Th17 based on the parent gate: non-naive CD4 T cells
    non_naive_T = ag.AGgate(effector_CD4pos() + central_memory_CD4pos() + effector_memory_CD4pos(), CD4pos, CCR7, CD45RA, name="non_naive_T")
    fcs.update(ag.AGgate(non_naive_T, CD4pos, xCol=CCR6, yCol=CD45RA, name="non_naive_T"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    xlim = ag.valleySeek(fcs, xCol=CCR6, interval=[100,1500], require_local_min=True, parentGate=CD4pos, scale='bilog', T=200)
    if xlim < 300:
        xlim = ag.valleySeek(fcs, xCol=CCR6, interval=[500,1500], require_local_min=True, parentGate=CD4pos, scale='bilog', T=200)
    elif xlim > 800:  
        Th_CCR6_interval = [800,5000]
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("Th_CCR6_interval"):
                interval_match = re.search(r'\[(\d+),\s*(\d+)\]', keys)
                if interval_match and ID in ids:
                    Th_CCR6_interval = list(map(int,interval_match.groups()))
                    break
        xlim = ag.valleySeek(fcs, xCol=CCR6, interval=Th_CCR6_interval, require_local_min=True, parentGate=CD4pos, scale='bilog', T=200)
    temp = ag.gateThreshold(fcs,xCol=CCR6, yCol=CD194,name="topTmp",parentGate=non_naive_T, thresh = xlim, population='upper',orientation='vertical',scale='bilog', T=200,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CD194, vI=temp(),scale="bilog", T=1000)
    if ID in corrections_dict_raw['CD194_sigma_1 -1.8']:
        CD194_sigma_1 = -1.8
    elif ID in corrections_dict_raw['CD194_sigma_1 -2.0']:
        CD194_sigma_1 = -2
    elif ID in corrections_dict_raw['CD194_sigma_1 -2.5']:
        CD194_sigma_1 = -2.5
    elif ID in corrections_dict_raw['CD194_sigma_1 5.0']:
        CD194_sigma_1 = 5
    else:
        CD194_sigma_1 = 0.5
    max_1 = ag.inverseTransformWrapper(vX=max+CD194_sigma_1*sigma, scale='bilog', T=1000)
    if max_1 > 3000:
        if ID in corrections_dict_raw['CD194_sigma_2 3.0']:
            CD194_sigma_2 = 3
        elif ID in corrections_dict_raw['CD194_sigma_2 4.0']:
            CD194_sigma_2 = 4
        elif ID in corrections_dict_raw['CD194_sigma_2 6.0']:
            CD194_sigma_2 = 6
        else:
            CD194_sigma_2 = 1.5
        max_1  = ag.inverseTransformWrapper(vX=max-CD194_sigma_2*sigma, scale='bilog', T=1000)
    ylim = max_1
    if ID in corrections_dict_raw['Th1_2_17_ylim_1000']:
        ylim = 1000

    Th1_Th2_ThG_part1,Th22_Th17,Th9,Th1_Th2_ThG_part2 = ag.quadGate(fcs=fcs, names=['Th1_Th2_ThG_part1','Th22_Th17','Th9','Th1_Th2_ThG_part2'], 
                                        xCol=CCR6, yCol=CD194, xThresh=xlim, yThresh=ylim, parentGate=non_naive_T, scale='bilog', T=1000, filePlot=None)
    Th1_Th2_ThG = ag.AGgate(Th1_Th2_ThG_part1() + Th1_Th2_ThG_part2(), non_naive_T, CCR6, CD194, name="Th1_Th2_ThG")
    fcs.update(ag.AGgate(Th1_Th2_ThG, non_naive_T, xCol=CCR6, yCol=CD194, name="h1_Th2_ThG"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list,xlim=[0, 3000], ylim=[-200, 1000], scale='bilog', T=1000)
    fcs.update(ag.AGgate(Th22_Th17, non_naive_T, xCol=CCR6, yCol=CD194, name="Th22_Th17"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th9, non_naive_T, xCol=CCR6, yCol=CD194, name="Th9"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    ############################################# the main part: Th1, Th2, ThG gating #############################################
    ## Th1, Th2 & ThG subpopulation
    Th1_Th2_ThG = ag.AGgate(Th1_Th2_ThG_part1() + Th1_Th2_ThG_part2(), non_naive_T, CCR6, CD194, name="Th1_Th2_ThG")
    fcs.update(ag.AGgate(Th1_Th2_ThG, non_naive_T, xCol=CCR6, yCol=CD194, name="non_naive_T"), QC=False, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    ### calculate the CCR10_A
    temphiCCR10 = ag.gateThreshold(fcs,xCol=CCR10_A, yCol=CXCR3,name="temphiCCR10",parentGate=Th1_Th2_ThG, thresh=0, population='upper',orientation='vertical',scale='bilog', T=500,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CCR10_A,vI=temphiCCR10(),scale="bilog", T=500)
    CCR10_lim = ag.inverseTransformWrapper(vX=max+CCR10_sigma*sigma, scale='bilog',T=500)
    # the default is 2.5
    ### calculate the CXCR3
    temphiCXCR3 = ag.gateThreshold(fcs,xCol=CCR10_A,yCol=CXCR3,name="temphiCXCR3",parentGate=Th1_Th2_ThG, thresh=0, population='upper',orientation='horizontal',scale='bilog', T=500,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CXCR3,vI=temphiCXCR3(),scale="bilog", T=500)
    sigma_converted= ag.inverseTransformWrapper(vX=sigma, scale='bilog',T=500)
    if sigma_converted > -400:
        CXCR3_lim = ag.inverseTransformWrapper(vX=max+Th1_2_G_sigma_1*sigma, scale='bilog',T=500)
        # the default is 0.7
    else:
        CXCR3_lim = ag.inverseTransformWrapper(vX=max+Th1_2_G_sigma_2*sigma, scale='bilog',T=500)
        # the default is 1.5
    fileName1="/home/wenxiaren/AGcorrector_T_panel/tmp"+Image_ID_prefix+"-Th_1_2_G.png"
    Th1,NS,ThG,Th2 = ag.quadGate(fcs=fcs, names=['Th1','NS','ThG','Th2'], xCol=CCR10_A, yCol=CXCR3, xThresh=CCR10_lim, yThresh=CXCR3_lim, parentGate=Th1_Th2_ThG, scale='bilog', T=500, filePlot=fileName1)
    return sigma_converted, fileName1

# Configure layout
st.set_page_config(layout="wide")

# --- SESSION STATE: to track which sample we're viewing ---
if "current_sample_index" not in st.session_state:
    st.session_state.current_sample_index = 0

# Sidebar: parameter controls
with st.sidebar:  
    st.markdown("### Parameters")  
    CCR10_sigma= st.number_input("CCR10_sigma", min_value=-10.0, max_value=10.0,step=1.,format="%.1f")
    Th1_2_G_sigma_1= st.number_input("Th1_2_G_sigma_1", min_value=-10.0, max_value=10.0,step=1.,format="%.1f")
    Th1_2_G_sigma_2= st.number_input("Th1_2_G_sigma_2", min_value=-10.0, max_value=10.0,step=1.,format="%.1f")
    st.markdown("---")
    # Parameter selection for Update/Confirm
    st.markdown("### Select Parameters to Apply")
    parameter_options = [
        "CCR10_sigma",
        "Th1_2_G_sigma_1",
        "Th1_2_G_sigma_2"
       ]
    selected_parameters = st.multiselect(
        "Choose which parameters to apply:",
        options=parameter_options,
        default=parameter_options,  # All selected by default
        help="Select only the parameters you want to update/save"
    )
    st.markdown("---")
    Update = st.button("Update Image")
    Confirm = st.button("Save Parameters")
    if st.button("Done! Check next !"):
        st.session_state.current_sample_index += 1  # Go to next sample
    if st.button("Go BacK! Check the previous !"):
        st.session_state.current_sample_index -= 1  # Go to the previous one

row1 = st.columns(1)
for col in row1:
    content_container = col.container(height=160)
    with content_container:
        col1, col2, col3 = st.columns([1,16,1])
        with col2:
            st.title("Upload Files")
            # update the ID file
            st.markdown("""**Please Upload the ID file**""")
            Smaple_IDs = st.file_uploader("Upload ID", type="txt")
            # update the PATH file
            st.markdown("""**Please Upload the PATH file**""")
            Path_file = st.file_uploader("Upload PATH", type="txt")
            # update the PNG file
            st.markdown("""**Please Upload the PNG File**""")
            Image_folder = st.file_uploader("Upload images", type="png",accept_multiple_files=True)
            # the st.file_uploader will return many attributes
# set an empty list
sample_ids = []
if Smaple_IDs:
    sample_ids = [line.decode('utf-8').strip() for line in Smaple_IDs.readlines()]
whole_path_file = []
if Path_file:
    whole_path_file = [line.decode('utf-8').strip() for line in Path_file.readlines()]

# For row2, there are two sub-rows, the first one the original image from the input 
row2 = st.columns(1)
Images = []
for col in row2:
    # Layout for image display
    content_container = st.container(height=800)
    with content_container:
        col1, col2, col3 = st.columns([1,16,1])
        with col2:
            st.header("Current Image")
            # Check if the input files are correctly uploaded
            if Smaple_IDs and Image_folder:
                current_idx = st.session_state.current_sample_index
                if current_idx < len(sample_ids):
                    current_id = sample_ids[current_idx]
                    # load the ID of the current sample
                    st.markdown(
    f"<div style='font-size:24px; font-weight:600;'>Sample {current_idx + 1}/{len(sample_ids)}: {current_id} </div>",
    unsafe_allow_html=True
)
                    # Extract images that match the current sample ID, then show the image
                    matched_images = []
                    for img in Image_folder:
                        # img.name function relates with the attributes of st.file_uploader
                        imd_name = img.name
                        # extract the name, the extracted name should be identical with name in the ID files
                        img_id = re.sub(r'-([^-\s]+)\.png$','',imd_name)
                        # collect all the images for the same sample in one folder
                        if current_id in img_id:
                            matched_images.append(img)

                    # tracing back to find the path fot the original fcs file
                    # get the fcs data in case the images need to be updated later when press "Update Image" button
                    modified_id  = current_id.replace("+","/")+".fcs"
                    modified_id = "/home/wenxiaren/cbio4/data/Aitzkoa/ImmSight/" + modified_id
                    for path in whole_path_file:
                        #f re.search(modified_id, path):
                        if path ==  modified_id:
                            st.session_state['current_sample_path'] = path
                            st.session_state['current_id'] = current_id
                            st.session_state['fcs'] = ag.loadFCS(path, return_type="agsample", flourochrome_area_filter=True)
                    # load the original images to check 
       
                    if matched_images:
                        cols = st.columns(4)
                        for i, img in enumerate(matched_images[:4]):
                            with cols[i]:
                                st.image(Image.open(img), use_container_width=True)
                    else:
                        st.warning(f"No images found for ID: {current_id}")
                else:
                    st.success("✅ All samples reviewed!")
       
            # if the images in the first column needs to be modified, then input the original fcs data 
            st.header("Updated Image")
            if Update:
                if not selected_parameters:
                    st.warning("⚠️ Please select at least one parameter to apply!")
                elif 'fcs' in st.session_state and st.session_state['fcs'] is not None:
                    fcs = st.session_state['fcs']
                    PATH = fcs.filePath
                    st.write(f"PATH:{PATH}")
                    
                    # Get default values for parameters not selected
                    default_values = {
                        "CCR10_sigma":2.5,
                        "Th1_2_G_sigma_1":0.7,
                        "Th1_2_G_sigma_2":1.5
                
                    }
                    # Create parameter dictionary with selected parameters or defaults
                    final_params = {}
                    for param in parameter_options:
                        if param in selected_parameters:
                            final_params[param] = locals()[param]  # Use current value
                        else:
                            final_params[param] = default_values[param]  # Use default value
                    
                    st.info(f"🔧 Applying selected parameters: {', '.join(selected_parameters)}")
                    
                    # call the function with final parameters
                    sigma_converted, fileName1 = gateGeneralDataSet(
                        fcs, 
                        final_params['CCR10_sigma'],
                        final_params['Th1_2_G_sigma_1'],
                        final_params['Th1_2_G_sigma_2']
                    )
                    st.write(f"Th1_2_G: {sigma_converted}") 
                    cols = st.columns(4)
                    with cols[0]:
                        st.image(fileName1, use_container_width=True)
            
  
            if Confirm:
                if not selected_parameters:
                    st.warning("⚠️ Please select at least one parameter to save!")
                elif 'current_id' in st.session_state:
                    current_id = st.session_state['current_id']
                    
                    # Only create dictionary entries for selected parameters
                    para_dict = {}
                    if "Th1_2_G_sigma_1" in selected_parameters:
                        para_dict[f"Th1_2_G_sigma_1 {Th1_2_G_sigma_1}"] = Th1_2_G_sigma_1
                    if "Th1_2_G_sigma_2" in selected_parameters:
                        para_dict[f"Th1_2_G_sigma_2 {Th1_2_G_sigma_2}"] = Th1_2_G_sigma_2
                    if "CCR10_sigma" in selected_parameters:
                        para_dict[f"CCR10_sigma {CCR10_sigma}"] = CCR10_sigma

                    YAML_PATH = "/home/wenxiaren/AGcorrector_T_panel/tmp/Th_sub.yaml"
                    if os.path.exists(YAML_PATH):
                        with open(YAML_PATH,"r") as in1:
                            correction_yaml = yaml.safe_load(in1) or {}
                    else:
                        correction_yaml = {}
                    
                    # Step 3: Update only selected parameters
                    for key, val in para_dict.items():
                        correction_yaml.setdefault(key, [])
                        if current_id not in correction_yaml[key]:
                            correction_yaml[key].append(current_id)
                    
                    # Step 4: Save
                    with open(YAML_PATH, "w") as out1:
                        yaml.dump(correction_yaml, out1, sort_keys=False)
                    
                    st.success(f"✅ Saved selected parameters ({', '.join(selected_parameters)}) for sample: {current_id}")
                else:
                    st.error("❌ No current sample ID found!")



        

            

    
