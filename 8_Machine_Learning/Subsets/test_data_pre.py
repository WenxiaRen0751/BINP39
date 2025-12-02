# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Author: Wenxia 
Date: 2025-10-01
Description: 
This script is to extract the test data. The classification of 7 classes
    # 0: CD19pos，1: CM, 2: Inter_Mono, 3: non_CM, 4: CD14neg, 5: CD45posCD3pos, 6: rest of singlets generated from the finalized gating strategy, served as the ground truth.
"""
# Import packages
import aligater as ag
import yaml,os，torch
from math import inf
import pandas as pd
import numpy as np
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
    # sample_ids are used to match wtih the id in the yaml file
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
    # Next: Singlet gating
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/singlets/"+date+"-"+plate+"-"+sampleName+"-singlets.png"
    # Setting the Angle
    if sample_ids in correction_14:
        adjustAngle_1 = 10
    elif sample_ids in correction_15:
        adjustAngle_1 = 5
    elif sample_ids in correction_17:
        adjustAngle_1 = 4
    elif sample_ids in correction_19:
        adjustAngle_1 = 14
    else:
        adjustAngle_1 = 1
    # Setting the height
    if sample_ids in correction_16 or sample_ids in correction_17:
        heightScale_1 = 3
    elif sample_ids in correction_20: 
        heightScale_1 = 0.5
    else:
        heightScale_1 = 4 
    # Setting width
    if sample_ids in correction_18:
        widthScale_1 = 4
    elif sample_ids in correction_21:
        widthScale_1 = 6
    else:
        widthScale_1 = 8

    singlets=ag.gatePC(fcs,xCol=FSC, yCol="FSC 488/10-H",name="singlets",center='density', adjustAngle= adjustAngle_1, widthScale=widthScale_1, 
                       heightScale=heightScale_1,parentGate=PBMC,filePlot=None)
    fcs.update(ag.AGgate(singlets,PBMC,FSC,"FSC 488/10-H","singlets"), QC=True, xlim=[25000,214000], ylim=[25000,214000],MFI=True, MFI_type='all', extra_MFI=None)
    # Next: CD45 gating
    ylim = ag.valleySeek(fcs, xCol=CD45,parentGate=singlets, interval=[800,2000],bins=100, scale='bilog', T=500)
 
    CD45_step1 = ag.gateThreshold(fcs, "CD45_step1", xCol=CD45, yCol=SSC, thresh=ylim, parentGate=singlets, orientation='vertical', population="upper")
    CD45_step2= ag.gateThreshold(fcs, name="CD45_step2", xCol=CD45, yCol=SSC, thresh=ylim,parentGate=singlets,
                                   orientation='vertical', population="lower")
    CD45_lim = ag.valleySeek(fcs, xCol=CD3,parentGate=CD45_step2, interval=[0,2000],bins=100, scale='bilog', T=500)
    mean, median, sigma, maxVal = ag.axisStats(fcs(), xCol=CD3, vI=CD45_step2(), bins=100, sigma=3, scale='bilog', T=500)
    CD3_lim = mean + 1.2*sigma
    CD3_lim = ag.inverseTransformWrapper(CD3_lim, scale = 'bilog', T = 500)
    if CD45_lim < CD3_lim:
        CD3_lim = CD45_lim
    xlim = ag.valleySeek(fcs, xCol=CD3,parentGate=CD45_step1, interval=[CD3_lim,2000],bins=100, scale='bilog', T=500)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/CD45/"+date+"-"+plate+"-"+sampleName+"-CD45pos.png"

    if sample_ids in corrections_dict_raw['CD45_CD19 Gating']['CUSTOM_y_to_2500']:
        ylim1=2500
    elif sample_ids in correction_23:
        ylim1=ylim + 100
    elif sample_ids in correction_24:
        ylim1=ylim + 200
    elif sample_ids in correction_25:
        ylim1=ylim + 500
    elif sample_ids in correction_26:
        ylim1=ylim + 1000
    elif sample_ids in correction_27:
        ylim1=ylim -100
    else:
        ylim1=ylim

    if sample_ids in correction_22:
        xlim1=2000
    elif sample_ids in correction_29:
        xlim1=xlim-100
    elif sample_ids in correction_30:
        xlim1=xlim-200
    else:
        xlim1=xlim
    CD45pos = ag.gateCorner(fcs, "CD45pos", xCol=CD3,yCol=CD45, parentGate=singlets,xThresh=xlim1,yThresh=ylim1,xOrientation='lower',
                            yOrientation='upper',scale='bilog',T=500, filePlot=None)
    fcs.update(ag.AGgate(CD45pos,singlets,CD3,CD45,"CD45pos"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)

    # Next: gating the CD45posCD3pos
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/CD45posCD3pos/"+date+"-"+plate+"-"+sampleName+"-CD45posCD3pos.png"
    CD45posCD3pos = ag.gateCorner(fcs, "CD45posCD3pos", xCol=CD3,yCol=CD45, parentGate=singlets,xThresh=xlim1,yThresh=ylim1,xOrientation='upper',
                                  yOrientation='upper',scale='bilog',T=500, filePlot=None)
    fcs.update(ag.AGgate(CD45posCD3pos,singlets,CD3,CD45,"CD45posCD3pos"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    # Next: gating NKT, the CD56pos under CD45posCD3pos
    temphiCD3pos = ag.gateThreshold(fcs,"temphiCD3pos", xCol=CD3, yCol=CD56, parentGate=CD45posCD3pos,orientation='vertical',population='upper', thresh=2500,scale='bilog',T=100,filePlot=None)
    ylim = ag.valleySeek(fcs,parentGate=temphiCD3pos, xCol=CD56, interval=[110,1100], scale='bilog',T=100)
    if sample_ids in corrections_dict_raw['NKT Gating']['CUSTOM_y_to_2000']:
        ylim = 2000
    elif sample_ids in corrections_dict_raw['NKT Gating']['y-200']:
        ylim = ylim-200
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/NKT/"+date+"-"+plate+"-"+sampleName+"-NKT.png"
    NKT=ag.gateThreshold(fcs, "NKT", xCol=CD3, yCol=CD56, parentGate=CD45posCD3pos,orientation='horizontal',population='upper', thresh=ylim,scale='bilog',T=100,filePlot=None)
    fcs.update(ag.AGgate(NKT,CD45posCD3pos,CD3,CD56,"NKT"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=100,MFI=True, MFI_type='all', extra_MFI=None)
    # Next: CD19 gate
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/CD19/"+date+"-"+plate+"-"+sampleName+"-CD19pos.png"
    xlim = ag.valleySeek(fcs,xCol=CD19, interval=[0,2000], parentGate=CD45pos, scale='bilog',T=1000)
    if sample_ids in corrections_dict_raw['CD19 Gating']['xlim=200']:
        xlim=200
    elif sample_ids in corrections_dict_raw['CD19 Gating']['xlim=400']:
        xlim=400
    else:
        xlim=xlim
    CD19neg=ag.gateThreshold(fcs, "CD19neg", xCol=CD19, yCol=CD45, parentGate=CD45pos, population='lower', thresh=xlim,scale='bilog',T=1000,filePlot=None)
    CD19pos=ag.gateThreshold(fcs, "CD19pos", xCol=CD19, yCol=CD45, parentGate=CD45pos, population='upper', thresh=xlim,scale='bilog',T=1000,filePlot=None)
    fcs.update(ag.AGgate(CD19pos,CD45pos,CD19,CD45,"CD19pos"), QC=True, xlim=[0,20000], ylim=[0,20000],scale='bilog', T=1000,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(CD19neg,CD45pos,CD19,CD45,"CD19neg"), QC=True, xlim=[0,20000], ylim=[0,20000],scale='bilog', T=1000,MFI=True, MFI_type='all', extra_MFI=None)                                                  
    # Next: Monocyte gate
    no_clutter_CD19neg = ag.gateThreshold(fcs,"tmp", xCol=CD14, parentGate=CD19neg, population='lower', thresh=10000,scale='bilog', T=500)
    tmp_low_CD16 = ag.gateThreshold(fcs, "tmplowCD16", xCol=CD16, parentGate=no_clutter_CD19neg, population='lower', thresh=800, scale='bilog', T=500)
    xlim_lo = ag.valleySeek(fcs, xCol=CD14, parentGate=tmp_low_CD16, interval=[200,1000], bins=100,sigma=1, scale='bilog', T=500)
    tmp_hi_CD14 = ag.gateThreshold(fcs, "tmphiCD14", xCol=CD14, parentGate=no_clutter_CD19neg, thresh=xlim_lo,population='upper', scale='bilog',T=500)
    ylim = ag.valleySeek(fcs, xCol=CD16, parentGate=tmp_hi_CD14, interval=[300,10000], scale='bilog',T=500)
    tmp_hi_CD16 = ag.gateThreshold(fcs, "tmplowCD16", xCol=CD16, parentGate=no_clutter_CD19neg, population='upper', thresh=ylim, scale='bilog', T=500)
    peak_density = ag.getHighestDensityPoint(fcs,xCol=CD14, yCol = CD16, parentGate=tmp_hi_CD16, bins=100, scale='bilog', T=500)
    peak_density = ag.inverseTransformWrapper(maxVal, scale='bilog',T=500)
    mean,sigma = ag.halfNormalDistribution(fcs,xCol=CD14,parentGate=tmp_hi_CD16, mean=peak_density,direction='left',scale='bilog',T=500)
    xlim_hi = ag.inverseTransformWrapper(vX=mean+3*sigma, scale='bilog',T=500)
    xlim_hi2 = ag.valleySeek(fcs, xCol=CD14, parentGate=tmp_hi_CD14, interval=[0,1000], scale='bilog',T=500)
    xlim_hi = (xlim_hi2+xlim_hi)/2
    # Correct the xlim-hi
    if sample_ids in corrections_dict_raw['Monocyte Gating']['hi-100']:
        xlim_hi = xlim_hi - 100
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['hi-200']:
        xlim_hi = xlim_hi - 200
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['hi+100']:
        xlim_hi = xlim_hi + 100
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['hi+200']:
        xlim_hi = xlim_hi + 200
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['hi+500']:
        xlim_hi = xlim_hi + 500
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_hi=1800']:
        xlim_hi = 1800
    else:
        xlim_hi = xlim_hi
        
    # Correct the xlim_lo
    if sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_lo=1800']:
        xlim_lo = 1800 
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_lo=800']:
       xlim_lo = 800
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_lo=1100']:
        xlim_lo = 1100
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_lo=1500']:
        xlim_lo = 1500
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['lo-200']:
        xlim_lo = xlim_lo - 200
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['lo+100']:
        xlim_lo = xlim_lo + 100
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['lo+200']:
        xlim_lo = xlim_lo + 200
    else:
        xlim_lo = xlim_lo
    # Correct the ylim
    if sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y+1000']:
        ylim = ylim + 1000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y-1000']:
        ylim = ylim - 1000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y-2000']:
        ylim = ylim - 2000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y-3000']:
        ylim = ylim - 3000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y-4000']:
        ylim = ylim - 4000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y=10000']:
        ylim = 10000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y=2000']:
        ylim = 2000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y=1000']:
        ylim = 1000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y=3000']:
        ylim = 3000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y=20000']:
        ylim = 20000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CUSTOM_y=8000']:
        ylim = 8000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['ylim+100']:
        ylim = ylim + 100
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['ylim-200']:
        ylim = ylim - 200
    else:
        ylim = ylim
   
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Monocytes/"+date+"-"+plate+"-"+sampleName+"-Non-classic-Monocytes.png"
    nonclassical_monocytes = ag.gateCorner(fcs, "NonClassicalMonocytes", xCol=CD14, yCol=CD16, parentGate=no_clutter_CD19neg, 
                                           xThresh=xlim_hi, yThresh=ylim, xOrientation='upper',yOrientation='upper',scale='bilog', T=500,filePlot=None)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Monocytes/"+date+"-"+plate+"-"+sampleName+"-CD14neg.png"    
    # Get the smallest number of the x values
    if xlim_lo < xlim_hi:
        vPL=[[xlim_lo,-5000],[xlim_lo, ylim],[xlim_hi,ylim],[xlim_hi, 215000]]
        vI = ag.gatePointList(fcs(), xCol=CD14, yCol=CD16, vPL=vPL, population='lower',vI=no_clutter_CD19neg(), bhorizontal=False, scale='bilog', T=500)
        CD14neg = ag.AGgate(vI, no_clutter_CD19neg, name="CD14neg", xCol=CD14, yCol=CD16)
        ag.backGate(fcs,xCol=CD14, yCol=CD16,population=CD14neg,background_population=no_clutter_CD19neg,scale='bilog',T=500,filePlot=fileName, markersize=0.1)
    elif xlim_lo > xlim_hi:
        vPL=[[xlim_hi, 215000],[xlim_hi,ylim],[xlim_lo, ylim],[xlim_lo,-5000]]
        vI = ag.gatePointList(fcs(), xCol=CD14, yCol=CD16, vPL=vPL, population='lower',vI=no_clutter_CD19neg(), bhorizontal=False, scale='bilog', T=500)
        CD14neg = ag.AGgate(vI, no_clutter_CD19neg, name="CD14neg", xCol=CD14, yCol=CD16)
        ag.backGate(fcs,xCol=CD14, yCol=CD16,population=CD14neg,background_population=no_clutter_CD19neg,scale='bilog',T=500,filePlot=fileName, markersize=0.1)
    elif xlim_lo == xlim_hi:
        CD14neg = ag.gateThreshold(fcs, 'CD14neg', xCol=CD14, yCol=CD16, parentGate=no_clutter_CD19neg, thresh=xlim, population='lower', scale='bilog',T=500,filePlot=fileName)

    if sample_ids in corrections_dict_raw['Monocyte Gating']['CD16_range_100_6000']:
        range = [100,6000]
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16_range_2000']:
        range = [2000,10000]
    else:
        range = [1000,10000]
    CD16low = ag.valleySeek(fcs, parentGate=CD14neg, xCol=CD16, interval=range, scale='bilog',T=100)
    if sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=200']:
        CD16low = 200
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=5000']:
        CD16low = 5000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=500']:
        CD16low = 500
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=1000']:
        CD16low = 1000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=2000']:
        CD16low = 2000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=3000']:
        CD16low = 3000
    elif sample_ids in corrections_dict_raw['Monocyte Gating']['CD16lim=4000']:
        CD16low = 4000
    # using gatebox to gate the intermediate monocytes (CD14++CD16+)
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Monocytes/"+date+"-"+plate+"-"+sampleName+"-Intermediate_monocytes.png"
    Intermediate_monocytes = ag.gateBox(fcs,'Intermediate_monocytes', xCol=CD14, yCol=CD16, xThreshRight=10000, yThreshTop=ylim, xThreshLeft=xlim_lo, yThreshBottom=CD16low, 
        Outer=False, parentGate=no_clutter_CD19neg, bins=300, scale='bilog',T=500, update=False,filePlot=fileName, QC=False)
    # using gatecorner to gate the Classical Monocytes(CD14++CD16-)
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Monocytes/"+date+"-"+plate+"-"+sampleName+"-Classical_Monocytes.png"
    classical_monocytes = ag.gateCorner(fcs,"classical_monocytes", xCol=CD14, yCol=CD16, parentGate=no_clutter_CD19neg, xThresh=xlim_lo, yThresh=CD16low, xOrientation='upper',yOrientation='lower',scale='bilog',T=500,filePlot=None)
    # the rest cells from the CD19neg
    Others_index = list(set(CD19neg())-set(classical_monocytes())-set(Intermediate_monocytes())-set(nonclassical_monocytes())-set(CD14neg()))
    Others_from_CD19neg = ag.AGgate(Others_index,CD19neg,CD14,CD16,"Others_from_CD19neg")
    # Now, update the date for each population
    fcs.update(ag.AGgate(nonclassical_monocytes,CD19neg,CD14,CD16,"nonclassical_monocytes"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(Intermediate_monocytes,CD19neg,CD14,CD16,"Intermediate_monocytes"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(classical_monocytes,CD19neg,CD14,CD16,"classical_monocytes"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(CD14neg,CD19neg,CD14,CD16,"CD14neg"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(Others_from_CD19neg,CD19neg,CD14,CD16,"Others_from_CD19neg"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    
    ############################################################################# extracting the ground truth ##########################################################
    # extract data
    df_all =fcs()
    singlets = singlets()
    CD45pos = CD45pos()
    # adding labels: 0~6
    # CD19pos+CM+Inter_Mono+non_CM+CD14neg = CD45pos
    # 0+1+2+3+4+5+6 = singlets
    CD19pos=CD19pos() # label as 0
    # the rest of them make it as CD19neg
    CM = classical_monocytes() # label as 1
    Inter_Mono = Intermediate_monocytes() # label as 2
    non_CM = nonclassical_monocytes() # label as 3
    CD14neg = CD14neg() # label as 4
    CD45posCD3pos = CD45posCD3pos() # label as 5
    Others_from_CD19 = Others_from_CD19neg()
    #est_of_singlets = list(set(singlets) - set(CD45pos) - set(CD45posCD3pos)|set(Others_from_CD19))# label as 6
    rest_of_singlets = list(set(singlets) - set(CD19pos)-set(CM) - set(Inter_Mono) - set(non_CM) - set(CD14neg) - set(CD45posCD3pos))
    # all the features are density markers
    features = ["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A"]
    scaler = StandardScaler()
    scaler.fit(df_all[["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A"]])
    df_all[["CD14","CD19","CD3","CD45","CD16"]] = scaler.transform(df_all[["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A"]])
    features =["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A","CD14","CD19","CD3","CD45","CD16"]
    logtrans_dict = {"BV786 CD14-A":"BV786 CD14-A_v2","BV605 CD19-A":"BV605 CD19-A_v2","Alexa 700 CD3-A":"Alexa 700 CD3-A_v2","APC-H7 CD45-A":"APC-H7 CD45-A_v2","APC CD16-A":"APC CD16-A_v2"}
    for col, new_col in logtrans_dict.items():
        if col == "BV605 CD19-A":
            df_all[new_col] = [ag.transformWrapper(x, 1000, "bilog") for x in  df_all[col]]
        else:
            df_all[new_col] = [ag.transformWrapper(x, 500, "bilog") for x in  df_all[col]]
    features =["BV786 CD14-A","BV605 CD19-A","Alexa 700 CD3-A","APC-H7 CD45-A","APC CD16-A","CD14","CD19","CD3","CD45","CD16",'BV786 CD14-A_v2',
            'BV605 CD19-A_v2',
            'Alexa 700 CD3-A_v2',
            'APC-H7 CD45-A_v2',
            'APC CD16-A_v2']
    # create the dataframe for model input
    df_all.loc[CD19pos, "label"] = 0
    df_all.loc[CM, "label"] = 1
    df_all.loc[Inter_Mono, "label"] = 2
    df_all.loc[non_CM, "label"] = 3
    df_all.loc[CD14neg, "label"] = 4
    df_all.loc[CD45posCD3pos, "label"] = 5
    df_all.loc[rest_of_singlets, "label"] = 6
    df_0 = df_all.loc[singlets, features + ["label"]].copy()
    # prepare the data for model
    x = df_0[features].values
    y = df_0["label"].values
    X = np.vstack(x)
    Y = np.hstack(y)
    X_tensor = torch.tensor(X, dtype=torch.float32)
    Y_tensor = torch.tensor(Y, dtype=torch.long)
    # set the output path
    out_dir = "/home/wenxiaren/Part_3_ML/Dataset/Multiple_Classes_CD3neg_monocytes/test_data"
    os.makedirs(out_dir, exist_ok=True)
    out_pt = os.path.join(out_dir, f"{sample_ids}.pt")
    torch.save({"X": X_tensor, "y": Y_tensor, "ID": sample_ids }, out_pt)
    return {
    "ID": sample_ids,
    "path": out_pt
    }
    
if __name__ == '__main__':
    results = []
    # Load test data paths
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
    
