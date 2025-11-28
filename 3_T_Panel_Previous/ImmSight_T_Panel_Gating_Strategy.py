# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Title:ImmSight_T_Panel_Gating_Strategy.py
Date: 2025-06-25
Author: Wenxia Ren
Description:
This program is used for gating the T cell panel in ImmSight project for data generated before 05/30/2025
Note:
1) To run this program, user need to provide the YAML file and the FCS file path. The YAML file contains the correction information for each sample.
2) For the CD39 gating, the threshold for the CD39 is adjusted in the following parentgate: Tregs, CD194posTregs. in the subpopulation of Treg, 
the threshold is as the same as the one in "Tregs". in the subpopulation of CD194posTregs, the threshold is as the same as the one in "CD194posTregs".(we are more interested in these population regarding the CD39).
The rest of the population use the fixed threshold, 1000.
3) For HLADR gating, we use the fixed threshold, 400.
4) The output images are saved in the folder: /home/wenxiaren/aligater/plots/phase_III/ImmSight/T_Panel/
5) The event number and MFT for each population are calculated and saved in txt files
"""

# Import packages
import aligater as ag
import yaml
import os
from math import inf
import re 

# ********************************************************************************************************************************************************************************************************
# Next: Creating output folders 
image_output_folder = "/home/wenxiaren/aligater/plots/phase_III/ImmSight/T_Panel"
if not os.path.exists(image_output_folder):
    os.makedirs(image_output_folder)
subdir = ["PBMCs","singlets","CD3","CD3pos_CD39pos","CD3neg_CD39pos","CD3pos_HLADRpos","CD3neg_HLADRpos","CD4_CD8","CD4_CD8_CD39pos","CD4_CD8_HLADRpos","CD8pos_N_CM_E_EM","CD8pos_CD39pos", 
 'CD8pos_HLADRpos',"CD4pos_CXCR5pos","CD4posCXCR5posCD39pos","CD4posCXCR5posHLADRpos",'CD4posCXCR5pos_Tfh',"Tfh_CD39pos","Tfh_HLADRpos","CD4pos_N_CM_E_EM","CD4pos_CD39pos","CD4pos_HLADRpos",
 "CD4pos_Tregs","CD4pos_Tregs_CD39","CD4pos_Tregs_HLADR","CD194_Tregs","CD4pos_CD194_Tregs_CD39","CD4pos_CD194_Tregs_HLADR","Treg_subpopulations","Treg_subpopulations_CD39","Treg_subpopulations_HLADR",
 "CD194posTreg_subpopulation","CD194treg_subpops_CD39","CD194treg_subpops_HLADR","Th_1_2_17","Th_1_2_17_CD39","Th_1_2_17_HLADR","Th_1_2_G","Th_1_2_G_CD39","Th_1_2_G_HLADR","Th22_Th17_subpopulation","Th17_Th22_CD39","Th17_Th22_HLADR"]
for dir in subdir:
    subdir_path = os.path.join(image_output_folder, dir)
    os.makedirs(subdir_path,exist_ok=True)

# Next: Handling YAML file
corrections_yaml_file = "/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_T_Panel_1006/Scripts/T_correction.yaml"
with open(corrections_yaml_file, 'r') as stream:
    corrections_dict_raw = yaml.safe_load(stream)
# ********************************************************************************************************************************************************************************************************
# Next: Defining the "gateGeneralDataSet" function
def gateGeneralDataSet(fcs, *args, **kwargs):
    # define the marker panel
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
    extra_MFI_list = [SSC, FSC, FSC_H, CD39, CD25, CCR10_A, CD127, CCR6, HLADR, CCR7, CXCR5, CXCR3, CD194, CD3, CD4, CD8, CD45RA]
    # handle the filename
    ag.agconf.ag_verbose=False
    base_path = "/home/wenxiaren/cbio4/data/Aitzkoa/ImmSight/"
    root_path = os.path.relpath(fcs.filePath,base_path)
    Image_ID_prefix_1 = root_path.replace("/","+")
    # Image_ID_prefix will be part of the name of the images
    Image_ID_prefix =Image_ID_prefix_1.replace(".fcs", "")
    # the "ID" here will be used for matching the ID in the YAML file
    ID =  Image_ID_prefix
    # define two functions, one returns a float number, one returns a range of value.
    # para_extracter returns a float number, if not matched, then return 0
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
# ********************************************************************************************************************************************************************************************************
    # Next: 1) PBMC gating
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
    if ID in corrections_dict_raw['direction_up_xboundaries_71000']:
        xboundaries_1 = [71000, 214000]
        direction_1 = 'up'
    elif ID in corrections_dict_raw['direction_up_xboundaries_41000']:
        direction_1 = 'up'
        xboundaries_1 = [41000, 214000]
    elif ID in corrections_dict_raw['direction_up_xboundaries_null']:
        direction_1 = 'up'
        xboundaries_1=None
    else:
        xboundaries_1 = None
        direction_1 = 'both'
    # // different startY value
    if ID in corrections_dict_raw['Another_Cell_Distribution True']:
        startY_1 = 30000
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
    if ID in corrections_dict_raw['X_Lim 65000']:
        cut_off_xlim = 65000
    elif ID in corrections_dict_raw['X_Lim 80000']:
        cut_off_xlim = 80000
    elif ID in corrections_dict_raw['X_Lim 85000']:
        cut_off_xlim = 85000
    else:
        cut_off_xlim = 75000
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/PBMCs/"+Image_ID_prefix+"-xlim.png"
    outline = ag.gateThreshold(fcs, name="outline", xCol=FSC, yCol=SSC, thresh=cut_off_xlim,parentGate=PBMC_step2,orientation='vertical', population="upper",filePlot=fileName)
    # // different cut_off_ylim value
    if ID in corrections_dict_raw['Y_Lim 70000']:
        cut_off_ylim = 70000
    elif ID in corrections_dict_raw['Y_Lim 90000']:
        cut_off_ylim = 90000
    elif ID in corrections_dict_raw['Y_Lim 100000']:
        cut_off_ylim = 100000
    elif ID in corrections_dict_raw['Y_Lim 120000']:
        cut_off_ylim = 120000
    else:
        cut_off_ylim = 110000
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/PBMCs/"+Image_ID_prefix+"-ylim.png"
    outline_1 = ag.gateThreshold(fcs, name="outline_1", xCol=FSC, yCol=SSC, thresh=cut_off_ylim,parentGate=outline,orientation='horizontal', population="lower",filePlot=fileName)
    # // save the image,PBMCs.png
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/PBMCs/"+Image_ID_prefix+"-PBMCs.png"
    PBMC=ag.gatePC(fcs,name="PBMC",xCol=FSC,yCol=SSC,center='centroid', adjustAngle=2, widthScale=8, heightScale = 6, parentGate=outline_1,filePlot=fileName)
    # // save the image,PBMCs_backgate.png
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/PBMCs/"+Image_ID_prefix+"-PBMCs_backgate.png"
    ag.backGate(fcs,population=PBMC,background_population=None,xCol=FSC, yCol=SSC,markersize=0.1, filePlot=fileName)
    fcs.update(ag.AGgate(PBMC,None,FSC,SSC,"PBMC"), QC=True, xlim=[0,215000], ylim=[0,215000],MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    
# ********************************************************************************************************************************************************************************************************
    # Next: 2)  Singlet gating
    # // setting the adjustAngle
    if ID in corrections_dict_raw['adjustAngle 4']:
        adjustAngle_1 = 4
    elif ID in corrections_dict_raw['adjustAngle 12']:
        adjustAngle_1 = 12
    else:
        adjustAngle_1 = 1
    # // setting the widthScale 
    if ID in corrections_dict_raw['widthScale 10']:
        widthScale_1 = 10
    else:
        widthScale_1 = 8
    # // setting the heightScale
    if ID in corrections_dict_raw['heightScale 2']:
        heightScale_1 = 2
    elif ID in corrections_dict_raw['heightScale 3']:
        heightScale_1 = 3
    else:
        heightScale_1 = 4 
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/singlets/"+Image_ID_prefix+"-singlets.png"
    singlets=ag.gatePC(fcs,xCol=FSC, yCol="FSC 488/10-H",name="singlets",center='density', adjustAngle= adjustAngle_1, widthScale=widthScale_1, heightScale=heightScale_1,parentGate=PBMC,filePlot=fileName)
    fcs.update(ag.AGgate(singlets,PBMC,FSC,"FSC 488/10-H","singlets"), QC=True,MFI=True, MFI_type='current', extra_MFI=extra_MFI_list) 

# ********************************************************************************************************************************************************************************************************
    # Next: 3)  CD3pos/neg gating
    # // setting the interval for CD3 
    if ID in corrections_dict_raw['CD3_interval_L']: 
        interval_1 = [800,4000]
    else:
        interval_1 = [300, 2000]
    xlim = ag.valleySeek(fcs=fcs, xCol=CD3, parentGate=singlets, interval=interval_1, require_local_min=True,
                        scale='bilog', T=200)
    if xlim == ag.np.inf:
        xlim = 400
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD3/"+Image_ID_prefix+"-CD3pos.png"
    CD3pos = ag.gateThreshold(fcs=fcs, name='CD3pos', xCol=CD3, yCol="FSC 488/10-H", parentGate=singlets, thresh=xlim, population='upper', scale='bilog', T=200, filePlot=fileName)
    CD3neg = ag.gateThreshold(fcs=fcs, name='CD3neg', xCol=CD3, yCol="FSC 488/10-H", parentGate=singlets, thresh=xlim,population='lower', scale='bilog', T=200)
    fcs.update(ag.AGgate(CD3pos, singlets, xCol=CD3, yCol="FSC 488/10-H", name="CD3pos"), QC=True, MFI=True,MFI_type='current', extra_MFI=extra_MFI_list,scale='bilog', T=200)
    fcs.update(ag.AGgate(CD3neg, singlets, xCol=CD3, yCol="FSC 488/10-H", name="CD3neg"), QC=True, MFI=True,MFI_type='current', extra_MFI=extra_MFI_list,scale='bilog', T=200)

    # Next:  3.0.1)  CD39pos in CD3pos/neg 
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD3pos_CD39pos/"+Image_ID_prefix+"-CD3posCD39pos.png"
    CD3pos_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD3pos_CD39pos", thresh=1000, parentGate=CD3pos,
                                     scale='bilog', T=1000, population='upper', update=False, filePlot=fileName)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD3neg_CD39pos/"+Image_ID_prefix+"-CD3negCD39pos.png"
    CD3neg_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD3negCD39pos", thresh=1000, parentGate=CD3neg,
                                     scale='bilog', T=1000, population='upper', update=False, filePlot=fileName)

    fcs.update(ag.AGgate(CD3pos_CD39pos, CD3pos, xCol=CD39, yCol="FSC 488/10-H", name="CD3pos_CD39pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD3neg_CD39pos, CD3neg, xCol=CD39, yCol="FSC 488/10-H", name="CD3neg_CD39pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)



    # Next:  3.0.2)  HLADR in CD3pos/neg
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD3pos_HLADRpos/"+Image_ID_prefix+"-CD3pos_HLADRpos.png"
    CD3pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="CD3pos_HLADRpos", parentGate=CD3pos, thresh=400,
                                       population='upper', orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD3neg_HLADRpos/"+Image_ID_prefix+"-CD3neg_HLADRpos.png"
    CD3neg_HLADRpos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD3neg_HLADRpos", thresh=400, parentGate=CD3neg,
                                     scale='bilog', T=400, population='upper', update=False, filePlot=fileName)

    fcs.update(ag.AGgate(CD3pos_HLADRpos, CD3pos, xCol=HLADR, yCol="FSC 488/10-H", name="CD3pos_HLADRpos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD3neg_HLADRpos, CD3neg, xCol=HLADR, yCol="FSC 488/10-H", name="CD3neg_HLADRpos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)

# ********************************************************************************************************************************************************************************************************
    # Next: 4)  CD4/CD8 gating
    # // Setting the interval of the CD4 
    if ID in corrections_dict_raw['CD4_interval [1000, 10000]']:
        interval_2 = [1000,10000]
    elif ID in corrections_dict_raw['CD4_interval [0, 3000]']:
        interval_2 = [0, 3000]
    elif ID in corrections_dict_raw['CD4_interval [2000, 10000]']:
        interval_2 = [2000, 10000]
    elif ID in corrections_dict_raw['CD4_interval [3000, 10000]']:
        interval_2 = [3000, 10000]
    elif ID in corrections_dict_raw['CD4_interval [5000, 10000]']:
        interval_2 = [5000, 10000]  
    else:
        interval_2 = [0, 10000]
    xlim = ag.valleySeek(fcs, xCol=CD8, parentGate=CD3pos, interval=interval_2, require_local_min=True, scale='bilog',T=200)
    if xlim == ag.np.inf:
        xlim = 1500
    CD8neg_tmp = ag.gateThreshold(fcs=fcs, name='CD8neg_temp', xCol=CD8, parentGate=CD3pos, thresh=xlim,population='lower', scale='bilog', T=200)
    # // Setting the interval of the CD8
    if ID in corrections_dict_raw['CD8_interval [0, 1000]']:
        interval_3 = [0, 1000]
    elif ID in corrections_dict_raw['CD8_interval [200, 1000]']:
        interval_3 = [200, 1000]
    elif ID in corrections_dict_raw['CD8_interval [300, 1000]']:
        interval_3 =  [300, 1000]
    elif ID in corrections_dict_raw['CD8_interval [450, 1000]']:
        interval_3 = [450, 1000]
    elif ID in corrections_dict_raw['CD8_interval [500, 1000]']:
        interval_3 = [500, 1000] 
    else:
        interval_3 = [50, 1000]
    ylim = ag.valleySeek(fcs, xCol=CD4, parentGate=CD8neg_tmp, interval=interval_3, require_local_min=True, scale='bilog', T=200)
    if ylim == ag.np.inf:
        ylim = 200
    if ylim < 100: 
        ylim = 200  
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4_CD8/"+Image_ID_prefix+"-CD4CD8quad.png"
    CD4pos, CD4CD8_doublepos, CD8pos, CD4CD8_doubleneg = ag.quadGate(fcs, names=["CD4pos", "CD4posCD8pos", "CD8pos","CD4negCD8neg"], xCol=CD8, yCol=CD4,
                                                                    parentGate=CD3pos, xThresh=xlim, yThresh=ylim,scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(CD4pos, CD3pos, xCol=CD8, yCol=CD4, name="CD4pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list, xlim=[-5000, 25000], ylim=[0, 1000], scale='bilog', T=200)
    fcs.update(ag.AGgate(CD8pos, CD3pos, xCol=CD8, yCol=CD4, name="CD8pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doublepos, CD3pos, xCol=CD8, yCol=CD4, name="CD4posCD8pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doubleneg, CD3pos, xCol=CD8, yCol=CD4, name="CD4negCD8neg"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:  4.0.1)  CD39pos in CD4/CD8 
    # // CD39pos in these four populations: CD4pos, CD4CD8_doublepos, CD8pos, CD4CD8_doubleneg
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_CD39pos/" + Image_ID_prefix + "-CD4posCD39pos.png"
    CD4pos_CD39pos = ag.gateThreshold(fcs, xCol=CD39,yCol="FSC 488/10-H", name="CD4pos_CD39pos", thresh=1000, parentGate=CD4pos, scale='bilog',T=1000, population='upper', update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_CD39pos/" + Image_ID_prefix + "-CD8posCD39pos.png"
    CD8pos_CD39pos = ag.gateThreshold(fcs, xCol=CD39,yCol="FSC 488/10-H", name="CD8pos_CD39pos", thresh=1000, parentGate=CD8pos, scale='bilog',T=1000, population='upper', update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_CD39pos/" + Image_ID_prefix + "-CD4CD8doubleposCD39pos.png"
    CD4CD8_doublepos_CD39pos = ag.gateThreshold(fcs, xCol=CD39,yCol="FSC 488/10-H", name="CD4CD8_doublepos_CD39pos", thresh=1000,parentGate=CD4CD8_doublepos, scale='bilog', T=1000, population='upper', update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_CD39pos/" + Image_ID_prefix + "-CD4CD8doublenegCD39pos.png"
    CD4CD8_doubleneg_CD39pos = ag.gateThreshold(fcs, xCol=CD39,yCol="FSC 488/10-H", name="CD4CD8_doubleneg_CD39pos", thresh=1000,parentGate=CD4CD8_doubleneg, scale='bilog', T=1000, population='upper',update=False, filePlot=fileName)

    fcs.update(ag.AGgate(CD4pos_CD39pos, CD4pos, xCol=CD39, yCol="FSC 488/10-H", name="CD4pos_CD39pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD8pos_CD39pos, CD8pos, xCol=CD39, yCol="FSC 488/10-H", name="CD8pos_CD39pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doublepos_CD39pos, CD4CD8_doublepos, xCol=CD39, yCol="FSC 488/10-H", name="CD4CD8_doublepos_CD39pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doubleneg_CD39pos, CD4CD8_doubleneg, xCol=CD39, yCol="FSC 488/10-H", name="CD4CD8_doubleneg_CD39pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:  4.0.2)  HLADR in CD4/CD8 
    # // HLADR in these four populations: CD4pos, CD4CD8_doublepos, CD8pos, CD4CD8_doubleneg
    
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_HLADRpos/" + Image_ID_prefix + "-CD4posHLADRpos.png"
    CD4pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="CD4pos_HLADRpos", parentGate=CD4pos, thresh=400,
                                       population='upper', orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_HLADRpos/" + Image_ID_prefix + "-CD4CD8_doublepos_HLADRpos.png"
    CD4CD8_doublepos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR,yCol="FSC 488/10-H", name="CD4CD8_doublepos_HLADRpos",parentGate=CD4CD8_doublepos, thresh=400, population='upper',
                                                 orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_HLADRpos/" + Image_ID_prefix  + "-CD8posHLADRpos.png"
    CD8pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR,yCol="FSC 488/10-H", name="CD8pos_HLADRpos", parentGate=CD8pos, thresh=400, population='upper', 
                                       orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4_CD8_HLADRpos/" + Image_ID_prefix  + "-CD4CD8_doubleneg_HLADRpos.png"
    CD4CD8_doubleneg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR,yCol="FSC 488/10-H",name="CD4CD8_doubleneg_HLADRpos",parentGate=CD4CD8_doubleneg, thresh=400, population='upper',
                                                 orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    # update the HLADR
    fcs.update(ag.AGgate(CD4pos_HLADRpos, CD4pos, xCol=HLADR, yCol="FSC 488/10-H", name="CD4pos_HLADRpos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doublepos_HLADRpos, CD4CD8_doublepos, xCol=HLADR, yCol="FSC 488/10-H", name="CD4CD8_doublepos_HLADRpos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD8pos_HLADRpos, CD8pos, xCol=HLADR, yCol="FSC 488/10-H", name="CD8pos_HLADRpos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD4CD8_doubleneg_HLADRpos, CD4CD8_doubleneg, xCol=HLADR, yCol="FSC 488/10-H", name="CD4CD8_doubleneg_HLADRpos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
# ********************************************************************************************************************************************************************************************************
    # Next: 5)  CD8pos_N_CM_E_EM
    # the naive, the central memory, the memoryeffector, the effector
    CD8_ylim = [50,400]
    output = para_extracter_1("CD8_ylim", r'\[(\d+),\s*(\d+)\]')
    if output:
        CD8_ylim = output
    ylim = ag.valleySeek(fcs, xCol=CD45RA, interval=CD8_ylim, require_local_min=True, parentGate=CD8pos, scale='bilog', T=200)
    if ylim == ag.np.inf:
        ylim=200
    if ylim < 200:
        ylim=200
    # // top_thresh
    if ID in corrections_dict_raw['Top_thresh 1000']:
        top_thresh = 1000
    else:
        top_thresh = 500
    tempCD8quad_top = ag.gateThreshold(fcs,xCol=CCR7, yCol=CD45RA,name="topTmp",parentGate=CD8pos, thresh =top_thresh, population='upper',orientation='horizontal',scale='bilog', T=200,update=False)
    # // CD8_xlim
    CD8_xlim = [100,800]
    output = para_extracter_1("CD8_xlim", r'\[(\d+),\s*(\d+)\]')
    if output:
        CD8_xlim = output
    mean, median, sigma, maxVal = ag.axisStats(fcsDF = fcs(), xCol = CCR7,vI = tempCD8quad_top(), scale = 'linear', T = 2000)
    
    xlim1 = ag.valleySeek(fcs,xCol=CCR7,interval=CD8_xlim,require_local_min=True, parentGate=tempCD8quad_top,scale='bilog', T=1000)
    if xlim1 == ag.np.inf:
        xlim1=600
    xlim = ag.densityDelimitation(fcs, xCol=CCR7, parentGate=tempCD8quad_top, interval=[maxVal,median+2*sigma], sigma=2, bins=300, limit_threshold=0.08, direction='right',scale='bilog', T=1000)
    if xlim1 < xlim:
        xlim = xlim1
    if ID in corrections_dict_raw['CD8 second correction']['xlim+50']:
        xlim +=50
    elif ID in corrections_dict_raw['CD8 second correction']['xlim=250']:
        xlim = 250
    elif ID in corrections_dict_raw['CD8 second correction']['xlim=300']:
        xlim = 300
    elif ID in corrections_dict_raw['CD8 second correction']['xlim=400']:
        xlim = 400
    elif ID in corrections_dict_raw['CD8 second correction']['xlim-50']:
        xlim -=50
    elif ID in corrections_dict_raw['CD8 second correction']['xlim-100']:
        xlim -=100
    elif ID in corrections_dict_raw['CD8 second correction']['xlim-150']:
        xlim -=150
    elif ID in corrections_dict_raw['CD8 second correction']['xlim-200']:
        xlim -=200
    elif ID in corrections_dict_raw['CD8 second correction']['xlim-300']:
        xlim -=300
    

    ################# add new correction here, move the xlim a little bit to the left, make the seperation based on the CD45RAhi part ############################
    if ID in corrections_dict_raw['CCR7_move_to_left']:
        CD45RAhi = ag.gateThreshold(fcs,xCol=CD45RA,name="CD45RAhi",parentGate=CD8pos,thresh=8000, population='upper',scale='bilog',T=1000,update=False)
        xlim = ag.valleySeek(fcs,xCol=CCR7,interval=[100,500],require_local_min=True, parentGate=CD45RAhi,scale='bilog', T=1000)
        xlim = xlim-100
    ##############################################################################################################################################################
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_N_CM_E_EM/"+Image_ID_prefix+"-CD8pos_N_CM_E_EM.png"
    effector_CD8pos, naive_CD8pos, central_memory_CD8pos, effector_memory_CD8pos = ag.quadGate(fcs=fcs, names=['effector_CD8pos', 'naive_CD8pos', 'central_memory_CD8pos', 'effector_memory_CD8pos'], 
                    xCol=CCR7, yCol=CD45RA, xThresh=xlim, yThresh=ylim, parentGate=CD8pos,scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(effector_CD8pos, CD8pos, xCol=CCR7, yCol=CD45RA, name="effector_CD8pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list,xlim=[0, 3000], ylim=[-200, 1000],scale='bilog', T=200)
    fcs.update(ag.AGgate(naive_CD8pos, CD8pos, xCol=CCR7, yCol=CD45RA, name="naive_CD8pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(central_memory_CD8pos, CD8pos, xCol=CCR7, yCol=CD45RA, name="central_memory_CD8pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(effector_memory_CD8pos, CD8pos, xCol=CCR7, yCol=CD45RA, name="effector_memory_CD8pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)

    # Next:  5.0.1)  CD39pos in CD8pos_N_CM_E_EM
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_CD39pos/"+Image_ID_prefix+"-CD39effector_CD8pos.png"    
    CD39effector_CD8pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39effector_CD8pos", thresh=1000,parentGate=effector_CD8pos, scale='bilog', T=1000, population='upper',update=False, filePlot=fileName)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_CD39pos/"+Image_ID_prefix+"-CD39naive_CD8pos.png"
    CD39naive_CD8pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39naive_CD8pos", thresh=1000, parentGate=naive_CD8pos,scale='bilog', T=1000, population='upper', update=False, filePlot=fileName)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_CD39pos/"+Image_ID_prefix+"-CD39effector_memory_CD8pos.png"
    CD39central_memory_CD8pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39central_memory_CD8pos", thresh=1000,parentGate=central_memory_CD8pos, scale='bilog', T=1000,population='upper',
                                                  update=False, filePlot=fileName)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_CD39pos/"+Image_ID_prefix+"-CD39effector_memory_CD8pos.png"
    CD39effector_memory_CD8pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39effector_memory_CD8pos", thresh=1000,parentGate=effector_memory_CD8pos, scale='bilog', 
                                                  T=1000,population='upper', update=False, filePlot=fileName)

    fcs.update(ag.AGgate(CD39effector_CD8pos, effector_CD8pos, xCol=CD39, yCol="FSC 488/10-H", name="CD39effector_CD8pos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39naive_CD8pos, naive_CD8pos, xCol=CD39, yCol="FSC 488/10-H", name="CD39naive_CD8pos"), QC=True,MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39central_memory_CD8pos, central_memory_CD8pos, xCol=CD39, yCol="FSC 488/10-H",name="CD39central_memory_CD8pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39effector_memory_CD8pos, effector_memory_CD8pos, xCol=CD39, yCol="FSC 488/10-H",name="CD39effector_memory_CD8pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)

    # Next:  5.0.2)  HLADR in CD8pos_N_CM_E_EM
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_HLADRpos/"+Image_ID_prefix+ "-effector_CD8pos_HLADRpos.png"
    effector_CD8pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, name="effector_CD8pos_HLADRpos",parentGate=effector_CD8pos, thresh=400, population='upper',
                                                orientation='vertical', scale='bilog', T=400, update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_HLADRpos/"+Image_ID_prefix + "-naive_CD8pos_HLADRpos.png"
    naive_CD8pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, name="naive_CD8pos_HLADRpos", parentGate=naive_CD8pos, thresh=400, population='upper', 
                                                orientation='vertical', scale='bilog',T=400, update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_HLADRpos/"+Image_ID_prefix + "-central_memory_CD8pos_HLADRpos.png"
    central_memory_CD8pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, name="central_memory_CD8pos_HLADRpos", parentGate=central_memory_CD8pos, thresh=400, population='upper',
                                                orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD8pos_HLADRpos/"+Image_ID_prefix + "-effector_memory_CD8pos_HLADRpos.png"
    effector_memory_CD8pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, name="effector_memory_CD8pos_HLADRpos",parentGate=effector_memory_CD8pos, thresh=400,
                                                       population='upper', orientation='vertical', scale='bilog', T=400,update=False, filePlot=fileName)
    fcs.update(ag.AGgate(effector_CD8pos_HLADRpos, effector_CD8pos, xCol=HLADR, yCol="FSC 488/10-H", name="effector_CD8pos_HLADRpos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(naive_CD8pos_HLADRpos, naive_CD8pos, xCol=HLADR, yCol="FSC 488/10-H", name="naive_CD8pos_HLADRpos"),QC=True, MFI=True, MFI_type='current',extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(central_memory_CD8pos_HLADRpos, central_memory_CD8pos, xCol=HLADR, yCol="FSC 488/10-H",name="central_memory_CD8pos_HLADRpos"), QC=True, MFI=True, MFI_type='current',extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(effector_memory_CD8pos_HLADRpos, effector_memory_CD8pos, xCol=HLADR, yCol="FSC 488/10-H",name="effector_memory_CD8pos_HLADRpos"), QC=True, MFI=True, MFI_type='current',extra_MFI=extra_MFI_list)
 # ********************************************************************************************************************************************************************************************************
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
    # Next:  6.1)  CD4+CXCR5+ gating
    # gating CD4+CXCR5+ CELLS USING "ylim" from above
    tmpHiCD45RA = ag.gateThreshold(fcs, xCol=CD45RA, name="tmpHiCD45RA", parentGate=CD4pos, thresh=ylim,population='upper', orientation='vertical', scale='bilog', T=200, update=False)
    # // setting the CXCR5_lim_threshold
    if ID in corrections_dict_raw['CXCR5_lim_threshold 0.4']:
        CXCR5_lim_threshold = 0.4
    elif ID in corrections_dict_raw['CXCR5_lim_threshold 0.3']:
        CXCR5_lim_threshold = 0.3
    else:
        CXCR5_lim_threshold = 0.2
    CXCR5_lim = ag.densityDelimitation(fcs, xCol=CXCR5, parentGate=tmpHiCD45RA, interval=[0, 2500], direction='right',limit_threshold=CXCR5_lim_threshold, scale='bilog', T=200)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CXCR5pos/"+Image_ID_prefix+"-CD4posCXCR5pos.png"
    CD4posCXCR5pos = ag.gateCorner(fcs, name='CD4posCXCR5pos', xCol=CD45RA, yCol=CXCR5, parentGate=CD4pos, xThresh=ylim, yThresh=CXCR5_lim, 
                             xOrientation='lower', yOrientation='upper', scale='bilog', T=200,filePlot=fileName)
    fcs.update(ag.AGgate(CD4posCXCR5pos, CD4pos, xCol=CD45RA, yCol=CXCR5, name="CD4posCXCR5pos"), QC=True, MFI=True, MFI_type='current',
               extra_MFI=extra_MFI_list, scale='bilog', T=200)
    # Next:   6.1.1)  CD39pos in CD4+CXCR5+
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/CD4posCXCR5posCD39pos/" + Image_ID_prefix + "-CD4posCXCR5posCD39pos.png"
    CD4posCXCR5posCD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD4posCXCR5posCD39pos", thresh=1000, parentGate=CD4posCXCR5pos, scale='bilog',T=1000, population='upper', update=False, filePlot=fileName)
    fcs.update(ag.AGgate(CD4posCXCR5posCD39pos,CD4posCXCR5pos, xCol=CD39, yCol="FSC 488/10-H", name="CD4posCXCR5posCD39pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:   6.1.2)  HLADR in CD4+CXCR5+
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4posCXCR5posHLADRpos/"+Image_ID_prefix+ "-CD4posCXCR5posHLADRpos.png"
    CD4posCXCR5posHLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="CD4posCXCR5posHLADRpos",parentGate=CD4posCXCR5pos, thresh=400, population='upper',
                                                orientation='vertical', scale='bilog', T=400, update=False, filePlot=fileName)
    fcs.update(ag.AGgate(CD4posCXCR5posHLADRpos, CD4posCXCR5pos, xCol=HLADR, yCol="FSC 488/10-H", name="CD4posCXCR5posHLADRpos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)

    # Next:  6.2)  CD4+CXCR5+ Tfh gating
    # gating follicular, small gate mainly to remove debris, gating Follicular_Th from CD4posCXCR5pos cells
    CXCR5pos_xlim = ag.densityDelimitation(fcs, xCol=CD25, parentGate=CD4posCXCR5pos, interval=[0, 3000], direction='right',limit_threshold=0.05, scale='bilog', T=200)
    if CXCR5pos_xlim == ag.np.inf:
        CXCR5pos_xlim = 500
    CXCR5pos_ylim = ag.densityDelimitation(fcs, xCol=CD127, parentGate=CD4posCXCR5pos, interval=[0, 2000], direction='left',limit_threshold=0.1, scale='bilog', T=200)
    if CXCR5pos_ylim == ag.np.inf:
        CXCR5pos_ylim = 300
    # Default xlim_bottom
    CXCR5pos_xlim_bottom = 300
    CXCR5pos_ylim_top = 1000
    # Custom CXCR5pos_xlim_bottom, CXCR5pos_ylim, CXCR5pos_xlim, CXCR5pos_ylim_top based on the yaml file
    # // in_1 is the increase or decrease value for CXCR5pos_xlim_bottom
    in_1 = para_extracter('in_1',r'-?\d+\.?\d*$')
    CXCR5pos_xlim_bottom = CXCR5pos_xlim_bottom + in_1
    # //n in_2 is the increase or decrease value for CXCR5pos_ylim 
    in_2 = para_extracter('in_2',r'-?\d+\.?\d*$')
    CXCR5pos_ylim = CXCR5pos_ylim + in_2
    # // in_3 is the increase or decrease value for CXCR5pos_xlim
    in_3 = para_extracter('in_3',r'-?\d+\.?\d*$')
    CXCR5pos_xlim = CXCR5pos_xlim + in_3
    # // in_4 is the increase or decrease value for XCR5pos_ylim_top
    in_4 = para_extracter('in_4',r'-?\d+\.?\d*$')
    CXCR5pos_ylim_top = CXCR5pos_ylim_top + in_4
    # Follicular T helper (Tfh)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4posCXCR5pos_Tfh/"+Image_ID_prefix+"-CD4posCXCR5pos_Tfh.png"
    # here, the xParam is the same as the yParam,meaning the seperating curve is a line
    CD4posCXCR5pos_Tfh = ag.gateBezier(fcs, name='CD4posCXCR5pos_Tfh', xCol=CD25, yCol=CD127, parentGate=CD4posCXCR5pos,endExtension='up', startExtension='left',
                                  points=[(CXCR5pos_xlim_bottom, CXCR5pos_ylim), (CXCR5pos_xlim, CXCR5pos_ylim_top)], population='upper', xParam=[0.8],yParam=[0.8], scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(CD4posCXCR5pos_Tfh, CD4posCXCR5pos, xCol=CD25, yCol=CD127, name="CD4posCXCR5pos_Tfh"), QC=True, MFI=True,
               MFI_type='current', extra_MFI=extra_MFI_list, xlim=[0, 2000], ylim=[0, 2000], scale='bilog',T=200)

    # Next:   6.2.1)  CD39pos in CD4+CXCR5+ Tfh
    fileName = ag.AGConfig.ag_home + "/plots/phase_III/ImmSight/T_Panel/Tfh_CD39pos/" + Image_ID_prefix + "-Tfh_CD39pos.png"
    Tfh_CD39pos = ag.gateThreshold(fcs, xCol=CD39,yCol="FSC 488/10-H", name="Tfh_CD39pos", thresh=1000, parentGate=CD4posCXCR5pos_Tfh, scale='bilog',T=1000, population='upper', update=False, filePlot=fileName)
    fcs.update(ag.AGgate(Tfh_CD39pos, CD4posCXCR5pos_Tfh, xCol=CD39, yCol="FSC 488/10-H", name="Tfh_CD39pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:   6.2.2)  HLADR in CD4+CXCR5+ Tfh
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Tfh_HLADRpos/"+Image_ID_prefix+ "-Tfh_HLADRpos.png"
    Tfh_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="Tfh_HLADRpos",parentGate=CD4posCXCR5pos_Tfh, thresh=400, population='upper',
                                                orientation='vertical', scale='bilog', T=400, update=False, filePlot=fileName)
    fcs.update(ag.AGgate(Tfh_HLADRpos, CD4posCXCR5pos_Tfh, xCol=HLADR, yCol="FSC 488/10-H", name="Tfh_HLADRpos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)

# ********************************************************************************************************************************************************************************************************
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
    effector_CD4pos, naive_CD4pos, central_memory_CD4pos, effector_memory_CD4pos = ag.quadGate(fcs=fcs, names=['effector_CD4pos', 'naive_CD4pos', 'central_memory_CD4pos', 'effector_memory_CD4pos'], 
                                                                                            xCol=CCR7, yCol=CD45RA, xThresh=xlim, yThresh=ylim, parentGate=CD4pos, 
                                                                                            scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(effector_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="effector_CD4pos"), QC=True, MFI=True,MFI_type='current', extra_MFI=extra_MFI_list, 
                        xlim=[0, 3000], ylim=[-200, 1000],scale='bilog', T=200)
    fcs.update(ag.AGgate(naive_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="naive_CD4pos"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(central_memory_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="central_memory_CD4pos"), QC=True,MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(effector_memory_CD4pos, CD4pos, xCol=CCR7, yCol=CD45RA, name="effector_memory_CD4pos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)

    # Next: 6.0.1)  CD39pos in CD4pos-N_CM_E_EM
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CD39pos/"+Image_ID_prefix+ "-effector_CD4pos_CD39pos.png"
    CD39effector_CD4pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD39effector_CD4pos", thresh=1000,parentGate=effector_CD4pos, scale='bilog', T=1000,
                                            population='upper',update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CD39pos/"+Image_ID_prefix+ "-naive_CD4pos_CD39pos.png"
    CD39naive_CD4pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39naive_CD4pos", thresh=1000, parentGate=naive_CD4pos,scale='bilog', T=1000, population='upper', 
                                        update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CD39pos/"+Image_ID_prefix+ "-central_memory_CD4pos_CD39pos.png"
    CD39central_memory_CD4pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39central_memory_CD4pos", thresh=1000,parentGate=central_memory_CD4pos, scale='bilog', 
                                                 T=1000,population='upper', update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CD39pos/"+Image_ID_prefix+ "-effector_memory_CD4pos_CD39pos.png"
    CD39effector_memory_CD4pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="CD39effector_memory_CD4pos", thresh=1000,parentGate=effector_memory_CD4pos, scale='bilog',
                                                   T=1000,population='upper', update=False, filePlot=fileName)

    fcs.update(ag.AGgate(CD39effector_CD4pos, effector_CD4pos, xCol=CD39, yCol="FSC 488/10-H", name="CD39effector_CD4pos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39naive_CD4pos, naive_CD4pos, xCol=CD39, yCol="FSC 488/10-H", name="CD39naive_CD4pos"), QC=True,MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39central_memory_CD4pos, central_memory_CD4pos, xCol=CD39, yCol="FSC 488/10-H",name="CD39central_memory_CD4pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39effector_memory_CD4pos, effector_memory_CD4pos, xCol=CD39, yCol="FSC 488/10-H",name="CD39effector_memory_CD4pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)

    # Next: 6.0.2)  HLADR in CD4pos-N_CM_E_EM
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_HLADRpos/"+Image_ID_prefix+ "-effector_CD4pos_HLADRpos.png"
    effector_CD4pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="effector_CD4pos_HLADRpos",parentGate=effector_CD4pos, thresh=400, population='upper',
                                                orientation='vertical', scale='bilog', T=400, update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_HLADRpos/"+Image_ID_prefix + "-naive_CD4pos_HLADRpos.png"
    naive_CD4pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="naive_CD4pos_HLADRpos", parentGate=naive_CD8pos, thresh=400, population='upper', 
                                                orientation='vertical', scale='bilog',T=400, update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_HLADRpos/"+Image_ID_prefix + "-central_memory_CD4pos_HLADRpos.png"
    central_memory_CD4pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="central_memory_CD4pos_HLADRpos", parentGate=central_memory_CD8pos, thresh=400, population='upper',
                                                orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_HLADRpos/"+Image_ID_prefix + "-effector_memory_CD4pos_HLADRpos.png"
    effector_memory_CD4pos_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="effector_memory_CD4pos_HLADRpos",parentGate=effector_memory_CD4pos, thresh=400,
                                                       population='upper', orientation='vertical', scale='bilog', T=400,update=False, filePlot=fileName)
    fcs.update(ag.AGgate(effector_CD4pos_HLADRpos, effector_CD4pos, xCol=HLADR, yCol="FSC 488/10-H", name="effector_CD4pos_HLADRpos"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(naive_CD4pos_HLADRpos, naive_CD4pos, xCol=HLADR, yCol="FSC 488/10-H", name="naive_CD4pos_HLADRpos"),QC=True, MFI=True, MFI_type='current',extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(central_memory_CD4pos_HLADRpos, central_memory_CD4pos, xCol=HLADR, yCol="FSC 488/10-H",name="central_memory_CD4pos_HLADRpos"), QC=True, MFI=True, MFI_type='current',extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(effector_memory_CD4pos_HLADRpos, effector_memory_CD4pos, xCol=HLADR, yCol="FSC 488/10-H",name="effector_memory_CD4pos_HLADRpos"), QC=True, MFI=True, MFI_type='current',extra_MFI=extra_MFI_list)
# ********************************************************************************************************************************************************************************************************
    # Next: 7)  CD4pos_Tregs
    # now, gating CD4pos_Treg based on the CD4 positive cells
    # there are two points, one is at lower-left (xmaxVal, ylim), the second one is at the top-right (xlim, ymaxVal)
    # if X param > Y param, the curve is like spoon
    # if X param < Y param, the curve is like the back of the spoon
    # if X param = Y param, then the curve is a direct line
    # the xlim is based on the CD25 density distribution under the parent gate of CD4pos
    CD4pos_hiCD25 = ag.gateThreshold(fcs, name="tmp", xCol=CD25, parentGate=CD4pos, thresh=-200, population='upper')
    mean, median, sigma, maxVal = ag.axisStats(fcsDF=fcs(), xCol=CD25, vI=CD4pos_hiCD25(), scale='bilog', T=200)
    lim_upper = ag.inverseTransformWrapper(maxVal + 0.5 * abs(sigma), scale='bilog', T=200)
    lim_lower = ag.inverseTransformWrapper(maxVal - 0.5 * abs(sigma), scale='bilog', T=200)
    # // setting the xmaxVal
    xmaxval_step = para_extracter('xparam_step',r'-?\d+\.?\d*$')
    if ID in corrections_dict_raw["xmaxval_step -500"]:
        xmaxval_step = -500
    elif ID in corrections_dict_raw["xmaxval_step 200"]:
        xmaxval_step = 200
    else:
        xmaxval_step = 0
    xmaxVal = lim_upper+ xmaxval_step
    # xmaxVal is pne of the two points, one is at lower-left (xmaxVal, ylim)
    tmp_maxCD25 = ag.gateThreshold(fcs, name="tmp", xCol=CD25, parentGate=CD4pos, thresh=lim_upper, population='lower')
    maxCD25 = ag.gateThreshold(fcs, name="tmp", xCol=CD25, parentGate=tmp_maxCD25, thresh=lim_lower, population='upper')
    # the ylim is based on the CD127 distribution under the parent gate of CD25
    ylim = ag.densityDelimitation(fcs=fcs, xCol=CD127, parentGate=maxCD25, interval=[-200, 1000], limit_threshold=0.05,direction='left', scale='bilog', T=200)
    # // setting the ylim
    ylim_step = para_extracter('ylim_step',r'-?\d+\.?\d*$')
    ylim = ylim+ylim_step
    # ylim is pne of the two points, one is at lower-left (xmaxVal, ylim)
    mean, median, sigma, maxVal = ag.axisStats(fcsDF=fcs(), xCol=CD127, vI=CD4pos(), scale='bilog', T=200)
    lim_upper = ag.inverseTransformWrapper(maxVal + 0.5 * abs(sigma), scale='bilog', T=200)
    lim_lower = ag.inverseTransformWrapper(maxVal - 0.5 * abs(sigma), scale='bilog', T=200)
    # // setting the ymaxVal
    ymaxval_step = para_extracter('ymaxval_step',r'-?\d+\.?\d*$')
    ymaxVal = ymaxval_step+lim_lower
    # ymaxVal is one of the two points, second one is at the top-right (xlim, ymaxVal)
    tmp_maxCD127 = ag.gateThreshold(fcs, name="tmp", xCol=CD127, parentGate=CD4pos, thresh=lim_upper,population='lower')
    maxCD127 = ag.gateThreshold(fcs, name="tmp", xCol=CD127, parentGate=tmp_maxCD127, thresh=lim_lower,population='upper')
    # the xlim is based on the CD25 distribution under the parent gate of CD127
    xlim = ag.densityDelimitation(fcs=fcs, xCol=CD25, parentGate=maxCD127, interval=[200, 10000], limit_threshold=0.05,
                                  direction='right', scale='bilog', T=200)
    # // setting the xlim
    xlim_step = para_extracter('xlim_step',r'-?\d+\.?\d*$')
    xlim = xlim + xlim_step
    # xlim is one of the two points, second one is at the top-right (xlim, ymaxVal)
    if ymaxVal < 0:
        ymaxVal = 0
    # here, the X param > Y param, the curve is like spoon
    bezierXParam = [0.7]
    bezierYParam = [0.3]
    # Take midpoint between ylim and ymaxVal (where ymaxval >> ylim) in transformed coordinates
    tyMid = (ag.transformWrapper([ymaxVal], T=200, scale='bilog')[0] + ag.transformWrapper([ylim], T=200, scale='bilog')[0]) / 2
    # Copy code from above to make a 'street' and calc density delim (we have already calced mean and sigma)
    lim_upper = ag.inverseTransformWrapper(tyMid + 0.5 * abs(sigma), scale='bilog', T=200)
    lim_lower = ag.inverseTransformWrapper(tyMid - 0.5 * abs(sigma), scale='bilog', T=200)
    tmp_midCD127 = ag.gateThreshold(fcs, name="tmp", xCol=CD127, parentGate=CD4pos, thresh=lim_upper,population='lower')
    midCD127 = ag.gateThreshold(fcs, name="tmp", xCol=CD127, parentGate=tmp_midCD127, thresh=lim_lower,population='upper')
    xMid = ag.densityDelimitation(fcs=fcs, xCol=CD25, parentGate=midCD127, interval=[-200, 10000], limit_threshold=0.3,
                                  direction='right', scale='bilog', T=200)
    # See which percent this xMid is at from relative to xlim and xmaxVal (xlim >> xmaxVal)
    fraction_traveled_x = (xMid - xmaxVal) / (xlim - xmaxVal)              
    if fraction_traveled_x < 0.30:
        bezierXParam = [0.3]
        bezierYParam = [0.7]
    elif fraction_traveled_x > 0.70:
        bezierXParam = [0.7]
        bezierYParam = [0.3]
    else:
        bezierXParam = [0.5]
        bezierYParam = [0.5]
    if ymaxVal < 0:
        ymaxVal = 0
    # // customize the xparam_step, then customzing the bezierXParam[0], bezierYParam[0]
    xparam_step = para_extracter('xparam_step', r'-?\d+\.?\d*$')
    bezierXParam[0] = bezierXParam[0]+xparam_step
    yparam_step = para_extracter('yparam_step', r'-?\d+\.?\d*$')
    bezierYParam[0] = bezierYParam[0]+yparam_step
    # bezierXParam[0] and bezierYParam[0] will determine the shape of the curve
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_Tregs/"+Image_ID_prefix+"-Tregs.png"
    Tregs = ag.gateBezier(fcs, name="Tregs", xCol=CD25, yCol=CD127, endExtension='right', startExtension='down', population='lower', parentGate=CD4pos, 
                          points=[(xmaxVal, ylim), (xlim, ymaxVal)], xParam=bezierXParam, yParam=bezierYParam, scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(Tregs, CD4pos, xCol=CD25, yCol=CD127, name="Tregs"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list, xlim=[0, 10000], ylim=[0, 2000], scale='bilog', T=200)
    # Next:  7.0.1)  CD39pos in Tregs
    tmp_lowCD39Tregs = ag.gateThreshold(fcs, name="tmp_lowCD39Tregs", xCol=CD39, parentGate=Tregs, thresh=1000, population='lower')
    mean, median, sigma, maxVal = ag.axisStats(fcsDF=fcs(), xCol=CD39, vI=tmp_lowCD39Tregs(), scale='bilog', T=200)
    lim_lower = ag.inverseTransformWrapper(mean, scale='bilog', T=200)
    Tregs_CD39pos_ylim_global = ag.valleySeek(fcs, xCol=CD39, interval=[lim_lower,1100], require_local_min=True, parentGate=Tregs, scale='bilog', T=1000)

    if ID in corrections_dict_raw["Tregs_CD39pos_ylim_global 400"]:
        Tregs_CD39pos_ylim_global = 400
    elif ID in corrections_dict_raw["Tregs_CD39pos_ylim_global 600"]:
        Tregs_CD39pos_ylim_global = 600
    elif ID in corrections_dict_raw["Tregs_CD39pos_ylim_global 700"]:
        Tregs_CD39pos_ylim_global = 700
    elif ID in corrections_dict_raw["Tregs_CD39pos_ylim_global 1000"]:
        Tregs_CD39pos_ylim_global = 1000 
    elif ID in corrections_dict_raw["Tregs_CD39pos_ylim_global 7000"]:
        Tregs_CD39pos_ylim_global = 7000

    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_Tregs_CD39/"+Image_ID_prefix+ "-CD39Tregs.png"
    CD39Tregs = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD39Tregs", thresh=Tregs_CD39pos_ylim_global, parentGate=Tregs, scale='bilog', T=1000,population='upper', update=False, filePlot=fileName)
    fcs.update(ag.AGgate(CD39Tregs, Tregs, xCol=CD39, yCol="FSC 488/10-H", name="CD39Tregs"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    # Next:  7.0.2)  HLADR in Tregs
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_Tregs_HLADR/"+Image_ID_prefix+ "-HLADRTregs.png"
    Tregs_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Tregs_HLADRpos", parentGate=Tregs, thresh=400, population='upper', orientation='vertical', scale='bilog',
                                       T=400, update=False,filePlot=fileName)
    fcs.update(ag.AGgate(Tregs_HLADRpos, Tregs, xCol=HLADR, yCol="FSC 488/10-H", name="Tregs_HLADRpos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    # Next:  7.1)  CD194+ Tregs
    CD194_Tregs_interval = [500, 5000]
    for keys, ids in corrections_dict_raw.items():
        if keys.startswith('CD194_Tregs_interval'):
            CD194_Tregs_match = re.search(r'\[(\d+),\s*(\d+)\]', keys)
            if CD194_Tregs_match and ID in ids:
                CD194_Tregs_interval = list(map(int,CD194_Tregs_match.groups()))
                break
    xlim = ag.valleySeek(fcs, xCol=CD194, parentGate=Tregs, interval=CD194_Tregs_interval, require_local_min=True, sigma=3, scale='bilog', T=200)
    if xlim == ag.np.inf:
        xlim = 1500
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194_Tregs/"+Image_ID_prefix+ "-CD194_Tregs.png"
    CD194_Tregs = ag.gateThreshold(fcs=fcs, name="CD194_Tregs", xCol=CD194, yCol=CD127, parentGate=Tregs, thresh=xlim,
                                   orientation='vertical', scale='bilog', T=200, population='upper', filePlot=fileName)
    fcs.update(ag.AGgate(CD194_Tregs, Tregs, xCol=CD194, yCol=CD127, name="CD194_Tregs"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list, 
               xlim=[0, 25000], ylim=[0, 700],scale='bilog', T=200)
    # Next:   7.1.1)  CD39pos in CD194+ Tregs
    tmp_lowCD39CD194posTregs = ag.gateThreshold(fcs, name="tmp_lowCD39CD194posTregs", xCol=CD39, parentGate=CD194_Tregs, thresh=1000, population='lower')
    mean, median, sigma, maxVal = ag.axisStats(fcsDF=fcs(), xCol=CD39, vI=tmp_lowCD39CD194posTregs(), scale='bilog', T=200)
    lim_lower = ag.inverseTransformWrapper(mean, scale='bilog', T=200)
    CD194_Tregs_CD39pos_ylim_global = ag.valleySeek(fcs, xCol=CD39, interval=[lim_lower,1100], require_local_min=True, parentGate=CD194_Tregs, scale='bilog', T=1000)
    if ID in corrections_dict_raw["CD194_Tregs_CD39pos_ylim_global 700"]:
        CD194_Tregs_CD39pos_ylim_global = 700
    elif ID in corrections_dict_raw["CD194_Tregs_CD39pos_ylim_global 7000"]:
        CD194_Tregs_CD39pos_ylim_global = 7000
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CD194_Tregs_CD39/"+Image_ID_prefix+ "-CD39posCD194posTregs.png"
    CD39posCD194posTregs = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD39posCD194posTregs", thresh=CD194_Tregs_CD39pos_ylim_global, parentGate=CD194_Tregs, scale='bilog', T=1000,population='upper', update=False, filePlot=fileName)
    fcs.update(ag.AGgate(CD39posCD194posTregs, CD194_Tregs, xCol=CD39, yCol="FSC 488/10-H", name="CD39posCD194posTregs"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    # Next:   7.1.2)  HLADR in CD194+ Tregs
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD4pos_CD194_Tregs_HLADR/"+Image_ID_prefix+ "-HLADRposCD194posTregs.png"
    HLADRposCD194posTregs = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="HLADRposCD194posTregs", parentGate=CD194_Tregs, thresh=400, population='upper', orientation='vertical', scale='bilog', 
                                             T=400, update=False,filePlot=fileName)
    fcs.update(ag.AGgate(HLADRposCD194posTregs, CD194_Tregs, xCol=HLADR, yCol="FSC 488/10-H", name="HLADRposCD194posTregs"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
# ********************************************************************************************************************************************************************************************************
    # Next:  7.2)  CD4pos_Tregs subpopulations 
    # activated, secreting, resting
    # First gate out the resting (top) by estimating half-normal distribution on bottom cluster.
    # Then look at the CD25 expression of these resting Tregs, and set the limit between secreting and activated Tregs
    # at some point of the right-tail of the CD25 distribution
    # First do arbitrary cut in the y-axis somewhere between 0-1000
    ylim = ag.valleySeek(fcs, xCol=CD45RA, parentGate=Tregs, interval=[0,1000], scale='bilog', T=200)
    tmpbottomTreg = ag.gateThreshold(fcs, xCol=CD25, yCol=CD45RA, name="bottomTreg", parentGate=Tregs, thresh=ylim,population='lower', orientation='horizontal', scale='bilog', T=200, update=False)
    # In this population, find the most dense point
    tmp, tmp1, tmp2, maxVal = ag.axisStats(fcs(), xCol=CD45RA, vI=tmpbottomTreg(), scale='bilog', T=200)
    maxVal = ag.inverseTransformWrapper(maxVal, scale='bilog', T=200)
    mean, sigma = ag.halfNormalDistribution(fcs, xCol=CD45RA, mean=maxVal, direction='left', parentGate=tmpbottomTreg,
                                            scale='bilog', T=200)
    # Now we can set the final ylim with big margin (4 std devs):
    ylim = ag.inverseTransformWrapper(mean + 4 * abs(sigma), scale='bilog', T=200)
    if ylim > 1000:
        ylim = 1000
    if ylim < 100:
        ylim = 100

    if ID in corrections_dict_raw['Tregs_ylim 1000']:
        ylim = 1000
    elif ID in corrections_dict_raw['Tregs_ylim 700']:
        ylim = 700
    # We use this to gate the activated Tregs (Dont save image tho, unneccesary)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations/"+Image_ID_prefix+"-restingTreg.png"
    restingTreg = ag.gateThreshold(fcs, xCol=CD25, yCol=CD45RA, name="restingTreg", thresh=ylim, parentGate=Tregs,
                                   scale='bilog', T=200, population='upper', orientation='horizontal', update=False,
                                   filePlot=fileName)
    # ********************************************************************************************************************************************************************************************************
    # THIS WAS THE P3 VERSION AS OF MAY 2021. Test to put in the P1 limit instead (mean + sigma, vs mean+2*sigma in p3)
    # We then investigate the CD25 distribution in this population and set the treshold at one standard dev on the right tail
    mean, median, sigma, maxVal = ag.axisStats(fcs(), xCol=CD25, vI=restingTreg(), scale='bilog', T=200)
    # With this we can set the final xlim
    xlim = mean + 2 * sigma
    xlim = ag.inverseTransformWrapper(xlim, scale='bilog', T=200)
    if xlim > 10000:
        xlim = 10000
    if xlim < 1000:
        xlim = 1000
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations/"+Image_ID_prefix+"-activatedTreg.png"
    activatedTreg = ag.gateCorner(fcs, name="activatedTreg", xCol=CD25, yCol=CD45RA, xThresh=xlim, yThresh=ylim,
                                  xOrientation="upper", yOrientation="lower", parentGate=Tregs, scale='bilog', T=200,
                                  filePlot=fileName)
    fileName =ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations/"+Image_ID_prefix+"-secretingTreg.png"
    secretingTreg = ag.gateCorner(fcs, name="secretingTreg", xCol=CD25, yCol=CD45RA, xThresh=xlim, yThresh=ylim,
                                  xOrientation="lower", yOrientation="lower", parentGate=Tregs, scale='bilog', T=200,
                                  filePlot=fileName)
    fcs.update(ag.AGgate(activatedTreg, Tregs, xCol=CD25, yCol=CD45RA, name="activatedTreg"), QC=True, MFI=True,
               MFI_type="current", extra_MFI=extra_MFI_list, xlim=[0, 3000], ylim=[0, 1000], scale='bilog',
               T=200)
    fcs.update(ag.AGgate(secretingTreg, Tregs, xCol=CD25, yCol=CD45RA, name="secretingTreg"), QC=True, MFI=True,
               MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(restingTreg, Tregs, xCol=CD25, yCol=CD45RA, name="restingTreg"), QC=True, MFI=True,
               MFI_type='current', extra_MFI=extra_MFI_list)

    # Next:   7.2.1)  CD39+ in Tregs subpopulations
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations_CD39/"+Image_ID_prefix+ "-restingTreg_CD39pos.png"
    restingTreg_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="restingTreg_CD39pos", thresh=Tregs_CD39pos_ylim_global, parentGate=restingTreg, scale='bilog', T=1000,population='upper', update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations_CD39/"+Image_ID_prefix+ "-activatedTreg_CD39pos.png"
    activatedTreg_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="activatedTreg_CD39pos", thresh=Tregs_CD39pos_ylim_global, parentGate=activatedTreg, scale='bilog', T=1000,population='upper', update=False, filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations_CD39/"+Image_ID_prefix+ "-secretingTreg_CD39pos.png"
    secretingTreg_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="secretingTreg_CD39pos", thresh=Tregs_CD39pos_ylim_global, parentGate=secretingTreg, scale='bilog', T=1000,population='upper', update=False, filePlot=fileName)
 
    fcs.update(ag.AGgate(restingTreg_CD39pos, restingTreg, xCol=CD39, yCol="FSC 488/10-H", name="restingTreg_CD39pos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(activatedTreg_CD39pos, activatedTreg, xCol=CD39, yCol="FSC 488/10-H", name="activatedTreg_CD39pos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(secretingTreg_CD39pos, secretingTreg, xCol=CD39, yCol="FSC 488/10-H", name="secretingTreg_CD39pos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)

    # Next:   7.2.2)  HLADR in Tregs subpopulations
    fileName =ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations_HLADR/"+Image_ID_prefix+"-restingTreg_HLADRpos.png"
    restingTreg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="restingTreg_HLADRpos", parentGate=restingTreg, thresh=400, population='upper', orientation='vertical', 
                                            scale='bilog',T=400, update=False,filePlot=fileName)
    fileName =ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations_HLADR/"+Image_ID_prefix+"-activatedTreg_HLADRpos.png"
    activatedTreg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="activatedTreg_HLADRpos", parentGate=activatedTreg,thresh=400, population='upper', orientation='vertical', 
                                              scale='bilog',T=400, update=False,filePlot=fileName)
    fileName =ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Treg_subpopulations_HLADR/"+Image_ID_prefix+"-secretingTreg_HLADRpos.png"
    secretingTreg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="secretingTreg_HLADRpos", parentGate=secretingTreg,thresh=400, population='upper', orientation='vertical', 
                                              scale='bilog',T=400, update=False,filePlot=fileName)

    fcs.update(ag.AGgate(restingTreg_HLADRpos, restingTreg, xCol=HLADR, yCol="FSC 488/10-H", name="restingTreg_HLADRpos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(activatedTreg_HLADRpos, activatedTreg, xCol=HLADR, yCol="FSC 488/10-H", name="activatedTreg_HLADRpos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(secretingTreg_HLADRpos, secretingTreg, xCol=HLADR, yCol="FSC 488/10-H", name="secretingTreg_HLADRpos"),QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)

    # ***************************************************************************************************************************************************************************************************************
    # Next:  7.3)  CD194+ Treg subpopulations 
    # add CD194treg subpopulations ('pure' treg + subpops + CD39/HLA-DRs), Tregs subpopulations - activated, secreting, resting
    # First gate out the resting (top) by estimating half-normal distribution on bottom cluster.
    # Then look at the CD25 expression of these resting Tregs, and set the limit between secreting and activated Tregs
    # at some point of the right-tail of the CD25 distribution
    # First do arbitrary cut in the y-axis somewhere between 0-1000
    ylim = ag.valleySeek(fcs, xCol=CD45RA, parentGate=Tregs, interval=[0, 1000], scale='bilog', T=200)
    tmpbottomTreg = ag.gateThreshold(fcs, xCol=CD25, yCol=CD45RA, name="bottomTreg", parentGate=CD194_Tregs,thresh=ylim, population='lower', orientation='horizontal', scale='bilog', T=200,update=False)
    # In this population, find the most dense point
    tmp, tmp1, tmp2, maxVal = ag.axisStats(fcs(), xCol=CD45RA, vI=tmpbottomTreg(), scale='bilog', T=200)
    maxVal = ag.inverseTransformWrapper(maxVal, scale='bilog', T=200)
    mean, sigma = ag.halfNormalDistribution(fcs, xCol=CD45RA, mean=maxVal, direction='left', parentGate=tmpbottomTreg,scale='bilog', T=200)
    # Now we can set the final ylim with big margin (4 std devs):
    ylim = ag.inverseTransformWrapper(mean + 4 * abs(sigma), scale='bilog', T=200)
    if ylim > 1000:
        ylim = 1000
    if ylim < 100:
        ylim = 100
    if ID in corrections_dict_raw['CD194posTregs_ylim 1000']:
        ylim = 1000
    elif ID in corrections_dict_raw['CD194posTregs_ylim 700']:
        ylim = 700
    elif ID in corrections_dict_raw['CD194posTregs_ylim 300']:
        ylim = 300
    # We use this to gate the activated Tregs (Dont save image tho, unneccesary)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194posTreg_subpopulation/"+Image_ID_prefix+"-restingCD194posTreg.png"
    restingCD194posTreg = ag.gateThreshold(fcs, xCol=CD25, yCol=CD45RA, name="restingCD194posTreg", thresh=ylim,parentGate=CD194_Tregs, scale='bilog', T=200, population='upper',
                                           orientation='horizontal', update=False, filePlot=fileName)
    # P1 version of secreting/activated Treg
    bottomTreg = ag.gateThreshold(fcs, xCol=CD25, yCol=CD45RA, name="tmp_bottom_CD194postreg", thresh=ylim,parentGate=CD194_Tregs, scale='bilog', T=200, population='lower',
                                  orientation='horizontal', update=False, filePlot=None)
    mean, median, sigma, maxVal = ag.axisStats(fcs(), xCol=CD25, vI=bottomTreg(), scale='bilog', T=200)
    # With this we can set the final xlim
    xlim = mean + sigma
    xlim = ag.inverseTransformWrapper(xlim, scale='bilog', T=200)
    if xlim > 10000:
        xlim = 10000
    if xlim < 1000:
        xlim = 1000
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194posTreg_subpopulation/"+Image_ID_prefix+"-activatedCD194posTreg.png"
    activatedCD194posTreg = ag.gateCorner(fcs, name="activatedCD194posTreg", xCol=CD25, yCol=CD45RA, xThresh=xlim,yThresh=ylim, xOrientation="upper", yOrientation="lower",
                                          parentGate=CD194_Tregs, scale='bilog', T=200, filePlot=fileName)
    fileName =ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194posTreg_subpopulation/"+Image_ID_prefix+"-secretingCD194posTreg.png"
    secretingCD194posTreg = ag.gateCorner(fcs, name="secretingCD194posTreg", xCol=CD25, yCol=CD45RA, xThresh=xlim,yThresh=ylim, xOrientation="lower", yOrientation="lower",
                                          parentGate=CD194_Tregs, scale='bilog', T=200, filePlot=fileName)
    # Update calls for CD194postregs
    fcs.update(ag.AGgate(activatedCD194posTreg, CD194_Tregs, xCol=CD25, yCol=CD45RA, name="activatedCD194posTreg"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list, 
               xlim=[0, 3000],ylim=[0, 1000],scale='bilog', T=200)
    fcs.update(ag.AGgate(secretingCD194posTreg, CD194_Tregs, xCol=CD25, yCol=CD45RA, name="secretingCD194posTreg"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(restingCD194posTreg, CD194_Tregs, xCol=CD25, yCol=CD45RA, name="restingCD194posTreg"),QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:   7.3.1)  CD39pos in CD194+ Treg subpopulations 
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194treg_subpops_CD39/"+Image_ID_prefix + "-CD39activatedCD194posTreg.png"
    CD39activatedCD194posTreg = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD39activatedCD194posTreg", thresh=CD194_Tregs_CD39pos_ylim_global,parentGate=activatedCD194posTreg, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194treg_subpops_CD39/"+Image_ID_prefix + "-CD39secretingCD194posTreg.png"
    CD39secretingCD194posTreg = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD39secretingCD194posTreg", thresh=CD194_Tregs_CD39pos_ylim_global,parentGate=secretingCD194posTreg, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194treg_subpops_CD39/"+Image_ID_prefix + "-CD39restingCD194posTreg.png"
    CD39restingCD194posTreg = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="CD39restingCD194posTreg", thresh=CD194_Tregs_CD39pos_ylim_global,parentGate=restingCD194posTreg, scale='bilog', T=1000, population='upper',
                                               update=False, filePlot=fileName)

    fcs.update(ag.AGgate(CD39activatedCD194posTreg, activatedCD194posTreg, xCol=CD39, yCol="FSC 488/10-H",name="CD39activatedCD194posTreg"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39secretingCD194posTreg, secretingCD194posTreg, xCol=CD39, yCol="FSC 488/10-H",name="CD39secretingCD194posTreg"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(CD39restingCD194posTreg, restingCD194posTreg, xCol=CD39, yCol="FSC 488/10-H", name="CD39restingCD194posTreg"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)
    # Next:   7.3.2)  HLADR in CD194+ Treg subpopulations 
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194treg_subpops_HLADR/"+Image_ID_prefix + "-restingCD194posTreg_HLADRpos.png"
    restingCD194posTreg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="restingCD194posTreg_HLADRpos",parentGate=restingCD194posTreg, thresh=400, population='upper',
                                                    orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194treg_subpops_HLADR/"+Image_ID_prefix + "-activatedCD194posTreg_HLADRpos.png"
    activatedCD194posTreg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H",name="activatedCD194posTreg_HLADRpos",parentGate=activatedCD194posTreg, thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/CD194treg_subpops_HLADR/"+Image_ID_prefix + "-secretingCD194posTreg_HLADRpos.png"                               
    secretingCD194posTreg_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="secretingCD194posTreg_HLADRpos",parentGate=secretingCD194posTreg, thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)

    fcs.update(ag.AGgate(restingCD194posTreg_HLADRpos, restingCD194posTreg, xCol=HLADR, yCol="FSC 488/10-H", name="restingCD194posTreg_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(activatedCD194posTreg_HLADRpos, activatedCD194posTreg, xCol=HLADR, yCol="FSC 488/10-H",name="activatedCD194posTreg_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(secretingCD194posTreg_HLADRpos, secretingCD194posTreg, xCol=HLADR, yCol="FSC 488/10-H",name="secretingCD194posTreg_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)

# ********************************************************************************************************************************************************************************************************
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
    if ID in corrections_dict_raw['CD194_sigma_1 1.0']:
        CD194_sigma_1 = 1
    elif ID in corrections_dict_raw['CD194_sigma_1 2.0']:
        CD194_sigma_1 = 2
    elif ID in corrections_dict_raw['CD194_sigma_1 -0.5']:
        CD194_sigma_1 = -0.5
    else:
        CD194_sigma_1 = 0.5
    max_1 = ag.inverseTransformWrapper(vX=max+CD194_sigma_1*sigma, scale='bilog', T=1000)
    if max_1 > 3000:
        if ID in corrections_dict_raw['CD194_sigma_2 0.5']:
            CD194_sigma_2 = 0.5
        else:
            CD194_sigma_2 = 1.5
        max_1  = ag.inverseTransformWrapper(vX=max-CD194_sigma_2*sigma, scale='bilog', T=1000)
    ylim = max_1
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17/"+Image_ID_prefix+"-Th_1_2_17.png"
    Th1_Th2_ThG_part1,Th22_Th17,Th9,Th1_Th2_ThG_part2 = ag.quadGate(fcs=fcs, names=['Th1_Th2_ThG_part1','Th22_Th17','Th9','Th1_Th2_ThG_part2'], 
                                        xCol=CCR6, yCol=CD194, xThresh=xlim, yThresh=ylim, parentGate=non_naive_T, scale='bilog', T=1000, filePlot=fileName)
    Th1_Th2_ThG = ag.AGgate(Th1_Th2_ThG_part1() + Th1_Th2_ThG_part2(), non_naive_T, CCR6, CD194, name="Th1_Th2_ThG")
    fcs.update(ag.AGgate(Th1_Th2_ThG, non_naive_T, xCol=CCR6, yCol=CD194, name="Th1_Th2_ThG"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list,xlim=[0, 3000], ylim=[-200, 1000], scale='bilog', T=1000)
    fcs.update(ag.AGgate(Th22_Th17, non_naive_T, xCol=CCR6, yCol=CD194, name="Th22_Th17"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th9, non_naive_T, xCol=CCR6, yCol=CD194, name="Th9"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:  8.0.1)  CD39 in Internediate_Th gating
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17_CD39/"+Image_ID_prefix + "-Th1_Th2_ThG_CD39pos.png"
    Th1_Th2_ThG_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="Th1_Th2_ThG_CD39pos", thresh=1000,parentGate=Th1_Th2_ThG, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17_CD39/"+Image_ID_prefix + "-Th22_Th17_CD39pos.png"
    Th22_Th17_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="Th22_Th17_CD39pos", thresh=1000,parentGate=Th22_Th17, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17_CD39/"+Image_ID_prefix + "-Th9_CD39pos.png"
    Th9_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H",name="Th9_CD39pos", thresh=1000,parentGate=Th9, scale='bilog', T=1000, population='upper',
                                               update=False, filePlot=fileName)

    fcs.update(ag.AGgate(Th1_Th2_ThG_CD39pos, Th1_Th2_ThG, xCol=CD39, yCol="FSC 488/10-H",name="Th1_Th2_ThG_CD39pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th22_Th17_CD39pos, Th22_Th17, xCol=CD39, yCol="FSC 488/10-H",name="Th22_Th17_CD39pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th9_CD39pos, Th9, xCol=CD39, yCol="FSC 488/10-H", name="Th9_CD39pos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)

    # Next:  8.0.2)  HLADR in Internediate_Th gating
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17_HLADR/"+Image_ID_prefix + "-Th1_Th2_ThG_HLADRpos.png"
    Th1_Th2_ThG_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th1_Th2_ThG_HLADRpos",parentGate=Th1_Th2_ThG, thresh=400, population='upper',
                                                    orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17_HLADR/"+Image_ID_prefix + "-Th22_Th17_HLADRpos.png"
    Th22_Th17_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th22_Th17_HLADRpos",parentGate=Th22_Th17, thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_17_HLADR/"+Image_ID_prefix + "-Th9_HLADRpos.png"                               
    Th9_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th9_HLADRpos",parentGate=Th9,thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)

    fcs.update(ag.AGgate(Th1_Th2_ThG_HLADRpos, Th1_Th2_ThG, xCol=HLADR, yCol="FSC 488/10-H", name="Th1_Th2_ThG_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th22_Th17_HLADRpos, Th22_Th17, xCol=HLADR, yCol="FSC 488/10-H",name="Th22_Th17_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th9_HLADRpos,Th9, xCol=HLADR, yCol="FSC 488/10-H",name="Th9_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)

# ********************************************************************************************************************************************************************************************************
    # Next: 9)  Th1_Th2_ThG subpopulation gating
    # calculate the CCR10_A
    temphiCCR10 = ag.gateThreshold(fcs,xCol=CCR10_A, yCol=CXCR3,name="temphiCCR10",parentGate=Th1_Th2_ThG, thresh=0, population='upper',orientation='vertical',scale='bilog', T=500,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CCR10_A,vI=temphiCCR10(),scale="bilog", T=500)
    if ID in corrections_dict_raw['CCR10_sigma 1.5']:
        CCR10_sigma = 1.5
    elif ID in corrections_dict_raw['CCR10_sigma 2.0']:
        CCR10_sigma = 2.0
    else:
        CCR10_sigma = 2.5
    CCR10_lim = ag.inverseTransformWrapper(vX=max+CCR10_sigma*sigma, scale='bilog',T=500)
    # calculate the CXCR3
    temphiCXCR3 = ag.gateThreshold(fcs,xCol=CCR10_A,yCol=CXCR3,name="temphiCXCR3",parentGate=Th1_Th2_ThG, thresh=0, population='upper',orientation='horizontal',scale='bilog', T=500,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CXCR3,vI=temphiCXCR3(),scale="bilog", T=500)
    sigma_converted= ag.inverseTransformWrapper(vX=sigma, scale='bilog',T=500)
    if sigma_converted > -400:
        Th1_2_G_sigma_1 = 0.7
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("Th1_2_G_sigma_1"):
                upper_sigma_match = re.search(r'-?\d+\.?\d*$', keys)
                if upper_sigma_match and ID in ids:
                    Th1_2_G_sigma_1 = float(upper_sigma_match.group())
                    break
        CXCR3_lim = ag.inverseTransformWrapper(vX=max+Th1_2_G_sigma_1*sigma, scale='bilog',T=500)
    else:
        Th1_2_G_sigma_2 = 1.5
        for keys, ids in corrections_dict_raw.items():
            if keys.startswith("Th1_2_G_sigma_1"):
                upper_sigma_match = re.search(r'-?\d+\.?\d*$', keys)
                if upper_sigma_match and ID in ids:
                    Th1_2_G_sigma_2 = float(upper_sigma_match.group())
                    break
        CXCR3_lim = ag.inverseTransformWrapper(vX=max+Th1_2_G_sigma_2*sigma, scale='bilog',T=500)
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G/"+Image_ID_prefix+"-Th_1_2_G.png"
    # NT means non-targeted population
    Th1,NT,ThG,Th2 = ag.quadGate(fcs=fcs, names=['Th1','NS','ThG','Th2'], xCol=CCR10_A, yCol=CXCR3, xThresh=CCR10_lim, yThresh=CXCR3_lim, parentGate=Th1_Th2_ThG, scale='bilog', T=500, filePlot=fileName)
    fcs.update(ag.AGgate(Th1, Th1_Th2_ThG, xCol=CCR10_A, yCol=CXCR3, name="Th1"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list, xlim=[0, 3000], ylim=[-200, 1000], scale='bilog', T=1000)
    fcs.update(ag.AGgate(ThG, Th1_Th2_ThG, xCol=CCR10_A, yCol=CXCR3, name="ThG"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th2, Th1_Th2_ThG, xCol=CCR10_A, yCol=CXCR3, name="Th2"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)

     # Next:  9.0.1)  CD39 in Th1_Th2_ThG subpopulation gating
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G_CD39/"+Image_ID_prefix + "-Th1_CD39pos.png"
    Th1_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="Th1_CD39pos", thresh=1000,parentGate=Th1, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G_CD39/"+Image_ID_prefix + "-Th2_CD39pos.png"
    Th2_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="Th2_CD39pos", thresh=1000,parentGate=Th2, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G_CD39/"+Image_ID_prefix + "-ThG_CD39pos.png"
    ThG_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="ThG_CD39pos", thresh=1000,parentGate=ThG, scale='bilog', T=1000, population='upper',
                                               update=False, filePlot=fileName)

    fcs.update(ag.AGgate(Th1_CD39pos, Th1, xCol=CD39, yCol="FSC 488/10-H",name="Th1_CD39pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th2_CD39pos, Th2, xCol=CD39, yCol="FSC 488/10-H",name="Th2_CD39pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(ThG_CD39pos, ThG, xCol=CD39, yCol="FSC 488/10-H", name="ThG_CD39pos"), QC=True, MFI=True, MFI_type="current", extra_MFI=extra_MFI_list)

    # Next:  9.0.2)  HLADR in Th1_Th2_ThG subpopulation gating
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G_HLADR/"+Image_ID_prefix + "-Th1_HLADRpos.png"
    Th1_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th1_HLADRpos",parentGate=Th1, thresh=400, population='upper',
                                                    orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G_HLADR/"+Image_ID_prefix + "-Th2_HLADRpos.png"
    Th2_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th2_HLADRpos",parentGate=Th2, thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th_1_2_G_HLADR/"+Image_ID_prefix + "-ThG_HLADRpos.png"                               
    ThG_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="ThG_HLADRpos",parentGate=ThG,thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)

    fcs.update(ag.AGgate(Th1_HLADRpos,Th1, xCol=HLADR, yCol="FSC 488/10-H", name="Th1_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th2_HLADRpos,Th2, xCol=HLADR, yCol="FSC 488/10-H",name="Th2_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(ThG_HLADRpos,ThG, xCol=HLADR, yCol="FSC 488/10-H",name="ThG_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)

# ********************************************************************************************************************************************************************************************************
    # Next: 10)  Th17_Th22 subpopulation gating
    # calculate the CCR10_A
    temphiCCR10 = ag.gateThreshold(fcs,xCol=CCR10_A, yCol=CXCR3,name="temphiCCR10",parentGate=Th22_Th17, thresh=0, population='upper',orientation='vertical',scale='bilog', T=500,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CCR10_A,vI=temphiCCR10(),scale="bilog", T=500)
    if ID in corrections_dict_raw['Th_17_22_CCR10_sigma 0.9']:
        Th_17_22_CCR10_sigma = 0.9
    elif ID in corrections_dict_raw['Th_17_22_CCR10_sigma 0.6']:
        Th_17_22_CCR10_sigma = 0.6
    elif ID in corrections_dict_raw['Th_17_22_CCR10_sigma 1.0']:
        Th_17_22_CCR10_sigma = 1.0
    else:
        Th_17_22_CCR10_sigma = 1.5
    CCR10_lim = ag.inverseTransformWrapper(vX=max+Th_17_22_CCR10_sigma*sigma, scale='bilog',T=500)
    # calculate the CXCR3
    temphiCXCR3 = ag.gateThreshold(fcs,xCol=CCR10_A,yCol=CXCR3,name="temphiCXCR3",parentGate=Th22_Th17, thresh=0, population='upper',orientation='horizontal',scale='bilog', T=500,update=False)
    mean,median,sigma,max = ag.axisStats(fcs(), xCol=CXCR3,vI=temphiCXCR3(),scale="bilog", T=500)
    sigma_converted = ag.inverseTransformWrapper(vX=sigma, scale='bilog',T=500)
    if sigma_converted > -400:
        if ID in corrections_dict_raw['Th17_22_sigma_1 1.3']:
            Th17_22_sigma_1 = 1.3
        elif ID in corrections_dict_raw['Th17_22_sigma_1 0.5']:
            Th17_22_sigma_1 = 0.5
        else:
            Th17_22_sigma_1 = 0.8
        CXCR3_lim = ag.inverseTransformWrapper(vX=max+Th17_22_sigma_1*sigma, scale='bilog',T=500)
    else:
        if ID in corrections_dict_raw['Th17_22_sigma_2 1.0']:
            Th17_22_sigma_2 = 1.0
        elif ID in corrections_dict_raw['Th17_22_sigma_2 0.8']:
            Th17_22_sigma_2 = 0.8
        else:
            Th17_22_sigma_2 = 1.5
        CXCR3_lim = ag.inverseTransformWrapper(vX=max+Th17_22_sigma_2*sigma, scale='bilog',T=500)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th22_Th17_subpopulation/"+Image_ID_prefix+"-Th22_Th17_subpopulation.png"
    # NT means non-targeted population
    NT_1,NT_2,Th22,Th17 = ag.quadGate(fcs=fcs, names=[' NS_1','NS_2','Th22','Th17'], xCol=CCR10_A, yCol=CXCR3, xThresh=CCR10_lim, yThresh=CXCR3_lim, parentGate=Th22_Th17, 
                                      scale='bilog', T=500, filePlot=fileName)
    fcs.update(ag.AGgate(Th22, Th22_Th17, xCol=CCR10_A, yCol=CXCR3, name="Th22"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list, xlim=[0, 3000], ylim=[-200, 1000], scale='bilog', T=1000)
    fcs.update(ag.AGgate(Th17, Th22_Th17, xCol=CCR10_A, yCol=CXCR3, name="Th17"), QC=True, MFI=True, MFI_type='current', extra_MFI=extra_MFI_list)
    # Next:  10.0.1)  CD39 in Th17_Th22 subpopulation gating
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th17_Th22_CD39/"+Image_ID_prefix + "-Th17_CD39pos.png"
    Th17_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="Th17_CD39pos", thresh=1000,parentGate=Th17, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)
    
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th17_Th22_CD39/"+Image_ID_prefix + "-Th22_CD39pos.png"
    Th22_CD39pos = ag.gateThreshold(fcs, xCol=CD39, yCol="FSC 488/10-H", name="Th22_CD39pos", thresh=1000,parentGate=Th22, scale='bilog', T=1000,
                                                 population='upper', update=False, filePlot=fileName)

    fcs.update(ag.AGgate(Th17_CD39pos, Th17, xCol=CD39, yCol="FSC 488/10-H",name="Th17_CD39pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th22_CD39pos, Th22, xCol=CD39, yCol="FSC 488/10-H",name="Th22_CD39pos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    
    # Next:  10.0.2)  HLADR in Th17_Th22 subpopulation gating
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th17_Th22_HLADR/"+Image_ID_prefix + "-Th17_HLADRpos.png"
    Th17_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th17_HLADRpos",parentGate=Th17, thresh=400, population='upper',
                                                    orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/T_Panel/Th17_Th22_HLADR/"+Image_ID_prefix + "-Th22_HLADRpos.png"
    Th22_HLADRpos = ag.gateThreshold(fcs, xCol=HLADR, yCol="FSC 488/10-H", name="Th22_HLADRpos",parentGate=Th22, thresh=400, population='upper',
                                                      orientation='vertical', scale='bilog', T=400, update=False,filePlot=fileName)

    fcs.update(ag.AGgate(Th17_HLADRpos,Th17, xCol=HLADR, yCol="FSC 488/10-H", name="Th17_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)
    fcs.update(ag.AGgate(Th22_HLADRpos,Th22, xCol=HLADR, yCol="FSC 488/10-H",name="Th22_HLADRpos"), QC=True, MFI=True, MFI_type="current",extra_MFI=extra_MFI_list)

    # Next: The End
    return fcs


if __name__ == '__main__':

    # Filter the inputfile
    path = "/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_T_Panel_1006/Scripts/old_data_filepath_T_panel_final_1007.csv"   
    all_files = []
    with open(path,"r") as f:
        for line in f:
            line = line.strip()
            all_files.append(line)
    random_inputfile = all_files
    #fcs_path = [x[0] for x in all_files]
    #Image_ID_prefix = [x[1] for x in all_files]
    Gen_exp=ag.AGExperiment(random_inputfile, filters=["T panel"], mask=["Comp"] , compList=[], experiment_name="ImmSight_T_1006-2025", 
                            flourochrome_area_filter=True, QC=True, QCbins=128)
    Gen_exp.apply(gateGeneralDataSet,n_processes=15)
    output = "/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_T_Panel_1006/ImmSight_T_Cell_10062025.txt"
    output_MFI  = "/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_T_Panel_1006/ImmSight_T_Cell__MFI_10062025.txt"
    Gen_exp.printExperiment(output,MFI_file=output_MFI)
    
    

