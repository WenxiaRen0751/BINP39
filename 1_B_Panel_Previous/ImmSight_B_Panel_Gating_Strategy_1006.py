# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Title:ImmSight_B_Panel_Gating_Strategy_1006.py
Date: 2025-05-01
Author: Wenxia Ren
Description:
This program is used for gating the B cell panel in ImmSight project for data generated before 05/30/2025
Note:
1) To run this program, user need to provide the YAML file and the FCS file path. The YAML file contains the correction information for each sample.
2) The output images are saved in the folder: /home/wenxiaren/aligater/plots/phase_III/ImmSight/
3) The event number and MFT for each population are calculated and saved in txt files.
"""

# Import packages
import aligater as ag
import yaml
import os
import random
from math import inf

# Define channel names
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
#corrections_yaml_file = "/home/wenxiaren/aligater/Scripts/ImmSight_B/correction.yaml"
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


# Create output folders
# Define base path and subdirectories
image_output_folder = "/home/wenxiaren/programs/aligater/plots/phase_III/ImmSight"
subdirs = [
    'PBMCs', 'singlets', 'CD45', 'CD19', 'Monocytes', 'NK-cells',
    'HLA-DRpos', 'Dendritic_cells', 'CD45posCD3pos', 'NKT', 'B_quad',
    'IgApos', 'plasmablasts', 'transitionals'
]
# Create base directory and all subdirectories
for subdir in [''] + subdirs:  
    path = os.path.join(image_output_folder, subdir) if subdir else image_output_folder
    os.makedirs(path, exist_ok=True)

# Define functions for gating strategies for Plasmablasts
def gatePlasmablasts(fcs, gate, lim):
    solutions = []
    for testThresh in range(1500, 5000, 500):
        plasmablasts, doublePos, tmp1, tmp2 = ag.quadGate(fcs, names=['plasmablasts', 'doublePos', 'tmp1', 'tmp2'],
                                                        xCol=CD24, yCol=CD38, xThresh=lim, yThresh=testThresh,
                                                        parentGate=gate, scale='bilog', T=200)
        solutions.append([(len(plasmablasts.current) - len(doublePos.current)), testThresh])
    maxVal = -ag.np.inf
    maxThresh = testThresh
    for solution in solutions:
        if solution[0] > maxVal:
            maxVal = solution[0]
            maxThresh = solution[1]
    return maxThresh

# Define functions for gating strategies for General DataSet
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
    
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/PBMCs/"+date+"-"+plate+"-"+sampleName+"-xlim.png"
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
                                    orientation='vertical', population="upper",filePlot=fileName)
    
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/PBMCs/"+date+"-"+plate+"-"+sampleName+"-ylim.png"

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
                                    orientation='horizontal', population="lower",filePlot=fileName)
    
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/PBMCs/"+date+"-"+plate+"-"+sampleName+"-PBMCs.png"
    PBMC=ag.gatePC(fcs,name="PBMC",xCol=FSC,yCol=SSC,center='centroid', adjustAngle=2, widthScale=8, heightScale = 6, parentGate=outline_1,filePlot=fileName)
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/PBMCs/"+date+"-"+plate+"-"+sampleName+"-PBMCs_backgate.png"
    ag.backGate(fcs,population=PBMC,background_population=None,xCol=FSC, yCol=SSC,markersize=0.1, filePlot=fileName)
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

    singlets=ag.gatePC(fcs,xCol=FSC, yCol="FSC 488/10-H",name="singlets",center='density', adjustAngle= adjustAngle_1, widthScale=widthScale_1, heightScale=heightScale_1,parentGate=PBMC,filePlot=fileName)
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

    CD45pos = ag.gateCorner(fcs, "CD45pos", xCol=CD3,yCol=CD45, parentGate=singlets,xThresh=xlim1,yThresh=ylim1,xOrientation='lower',yOrientation='upper',scale='bilog',T=500, filePlot=fileName)
    fcs.update(ag.AGgate(CD45pos,singlets,CD3,CD45,"CD45pos"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)

    # Next: gating the CD45posCD3pos
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/CD45posCD3pos/"+date+"-"+plate+"-"+sampleName+"-CD45posCD3pos.png"
    CD45posCD3pos = ag.gateCorner(fcs, "CD45posCD3pos", xCol=CD3,yCol=CD45, parentGate=singlets,xThresh=xlim1,yThresh=ylim1,xOrientation='upper',yOrientation='upper',scale='bilog',T=500, filePlot=fileName)
    fcs.update(ag.AGgate(CD45posCD3pos,singlets,CD3,CD45,"CD45posCD3pos"), QC=True, xlim=[-5000,100000], ylim=[-5000,100000], scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    # Next: gating NKT, the CD56pos under CD45posCD3pos
    temphiCD3pos = ag.gateThreshold(fcs,"temphiCD3pos", xCol=CD3, yCol=CD56, parentGate=CD45posCD3pos,orientation='vertical',population='upper', thresh=2500,scale='bilog',T=100,filePlot=None)
    ylim = ag.valleySeek(fcs,parentGate=temphiCD3pos, xCol=CD56, interval=[110,1100], scale='bilog',T=100)
    if sample_ids in corrections_dict_raw['NKT Gating']['CUSTOM_y_to_2000']:
        ylim = 2000
    elif sample_ids in corrections_dict_raw['NKT Gating']['y-200']:
        ylim = ylim-200
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/NKT/"+date+"-"+plate+"-"+sampleName+"-NKT.png"
    NKT=ag.gateThreshold(fcs, "NKT", xCol=CD3, yCol=CD56, parentGate=CD45posCD3pos,orientation='horizontal',population='upper', thresh=ylim,scale='bilog',T=100,filePlot=fileName)
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
    CD19pos=ag.gateThreshold(fcs, "CD19pos", xCol=CD19, yCol=CD45, parentGate=CD45pos, population='upper', thresh=xlim,scale='bilog',T=1000,filePlot=fileName)
    fcs.update(ag.AGgate(CD19pos,CD45pos,CD19,CD45,"CD19pos"), QC=True, xlim=[0,20000], ylim=[0,20000],scale='bilog', T=1000,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(CD19neg,CD45pos,CD19,CD45,"CD19neg"), QC=True, xlim=[0,20000], ylim=[0,20000],scale='bilog', T=1000,MFI=True, MFI_type='all', extra_MFI=None)

    # Next: the QUADGate: switchB, preSwitchB, naiveB, exhaustedB (xCol=IgD, yCol=CD27)
    no_clutter_CD19pos = ag.gateThreshold(fcs, "tmp", xCol=IgD, yCol=CD27, thresh=-500, population='upper',
                                          parentGate=CD19pos, scale='bilog', T=500)
    xlim = ag.valleySeek(fcs, xCol=IgD, parentGate=no_clutter_CD19pos, interval=[500, 2500], bins='auto', sigma=2,
                         scale='bilog', T=500)

    rightQuad = ag.gateThreshold(fcs, "tmp", xCol=IgD, yCol=CD27, thresh=xlim, parentGate=no_clutter_CD19pos,
                                 scale='bilog', T=500)

    mean, median, sigma, maxval = ag.axisStats(fcs(), xCol=CD27, vI=rightQuad(), scale='bilog', T=500)
    gaussian_center = ag.inverseTransformWrapper(vX=maxval, T=500, scale='bilog')
    if gaussian_center > 800:  # Cover the case if upper cluster is more dense than lower
        gaussian_center = 0
    mean, sigma = ag.halfNormalDistribution(fcs, xCol=CD27, mean=gaussian_center, direction='left',
                                            parentGate=rightQuad, scale='bilog', T=200)
    ylim = ag.inverseTransformWrapper(vX=mean + 4 * abs(sigma), T=200, scale='bilog')
    if ylim < 300:
        ylim = 300
    if ylim > 500:
        ylim = 500

    bottom_IgD = ag.gateThreshold(fcs, "tmp", xCol=IgD, yCol=CD27, parentGate=no_clutter_CD19pos,
                                  orientation='horizontal', population='lower', thresh=ylim, scale='bilog', T=500)
    # this is for the lower part
    xlim = ag.valleySeek(fcs, xCol=IgD, parentGate=bottom_IgD, sigma=1, interval=[0, 4000], require_local_min=True,
                         scale='bilog', T=500, bins='auto')
    if xlim == inf:
        xlim = 500
    if xlim < 500:
        xlim = 500
    if xlim > 800:
        xlim = 800
 
    top_IgD = ag.gateThreshold(fcs, "tmp", xCol=IgD, yCol=CD27, parentGate=no_clutter_CD19pos, orientation='horizontal',
                               population='upper', thresh=ylim, scale='bilog', T=500)

    # Here is the valley for the lower part
    upper_xlim = ag.valleySeek(fcs, xCol=IgD, sigma=1, parentGate=top_IgD, interval=[100, 1500], require_local_min=True,
                               scale='bilog', T=500, bins='auto')
    if upper_xlim == inf:
        mean, median, sigma, maxval = ag.axisStats(fcs(), xCol=IgD, vI=top_IgD(), scale='bilog', T=500)
        upper_xlim = ag.inverseTransformWrapper(vX=median+1*sigma, T=500, scale='bilog')
    if upper_xlim < xlim:
        mean, median, sigma, maxval = ag.axisStats(fcs(), xCol=IgD, vI=top_IgD(), scale='bilog', T=500)
        upper_xlim = ag.inverseTransformWrapper(vX=mean+1*sigma, T=500, scale='bilog')
    

    if sample_ids in corrections_dict_raw['B_sub Gating']['CD27_ylim_100']:
        ylim = 100
    elif sample_ids in corrections_dict_raw['B_sub Gating']['CD27_ylim_1000']:
        ylim = 1000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['CD27_ylim_2000']:
        ylim = 2000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['CD27_ylim_500']:
        ylim = 500
    
    if sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_600']:
        upper_xlim = 600
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_700']:
        upper_xlim = 700
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_1000']:
        upper_xlim = 1000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_1500']:
        upper_xlim = 1500
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_2000']:
        upper_xlim = 2000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_10000']:
        upper_xlim = 10000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_500']:
        upper_xlim = 500
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_2500']:
        upper_xlim = 2500
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_3000']:
        upper_xlim = 3000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_upper_xlim_5000']:
        upper_xlim = 5000

    if sample_ids in corrections_dict_raw['B_sub Gating']['IgD_xlim_300']:
        xlim = 300
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_xlim_1000']:
        xlim = 1000
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_xlim_1500']:
        xlim = 1500
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_xlim_2500']:
        xlim = 2500
    elif sample_ids in corrections_dict_raw['B_sub Gating']['IgD_xlim_100']:
        xlim = 100
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/B_quad/"+date+"-"+plate+"-"+sampleName+"-quadgate.png"
    switchB, preSwitchB, naiveB, exhaustedB = ag.customQuadGate(fcs, ['switchB', 'preSwitchB', 'naiveB', 'exhaustedB'],
                                                                xCol=IgD, yCol=CD27, threshList=[xlim, upper_xlim, ylim, ylim],
                                                                parentGate=no_clutter_CD19pos, scale='bilog', T=500,filePlot=fileName)
    fcs.update(ag.AGgate(switchB, CD19pos, IgD, CD27, "switchB"), QC=True, MFI=True, extra_MFI=None, MFI_type="all",
               xlim=[-2000, 100000], ylim=[-2000, 70000], scale='bilog', T=500)
    fcs.update(ag.AGgate(preSwitchB, CD19pos, IgD, CD27, "preSwitchB"), QC=True, MFI=True, extra_MFI=None,MFI_type="all"
        ,xlim=[-2000, 100000], ylim=[-2000, 70000], scale='bilog', T=500)
    fcs.update(ag.AGgate(naiveB, CD19pos, IgD, CD27, "naiveB"), QC=True, MFI=True, extra_MFI=None, MFI_type="all"
        ,xlim=[-2000, 100000], ylim=[-2000, 70000], scale='bilog', T=500)
    fcs.update(ag.AGgate(exhaustedB, CD19pos, IgD, CD27, "exhaustedB"), QC=True, MFI=True, extra_MFI=None,MFI_type="all"
        ,xlim=[-2000, 100000], ylim=[-2000, 70000], scale='bilog', T=500)

    # Next: IgApos
    lim = ag.valleySeek(fcs, xCol=IgA, parentGate=switchB, interval=[0, 1750], require_local_min=True, bins='auto',
                    scale='bilog', T=500)
    if lim == ag.np.inf:
        lim = 500
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/IgApos/"+date+"-"+plate+"-"+sampleName+"-IgApos.png"

    if sample_ids in corrections_dict_raw['IgA']['IgA_400']:
        lim = 400
    elif sample_ids in corrections_dict_raw['IgA']['IgA_2000']:
        lim = 2000
    IgApos = ag.gateThreshold(fcs, name="IgAPos", xCol=IgA, yCol=CD19, orientation='vertical', parentGate=switchB,
                            thresh=lim, scale='bilog', T=500, filePlot=fileName)
    fcs.update(ag.AGgate(IgApos, switchB, IgA, CD19, "IgApos"), QC=True, MFI=True, extra_MFI=None, MFI_type="all",
               xlim=[-2000, 100000], ylim=[0, 100000], scale='bilog', T=200)

    # Next: plasmablasts
    plasmablastStep1 = ag.gateThreshold(fcs=fcs, name="plasmatmp", xCol=CD24, yCol=CD38, thresh=5000,
                                        orientation='horizontal', parentGate=switchB, scale='bilog', T=200)
    if not len(plasmablastStep1()) == 0:
        xlim = max(ag.getGatedVector(fcsDF=fcs(), gate=CD24, vI=plasmablastStep1())) + 50
        if xlim <= 150 or xlim >= 750:
            xlim = 500
    else:
        xlim = 500
    CD38thresh = gatePlasmablasts(fcs, switchB, xlim)
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/plasmablasts/"+date+"-"+plate+"-"+sampleName+"-plasmablasts.png"

    if sample_ids in corrections_dict_raw['plasmablasts']['CD24_200']:
        xlim= 200
    elif sample_ids in corrections_dict_raw['plasmablasts']['CD24_1000']:
        xlim = 1000
    elif sample_ids in corrections_dict_raw['plasmablasts']['CD24_2000']:
        xlim = 2000

    if sample_ids in corrections_dict_raw['plasmablasts']['CD38_1000']:
        CD38thresh= 1000
    elif sample_ids in corrections_dict_raw['plasmablasts']['CD38_4000']:
        CD38thresh = 4000
    elif sample_ids in corrections_dict_raw['plasmablasts']['CD38_6000']:
        CD38thresh = 6000
 
    plasmablasts = ag.gateCorner(fcs=fcs, name="plasmablasts", xCol=CD24, yCol=CD38, xThresh=xlim, yThresh=CD38thresh,
                                xOrientation="lower", yOrientation="upper", parentGate=switchB, scale='bilog', T=200,
                                filePlot=fileName)
    fcs.update(ag.AGgate(plasmablasts, switchB, CD24, CD38, "plasmablasts"), QC=True, MFI=True, extra_MFI=None,
               MFI_type="all", xlim=[-2000, 50000], ylim=[-2000, 70000], scale='bilog', T=200)

    # Next: transitionals
    transitionalStep1 = ag.gateThreshold(fcs=fcs, name="transitionals_temp", xCol=CD24, yCol=CD38, thresh=600,
                                     orientation='vertical', population='lower', parentGate=naiveB, scale='bilog',
                                     T=200)
    xmean, xmedian, xsigma, xmaxVal = ag.axisStats(fcs(), xCol=CD24, vI=naiveB(), scale='bilog', T=200)
    xmaxVal = ag.inverseTransformWrapper(xmaxVal, scale='bilog', T=200)
    if xmaxVal < 200:
        xmaxVal = 200
    ylim = ag.densityDelimitation(fcs, xCol=CD38, parentGate=transitionalStep1, bins='auto', interval=[0, 30000],
                                direction='right', limit_threshold=0.10, scale='bilog', sigma=1, T=1000)
    if ylim == ag.np.inf or ylim < 400:
        ylim = 400  # No FMO available
    fileName = ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/transitionals/"+date+"-"+plate+"-"+sampleName+"-transitionals.png"

    if sample_ids in corrections_dict_raw['transitional']['CD38_1500']:
        ylim = 1500
    elif sample_ids in corrections_dict_raw['transitional']['CD38_2000']:
        ylim = 2000
    elif sample_ids in corrections_dict_raw['transitional']['CD38_3000']:
        ylim = 3000

    if sample_ids in corrections_dict_raw['transitional']['CD24_400']:
        xmaxVal = 400

    transitionals = ag.gateCorner(fcs, name="transitionals", xCol=CD24, yCol=CD38, xThresh=xmaxVal, yThresh=ylim,
                                parentGate=naiveB, scale='bilog', T=200, filePlot=fileName)
    fcs.update(ag.AGgate(transitionals, naiveB, CD24, CD38, "transitionals"), QC=True, MFI=True, extra_MFI=None,
               MFI_type="all", xlim=[-2000, 40000], ylim=[-2000, 40000], scale='bilog', T=200)
                                                             
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
    nonclassical_monocytes = ag.gateCorner(fcs, "NonClassicalMonocytes", xCol=CD14, yCol=CD16, parentGate=no_clutter_CD19neg, xThresh=xlim_hi, yThresh=ylim, xOrientation='upper',yOrientation='upper',scale='bilog', T=500,filePlot=fileName)

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

    # New Gating for the monocytes: Classical Monocytes(CD14++CD16-), non-classical monocytes (CD14+CD16++), and intermediate monocytes (CD14++CD16+)
    # keep the original nonclassical_monocytes. now they are considered as the non-classical monocytes (CD14+CD16++)
    # now, find the threshold of CD16neg and CD16pos based on the CD14neg population, use that threshold to seperate the classical_monocytes, the lower part is the new classical_monocytes, the upper part is new intermediate monocytes
    # find the threshold
    # finding the cutting line based on the CD14neg population\
    # trying new range here 
    # CD16low = ag.valleySeek(fcs, parentGate=CD14neg, xCol=CD16, interval=[400,10000], scale='bilog',T=100)
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
    classical_monocytes = ag.gateCorner(fcs,"classical_monocytes", xCol=CD14, yCol=CD16, parentGate=no_clutter_CD19neg, xThresh=xlim_lo, yThresh=CD16low, xOrientation='upper',yOrientation='lower',scale='bilog',T=500,filePlot=fileName)
    # the rest cells from the CD19neg
    Others_index = list(set(CD19neg())-set(classical_monocytes())-set(Intermediate_monocytes())-set(nonclassical_monocytes())-set(CD14neg()))
    Others_from_CD19neg = ag.AGgate(Others_index,CD19neg,CD14,CD16,"Others_from_CD19neg")
    # Now, update the date for each population
    fcs.update(ag.AGgate(nonclassical_monocytes,CD19neg,CD14,CD16,"nonclassical_monocytes"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(Intermediate_monocytes,CD19neg,CD14,CD16,"Intermediate_monocytes"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(classical_monocytes,CD19neg,CD14,CD16,"classical_monocytes"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(CD14neg,CD19neg,CD14,CD16,"CD14neg"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    fcs.update(ag.AGgate(Others_from_CD19neg,CD19neg,CD14,CD16,"Others_from_CD19neg"), QC=True, xlim=[0,10000], ylim=[0,100000],scale='bilog',T=500,MFI=True, MFI_type='all', extra_MFI=None)
    
    # Next: NKCells gate
    # Find y separation between the two rightmost clusters
    tmpHiCD16 = ag.gateThreshold(fcs, "tmpHiCD16", xCol=CD16, parentGate=CD14neg, thresh=1000, population='upper', scale='bilog', T=100)
    ylim_1 = ag.valleySeek(fcs, parentGate=tmpHiCD16, xCol=CD56, interval=[100,4000], scale='bilog',T=100)
    #### Correct the ylim_1
    if sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_y1=1000']:
        ylim_1 = 1000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_y1=800']:
        ylim_1 = 800
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_y1=200']:
        ylim_1 = 200
    elif sample_ids in corrections_dict_raw['NK Gating']['y1+50%']:
        ylim_1 =  ylim_1 *1.5
    elif sample_ids in corrections_dict_raw['NK Gating']['y1-50%']:
        ylim_1 =  ylim_1 *0.5
    else:
        ylim_1 =  ylim_1
    # Find xlim separation between the upper right and lower left clusters
    tmpHiCD56 = ag.gateThreshold(fcs, "tmpHiCD16", xCol=CD56, parentGate=CD14neg, thresh=0, population='upper', scale='bilog', T=100)
    hixlim = ag.valleySeek(fcs, parentGate=tmpHiCD56, xCol=CD16, interval=[1000,10000],scale='bilog',T=100,sigma=3)
    # preliminary cut there
    doublePosNKStep1 = ag.gateCorner(fcs, "doublePosNKStep1", parentGate=CD14neg, xCol=CD16, yCol=CD56, xThresh=hixlim, yThresh=ylim, yOrientation='upper', xOrientation='upper', scale='bilog',T=100)
    # Refine this cut based on top right cluster size
    mean,median,sigma,maxVal = ag.axisStats(fcs(),xCol=CD16, vI=doublePosNKStep1(),scale='bilog',T=100)

    #### Correct the sigma number
    if sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_sigma+4']:
        sigma_n = 4
    elif sample_ids in corrections_dict_raw['NK Gating']['sigma_1']:
        sigma_n = 1
    elif sample_ids in corrections_dict_raw['NK Gating']['sigma_1.5']:
        sigma_n = 1.5
    elif sample_ids in corrections_dict_raw['NK Gating']['sigma_2']:
        sigma_n = 2
    elif sample_ids in corrections_dict_raw['NK Gating']['sigma_2.5']:
        sigma_n = 2.5
    elif sample_ids in corrections_dict_raw['NK Gating']['sigma_0.5']:
        sigma_n = 0.5
    elif sample_ids in corrections_dict_raw['NK Gating']['sigma_3.5']:
        sigma_n = 3.5
    else:
        sigma_n = 3

    bilogXlim = mean- sigma*sigma_n
    hixlim = ag.inverseTransformWrapper(vX=[bilogXlim], scale='bilog', T=100)[0]

    if sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_x=1000']:
        hixlim = 1000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_x=2000']:
        hixlim = 2000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_x=300']:
        hixlim = 300
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_x=3000']:
        hixlim = 3000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_x=4000']:
        hixlim = 4000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_x=800']:
        hixlim = 800
    elif sample_ids in corrections_dict_raw['NK Gating']['CD16=5000']:
        hixlim = 5000
    elif sample_ids in corrections_dict_raw['NK Gating']['CD16=10000']:
        hixlim = 10000
    else:
        hixlim = hixlim 
    
    # Next: CD56+CD16+ NK
    # doublePosNK = ag.gateCorner(fcs, "doublePosNK", parentGate=CD14neg, xCol=CD16, yCol=CD56, xThresh=hixlim, yThresh=ylim, yOrientation='upper', xOrientation='upper', scale='bilog',T=100)
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/NK-cells/"+date+"-"+plate+"-"+sampleName+"-doublePosNK.png"
    doublePosNK = ag.gateCorner(fcs, "CD56posCD16posNK", parentGate=CD14neg, xCol=CD16, yCol=CD56, xThresh=hixlim, yThresh=ylim_1, yOrientation='upper', xOrientation='upper', scale='bilog',T=100, filePlot=fileName)
    fcs.update(ag.AGgate(doublePosNK,CD14neg,CD16,CD56,"doublePosNK"), QC=True, xlim=[-1000,120000], ylim=[-1000,120000], scale='bilog', T=100,MFI=True, MFI_type='all', extra_MFI=None)
    # Next: CD56-CD16+NK
    xThresh2 = 10000
    if sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=3000']:
        xThresh2 = 3000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=4000']:
        xThresh2 = 4000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=5000']:
        xThresh2 = 5000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=6000']:
        xThresh2 = 6000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=7000']:
        xThresh2 = 7000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=8000']:
        xThresh2 = 8000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=1000']:
        xThresh2 = 1000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=10000']:
        xThresh2 = 10000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=15000']:
        xThresh2 = 15000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=20000']:
        xThresh2 = 20000
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['xThresh2=25000']:
        xThresh2 = 25000
    
    adjustAngle_2 = 0
    if sample_ids in corrections_dict_raw['CD56-CD16+NK']['adjustAngle_2=-2']:
        adjustAngle_2 = -2
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['adjustAngle_2=0']:
        adjustAngle_2 = 0
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['adjustAngle_2=2']:
        adjustAngle_2 = 2
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['adjustAngle_2=-1']:
        adjustAngle_2 = -1
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['adjustAngle_2=10']:
        adjustAngle_2 = 10
    
    widthScale_2 = 8
    if sample_ids in corrections_dict_raw['CD56-CD16+NK']['widthScale_2=8']:
        widthScale_2 = 8 
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['widthScale_2=4']:
        widthScale_2 = 4
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['widthScale_2=10']:
        widthScale_2 = 10
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['widthScale_2=6']:
        widthScale_2 = 6
  
    heightScale_2 = 6
    if sample_ids in corrections_dict_raw['CD56-CD16+NK']['heightScale_2=4']:
        heightScale_2 = 4
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['heightScale_2=3']:
        heightScale_2  = 3
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['heightScale_2=2']:
        heightScale_2  = 2
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['heightScale_2=1']:
        heightScale_2  = 1
  
    Thresh2 = -100
    if sample_ids in corrections_dict_raw['CD56-CD16+NK']['Thresh2=0']:
        Thresh2 = 0
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['Thresh2=-50']:
        Thresh2 = -50

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/NK-cells/"+date+"-"+plate+"-"+sampleName+"-CD56negCD16posNK.png"
    CD56negCD16posNK = ag.gateCorner(fcs, "CD56negCD16posNK", parentGate=CD14neg, xCol=CD16, yCol=CD56, xThresh=xThresh2, yThresh=ylim_1, yOrientation='lower', xOrientation='upper', scale='bilog',T=100, filePlot=None)
    CD56low = ag.gateThreshold(fcs, "CD56low", xCol=CD56, parentGate=CD56negCD16posNK, thresh=Thresh2, population='lower', scale='bilog', T=100)

    if sample_ids in corrections_dict_raw['CD56-CD16+NK']['custom']:
        CD56negCD16posNK = ag.gatePC(fcs,xCol=CD16, yCol=CD56,name="CD56negCD16posNK",center='custom', customCenter=[25000,-500],adjustAngle=0, widthScale=6, heightScale=1,parentGate=CD56low,filePlot=None)
    elif sample_ids in corrections_dict_raw['CD56-CD16+NK']['centroid']:  
        CD56negCD16posNK =ag.gatePC(fcs,xCol=CD16, yCol=CD56,name="CD56negCD16posNK",center='centroid',adjustAngle=0, widthScale=8, heightScale=3,parentGate=CD56low,filePlot=None)
    else:
        CD56negCD16posNK=ag.gatePC(fcs,xCol=CD16, yCol=CD56,name="CD56negCD16posNK",center='density', adjustAngle=adjustAngle_2, widthScale=widthScale_2, heightScale=heightScale_2,parentGate=CD56low,filePlot=None)

    ag.backGate(fcs,population=CD56negCD16posNK,background_population=CD14neg,xCol=CD16, yCol=CD56,markersize=0.5, filePlot=fileName,scale='bilog', T=100)
    fcs.update(ag.AGgate(CD56negCD16posNK,CD14neg,CD16,CD56,"CD56negCD16posNK"), QC=True, xlim=[-1000,120000], ylim=[-1000,120000],MFI=True, MFI_type='all', extra_MFI=None)
    # Next: CD56+CD16- NK
    tmplowCD16 = ag.gateThreshold(fcs, "tmplowCD16", xCol=CD16, parentGate=CD14neg, population='lower', thresh=hixlim, scale='bilog', T=100)
    ylim_2=ag.valleySeek(fcs, xCol=CD56, parentGate=tmplowCD16, interval=[600,10000], scale='bilog',T=100)
    # Correct the ylim-2
    if sample_ids in corrections_dict_raw['NK Gating']['y2-20%']:
        ylim_2 = ylim_2 *0.8
    elif sample_ids in corrections_dict_raw['NK Gating']['y2-50%']:
        ylim_2 = ylim_2 *0.5
    elif sample_ids in corrections_dict_raw['NK Gating']['y2-70%']:
        ylim_2 =  ylim_2 *0.3
    elif sample_ids in corrections_dict_raw['NK Gating']['y2-10%']:
        ylim_2 =  ylim_2 *0.9
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_y2=1000']:
        ylim_2 = 1000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_y2=2000']:
        ylim_2 = 2000
    elif sample_ids in corrections_dict_raw['NK Gating']['CUSTOM_y2=3000']:
        ylim_2 = 3000
    else:
        ylim_2 = ylim_2
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/NK-cells/"+date+"-"+plate+"-"+sampleName+"-CD16negCD56posNK.png"
    CD16negCD56posNK = ag.gateCorner(fcs, "CD16negCD56posNK", xCol=CD16, yCol=CD56, parentGate=CD14neg, xOrientation='lower',yOrientation='upper', xThresh=hixlim, yThresh=ylim_2, scale='bilog',T=100,filePlot=fileName)
    fcs.update(ag.AGgate(CD16negCD56posNK,CD14neg,CD16,CD56,"CD16negCD56posNK"), QC=True, xlim=[-1000,120000], ylim=[-1000,120000],MFI=True, MFI_type='all', extra_MFI=None)
   
    # Next: CD56-CD16- NK
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/NK-cells/"+date+"-"+plate+"-"+sampleName+"-doubleNegNK.png"
    if sample_ids in corrections_dict_raw['DNNK']['CD16=1000']:
        hixlim = 1000
    elif sample_ids in corrections_dict_raw['DNNK']['CD16=1500']:
        hixlim = 1500
    elif sample_ids in corrections_dict_raw['DNNK']['CD16=2000']:
        hixlim = 2000
    elif sample_ids in corrections_dict_raw['DNNK']['CD16=2500']:
        hixlim = 2500
    elif sample_ids in corrections_dict_raw['DNNK']['CD16=6000']:
        hixlim = 6000
    elif sample_ids in corrections_dict_raw['DNNK']['CD16=8000']:
        hixlim = 8000
    elif sample_ids in corrections_dict_raw['DNNK']['CD16=4000']:
        hixlim = 4000
    
    if sample_ids in corrections_dict_raw['DNNK']['CD56=1000']:
        ylim_2 = 1000

    doubleNegNK = ag.gateCorner(fcs, "doubleNegNK", xCol=CD16, yCol=CD56, parentGate=CD14neg, xOrientation='lower',yOrientation='lower', xThresh=hixlim, yThresh=ylim_2, scale='bilog',T=100,filePlot=fileName)
    fcs.update(ag.AGgate(doubleNegNK,CD14neg,CD16,CD56,"doubleNegNK"), QC=True, xlim=[-1000,120000], ylim=[-1000,120000],MFI=True, MFI_type='all', extra_MFI=None)  
    # Next:  the rest cells from the CD14neg
    Others_from_CD14neg_index = list(set(CD14neg())-set(CD56negCD16posNK())-set(CD16negCD56posNK())-set(doubleNegNK())-set(doublePosNK()))
    Others_from_CD14neg = ag.AGgate(Others_from_CD14neg_index,CD14neg,CD16,CD56,"Others_from_CD14neg") # the rest cells from the CD19neg
    fcs.update(ag.AGgate(Others_from_CD14neg,CD14neg,CD16,CD56,"Others_from_CD14neg"), QC=True, xlim=[-1000,120000], ylim=[-1000,120000],MFI=True, MFI_type='all', extra_MFI=None) 

    # Next: Dendritic cells
    mean,median,sigma,maxVal=ag.axisStats(fcs(),xCol=HLADR,vI=doubleNegNK(),scale='bilog',T=500)
    maxVal = ag.inverseTransformWrapper(maxVal, scale='bilog',T=500)
    mean,sigma=ag.halfNormalDistribution(fcs,xCol=HLADR,parentGate=doubleNegNK, direction='left', mean=maxVal)
    #xlim = ag.valleySeek(fcs,xCol=HLADR, interval=[10,8000], parentGate=doubleNegNK, scale='bilog',T=1000)
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/HLA-DRpos/"+date+"-"+plate+"-"+sampleName+"-HLADRpos.png"
    if sample_ids in corrections_dict_raw['DC Gating']['sigma+0.5']:
        sigma_n = 3.5
    elif sample_ids in corrections_dict_raw['DC Gating']['sigma+1']:
        sigma_n =  4
    elif sample_ids in corrections_dict_raw['DC Gating']['sigma+1.5']:
        sigma_n =  4.5
    elif sample_ids in corrections_dict_raw['DC Gating']['sigma-0.5']:
        sigma_n =  2.5
    elif sample_ids in corrections_dict_raw['DC Gating']['sigma+2.0']:
        sigma_n = 5
    elif sample_ids in corrections_dict_raw['DC Gating']['sigma+2.5']:
        sigma_n = 5.5
    elif sample_ids in corrections_dict_raw['DC Gating']['sigma+3.5']:
        sigma_n = 6.5
    elif sample_ids in corrections_dict_raw['DC Gating']['xlim=1100']:
        sigma_n = -1
    else:
        sigma_n = 3
    
    if sample_ids in corrections_dict_raw['Second DC Gating']['HLADR']:
        HLADRpos=ag.gateThreshold(fcs, "HLADRpos", xCol=HLADR, yCol=CD56, parentGate=doubleNegNK, thresh=1200, population='upper', scale='bilog', T=500,filePlot=fileName)
    else:
        HLADRpos=ag.gateThreshold(fcs, "HLADRpos", xCol=HLADR, yCol=CD56, parentGate=doubleNegNK, thresh=(mean+sigma_n*sigma), population='upper', scale='bilog', T=500,filePlot=fileName)
    #HLADRpos=ag.gateThreshold(fcs, "HLADRpos", xCol=HLADR, yCol=CD56, parentGate=doubleNegNK, thresh=xlim, population='upper', scale='bilog', T=500,filePlot=fileName)
    fcs.update(ag.AGgate(HLADRpos,doubleNegNK,HLADR,CD56,"HLADRpos"), QC=True, xlim=[0,10000], ylim=[0,10000], scale='bilog', T=500,MFI=True, MFI_type='all', extra_MFI=None)

    # Next: Dendritic subtypes
    ylim = ag.valleySeek(fcs, xCol=CD123, parentGate=HLADRpos, interval=[200,8000], scale='bilog',T=500)
    CD123neg = ag.gateThreshold(fcs, "CD123neg", xCol=CD141, yCol=CD123, parentGate=HLADRpos, thresh=ylim, orientation='horizontal', population='lower', scale='bilog', T=500)
    xlim = ag.valleySeek(fcs, xCol=CD141, parentGate=CD123neg, interval=[500,30000], scale='bilog',T=500)
 
    if sample_ids in corrections_dict_raw['DC_sub Gating']['xlim+10%']:
        xlim = xlim*1.1
    elif sample_ids in corrections_dict_raw['DC_sub Gating']['CUSTOM_x=10000']:
        xlim = 10000
    elif sample_ids in corrections_dict_raw['DC_sub Gating']['xlim+20%']:
        xlim = xlim*1.2
    elif sample_ids in corrections_dict_raw['DC_sub Gating']['xlim+30%']:
        xlim = xlim*1.3
    elif sample_ids in corrections_dict_raw['DC_sub Gating']['xlim+50%']:
        xlim = xlim*1.5
    else:
        xlim = xlim

    if sample_ids in corrections_dict_raw['DC_sub Gating']['CUSTOM_y=1000']:
        ylim = 1000
    else:
        ylim = ylim 
    
    if sample_ids in corrections_dict_raw['second DC gating']['CD141=8000']:
        xlim =  8000
    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Dendritic_cells/"+date+"-"+plate+"-"+sampleName+"-pDCs.png"
    pDC = ag.gateCorner(fcs, "pDC", parentGate=HLADRpos, xCol=CD141, yCol=CD123, xThresh=xlim, yThresh=ylim, xOrientation='lower', yOrientation='upper', scale='bilog', T=500,filePlot=fileName)
    fcs.update(ag.AGgate(pDC,doubleNegNK,CD141,CD123,"pDC"), QC=True, xlim=[0,20000], ylim=[0,20000], scale='bilog', T=500,MFI=True, MFI_type='all', extra_MFI=None)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Dendritic_cells/"+date+"-"+plate+"-"+sampleName+"-cDC1s.png"
    cDC1 = ag.gateCorner(fcs, "cDC1", parentGate=HLADRpos, xCol=CD141, yCol=CD123, xThresh=xlim, yThresh=ylim, xOrientation='upper', yOrientation='lower', scale='bilog', T=500,filePlot=fileName)
    fcs.update(ag.AGgate(cDC1,doubleNegNK,CD141,CD123,"cDC1"), QC=True, xlim=[0,20000], ylim=[0,20000],MFI=True, MFI_type='all', extra_MFI=None)

    fileName=ag.AGConfig.ag_home+"/plots/phase_III/ImmSight/Dendritic_cells/"+date+"-"+plate+"-"+sampleName+"-cDC2s.png"
    cDC2 = ag.gateCorner(fcs, "cDC2", parentGate=HLADRpos, xCol=CD141, yCol=CD123, xThresh=xlim, yThresh=ylim, xOrientation='lower', yOrientation='lower', scale='bilog', T=500,filePlot=fileName)
    fcs.update(ag.AGgate(cDC2,doubleNegNK,CD141,CD123,"cDC2"), QC=True, xlim=[0,20000], ylim=[0,20000],MFI=True, MFI_type='all', extra_MFI=None)

    return fcs

    
if __name__ == '__main__':

    # Filter the inputfile
    path = "/home/wenxiaren/aligater/Scripts/ImmSight_B/Scripts/old_data_path_B_panel_final_1007.csv"   
    all_files = []
    with open(path,"r") as f:
        for line in f:
            line = line.strip()
            all_files.append(line)
    random_inputfile = all_files
    Gen_exp=ag.AGExperiment(random_inputfile, filters=["B panel"], mask=["Comp"] , compList=[], 
        experiment_name="ImmSight_B_Cell_10062025", flourochrome_area_filter=True, QC=True, QCbins=128)
    Gen_exp.apply(gateGeneralDataSet,n_processes=15)
    output = "/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_B_Panel_1006/ImmSight_B_Cell_10062025.txt"
    output_MFI  = "/home/wenxiaren/cbio4/projects/Wenxia/ImmSight_B_Panel_1006/ImmSight_B_Cell_MFI_10062025.txt"
    Gen_exp.printExperiment(output,MFI_file=output_MFI)
