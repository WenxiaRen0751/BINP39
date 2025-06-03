# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Title:2_PBMC.py
Date: 2025-05-16
Author: Wenxia Ren
"""

from curses import BUTTON1_RELEASED
from matplotlib import pyplot as plt
import streamlit as st
import aligater as ag
from math import inf
import os
from PIL import Image
import re
import yaml
import json
# Defining the gating startegy
# From the past experience, in the PBMC gating strategies, these four parameters may significantly affect the final gating result.
def PBMC_Gating(fcs, direction, xboundaries,no_clutter_thresh, cut_off_xlim, cut_off_ylim):

    SSC = "SSC 488/10-A"
    FSC = "FSC 488/10-A"
    FSC_H = "FSC 488/10-H"
    CD39 = "BB515 CD39-A"
    CCR10_A = "PerCP-Cy5.5 CCR10_A"
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
    date_plate = ag.getFileName(ag.getParent(fcs.filePath))
    date = date_plate.split("T panel")[0]
    plate = date_plate.split("T panel")[1]
    sampleName=ag.getFileName(fcs.filePath)


    # PBMC gating
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
    
    ### direction_1, 
    startY_1= ylim_bot+10000
    xboundaries_1 = xboundaries
    direction_1 = direction
    PBMC_step2 = ag.horizontalPath(fcs, name="PBMC",
                    xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                    startY=startY_1, endY=ylim_top, xboundaries=xboundaries_1,
                    yboundaries=[ylim_bot-5000,ylim_top+5000],
                    leftRight=True , direction=direction_1,
                    maxStep=2, phi=0.1, bins=100, sigma=1,
                    scale='linear', parentGate=PBMC_step1)
    
    # set the path for all the images which need to be checked later
    # in streamlit, use st.image to load the image, st.image takes the image path as input
    # all the image in this function will be saved to a tmp dir which is as the same path as the APP
    fileName_2 ="/home/wenxiaren/APP/tmp/"+date+"-"+plate+"-"+sampleName+"- xlim.png"
    outline = ag.gateThreshold(fcs, name="outline", xCol=FSC, yCol=SSC, thresh=cut_off_xlim,parentGate=PBMC_step2,
                                    orientation='vertical', population="upper",filePlot=fileName_2)
    fileName_3 ="/home/wenxiaren/APP/tmp/"+date+"-"+plate+"-"+sampleName+"- ylim.png"
    outline_1 = ag.gateThreshold(fcs, name="outline_1", xCol=FSC, yCol=SSC, thresh=cut_off_ylim,parentGate=outline,
                                    orientation='horizontal', population="lower",filePlot=fileName_3)
    fileName_4="/home/wenxiaren/APP/tmp/"+date+"-"+plate+"-"+sampleName+"- PBMC.png"
    PBMC=ag.gatePC(fcs,name="PBMC",xCol=FSC,yCol=SSC,center='centroid', adjustAngle=2, widthScale=8, heightScale = 6, parentGate=outline_1,filePlot=fileName_4)
    fileName_1 ="/home/wenxiaren/APP/tmp/"+date+"-"+plate+"-"+sampleName+"- PBMC_back.png"
    ag.backGate(fcs,population=PBMC,background_population=None,xCol=FSC, yCol=SSC,markersize=0.1, filePlot=fileName_1)
    return fileName_1, fileName_4, fileName_3,fileName_2

# Configure layout
st.set_page_config(layout="wide")

# --- SESSION STATE: to track which sample we're viewing ---
if "current_sample_index" not in st.session_state:
    st.session_state.current_sample_index = 0

# Sidebar: parameter controls
with st.sidebar:
    st.markdown("### Parameters")
    direction = st.selectbox("Select direction", ["both","up"], index=0)
    xboundaries =  st.slider("Select a range", min_value=0, max_value=214000, value=(65000, 214000), step=1000)
    if xboundaries[0] == xboundaries[1]:
        xboundaries = None
    else:
        xboundaries = list(xboundaries)
    no_clutter_thresh = st.selectbox("No_clutter_thresh", (214000, 220000, 230000, 250000))
    cut_off_xlim = st.number_input("cut_off_xlim", min_value=0, max_value=250000, value=75000, step=500)
    cut_off_ylim = st.number_input("cut_off_ylim", min_value=0, max_value=250000, value=110000, step=500)
    Another_Cell_Distribution = st.checkbox("Another_Cell_Distribution",value=False)
    Update = st.button("Update Image")
    Confirm = st.button("Save Parameters")
    if st.button("Done! Check next !"):
        st.session_state.current_sample_index += 1  # Go to next sample
    if st.button("Go BacK! Check the previous !"):
        st.session_state.current_sample_index -= 1  # Go to the previous one
    
# For row1, the input is the IDs list, the absolute path for all the fcs file and the image folder.
# The layout of the first row, just to input all the files needed for the next process.
row1 = st.columns(1)
for col in row1:
    content_container = col.container(height=160)
    with content_container:
        col1, col2, col3 = st.columns([1,16,1])
        with col2:
            st.title("Customed Gating")
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
                    # modified_id = re.findall(r"-\d{8}-\d{4}-[a-zA-Z]\d",current_id)
                    modified_id  = current_id.replace("+","/")+".fcs"
                    modified_id = "/home/wenxiaren/cbio4/data/Aitzkoa/ImmSight/" + modified_id
                    #    Aitzkoa-2023-05-30-- with CCR10-ImmSight-Cat-20230530-1657 (1)-A9 7913-3
                    #    ID problem ????????????
        
                    # modified_id  = modified_id[0]
                    # last_dash_index = modified_id.rfind('-')
                    # s_modified = modified_id[:last_dash_index] + '/' + modified_id[last_dash_index+1:]
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
                if 'fcs' in st.session_state and st.session_state['fcs'] is not None:
                    fcs = st.session_state['fcs']
                    PATH = fcs.filePath
                    st.write(f"PATH:{PATH}")
                    # call the function, it will return a tuple which contains the absolute imgae path.
                    fig_path = PBMC_Gating(fcs, direction,xboundaries,no_clutter_thresh, cut_off_xlim, cut_off_ylim)
                    #set the layput of the updated images
                    cols = st.columns(4)
                    clo_i = 0
                    for i  in fig_path:
                        with cols[clo_i]:
                            st.image(i, use_container_width=True)
                            clo_i+=1 
            if Confirm:
                if 'current_id' in st.session_state:
                    current_id = st.session_state['current_id']
                    combo_key = f"direction_{direction}_xboundaries_{json.dumps(xboundaries)}"
                    para_dict = {
                        f"No_Clutter_Threshold {no_clutter_thresh}": no_clutter_thresh,
                        f"X_Lim {cut_off_xlim}": cut_off_xlim,
                        f"Y_Lim {cut_off_ylim}": cut_off_ylim,
                        f"Another_Cell_Distribution {Another_Cell_Distribution}":Another_Cell_Distribution,
                    }
                YAML_PATH = "/home/wenxiaren/APP/tmp/T_PBMC_correction.yaml"
                if os.path.exists(YAML_PATH):
                    with open(YAML_PATH,"r") as in1:
                        correction_yaml = yaml.safe_load(in1) or {}
                else:
                    correction_yaml = {}
                # Step 2: Update combo_key entry
                correction_yaml.setdefault(combo_key, [])
                if current_id not in correction_yaml[combo_key]:
                    correction_yaml[combo_key].append(current_id)
                # Step 3: Update other individual parameters
                for key, val in para_dict.items():
                    correction_yaml.setdefault(key, [])
                    if current_id not in correction_yaml[key]:
                        correction_yaml[key].append(current_id)
                # Step 4: Save
                with open(YAML_PATH, "w") as out1:
                    yaml.dump(correction_yaml, out1, sort_keys=False)




        

            

    