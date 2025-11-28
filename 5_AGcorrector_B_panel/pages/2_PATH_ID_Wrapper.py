# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Title:2_PATH_ID_Wrapper.py
Date: 2025-05-16
Author: Wenxia Ren
Description:
This program is used to generate the absolute path for the raw .fcs files and the image IDs for the ImmSight B panel data.
The image IDs will be used to match with the images in the later steps.
The program allows users to filter the .fcs files based on a specified cutoff date.
Note:
1) User need to provide the path where the raw .fcs files are located.
2) User can select a cutoff date and filtering condition (earlier or later than the cutoff
date) to filter the .fcs files.
"""
import streamlit as st
import os,re
import aligater as ag
from datetime import datetime
# Configure layout
st.set_page_config(layout="wide")
content_container = st.container()
with content_container:
    col1, col2, col3 = st.columns([1, 7, 1])
    with col2:
        st.title("Generating FCM data path & Image_IDs")
        st.markdown("") 
        st.markdown("""**Please enter the path where the raw data are, such as: /home/wenxiaren/cbio4/data/Aitzkoa/ImmSight**""")
        st.markdown("""**It will output two files**""")
        st.markdown("""File 1: the absolute path for the raw .fcs file""")
        st.markdown("""File 2: the image ID which will be used to match with the images in the later steps""")
        # Date filter options
        st.markdown("**Selecting a cutoff sampling date to filter the .fcs files:**")
        cutoff_date_ui = st.date_input(
            "**Select cutoff date:**",
            value=datetime(2025, 5, 13)
        )
        date_condition = st.radio(
            "Select filtering condition:",
            ("Later than this date", "Earlier than this date")
        )
        # Convert UI date to datetime object
        cutoff_date = datetime.combine(cutoff_date_ui, datetime.min.time())
        # File uploader widget, accepts only .txt files
        FCS_PATH = st.text_input("**Enter the path where the fcs files are:**")
        st.markdown("such as:/home/wenxiaren/cbio4/data/Aitzkoa/ImmSight")
        if FCS_PATH is not None:
            output_path_1 = "/home/wenxiaren/AGcorrector_B_panel/tmp/fsc_path.txt"
            output_path_2 = "/home/wenxiaren/AGcorrector_B_panel/tmp/Image_ID_Prefix.txt"
            st.markdown("Presse the button below to extract the fcs paths and IDs:") 
            if st.button("Extract IDs"):
                with open(output_path_1 , "w") as output_1, open(output_path_2 , "w") as output_2:
                    for item in os.listdir(FCS_PATH):
                        item_path = os.path.join(FCS_PATH, item)
                        base_path_1 = os.path.dirname(item_path)
                        root_path_1 = os.path.relpath(item_path,base_path_1)
                        if os.path.isdir(item_path) and 'B panel' in item and "Comp" not in item:
                            for subitem in os.listdir(item_path):
                                subitem_path = os.path.join(item_path, subitem)
                                if os.path.isfile(subitem_path) and subitem_path.endswith(".fcs"):
                                    matched_date = re.search(r"(\d{8})",subitem_path)
                                    if matched_date:
                                        matched_date = str(matched_date.group(1))
                                        # Convert to datetime object
                                        date_obj = datetime.strptime(matched_date, "%Y%m%d")
                                        if date_condition == "Earlier than this date":
                                            if date_obj <= cutoff_date:
                                                date_plate = ag.getFileName(ag.getParent(subitem_path))
                                                date = date_plate.split("B panel")[0]
                                                plate = date_plate.split("B panel")[1]
                                                sampleName=ag.getFileName(subitem_path)
                                                # sample_ids are used to match wtih the id in the yaml file
                                                sample_ids = date+"-"+plate+"-"+sampleName
                                                output_1.write(subitem_path+'\n')
                                                output_2.write(sample_ids+'\n')
                                        else:  
                                            if date_obj >= cutoff_date:
                                                date_plate = ag.getFileName(ag.getParent(subitem_path))
                                                date = date_plate.split("B panel")[0]
                                                plate = date_plate.split("B panel")[1]
                                                sampleName=ag.getFileName(subitem_path)
                                                # sample_ids are used to match wtih the id in the yaml file
                                                sample_ids = date+"-"+plate+"-"+sampleName
                                                output_1.write(subitem_path+'\n')
                                                output_2.write(sample_ids+'\n')

                                    
                st.success("Extraction completed!")
                st.markdown("Please download the two files below:")
                st.markdown("**File 1:** fcs_path.txt")
                st.markdown("**File 2:** Image_ID_Prefix.txt")
                # Provide a download button
                if st.button("Download FCM data path & Image_IDs"):
                    st.success("Files are ready for download!")
                    with open(output_path_1, "rb") as file1:
                        st.download_button(label="Download fcs path", data=file1, file_name="fsc_path.txt")
                    with open(output_path_2, "rb") as file2:
                        st.download_button(label="Download ID path", data=file2, file_name="Image_ID_Prefix.txt")  



                       
                                        
               
 


