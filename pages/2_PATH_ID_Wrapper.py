# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Title:2_PATH_ID_Wrapper.py
Date: 2025-05-16
Author: Wenxia Ren
"""
import streamlit as st
import os
import aligater as ag


# Configure layout
st.set_page_config(layout="wide")
content_container = st.container()
with content_container:
    col1, col2, col3 = st.columns([1, 7, 1])
    with col2:
        st.title("ID_Wrapper")
        st.markdown("") 
        st.markdown("""**Need To Know: This step may take some time, prepare your 🍵 and 🍰**""")
        # File uploader widget, accepts only .txt files
        FCS_PATH = st.text_input("Enter the path for the folder where the fcs files are:")
        if FCS_PATH is not None:
            output_path_1 = "/home/wenxiaren/APP/tmp/fsc_path.txt"
            output_path_2 = "/home/wenxiaren/APP/tmp/Image_ID_Prefix.txt"
            if st.button("Extract IDs"):
                with open(output_path_1 , "w") as output_1, open(output_path_2 , "w") as output_2:
                    for item in os.listdir(FCS_PATH):
                        item_path = os.path.join(FCS_PATH, item)
                        base_path_1 = os.path.dirname(item_path)
                        root_path_1 = os.path.relpath(item_path,base_path_1)
                        if os.path.isdir(item_path) and 'T panel' in item and "Comp" not in item:
                            for subitem in os.listdir(item_path):
                                subitem_path = os.path.join(item_path, subitem)
                                if os.path.isfile(subitem_path) and subitem_path.endswith(".fcs"):
                                    output_1.write(subitem_path+'\n')
                                    base_path_2 = os.path.dirname(subitem_path)
                                    root_path_2 = os.path.relpath(subitem_path,base_path_2)
                                    Image_ID_prefix = root_path_1 + "+" + root_path_2.replace(".fcs", "")
                                    output_2.write(Image_ID_prefix+'\n')

                # Provide a download button
                if st.button("Generate PATH & IDs"):
                    with open(output_path_1, "rb") as file1:
                        st.download_button(label="Download fcs path", data=file1, file_name="fsc_path.txt")
                    with open(output_path_2, "rb") as file2:
                        st.download_button(label="Download ID path", data=file2, file_name="Image_ID_Prefix.txt")  
 


