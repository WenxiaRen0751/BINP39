# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Title:1_📁_Homepage.py
Date: 2025-06-02
Author: Wenxia Ren

Description:
AGcorrector is a web-based application that extends the AliGater framework to provide an interactive environment for the quality control and refinement of flow cytometry gating outputs. 
Built on Python 3.12.2 and Streamlit 1.43.2, AGcorrector enables researchers to inspect gating results and dynamically adjust parameters, including thresholds, filters, 
and smoothing with immediate visual feedback on resulting images and metrics. Adjustments can be applied globally or to individual samples, and the tool supports reproducibility through versioned YAML presets, 
side-by-side strategy comparisons, and one-click export of publication-quality figures, quantitative metrics, and the exact configuration used.

"""
import streamlit as st

# Set up the Streamlit page configuration
st.set_page_config(
    page_title="AGcorrector (T_panel)",
    page_icon="🔴",
    layout="wide"
)
# Create a container to hold the page content
content_container = st.container()
with content_container:
    # Define a three-column layout where the middle column (col2) is wider
    col1, col2, col3 = st.columns([1, 8, 1])
    # Place main content in the center column
    with col2:
        # Display the main title
        st.title("AGcorrector (T_panel)")
        # Add a new line for spacing
        st.markdown("") 
        st.markdown("**🔈AGcorrector stands for correcting gating parameters for AliGater, AGcorrector is written in Python (v3.12.2). It depends on python packages pandas (v 2.2.3) and streamlit (v1.43.2)**")
        st.markdown("AliGater: https://github.com/LudvigEk/aligater.git") 
        st.markdown("")
        # Project Introduction
        st.markdown("**🔎 Research Project Title:**")
        st.markdown("""Extending AliGater: Advancing Gating Strategies and Machine Learning for Flow Cytometry Analysis
        """)
        st.markdown("**🔎 What does AGcorrector do:**")
    
        st.write("""
        Enable researchers to inspect gating results and dynamically adjust parameters, including thresholds, filters, with immediate visual feedback on resulting images and metrics. 
        Adjustments can be applied to individual samples later with the combination of AliGater.

        """)
        # Add a new line for spacing
        st.markdown("") 
        # Dataset introduction
        st.markdown("**📚Input for AGcorrector:**")

        st.write("""
        📂 Images generated using AliGater      
        📂 The gating strategy used in the current gating steps         
        📂 The raw .fcs file
                 
        """)

        # Add a new line for spacing
        st.markdown("") 
        # Dataset introduction
        st.markdown("**📚Output for AGcorrector:**")

        st.write("""
        📋 YAML file containing sample identifiers and corresponding optimized values""")
        # Add a new line for spacing
        st.markdown("") 
        st.markdown("**🔄The workflow of AGcorrector:**")
        # Workflow introduction
        col1, col2,col3 = st.columns([2,0.3,1])
        with col1:
            st.image("/home/wenxiaren/ImmSight_B_visualization/WORKFLOW.png", caption="**Fig.1. The workflow of AGcorrector**", width=700)
            st.write("""
                The AGcorrector workflow includes sample ID and path wrapper, QC & refinement to fine-tune parameters
                and generate YAML files which wiil be merged to the original gating strategy and then rerun AliGater
                for future Expert Review to verify improvements.
        """)
        with col2:
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("")
            st.markdown("<h1 style='text-align: center;'>⏭</h1>", unsafe_allow_html=True)
        with col3:
            st.image("/home/wenxiaren/Desktop/workflow.png", caption="**Fig.2. The combination of AliGater and AGcorrector**", width=400)
            st.write("""
            Ther adjustments from AGcorrector can be applied to individual samples with the combination of AliGater.
        """)
          
        
