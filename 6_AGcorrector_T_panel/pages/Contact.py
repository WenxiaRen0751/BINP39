
import streamlit as st
import pandas as pd

st.set_page_config(layout="wide")
content_container = st.container()
with content_container:
    col1, col2, col3 = st.columns([1, 7, 1])
    with col2:
        st.title("Contacts")

        with st.expander("Contact Information", expanded=True):
            st.write(f"""
            **Author:** Wenxia Ren
                     
            **Program:** Master in Bioinformatics
                     
            **Department:** Department of Biology
                     
            **Faculty:** Faculty of Science
                     
            **Email:** we0751re-s@student.lu.se  
                     
            **Address:** Sölvegatan 30, 223 62 Lund, Sweden 
                           
            **Working Hour:** 24 hours
                     
            """)
        with st.expander("Contact Form", expanded=True):
            data = {
                "Field": ["First_Name", "Last_Name","Email","Field","Question"]
            }
            df = pd.DataFrame(data)
            st.table(df)
            with st.form("contact_form"):
                First_name = st.text_input("First_Name")
                Last_name = st.text_input("Last_Name")
                Email = st.text_input("Email")
                Field = st.text_input("Feild")
                Question = st.text_area("Question")
                submitted = st.form_submit_button("Submit")
                if submitted:
                    st.write(f"Name: {First_name +' '+ Last_name}")
                    st.write(f"Email: {Email}")
                    st.write(f"Question: {Question}")
                    st.write("Thank you for your question, we will answer back ASAP !")