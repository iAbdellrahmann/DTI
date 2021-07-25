# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 06:37:06 2021

@author: Abdel-Rahman
"""

import streamlit as st
import fitz





def write():
    st.title('Documentation')
    st.write("[Download documentation as pdf](https://github.com/iAbdellrahmann/DTI/blob/master/d.pdf)")
    f='D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Frontend/d.pdf'
    if f is not None:
        with fitz.open(f) as doc:
            text = ""
            for page in doc:
                text += page.getText()
            st.write(text) 
