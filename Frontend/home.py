# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 07:04:24 2021

@author: Abdel-Rahman
"""

import streamlit as st

def write():
    st.title('Drug Target Interaction Prediction')
    st.markdown("""
               * Motivated by the limitations of the exiting methods for the prediction of
    the potential DTIs, the types of the effects a drug exert on the target,
    and aiming to further improve their prediction accuracy, we present a
    method that utilize a heterogeneous drug-target graph that contains
    information about DTIs, effect types of DTIs,
    as well as multiple similarities between drugs and multiple similarities between target proteins.
                """)
    st.write("""
    ## Drug Target Interaction Visualization
             """)
    st.video(
                "D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Frontend/video.mp4"
            )
    