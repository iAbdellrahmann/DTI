# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 06:31:05 2021

@author: Abdel-Rahman
"""
import streamlit as st

import awesome_streamlit as ast
import home
import prediction_tool
import documentation


ast.core.services.other.set_logging_format()
st.markdown(""" <style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style> """, unsafe_allow_html=True)





PAGES = {
    "Home":home,
    "Prediction tool": prediction_tool,
    "Documentation": documentation,
}


def main():
    """Main function of the App"""
    st.sidebar.title("Navigation")
    selection = st.sidebar.radio("Go to", list(PAGES.keys()))

    page = PAGES[selection]

    with st.spinner(f"Loading {selection} ..."):
        ast.shared.components.write_page(page)
    
    st.sidebar.title("Contribute")
    st.sidebar.info(
        "This an open source project and you are very welcome to **contribute** your awesome "
        "comments, questions, resources and apps "
        "to the [source code](https://github.com/iAbdellrahmann/DTI). "
    )
    st.sidebar.title("About")
    st.sidebar.info(
        """
        This app is maintained by Abdelrahman Moustafa, Khaled Waheed, Zyad Ayman and Ahmed Alaa.
"""
    )


if __name__ == "__main__":
    main()
