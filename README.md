# MMRF-WGCNA-ANALYSIS
This project focuses on understanding the molecular mechanism behind early multiple myeloma (MM)progression in newly diagnosed 270 Multiple Myeloma patients (downloaded from MMRF COMMPASS study) who were treated early with MM standard therapies:lenalidomide, bortezomide and dexamethasone (RVD).

This project is with the assuption that you have cleaned your expression data; otherwise, this code did not cover data prep. 
For data prep you will need another set of code for that purpose. 

# Project Ojective
The objectives of this research is to identify some of the molecular mechanism using:
Weighted gene co-expression network analysis  (WGCNA) to determine co-expressed gene networks associated with clinical phenotypes of MM and disease progression.

GO-Elite to determine the biological functions and cellular component of the clinically relevant gene clusters.

Receiver Operator characteristic Curve (ROC) to predict if identified genes are biomarkers for MM progression.

Confirm biological relevance of HUB genes using cell-based assays (proliferation, survival, and adhesion)



# Partners
* Morehouse School of Medicine depatment of Microbiology, Biochemistry, and Immunology
* Emory School of Medicine, Winship cancer institute, Clifton Road NE, Atlanta, GA.
* Emory University School of Medicine, Center for Neurodegenerative Diseases, Atlanta, GA 30322, USA
* Georgia Institute of Technology department of applied system engineering, Atlanta, GA.

Methods Used
* Inferential Statistics
* Data Visulization

# Technologies
* R
* Python
* kmploter
* ROC
* Cytoscape

# System Requirements 
* Most of the current new OS will work
* Minimum RAM required 16GB

# Getting Started
* Dwnload a secondary expression data for MM samples and prep it before using this WGCNA code
* Downloa the clinical data for the samples
