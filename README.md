# cyp2d6bias
Bias analysis for effect of CYP2D6 on Breast Cancer recurrence and/or death

These files accompany the paper "CYP2D6 Phenotype and Breast Cancer Recurrence and/or death: A Bias Analysis and Meta-Analysis."

meta_data.xls is the dataset necessary to run the analysis
cyp2d6 bias analysis v2.r is the main program to replicate results

It calls the following programs:
cyp2d6 data.r reads in data and formats it appropriately
cyp2d6 bias functions.r are the functions that implement the bias analysis
genetic fxns.r are functions computing activity scores and phenotypes
xtra fxns.r contain various functions for inverting matrices and computing effects
