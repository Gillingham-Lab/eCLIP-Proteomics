import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import os
from _plotly_utils.utils import PlotlyJSONEncoder
from PIL import Image
import io
import requests
import json

#---DEFINE PATH & FILES---
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'Data/')

MV411WT_SC026_CV070 = "P284-SC-TMT-red-all-20180709.xlsx"
MV411WT_SC063_Indi = "P284_SC_TMT_blue_all_20180709.xlsx"
HCT116WT_Indi_E7820 = "SQ_P284_AS_SC_SC240243__TMT9plex_240_PROTEIN.xlsx"
HCT116WT_CQS_Tas = "SQ_P284_AS_SC_SC240243__TMT9plex_241_PROTEIN.xlsx"
HCT116WT_PR041_PR064 = "SQ_P284_AS_SC_SC240243__TMT9plex_243_PROTEIN.xlsx"
HCT116KO_Indi_E7820_CQS_Tas = "P284_AS_SC_SC269TMT_15plex_PROTEIN.xlsx"


#---DEFINE DATAFRAMES---
df1 = pd.read_excel(filename+MV411WT_SC063_Indi, usecols=["geneName", "log2ratio_Condition3", "qValue_Condition3"])
df2 = pd.read_excel(filename+'RBM39.Wang_el_al.peaks.xlsx', usecols=["gene_name", "scores_1", "scores_2", "scores_3", "scores_4", "scores_5", "scores_6"]) # for individual scores : CASE 1
#df2 = pd.read_excel(filename+'RBM39.Wang_el_al.peaks.xlsx', usecols=["gene_name", "total_score"]) #for total_score : CASE 2

#---FOR HIGHEST INDIVIDUAL SCORES ## Activate for CASE 1 / Deactivate for CASE 2
total_score = np.amax(df2, axis=1)
df2 = df2.drop(["scores_1", "scores_2", "scores_3", "scores_4", "scores_5", "scores_6"], axis=1)
df2["total_score"] = total_score
df2 = df2[["total_score", "gene_name"]]

#keep next 2 lines
df2 = df2.sort_values(by=['total_score'], ascending=False) #important to sort for FILTER 1
df2 = df2.dropna()
#df2 = df2.head(2000)


#---GENE DICTIONARY - List of genes which are in ALL Proteomics datasets---
#df3 = pd.read_excel(filename+MV411WT_SC063_Indi, sep=";", error_bad_lines=False, usecols=["geneName"])
#df4 = pd.read_excel(filename+HCT116WT_Indi_E7820, sep=";", error_bad_lines=False, usecols=["geneName"])
#df5 = pd.read_excel(filename+HCT116WT_CQS_Tas, sep=";", error_bad_lines=False, usecols=["geneName"])
#df6 = pd.read_excel(filename+HCT116KO_Indi_E7820_CQS_Tas, sep=";", error_bad_lines=False, usecols=["geneName"])

#geneDict = df3.merge(df4, left_on="geneName", right_on="geneName", how="inner")
#geneDict = geneDict.merge(df5, left_on="geneName", right_on="geneName", how="inner")
#geneDict = geneDict.merge(df6, left_on="geneName", right_on="geneName", how="inner")
#geneDict = geneDict.set_index("geneName", drop= False)
#geneDict.columns = ["geneLib"]
#geneDict = geneDict.loc[~geneDict.index.duplicated(keep='first')] #FILTER 1.1 - Merge Gene List from all datasets


#---MERGE OF ECLIP AND PROTEOMCIS---
#df1 = geneDict.merge(df1, right_on="geneName", left_on="gene_name", how="inner")  #FILTER 1.2 - Only if GENE DICTIONARY WAS USED
newDataframe = df1.merge(df2, left_on="geneName", right_on="gene_name", how="inner") 
newDataframe = newDataframe.sort_values(by=['total_score'], ascending=False) #sorting is important for percentile selction
newDataframe.columns = ["geneName","log2ratio","q_p","total_score","gene_name"]
#newDataframe = newDataframe.head(1000)
#newDataframe = newDataframe[400:1000]

newDataframe = newDataframe.set_index("geneName", drop=False)
newDataframe = newDataframe.loc[~newDataframe.index.duplicated(keep='first')] #FILTER 2 - remove clusters beloning to the same gene from eCLIP, keeping the one with the highest score
newDataframe["total_score"].fillna(0, inplace=True) # fill NA with 0


#---SELECT PERCENTILE OF SCORES---
num_rows = newDataframe.shape[0]
num_rows = num_rows*0.001
newDataframe = newDataframe[int(num_rows):] #rmeoves top percentile (as defined in previous line) scores and keeps rest
print(newDataframe)


#---DISTRIBUTION OF SCORES USED---
dist = newDataframe['total_score']
dist = pd.concat([dist,newDataframe['total_score']], axis=1, ignore_index=True)


#---GRAPH VISUALIZATION---
newDataframe["q_p"] = -1*np.log10(newDataframe["q_p"])
fig = px.scatter(newDataframe, y="q_p", x="log2ratio", color='total_score', hover_name="geneName", color_continuous_scale=["#d8e7f5", "#083573"]) #, text="geneName"
#color_continuous_scale=px.colors.sequential.Cividis_r
#fig.update_layout(yaxis_type="log")


#---MARK A CERTAIN GENE---
geneOfInterest = "KIF20B"
#geneOfInterest2 = "ZNF791"
x_1 = newDataframe.loc[geneOfInterest,"log2ratio"]
#x_1 = x_1.loc[~x_1.index.duplicated(keep='first')]
#x_1 = x_1[0]
#x_2 = newDataframe.loc[geneOfInterest2,"log2ratio"]
y_1 = newDataframe.loc[geneOfInterest,"q_p"]
#y_1 = y_1.loc[~y_1.index.duplicated(keep='first')]
#y_1 = y_1[0]
#y_2 = newDataframe.loc[geneOfInterest2,"q_p"]


#---PLOT LAYOUT AND OUTPUT---
fig.update_layout(yaxis_title="-log10_pValue")
fig.update_layout(xaxis_title="log2ratio_Indisulam_DMSO_Mv411", font=dict(family="Arial"), plot_bgcolor="#ffffff")
fig.update_traces(textposition='top center', marker=dict(opacity=0.5, size=10))
fig.update_xaxes(ticks="outside")
fig.update_yaxes(ticks="outside")
fig.write_image(filename+"fig1.jpeg")
fig.write_image(filename+"fig1.svg")


#---SHOW MARKED GENE ON PLOT---
#fig.add_annotation(x=x_1, y=y_1, text=geneOfInterest)
#fig.add_annotation(x=x_2, y=y_2, text=geneOfInterest2)
#fig.update_annotations(dict(xref="x", yref="y", showarrow=True,arrowhead=0,ax=0,ay=-40))

fig.show()


#---SHOW DISTRIBUTION OF SCORES---
#dist.columns = ["eCLIP", "Merge", "Top 1000"]
dist.plot.hist(bins=200, alpha=0.2)
plt.show()


#---EXPORT---
#export_csv = newDataframe.to_csv(r"C:\Users\Richi\MEGA\Work\Python\eCLip\eCLIP_Indisulam_HCT116.csv", index = None, header=True)


#---OBSOLETE---
#---DATA MAP ON ECLIP---
#dict1 = pd.Series(df2.score.values, index=df2.gene).to_dict()
#df1["eClip"] = df1["geneName"].map(dict1)
#df1["eClip"].fillna(0, inplace=True)

#newDataframe.loc[:,"qValue_Condition2"] *= -1