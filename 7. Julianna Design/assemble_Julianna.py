import pandas as pd 
import numpy as np


# July 2025
# Assembling Julianna: Mouse Genome-wide Cas9 CRISPRko library 

# Designed with CRISPick: targeting GRCm38 with all protein coding genes recognized by Ensembl v.102 annotations at 4 guides/gene
# On-target efficacy score reflects RS3seq-Chen2013+RS3target
keepcols= ["sgRNA Sequence","Target Gene Symbol","Target Gene ID","Target Transcript","Reference Sequence","Strand of sgRNA","sgRNA Cut Position (1-based)","Exon Number","Aggregate CFD Score","Off-Target CFD100 Hits","Off-Target Tier I CFD100 Hits","On-Target Efficacy Score","Pick Order"]
julianna=pd.read_table("julianna-sgrna-designs.txt",usecols=keepcols)
julianna_negative_controls=pd.read_table("julianna_negative_controls-sgrna-designs.txt",usecols=keepcols+["Other Target Matches"])
julianna_negative_controls["Target Gene Symbol"]=julianna_negative_controls["Other Target Matches"]
julianna_negative_controls["Target Gene ID"]=julianna_negative_controls["Other Target Matches"]
julianna_negative_controls=julianna_negative_controls.drop("Other Target Matches",axis=1)


#exclude guides with potential for sufficient cumulative off-target activity for antiproliferative DNA damage response
julianna=julianna[(julianna["Aggregate CFD Score"]<=4.8)]


#identify the set of all Ensembl v.102 genes to which each guide maps 
guide_mappings=pd.read_csv("adhoc_sgRNA_disco_GRCm38_Ensembl_SpyoCas9Ko_strict.csv",usecols=["Target Sequence", "On-target Gene IDs","On-target Gene Symbols"]) 
guide_mappings["On-target Gene IDs"]=guide_mappings["On-target Gene IDs"].str.split(",")
guide_mappings["On-target Gene Symbols"]=guide_mappings["On-target Gene Symbols"].str.split(",")


#Julianna has 3 guides per gene, but providing guides with this picking scheme with up to 4 for those interested
julianna_quota4=julianna.copy()
julianna_quota4=julianna_quota4.drop("Pick Order",axis=1)

# Use of target local mode allows the selection of the same guide for multiple genes
# collapse guides targeting the same gene
julianna_quota4=julianna_quota4.groupby("sgRNA Sequence").agg(list).reset_index()

#also report any other genes to which the guide maps (even if not picked as a target for the other gene)
julianna_quota4=pd.merge(julianna_quota4,guide_mappings,left_on="sgRNA Sequence", right_on="Target Sequence")

julianna_quota4["Target Gene ID"]=julianna_quota4.apply(lambda x:x["Target Gene ID"]+[gene for gene in x["On-target Gene IDs"] if gene not in x["Target Gene ID"]],axis=1)
julianna_quota4["Target Gene ID"]=julianna_quota4["Target Gene ID"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Target Gene Symbol"]=julianna_quota4.apply(lambda x:x["Target Gene Symbol"]+[gene for gene in x["On-target Gene Symbols"] if gene not in x["Target Gene Symbol"]],axis=1)
julianna_quota4["Target Gene Symbol"]=julianna_quota4["Target Gene Symbol"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Target Transcript"]=julianna_quota4["Target Transcript"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Reference Sequence"]=julianna_quota4["Reference Sequence"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Strand of sgRNA"]=julianna_quota4["Strand of sgRNA"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["sgRNA Cut Position (1-based)"]=julianna_quota4["sgRNA Cut Position (1-based)"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Exon Number"]=julianna_quota4["Exon Number"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Aggregate CFD Score"]=julianna_quota4["Aggregate CFD Score"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Off-Target CFD100 Hits"]=julianna_quota4["Off-Target CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["Off-Target Tier I CFD100 Hits"]=julianna_quota4["Off-Target Tier I CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota4["On-Target Efficacy Score"]=julianna_quota4["On-Target Efficacy Score"].apply(lambda x:"|".join([str(i) for i in x]))

julianna_quota4=julianna_quota4.drop(["Target Sequence", "On-target Gene IDs","On-target Gene Symbols"],axis=1)

julianna_quota4=julianna_quota4.sort_values(by="Target Gene ID")
julianna_quota4=pd.concat([julianna_quota4,julianna_negative_controls])


julianna_quota4.to_csv("Julianna-GRCm38-Ensembl102-CRISPRkoSpCas9-Designs-Quota4.csv",index=False)

# Now generate report for Jacquere with 3 guides/gene 
julianna_quota3=julianna[(julianna["Pick Order"]!=4)]
julianna_quota3=julianna_quota3.drop("Pick Order",axis=1)

julianna_quota3=julianna_quota3.groupby("sgRNA Sequence").agg(list).reset_index()

julianna_quota3=pd.merge(julianna_quota3,guide_mappings,left_on="sgRNA Sequence", right_on="Target Sequence")

julianna_quota3["Target Gene Symbol"]=julianna_quota3.apply(lambda x:x["Target Gene Symbol"]+[gene for gene in x["On-target Gene Symbols"] if gene not in x["Target Gene Symbol"]],axis=1)
julianna_quota3["Target Gene Symbol"]=julianna_quota3["Target Gene Symbol"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Target Gene ID"]=julianna_quota3.apply(lambda x:x["Target Gene ID"]+[gene for gene in x["On-target Gene IDs"] if gene not in x["Target Gene ID"]],axis=1)
julianna_quota3["Target Gene ID"]=julianna_quota3["Target Gene ID"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Target Transcript"]=julianna_quota3["Target Transcript"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Reference Sequence"]=julianna_quota3["Reference Sequence"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Strand of sgRNA"]=julianna_quota3["Strand of sgRNA"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["sgRNA Cut Position (1-based)"]=julianna_quota3["sgRNA Cut Position (1-based)"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Exon Number"]=julianna_quota3["Exon Number"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Aggregate CFD Score"]=julianna_quota3["Aggregate CFD Score"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Off-Target CFD100 Hits"]=julianna_quota3["Off-Target CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["Off-Target Tier I CFD100 Hits"]=julianna_quota3["Off-Target Tier I CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
julianna_quota3["On-Target Efficacy Score"]=julianna_quota3["On-Target Efficacy Score"].apply(lambda x:"|".join([str(i) for i in x]))

julianna_quota3=julianna_quota3.drop(["Target Sequence", "On-target Gene IDs","On-target Gene Symbols"],axis=1)

julianna_quota3=julianna_quota3.sort_values(by="Target Gene ID")

julianna_quota3=pd.concat([julianna_quota3,julianna_negative_controls])

julianna_quota3.to_csv("Julianna-GRCm38-Ensembl102-CRISPRkoSpCas9-Designs-Quota3.csv",index=False)
