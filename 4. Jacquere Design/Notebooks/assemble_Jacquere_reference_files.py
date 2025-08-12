import pandas as pd 

# Assembling jacquere: Human Genome-wide Cas9 CRISPRko library 

# Designed with CRISPick: targeting GRCh38 with all protein coding genes recognized by Ensembl v.113 annotations at 4 guides/gene
# On-target efficacy score reflects RS3seq-Chen2013+RS3target
keepcols= ["Target Gene ID","Target Gene Symbol","Target Transcript","Reference Sequence","Strand of sgRNA","sgRNA Cut Position (1-based)","sgRNA Sequence","Exon Number","Aggregate CFD Score","Off-Target CFD100 Hits","Off-Target Tier I CFD100 Hits","On-Target Efficacy Score","Pick Order"]
jacquere=pd.read_csv("../Data/jacquere_assembled_crispick.csv",usecols=keepcols)

#exclude guides with potential for sufficient cumulative off-target activity for antiproliferative DNA damage response
jacquere=jacquere[(jacquere["Aggregate CFD Score"]<=4.8)]

#jacquere has 3 guides per gene, but providing guides with this picking scheme with up to 4 for those interested
jacquere_quota4=jacquere.copy()

# Use of target local mode allows the selection of the same guide for multiple genes
# collapse guides targeting the same gene
jacquere_quota4=jacquere_quota4.groupby("sgRNA Sequence").agg(list).reset_index()


jacquere_quota4["Target Gene ID"]=jacquere_quota4["Target Gene ID"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Target Gene Symbol"]=jacquere_quota4["Target Gene Symbol"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Target Transcript"]=jacquere_quota4["Target Transcript"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Reference Sequence"]=jacquere_quota4["Reference Sequence"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Strand of sgRNA"]=jacquere_quota4["Strand of sgRNA"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["sgRNA Cut Position (1-based)"]=jacquere_quota4["sgRNA Cut Position (1-based)"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Exon Number"]=jacquere_quota4["Exon Number"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Aggregate CFD Score"]=jacquere_quota4["Aggregate CFD Score"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Off-Target CFD100 Hits"]=jacquere_quota4["Off-Target CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Off-Target Tier I CFD100 Hits"]=jacquere_quota4["Off-Target Tier I CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["On-Target Efficacy Score"]=jacquere_quota4["On-Target Efficacy Score"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota4["Pick Order"]=jacquere_quota4["Pick Order"].apply(lambda x:"|".join([str(i) for i in x]))

jacquere_quota4=jacquere_quota4.sort_values(by="Target Gene Symbol")

jacquere_quota4.to_csv("../Figures/jacquere-GRCh38-Ensembl113-CRISPRkoSpCas9-Designs-Quota4.csv",index=False)

jacquere_quota3=jacquere[(jacquere["Pick Order"]<4)]

jacquere_quota3=jacquere_quota3.groupby("sgRNA Sequence").agg(list).reset_index()

jacquere_quota3["Target Gene ID"]=jacquere_quota3["Target Gene ID"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Target Gene Symbol"]=jacquere_quota3["Target Gene Symbol"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Target Transcript"]=jacquere_quota3["Target Transcript"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Reference Sequence"]=jacquere_quota3["Reference Sequence"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Strand of sgRNA"]=jacquere_quota3["Strand of sgRNA"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["sgRNA Cut Position (1-based)"]=jacquere_quota3["sgRNA Cut Position (1-based)"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Exon Number"]=jacquere_quota3["Exon Number"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Aggregate CFD Score"]=jacquere_quota3["Aggregate CFD Score"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Off-Target CFD100 Hits"]=jacquere_quota3["Off-Target CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Off-Target Tier I CFD100 Hits"]=jacquere_quota3["Off-Target Tier I CFD100 Hits"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["On-Target Efficacy Score"]=jacquere_quota3["On-Target Efficacy Score"].apply(lambda x:"|".join([str(i) for i in x]))
jacquere_quota3["Pick Order"]=jacquere_quota3["Pick Order"].apply(lambda x:"|".join([str(i) for i in x]))

jacquere_quota3=jacquere_quota3.sort_values(by="Target Gene Symbol")

jacquere_quota3.to_csv("../Figures/jacquere-GRCh38-Ensembl113-CRISPRkoSpCas9-Designs-Quota3.csv",index=False)

