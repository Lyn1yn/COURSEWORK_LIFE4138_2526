#input data(put path of tsv files in brackets)
import pandas as pd
AD=pd.read_csv("/Users/linl/COURSEWORK_LIFE4138_2526/python_coursework/A_vs_D.deseq2.results.tsv", sep='\t')
GD=pd.read_csv("/Users/linl/COURSEWORK_LIFE4138_2526/python_coursework/G_vs_D.deseq2.results.tsv", sep='\t')
AD.head()
GD.head()

# Summary Statistics:
#number of significantly upregulated(log2FoldChange > 1) and downregulated(log2FoldChange < -1) genes(p<0.05)
#A vs D: upregulated: 365; down regulated: 653
#G vs D: upregulated: 833; down regulated: 1202
AD_up = AD[(AD['pvalue'] < 0.05) & (AD['log2FoldChange'] > 1)] #select significant upregulated data in AD
AD_down = AD[(AD['pvalue'] < 0.05) & (AD['log2FoldChange'] < -1)] #select significant downregulated data in AD

GD_up = GD[(GD['pvalue'] < 0.05) & (GD['log2FoldChange'] > 1)] #select significant upregulated data in GD
GD_down = GD[(GD['pvalue'] < 0.05) & (GD['log2FoldChange'] < -1)] #select significant downregulated data in GD

len(AD_up)
len(AD_down)
len(GD_up)
len(GD_down)
#summary 
AD[['log2FoldChange', 'pvalue']].describe() #Generate descriptive statistics
GD[['log2FoldChange', 'pvalue']].describe() #Generate descriptive statistics


# Plots:
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#data cleaning
AD.isnull().sum() #Check for missing values in AD
AD.dropna(subset=['pvalue', 'padj'], inplace=True) #Remove missing values in AD

GD.isnull().sum() #Check for missing values in GD
GD.dropna(subset=['pvalue', 'padj'], inplace=True) #Remove missing values in GD

#Volcano plot
#plot A vs D
AD['-log10padj'] = -np.log10(AD['padj']) #Convert the adjusted p-values ​​into more intuitive numerical values ​​to visualize significance.
plot1=sns.scatterplot(data=AD,x="log2FoldChange", y="-log10padj") #Use Seaborn to draw scatter plots
plt.title('The relationship between log2foldchange and adjusted p value in A vs D result')
plt.show()
#plot G vs D
GD['-log10padj'] = -np.log10(GD['padj']) #Convert the adjusted p-values ​​into more intuitive numerical values ​​to visualize significance.
plot2=sns.scatterplot(data=GD,x="log2FoldChange", y="-log10padj") #Use Seaborn to draw scatter plots
plt.title('The relationship between log2foldchange and adjusted p value in G vs D result')
plt.show()

#MA plot
#plot A vs D
from bioinfokit import visuz #Import Module
visuz.GeneExpression.ma(df=AD, #Input data
                        lfc='log2FoldChange',  #y axis
                        ct_count='baseMean', #x axis
                        st_count='baseMean', #x axis
                        pv='padj', #p value
                         title='the relationship between the log fold change and mean expression in A and D comparison',
                         axxlabel='mean expression(baseMean)', 
                         axylabel='log 2 fold change')
#plot G vs D
visuz.GeneExpression.ma(df=GD, lfc='log2FoldChange', ct_count='baseMean', st_count='baseMean', pv='padj',
                        title='the relationship between the log fold change and mean expression in G and D comparison',
                        axxlabel='mean expression(baseMean)', 
                        axylabel='log 2 fold change')

#Histogram
#plot A vs D
AD['pvalue'].plot.hist() #select the pvalue column from the A vs D, use Pandas' built-in plotting features to create a histogram.
plt.title('distribution of statistical significance in A vs D comparison')
plt.xlabel('p value')
plt.ylabel('frequency')
plt.show()
#plot G vs D
GD['pvalue'].plot.hist() #select the pvalue column from the G vs D, use Pandas' built-in plotting features to create a histogram.
plt.title('distribution of statistical significance in G vs D comparison')
plt.xlabel('p value')
plt.ylabel('frequency')
plt.show()

#Heatmap
#plot A vs D
top_gene1 = AD.sort_values('log2FoldChange', key=lambda x: x.abs(), ascending=False).head(30) #select the top30 expressed genes
top_gene_2d_1 = top_gene1.set_index('gene_id')[['log2FoldChange']] #select only gene id and log2fold change value to do plot
plot = sns.heatmap(top_gene_2d_1, cmap="viridis")
plt.title('The gene expression patterns in A vs D results')
plt.show()
#plot G vs D
top_gene2 = GD.sort_values('log2FoldChange', key=lambda x: x.abs(), ascending=False).head(30) #select the top30 expressed genes
top_gene_2d_2 = top_gene2.set_index('gene_id')[['log2FoldChange']] #select only gene id and log2fold change value to do plot
plot = sns.heatmap(top_gene_2d_2, cmap="viridis")
plt.title('The gene expression patterns in G vs D results')
plt.show()


#Tables:
table_AD = pd.concat([AD_up[['gene_id', 'log2FoldChange', 'pvalue', 'padj']], AD_down[['gene_id', 'log2FoldChange', 'pvalue', 'padj']]]) #combine upregulated and downregulated DataFrames
print(table_AD)
table_AD.to_csv("table_AD.csv", index=False) #Save the table as a CSV file


table_GD = pd.concat([GD_up[['gene_id', 'log2FoldChange', 'pvalue', 'padj']], GD_down[['gene_id', 'log2FoldChange', 'pvalue', 'padj']]]) #combine upregulated and downregulated DataFrames
print(table_GD)
table_GD.to_csv("table_GD.csv", index=False)


#Additional analysis:
#clustering
top_gene_2d_1index = top_gene1.set_index("gene_id")[["log2FoldChange"]]#make A vs D gene id into index, during the subsequent concat operation, genes with the same gene_id will be aligned.
top_gene_2d_2index = top_gene2.set_index("gene_id")[["log2FoldChange"]]#make G vs D gene id into index, during the subsequent concat operation, genes with the same gene_id will be aligned.
common_genes=pd.concat([top_gene_2d_1index,top_gene_2d_2index], axis=1, join="inner")#combine the same gene in A vs D and G vs D together
common_genes.columns = ["AvsD", "GvsD"]#change the column name
sns.clustermap(common_genes, cmap="vlag")
plt.title("Clustering of top expressed genes (based on log2FC)")
plt.show()

#potential pathway (Gene Ontology Biological Processes)
import gseapy as gp
gp.enrichr(
    gene_list=AD,
    gene_sets="GO_Biological_Process_2023",
    organism="Yeast",
)
