import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv("data_mutations.txt", sep='\t', comment='#', low_memory=False)

# نطلع أعلى 10 جينات حصل فيها طفرات في المرضى
top_genes = df['Hugo_Symbol'].value_counts().head(10)
print("Top 10 mutated genes in skin cancer",top_genes)

# رسم بياني سريع بشكل شيك
plt.figure(figsize=(10, 5))
sns.barplot(x=top_genes.values, y=top_genes.index, palette="viridis")
plt.title('Top 10 Mutated Genes in Skin Cutaneous Melanoma (TCGA)')
plt.xlabel('Number of Mutations')
plt.ylabel('Gene')
plt.tight_layout()
plt.show()








immune_genes = ["CTLA4","PDCD1","CD274","STAT3","FOXP3","IL7R","IFNG","JAK1","JAK2","HLA-A","HLA-B","HLA-C","B2M"]

immune_counts = df[df['Hugo_Symbol'].isin(immune_genes)]['Hugo_Symbol'].value_counts()
print("Mutation Counts in Immune Genes (SKCM)",immune_counts)

plt.figure(figsize=(8,5))
sns.barplot(x=immune_counts.values, y=immune_counts.index, palette="magma")
plt.title("Mutation Counts in Immune Genes (SKCM)")
plt.xlabel("Number of Mutations")
plt.ylabel("Immune Gene")
plt.tight_layout()
plt.show()








variant_counts = df[df["Hugo_Symbol"].isin(immune_genes)]["Variant_Classification"].value_counts()


print("Distribution of mutation types in immune genes:", variant_counts)

# رسم بياني
plt.figure(figsize=(8,5))
sns.barplot(x=variant_counts.values, y=variant_counts.index, palette="coolwarm")
plt.title("Variant Classification in Immune Genes (SKCM)")
plt.xlabel("Number of Mutations")
plt.ylabel("Variant Type")
plt.tight_layout()
plt.show()








# حساب عدد المرضى المطفرين لكل جين
mutated_patients = df.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique()

# فلترة على الجينات المناعية فقط
immune_mutated_patients = mutated_patients[mutated_patients.index.isin(immune_genes)]

print("Number of mutated patients per immune gene:")
print(immune_mutated_patients.sort_values(ascending=False))









# إجمالي عدد المرضى في الكوهورت (مثلاً 470 في TCGA-SKCM)
total_patients = df["Tumor_Sample_Barcode"].nunique()

# عدد المرضى المطفرين في كل جين
mutated_patients = df.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique()

# فلترة على الجينات المناعية فقط
immune_mutated_patients = mutated_patients[mutated_patients.index.isin(immune_genes)]

# حساب النسبة المئوية
mutation_freq = (immune_mutated_patients / total_patients) * 100

# جدول مرتب
mutation_table = mutation_freq.sort_values(ascending=False).round(2)

print("Mutation Frequency (%) per immune gene:")
print(mutation_table)

# رسم بياني
plt.figure(figsize=(8,6))
sns.barplot(x=mutation_table.values, y=mutation_table.index, palette="plasma" , legend=False  )
plt.title("Mutation Frequency (%) in Immune Genes (SKCM)")
plt.xlabel("Frequency (%)")
plt.ylabel("Immune Gene")
plt.tight_layout()
plt.show()









antigen_genes = ["HLA-A","HLA-B","HLA-C","B2M"]
cytokine_genes = ["JAK1","JAK2","STAT3","IFNG"]

antigen_count = df[df['Hugo_Symbol'].isin(antigen_genes)].shape[0]
cytokine_count = df[df['Hugo_Symbol'].isin(cytokine_genes)].shape[0]

print("Antigen presentation mutations:", antigen_count)
print("Cytokine signaling mutations:", cytokine_count)





















