import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



df = pd.read_csv("data_mutations.txt", sep='\t', comment='#', low_memory=False)






# فلترة الطفرات الخاصة بجين IL7R
il7r_mutations = df[df["Hugo_Symbol"]=="IL7R"]["Variant_Classification"].value_counts()

print("Mutation Spectrum only for the IL7R gene",il7r_mutations)

# رسم Bar Plot يوضح أنواع الطفرات في IL7R
il7r_mutations.plot(
    kind="bar",
    figsize=(8,5),
    color="skyblue"
)

plt.title("Mutation Spectrum for IL7R (SKCM)")
plt.xlabel("Variant Classification")
plt.ylabel("Number of Mutations")
plt.tight_layout()
plt.show()











# مثال بيانات الطفرات لكل مريض
mut_matrix = pd.DataFrame({
    "Patient": ["Patient1","Patient2","Patient3","Patient4"],
    "IL7R": [1,0,1,0],
    "JAK1": [0,1,0,1],
    "STAT3": [1,0,0,1]
}).set_index("Patient")

# حساب التشارك (co-mutation)
co_mutations = {}
for gene in ["JAK1","STAT3"]:
    co_mutations[f"IL7R + {gene}"] = ((mut_matrix["IL7R"]==1) & (mut_matrix[gene]==1)).sum()

# عرض النتائج
co_df = pd.DataFrame(list(co_mutations.items()), columns=["Gene Pair","Patients"])
print("Co-mutation Analysis",co_df)










# عمل mapping: أول 12 حرف من الـ Tumor_Sample_Barcode = Patient_ID
df["Patient_ID"] = df["Tumor_Sample_Barcode"].str[:12]

# استخراج المرضى اللي عندهم طفرة في IL7R
il7r_patients = df[df["Hugo_Symbol"]=="IL7R"]["Patient_ID"].unique()

# قائمة الجينات اللي عايز تعمل لها Co-mutation
genes_of_interest = ["BRAF","NRAS","TP53"]

# جدول لتخزين النتائج
co_mutation_results = []

for gene in genes_of_interest:
    patients_with_gene = df[df["Hugo_Symbol"]==gene]["Patient_ID"].unique()
    
    overlap = len(set(il7r_patients) & set(patients_with_gene))
    il7r_only = len(il7r_patients) - overlap
    gene_only = len(patients_with_gene) - overlap
    
    co_mutation_results.append({
        "Gene": gene,
        "IL7R+Gene": overlap,
        "IL7R_only": il7r_only,
        "Gene_only": gene_only
    })

# تحويل النتائج إلى DataFrame مرتب
co_mutation_df = pd.DataFrame(co_mutation_results)
print(co_mutation_df)













import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency

# عمل mapping للـ Patient_ID
df["Patient_ID"] = df["Tumor_Sample_Barcode"].str[:12]

# استخراج المرضى المطفرين في IL7R و BRAF
il7r_patients = set(df[df["Hugo_Symbol"]=="IL7R"]["Patient_ID"].unique())
braf_patients = set(df[df["Hugo_Symbol"]=="BRAF"]["Patient_ID"].unique())

# إجمالي عدد المرضى
total_patients = df["Patient_ID"].nunique()

# حساب الأعداد
IL7R_and_BRAF = len(il7r_patients & braf_patients)
IL7R_only = len(il7r_patients - braf_patients)
BRAF_only = len(braf_patients - il7r_patients)
Neither = total_patients - (IL7R_and_BRAF + IL7R_only + BRAF_only)

# بناء جدول 2x2
contingency_table = [
    [IL7R_and_BRAF, IL7R_only],
    [BRAF_only, Neither]
]

print("Contingency Table (IL7R vs BRAF):")
print(pd.DataFrame(contingency_table,
                   index=["BRAF+","BRAF-"],
                   columns=["IL7R+","IL7R-"]))

# Fisher Exact Test
oddsratio, pval_fisher = fisher_exact(contingency_table)
print("Fisher Exact Test p-value:", pval_fisher)

# Chi-square Test
chi2, pval_chi, dof, expected = chi2_contingency(contingency_table)
print("Chi-square Test p-value:", pval_chi)


















data = {
    "Patient": ["Patient1", "Patient2", "Patient3"],
    "IL7R": [1, 0, 1],
    "JAK1": [0, 1, 0],
    "JAK2": [0, 0, 1],
    "STAT3": [1, 0, 0]
}
df = pd.DataFrame(data).set_index("Patient")

print("Oncoplot",df)

# Heatmap
plt.figure(figsize=(6,4))
sns.heatmap(df, cmap="Reds", cbar=False, linewidths=0.5, linecolor="black")
plt.title("Oncoplot (Immune Gene Mutations)")
plt.show()