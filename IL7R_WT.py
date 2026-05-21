
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



df = pd.read_csv("data_mutations.txt", sep='\t', comment='#', low_memory=False)





from scipy.stats import mannwhitneyu

# قراءة ملف الطفرات
mut_df = pd.read_csv("data_mutations.txt", sep="\t")

# تحديد العمود اللي فيه الـ Sample IDs
possible_cols = ["Tumor_Sample_Barcode", "Sample_ID", "PATIENT_ID", "Tumor_Sample_UUID"]
sample_col = None
for col in possible_cols:
    if col in mut_df.columns:
        sample_col = col
        break

if sample_col is None:
    raise ValueError("No patient/sample ID column found in mutation file")

# حساب عدد الطفرات لكل مريض (TMB)
tmb = mut_df.groupby(sample_col).size().reset_index(name="Mutation_Count")

# تحديد المرضى اللي عندهم طفرة في IL7R
il7r_mutants = mut_df[mut_df["Hugo_Symbol"]=="IL7R"][sample_col].unique()

# إضافة حالة IL7R لكل مريض
tmb["IL7R_status"] = tmb[sample_col].apply(lambda x: "Mutant" if x in il7r_mutants else "WT")

print("Tumor Mutation Burden (TMB) in IL7R Mutant vs WT",tmb)
# مقارنة التوزيع
sns.boxplot(x="IL7R_status", y="Mutation_Count", data=tmb, palette="Set2")
plt.title("Tumor Mutation Burden (TMB) in IL7R Mutant vs WT")
plt.show()

# اختبار إحصائي (Mann-Whitney U)
mutant_counts = tmb[tmb["IL7R_status"]=="Mutant"]["Mutation_Count"]
wt_counts = tmb[tmb["IL7R_status"]=="WT"]["Mutation_Count"]

stat, pval = mannwhitneyu(mutant_counts, wt_counts, alternative="two-sided")
print("Mann-Whitney U test p-value Tumor Mutation Burden (TMB) in IL7R Mutant vs WT:", pval)









import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. قراءة بيانات RNA-seq
expr = pd.read_csv("data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)

# 2. قراءة بيانات الطفرات
mut_df = pd.read_csv("data_mutations.txt", sep="\t")

# تحديد العمود اللي فيه الـ Sample IDs
possible_cols = ["Tumor_Sample_Barcode", "Sample_ID", "PATIENT_ID", "Tumor_Sample_UUID"]
sample_col = next((col for col in possible_cols if col in mut_df.columns), None)

# 3. تحديد المرضى IL7R mutant
il7r_mutants = mut_df[mut_df["Hugo_Symbol"]=="IL7R"][sample_col].unique()

# 4. تحويل الـ RNA-seq بحيث المرضى في الصفوف
expr_T = expr.T

# إضافة حالة IL7R لكل عينة
expr_T["IL7R_status"] = ["Mutant" if sample in il7r_mutants else "WT" for sample in expr_T.index]

# 5. الجينات المناعية المهمة
immune_genes = ["CD8A","GZMB","PRF1","IFNG","PDCD1","LAG3","TIGIT","CD274","CXCL9","CXCL10"]

# 6. استخراج التعبير للجينات المطلوبة فقط
subset = expr_T[immune_genes + ["IL7R_status"]]

# 7. حساب مجموع التعبير لكل مريض
subset["Immune_sum"] = subset[immune_genes].sum(axis=1)

# 8. اختيار أعلى 10 مرضى
top10 = subset.sort_values("Immune_sum", ascending=False).head(10)

print("Top 10 patients by immune gene expression:")
print(top10[["Immune_sum","IL7R_status"]])

# 9. رسم Heatmap لأعلى 10 مرضى
plt.figure(figsize=(10,6))
sns.heatmap(top10[immune_genes].T, cmap="viridis", cbar=True,
            xticklabels=top10.index, yticklabels=immune_genes)

plt.title("Top 10 Patients - Immune Gene Expression")
plt.xlabel("Patients")
plt.ylabel("Immune Genes")
plt.show()















from scipy.stats import mannwhitneyu
import numpy as np

# 1. قراءة بيانات RNA-seq
expr = pd.read_csv("data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)

# 2. قراءة بيانات الطفرات
mut_df = pd.read_csv("data_mutations.txt", sep="\t")

# تحديد العمود اللي فيه الـ Sample IDs
possible_cols = ["Tumor_Sample_Barcode", "Sample_ID", "PATIENT_ID", "Tumor_Sample_UUID"]
sample_col = None
for col in possible_cols:
    if col in mut_df.columns:
        sample_col = col
        break

# 3. تحديد المرضى IL7R mutant
il7r_mutants = mut_df[mut_df["Hugo_Symbol"]=="IL7R"][sample_col].unique()

# 4. تحويل الـ RNA-seq بحيث المرضى في الصفوف
expr_T = expr.T

# إضافة حالة IL7R لكل عينة
expr_T["IL7R_status"] = ["Mutant" if sample in il7r_mutants else "WT" for sample in expr_T.index]

# 5. حساب الـ scores
expr_T["CD8_score"] = expr_T["CD8A"] if "CD8A" in expr.index else np.nan
expr_T["CYT_score"] = np.sqrt(expr_T["GZMB"] * expr_T["PRF1"])

# 6. Boxplots للمقارنة
for score in ["CD8_score","CYT_score"]:
    sns.boxplot(x="IL7R_status", y=score, data=expr_T, palette="Set2")
    plt.title(f"{score} in IL7R Mutant vs WT")
    plt.show()

    # اختبار إحصائي
    mutant_vals = expr_T[expr_T["IL7R_status"]=="Mutant"][score].dropna()
    wt_vals = expr_T[expr_T["IL7R_status"]=="WT"][score].dropna()
    stat, pval = mannwhitneyu(mutant_vals, wt_vals, alternative="two-sided")
    print(f"{score}: Mann-Whitney p-value  Cytolytic score = {pval}")










import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# 1. قراءة بيانات RNA-seq
expr = pd.read_csv("data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)

# 2. قراءة بيانات الطفرات
mut_df = pd.read_csv("data_mutations.txt", sep="\t")

# تحديد العمود اللي فيه الـ Sample IDs
possible_cols = ["Tumor_Sample_Barcode", "Sample_ID", "PATIENT_ID", "Tumor_Sample_UUID"]
sample_col = next((col for col in possible_cols if col in mut_df.columns), None)

# 3. تحديد المرضى IL7R mutant
il7r_mutants = mut_df[mut_df["Hugo_Symbol"]=="IL7R"][sample_col].unique()

# 4. تحويل الـ RNA-seq بحيث المرضى في الصفوف
expr_T = expr.T

# إضافة حالة IL7R لكل عينة
expr_T["IL7R_status"] = ["Mutant" if sample in il7r_mutants else "WT" for sample in expr_T.index]

# 5. الجينات المطلوبة
genes = ["CXCL9","CXCL10"]

# 6. مقارنة التعبير لكل جين
for gene in genes:
    if gene in expr.index:
        sns.boxplot(x="IL7R_status", y=gene, data=expr_T, palette="Set2")
        plt.title(f"{gene} expression in IL7R Mutant vs WT")
        plt.show()

        # اختبار إحصائي
        mutant_vals = expr_T[expr_T["IL7R_status"]=="Mutant"][gene].dropna()
        wt_vals = expr_T[expr_T["IL7R_status"]=="WT"][gene].dropna()
        stat, pval = mannwhitneyu(mutant_vals, wt_vals, alternative="two-sided")
        print(f"{gene}: Mann-Whitney p-value  Chemokine Score= {pval}")











import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# 1. قراءة بيانات RNA-seq
expr = pd.read_csv("data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0)

# 2. قراءة بيانات الطفرات
mut_df = pd.read_csv("data_mutations.txt", sep="\t")

# تحديد العمود اللي فيه الـ Sample IDs
possible_cols = ["Tumor_Sample_Barcode", "Sample_ID", "PATIENT_ID", "Tumor_Sample_UUID"]
sample_col = next((col for col in possible_cols if col in mut_df.columns), None)

# 3. تحديد المرضى IL7R mutant
il7r_mutants = mut_df[mut_df["Hugo_Symbol"]=="IL7R"][sample_col].unique()

# 4. تحويل الـ RNA-seq بحيث المرضى في الصفوف
expr_T = expr.T

# إضافة حالة IL7R لكل عينة
expr_T["IL7R_status"] = ["Mutant" if sample in il7r_mutants else "WT" for sample in expr_T.index]

# 5. الجينات المطلوبة (immune checkpoints)
genes = ["CD274","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2"]

# 6. استخراج التعبير للجينات المطلوبة فقط
subset = expr_T[genes + ["IL7R_status"]]

# 7. رسم Heatmap مجمّع
plt.figure(figsize=(10,6))
sns.heatmap(subset.sort_values("IL7R_status")[genes].T, cmap="mako", cbar=True,
            xticklabels=False, yticklabels=genes)

plt.title("Immune Checkpoint Gene Expression (IL7R Mutant vs WT)")
plt.xlabel("Patients (sorted by IL7R status)")
plt.ylabel("Immune Checkpoints")
plt.show()

# 8. اختبار إحصائي لكل جين
for gene in genes:
    if gene in expr.index:
        mutant_vals = expr_T[expr_T["IL7R_status"]=="Mutant"][gene].dropna()
        wt_vals = expr_T[expr_T["IL7R_status"]=="WT"][gene].dropna()
        stat, pval = mannwhitneyu(mutant_vals, wt_vals, alternative="two-sided")
        print(f"{gene}: Mann-Whitney p-value = Checkpoint {pval}")










import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt

# 1. قراءة بيانات الطفرات
mut_df = pd.read_csv("data_mutations.txt", sep="\t")

# تحديد العمود اللي فيه الـ Sample IDs
possible_cols = ["Tumor_Sample_Barcode", "Sample_ID", "PATIENT_ID", "Tumor_Sample_UUID"]
sample_col = next((col for col in possible_cols if col in mut_df.columns), None)

# 2. تحديد المرضى IL7R mutant
il7r_mutants = mut_df[mut_df["Hugo_Symbol"]=="IL7R"][sample_col].unique()

# 3. حساب عدد الطفرات لكل مريض في BRAF, NRAS, NF1
target_genes = ["BRAF","NRAS","NF1"]
mut_counts = mut_df[mut_df["Hugo_Symbol"].isin(target_genes)].groupby([sample_col,"Hugo_Symbol"]).size().unstack(fill_value=0)

# 4. إضافة حالة IL7R لكل مريض
mut_counts["IL7R_status"] = ["Mutant" if sample in il7r_mutants else "WT" for sample in mut_counts.index]

# 5. Boxplots للمقارنة
for gene in target_genes:
    sns.boxplot(x="IL7R_status", y=gene, data=mut_counts, palette="Set2")
    plt.title(f"{gene} mutation count in IL7R Mutant vs WT")
    plt.show()

    # اختبار إحصائي
    mutant_vals = mut_counts[mut_counts["IL7R_status"]=="Mutant"][gene]
    wt_vals = mut_counts[mut_counts["IL7R_status"]=="WT"][gene]
    stat, pval = mannwhitneyu(mutant_vals, wt_vals, alternative="two-sided")
    print(f"{gene}: Mann-Whitney p-value = {pval}")
