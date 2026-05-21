
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# قائمة الجينات المناعية
immune_genes = ["CTLA4","PDCD1","CD274","STAT3","FOXP3","IL7R","IFNG","JAK1","JAK2","HLA-A","HLA-B","HLA-C","B2M"]

# قراءة ملف الطفرات (MAF)
df = pd.read_csv("data_mutations.txt", sep='\t', comment='#', low_memory=False)

# عمل mapping: أول 12 حرف من الـ Tumor_Sample_Barcode = Patient_ID
df["Patient_ID"] = df["Tumor_Sample_Barcode"].str[:12]

# قراءة ملف الـ clinical
clinical = pd.read_csv("data_clinical_patient.txt", sep='\t')

# تنظيف الأعمدة وتحويلها لأرقام
clinical = clinical.dropna(subset=["Overall Survival (Months)", "Overall Survival Status"])
clinical["OS_time"] = pd.to_numeric(clinical["Overall Survival (Months)"], errors="coerce")
clinical["OS_status"] = clinical["Overall Survival Status"].map({"DECEASED":1, "LIVING":0})

# جدول لتخزين النتائج
results_table = []

for gene in immune_genes:
    # المرضى المطفرين في الجين
    patients_with_gene = df[df["Hugo_Symbol"]==gene]["Patient_ID"].unique()
    
    # إضافة عمود يحدد حالة الطفرة
    clinical[gene+"_mut"] = clinical["#Patient Identifier"].isin(patients_with_gene)
    
    # تقسيم المجموعتين
    mask_mut = clinical[gene+"_mut"]==True
    mask_wt = clinical[gene+"_mut"]==False
    
    # لو فيه بيانات كافية (على الأقل مريض واحد في كل مجموعة)
    if mask_mut.sum()>0 and mask_wt.sum()>0:
        test = logrank_test(
            clinical.loc[mask_mut,"OS_time"], clinical.loc[mask_wt,"OS_time"],
            event_observed_A=clinical.loc[mask_mut,"OS_status"],
            event_observed_B=clinical.loc[mask_wt,"OS_status"]
        )
        pval = test.p_value
    else:
        pval = None
    
    results_table.append({
        "Gene":gene,
        "Mutated Patients":mask_mut.sum(),
        "Wildtype Patients":mask_wt.sum(),
        "Log-rank p-value":pval
    })

# تحويل النتائج إلى DataFrame مرتب
results_df = pd.DataFrame(results_table)
print(results_df.sort_values("Log-rank p-value"))

# رسم Barplot لعدد المرضى المطفرين
plt.figure(figsize=(8,6))
sns.barplot(
    x=results_df["Mutated Patients"],
    y=results_df["Gene"],
    hue=results_df["Gene"],
    palette="plasma",
    legend=False
)
plt.title("Number of Mutated Patients per Immune Gene (SKCM)")
plt.xlabel("Mutated Patients")
plt.ylabel("Immune Gene")
plt.tight_layout()
plt.show()
