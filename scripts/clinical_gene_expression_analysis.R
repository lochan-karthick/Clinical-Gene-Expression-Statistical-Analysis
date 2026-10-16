#importing libraries:
library(ggplot2)
library(tibble)

data_dir <- Sys.getenv("DATA_DIR", unset = "data")

#loading the tables:

sample_info_table = read.table(file.path(data_dir, "Sample_Information.csv"), header = TRUE, 
                               row.names = 1, sep = "\t")
DE_SA_vs_HC_table = read.table(file.path(data_dir, "DE_SA_vs_HC.csv"), header = TRUE, 
                         row.names = 1, sep = "\t")
DE_GOUT_vs_HC_table = read.table(file.path(data_dir, "DE_GOUT_vs_HC.csv"), header = TRUE,
                           row.names = 1, sep = "\t")
expression_table = read.table(file.path(data_dir, "Expression_Table.csv"), header = TRUE, 
                              row.names = 1, sep = "\t")
annotations_table = read.table(file.path(data_dir, "Annotations.csv"), header = TRUE,
                               row.names = 1, sep = "\t")


#1.Clinical Matching:

summary(sample_info_table)

neutrophils_summary = aggregate(NEUTROPHILS ~ SAMPLE_GROUP, data = sample_info_table, 
                                 FUN = function(x) c(mean = mean(x), sd = sd(x)))

sex_summary = aggregate(SEX ~ SAMPLE_GROUP, data = sample_info_table, 
                         FUN = function(x) sum(x == "M") / length(x) * 100)

colnames(sex_summary) = c("SAMPLE_GROUP", "Percentage_Male")

summary_table = merge(sex_summary, neutrophils_summary, by = "SAMPLE_GROUP")

# Separate the aggregated neutrophil statistics into columns
summary_table = data.frame(
  Group = summary_table$SAMPLE_GROUP,
  Percentage_Male = summary_table$Percentage_Male,
  Neutrophils_Mean = summary_table$NEUTROPHILS[, "mean"],
  Neutrophils_SD = summary_table$NEUTROPHILS[, "sd"]
)

print(summary_table)

#To check if clinical data is well matched:
#We test the categorical variables which in our set is only sex
# Using Fisher's exact test for sex distribution

fisher_test_sex = fisher.test(table(sample_info_table$SEX, sample_info_table$SAMPLE_GROUP))
print(fisher_test_sex)
fisher_test_sex$p.value

#We test the continuous variable which in our set is only neutrophils
#Using the anova test to check the continuous variable neutrophils:

anova_neutrophil = aov(NEUTROPHILS ~ SAMPLE_GROUP, data = sample_info_table)
summary(anova_neutrophil) #we obtained a p-value of 3.25e-05 which implies there is a significant statistical difference in neutrophil counts among the groups




#2.Exploratory data analysis:

#To see how the grouped genes are clustered we use Principle component analysis:

#log transforming the expression data to manage variance:

expression_data = log2(expression_table + 1)
pca = prcomp(t(expression_data))
pca_df = data.frame(pca$x)
pca_df$SAMPLE_GROUP = sample_info_table$SAMPLE_GROUP

#making a scatterplot for the principle components obtained from the above step:

ggp=ggplot(pca_df, aes(x=PC1, y=PC2, color=SAMPLE_GROUP)) + 
    geom_point(size=3) +
    labs(x = "PC1", y = "PC2", title = "Principle components of Expression Data") 
# Plot saved below
#3. Analyzing and Visualizing the Differential Expression data of HC vs Gout and HC vs SA:

#selecting only the significantly expressed genes based off the adjusted p-value:

sig_genes_hc_gout = subset(DE_GOUT_vs_HC_table, p.adj<0.05)
sig_genes_hc_SA = subset(DE_SA_vs_HC_table, p.adj<0.05)

#scatterplot for HC vs gout:

ggp = ggplot(sig_genes_hc_gout, aes(x = log2Fold, y = -log10(p.adj))) + 
  geom_point(alpha = 0.5) +
  labs(x = "Log2Fold change", y = "-log 10 p.adj", title = "scatterplot: HC vs Gout")
# Plot saved below
#scatterplot for HC vs SA

ggp = ggplot(sig_genes_hc_SA, aes(x = log2Fold, y = -log10(p.adj))) + 
  geom_point(alpha = 0.5) +
  labs(x = "Log2Fold change", y = "-log 10 p.adj", title = "scatterplot: HC vs SA")
# Plot saved below
#4. Finding Correlations between the differential expression and clinical measurements:

sig_genes_hc_gout = subset(DE_GOUT_vs_HC_table, p.adj<0.05)
sig_genes_hc_SA = subset(DE_SA_vs_HC_table, p.adj<0.05)

#arranging the significant genes so the top 10 are at the top and bringing the gene_id's into the dataframe as a column

sig_genes_hc_gout = sig_genes_hc_gout[order(sig_genes_hc_gout$p.adj), ]
sig_genes_hc_SA = sig_genes_hc_SA[order(sig_genes_hc_SA$p.adj), ]

sig_genes_hc_gout = rownames_to_column(sig_genes_hc_gout, var = "Gene")
sig_genes_hc_SA = rownames_to_column(sig_genes_hc_SA, var = "Gene")

#selecting the top 10 genes:
top_genes_gout_HC = sig_genes_hc_gout$"Gene"[1:10]
top_genes_SA_HC = sig_genes_hc_SA$"Gene"[1:10]

#Correlating between expression and the other variables:

#Preparing the data for testing:

#obtaining expression data for the top ten HC vs Gout genes
top_genes_data_gout = intersect(row.names(expression_table),top_genes_gout_HC)

#creating an expression table of only the top 10 significant genes in gout
filtered_expression_table_gout = expression_table[top_genes_data_gout,]

#transposing the filtered table to get the gene_id as columns
transposed_filtered_expression_gout = t(filtered_expression_table_gout)

#merging the filtered table with the sample_info table to compare expression with the other variables:
transposed_filtered_expression_gout_correlation = merge(transposed_filtered_expression_gout, sample_info_table, by.x = 0, by.y = 0)
names(transposed_filtered_expression_gout_correlation)[1] = "Sample_names" 

#obtaining expression data for the top ten HC vs SA genes
top_genes_data_SA = intersect(row.names(expression_table),top_genes_SA_HC)

#creating an expression table of only the top 10 significant genes in SA
filtered_expression_table_SA = expression_table[top_genes_data_SA,]

#transposing the filtered table to get the gene_id as columns
transposed_filtered_expression_SA = t(filtered_expression_table_SA)

#merging the filtered table with the sample_info table to compare expression with the other variables:
transposed_filtered_expression_SA_correlation = merge(transposed_filtered_expression_SA, sample_info_table, by.x = 0, by.y = 0)
names(transposed_filtered_expression_SA_correlation)[1] = "Sample_names"

#Conducting t-test for categorical variable sex against the top 10 gout genes:

#creating an empty data frame for the t-test results:
t_test_results = data.frame(Gene = character(),
                         t_statistic = numeric(),
                         p = numeric())

#iterating over the top 10 significant genes in gout:
for (gene in top_genes_gout_HC) {
  
  #conducting t-tests between the gene expression and gender:
  t_test = t.test(transposed_filtered_expression_gout_correlation[[gene]]~transposed_filtered_expression_gout_correlation$SEX)
  
  #appending the results:
  t_test_results = rbind(t_test_results, data.frame(
    Gene = gene,
    t_statistic = unname(t_test$statistic),
    p = t_test$p.value))
  
}
print(t_test_results)#p-value for all genes across sex are above 0.005 implying that sex does not play a role.

#Conducting correlation tests to find relationships between the top 10 gout genes against neutrophil counts:

cor_results = data.frame(Gene = character(),
                            Correlation = numeric(),
                            p = numeric())

#iterating over the top 10 significant genes in gout:
for (gene in top_genes_gout_HC) {
  
  #conducting t-tests between the gene expression and gender:
  cor_test = cor.test(transposed_filtered_expression_gout_correlation[[gene]], transposed_filtered_expression_gout_correlation$NEUTROPHILS, method = "pearson")
  
  #appending the results:
  cor_results = rbind(cor_results, data.frame(
    Gene = gene,
    Correlation = cor_test$estimate,
    p = cor_test$p.value))
  
}
print(cor_results)

#Conducting t-test for categorical variable sex against top 10 SA genes to check for relationships:

#creating an empty data frame for the t-test results:
t_test_results = data.frame(Gene = character(),
                            t_statistic = numeric(),
                            p = numeric())

#iterating over the top 10 significant genes in SA:
for (gene in top_genes_SA_HC) {
  
  #conducting t-tests between the gene expression and gender:
  t_test = t.test(transposed_filtered_expression_SA_correlation[[gene]]~transposed_filtered_expression_SA_correlation$SEX)
  
  #appending the results:
  t_test_results = rbind(t_test_results, data.frame(
    Gene = gene,
    t_statistic = unname(t_test$statistic),
    p = t_test$p.value))
  
}
print(t_test_results)#p-value for all genes across sex are above 0.005 implying that sex does not play a role.

#Conducting correlation tests to find relationships between the top 10 SA genes against neutrophil counts:

cor_results = data.frame(Gene = character(),
                         Correlation = numeric(),
                         p = numeric())

#iterating over the top 10 significant genes in SA:
for (gene in top_genes_SA_HC) {
  
  #conducting t-tests between the gene expression and gender:
  cor_test = cor.test(transposed_filtered_expression_SA_correlation[[gene]], transposed_filtered_expression_SA_correlation$NEUTROPHILS, method = "pearson")
  
  #appending the results:
  cor_results = rbind(cor_results, data.frame(
    Gene = gene,
    Correlation = cor_test$estimate,
    p = cor_test$p.value))
  
}
print(cor_results) 
    

#5. Differentiating Gout from SA
#identifying genes that differ significantly between Gout and SA:

#preparing the data:

#subsetting the expression data for only Gout and SA:

gout_SA_samples = rownames(sample_info_table[sample_info_table$SAMPLE_GROUP %in% c("GOUT","SEPSIS"), ])
expression_gout_SA = log2(expression_table[, gout_SA_samples, drop = FALSE] + 1)

#obtaining the labels:
sample_groups = sample_info_table[gout_SA_samples,"SAMPLE_GROUP"]
gout_samples = rownames(sample_info_table[sample_info_table$SAMPLE_GROUP == "GOUT", ])
sepsis_samples = rownames(sample_info_table[sample_info_table$SAMPLE_GROUP == "SEPSIS", ])
#Performing the tests:

#vector to store p-values:
p_values = numeric(nrow(expression_gout_SA))

#looping through every gene and performing a t-test between gout and SA:
for(i in 1:nrow(expression_gout_SA)){
  
  #getting the expression data for the genes across gout and SA samples
  gout_expr = as.numeric(expression_gout_SA[i, gout_samples])
  SA_expr = as.numeric(expression_gout_SA[i, sepsis_samples])
  
  #performing t_tests between gout and sepsis groups for each gene
  test_result = t.test(gout_expr, SA_expr)
  p_values[i] = test_result$p.value
}
cat("Genes tested:", length(p_values), "\n")

#adjusting the p_value:
adjusted_p_values = p.adjust(p_values, method = "BH")

#identifying the significant genes:
significant_genes_Gout_SA = rownames(expression_gout_SA)[which(!is.na(adjusted_p_values) & adjusted_p_values < 0.05)]


# 6. Save selected figures

dir.create("results", showWarnings = FALSE, recursive = TRUE)

pca_plot = ggplot(pca_df, aes(PC1, PC2, color = SAMPLE_GROUP)) +
  geom_point(size = 3) +
  labs(title = "PCA of gene expression",
       x = "PC1", y = "PC2", color = "Group") +
  theme_minimal()

ggsave("results/pca.png", pca_plot, width = 8, height = 5, dpi = 200)

save_de_plot = function(data, title, filename) {
  selected = subset(data,
    !is.na(p.adj) & p.adj > 0 & p.adj < 0.05 &
    is.finite(log2Fold))

  plot = ggplot(selected, aes(log2Fold, -log10(p.adj))) +
    geom_point(alpha = 0.5, size = 1) +
    labs(title = title,
         x = "Log2 fold change",
         y = "-log10 adjusted p-value") +
    theme_minimal()

  ggsave(filename, plot, width = 8, height = 5, dpi = 200)
}

save_de_plot(
  DE_GOUT_vs_HC_table,
  "Significant genes: GOUT vs HC",
  "results/gout_vs_hc.png"
)

save_de_plot(
  DE_SA_vs_HC_table,
  "Significant genes: SEPSIS vs HC",
  "results/sepsis_vs_hc.png"
)

# Visualise the six genes with the lowest adjusted p-values
valid = which(!is.na(adjusted_p_values) &
              adjusted_p_values < 0.05)

if (length(valid) > 0) {
  top = head(valid[order(adjusted_p_values[valid])], 6)
  selected = expression_gout_SA[top, , drop = FALSE]

  plot_data = data.frame(
    Gene = rep(rownames(selected), times = ncol(selected)),
    Sample = rep(colnames(selected), each = nrow(selected)),
    Expression = as.vector(as.matrix(selected))
  )

  plot_data$Sample_Group =
    sample_info_table[plot_data$Sample, "SAMPLE_GROUP"]

  plot = ggplot(
    plot_data,
    aes(Sample_Group, Expression, fill = Sample_Group)
  ) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.5) +
    facet_wrap(~Gene, scales = "free_y") +
    labs(
      title = "Top genes: GOUT vs SEPSIS",
      x = "Sample group",
      y = "Log2(expression + 1)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  ggsave(
    "results/gout_vs_sepsis_top_genes.png",
    plot, width = 11, height = 7, dpi = 200
  )
}

cat("\nFigures saved to results/\n")
