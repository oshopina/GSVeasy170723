library(gsveasyr)
#library(stringr)

env.all = read.table('Data/mag_env_for_shotgun_samples.csv',
                      header=TRUE, sep = ",")

rownames(env.all) = env.all$SampleID

# # Check if number of samples in one of the tables is not correct
# gene.list = read.table("Data/N_cycle/graftm_gene_counts.tsv", header = T)
# gene.list$SampleId = str_extract(gene.list$SampleId, "(?<=_)[A-Za-z0-9]+")
# setdiff(rownames(gene.list), rownames(env.all))

C_folder = 'Data/C_cycle/'
N_folder = 'Data/N_cycle/'
S_folder = 'Data/S_cycle/'


load_all_gsv_output(C_folder, N_folder, S_folder,
                    min_num_reads = 1000000, clustlvls=c(65,75,85),
                    changeSampleName = T, env_df = env.all,
                    refColumn = env.all$gsveasy_sample)

GSV_gene_count_list$env$pH.group.label = as.factor(GSV_gene_count_list$env$pH.group.label)

auto_uniplot(GSV_gene_count_list$C_per_million, GSV_gene_count_list$env, "pH.group.label")

df.tmp = data_frame(aprA = GSV_gene_count_list$S_per_million$aprA,
                    pH_group = GSV_gene_count_list$env$pH.group.label)

posthoc.res = cld(lsmeans(aov(aprA ~ pH_group, data = df.tmp),
                          ~ pH_group, Letters = letters,
                          adjust = 'sidak'))

auto_barplot_posthoc(GSV_gene_count_list$C_per_million, 'pH.group.label',
                     GSV_gene_count_list$env)

Data = GSV_alpha_div_cluster65$C_sobs$ackA
anova(lm(Data ~ pH, data = env.all))

df = GSV_alpha_div_cluster65$C_sobs
aov_results = auto_aov_fixed(df, ~ pH, env.all)

adonis_res = vegan::adonis2(GSV_gene_count_list$C_per_million ~ pH,
                     data = GSV_gene_count_list$env,
                     method = 'euc', na.rm = T)

auto_adonis_fixed(GSV_rarefied_gsvtables65, ~ pH, env.all)
