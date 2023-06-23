library(gsveasyr)
#library(stringr)

env.all = read.table('Data/mag_env_for_shotgun_samples.csv',
                      header=TRUE, sep = ",")

# # Check if number of samples in one of the tables is not correct
# gene.list = read.table("Data/N_cycle/graftm_gene_counts.tsv", header = T)
# gene.list$SampleId = str_extract(gene.list$SampleId, "(?<=_)[A-Za-z0-9]+")
# setdiff(rownames(gene.list), rownames(env.all))


rownames(env.all) = env.all$SampleID

C_folder = 'Data/C_cycle/'
N_folder = 'Data/N_cycle/'
S_folder = 'Data/S_cycle/'

load_all_gsv_output(c_folder = C_folder, n_folder = N_folder, s_folder = S_folder,
                    min_num_reads = 1000000, clustlvls=c(65,75,85),
                    changeSampleName = T, env_df = env.all,
                    refColumn = env.all$SampleID)

bargraph.CI(GSV_gene_count_list$env$pH,
            GSV_gene_count_list$C_raw$ackA,
            legend = T, xlab = '', ylab = 'Gene count per million reads',
            main = 'aprA')
