library(gsveasyr)
library(vegan)
library(dplyr)
#library(stringr)

# Read the environmental data
env.all <- read.table('Data/mag_env_for_shotgun_samples.csv', header=TRUE, sep = ",")
rownames(env.all) <- env.all$SampleID

# # Check if number of samples in one of the tables is not correct
# gene.list = read.table("Data/N_cycle/graftm_gene_counts.tsv", header = T)
# gene.list$SampleId = str_extract(gene.list$SampleId, "(?<=_)[A-Za-z0-9]+")
# setdiff(rownames(gene.list), rownames(env.all))

# Define folder paths
C_folder <- 'Data/C_cycle/'
N_folder <- 'Data/N_cycle/'
S_folder <- 'Data/S_cycle/'

# Load GSV data and process
load_all_gsv_output(C_folder, N_folder, S_folder,
                    min_num_reads = 1000000, clustlvls=c(65,75,85),
                    changeSampleName = TRUE, env_df = env.all,
                    refColumn = env.all$gsveasy_sample)

# Prepare pH group label as a factor
GSV_gene_count_list$env$pH.group.label <- as.factor(GSV_gene_count_list$env$pH.group.label)

rm(C_folder, N_folder, S_folder)

## Alpha-diversity

# Calculate alpha diversity metrics
alphas <- list()
for (i in c('C_per_million', 'N_per_million', 'S_per_million')) {
  abundance_table <- GSV_gene_count_list[[i]]
  shannon <- diversity(abundance_table, index = "shannon")
  sobs <- rowSums(decostand(abundance_table, method = 'pa'))
  chao <- estimateR(round(abundance_table))[2,]
  alpha <- data.frame(shannon, sobs, chao, GSV_gene_count_list$env$pH, samples = names(shannon))
  colnames(alpha) <- c('shannon', 'sobs', 'chao', 'pH', 'samples')
  alphas[[i]] <- alpha
}

rm(abundance_table, alpha, chao, i, shannon, sobs)

library(ggplot2)
library(ggpmisc)

alpha_plots <- list()
for (i in c('C_per_million', 'N_per_million', 'S_per_million')) {
  df <- alphas[[i]]
  df_filtered <- na.omit(df)
  plots <- list()
  for (j in c('shannon', 'sobs', 'chao')) {
    model <- lm(paste(j, "~ pH"), data = df)
    p_value <- summary(model)$coefficients["pH", "Pr(>|t|)"]
    color <- ifelse(coef(model)["pH"] > 0, "red", "blue")
    if (p_value > 0.05) {
      color <- "gray"
    }
    # Create the scatter plot with trend line using the lm method
    plot <- ggplot(df_filtered, aes(x = pH, y = !!sym(j))) +
      geom_point() +
      geom_smooth(method = 'loess', color = color) +
      stat_poly_eq(use_label(c("R2", "p"))) +
      geom_point() +
      ggtitle(i)
    plots[[j]] <- plot
  }
  alpha_plots[[i]] <- plots
}
 rm(df, df_filtered, model, plot, plots, color, i, j, p_value)

# pdf("Figures/alpha_plots.pdf")
# for (plot in alpha_plots) {
#   print(plot)
# }
# dev.off()

## Gene abundance 
gene_plots <- list()
anova_results = list()
for (i in c('C_per_million', 'N_per_million', 'S_per_million')) {
  anova <- auto_aov_fixed(GSV_gene_count_list[[i]], ~ pH, GSV_gene_count_list$env)
  anova_results[[i]] = anova$Results
  genes_list <- unique(anova$Results$Data)
  plots <- list()
  for (j in genes_list) {
    df <- cbind(GSV_gene_count_list[[i]][[j]], GSV_gene_count_list$env$pH) %>% as.data.frame()
    model <- lm(V2 ~ V1, data = df)
    p_value <- summary(model)$coefficients["V1", "Pr(>|t|)"]
    color <- ifelse(coef(model)["V1"] > 0, "red", "blue")
    if (p_value > 0.05) {
      color <- "gray"
    }
    # Create the scatter plot with trend line using the lm method
    plot <- ggplot(df, aes(x = V2, y = V1)) +
      geom_point() +
      geom_smooth(method = "loess", color = color) +
      stat_poly_eq(use_label(c("R2", "p"))) +
      xlab('pH') +
      ylab(j) +
      ggtitle(i)
    plots[[j]] <- plot
  }
  gene_plots[[i]] <- plots
}

rm(anova, df, model, plot, plots, color, genes_list, i, j, p_value)

# pdf("Figures/gene_plots.pdf")
# for (plot in gene_plots) {
#   print(plot)
# }
# dev.off()

