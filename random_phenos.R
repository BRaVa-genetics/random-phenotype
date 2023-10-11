#!/usr/bin/env Rscript

# Load required packages
packages = c('dplyr', 'optparse', 'tidyverse', 'Matrix', 'sparseMVN', 'spdep', 'qlcMatrix', 'ggplot2', 'stringr')
for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p, repos = "http://cran.us.r-project.org")
  }
}

# Define input arguments
option_list = list(
  make_option(c("-g", "--grm"), type="character", default=NULL,
              help="path to the sparse GRM (.mtx file)", metavar="character"),
  make_option(c("--heritabilities"), type="character", default='0.5',
              help="heritability of the phenotypes (a list of comma split values)", metavar="character"),
  make_option(c("-n", "--n"), type="integer", default=1000,
              help="number of raw phenotypes needed", metavar="integer"),
  make_option(c("-p", "--prevalence"), type="character", default='100,300,1000',
              help="number of cases (a list of comma split integers) or prevalences (a list of comma split value between 0 and 1)", metavar="character"),
  make_option(c("-r", "--n_reps"), type="integer", default=1,
              help="number of draws to run", metavar="integer"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="path to the sample ID file", metavar="character"),
  make_option(c("-i", "--seed"), type="double", default=663,
              help="set seed for random sampling", metavar="double"),
  make_option(c("-t", "--trait_type"), type="character", default='both',
              help="the type of phenotypes to simulate, (choose from `continuous`, `binary` or `both`)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="random_phenos",
              help="output file name [default= %default]", metavar="character")
)
opt_parser = OptionParser(add_help_option=TRUE, option_list=option_list)
opt = parse_args(opt_parser);

# Confirm that a valid GRM is provided
if (is.null(opt$grm)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Define variables
data = readMM(opt$grm) # GRM matrix
samples = read.table(opt$samples) # a list of sample ID
n_samples = dim(data)[1] # Number of samples
print(paste0('N samples:', n_samples))

n_phenos_to_simulate = opt$n # Number of raw phenotypes to simulate
heritabilities = as.numeric(unlist(strsplit(opt$heritabilities, ','))) # List of heritabilities to use for the random phenotypes
n_heritabilities = length(heritabilities) # Number of heritabilities to use for the random phenotypes
binary_info = as.numeric(unlist(strsplit(opt$prevalence, ',')))
if(mean(binary_info) > 1){
  n_cases = binary_info # List of n_cases to use for the random phenotypes
  prevalences = n_cases/n_samples # Compute the prevalences from the n_cases
}else{
  prevalences = binary_info
  n_cases = prevalences*n_samples
}
n_prevalences = length(prevalences) # Number of prevalences to use for the random phenotypes
qnorms = qnorm(prevalences) # Compute the corresponding quantile of each prevalence in a standard normal distribution
n_draws = opt$n_reps # Number of draws to run per heritbility per n_case
n_phenos = ((n_prevalences + 1)*n_draws) # Number of phenotypes to draw per heritability

# Compute the Cholesky decomposition of the GRM
decomp = Cholesky(data)

# Build `n_phenos_to_simulate` raw random phenotypes from multivariate normal distributions
# with mean 0 and variance information from decomposition computed in the previous step
set.seed(opt$seed)
phenos = rmvn.sparse(n_phenos_to_simulate, rep(0, n_samples), decomp, prec=FALSE)

# Pivot the values in GRM to a data.frame
orig = data.frame(x=data@i + 1, y=data@j + 1, orig=data@x) %>% filter(x != y)

# Randomly sample a same number of uncorrelated pairs of samples (genetic correlation = 0)
# from a uniform distribution with min 0 and max = `n_samples`
# and merge it with the original data.frame
set.seed(opt$seed)
orig_uncorr = data.frame(x=round(runif(nrow(orig), max=n_samples)),
                         y=round(runif(nrow(orig), max=n_samples)),
                         orig=0)
orig = orig %>%
  union_all(orig_uncorr %>%
              anti_join(orig, by=c('x', 'y'))) %>%
  filter(x > 0 & y > 0)

# Annotate the corresponding phenotypic correlation from the previously simulated phenotypes
output = orig %>%
  rowwise() %>%
  mutate(pheno_corr=cor(phenos[,x], phenos[,y]))

# Plot and save a figure of the relationship between genetic and phenotypic correlation
output %>%
  ggplot + aes(x = orig, y = pheno_corr) +
  geom_point() +
  labs(x = 'GRM', y = 'Random pheno correlation')
ggsave(paste0(opt$out, '.png'))

# Create a data.frame with only userIds
output_df = samples %>%
  transmute(userId=V1)

set.seed(opt$seed)
for (j in 1:n_heritabilities) {
  heritability = heritabilities[j]
  # Subset to the `n_phenos` unused raw phenotypes for this round
  start = (j - 1) * n_phenos + 1
  end = j * n_phenos
  use_phenos = t(phenos[start:end,])

  # Create `n_phenos` noise phenotypes sampled from the standard normal distribution
  noise_pheno = matrix(rnorm(length(use_phenos)), nrow=nrow(use_phenos))

  # Create `n_phenos` phenotypes with the current heritability from the selected raw and noise phenotypes
  use_phenos = use_phenos*sqrt(heritability) + noise_pheno*sqrt(1 - heritability)

  # Pivot the phenotypes to a long table and compute the phenotypic correlation and print the info
  test = orig %>%
    rowwise() %>%
    mutate(pheno_corr=cor(use_phenos[x,], use_phenos[y,]))
  print(paste('Heritability:', heritability,
              '; Slope:', round(lm(pheno_corr ~ orig, test)$coefficients[[2]], 4),
              '; Correlation:', round(cor(test$pheno_corr, test$orig), 3)))

  # Create `n_prevalence` binary phenotypes with the current heritability
  if(opt$trait_type %in% c('binary', 'both') | !is.null(n_cases)){
      for (i in 1:n_prevalences) {
        start = (i - 1) * n_draws + 1
        end = i * n_draws
        out = data.frame(use_phenos[,start:end] < qnorms[i])
        output_df = output_df %>%
          bind_cols(out %>%
                      rename_with(function(x) {if(n_draws==1){paste0('p_', heritability, '_', binary_info[i], '_1')}else{paste0('p_', heritability, '_', binary_info[i], '_', gsub('X', '', x))}}) %>%
                      mutate_all(as.numeric))
      }
  }

  # Create `n_heritabilities` continuous phenotypes
  if(opt$trait_type %in% c('continuous', 'both')){
    continuous_df = data.frame(use_phenos[,(n_prevalences*n_draws + 1):n_phenos])
    output_df = output_df %>%
      bind_cols(continuous_df %>%
                  rename_with(function(x) {if(n_draws==1){paste0('p_', heritability, '_continuous_1')}else{paste0('p_', heritability, '_continuous_', gsub('X', '', x))}}))
  }

}

# Print a summary of the output random phenotype table
output_df %>% summary

# Write a tsv file of the random phenotypes
output_df %>%
  write_tsv(paste0(opt$out, '.tsv'))
