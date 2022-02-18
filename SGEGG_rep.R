# Packages
suppressMessages(require(data.table))
suppressMessages(require(parallel))
suppressMessages(library(optparse))
suppressMessages(library(BEDMatrix))
suppressMessages(library(lme4))


options(stringsAsFactors = F)
option_list = list(
  make_option("--pheno", action = "store", default = NA, type = "character", help = "The path of the phenotype file"),
  make_option("--geno", action = "store", default = NA, type = "character", help = "The path of the genotype file"),
  make_option("--interaction", action = "store", default = NA, type = "character", help = "The path to the variable of interaction file"),
  make_option("--covar", action = "store", default = NA, type = "character", help = "The path of the covariate file"),
  make_option("--output", action = "store", default = NA, type = "character", help = "The path of the output file"),
  make_option("--num_cores", action = "store", default = NA, type = "integer", help = "Numbers of cores for parellel computing"),
  make_option("--type", action = "store", default = NA, type = "character", help = "GxE or GxG analysis (This will only change the columns names for output)"),
  make_option("--start", action = "store", default = NA, type = "integer", help = "Index of SNP that starts computing."),
  make_option("--end", action = "store", default = NA, type = "integer", help = "Index of SNP that ends computing.")
)

opt = parse_args(OptionParser(option_list=option_list))

# --- 0. I/O check
if (is.na(opt$pheno) | is.na(opt$geno) |  is.na(opt$interaction) | is.na(opt$covar) |is.na(opt$output) | is.na(opt$num_cores)) {
    cat("ERROR: Missing essential inputs.\n")
    q("no")
}


# Input parameters
phenotype <- opt$pheno
genotype <- opt$geno
covariate <- opt$covar
interaction_file <- opt$interaction
output <- opt$output
num_cores <- opt$num_cores
type <- opt$type
snp_start <- opt$start
snp_end <- opt$end

# phenotype <- "/z/Comp/lu_group/Members/jmiao24/Collaboration/Lauren/Lauren_LHS/10_SGEGG_Rep/data/pheno/bmi_rep_pheno.txt"
# genotype <- "/z/Comp/lu_group/Members/jmiao24/Collaboration/Lauren/Lauren_LHS/8_run_GxE/data/geno/lunghealth_imputed_hm3_qced"
# covariate <- "/z/Comp/lu_group/Members/jmiao24/Collaboration/Lauren/Lauren_LHS/8_run_GxE/data/pheno/covar.txt"
# interaction_file <- "/z/Comp/lu_group/Members/jmiao24/Collaboration/Lauren/Lauren_LHS/8_run_GxE/data/pheno/E.txt"
# output <- "/z/Comp/lu_group/Members/jmiao24/Collaboration/Lauren/Lauren_LHS/8_run_GxE/results/sumstats_tmp/bmi/sumstats_1.txt"
# num_cores <- 1
# type <- "GxE"
# start <- 1
# end <- 10

geno_bed <- suppressMessages(BEDMatrix((paste0(genotype, ".bed"))))
geno_bim <- fread(paste0(genotype, ".bim"))

if(is.na(opt$start) | is.na(opt$end)){
    snp_start <- 1
    snp_end <- nrow(geno_bim)
    cat("No index of SNP that starts or ends computing.\nSGEGG will fit Genome-wide interaction analysis for all SNPs.\n")
}


cat("Preparing the data:\n")
pheno <- as.data.frame(fread(phenotype, stringsAsFactors = F))
covar <- as.data.frame(fread(covariate, data.table = F, stringsAsFactors = F))
interaction <- as.data.frame(fread(interaction_file, data.table = F, stringsAsFactors = F))

# align the  phenotype, covaraite, genotype files
IID <- gsub(pattern = paste0(".*_(.*)"), rownames(geno_bed), replacement = "\\1")
IID_overlap <- intersect(pheno$IID, IID)
IID_overlap <- intersect(covar$IID, IID_overlap)
IID_overlap <- intersect(interaction$IID, IID_overlap)
pheno_lm <- pheno[pheno$IID %in% IID_overlap, ]

index_covar <- match(pheno_lm$IID, covar$IID)
index_interaction <- match(pheno_lm$IID, interaction$IID)
index_geno <- match(pheno_lm$IID, IID)
covar_lm <- as.data.frame(covar[index_covar, 3:ncol(covar)])
interaction_lm <- as.data.frame(interaction[index_interaction, 3:3])
colnames(interaction_lm) <- "I"

pheno_names <- colnames(pheno_lm)[3]
covar_names <- colnames(covar_lm)
interaction_names <- colnames(interaction_lm)
df_cb <- cbind(pheno_lm, covar_lm, interaction_lm)

# Progress bar of the mclapply; Source: https://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply/26892969#26892969
mclapply2 <- function(X, FUN, ..., 
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
    mc.cleanup = TRUE, mc.allow.recursive = TRUE,
    mc.progress=TRUE, mc.style=3) {
    if (!is.vector(X) || is.object(X)) X <- as.list(X)

    if (mc.progress) {
        f <- fifo(tempfile(), open="w+b", blocking=T)
        p <- parallel:::mcfork()
        pb <- txtProgressBar(0, length(X), style=mc.style)
        setTxtProgressBar(pb, 0) 
        progress <- 0
        if (inherits(p, "masterProcess")) {
            while (progress < length(X)) {
                readBin(f, "double")
                progress <- progress + 1
                setTxtProgressBar(pb, progress) 
            }
            cat("\n")
            parallel:::mcexit()
        }
    }
    tryCatch({
        result <- mclapply(X, function(...) {
                res <- FUN(...)
                if (mc.progress) writeBin(1, f)
                res
            }, 
            mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
            mc.silent = mc.silent, mc.cores = mc.cores,
            mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
        )

    }, finally = {
        if (mc.progress) close(f)
    })
    result
}

# SGEGG function 
SGEGG <- function(i){
    # Obtain the SNP information
    SNP <- geno_bed[index_geno, i]
    snp_name <- geno_bim$V2[i]
    chr <- geno_bim$V1[i]
    bp <- geno_bim$V4[i]
    a1 <- geno_bim$V5[i]
    a2 <- geno_bim$V6[i]

    ## Keep the non-NA in SNP
    df_cb$SNP <- SNP
    index_non_NA <- which(!is.na(SNP))
    df_cb <- df_cb[index_non_NA, ]
    
    # Calculate MAF
    maf <- sqrt(sum(df_cb$SNP == 2, na.rm=T)/sum(!is.na(df_cb$SNP)))

    # Add E*covar interaction
    predictors_I <- paste0(c("SNP", covar_names), "*I")
    f1 <- as.formula(paste(paste(pheno_names, paste("(1|IID)+" , paste(predictors_I, collapse = " + ")), sep = " ~ ")))

    fit <- lmer(f1, data = df_cb)
    # coeff <- summary(lm(f1, data = df_cb))$coefficients 
    summ <-  summary(fit)
    coeff <- summ$coefficients 
    ind1 <- which(rownames(coeff) == "SNP")
    ind2 <- which(rownames(coeff) == "I")
    ind3 <- which(rownames(coeff) == "SNP:I")
    SGEGG_results <- data.frame(chr, snp_name, bp, a1, a2, maf, 
    beta_M = coeff[ind1, 1], SE_M = coeff[ind1, 2], P_M = 2*pnorm(abs(coeff[ind1, 1]/coeff[ind1, 2]), 0, 1, lower.tail = F),
    beta_E = coeff[ind2, 1],SE_E = coeff[ind2, 2], P_E = 2*pnorm(abs(coeff[ind2, 1]/coeff[ind2, 2]), 0, 1, lower.tail = F),
    beta_I = coeff[ind3, 1], SE_I = coeff[ind3, 2], P_I = 2*pnorm(abs(coeff[ind3, 1]/coeff[ind3, 2]), 0, 1, lower.tail = F),
    N_rep = length(summ$residuals), N = summ$ngrps)
    return(SGEGG_results)
}


# Parallel of SGEGG
Fit_SGEGG <- function(start = snp_start, end = snp_end){
    df_out <- mclapply2(start:end, SGEGG, mc.cores = num_cores)
    df_out <- do.call(rbind, df_out)
    df_out <- as.data.frame(df_out)
    if (type == "GxE"){
        colnames(df_out) <-  c('CHR', 'SNP', 'BP', 'A1', 'A2', 'FREQ', 
        'BETA_M','SE_M','P_M', 
        'BETA_E','SE_E','P_E', 
        'BETA_I','SE_I','P_I', 'N')
    }else if (type == "GxG") {
       colnames(df_out) <-  c('CHR', 'SNP', 'BP', 'A1', 'A2', 'FREQ', 
        'BETA_G1','SE_G1','P_G1', 
        'BETA_G2','SE_G2','P_G2', 
        'BETA_I','SE_I','P_I', 'N')
    }

    return(df_out)
}

cat("Begin running Genome-wide interaction analysis:\n")

# To suppress some Warnings that are harmless: https://github.com/HenrikBengtsson/future/issues/218
df_all_snp <- suppressWarnings(Fit_SGEGG(start=snp_start, end =snp_end))

# Write out the results
cat("Write out the summary statistics.\n")
fwrite(df_all_snp, output, col.names = T, row.names = F, quote = F, sep = "\t", na = "NA")
cat("Finish!\n")
