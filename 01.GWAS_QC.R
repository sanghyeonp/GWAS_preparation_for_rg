####
# - ANNOVAR if needed
# - Filter HapMap3 SNPs
####
source("fnc.GWAS_QC.R")
library(data.table)
library(dplyr)

metadata_prefix <- "metadata_for_external_trait_rg"
nthread <- 20


df.metadata <- fread(paste0(metadata_prefix, ".csv"), data.table=F)

file.hm3 <- "/data1/sanghyeon/wonlab_contribute/combined/software/ldsc/data/w_hm3.snplist"
snplist.hm3 <- fread(file.hm3, data.table=F, select=c("SNP")) %>% pull(SNP)

dir.create("./qced_gwas", showWarnings=F)
dir.create("./log.GWAS_QC", showWarnings=F)

df.qced_path <- data.frame()
for (idx1 in 1:nrow(df.metadata)){
    prefix <- df.metadata[idx1,]$Trait_name_unique

    f.out <- paste0("/data1/sanghyeon/data/GWAS_for_rg/qced_gwas/", prefix, ".tsv")
    if(file.exists(f.out)) {df.qced_path <- rbind(df.qced_path, data.frame(Prefix=prefix, Path_QCed_GWAS=f.out)); next}

    file.gwas <- df.metadata[idx1,]$Path_raw_GWAS
    col_snp <- df.metadata[idx1,]$col_SNP; col_chr <- df.metadata[idx1,]$col_CHR; col_pos <- df.metadata[idx1,]$col_POS
    col_a1 <- df.metadata[idx1,]$col_EA; col_a2 <- df.metadata[idx1,]$col_OA; parse_OA <- grepl(";", col_a2)
    col_b <- df.metadata[idx1,]$col_BETA; is_or <- as.logical(df.metadata[idx1,]$is_OR); col_se <- df.metadata[idx1,]$col_SE
    col_p <- df.metadata[idx1,]$col_P; col_n <- df.metadata[idx1,]$col_N

    colnames.given <- c(col_snp, col_chr, col_pos, col_a1, col_b, col_se, col_p)
    if(!grepl(";", col_a2)){colnames.given <- c(colnames.given, col_a2)}
    parse_A2_from_alleles <- FALSE; if(grepl(";", col_a2)){parse_A2_from_alleles <- TRUE}
    if(!is.na(col_n) & !grepl(";", col_n)) {colnames.given <- c(colnames.given, col_n)}
    parse_Nobs_from_case_control <- FALSE; if(grepl(";", col_n)) {parse_Nobs_from_case_control <- TRUE}
    
    ### Drop any NA columns
    colnames.given <- colnames.given[which(!is.na(colnames.given))]

    Nobs <- df.metadata[idx1,]$Nobs_input; annovar <- as.logical(df.metadata[idx1,]$ANNOVAR); gbuild <- df.metadata[idx1,]$genome_build

    log.file <- file(paste0("log.GWAS_QC/GWAS_QC.", prefix, ".log"),open="wt"); .LOG(":: GWAS QC ::\n\t- Prefix: ", prefix, file=log.file)

    ## Download if GWAS is absent
    if (!file.exists(file.gwas)) {download_GWAS(cmd_wget=df.metadata[idx1,]$Download_wget); .LOG("\n[Downloading GWAS]", file=log.file)}

    ## Read GWAS
    df.gwas <- fread(file.gwas, data.table=F, fill=TRUE)
    colnames.gwas <- colnames(df.gwas)

    ## Sanity check
    .LOG("\n[Sanity check]", file=log.file); sanity <- TRUE
    if (all(colnames.given %in% colnames.gwas)!=TRUE) {.LOG("\tERROR: Specified column name(s) is/are missing. [", paste(setdiff(colnames.given, colnames.gwas), collapse=", "), "]", file=log.file); sanity <- FALSE}
    if (annovar & is.na(gbuild)) {.LOG("\tERROR: Need genome build to run ANNOVAR.", file=log.file); sanity <- FALSE}
    if (annovar & (is.na(col_chr) | is.na(col_pos))) {.LOG("\t:ERROR: Need chromosome and base position to run ANNOVAR.", file=log.file); sanity <- FALSE}
    if (!sanity) {next}
    .LOG("Passed all sanity checks!", file=log.file)

    ## QC in action
    .LOG("\n[QC arguments]", file=log.file)
    .LOG("\t- Specified columns: ", paste(colnames.given, collapse=", "), file=log.file)
    .LOG("\t- rsID mapping by ANNOVAR? ", annovar, file=log.file)
    .LOG("\t- Have N column? ", !is.na(col_n), file=log.file)
    .LOG("\t\t- Compute total sample size from Ncase and Ncontrol? ", parse_Nobs_from_case_control, file=log.file)
    col_ncase <- NA; col_ncontrol <- NA
    if(parse_Nobs_from_case_control){
        col_ncase <- strsplit(col_n, ";")[[1]][1]; col_ncontrol <- strsplit(col_n, ";")[[1]][2]
        if (!(all(c(col_ncase, col_ncontrol) %in% colnames.gwas))){.LOG("\tERROR: Columns specified in N is not present.", file=.log.file)}
        colnames.given <- c(colnames.given, col_ncase, col_ncontrol); col_n <- NA
    }
    .LOG("\t\t- Nobs is specified? ", !is.na(Nobs), " (", ifelse(!is.na(Nobs), format(Nobs, big.mark=","), "NA"), ")", file=log.file)
    .LOG("\t- Need OR -> BETA conversion? ", is_or, file=log.file)
    .LOG("\t- Need to parse OA from EA, Allele1, and Allele2? ", parse_OA, file=log.file)
    col_allele1 <- NA; col_allele2 <- NA
    if (parse_OA){
        col_allele1 <- strsplit(col_a2, ";")[[1]][1]; col_allele2 <- strsplit(col_a2, ";")[[1]][2]
        if (!(all(c(col_allele1, col_allele2) %in% colnames.gwas))) {.LOG("\tERROR: Columns specified in OA is not present.", file=log.file)}
        .LOG("\t\t- Allele columns being used: ", col_allele1, ", ", col_allele2, file=log.file)
        colnames.given <- c(colnames.given, col_allele1, col_allele2); col_a2 <- NA
    }
    
    lookup <- renamming_vector(col_snp, col_chr, col_pos, col_a1, col_a2, col_b, col_se, col_p, col_n, col_allele1, col_allele2, col_ncase, col_ncontrol)
    df.gwas.subset <- df.gwas %>%
        dplyr::select(all_of(lookup))

    if (parse_A2_from_alleles){df.gwas.subset$A2 <- apply(df.gwas.subset, 1, function(row){ifelse(row[["A1"]] == row[["Allele1"]], row[["Allele2"]], row[["Allele1.a"]])})}
    df.gwas.subset <- df.gwas.subset %>%
        mutate(A1=toupper(A1), A2=toupper(A2))

    if (annovar) {df.gwas.subset <- rsid_mapping(df.gwas.subset, gbuild, nthread)}

    df.gwas.subset <- df.gwas.subset %>% filter(SNP %in% snplist.hm3)

    if (parse_Nobs_from_case_control){df.gwas.subset$N <- apply(df.gwas.subset, 1, function(row){as.numeric(row[["Ncase"]]) + as.numeric(row[["Ncontrol"]])})}
    if (!is.na(Nobs)){df.gwas.subset$N <- Nobs}

    if (is_or){df.gwas.subset$BETA <- apply(df.gwas.subset, 1, function(row) log(as.numeric(row[["BETA"]])))}

    df.gwas.final <- df.gwas.subset %>%
        dplyr::select(SNP, A1, A2, BETA, SE, P, N)
    
    write.table(df.gwas.final, f.out, sep="\t", row.names=F, quote=F)
    df.qced_path <- rbind(df.qced_path, data.frame(Prefix=prefix, Path_QCed_GWAS=f.out))
    
    .LOG("QCed GWAS written: ", f.out, "\n\n", file=log.file)
    flush(log.file); close(log.file)
}

df.metadata <- df.metadata %>%
    left_join(df.qced_path, by=c("Trait_name_unique"="Prefix"))

write.table(df.metadata, 
            paste0(metadata_prefix, ".qced_gwas_added.csv"),
            sep=",", row.names=F, quote=T)
