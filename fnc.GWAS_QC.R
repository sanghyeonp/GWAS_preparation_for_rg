

.LOG <- function(..., file, print = TRUE) {
    msg <- paste0(..., "\n")
    if (print) cat(msg)
    cat(msg, file = file, append = TRUE)
}



download_GWAS <- function(cmd_wget){
    cmd_wget <- paste0("cd /data1/sanghyeon/data/GWAS_for_rg/raw_gwas; ", cmd_wget)
    system(cmd_wget, wait=T, ignore.stdout=T)
}


renamming_vector <- function(col_snp, col_chr, col_pos, col_a1, col_a2, col_b, col_se, col_p, col_n, col_allele1, col_allele2, col_ncase, col_ncontrol){
    initial_vec <- c(col_snp, col_chr, col_pos, col_a1, col_a2, col_b, col_se, col_p, col_n, col_allele1, col_allele2, col_ncase, col_ncontrol)
    idx <- which(!is.na(initial_vec))
    rename_vec <- c("SNP", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "N", "Allele1", "Allele2", "Ncase", "Ncontrol")
    vec4rename <- initial_vec[idx]
    names(vec4rename) <- rename_vec[idx]
    return(vec4rename)
}

rsid_mapping <- function(df, gbuild, nthread=1){
    write.table(df %>% dplyr::select(CHR, POS, A1, A2), 
            "tmp.annovar_in.txt", sep="\t", row.names=F, quote=F)

    src <- "/data1/sanghyeon/wonlab_contribute/combined/src/annovar/annovar_rsid_map.R"
    cmd_annovar <- paste0("/data1/software/anaconda3/envs/R_SP_4.4/bin/Rscript ", src,
                        " --gwas tmp.annovar_in.txt --delim-in tab",
                        " --chr-col CHR --pos-col POS --ref-col A1 --alt-col A2",
                        " --genome-build ", ifelse(gbuild == "GRCh37", 19, 38),
                        " --dbgap-build 150", " --nthread ", nthread, " --outf tmp.annovar_out")
    system(cmd_annovar, wait=T, ignore.stdout=T)

    df.map <- fread("tmp.annovar_out.rsid_annotated", data.table=F, nThread=nthread) %>%
        dplyr::select(-flipped) %>% rename(SNP=rsid)
    if("SNP" %in% colnames(df)){df <- dplyr::select(df, -SNP)}
    df <- df %>%
        left_join(df.map, by=c("CHR", "POS", "A1", "A2")) %>%
        # Since LDSC uses rsID only. Drop all.
        filter(!is.na(SNP))

    system("rm tmp.annovar_in.txt tmp.annovar_out.annovin.hg19_avsnp150_dropped.merged tmp.annovar_out.rsid_annotated",
        wait=T, ignore.stdout=T)
    return(df)
}