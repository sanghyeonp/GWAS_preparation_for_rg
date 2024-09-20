library(data.table)
library(dplyr)

metadata_prefix <- "metadata_for_external_trait_rg"

### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ###
### munged path는 metadata 파일에 추가 and 저장
df.metadata <- fread(paste0(metadata_prefix, ".qced_gwas_added.csv"), data.table=F)
unique_trait.list <- df.metadata$Trait_name_unique
df.munged <- data.frame()
for (unique_trait in unique_trait.list){
    f.munged <- normalizePath(paste0("./munged/", unique_trait, ".sumstats.gz"))
    if (!file.exists(f.munged)){f.munged <- NA}
    df.munged <- rbind(df.munged, data.frame(Trait_name_unique=unique_trait, Path_Munged=f.munged))
}

df.metadata <- df.metadata %>%
    left_join(df.munged, by="Trait_name_unique")
write.table(df.metadata, paste0(metadata_prefix, ".qced_gwas_munged_added.csv"),
            sep=",", row.names=F, quote=T)

### Munged 정보 따로 저장.
df.munged_parse <- data.frame()
for (unique_trait in unique_trait.list){
    f.munged <- normalizePath(paste0("./munged/", unique_trait, ".sumstats.gz"))
    nsnp <- NA; nsnp.nonmissing <- NA; mean_chisq <- NA; lambdaGC <- NA; max_chisq <- NA; nsnp_GWS <- NA
    if (file.exists(f.munged)){
        lines <- readLines(paste0("./munged/", unique_trait, ".log"))
        for (line in lines){
            if (grepl("Writing summary statistics for", line)){
                nsnp <- as.integer(strsplit(line, " ")[[1]][5])
                nsnp.nonmissing <- as.integer(gsub("\\(", "", strsplit(line, " ")[[1]][7]))
            }
            if (grepl("Mean chi", line)){mean_chisq <- as.numeric(strsplit(line, " ")[[1]][4])}
            if (grepl("Lambda GC", line)){lambdaGC <- as.numeric(strsplit(line, " ")[[1]][4])}
            if (grepl("Max chi", line)){max_chisq <- as.numeric(strsplit(line, " ")[[1]][4])}
            if (grepl("Genome-wide significant SNPs", line)){nsnp_GWS <- as.numeric(strsplit(line, " ")[[1]][1])}
        }
    }
    df.munged_parse <- rbind(df.munged_parse,
                        data.frame(Trait_name_unique=unique_trait,
                                Nsnp = nsnp, Nsnp_nonmissing=nsnp.nonmissing,
                                Mean_chisq=mean_chisq, lambdaGC=lambdaGC, Max_chisq=max_chisq,
                                Nsnp_GWS=nsnp_GWS, check.names=F))
}
write.table(df.munged_parse, paste0(metadata_prefix, ".munging_parsed.csv"),
            sep=",", row.names=F, quote=F)
