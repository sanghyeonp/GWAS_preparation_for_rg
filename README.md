# Preparation of external traits for genetic correlation analysis using LDSC

```
Working directory: /data1/sanghyeon/data/GWAS_for_rg/
```

Acknowledgment: The origin of the metadata being utilized is from Marios Georgakis tweet (https://x.com/MariosGeorgakis/status/1825216813342605421).

## Step 1: Preparse metadata

Specify the metadata (csv) file name prefix to variable, `metadata_prefix` and number of threads to use at `nthread`.

### Metadata (csv) structure

Example metadata: `metadata_for_external_trait_rg.csv`

*Requires following information:*
- `Trait_name_unique`: Any UNIQUE name for a trait-source
- `Abbreviation` (optional): Abbreviation for the trait (does not have to be unique)
- `Source` (optional): \<Author>\<Year> format to identify source of the GWAS summary statistics
- `Category` (optional): Category that a trait belongs to
- `Trait_name_figure` (optional): Name of the trait that will be used in the figure
- `Trait_type`: Either "continuous" or "binary"
- `Ancestry` (optional): 
- `Path_raw_GWAS`: The absolute path to the GWAS summary statistics
- `col_SNP`: SNP column name (if not present, NA)
- `ANNOVAR`: Specify to perform ANNOVAR to map rsID (either F or T)
- `genome_build`: Genome-build for the variant genomic position (either GRCh36, GRCh37, or GRCh38)
- `col_CHR` (optional if not running ANNOVAR): Chromosome column name (if not present, NA)
- `col_POS` (optional if not running ANNOVAR): Base position column name (if not present, NA)
- `col_EA`: Effect allele column name
- `col_OA`: Other allele column name
    - Sometimes OA needs to be inferred from given EA column and Allele1/2 columns. Then, specify the name of the columns for the two alleles with `;` separator to parse OA automatically. (e.g., if `col_EA` is "Effect_allele" and other two allele columns are "ALT" and "REF", where OA needs to be inferred based on "Effect_allele", "ALT", and "REF", then specify "ALT;REF" for `col_OA`)
- `col_BETA`: Effect size column name
- `is_OR`: Specify whether `col_BETA` is odds ratio or not (either F or T)
- `col_SE`: Specify standard error column name
- `col_P`: Specify P-value column name
- `col_N` (optional if `Nobs_input` is specified): Specify N column name if present
- `Nobs_input` (optional if `col_N` is specified): Specify N if `col_N` is not present (For binary traits, specify sum of Ncase and Ncontrol)
- `Source_Nobs` (optional): Note on where Nobs was found
- `Download_wget`: wget command to download the GWAS summary statistics
- `Reference`: Reference for the journal article
- `Link_to_data`: Link to the GWAS summary statistics



## Step 2: Perform GWAS QC

Use `01.GWAS_QC.R` to perform QC on the raw GWAS summary statistics.

The QC process goes as follows:
1. Sanity check for the specified column names specified, sample size, parsing OA and so on. 
2. Prase OA if needed.
3. Run ANNOVAR to map rsID if needed.
4. Subset only HapMap3 variants that are being used in LDSC munge.
5. Handle N column.
6. Conversion of OR to BETA if needed.
7. Save the QCed-GWAS summary statistics.
8. Generate new metadata with additional column, `Path_QCed_GWAS`, with the following suffix `qced_gwas_added` added to the metadata file name prefix.

The output QCed GWAS will be saved in `./qced_gwas/` with the following file name format, `<unique_name_given>.tsv`.

The log file for each GWAS being QCed will be saved in `./log.GWAS_QC/` with the following file name format, `GWAS_QC.<unique_name_given>.log`.
The content of the log will contain the information like the example below:
```
:: GWAS QC ::
	- Prefix: AAA.Klarin2020

[Sanity check]
Passed all sanity checks!

[QC arguments]
	- Specified columns: SNP_ID, CHROM, POS, EA, OR, SE, P, N
	- rsID mapping by ANNOVAR? FALSE
	- Have N column? TRUE
		- Compute total sample size from Ncase and Ncontrol? FALSE
		- Nobs is specified? FALSE (NA)
	- Need OR -> BETA conversion? TRUE
	- Need to parse OA from EA, Allele1, and Allele2? TRUE
		- Allele columns being used: Allele1, Allele2
QCed GWAS written: /data1/sanghyeon/data/GWAS_for_rg/qced_gwas/AAA.Klarin2020.tsv
```

**Please search the content of the log with a term "ERROR:" to check if there were any error while performing the QC. QC process will not stop and continue to the next GWAS in queue although an error has occurred.**


## Step 3. Munging

Use `02_1.munge.sh` to run munging using LDSC and within the script, specify variable `metadata_prefix`.
Then, use `02_2.parse_munge.R` to make metadata with `Path_Munged` column added and saved with the following suffix `qced_gwas_munged_added` added to the metadata file name prefix. Also, a separate file with parsed information from munged logs which will be saved with the following suffix `munging_parsed`.

The `Path_Munged column` in the metadata output `<metadata_prefix>.qced_gwas_munged_added.csv` can be used for further genetic correlation analysis with LDSC.
