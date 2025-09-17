library(optparse)
library(data.table)
library(dplyr)
library(coloc)



parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--eqtl",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--eqtl_data_from",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--pph",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--snp_level",
              action = "store",
              default = NA,
              type = 'character')
)

if (interactive()) {
  opt <- list(gwas = "data/COL27A1_region.txt",
              eqtl="results/eQTL_catalogue/temp/GTEx_tx_nerve_tibial_chr9.txt",
              eqtl_data_from="GTEx_tx_skin_not_sun_exposed_chr9.txt",
              pph="results/eQTL_catalogue/",
              snp_level="results/eQTL_catalogue/"
             )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}





# -----------------------------
# Step 1: Load your datasets
# -----------------------------
## gwas summary statistics
gwas <- fread(opt$gwas,header=TRUE) 
print(nrow(gwas))
head(gwas)
#beta is effect for A1 allele



# eQTL summary statistics
if (file.info(opt$eqtl)$size > 0) {

eqtl <- fread(opt$eqtl,header = FALSE,col.names = c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid"))
head(eqtl)
## beta is effect of alernative allele
                

## eqtl data origin
eqtl_data_from=sub("_chr9\\.txt$", "", opt$eqtl_data_from)

# -----------------------------
# Step 2: flip SNP effect to match and handel duplocated SNP in eQTL data
# -----------------------------
# filter to the same SNPs
common_snps <- intersect(gwas$ID, eqtl$rsid)
gwas <- gwas[ID %in% common_snps] %>%
         mutate(sdY_est = sqrt(SE^2 * 2 * N * MAF * (1 - MAF)))

eqtl <- eqtl[rsid %in% common_snps]


## match effect allele
eqtl2 <- merge(eqtl,gwas[,c("ID","A1","A2")],by.x="rsid",by.y="ID")
eqtl2 <- eqtl2 %>%
        mutate(BETA = case_when(alt==A1 & ref==A2 ~ beta,
		alt==A2 & ref==A1 ~ -beta,
		TRUE ~ NA)
        ) %>%
	filter(!is.na(BETA)) %>%
        mutate(sdY_est = sqrt(se^2 * 2 * (an/2) * maf * (1 - maf)))
#so BETA is for A1

## Keep one record per SNP-gene, using smallest p-value
#in eQTL data, the same SNP can appear for different genes. This is normal because:
#	1.	A SNP can be a cis-eQTL for multiple nearby genes.
#	2.	Some SNPs may influence transcripts of different genes.
#	3.	eQTL datasets often include all significant SNP-gene pairs, so duplication by SNP ID across genes is expected.

eqtl_unique <- eqtl2 %>%
  group_by(rsid, gene_id) %>%
  slice_min(order_by = pvalue, n = 1) %>%  # keep row with smallest p-value
  ungroup()




# -----------------------------
# Step 3: Format gwas datasets for coloc
# -----------------------------
#
dataset1 <- list(
  beta    = gwas$BETA,
  varbeta = gwas$SE^2,
  type    = "quant",          # quantitative trait; use "cc" for case-control
  snp     = gwas$ID,
  MAF     = gwas$MAF,
  N       = gwas$N,
  sdY = gwas$sdY_est  
)


# -----------------------------
# Step 4: Format qQTL datasets for coloc, and run for each gene
# -----------------------------


#For coloc.abf():
#	•	You run colocalization per gene.
## Get list of genes and Initialize result lists
genes <- unique(eqtl_unique$gene_id)
gene_results_list <- list()
snp_results_list <- list()

## Run coloc per gene
for(g in genes) {

    eqtl_gene <- eqtl_unique %>% 
        filter(gene_id == g)
    cat("number of lines for gene ",g, ": ",nrow(eqtl_gene),"\n")
    
    
    ## remove duplicated SNPs
    eqtl_gene <- eqtl_gene %>%
        group_by(rsid) %>%
        slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
        ungroup()
   cat("number of lines for gene ",g, " after remove duplicated SNPs: ",nrow(eqtl_gene),"\n")
   
  
    #for eQTL data for coloc
    dataset2 <- list(
        beta    = eqtl_gene$BETA,
        varbeta = eqtl_gene$se^2,
        type    = "quant",
        snp     = eqtl_gene$rsid,
        MAF     = eqtl_gene$maf,
        N       = eqtl_gene$an /2,
        sdY  = eqtl_gene$sdY_est
    )

  
  # Run colocalization
  res <- coloc.abf(dataset1, dataset2)
  
  # Return summary
    # Example output columns:
    # PP.H0  PP.H1  PP.H2  PP.H3  PP.H4
    # posterior probability for hypotheses:
    # H0: neither trait associated
    # H1: trait1 only
    # H2: trait2 only
    # H3: both associated, different causal variants
    # H4: both associated, shared causal variant (colocalized)

    gene_results_list[[g]] <- data.frame(
    eqtl_data_from = eqtl_data_from,
    gene_id = g,
    PP.H0 = res$summary["PP.H0.abf"],
    PP.H1 = res$summary["PP.H1.abf"],
    PP.H2 = res$summary["PP.H2.abf"],
    PP.H3 = res$summary["PP.H3.abf"],
    PP.H4 = res$summary["PP.H4.abf"]
  )
    snp_level <- res$results
    snp_level$gene_id <- g
    snp_level$eqtl_data_from <- eqtl_data_from
    snp_results_list[[g]] <- snp_level
}

## Combine results into a data frame
gene_summary <- do.call(rbind, gene_results_list)
snp_summary  <- do.call(rbind, snp_results_list)


} else {
  message("File is empty: ", opt$eqtl)
  eqtl <- NULL
	gene_summary <- data.frame(eqtl_data_from=character(),gene_id=character(),PP.H0 = numeric(),  PP.H1 = numeric(),  PP.H2 = numeric(),PP.H3 = numeric(),  PP.H4 = numeric(),  stringsAsFactors = FALSE)
	snp_summary  <- data.frame(snp=character(),V.df1=numeric(),z.df1 = numeric(),  r.df1 = numeric(),  lABF.df1= numeric(),V.df2 = numeric(), z.df2 = numeric(),r.df2=numeric(),lABF.df2= numeric(),internal.sum.lABF=numeric(),SNP.PP.H4=numeric(),gene_id=character(),eqtl_data_from=character() ,stringsAsFactors = FALSE)
}

write.table(gene_summary,opt$pph,row.names=F,quote=F,sep = "\t",col.names=T )
write.table(snp_summary,gzfile(opt$snp_level),row.names=F,quote=F,sep = "\t", col.names=T)
