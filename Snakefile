## colocalization for COL27A1 ---------------------------
## @Xiaoping WU
## 2025-08-27


## define the main directory
#configfile: "workflow/config.yaml"
workdir: "/mnt/scratch/xiaoping/COL27A1/"


## define wildcards
import pandas as pd

# Read the manifest
manifest="resources/eQTL_Catalogue_manifest.txt"
eQTL_catalogue_file = pd.read_csv(manifest, sep='\t', header=0)

# Create unique study names
eQTL_catalogue_list = (eQTL_catalogue_file['study_label'] + '_' + eQTL_catalogue_file['quant_method'] + '_' + eQTL_catalogue_file['sample_group']).unique().tolist()



rule all:
	input:
		"results/colocalization/eQTL_catalogue_pph.tsv",
		"results/colocalization/eQTL_catalogue_snp_level.tsv"

## define target rules
rule coloc_eqtl:
	input:
		expand("results/colocalization/eQTL_catalogue/{eqtl_catalogue}_pph.tsv",eqtl_catalogue=eQTL_catalogue_list),
		expand("results/colocalization/eQTL_catalogue/{eqtl_catalogue}_snp_level.tsv.gz",eqtl_catalogue=eQTL_catalogue_list)


rule download:
	input:
		expand("results/eQTL_catalogue/temp/{eqtl_catalogue}_chr9.txt",eqtl_catalogue=eQTL_catalogue_list)


## define rules

## download eQTL data
rule download_and_query_tabix:
	input:
		manifest="resources/eQTL_Catalogue_manifest.txt"
	output:
		"results/eQTL_catalogue/temp/{eqtl_catalogue}_chr9.txt"
	params:
		region="9:"  # region to query
	log:
		"log/download_{eqtl_catalogue}.log"
	run:
		import pandas as pd
		import time
		from snakemake.shell import shell

		# Read manifest and select the file
		d = pd.read_csv(input.manifest, sep='\t', header=0)
		d['study'] = d.study_label + '_' + d.quant_method + '_' + d.sample_group
		d = d.loc[d.study == wildcards.eqtl_catalogue, :]
		if d.empty:
			raise ValueError(f"No matching entry for {wildcards.eqtl_catalogue}")
            
		# Convert ftp -> http
		url = 'https://' + d.ftp_path.values[0].replace('ftp://', '')
		time.sleep(1)  # optional delay
        
		# Local paths
		local_file = str(output[0]).replace("_chr9.txt", ".tsv.gz")
		index_file = local_file + ".tbi"

		# Download  desired region and save to output
		print('Downloading from the following url: ' + url)
		shell(f"tabix {url} {params.region} > {output[0]}")
		tbi= url.split('/')[-1] + '.tbi'
		shell(f"rm {tbi}")
        
 
    
    
## colocalization
rule GTEx_eQTL_colocalization:
	input:
		Rcode="workflow/scripts/coloc_eQTL.R",
		gwas="data/COL27A1_region.txt",
		eqtl="results/eQTL_catalogue/temp/{eqtl_catalogue}_chr9.txt"
	params:
		eqtl_data_from="{eqtl_catalogue}_chr9.txt"
	output:
		pph="results/colocalization/eQTL_catalogue/{eqtl_catalogue}_pph.tsv",
		snp_level="results/colocalization/eQTL_catalogue/{eqtl_catalogue}_snp_level.tsv.gz"
	log:
		log="log/coloc_{eqtl_catalogue}.log"
	shell:
		"""
		Rscript --vanilla {input.Rcode} \
            		--gwas {input.gwas} \
            	--eqtl {input.eqtl} \
            	--eqtl_data_from {params.eqtl_data_from} \
            	--pph {output.pph} \
            	--snp_level {output.snp_level} \
            	2>&1 | tee {log}
		"""
        


## summerize colocalization results
rule summary_eqtl_coloc:
	input:
		files=expand("results/colocalization/eQTL_catalogue/{eqtl_catalogue}_pph.tsv",eqtl_catalogue=eQTL_catalogue_list)
	output:
		"results/colocalization/eQTL_catalogue_pph.tsv"
	params:
        	"results/colocalization/eQTL_catalogue"
	log:
		"log/coloc_eqtl_merge_filter.log"
	shell:
		"""
		chmod +x ./workflow/scripts/coloc_eqtl_filter.sh
		./workflow/scripts/coloc_eqtl_filter.sh {params} {output}  &> {log}
		"""
        
#	script:
#		"workflow/scripts/coloc_eqtl_merge_and_filter.py"
#

## summerize SNP-level posterior probability 
rule summary_eqtl_coloc_snp_level:
        input:
                files=expand("results/colocalization/eQTL_catalogue/{eqtl_catalogue}_snp_level.tsv.gz",eqtl_catalogue=eQTL_catalogue_list)
        output:
                "results/colocalization/eQTL_catalogue_snp_level.tsv"
        params:
                "results/colocalization/eQTL_catalogue"
        log:
                "log/coloc_eqtl_merge_filter_snp_level.log"
        shell:
                """
                chmod +x ./workflow/scripts/coloc_eqtl_filter_snp_level.sh
                ./workflow/scripts/coloc_eqtl_filter_snp_level.sh {params} {output}  &> {log}
                """
