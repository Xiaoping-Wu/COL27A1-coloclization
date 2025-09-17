rule tabix_eQTL_Catalog:
	'Download summary stats from eQTL Catalog using tabix, for all variants in chromosome 9.'
    'Select 150kb region abound SNP rs7023208, which has high LD with SNP rs72760655'
	priority: 1
	input:
		'resources/eQTL_Catalogue_manifest.txt'
	output:
		temp('results/eQTL_catalogue/temp/{eqtl_catalogue}.txt')
	threads: workflow.cores * 0.25
	run:
		d= pd.read_csv(input[0], sep= '\t', header= 0)
		d['study']= d.study_label + '_' + d.quant_method + '_' + d.sample_group
		d= d.loc[d.study == wildcards.eqtl_catalogue, :]
		url= 'http://' + d.ftp_path.values[0].replace('ftp://', '')
		time.sleep(2)
		print('Downloading from the following url: ' + url)
		shell("tabix {url} 9: > {output[0]}")
		tbi= url.split('/')[-1] + '.tbi'
		shell("rm {tbi}")


rule format_eQTL_catalog:
        'Format eQTL Catalogue data.'
	priority: 10
	input:
                'results/eQTL_catalogue/temp/{eqtl_catalogue}.txt'
        output:
                'results/eQTL_catalogue/formatted/{eqtl_catalogue}.txt.gz'
        run:
                cols= ['molecular_trait_id', 'CHR', 'BP', 'ref', 'alt', 'variant', 'ma_samples', 'maf', 'pvalue', 'beta', 'se', 'type', 'ac', 'an', 'r2', 'molecular_trait_object_id', 'gene_id', 'median_tpm', 'rsid']
                d= pd.read_csv(input[0], sep= '\t', header= None, names= cols)
                d['ref']= np.where(d.ref.str.len() > d.alt.str.len(), 'I', d.ref)
                d['alt']= np.where(d.ref.str.len() < d.alt.str.len(), 'I', d.alt)
                d['ref']= np.where(d.alt== 'I', 'D', d.ref)
                d['alt']= np.where(d.ref== 'I', 'D', d.alt)
                d['beta']= np.where(d.ref> d.alt, -1 * d.beta, d.beta)
                d.drop_duplicates(subset= ['rsid', 'gene_id'], inplace= True, keep= 'first')
                d['N']= d.an / 2
                d= d[['gene_id', 'rsid', 'BP', 'maf', 'N', 'ref', 'alt', 'pvalue', 'beta', 'se']]
                d.to_csv(output[0], sep= '\t', header= True, index= False, compression= 'gzip')


rule liftOver:
	'LiftOver from hg19 to hg38 using UCSC liftOver.'
	rule Osmolality_lifted_hg38:
	input:
		Rcode="workflow/scripts/liftover_to_hg38.R",
		chain_path="resources/hg19ToHg38.over.chain.tab",
		gwas="data/COL27A1_region.txt" 
	output:
		"results/d20230523/Bogers/Osmolality_lifted_to_hg38.tsv.gz"
	log:
		log="log/lifted_to_hg38.log"
	threads: 1
	resources:
		mem_gb = 10,
		runtime = "01:00:00"
	params:
		col_chr = "CHR",
		col_pos = "GENPOS",
		from_build = "hg19",
		to_build = "hg38"
	shell:
		"Rscript --vanilla {input.Rcode} --chain_path {input.chain_path} "
		"--input {input.gwas} --output {output} "
		"--chr_column {params.col_chr} --pos_column {params.col_pos} "
		"--from_build {params.from_build} --to_build {params.to_build} 2>&1 | tee {log.log}"


rule coloc_eQTL_catalog:
        'Colocalization COL27A1_region with eQTL data form eQTL Catalogue.'
        input:
                gwas="results/d20230523/Bogers/Osmolality_lifted_to_hg38.tsv.gz",
                eqtl="results/eQTL_catalogue/formatted/{eqtl_catalogue}.txt.gz"
        output:
                temp('results/eQTL_catalogue/temp/pph-{eqtl_catalogue}.txt'),
                temp('results/eQTL_catalogue/temp/SNP-{eqtl_catalogue}.txt'),
        conda:
                'workflow/envs/coloc.yml'
	threads: 4
        script:
                'workflow/scripts/coloc_eqtl.R'

rule concat_eQTL_catalog_coloc_pph:
	'concat results from Colocalization analysis with eQTL catalogue.'
	input:
		expand('results/eQTL_catalogue/temp/pph-{eqtl_catalogue}.txt', eqtl_catalogue= eQTL_catalogue_conditions)
	output:
		'results/eQTL_catalogue/pph_eqtl.txt'
	shell:
		head -1 {input[0]} > {output[0]}
		tail -n +2 -q {input} >> {output[0]}


rule concat_eQTL_catalog_coloc_results:
        'concat results from Colocalization analysis with eQTL catalogue.'
        input:
                expand('results/eQTL_catalogue/temp/SNP-{eqtl_catalogue}.txt', eqtl_catalogue= eQTL_catalogue_conditions)
        output:
                'results/eQTL_catalogue/SNP_eqtl.txt'
        shell:
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}