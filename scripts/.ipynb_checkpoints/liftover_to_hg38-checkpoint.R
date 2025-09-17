cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Initiating liftover...\n"))

### Packages and options ---------------
library(optparse)
library(rtracklayer)
library(dplyr)
library(stringr)
library(data.table)
library(readr)

parsed_opts <- list(
  make_option("--input",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--chain_path",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--output",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--from_build",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--to_build",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--chr_column",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--pos_column",
              action = "store",
              default = NA,
              type = 'character')
)

if (interactive()) {
  opt <- list(input = "data/COL27A1_region.txt",
              output = "results/COL27A1_region_hg38.txt",
              chain_path="resources/hg19ToHg38.over.chain.tab",
              chr_column = "CHR",
              pos_column = "GENPOS",
              from_build = "hg19",
              to_build = "hg38")
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}


### Import and Loading ---------------
cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Reading input file ", opt$input," ..."))

raw_input <- data.table::fread(opt$input, data.table = F)

raw_input <- raw_input %>% 
        mutate(CHR=paste("chr",CHROM,sep=""))

head(raw_input)
raw_input_columns <- colnames(raw_input)
raw_without_pos <- raw_input_columns[-which(raw_input_columns == opt$pos_col | raw_input_columns == opt$chr_col)]

unique_out_columns <- make.names(c("seqnames", "start", "end", "width", "strand", raw_without_pos),
                                 unique = T)

processed_column_order <- c(opt$chr_column, 
                            opt$pos_column, 
                            unique_out_columns[6:length(unique_out_columns)])

raw_input$input_row_index <- seq_len(nrow(raw_input))



cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Running liftover..."))

liftover_chain <- import.chain(opt$chain_path)

input_granges <- 
  GenomicRanges::makeGRangesFromDataFrame(raw_input,
                                          seqnames.field = "CHR",
                                          start.field = opt$pos_column,
                                          end.field = opt$pos_column,
                                          ignore.strand = TRUE,
                                          keep.extra.columns = TRUE)

gr_lifted <- liftOver(x = input_granges, chain = liftover_chain) %>% 
  unlist() %>% 
  data.frame()


cat("Complete!\n")

cat(paste0("[", format(Sys.time(), "%X"), "] ",
           "Writing output files..."))

colnames(gr_lifted)[1] <- opt$chr_column
colnames(gr_lifted)[2] <- opt$pos_column


gr_lifted %>% 
  select(all_of(processed_column_order)) %>% 
  readr::write_tsv(file = opt$output)

cat("Complete!\n")