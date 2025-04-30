# Load necessary libraries
library(ggplot2)
library(VariantAnnotation)
library(GenomicRanges)
library(reshape2)
library(ComplexHeatmap)

# Define file paths
control_vcf <- 'control_sample.vcf'
treated_vcf <- 'treated_sample.vcf'

# Read VCF files
control <- readVcf(control_vcf, 'hg19')
treated <- readVcf(treated_vcf, 'hg19')

# SNP Analysis Pipeline
control_snps <- rowRanges(control)[rowRanges(control)$FILTER == 'PASS' & isSNV(control)]
treated_snps <- rowRanges(treated)[rowRanges(treated)$FILTER == 'PASS' & isSNV(treated)]

unique_control_snps <- setdiff(control_snps, treated_snps)
unique_treated_snps <- setdiff(treated_snps, control_snps)


## SNP Comparison Barplot ##
snps_df <- data.frame(
  Sample = c(rep('Control', length(unique_control_snps)), 
             rep('Treated', length(unique_treated_snps))),
  SNP = c(names(unique_control_snps), names(unique_treated_snps))
)

ggplot(snps_df, aes(x = Sample, fill = Sample)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Unique SNP Distribution", 
       y = "Count", 
       subtitle = "Comparative Analysis of Control vs Treated Samples")

## Genomic Feature Heatmap ##
snp_matrix <- matrix(sample(0:1, 20, replace = TRUE), nrow = 4)
rownames(snp_matrix) <- paste0('SNP', 1:4)
colnames(snp_matrix) <- c('Control', 'Treated', 'Sample3', 'Sample4')

Heatmap(snp_matrix, 
        name = 'SNP Presence',
        column_title = "Sample Comparison",
        row_title = "Genetic Variants")

## Structural Variant Analysis ##
analyze_structural_variants <- function(vcf, genome = "hg19") {
  library(VariantAnnotation)
  library(GenomicRanges)
  library(ggplot2)
  library(karyoploteR)
  
  # Read VCF if input is file path
  if(is.character(vcf)) {
    vcf <- readVcf(vcf, genome)
  }
  
  # Extract structural variants
  sv_gr <- rowRanges(vcf)
  sv_types <- info(vcf)$SVTYPE
  sv_gr <- sv_gr[!is.na(sv_types) & sv_types %in% c("DEL", "DUP", "INV", "INS", "BND")]
  
  # Add metadata columns
  mcols(sv_gr)$SVTYPE <- sv_types[!is.na(sv_types) & sv_types %in% c("DEL", "DUP", "INV", "INS", "BND")]
  mcols(sv_gr)$SVLEN <- info(vcf)$SVLEN[match(names(sv_gr), rownames(vcf))]
  
  # Create summary statistics
  sv_summary <- data.frame(
    Type = names(table(mcols(sv_gr)$SVTYPE)),
    Count = as.numeric(table(mcols(sv_gr)$SVTYPE))
  )
  
  # Visualization: SV type distribution
  p1 <- ggplot(sv_summary, aes(x = Type, y = Count, fill = Type)) +
    geom_col() +
    theme_minimal() +
    labs(title = "Structural Variant Distribution", 
         subtitle = "By variant type")
  
  # Visualization: Genome-wide SV distribution
  p2 <- plotKaryotype(genome = genome) +
    kpPlotDensity(sv_gr[mcols(sv_gr)$SVTYPE == "DEL"], col = "red") +
    kpPlotDensity(sv_gr[mcols(sv_gr)$SVTYPE == "DUP"], col = "blue") +
    kpAddLegend(text = c("Deletions", "Duplications"), 
                col = c("red", "blue"))
  
  return(list(
    structural_variants = sv_gr,
    summary = sv_summary,
    plots = list(type_distribution = p1, genome_view = p2)
  ))
}


## Plasmid Sequence Identification ##
detect_plasmid_sequences <- function(assembly_file, 
                                     db_path = "plasmid_db/PlasmidFinder.fasta",
                                     min_identity = 95,
                                     min_coverage = 90) {
  # Load required packages
  if (!requireNamespace("blastr", quietly = TRUE)) {
    BiocManager::install("blastr")
  }
  library(blastr)
  library(ggplot2)
  library(knitr)
  
  # Validate inputs
  if (!file.exists(assembly_file)) {
    stop("Assembly file not found: ", assembly_file)
  }
  
  if (!file.exists(db_path)) {
    stop("Plasmid database not found. Download from:\n",
         "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/")
  }
  
  # Create BLAST database (if not existing)
  if (!file.exists(paste0(db_path, ".nin"))) {
    makeblastdb(db_path, db_type = "nucl")
  }
  
  # Run BLASTn search
  blast_results <- blastn(
    query = assembly_file,
    db = db_path,
    outfmt = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs",
    evalue = 1e-10
  ) %>% 
    filter(as.numeric(pident) >= min_identity,
           as.numeric(qcovs) >= min_coverage)
  
  # Annotate significant hits
  plasmid_hits <- blast_results %>%
    group_by(qseqid) %>%
    arrange(desc(bitscore)) %>%
    slice_head(n = 1) %>%
    mutate(
      Plasmid_Type = gsub("_.*", "", sseqid),
      Confidence = case_when(
        pident >= 98 & qcovs >= 98 ~ "High",
        pident >= 95 & qcovs >= 95 ~ "Medium",
        TRUE ~ "Low"
      )
    )
  
  # Visualization: Top plasmid hits
  p <- ggplot(plasmid_hits, aes(x = Plasmid_Type, y = pident, fill = Confidence)) +
    geom_col() +
    theme_minimal() +
    labs(title = "Plasmid Detection Summary",
         x = "Plasmid Type",
         y = "Percentage Identity (%)",
         caption = paste("Assembly:", assembly_file))
  
  # Generate report table
  report_table <- plasmid_hits %>%
    select(Contig = qseqid, Plasmid_Type, 
           Identity = pident, Coverage = qcovs, 
           Confidence, Evalue = evalue)
  
  return(list(
    plot = p,
    table = report_table,
    raw_data = blast_results
  ))
}



