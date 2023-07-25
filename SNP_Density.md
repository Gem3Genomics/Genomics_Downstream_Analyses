**SNP density** consists of the number of SNPs over a certain length of the genome, often ranging from 100 kilbases to 1 megabase. SNP density can be used as a way of visualizing genetic diversity across the genome via heatmaps. 

We will calculate SNP density with VCFtools, which is a package designed to work with VCF files and carry out the following: filter out specific variants, compare files, summarize variants, convert to different file types, validate and merge files, and create intersections and subsets of variants. To learn more about VCFtools or to download, visit [VCFtools](https://vcftools.github.io/index.html).

To apply VCFtools to your data to create SNP density heatmaps, do so in the following steps:

1) Analyze SNP Density with VCFtools

  1a) Create new VCF files isolating only heterozygous sites

  1b) Run VCFtools on data

  1c) Add sample names to each table

  1d) Concatenate

2) Plot in R


## 1) Analyze SNP Density with VCFtools

Steps 1-4 can be applied in one job script. Be sure modules are available and loaded on your HPC, and submit a job with the following code:

This function takes a VCF file as input and outputs a table with the SNP count for each region interval.

```Bash
#Create new VCF files isolating only heterozygous sites
vcftools --gzvcf individual_1.vcf.gz --recode --out individual_1_hetsites --maf 0.1

vcftools --gzvcf individual_2.vcf.gz --recode --out individual_2_hetsites --maf 0.1

vcftools --gzvcf individual_3.vcf.gz --recode --out individual_3_hetsites --maf 0.1

# Explanation
#gzvcf input.vcf.gz: specifies input as a VCF file in compressed gzip format
#recode: indicates that you want to create a new VCF file with the filtered results
#maf 0.1: Sets a minimum allele frequency (MAF) filter of 0.1, meaning that only variants with a MAF of at least 10% will be included in the output file

#Run VCFtools on data
vcftools --vcf individual_1_hetsites.recode.vcf --SNPdensity 1000000 --out individual_1_hetsites

vcftools --vcf individual_2_hetsites.recode.vcf --SNPdensity 1000000 --out individual_2_hetsites

vcftools --vcf individual_3.recode.vcf --SNPdensity 1000000 --out individual_3_hetsites

#Explanation
SNPdensity 1000000: calculates the SNP density in non-overlapping windows of 1,000,000 base pairs (1 Mb) across the genome. Note that SNP density is the number of SNPs per window.

#Add sample names to each table
awk -v sample="NN114296" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' Individual_1_hetsites.snpden > Individual_1_hetsites_id.snpden

awk -v sample="NN114296" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' Individual_2_hetsites.snpden > Individual_2_hetsites_id.snpden

awk -v sample="NN114296" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' Individual_3_hetsites.snpden > Individual_3_hetsites_id.snpden

#Explanation
'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}': an awk script that processes the input file line by line, checking if the current line number (NR) is equal to 1, and if it is, it prints the entire line ($0) followed by a tab and the string "Indiv". For all other lines (NR > 1), it prints the entire line followed by a tab and the value of the sample variable.
> Individual_1_hetsites_id.snpden: Redirects the output to a new file named "Individual_1_hetsites_id.snpden" in the current directory.

#Concatenate files
tail -q -n +2 *_id.snpden > _hetsites.snpden
```

## 2) Plot in R

```R
# Load the required packages
library(tidyverse)
library(gdata)

# Read the SNP density data file
snpden <- read.table("_hetsites.snpden.txt", header = T)

# Define the order of the scaffolds to be used in the visualization
target <- c("CM051599.1", "CM051600.1", "CM051601.1", "CM051602.1", "CM051603.1", "CM051604.1", 
            "CM051605.1", "CM051606.1", "CM051607.1", "CM051608.1", "CM051609.1", "CM051610.1", 
            "CM051611.1", "CM051612.1", "CM051613.1", "CM051614.1", "CM051615.1", "CM051616.1", 
            "CM051617.1")

#Define the order of the chromosomes to be used in the visualization
chr <-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chrX')

snpden.master <- snpden

# Reorder the chromosome column of the data frame according to the target order
snpden.master$CHROM <- reorder.factor(snpden.master$CHROM, new.order = target)

# Subset data from chromosomes that are not "NA"
snpden.master <-subset(snpden.master, snpden.master$CHROM!='NA')

snpden.master$groups <- cut(as.numeric(snpden.master$VARIANTS.KB), 
                            c(0,0.05,0.1,0.15,0.20,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,
                              3,3.25,3.5,3.75,4,4.25,4.5,4.75,5),
                            include.lowest = TRUE, labels = c("0", "0.05-0.1", "0.1-0.15", "0.15-0.2", 
                                                              "0.2-0.25", "0.25-0.5", "0.5-0.75", "0.75-1", 
                                                              "1-1.25", "1.25-1.5", "1.5-1.75", "1.75-2", 
                                                              "2-2.25", "2.25-2.5", "2.5-2.75", "2.75-3", 
                                                              "3-3.25", "3.25-3.5", "3.5-3.75", "3.75-4", 
                                                              "4-4.25", "4.25-4.5", "4.5-4.75", "4.75-5"))
snpden.master$groups[snpden.master$VARIANTS.KB == 0] <- "0"

# Rename CHROM levels
levels(snpden.master$CHROM) <-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chrX')

snpden.master$BIN_START <- as.numeric(as.character(snpden.master$BIN_START))

names_vec <- c("Individual_1", "Individual_2", "Individual_3")

for (individual in unique(snpden.master$Indiv)) {

# Subset the data for the current chromosome
snpden.chr <- subset(snpden.master, snpden.master$Indiv == individual)
  
# Define title
title<-expression(paste(italic("Title")))
  
# Define title
title<-expression(paste(italic("Title")))
  
# Create ggplot object 
  snpden_plot <- snpden.chr %>% 
    mutate(Indiv = factor(Indiv, levels = c("Individual_1", "Individual_2", "Individual_3"))) %>%
    ggplot(aes(x=BIN_START, y=1)) + 
    geom_tile(aes(fill=groups)) +
    facet_grid(CHROM ~ ., switch='y') +
    labs(x = 'Chromosome Length' , 
         y = 'Scaffold Number' , 
         title = expression(paste(italic("Title"))), 
         subtitle = paste0(individual, " heterozygous SNP densities")) + 
    theme_minimal() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.text.y.left = element_text(angle = 0, size=8),
          panel.spacing.y=unit(0.15, "lines"),
          plot.title = element_text(hjust = .5, size = 15),
          plot.subtitle = element_text(hjust = .5, size = 13, color = "gray")) +
    scale_fill_manual(values = c("#000081", "#0000f3", "#004dff", "#00b3ff", "#29ffce", 
                                          "#7bff7b", "#ceff29", "#ffc600", "#ff6800", "#f30900", "brown","#800000","black","pink","cyan","purple","magenta"),
                                          name = "Variants/kb",
                      labels = c("<0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25",
                                 "0.25-0.50","0.50-0.75","0.75-1.0","1.0-1.25","1.25-1.5",
                                 "1.5-1.75","1.75-2.0","2.0-2.25","2.25-2.5")) +  
    scale_x_continuous(name='Chromosome length', labels = c('0Mb',"50Mb", '100Mb', "150Mb", '200Mb','250Mb'),
                       breaks = c(0,50000000, 100000000, 150000000, 200000000,250000000), expand = c(0,0))
  
  ggsave(filename = paste0('Neofelis_',individual, '.1Mb.snpden.png'), plot = snpden_plot, device = 'png',
         dpi = 600, units = c('cm'), width = 28, height = 18, path = "plots/", bg = "white")}  
```
