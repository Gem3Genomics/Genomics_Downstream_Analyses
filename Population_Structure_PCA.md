# Population Structure via Principal Component Analysis

A **principal component analysis**, or **PCA** is a statistical procedure that allows you to summarize information in large datasets by means of smaller summary indices that can be visualized. PCAs identify main axes of variation in a dataset, and in genomics more specifically, identify main axes of allelic frequency. Individuals are then assigned coordinates and plotted along these axes, giving insight on population structure.

**PLINK** is a whole genome association toolset that includes a couple of options for dimension reduction, PCA and multidimensional scaling (MDS), along with many other tools and applications. To read more about PLINK's capabilities, or to download, visit [PLINK](https://zzz.bwh.harvard.edu/plink/tutorial.shtml). Apply PLINK in the following steps: 

1) Prepare data
2) Run PLINK
3) Visualize results in R

## 1) Prepare Data

Prepare data so that all assumptions within PCAs are met. A major assumption is that SNPs are independent, which is not the case with genomic data given that allelic frequencies are correlated due to linkage. Therefore, we need to remove variants that are in physical linkage.

Run the following code in an interactive or computational node:

```Bash
#
plink --vcf /path/to/vcf/file/file.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out file_plink
#
#Explanation
#--vcf: Path to the input VCF file.
#--double-id: this tells PLINK to set both the family ID and the within-family ID to be set to the sample ID
#--allow-extra-chr: allow additional chromosomes beyond the human chromosome set. PLINK by default expectes that the data is human
#--set-missing-var-ids: Set a variant ID for the SNPs. Human datas sets often have annotated SNP names, so plink will look for these. If not working wtih a human dataset, set to default to chromosome:position which can be achieved in plink by setting the option @:#.
#--indep-pairwise: This command performs the linkage pruning: "50" denotes a window of 50 Kb, "10" is the window step size - meaning we move 10 bp each time we calculate linkage. The third argument is an r2 threshold, or the threshold of linkage we are willing to tolerate. For example, variables that show an r2 of greater than 0.1 are pruned here.
```

## 2) Run PLINK

Run the following code in an interactive or computational node, or submit in a job script:

```Bash
#
plink --vcf /path/to/vcf/file/file.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract file.prune.in \
--make-bed --pca --out file_pca
#
#Explanation:
#--extract: this tells PLINK to extract only these positions from our .vcf file
#--make-bed: produces the output file can be used in downstream admixture analysis
#--pca: command to calculate a principal components analysis
```

Download the following files for the next step:

* file_pca.eigenval  
* file_pca.eigenvec

## 3) Visualize results in R

Plot PCA results with R and RStudio.

Create an R project and open a new Rmarkdown document or R script. A benefit to using markdown is that you can create a PDF or HTML document with your code to distribute to collaborators or colleagues. For more information on how to do this, visit [Rmarkdown](https://rmarkdown.rstudio.com). For a helpful cheat sheet on how to annotate your code so it "knits" nicely into a document, visit [RMarkdown Cheat Sheet](https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf).

```R
# Install tidyverse from CRAN. Tidyverse is a collection of R packages designed for data science and manipulation.
install.packages("tidyverse")

#load tidyverse library
library (tidyverse)

#Read our data using the read table function
pca <- read_table2("file_pca.eigenvec", col_names = FALSE)
eigenval <- scan("file_pca.eigenval")

#Set names for individuals
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

#Organize each individual by its designated name or ID, population, and species.Let's now organize each individual by its name and its population/species. For this you will create some vectors that will hold that information (i.e. nm for name and pop for population).

#Organize by name
nm <- rep(NA, length(pca$ind))
nm[grep("ID_001", pca$ind)] <- "name_1"
nm[grep("ID_002", pca$ind)] <- "name_2"
nm[grep("ID_003", pca$ind)] <- "name_3"

#Organize by Population
pop <- rep(NA, length(pca$ind))
pop[grep("ID_001", pca$ind)] <- "cl"
pop[grep("ID_002", pca$ind)] <- "cl"
pop[grep("ID_003", pca$ind)] <- "cl"

#Optional, organize by both color and population by creating a third vector:
nm_pop <- paste0(nm, "_", pop)

#Reorganize dataframe to a new one
pca2 <- as.tibble(data.frame(pca, nm, pop, nm_pop))

#Plot the percentage of variance by each principal component, convert your eigenvalues to percentages and create a new dataframe
pve <- data.frame(PC = 1:5, pve = eigenval/sum(eigenval)*100)

#Plot eignvalues and PCA after organizing
plot_pve <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
plot_pve <- plot_pve  + ylab("%  of variance expl") +
theme_classic()
plot_pve

#Calculate what percentage each principal component explains of variance:
cumsum_cal<-cumsum(pve$pve)
cumsum_cal

#Plot principal components
pca_plot <- ggplot(pca2, aes(PC1, PC2, col = nm, shape = pop)) + geom_point(size = 3)
pca_plot <- pca_plot + scale_colour_manual(values = c("red", "blue", "green", "purple","black"))
pca_plot <- pca_plot + coord_equal() + theme_light()
pca_plot <- pca_plot + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
pca_plot
#A simple way to interpret this is that the closer the data point, the more similar the individuals

```
