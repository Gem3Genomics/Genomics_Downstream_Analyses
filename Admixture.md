
# Admixture

Admixture is a program that calculates individual and population ancestries through maximum likelihood estimation from multilocus SNP genotype datasets. To read more about this program, visit [ADMIXTURE](https://dalexander.github.io/admixture/)

1) Prepare data
2) Run ADMIXTURE
3) Visualize Results in R

## 1) Prepare data, if applicable
Admixture can use input files from PLINK (see "Population Structure via PCA" files), including .bim and .bed files. However, PLINK was designed to work with human genomes, so if you are processing non-human data, you'll have to reorganize your data.

Manipulate this data on your command line interface:

```Bash
#Change the first column of the .bim file with a "0".
awk '{$1=0;print $0}' file.bim > file.bim.tmp

#Rename .bim file to match the name on the .bed file
mv file.bim.tmp file.bim
```

## 2) Run ADMIXTURE

Be sure programs are available and loaded on your HPC, and run the following code in an interactive or computational node:


```Bash
#tests ancestry assuming 2 genomic clusters
admixture --cv file.bed 2 > log2.out

#test admixture and ancestry on a range of K values to evaluate which K-value best fits the data; can be run in a loop
for i in {2..15}
do
admixture --cv efish.bed $i > log${i}.out
done
```

The output should be two files: a .Q file, which contains cluster assignments for each individual and a .P file which contains the population allele frequencies for each SNP. 

Identify the best K-value and cross-validation (cv) errors. There are three options of code you can apply:

Option 1:
```Bash
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'  > file.cv.error
```

Option 2:
```Bash
#grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > file.cv.error
```

Option 3:
```Bash
#awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > file.cv.error
```

Recommended step: extract a list of individuals from the analysis to best organize results.

```Bash
awk '{split($1,name,"."); print $1,name[2]}' file.nosex > file.list
```

## Visualize Results in R

```R
# read .Q file
tbl<-read.table("file.10.Q")

#Simple plot
barplot(t(as.matrix(tbl)), col=rainbow(10),
xlab="Individual #", ylab="Ancestry", border=NA)
```
