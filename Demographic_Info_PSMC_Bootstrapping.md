
# Demographic History via PSMC

Demographic history in often investigated by applying PSMC (Pairwise Sequential Markovian Coalescent) analysis. To read more about the PSMC software we will use or to download, visit [PSMC](https://github.com/lh3/psmc/tree/master).

Apply PSMC through the following steps:

1) Consensus building
2) Run PSMC
3) Optional: use PSMC to run simulations based on results with MS software
4) Plot results

## 1) Consensus Building

Create the input file for PSMC by generating a consensus sequence from the .bam files.

Write and submit a job with the following code:
```Bash
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
 #
bcftools mpileup -Ou -f <reference_genome.fasta> <file.bam> | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > <output.fq.gz
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
#Explanation
#bcftools mpileup -C50 -uf <reference_genome> <bam_file>: generates a textual pileup format of the input BAM file using the given reference genome. 
#C50: applies a coefficient to adjust the base alignment quality
#-u : outputs the results in the uncompressed BCF format; required for piping to bcftools
#-f : specifies the reference genome file
#bcftools call -c: performs variant calling on the input data received from the bcftools mpileup command (indicated by `` as input). 
#-c: uses the consensus caller, which is suitable for calling a diploid consensus sequence
#vcfutils.pl vcf2fq -d 10 -D 100: converts the output from bcftools call in vcf format to a FastQ format
#d 10 and D 100: set the minimum and maximum depth thresholds for filtering variants
#gzip > <output.fq.gz>: compresses the final output using gzip and saves it as a .fq.gz file 
```

## 2) Run PSMC

First, we need to take a q-zipped fastq file and create a PSMC input file. You may be able to run this on an interactive or computational node. Otherwise, write and submit a job with the following code: 

```Bash

fq2psmcfa -q20 <file.fq.gz> > <output.psmcfa>

#Explanation
#-q20: mask bases with quality score below 20
```

Next, run PSMC on your data. Write and submit a job with the following code:
```Bash
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o <output.psmc> <input_file.psmcfa>
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
#Explanation
#N25: sets the effective population size (N) to 25. The effective population size is a measure of the genetic diversity in a population and is used to calculate the time to the most recent common ancestor of the population.
#t15: sets the scaled mutation rate per generation (t) to 15. The scaled mutation rate is the product of the mutation rate per base pair per generation and the effective population size.
#r5: sets the scaled recombination rate per generation (r) to 5. The scaled recombination rate is the product of the recombination rate per base pair per generation and the effective population size.
#p "4+25*2+4+6": This flag sets the time intervals (p) for the PSMC model. The specified pattern, "4+25*2+4+6", means that there are 4 intervals of equal size at the start, followed by 25 intervals with twice the size of the previous intervals, and then 4 more intervals of equal size, and finally 6 more intervals of increasing size. This allows the model to have higher time resolution near the present and lower resolution in the more distant past.
```

## 3) Optional: use PSMC to run simulations based on results with MS software

You can simulate genetic data using MS software, a tool for simulating coalescent-based gene genealogies and genetic data under various demographic scenarios.
Be sure modules are available and loaded on your HPC, and submit a job with the following code:
```Bash
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
psmc2history.pl <file.psmc> | history2ms.pl > <output_ms-cmd.sh>
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
#Explanation
#utils/psmc2history.pl <input.psmc>: This script takes the PSMC output file (<input.psmc>) and converts it into a simple history format. This format represents the inferred demographic history of the population.
#utils/history2ms.pl > <output_ms-cmd.sh>: This script takes the output from the previous script (the simple history format) and generates an ms command. The ms command is then saved in an output shell script file (<output_ms-cmd.sh>).
```

## 4) Plot results

Next, proceed to plotting the results. Consider generation time and mutation rate designated in the code below. Note that researchers may need to rely on estimates because obtaining this information can be challenging, especially for non-model species.

Be sure modules are available and loaded on your HPC, and submit a job with the following code:

```Bash
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
psmc_plot.pl -g 7 -u 1e-8 -p <input_name> </path/to/file/psmc/input.psmc>
#
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME
#
#Explanation
#-g 7: Maximum generation time in years (tmax) that will be used for the analysis; In this case, the value is 7 years
#-u 1e-8: Mutation rate per site per generation (u) used in the analysis; The value is set to 1x10^-8
#-p: save plot in pdf format
```
