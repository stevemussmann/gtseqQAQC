# gtseqQAQC
Compare genotypes from two GTseq sequencing runs

## Dependencies
No dependencies (other than perl). 

## Instructions
Either install the script in your $PATH or place it in the same folder with your original GTseq genotypes file for a project and the QAQC genotypes file. The QAQC and original genotypes files should be the .csv outputs of Nate Campbell's GTseq_GenoCompile_v3.pl script (https://github.com/GTseq/GTseq-Pipeline). 

## Example
The following is an example command, assuming you have two files containing the samples from your original GTseq run (`original_genotypes.csv`) and your QAQC sequencing run (`qaqc_genotypes.csv`) which is expected to contain a subsample of the individuals in the `original_genotypes.csv` file.

```
./gtseqQAQC.pl -f original_genotypes.csv -q qaqc_genotypes.csv -o GT0012
```

Command explanations:
-f is used to input your original genotypes file
-q is used to specify your qaqc genotypes file
-o is used to specify a prefix that the script will use for naming the output files (e.g., `GT0012`).

The program will print to the screen a list of summary statistics for the amount of missing data in the two files, the proportion of detected genotype mismatches for individuals, and the proportion of detected genotype mismatches for loci. The script reports both 'raw' and 'adjusted' numbers of mismatches. The 'adjusted' values exclude missing data from their calculations. 

The program also produces two output files that are tab-delimited and can be opened in Excel. These files have the following columns:

`prefix.individuals.txt`
* Individual = sample
* Mismatch_Loci = number of loci for this individual that have mismatching genotypes
* Orig_Missing = number of loci for which comparisons could not be made because of missing data in this individual from the original GTseq run
* QAQC_Missing = number of loci for which comparisons could not be made because of missing data in this individual from the QAQC GTseq run
* Both_Missing = number of loci for which comparisons could not be made due to missing data in both runs
* Adjusted_Loci = number of loci successfully genotyped for this individual in both the original and QAQC runs
* Raw_Loci = number of loci common to both the original and QAQC runs
* Adjusted_%_Missing = (Mismatch_Loci / Adjusted_Loci) * 100
* Raw_%_Missing = (Mismatch_Loci / Raw_Loci) * 100


`prefix.loci.txt`
* Locus = GTseq locus
* Mismatch_Individuals = number of individual samples with mismatching genotypes
* Orig_Missing = number of individuals for which comparisons could not be made because of missing data at this locus from the original GTseq run
* QAQC_Missing = number of individuals for which comparisons could not be made because of missing data at this locus from the QAQC GTseq run
* Both_Missing = number of individuals for which comparisons could not be made due to missing data in both runs
* Adjusted_Individuals = number of individuals successfully genotyped at this locus in both the original and QAQC runs
* Raw_Individuals = number of individuals common to both the original and QAQC runs
* Adjusted_%_Missing = (Mismatch_Individuals / Adjusted_Individuals) * 100
* Raw_%_Missing = (Mismatch_Individuals / Raw_Individuals) * 100
