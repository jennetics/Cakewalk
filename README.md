# Conditional Analysis Cakewalk

Linkage disequilibrium (LD), the nonrandom association of alleles at different loci (Maruki & Lynch, 2014), is an important consideration when performing association analysis as it can bias the results. LD can cause some variants to appear significantly associated with a trait that conditional analysis will reveal to have no association. Therefore, ignoring LD can lead to a misinterpretation of genome wide association studies (GWAS), impacting conclusions or future studies derived from the results (Christoforou et al., 2012). As high throughput sequencing becomes more widely used, the challenge of managing the LD is growing. Conditional analysis can identify which variants are independently associated with an outcome (Knight et al., 2012), but running sequential conditional analyses can be tedious and time-consuming. To meet this challenge, we have developed the Conditional Analysis Cakewalk tool for automating this important process.

The Conditional Analysis Cakewalk begins with a section requiring input from the user. The user will enter:
- the desired trait
-	chromosome number
-	start/end positions for the region to be scanned
-	covariates 
- desired suffix for output files
- minor allele frequency (MAF) cutoff
-	four thresholds that can be adjusted:
    -	maximum p-values for including variants in both the initial and subsequent LD calculations 
    -	the p-value for establishing significance of variants as the termination condition 
    -	a p-value for limiting which variants are included in each association analysis, which was included to improve processing times 
-	the y-axis range and placement of significance lines within the LocusZoom 

Aditionally, the user will have the option of manually selecting the top variant for the initial analysis. IF this option is chosen, the top 10 variants will be printed to the output screen for selection. Finally, the results will be routed to a user-specified directory with an option for removing all other files to conserve storage space. 

The Conditional Analysis Cakewalk is available in two formats, utilizing either [MMAP](https://github.com/MMAP/MMAP-releases-issues-Q-and-A/releases/tag/mmap.2018_04_07_13_28.intel) or [PLINK](https://www.cog-genomics.org/plink2/) to conduct the initial association analysis between the specified trait and genomic region (O’Connell). The output file is then sorted and processed into a marker list for LD calculations, which are also conducted using MMAP and PLINK. The LD results are again sorted and processed into files formatted specifically for [LocusZoom](https://github.com/statgen/locuszoom-standalone). The LocusZoom plot includes all variants in the specified region, with a green vertical line in the location of the top variant. A boxplot will also be created with each iteration of the Cakewalk. The Conditional Analysis Cakewalk will continue to iterate through the MMAP analysis, LD calculations, LocusZoom and boxplot plotting until the p-value of the top variant is greater than the cutoff p-value specified in the beginning of the script.

## Running the Cakewalk
Instructions for running the Conditional Analysis Cakewalk using the simulated data are included below.
#### PLINK version
1. Download: 
    - PLINK_CAKEWALK.sh
    - simulated data PLINK.zip
2. Input information at the beginning of PLINK_CAKEWALK.sh file
3. Download PLINK and LocusZoom if needed at the links provided above or within PLINK_CAKEWALK.sh
4. If using a Unix machine, type "sh PLINK_CAKEWALK.sh" into the command line

#### MMAP version
1. Download:
    - MMAP_CAKEWALK.sh
    - simulated data MMAP.zip
    - dose_to_genotype.pl
2. Input information at the beginning of MMAP_CAKEWALK.sh file
3. Download PLINK, MMAP and LocusZoom if needed at the links provided above or within MMAP_CAKEWALK.sh
4. If using a Unix machine, type "sh MMAP_CAKEWALK.sh" into the command line


## References

Christoforou, A., Dondrup, M., Mattingsdal, M., Mattheisen, M., Giddaluru, S., Nothen, M.M., Rietschel, M., Cichon, S., Djurovic, S., Andreassen, O.A., Jonassen, I., Steen, V.M., Puntervoll, P., & Le Hellard, S. (2012). Linkage-disequilibrium-based binning affects the interpretation of GWASs. American Journal of Human Genetics, 90(4), 727-733.

Knight, J., Spain, S.L., Capon, F., Hayday, A., Nestle, F.O., Clop, A., Wellcome Trust Case Control Consortium, Genetic Analysis of Psoriasis Consortium, I-chip for Psoriasis Consortium, Barker, J.N., Weale, M.E., & Trembath, R.C. (2012). Conditional analysis identifies three novel major histocompatibility complex loci associated with psoriasis. Human Molecular Genetics, 21(23), 5185-5192. 

Maruki, T., & Lynch, M. (2014). Genome e-wide estimation of linkage disequilibrium from population-level high-throughput sequencing data. Genetics, 197(4), 1303-1313. 

O’Connell, J.R. MMAP: Mixed Model Analysis for Pedigrees and Populations. Available from: https://mmap.github.io/

Purcell, S, Neale, B, Todd-Brown, K, Thomas, L, Ferreira, MAR, Bender, D, Maller, J, Sklar, P, de Bakker, PIW, Daly, MJ, Sham PC. (2007). PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81. Available from https://www.cog-genomics.org/plink/1.9/#test_note. 

Welch, R., Pruim, R. (2010). LocusZoom standalone. Available from: https://github.com/statgen/locuszoom-standalone. 


