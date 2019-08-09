########################################################################
## Results can be found in file TOP_SNPS.csv. LocusZoom plots can be  ##
## found in pdf files listed in TOP_SNPS.csv. Please begin by         ##
## inputting data below.      										  ##
##########################input data here###############################

#Select suffix, trait, chromosome, minimum and maximum region
suffix=run1
trait=LDL
chr=19
min=11100000 
max=11300000

#Set location for genomic data
genomxs=/path/to/genomic/data

#Set location for binary data
bin_input_file=/path/to/binary/data

#Set location for covariate data
covfile=/path/to/covariate/data

#Does the data include a kinship file (0=yes, 1=no)?
kinship=1

#Set location for kinship file
kinbin=

#Set location for pedigree data
ped=/path/to/pedigree/data

#Set location for phenotype data - must be csv format
pheno=/path/to/phenotype/data

#Set location for MMAP - download at https://mmap.github.io/
myMMAP=/path/to/mmap

#Set location for PLINK - download at https://www.cog-genomics.org/plink2 
myPLINK=/path/to/plink

#Set location for LocusZoom - download at https://github.com/statgen/locuszoom-standalone
locuszoom=/path/to/locuszoom

#Set location for perl - download at https://www.perl.org/
perl=/path/to/perl

#Set up covariates: separate covariates with a "," 
covariates_list="AGE, SEX"

#Set minor allele frequency
MAF=0.01

#Set p-value for creating first markers list for LD calculation
pmark=0.01

#Set p-value for ending analysis (exponential form)
max_pvala=1E-4

#Set p-value for creating markers list in loop for LD calculations
pmark_loop=0.05

#Set p-value for cutting down SNP marker set after each MMAP analysis
pcut=0.10

#Set max value for y axis in LocusZoom plot
ymax=200

#Set significance lines for LocusZoom plot. Format "<upper_line>,<lower_line>"
sigline="6,8"

#Set name for results directory
dir=Results

#Set memory amount
mem=5000

#Is the variant name in chr:pos:ref:alt format (0=yes; 1=no)?
format=1

#Do you want the output other than results deleted? (0=yes; 1=no)
debug=0

#Do you want the top SNP to be automatically selected as the refSNP? (0=yes; 1=no; if selecting no, top SNPs will be output for you to select the refSNP for analysis)
select_SNP=0

#Are you using imputed genotype data (0=yes; 1=no)?
impute=1

###################beginning of script###################################
#create covariate file
echo $covariates_list > cov_col.csv
sed 's/\,/\t/g' cov_col.csv > covar.txt
covariates=`awk -F"\t" 'NR==1''{print; exit}' covar.txt` 

if [ $kinship -eq 0 ]
then
	$myMMAP --ped $ped --read_binary_covariance_file $kinbin --phenotype_filename $pheno --trait $trait --covariates $covariates --file_suffix $suffix --binary_genotype_filename $genomxs --covariate_filename $covfile --min_minor_allele_frequency $MAF --binary_covariate_filename $genomxs --genomic_region $chr $min $max --model add --all_output > $suffix.mmap.log
	
	#remove unnecessary columns
	awk -F"," '{OFS=","; print $1,$3,$4,$5,$6,$13,$17}' $trait.$suffix.add.mle.pval.csv > $trait.$suffix.csv
	
	#sort the add.mle.pval.csv file by pvalue	
	(head -n 1 $trait.$suffix.csv && tail -n +2 $trait.$suffix.csv | sort -gt, -k7,7) > sorted_$trait.$suffix.csv
	
	#create file for LD analysis
	awk -F"," '{OFS=","; print "chr"$2":"$3,$2,$3,$4,$5,$6,$7}' sorted_$trait.$suffix.csv > LD.$trait.$suffix.csv
	
	tail -n +2 LD.$trait.$suffix.csv > tail.$trait.$suffix.csv
	
	echo "SNPNAME,CHR,POS,NON_CODED_ALLELE,EFFECT_ALLELE,BETA_SNP,SNP_TSCORE_PVAL" > head.LD.csv
	
	cat head.LD.csv tail.$trait.$suffix.csv > LD.$trait.$suffix.csv
else
	$myMMAP --ped $ped --linear_regression --single_pedigree --phenotype_filename $pheno --trait $trait --covariates $covariates --file_suffix $suffix --binary_genotype_filename $genomxs --covariate_filename $covfile --min_minor_allele_frequency $MAF --binary_covariate_filename $genomxs --genomic_region $chr $min $max --model add --all_output > $suffix.mmap.log

	#remove unnecessary columns
	awk -F"," '{OFS=","; print $1,$3,$4,$5,$6,$27,$31}' $trait.$suffix.add.lm.pval.csv > $trait.$suffix.csv
	
	#sort the add.mle.pval.csv file by pvalue	
	(head -n 1 $trait.$suffix.csv && tail -n +2 $trait.$suffix.csv | sort -gt, -k7,7) > sorted_$trait.$suffix.csv
	
	#create file for LD analysis
	awk -F"," '{OFS=","; print "chr"$2":"$3,$2,$3,$4,$5,$6,$7}' sorted_$trait.$suffix.csv > LD.$trait.$suffix.csv
	
	tail -n +2 LD.$trait.$suffix.csv > tail.$trait.$suffix.csv
	
	echo "SNPNAME,CHR,POS,NON_CODED_ALLELE,EFFECT_ALLELE,BETA_SNP,SNP_TSCORE_PVAL" > head.LD.csv
	
	cat head.LD.csv tail.$trait.$suffix.csv > LD.$trait.$suffix.csv
fi

#loop based on format answer
if [ $format -eq 0 ]
then
	## print top hit SNPNAME - LZ
	if [ $select_SNP -eq 0 ]
	then
		## print top hit SNPNAME - LZ
		awk -F"," 'NR==2 {OFS=","; print $1}' sorted_$trait.$suffix.csv > pval_chrLZ.csv
	else
		list=`awk -F, '{print $1,$7,$6}' sorted_$trait.$suffix.csv | tail -n + 2 | head -n 10 | cat --number`
		echo "Enter the line number of the desired variant to select as the reference variant (columns are SNP, p-value, beta):
		$list"
		read line_num
		awk -F"," -v line=$line_num 'NR==line {OFS=","; print $1}' sorted_$trait.$suffix.csv > pval_chrLZ.csv
	fi
	####add pvalue and beta - LZ 	
	if [ $select_SNP -eq 0 ]
	then	
		cpos=`awk -F"," 'NR==2 {OFS=","; print $1}' LD.$trait.$suffix.csv`
		awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' pval_chrLZ.csv > pval_chrLZ2.csv
		pv=`awk -F"," 'NR==2 {OFS=","; print $7}' sorted_$trait.$suffix.csv`
		awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' pval_chrLZ2.csv > p_valLZ.csv
		beta=`awk -F"," 'NR==2 {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
		awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' p_valLZ.csv > pvalLZ.csv
	else
		cpos=`awk -F"," -v line=$line_num 'NR==line {OFS=","; print $1}' LD.$trait.$suffix.csv`
		awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' pval_chrLZ.csv > pval_chrLZ2.csv
		pv=`awk -F"," -v line=$line_num 'NR==line {OFS=","; print $7}' sorted_$trait.$suffix.csv`
		awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' pval_chrLZ.csv > p_valLZ.csv
		beta=`awk -F"," -v line=$line_num 'NR==line {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
		awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' p_valLZ.csv > pvalLZ.csv
	fi
	
	#add header to pval.csv - LZ
	echo "SNPNAME, CHR:POS, PVALUE, BETA" > pval.header
	cat pval.header pvalLZ.csv > TOP_SNPS.csv

else
	###TOP_SNPS for LD calc#####
	if [ $select_SNP -eq 0 ]
	then
		awk -F"," 'NR==2 {OFS=","; print $1,"chr"$2":"$3,$7}' sorted_$trait.$suffix.csv > pval_LD.csv
		beta=`awk -F"," 'NR==2 {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
		awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' pval_LD.csv > pvalLD.csv
	else
		awk -F"," -v line=$line_num 'NR==line {OFS=","; print $1,chr$2":"$3,$7}' sorted_$trait.$suffix.csv > pval_LD.csv
		beta=`awk -F"," -v line=$line_num 'NR==line {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
		awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' pval_LD.csv > pvalLD.csv
	fi

	#add header#
	echo "SNPNAME, CHR:POS, PVALUE, BETA" > pval.header
	cat pval.header pvalLD.csv > TOP_SNPS.csv
fi	

###create markers.list.text
awk -F, -v pmark=$pmark '$7 < pmark { print $1}' sorted_$trait.$suffix.csv > test1.marker.list.txt

############ use MMAP to calculate LD with top hit as ref variant#################
name=test1
listLD=$name.marker.list.txt
calcLD_refsnp=`awk -F, 'NR==2 {print $1}' TOP_SNPS.csv`

# BASIC DEFINITIONS:

datasetLD=$name.subset
resultsLD=$name.results
plinkOut1=$name.plinkOut1.log
mmapOut1=$name.mmapOut1.log
mmapOut2=$name.mmapOut2.log
mmapOut3=$name.mmapOut3.log
mmapOut4=$name.mmapOut4.log
skip_fields=8  # number of columns to skip before individual genotypes begin

#loop based on impute answer
if [ $impute -eq 0 ];
then	
	# STEP1: subset the genotype file with a marker_set and write out as a csv file
	$myMMAP --marker_by_subject_mmap2csv --binary_input_filename $bin_input_file --csv_output_filename $datasetLD.csv --min_minor_allele_frequency $MAF --marker_set $listLD --file_suffix $name.1 > $mmapOut1 2>&1
	
	# STEP2: change imputed dosages to genotype values of 0, 1 or 2
	#        This step can be skipped if genotypes are not imputed dosages
	$perl dose_to_genotype.pl $datasetLD.csv $skip_fields
	cp $datasetLD.csv $datasetLD.original.csv
	mv $datasetLD.csv.genotypes $datasetLD.csv

	# STEP3: convert csv to mmap MxS.bin
	$myMMAP --write_binary_genotype_file --csv_input_filename $datasetLD.csv --binary_output_filename $datasetLD.MxS.bin --num_skip_fields $skip_fields --file_suffix $name.2 > $mmapOut2 2>&1

	# STEP4: transpose MxS.bin to SxM.bin
	$myMMAP --transpose_binary_genotype_file --binary_input_filename $datasetLD.MxS.bin --binary_output_filename $datasetLD.SxM.bin --file_suffix $name.3 > $mmapOut3 2>&1

	# STEP5: create plink dataset from the SxM.bin
	$myMMAP --subject_by_marker_mmap2plink --binary_input_filename $datasetLD.SxM.bin --plink_output_prefix $datasetLD.plink --use_snpname --file_suffix $name.4 > $mmapOut4 2>&1

	# STEP6: run Plink to perform the LD calculations on everything in the "subset" dataset
	$myPLINK --file $datasetLD.plink --r2 dprime --ld-snp $calcLD_refsnp --make-founders --maf $MAF --out $resultsLD --ld-window 999999999 --ld-window-kb 999999999 --ld-window-r2 0 > $plinkOut1 2>&1
else
	# STEP4: transpose MxS.bin to SxM.bin
	$myMMAP --transpose_binary_genotype_file --binary_input_filename $bin_input_file --min_minor_allele_frequency $MAF --binary_output_filename $datasetLD.SxM.bin --file_suffix $name.3 > $mmapOut3 2>&1

	# STEP5: create plink dataset from the SxM.bin
	$myMMAP --subject_by_marker_mmap2plink --binary_input_filename $datasetLD.SxM.bin --plink_output_prefix $datasetLD.plink --use_snpname --file_suffix $name.4 > $mmapOut4 2>&1

	# STEP6: run Plink to perform the LD calculations on everything in the "subset" dataset
	$myPLINK --file $datasetLD.plink --r2 dprime --ld-snp $calcLD_refsnp --make-founders --maf $MAF --extract $listLD --out $resultsLD --ld-window 999999999 --ld-window-kb 999999999 --ld-window-r2 0 > $plinkOut1 2>&1
fi

#sort result file
(head -n 1 $name.results.ld && tail -n +2 $name.results.ld | sort -gr -k7,7) > sorted.$name.results.ld

#convert pheno file from csv to txt
cat $pheno | tr "," "\\t" > pheno.txt

#create boxplot
calcLD_refsnp=`awk -F, 'NR==2 {print $1}' TOP_SNPS.csv`

#create genotype file
$myPLINK --file $datasetLD.plink --out snp_genotypes --recode A --snp $calcLD_refsnp

ref=`awk -F"," 'NR==2 {print $4}' sorted_$trait.$suffix.csv`
alt=`awk -F"," 'NR==2 {print $5}' sorted_$trait.$suffix.csv`

ref_l=`echo -n $ref | wc -c`
alt_l=`echo -n $alt | wc -c`

if [ $ref_l -gt 1 ]
then 
	ref="I"
	alt="D"
fi
if [ $alt_l -gt 1 ] 
then 
	ref="D"
	alt="I"
fi

ref_homo=$ref$ref
hetero=$ref$alt
alt_homo=$alt$alt

#create sample size variables
#print column with genotypes
awk '{print $7}' snp_genotypes.raw | tail -n +2 > snp.raw
suba=`grep "0" snp.raw | wc | awk '{print $1}'`
subb=`grep "1" snp.raw | wc | awk '{print $1}'`
subc=`grep "2" snp.raw | wc | awk '{print $1}'`

#check for homozygotes for alt allele
awk '{print $7}' snp_genotypes.raw | tail -n +2 > col7.txt

grep "2" col7.txt > althomo.txt

calcLD_refsnp=`awk -F, 'NR==2 {print $1}' TOP_SNPS.csv`

if [ -s althomo.txt ]
then
	#make boxplot
	echo "geno1 = read.table('"snp_genotypes.raw"',header=T, fill = TRUE)
	geno2 = geno1[order(geno1$IID),]
	pheno1 = read.table('"pheno.txt"',header=T, fill = TRUE)
		
	geno_id=geno2[,2]
	pheno_id=pheno1[,1]
		
	common_subjects = intersect(geno_id,pheno_id)
		
	geno = geno2[geno_id %in% common_subjects,]
	pheno = pheno1[pheno_id %in% common_subjects,]
		
	colnames(geno)[7]='"genotype"'
	
	pdf('"$trait.$name.pdf"')
	par(mgp=c(3,1.5,0))
	boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects','"$alt_homo" \n "$subc" subjects'))
	stripchart(pheno$"$trait"~geno$"genotype", vertical=TRUE, add=TRUE, col='"gray39"', pch=1)
	status=dev.off()" > test.R
else 
	echo "geno1 = read.table('"snp_genotypes.raw"',header=T, fill = TRUE)
	geno2 = geno1[order(geno1$IID),]
	pheno1 = read.table('"pheno.txt"',header=T, fill = TRUE)
		
	geno_id=geno2[,2]
	pheno_id=pheno1[,1]
		
	common_subjects = intersect(geno_id,pheno_id)
		
	geno = geno2[geno_id %in% common_subjects,]
	pheno = pheno1[pheno_id %in% common_subjects,]
		
	colnames(geno)[7]='"genotype"'
	
	pdf('"$trait.$name.pdf"')
	par(mgp=c(3,1.5,0))
	boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects'))
	stripchart(pheno$"$trait"~geno$"genotype", vertical=TRUE, add=TRUE, col='"gray39"', pch=1)
	status=dev.off()" > test.R
fi

Rscript test.R 

echo $calcLD_refsnp > allrefs.txt

####create epacts & ld files####
if [ $format -eq 0 ]
then
	# remove letters 
	f=sorted_$trait.$suffix.csv
	awk -F, '{OFS="\t"; print $2,$3,$3,$1,$7}' $f | tail -n +1 > epacts.INPUT.txt
	sed 's/[A-Z]*//g' epacts.INPUT.txt > epacts.INPUT_cut.txt
	sed -i 's/:://g' epacts.INPUT_cut.txt
	sed -i '1d' epacts.INPUT_cut.txt

	#add prefix 
	prefix="chr"
	awk -v prefix=$prefix '{print $1,$2,$3,prefix $4,$5}' epacts.INPUT_cut.txt > epacts.INPUT_chr.txt
	
	#add in headers
	echo "#CHROM	BEGIN	END	MARKER_ID	PVALUE" > epacts.header
	cat epacts.header epacts.INPUT_chr.txt > epac.$name.INPUT_final.txt.epacts

	####create ld file####
	ldfile=sorted.$name.results.ld 
	awk '{OFS="\t"; print $3,$6,$8,$7}' $ldfile > locuszoom.INPUT.txt.ld

	# remove letters
	sed 's/[A-Z]*//g' locuszoom.INPUT.txt.ld > locuszoom.INPUT_cut.txt.ld
	sed -i 's/:://g' locuszoom.INPUT_cut.txt.ld
	sed -i '1d' locuszoom.INPUT_cut.txt.ld

	#add prefix
	prefix="chr"
	awk -v prefix=$prefix '{print prefix $1,prefix $2,$3,$4}' locuszoom.INPUT_cut.txt.ld > locuszoom.INPUT_chr.txt.ld

	#add in headers
	echo "snp1	snp2	dprime	rsquare" > ld.header
	cat ld.header locuszoom.INPUT_chr.txt.ld > locuszoom.$name.INPUT_final.txt.ld
else
	#create epacts file
	f=sorted_$trait.$suffix.csv
	awk -F, '{OFS="\t"; print $2,$3,$3,"chr"$2":"$3,$7}' $f | tail -n +2 > epacts.INPUT.txt
	
	#add headers
	echo "#CHROM	BEGIN	END	MARKER_ID	PVALUE" > epacts.header
	cat epacts.header epacts.INPUT.txt > epac.$name.INPUT_final.txt.epacts
	
	####create ld file####
	ldfile=sorted.$name.results.ld 
	awk '{OFS="\t"; print "chr"$1":"$2,"chr"$4":"$5,$8,$7}' $ldfile | tail -n +2 > locuszoom.INPUT.txt.ld
	
	#add in headers
	echo "snp1	snp2	dprime	rsquare" > ld.header
	cat ld.header locuszoom.INPUT.txt.ld > locuszoom.$name.INPUT_final.txt.ld
fi

##make files tab delimited###
sed 's/ /\t/g' epac.test1.INPUT_final.txt.epacts > ts.epac.test1.INPUT_final.txt.epacts
sed 's/ /\t/g' locuszoom.$name.INPUT_final.txt.ld > ts.locuszoom.test1.INPUT_final.txt.ld

## create variable for vertical line in LocusZoom plot
vertline=`awk -F, 'NR==2 {print $3}' sorted_$trait.$suffix.csv`

exit

########################## generate LocusZoom plot #######################################
LZ_refsnp=`awk -F, 'NR==2 {print $2}' TOP_SNPS.csv`

$locuszoom  --build hg38 --epacts ts.epac.test1.INPUT_final.txt.epacts --ld ts.locuszoom.test1.INPUT_final.txt.ld --refsnp $LZ_refsnp --chr $chr --start $min --end $max geneFontSize=.5 ymax=$ymax vertline=$vertline vertLineColor="green4" vertLineType="1" vertLineWidth="1" signifLine=$sigline signifLineColor="red,royalblue" title="$trait on Chr$chr" titleCex=1.3 title2="Top Variant: $LZ_refsnp Covariates: $covariates" title2Cex=1.0 title2Color="gray30" showRecomb=TRUE recombColor="gray80" recombAxisColor="gray60" --plotonly --cache None --no-date --prefix "LZ_$trait.$name.ImpRgn3TM5"  >> LZdialog.log 2>&1 

##move file name onto doc##
awk -F"," 'NR==2 {OFS=","; print}' TOP_SNPS.csv > pvalLD.csv

lzname="LZ_$trait.$name_$calcLD_refsnp.pdf"
awk -F"," -v lz=$lzname '{OFS=","} {$5=lz; print}' pvalLD.csv > pval_withLD.csv
bpname="$trait.$name.pdf"
awk -F"," -v bp=$bpname '{OFS=","; $6=bp; print}' pval_withLD.csv > pval_withFiles.csv

#add header#
echo "SNPNAME, CHR:POS, PVALUE, BETA, LOCUSZOOM, BOXPLOT" > pval_withFiles.header
cat pval_withFiles.header pval_withFiles.csv > TOP_SNPS.csv

#cut down markers for inclusion in further analysis
awk -F, -v pcut=$pcut '$7 < $pcut { print $1}' sorted_$trait.$suffix.csv > markers.csv

#################################################################################################################################
############################################### start of loop ###################################################################
#################################################################################################################################

a=2

#switch row to column
sed 's/ \+/,/g' covar.txt > top_row.csv
tr -s ',' '\n' < top_row.csv | tr -d '"' > top_col.csv

cp top_col.csv col_LZ.csv

#add to cov list
awk -v a=$a -F"," 'NR==a {print $1}' TOP_SNPS.csv >> top_col.csv
awk 'BEGIN {ORS = ","}{print}' top_col.csv > covar_new.csv
sed 's/\,/\t/g' covar_new.csv > cov_new.txt
sed 's/ /\	/g' cov_new.txt > covar.txt

#make LZ covar file for title 
LZ_refsnp=`awk -F, 'NR==2 {print $2}' TOP_SNPS.csv`
echo $LZ_refsnp >> col_LZ.csv
awk 'BEGIN {ORS=","}{ print}' col_LZ.csv > cov_LZ.csv
sed 's/\,/\ /g' cov_LZ.csv > LZ_cov.txt
	
##while p-value < max 
pval_dec=`awk -F, 'NR==2 {print $3}' TOP_SNPS.csv `

#cut decimals off of pvals
pval=`printf "%.0g" "$pval_dec"`
max_pval=`printf "%.0g" "$max_pvala"`

r=$(awk 'BEGIN {print ('${pval^^}'<'${max_pval^^}')?1:0}')

while [ "$r" == 1 ]
do
	#run analysis with SNP added to covariates
	markers=markers.csv
	covariates=`awk -F"\t" 'NR==1''{print; exit}' covar.txt`
	suffix=run$a
	if [ $kinship -eq 0 ]
	then
		$myMMAP --ped $ped --read_binary_covariance_file $kinbin --phenotype_filename $pheno --trait $trait --covariates $covariates --file_suffix $suffix --binary_genotype_filename $genomxs --covariate_filename $covfile --min_minor_allele_frequency $MAF --binary_covariate_filename $genomxs --genomic_region $chr $min $max --model add --all_output > $suffix.mmap.log
	
		#remove unnecessary columns
		awk -F"," '{OFS=","; print $1,$3,$4,$5,$6,$13,$17}' $trait.$suffix.add.mle.pval.csv > $trait.$suffix.csv
	
		#remove NA, nan, and -nan rows
		awk -F"," '$7 != "NA" && $7 != "nan" && $7 != "-nan"' $trait.$suffix.csv > $trait.$suffix.clean.csv 
		
		#sort the add.mle.pval.csv file by pvalue	
		(head -n 1 $trait.$suffix.csv && tail -n +2 $trait.$suffix.clean.csv | sort -gt, -k7,7) > sorted_$trait.$suffix.csv
		
		#create file for LD analysis
		awk -F"," '{OFS=","; print "chr"$2":"$3,$2,$3,$4,$5,$6,$7}' sorted_$trait.$suffix.csv > LD.$trait.$suffix.csv
	
		tail -n +2 LD.$trait.$suffix.csv > tail.$trait.$suffix.csv
	
		echo "SNPNAME,CHR,POS,NON_CODED_ALLELE,EFFECT_ALLELE,BETA_SNP,SNP_TSCORE_PVAL" > head.LD.csv
	
		cat head.LD.csv tail.$trait.$suffix.csv > LD.$trait.$suffix.csv
	else
		$myMMAP --ped $ped --linear_regression --single_pedigree --phenotype_filename $pheno --trait $trait --covariates $covariates --file_suffix $suffix --binary_genotype_filename $genomxs --covariate_filename $covfile --min_minor_allele_frequency $MAF --binary_covariate_filename $genomxs --genomic_region $chr $min $max --model add --all_output > $suffix.mmap.log

		#remove unnecessary columns
		awk -F"," '{OFS=","; print $1,$3,$4,$5,$6,$27,$31}' $trait.$suffix.add.lm.pval.csv > $trait.$suffix.csv
	
		#remove NA, nan, and -nan rows
		awk -F"," '$7 != "NA" && $7 != "nan" && $7 != "-nan"' $trait.$suffix.csv > $trait.$suffix.clean.csv 
		
		#sort the add.mle.pval.csv file by pvalue	
		(head -n 1 $trait.$suffix.csv && tail -n +2 $trait.$suffix.clean.csv | sort -gt, -k7,7) > sorted_$trait.$suffix.csv
		#create file for LD analysis
		awk -F"," '{OFS=","; print "chr"$2":"$3,$2,$3,$4,$5,$6,$7}' sorted_$trait.$suffix.csv > LD.$trait.$suffix.csv
	
		tail -n +2 LD.$trait.$suffix.csv > tail.$trait.$suffix.csv
	
		echo "SNPNAME,CHR,POS,NON_CODED_ALLELE,EFFECT_ALLELE,BETA_SNP,SNP_TSCORE_PVAL" > head.LD.csv
	
		cat head.LD.csv tail.$trait.$suffix.csv > LD.$trait.$suffix.csv
	fi
	
	#loop based on format
	if [ $format -eq 0 ]
	then
		## print top hit SNPNAME - LZ
		awk -F"," 'NR==2 {OFS=","; print $2}' sorted_$trait.$suffix.csv > newp_valueLZ.csv
	
		# remove letters
		sed 's/:[^:]*//2g' newp_valueLZ.csv > newpval_cutLZ.csv

		#add prefix
		prefix="chr"
		awk -F, -v prefix=$prefix '{print prefix $1}' newpval_cutLZ.csv > newpval_chrLZ.csv
	
		####add pvalue and beta - LZ
		cpos=`awk -F"," 'NR==2 {OFS=","; print $1}' LD.$trait.$suffix.csv`
		awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' newpval_chrLZ.csv > newpval_chrLZ2.csv
		pv=`awk -F"," 'NR==2 {OFS=","; print $7}' sorted_$trait.$suffix.csv`
		awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' newpval_chrLZ2.csv > newp_valLZ.csv
		beta=`awk -F"," 'NR==2 {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
		awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' newp_valLZ.csv > newpvalLZ.csv
		awk -F"," 'NR=1 {print; exit}' newpvalLZ.csv >> TOP_SNPS.csv
	else
		#extract SNPname and pval
		awk -F"," 'NR==2 {OFS=","; print $1,"chr"$2":"$3,$7}' sorted_$trait.$suffix.csv > pval_new.csv
		beta=`awk -F"," 'NR==2 {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
		awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' pval_new.csv > pval_new2.csv
		awk -F"," 'NR==1 {print; exit}' pval_new2.csv >> TOP_SNPS.csv
	fi
	
	###create markers.list.text  
	name=test$a
	awk -F, -v pmark_loop=$pmark_loop '$7 < pmark_loop {print $1}' sorted_$trait.$suffix.csv > $name.marker.list.txt
	
	#End script if file size for markers list is 0 
	if [ `du -b $name.marker.list.txt | awk '{print $1}'` -eq 0 ];
	then
		sed -i '$d' TOP_SNPS.csv
		sed -i '$d' TOP_SNPS_LZ.csv
		break
	fi

	#End script if new pval is bigger than cutoff
	pval_newdec=`awk -F, 'END {print $3}' TOP_SNPS.csv`
	pval_new=`printf "%.0g" "$pval_newdec"`
	t=$(awk 'BEGIN {print ('${pval_new^^}'<'${max_pval^^}')?1:0}')
	
	if [ "$t" != 1 ];
	then
		sed -i '$d' TOP_SNPS.csv
		sed -i '$d' TOP_SNPS_LZ.csv
		break
	fi
	
	#create b variable
	b=`expr $a + 1`
	
	#stop if pval=0
	zero_test=`awk -F"," -v b=$b 'NR==b {print $3}' TOP_SNPS.csv`
	if [ $zero_test == 0 ];
	then 
		sed -i '$d' TOP_SNPS.csv
		sed -i '$d' TOP_SNPS_LZ.csv
		break
	fi		
	
	##########################calculate LD####################################
	listLD=$name.marker.list.txt
	calcLD_refsnp=`awk -F, 'NR==2 {print $1}' sorted_$trait.$suffix.csv`

	# BASIC DEFINITIONS:
	datasetLD=$name.subset
	resultsLD=$name.results
	plinkOut1=$name.plinkOut1.log
	mmapOut1=$name.mmapOut1.log
	mmapOut2=$name.mmapOut2.log
	mmapOut3=$name.mmapOut3.log
	mmapOut4=$name.mmapOut4.log
	skip_fields=8  # number of columns to skip before individual genotypes begin
	
	#loop based on impute answer
	if [ $impute -eq 0 ];
	then	
		# STEP1: subset the genotype file with a marker_set and write out as a csv file
		$myMMAP --marker_by_subject_mmap2csv --binary_input_filename $bin_input_file --csv_output_filename $datasetLD.csv --min_minor_allele_frequency $MAF --marker_set $listLD --file_suffix $name.1 > $mmapOut1 2>&1
	
		# STEP2: change imputed dosages to genotype values of 0, 1 or 2
		#        This step can be skipped if genotypes are not imputed dosages
		$perl dose_to_geno.pl $datasetLD.csv $skip_fields
		cp $datasetLD.csv $datasetLD.original.csv
		mv $datasetLD.csv.genotypes $datasetLD.csv

		# STEP3: convert csv to mmap MxS.bin
		$myMMAP --write_binary_genotype_file --csv_input_filename $datasetLD.csv --binary_output_filename $datasetLD.MxS.bin --num_skip_fields $skip_fields --file_suffix $name.2 > $mmapOut2 2>&1

		# STEP4: transpose MxS.bin to SxM.bin
		$myMMAP --transpose_binary_genotype_file --binary_input_filename $datasetLD.MxS.bin --binary_output_filename $datasetLD.SxM.bin --file_suffix $name.3 > $mmapOut3 2>&1

		# STEP5: create plink dataset from the SxM.bin
		$myMMAP --subject_by_marker_mmap2plink --binary_input_filename $datasetLD.SxM.bin --plink_output_prefix $datasetLD.plink --use_snpname --file_suffix $name.4 > $mmapOut4 2>&1

		# STEP6: run Plink to perform the LD calculations on everything in the "subset" dataset
		$myPLINK --file $datasetLD.plink --r2 dprime --ld-snp $calcLD_refsnp --make-founders --maf $MAF --out $resultsLD --ld-window 999999999 --ld-window-kb 999999999 --ld-window-r2 0 > $plinkOut1 2>&1
	else
		# STEP4: transpose MxS.bin to SxM.bin
		$myMMAP --transpose_binary_genotype_file --binary_input_filename $bin_input_file --min_minor_allele_frequency $MAF --binary_output_filename $datasetLD.SxM.bin --file_suffix $name.3 > $mmapOut3 2>&1

		# STEP5: create plink dataset from the SxM.bin
		$myMMAP --subject_by_marker_mmap2plink --binary_input_filename $datasetLD.SxM.bin --plink_output_prefix $datasetLD.plink --use_snpname --file_suffix $name.4 > $mmapOut4 2>&1

		# STEP6: run Plink to perform the LD calculations on everything in the "subset" dataset
		$myPLINK --file $datasetLD.plink --r2 dprime --ld-snp $calcLD_refsnp --make-founders --maf $MAF --extract $listLD --out $resultsLD --ld-window 999999999 --ld-window-kb 999999999 --ld-window-r2 0 > $plinkOut1 2>&1
	fi
	
	#sort result file
	(head -n 1 $name.results.ld && tail -n +2 $name.results.ld | sort -gr -k7,7) > sorted.$name.results.ld
	
	#convert pheno file from csv to txt
	cat $pheno | tr "," "\\t" > pheno.txt

	#create boxplot
	#create genotype file
	$myPLINK --file $datasetLD.plink --out snp_genotypes.$name --recode A --snp $calcLD_refsnp

	ref=`awk -F"," 'NR==2 {print $4}' sorted_$trait.$suffix.csv`
	alt=`awk -F"," 'NR==2 {print $5}' sorted_$trait.$suffix.csv`

	
	ref_l=`echo -n $ref | wc -c`
	alt_l=`echo -n $alt | wc -c`

	if [ $ref_l -gt 1 ]
	then 
		ref="I"
		alt="D"
	fi
	if [ $alt_l -gt 1 ] 
	then 
		ref="D"
		alt="I"
	fi

	ref_homo=$ref$ref
	hetero=$ref$alt
	alt_homo=$alt$alt

	#create sample size variables
	#print column with genotypes
	awk '{print $7}' snp_genotypes.$name.raw | tail -n +2 > $name.snp.raw
	suba=`grep "0" $name.snp.raw | wc | awk '{print $1}'`
	subb=`grep "1" $name.snp.raw | wc | awk '{print $1}'`
	subc=`grep "2" $name.snp.raw | wc | awk '{print $1}'`

	#check for homozygotes for alt allele
	awk '{print $7}' snp_genotypes.$name.raw | tail -n +2 > col7.$name.txt

	grep "2" col7.$name.txt > althomo.$name.txt

	calcLD_refsnp=`awk -F, 'NR==2 {print $1}' TOP_SNPS.csv`

	if [ -s althomo.$name.txt ]
	then
		#make boxplot
		echo "geno1 = read.table('"snp_genotypes.$name.raw"',header=T, fill = TRUE)
		geno2 = geno1[order(geno1$IID),]
		pheno1 = read.table('"pheno.txt"',header=T, fill = TRUE)
		
		geno_id=geno2[,2]
		pheno_id=pheno1[,1]
		
		common_subjects = intersect(geno_id,pheno_id)
		
		geno = geno2[geno_id %in% common_subjects,]
		pheno = pheno1[pheno_id %in% common_subjects,]
		
		colnames(geno)[7]='"genotype"'
	
		pdf('"$trait.$name.pdf"')
		par(mgp=c(3,1.5,0))
		boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects','"$alt_homo" \n "$subc" subjects'))
		stripchart(pheno$"$trait"~geno$"genotype", vertical=TRUE, add=TRUE, col='"gray39"', pch=1)
		status=dev.off()" > test.$name.R
	else 
		echo "geno1 = read.table('"snp_genotypes.$name.raw"',header=T, fill = TRUE)
		geno2 = geno1[order(geno1$IID),]
		pheno1 = read.table('"pheno.txt"',header=T, fill = TRUE)
			
		geno_id=geno2[,2]
		pheno_id=pheno1[,1]
		
		common_subjects = intersect(geno_id,pheno_id)
		
		geno = geno2[geno_id %in% common_subjects,]
		pheno = pheno1[pheno_id %in% common_subjects,]
		
		colnames(geno)[7]='"genotype"'
	
		pdf('"$trait.$name.pdf"')
		par(mgp=c(3,1.5,0))
		boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects'))
		stripchart(pheno$"$trait"~geno$"genotype", vertical=TRUE, add=TRUE, col='"gray39"', pch=1)
		status=dev.off()" > test.$name.R
	fi
	
	Rscript test.$name.R 

	echo $calcLD_refsnp >> allrefs.txt
	
	####create epacts file####
	if [ $format -eq 0 ]
	then
		# remove letters
		f=sorted_$trait.$suffix.csv
		awk -F, '{OFS="\t"; print $2,$3,$3,$1,$7}' $f | tail -n +1 > epacts.INPUT.txt
		sed 's/[A-Z]*//g' epacts.INPUT.txt > epacts.INPUT_cut.txt
		sed -i 's/:://g' epacts.INPUT_cut.txt
		sed -i '1d' epacts.INPUT_cut.txt

		#add prefix
		prefix="chr"
		awk -v prefix=$prefix '{print $1,$2,$3,prefix $4,$5}' epacts.INPUT_cut.txt > epacts.INPUT_chr.txt

		#add in headers
		echo "#CHROM	BEGIN	END	MARKER_ID	PVALUE" > epacts.header
		cat epacts.header epacts.INPUT_chr.txt > epac.$name.INPUT_final.txt.epacts

		####create ld file####
		ldfile=sorted.$name.results.ld 
		awk '{OFS="\t"; print $3,$6,$8,$7}' $ldfile > locuszoom.INPUT.txt.ld

		# remove letters
		sed 's/[A-Z]*//g' locuszoom.INPUT.txt.ld > locuszoom.INPUT_cut.txt.ld
		sed -i 's/:://g' locuszoom.INPUT_cut.txt.ld
		sed -i '1d' locuszoom.INPUT_cut.txt.ld

		#add prefix
		prefix="chr"
		awk -v prefix=$prefix '{print prefix $1,prefix $2,$3,$4}' locuszoom.INPUT_cut.txt.ld > locuszoom.INPUT_chr.txt.ld

		#add in headers
		echo "snp1	snp2	dprime	rsquare" > ld.header
		cat ld.header locuszoom.INPUT_chr.txt.ld > locuszoom.$name.INPUT_final.txt.ld
		
		#format covariates file LZ- switch from row to columns
		sed 's/ \+/,/g' covar.txt > top_row_LZ.csv
		tr -s ',' '\n' < top_row_LZ.csv | tr -d '"' > top_col_LZ.csv
	else
		#create epacts file
		f=sorted_$trait.$suffix.csv
		awk -F, '{OFS="\t"; print $2,$3,$3,"chr"$2":"$3,$7}' $f | tail -n +2 > epacts.INPUT.txt
		
		#add headers
		echo "#CHROM	BEGIN	END	MARKER_ID	PVALUE" > epacts.header
		cat epacts.header epacts.INPUT.txt > epac.$name.INPUT_final.txt.epacts
	
		####create ld file####
		ldfile=sorted.$name.results.ld 
		awk '{OFS="\t"; print "chr"$1":"$2,"chr"$4":"$5,$8,$7}' $ldfile | tail -n +2 > locuszoom.INPUT.txt.ld
	
		#add in headers
		echo "snp1	snp2	dprime	rsquare" > ld.header
		cat ld.header locuszoom.INPUT.txt.ld > locuszoom.$name.INPUT_final.txt.ld
	fi
	
	##make files tab delimited###
	sed 's/ /\t/g' epac.$name.INPUT_final.txt.epacts > ts.epac.$name.INPUT_final.txt.epacts
	sed 's/ /\t/g' locuszoom.$name.INPUT_final.txt.ld > ts.locuszoom.$name.INPUT_final.txt.ld
	
	covariates_LZ=`awk -Ft 'NR==1 {print}' LZ_cov.txt`
	
	#generate LocusZoom
	epacts=ts.epac.$name.INPUT_final.txt.epacts
	txtld=ts.locuszoom.$name.INPUT_final.txt.ld
	LZ_refsnp=`awk -F, -v b=$b 'NR==b {print $2}' TOP_SNPS.csv`
	
	$locuszoom --build hg38 --epacts $epacts --ld $txtld --refsnp $LZ_refsnp --chr $chr --start $min --end $max geneFontSize=.5 ymax=$ymax vertline="$vertline" vertLineColor="green4" vertLineType="1" vertLineWidth="1" signifLine=$sigline signifLineColor="red,royalblue" title="$trait on Chr$chr" titleCex=1.3 title2="Top Variant: $LZ_refsnp Covariates: $covariates_LZ" title2Cex=1.0 title2Color="gray30" showRecomb=TRUE recombColor="gray80" recombAxisColor="gray60" --plotonly --cache None --no-date --prefix "LZ_$trait.$name.ImpRgn3TM5"  >> LZdialog.log 2>&1 
	
	##remove line from TOP_SNPS.csv
	sed -i '$d' TOP_SNPS.csv

	##add LZ file name to doc
	awk -F"," 'NR==2 {OFS=","; print $1,"chr"$2":"$3,$7}' sorted_$trait.$suffix.csv > newpval_LD.csv
	beta=`awk -F"," 'NR==2 {print $6}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
	awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' newpval_LD.csv > newpvalLD.csv
	
	lzname="LZ_$trait.$name_$calcLD_refsnp.pdf"
	awk -F"," -v lzname=$lzname '{OFS=","} {$5=lzname; print}' newpvalLD.csv > newpval_withLZ.csv
	
	bpname="$trait.$name.pdf"
	awk -F"," -v bp=$bpname '{OFS=","; $6=bp; print}' newpval_withLZ.csv > newpval_withFiles.csv
	
	awk -F"," 'NR==1 {print; exit}' newpval_withFiles.csv >> TOP_SNPS.csv

	#cut down markers for inclusion in further analysis
	awk -F, -v pcut=$pcut '$7 < $pcut { print $1}' sorted_$trait.$suffix.csv > markers.csv
	
	#add one to variables
	((a+=1))
	
	#format covariates file - switch from row to columns
	sed 's/ \+/,/g' covar.txt > top_row.csv
	tr -s ',' '\n' < top_row.csv | tr -d '"' > top_col.csv
	
	#add new SNP to covariates
	awk -v a=$b -F"," 'NR==a {print $1}' TOP_SNPS.csv >> top_col.csv
	awk 'BEGIN {ORS = ","}{print}' top_col.csv > covar_new.csv
	sed 's/\,/\t/g' covar_new.csv > cov_new.txt
	sed 's/ /\	/g' cov_new.txt > covar.txt
	
	#make LZ covar file
	sed 's/ \+/,/g' LZ_cov.txt > top_row.csv
	tr -s ',' '\n' < top_row.csv | tr -d '"' > top_col.csv
	LZ_refsnp=`awk -F, -v b=$b 'NR==b {print $2}' TOP_SNPS.csv`
	echo $LZ_refsnp >> col_LZ.csv
	awk 'BEGIN {ORS=","}{ print}' col_LZ.csv > cov_LZ.csv
	sed 's/\,/\ /g' cov_LZ.csv > LZ_cov.txt

	#set new pval for loop
	pval_dec=`awk -F, -v b=$b 'NR==b {print $3}' TOP_SNPS.csv`
	pval=`printf "%.0g" "$pval_dec"`

	r=$(awk 'BEGIN {print ('${pval^^}'<'${max_pval^^}')?1:0}')
	
	if [ $a -gt 9 ]
	then	
		break
	fi
done

#make final boxplot
paste col7*.txt > geno.raw

echo "rs123456789" > sums.raw

a=1
wc snp_genotypes.raw > wc.txt
line_count=`awk '{print $1}' wc.txt`

while [ $a -lt $line_count ] 
do
	awk -v a=$a '{if (NR==a) print $0}' geno.raw > line.raw
	#transpose
	awk '
	{
		for (i = 1; i <= NF; i++) {
			if(NR == 1) {
				s[i] = $i;
			} else {
				s[i] = s[i] " " $i;
			}
		}
	}
	END {
		for (i = 1; s[i] != ""; i++) {
			print s[i];
		}
	}' line.raw > snp.trans.raw
	
	#sum column
	awk '{s+=$1}END{print s}' snp.trans.raw >> sums.raw
	
	((a+=1))
	
done

paste snp_genotypes.raw sums.raw > snp.sums.raw

cat snp.sums.raw | tr "\\t" " " > final.snp.raw

#print column with genotypes
awk '{print $8}' final.snp.raw | tail -n +2 > col8.raw

echo "Allele_count,Subject_count" > final_count.csv

#get sample sizes
for all_num in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do	
	sub=`grep $all_num col8.raw | wc | awk '{print $1}'`
	
	echo $all_num,$sub > count.csv
	
	sub_num=`awk -F"," '{print $2}' count.csv`
	
	if [ $sub_num -gt 0 ]
	then
		awk -F"," 'NR==1 {print; exit}' count.csv >> final_count.csv
	fi
done

#add sample sizes for final boxplot to TOP_SNPS file
echo "  " > break.csv
echo "Statistics for final boxplot" > text.csv
cp TOP_SNPS.csv file1.csv
cat file1.csv break.csv > file2.csv
cat file2.csv text.csv > file3.csv
cat file3.csv final_count.csv > TOP_SNPS.csv

#create boxplot
echo "geno1 = read.table('"final.snp.raw"',header=T, fill = TRUE)
geno2 = geno1[order(geno1$IID),]
pheno1 = read.table('"pheno.txt"',header=T, fill = TRUE)
		
geno_id=geno2[,2]
pheno_id=pheno1[,1]
		
common_subjects = intersect(geno_id,pheno_id)
		
geno = geno2[geno_id %in% common_subjects,]
pheno = pheno1[pheno_id %in% common_subjects,]

pdf('"$trait.final.pdf"')
boxplot(pheno$"$trait"~geno$"rs123456789",ylab='"$trait"',xlab='"Allele Count"',col=c('darkolivegreen3','plum3','royalblue1','goldenrod1','indianred2','cadetblue1','slateblue1','gold','hotpink','thistle','steelblue2'))
stripchart(pheno$"$trait"~geno$"rs123456789", vertical=TRUE, add=TRUE, col='"gray39"', pch=1)
status=dev.off()" > test.final.R

Rscript test.final.R 

#organize/make directory for results
mkdir -p $dir
mv *pdf Results/
mv TOP_SNPS.csv Results/

#cleaning
if [ $debug == 0 ]
then 
	rm *epacts *header *ld *log *nosex *map *bin *ped *.o* *raw *.R *.linear *.snplist wc.txt *test*.marker.list.txt cov*.txt alleles*.txt allrefs.txt epacts*.txt markers.csv cov*.csv pval*.csv sorted*.csv LDL.run*.csv top*.csv run*.csv new*.csv snp*.csv tail* *.removed.subjects* LD* pheno.txt head.LD.csv *pedigree.err.txt althomo* col7* final_count.csv file*.csv break.csv count.csv text.csv p_valLZ.csv test*.csv LZ_cov.txt col_LZ.csv
fi

