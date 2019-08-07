########################################################################
## Results can be found in file TOP_SNPS.csv. LocusZoom plots and     ##
## boxplots can be found in pdf files listed in TOP_SNPS.csv. Please  ##
## begin by inputting data below.   				      			  ##
##########################input data here###############################

#Select suffix, trait, chromosome, minimum and maximum region
suffix=run1
trait=LDL
chr=19
min=8000000 
max=18000000

#Set location for covariate data
covfile=/path/to/covariate/data

#Set location for phenotype data
pheno=/path/to/phenotype/data

#Set location for genomic data
bfile=/path/to/geno/data

#Specify bim file location
bim=/path/to/bim

#Set location for PLINK - download at https://www.cog-genomics.org/plink2 
myPLINK=/path/to/plink

#Set location for LocusZoom - download at https://github.com/statgen/locuszoom-standalone
locuszoom=/path/to/locuszoom

#Set location for perl - download at https://www.perl.org/
perl=/path/to/locuszoom

#Set up covariates: separate covariates with a "," 
covariates_list="AGE, SEX"

#Set minor allele frequency
MAF=0.01

#Set p-value for creating first markers list for LD calculation
pmark=0.01

#Set p-value for ending analysis (exponential form)
max_pvala=1E-6

#Set p-value for creating markers list in loop for LD calculations
pmark_loop=0.05

#Set p-value for cutting down SNP marker set after each MMAP analysis
pcut=0.10

#Set max value for y axis in LocusZoom plot
ymax=15

#Set significance lines for LocusZoom plot. Format "<upper_line>,<lower_line>"
sigline="6,8"

#Set name for results directory
dir=Results

#Set memory amount
mem=5000

#Is the variant name in chr:pos:ref:alt format (0=yes; 1=no)?
format=1

#Do you want the output other than results deleted? (0=yes; 1=no)
debug=1

#Do you want the top SNP to be automatically selected as the refSNP? (0=yes; 1=no; if selecting no, top SNPs will be output for you to select the refSNP for analysis)
select_SNP=0

##############################################beginning of script#####################################################

#check for duplicates
$myPLINK --bfile $bfile --write-snplist --memory $mem --out all_snps
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist

#create covariate file
echo $covariates_list > cov_col.csv
sed 's/\,/\t/g' cov_col.csv > covar.txt
covariates=`awk -F"\t" 'NR==1''{print; exit}' covar.txt` 

$myPLINK --bfile $bfile --linear --pheno $pheno --covar $covfile --covar-name $covariates --maf $MAF --chr $chr --from-bp $min --to-bp $max --exclude duplicated_snps.snplist --allow-no-sex --keep-allele-order --memory $mem --out $trait.$suffix 

#run plink just without duplicated snps
$myPLINK --bfile $bfile --exclude duplicated_snps.snplist --keep-allele-order --make-bed --out nodup.$trait.$suffix

#make file to rename from rs to chr:pos (needed for LZ)
awk '{OFS=" "; print $2,"chr"$1":"$4}' nodup.$trait.$suffix.bim >> rename.txt

$myPLINK --bfile nodup.$trait.$suffix --update-name rename.txt 2 1 --keep-allele-order --make-bed --out update.name

$myPLINK --bfile update.name --linear --pheno $pheno --covar $covfile --covar-name $covariates --maf $MAF --chr $chr --from-bp $min --to-bp $max  --allow-no-sex --keep-allele-order --memory $mem --out $trait.$suffix.rename 

#convert output to csv
sed -r 's/\s+/,/g' $trait.$suffix.assoc.linear | sed -r 's/^,//g' > $suffix.nosort.csv

#take first header
awk -F"," 'NR==1 {print}' $suffix.nosort.csv > $suffix.csv

#keep add model only
grep "ADD" $suffix.nosort.csv >> $suffix.csv

##sort the add.mle.pval.csv file by pvalue	
awk -F"," '$9 != "NA"' $suffix.csv > $trait.$suffix.csv
(head -n 1 $trait.$suffix.csv && tail -n +2 $trait.$suffix.csv | sort -gt, -k 9,9) > rs_$trait.$suffix.csv

awk -F"," '{OFS=","; print $1,"chr"$1":"$3,$3,$4,$5,$6,$7,$8,$9}' rs_$trait.$suffix.csv > sorted_$trait.$suffix.csv

rm nodup*

## print top hit SNPNAME - LZ
if [ $select_SNP -eq 0 ]
then
	awk -F"," 'NR==2 {OFS=","; print $2}' rs_$trait.$suffix.csv > pval_chrLZ.csv
else
	list=`awk -F, '{print $2,$9,$7}' rs_$trait.$suffix.csv | tail -n +2 | head -n 11 | cat --number`
	echo -e "Enter the line number of the desired variant to select as the reference variant (columns are SNP, p-value, beta): \n $list"
	read line_num
	awk -F"," -v line=$line_num 'NR==line {OFS=","; print $2}' rs_$trait.$suffix.csv > pval_chrLZ.csv
fi
####add pvalue and beta - LZ - can be skipped if not in 20:10331:A:G format
if [ $select_SNP -eq 0 ]
then
	cpos=`awk -F"," 'NR==2 {OFS=","; print $2}' sorted_$trait.$suffix.csv`
	awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' pval_chrLZ.csv > p_valLZ2.csv
	pv=`awk -F"," 'NR==2 {OFS=","; print $9}' sorted_$trait.$suffix.csv`
	awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' p_valLZ2.csv > p_valLZ.csv
	beta=`awk -F"," 'NR==2 {print $7}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
	awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' p_valLZ.csv > pvalLZ.csv
else
	cpos=`awk -F"," -v line=$line_num 'NR==line {OFS=","; print $2}' sorted_$trait.$suffix.csv > p_valueLZ.csv`
	awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' pval_chrLZ.csv > p_valLZ2.csv
	pv=`awk -F"," -v line=$line_num 'NR==line {OFS=","; print $9}' sorted_$trait.$suffix.csv`
	awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' p_valLZ2.csv > p_valLZ.csv
	beta=`awk -F"," -v line=$line_num 'NR==line {print $7}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
	awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' p_valLZ.csv > pvalLZ.csv
fi

#add header to pval.csv - LZ
echo "SNPNAME, CHR:POS, PVALUE, BETA" > pval.header
cat pval.header pvalLZ.csv > TOP_SNPS.csv	

###create markers.list.text
awk -F, -v pmark=$pmark '$9 < pmark { print $1}' sorted_$trait.$suffix.csv > test1.markera.list.txt
(awk -F, '4 != "<"' test1.markera.list.txt) > test1.marker.list.txt

#create boxplot
name=test1
calcLD_refsnp=`awk -F, 'NR==2 {print $1}' TOP_SNPS.csv`

#create genotype file
$myPLINK --bfile $bfile --out snp_genotypes --recode A --snp $calcLD_refsnp

#get ref/alt alleles
grep "$calcLD_refsnp" $bim > alleles1.txt
	
pos=`awk -F, 'NR==2 {print $3}' rs_$trait.$suffix.csv`
	
grep "$pos" alleles1.txt > alleles.txt

ref=`awk '{print $6}' alleles.txt`
alt=`awk '{print $5}' alleles.txt`

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

if [ -s althomo.txt ]
then
	#make boxplot
	echo "geno1 = read.table('"snp_genotypes.raw"',header=T, fill = TRUE)
	geno2 = geno1[order(geno1$"IID"),]
	pheno1 = read.table('"$pheno"',header=T, fill = TRUE)
		
	geno_id=geno2[,2]
	pheno_id=pheno1[,1]
		
	common_subjects = intersect(geno_id,pheno_id)
		
	geno = geno2[geno_id %in% common_subjects,]
	pheno = pheno1[pheno_id %in% common_subjects,]
		
	colnames(geno)[7]='"genotype"'
	
	pdf('"$trait.$name.pdf"')
	par(mgp=c(3,1.5,0))
	boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects','"$alt_homo" \n "$subc" subjects'))
	status=dev.off()" > test.R
else 
	echo "geno1 = read.table('"snp_genotypes.raw"',header=T, fill = TRUE)
	geno2 = geno1[order(geno1$"IID"),]
	pheno1 = read.table('"$pheno"',header=T, fill = TRUE)
		
	geno_id=geno2[,2]
	pheno_id=pheno1[,1]
		
	common_subjects = intersect(geno_id,pheno_id)
		
	geno = geno2[geno_id %in% common_subjects,]
	pheno = pheno1[pheno_id %in% common_subjects,]
		
	colnames(geno)[7]='"genotype"'
	
	pdf('"$trait.$name.pdf"')
	par(mgp=c(3,1.5,0))
	boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects'))
	status=dev.off()" > test.R
fi

Rscript test.R 

echo $calcLD_refsnp > allrefs.txt

############ calculate LD with top hit as ref variant#################
name=test1
listLD=$name.marker.list.txt
calcLD_refsnp=`awk -F, 'NR==2 {print $2}' TOP_SNPS.csv`

# BASIC DEFINITIONS:
datasetLD=$name.subset
resultsLD=$name.results
plinkOut1=$name.plinkOut1.log
skip_fields=8  # number of columns to skip before individual genotypes begin

# Run Plink to perform the LD calculations on everything in the "subset" dataset
$myPLINK --bfile update.name --r2 dprime --ld-snp $calcLD_refsnp --make-founders --maf $MAF --out $resultsLD --keep-allele-order --ld-window 999999999 --ld-window-kb 999999999 --ld-window-r2 0 > $plinkOut1 2>&1

#sort result file
(head -n 1 $name.results.ld && tail -n +2 $name.results.ld | sort -gr -k7,7) > sorted.$name.results.ld

####create epacts & ld file####
f=sorted_$trait.$suffix.csv
awk -F, '{OFS="\t"; print $1,$3,$3,$2,$9}' $f | tail -n +2 > epacts.INPUT.txt
	
#add headers
echo "#CHROM	BEGIN	END	MARKER_ID	PVALUE" > epacts.header
cat epacts.header epacts.INPUT.txt > epac.$name.INPUT_final.txt.epacts
	
####create ld file####
ldfile=sorted.$name.results.ld 
awk '{OFS="\t"; print $3,$6,$8,$7}' $ldfile | tail -n +2 > locuszoom.INPUT.txt.ld
	
#add in headers
echo "snp1	snp2	dprime	rsquare" > ld.header
cat ld.header locuszoom.INPUT.txt.ld > locuszoom.$name.INPUT_final.txt.ld

##make files tab delimited###
sed 's/ /\t/g' epac.$name.INPUT_final.txt.epacts > ts.epac.$name.INPUT_final.txt.epacts
sed 's/ /\t/g' locuszoom.$name.INPUT_final.txt.ld > ts.locuszoom.$name.INPUT_final.txt.ld

#variables for LocusZoom
vertline=`awk -F, 'NR==2 {print $3}' sorted_$trait.$suffix.csv`
LZ_refsnp=`awk -F, 'NR==2 {print $2}' TOP_SNPS.csv`

########################## generate LocusZoom plot #######################################

$locuszoom  --build hg38 --epacts ts.epac.$name.INPUT_final.txt.epacts --ld ts.locuszoom.$name.INPUT_final.txt.ld --refsnp $LZ_refsnp --chr $chr --start $min --end $max geneFontSize=.5 ymax=$ymax vertline="$vertline" vertLineColor="green4" vertLineType="1" vertLineWidth="1" signifLine=$sigline signifLineColor="red,royalblue" title="$trait on Chr$chr" titleCex=1.3 title2="Top Variant: $calcLD_refsnp Covariates: $covariates" title2Cex=1.0 title2Color="gray30" showRecomb=TRUE recombColor="gray80" recombAxisColor="gray60" --plotonly --cache None --no-date --prefix "LZ_$trait.$name"  >> LZdialog.log 2>&1 

##move file name onto doc##
awk -F"," 'NR==2 {OFS=","; print}' TOP_SNPS.csv > pvalLD.csv

lzname="LZ_$trait.$name_$calcLD_refsnp.pdf"
awk -F"," -v lzname=$lzname '{OFS=","} {$5=lzname; print}' pvalLD.csv > pval_withLD.csv
bpname="$trait.$name.pdf"
awk -F"," -v bp=$bpname '{OFS=","; $6=bp; print}' pval_withLD.csv > pval_withFiles.csv

#add header#
echo "SNPNAME, CHR:POS, PVALUE, BETA, LOCUSZOOM, BOXPLOT" > pval_withFiles.header
cat pval_withFiles.header pval_withFiles.csv > TOP_SNPS.csv

#cut down markers for inclusion in further analysis
awk -F, -v pcut=$pcut '$9 < $pcut { print $1}' sorted_$trait.$suffix.csv > markers.csv

rm update.name nodup.$trait.$suffix* $trait.$suffix.rename*


#################################################################################################################################
############################################### start of loop ###################################################################
#################################################################################################################################

a=2

#switch row to column
sed 's/ \+/,/g' covar.txt > top_row.csv
tr -s ',' '\n' < top_row.csv | tr -d '"' > top_col.csv

#add to cov list
calcLD_refsnp=`awk -F, 'NR==2 {print $1}' TOP_SNPS.csv`
alt=`awk '{print $5}' alleles.txt`
new=$calcLD_refsnp"_"$alt
echo $new >> top_col.csv
awk 'BEGIN {ORS = ","}{print}' top_col.csv > covar_new.csv
sed 's/\,/\t/g' covar_new.csv > covar_LZ.txt

#cut first two columns of covar_LZ.txt
cut -f 3- $covfile > cov2.txt

#print genotypes
awk '{print $1,$2,$7}' snp_genotypes.raw > snp2.txt

#combine covar and geno
paste snp2.txt cov2.txt > newcovar.txt

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
	covariates=`awk 'NR==1 {print; exit}' covar_LZ.txt`
	suffix=run$a

	$myPLINK --bfile $bfile --linear --pheno $pheno --covar newcovar.txt --covar-name $covariates --maf $MAF --chr $chr --from-bp $min --to-bp $max --exclude duplicated_snps.snplist --allow-no-sex --keep-allele-order --memory $mem --out $trait.$suffix 
	
	#run plink just without duplicated snps
	$myPLINK --bfile $bfile --exclude duplicated_snps.snplist --keep-allele-order --make-bed --out nodup.$trait.$suffix

	#make file to rename from rs to chr:pos (needed for LZ)
	awk '{OFS=" "; print $2,"chr"$1":"$4}' nodup.$trait.$suffix.bim >> rename.$suffix.txt

	$myPLINK --bfile nodup.$trait.$suffix --update-name rename.txt 2 1 --keep-allele-order --make-bed --out update.name

	$myPLINK --bfile update.name --linear --pheno $pheno --covar newcovar.txt --covar-name $covariates --maf $MAF --chr $chr --from-bp $min --to-bp $max  --allow-no-sex --memory $mem --keep-allele-order --out $trait.$suffix.rename 

	#convert output to csv
	sed -r 's/\s+/,/g' $trait.$suffix.assoc.linear | sed -r 's/^,//g' > $suffix.nosort.csv

	#take first header
	awk -F"," 'NR==1 {print}' $suffix.nosort.csv > $suffix.csv

	#keep add model only
	grep "ADD" $suffix.nosort.csv >> $suffix.csv

	##sort  by pvalue	
	awk -F"," '$9 != "NA" && $4 != "<"' $suffix.csv > $trait.$suffix.csv
	(head -n 1 $trait.$suffix.csv && tail -n +2 $trait.$suffix.csv | sort -gt, -k 9,9) > rs_$trait.$suffix.csv

	awk -F"," '{OFS=","; print $1,"chr"$1":"$3,$3,$4,$5,$6,$7,$8,$9}' rs_$trait.$suffix.csv > sorted_$trait.$suffix.csv

	rm nodup*

	## print top hit SNPNAME - LZ
	awk -F"," 'NR==2 {OFS=","; print $2}' rs_$trait.$suffix.csv > new.csv
		
	####add pvalue and beta - LZ
	cpos=`awk -F"," 'NR==2 {OFS=","; print $2}' sorted_$trait.$suffix.csv`
	awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' new.csv > new2.csv
	pv=`awk -F"," 'NR==2 {OFS=","; print $9}' sorted_$trait.$suffix.csv`
	awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' new2.csv > new3.csv
	beta=`awk -F"," 'NR==2 {print $7}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
	awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' new3.csv > new4.csv
	awk -F"," 'NR=1 {print; exit}' new4.csv >> TOP_SNPS.csv
	
	###create markers.list.text  
	name=test$a
	awk -F, -v pmark_loop=$pmark_loop '$9 < pmark_loop {print $1}' sorted_$trait.$suffix.csv > $name.marker.list.txt
	
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
	
	#create boxplot
	calcLD_refsnp=`awk -F, 'NR==2 {print $2}' rs_$trait.$suffix.csv`
	
	#create genotype file
	$myPLINK --bfile $bfile --out $suffix.snp_genotypes --recode A --snp $calcLD_refsnp
	
	grep "$calcLD_refsnp" $bim > alleles1.$suffix.txt
	
	pos=`awk -F, 'NR==2 {print $3}' rs_$trait.$suffix.csv`
	
	grep "$pos" alleles1.$suffix.txt > alleles.$suffix.txt

	ref=`awk -F"\t" '{print $6}' alleles.$suffix.txt`
	alt=`awk -F"\t" '{print $5}' alleles.$suffix.txt`
	
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
	awk '{print $7}' $suffix.snp_genotypes.raw | tail -n +2 > $suffix.snp.raw
	suba=`grep "0" $suffix.snp.raw | wc | awk '{print $1}'`
	subb=`grep "1" $suffix.snp.raw | wc | awk '{print $1}'`
	subc=`grep "2" $suffix.snp.raw | wc | awk '{print $1}'`
	
	#check for homozygotes for alt allele
	awk '{print $7}' $suffix.snp_genotypes.raw | tail -n +2 > col7.$a.txt

	grep "2" col7.$a.txt > althomo.$name.txt

	if [ -s althomo.$name.txt ]
	then
		#make boxplot
		echo "geno1 = read.table('"$suffix.snp_genotypes.raw"',header=T, fill = TRUE)
		geno2 = geno1[order(geno1$"IID"),]
		pheno1 = read.table('"$pheno"',header=T, fill = TRUE)
		
		geno_id=geno2[,2]
		pheno_id=pheno1[,1]
		
		common_subjects = intersect(geno_id,pheno_id)
		
		geno = geno2[geno_id %in% common_subjects,]
		pheno = pheno1[pheno_id %in% common_subjects,]
		
		colnames(geno)[7]='"genotype"'
	
		pdf('"$trait.$name.pdf"')
		par(mgp=c(3,1.5,0))
		boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects','"$alt_homo" \n "$subc" subjects'))
		status=dev.off()" > test.$name.R
	else 
		echo "geno1 = read.table('"$suffix.snp_genotypes.raw"',header=T, fill = TRUE)
		geno2 = geno1[order(geno1$"IID"),]
		pheno1 = read.table('"$pheno"',header=T, fill = TRUE)
			
		geno_id=geno2[,2]
		pheno_id=pheno1[,1]
		
		common_subjects = intersect(geno_id,pheno_id)
		
		geno = geno2[geno_id %in% common_subjects,]
		pheno = pheno1[pheno_id %in% common_subjects,]
		
		colnames(geno)[7]='"genotype"'
	
		pdf('"$trait.$name.pdf"')
		par(mgp=c(3,1.5,0))
		boxplot(pheno$"$trait"~geno$"genotype",ylab='"$trait"',xlab='"$calcLD_refsnp"',col=c('darkolivegreen3','plum3','royalblue1'), names=c('"$ref_homo" \n "$suba" subjects', '"$hetero" \n "$subb" subjects'))
		status=dev.off()" > test.$name.R
	fi

	Rscript test.$name.R 

	echo $calcLD_refsnp >> allrefs.txt

	##########################calculate LD####################################
	listLD=$name.marker.list.txt
	
	# BASIC DEFINITIONS:
	datasetLD=$name.subset
	resultsLD=$name.results
	plinkOut1=$name.plinkOut1.log
	skip_fields=8  # number of columns to skip before individual genotypes begin
	calcLD_refsnp=`awk -F, 'NR==2 {print $2}' sorted_$trait.$suffix.csv`

	# Run Plink to perform the LD calculations on everything in the "subset" dataset
	$myPLINK --bfile update.name --r2 dprime --ld-snp $calcLD_refsnp --make-founders --maf $MAF --out $resultsLD --keep-allele-order --ld-window 999999999 --ld-window-kb 999999999 --ld-window-r2 0 > $plinkOut1 2>&1

	#sort result file
	(head -n 1 $name.results.ld && tail -n +2 $name.results.ld | sort -gr -k7,7) > sorted.$name.results.ld
	
	#create epacts file
	f=sorted_$trait.$suffix.csv
	awk -F, '{OFS="\t"; print $1,$3,$3,$2,$9}' $f | tail -n +2 > epacts.INPUT.txt
		
	#add headers
	echo "#CHROM	BEGIN	END	MARKER_ID	PVALUE" > epacts.header
	cat epacts.header epacts.INPUT.txt > epac.$name.INPUT_final.txt.epacts
	
	####create ld file####
	ldfile=sorted.$name.results.ld 
	awk '{OFS="\t"; print $3,$6,$8,$7}' $ldfile | tail -n +2 > locuszoom.INPUT.txt.ld
	
	#add in headers
	echo "snp1	snp2	dprime	rsquare" > ld.header
	cat ld.header locuszoom.INPUT.txt.ld > locuszoom.$name.INPUT_final.txt.ld
	
	##make files tab delimited###
	sed 's/ /\t/g' epac.$name.INPUT_final.txt.epacts > ts.epac.$name.INPUT_final.txt.epacts
	sed 's/ /\t/g' locuszoom.$name.INPUT_final.txt.ld > ts.locuszoom.$name.INPUT_final.txt.ld
	
	covariates_LZ=`awk -Ft 'NR==1 {print}' covar_LZ.txt`
	LZ_refsnp=`awk -F, 'NR==2 {print $2}' sorted_$trait.$suffix.csv`
	
	#generate LocusZoom
	$locuszoom  --build hg38 --epacts ts.epac.$name.INPUT_final.txt.epacts --ld ts.locuszoom.$name.INPUT_final.txt.ld --refsnp $LZ_refsnp --chr $chr --start $min --end $max geneFontSize=.5 ymax=$ymax vertline="$vertline" vertLineColor="green4" vertLineType="1" vertLineWidth="1" signifLine=$sigline signifLineColor="red,royalblue" title="$trait on Chr$chr" titleCex=1.3 title2="Top Variant: $calcLD_refsnp Covariates: $covariates_LZ" title2Cex=1.0 title2Color="gray30" showRecomb=TRUE recombColor="gray80" recombAxisColor="gray60" --plotonly --cache None --no-date --prefix "LZ_$trait.$name.ImpRgn3TM5"  >> LZdialog.log 2>&1 

	##remove line from TOP_SNPS.csv
	sed -i '$d' TOP_SNPS.csv

	## print top hit SNPNAME - LZ
	awk -F"," 'NR==2 {OFS=","; print $2}' rs_$trait.$suffix.csv > newp_valueLZ.csv
		
	####add pvalue and beta - LZ
	cpos=`awk -F"," 'NR==2 {OFS=","; print $2}' sorted_$trait.$suffix.csv`
	awk -F"," -v cp=$cpos 'BEGIN {OFS=","} {$2=cp; print}' newp_valueLZ.csv > p_valLZ2.csv
	pv=`awk -F"," 'NR==2 {OFS=","; print $9}' sorted_$trait.$suffix.csv`
	awk -F"," -v pv=$pv 'BEGIN {OFS=","} {$3=pv; print}' p_valLZ2.csv > newp_valLZ.csv
	beta=`awk -F"," 'NR==2 {print $7}' sorted_$trait.$suffix.csv | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'`
	awk -F"," -v b=$beta 'BEGIN {OFS=","} {$4=b; print}' newp_valLZ.csv > newpvalLZ.csv
	lzname="LZ_$trait.$name_$calcLD_refsnp.pdf"
	awk -F"," -v lzname=$lzname '{OFS=","} {$5=lzname; print}' newpvalLZ.csv > newpval_withLZ.csv
	bpname="$trait.$name.pdf"
	awk -F"," -v bp=$bpname '{OFS=","; $6=bp; print}' newpval_withLZ.csv > newpval_withFiles.csv
	awk -F"," 'NR==1 {print; exit}' newpval_withFiles.csv >> TOP_SNPS.csv

	#cut down markers for inclusion in further analysis
	awk -F, -v pcut=$pcut '$9 < $pcut { print $1}' sorted_$trait.$suffix.csv > markers.csv
	
	#add one to variables
	((a+=1))
	
	#format covariates file - switch from row to columns
	sed 's/ \+/,/g' covar_LZ.txt > top_row.csv
	tr -s ',' '\n' < top_row.csv | tr -d '"' > top_col.csv
	
	#add new SNP to covariates
	calcLD_refsnp=`awk -F, 'NR==2 {print $2}' rs_$trait.$suffix.csv`
	alt=`awk '{print $5}' alleles.$suffix.txt`
	new=$calcLD_refsnp"_"$alt
	echo $new >> top_col.csv
	awk 'BEGIN {ORS = ","}{print}' top_col.csv > covar_new.csv
	sed 's/\,/\t/g' covar_new.csv > covar_LZ.txt
	
	#cut first two columns of covar_LZ.txt
	cut -d " " -f 3- newcovar.txt > cov$a.txt

	#print genotypes
	awk '{print $1,$2,$7}' $suffix.snp_genotypes.raw > snp$a.txt

	#combine covar and geno
	paste snp$a.txt cov$a.txt > newcovar.txt

	#set new pval for loop
	pval_dec=`awk -F, -v b=$b 'NR==b {print $2}' TOP_SNPS.csv`
	pval=`printf "%.0g" "$pval_dec"`

	r=$(awk 'BEGIN {print ('${pval^^}'<'${max_pval^^}')?1:0}')
	
	rm update.name nodup.$trait.$suffix* $trait.$suffix.rename*
	
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

echo "Allele_count,Subjects" > final_count.csv

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
echo "Subject numbers for final boxplot" > text.csv
cp TOP_SNPS.csv file1.csv
cat file1.csv break.csv > file2.csv
cat file2.csv text.csv > file3.csv
cat file3.csv final_count.csv > TOP_SNPS.csv

#create boxplot

echo "geno1 = read.table('"final.snp.raw"',header=T, fill = TRUE)
geno2 = geno1[order(geno1$"IID"),]
pheno1 = read.table('"$pheno"',header=T, fill = TRUE)
		
geno_id=geno2[,2]
pheno_id=pheno1[,1]
		
common_subjects = intersect(geno_id,pheno_id)
		
geno = geno2[geno_id %in% common_subjects,]
pheno = pheno1[pheno_id %in% common_subjects,]

pdf('"$trait.final.pdf"')
boxplot(pheno$"$trait"~geno$"rs123456789",ylab='"$trait"',xlab='"Allele Count"',col=c('darkolivegreen3','plum3','royalblue1','goldenrod1','indianred2','cadetblue1','slateblue1','gold','hotpink','thistle','steelblue2'))
status=dev.off()" > test.final.R

Rscript test.final.R 

mv Rplots.pdf $trait.final.pdf

#organize/make directory for results
mkdir -p $dir
mv *pdf Results/
mv TOP_SNPS.csv Results/

#cleaning
if [ $debug == 0 ]
then 
	rm *epacts *header *ld *log *nosex *map *bin *ped *.o* *raw *.R *.linear *.snplist wc.txt *test*.marker.list.txt cov*.txt alleles*.txt allrefs.txt epacts*.txt markers.csv cov*.csv pval*.csv sorted*.csv LDL.run*.csv top*.csv run*.csv new*.csv snp*.csv update* rename* rs* col7*.txt snp*.txt p_val* althomo.txt *csv newcovar.txt althomo*.txt test*.txt
fi