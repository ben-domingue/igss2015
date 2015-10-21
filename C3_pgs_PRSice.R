cd ~/igss

#first do some fine-tuning just to get the data all as needed.
#.assoc: snp, a1, p beta
awk '{print $1, $2, $7, $5}' ~/igss/SNP_gwas_mc_merge_nogc.tbl.uniq > ~/igss/bmi-gwas-prsice.assoc #bash
sed 's/p/P/' bmi-gwas-prsice.assoc | sed 's/b/BETA/' > bmi-gwas-prsice.assoc2 #bash
#.phen file
awk '{print $2, $3}' ~/igss/bmi.phen > bmi-prsice.phen #bash


#next block in bash
R -q --file=/home/ubuntu/apps/PRSice_v1.23/PRSice_v1.23.R --args \
base bmi-gwas-prsice.assoc2 \
target ceu-qc-small \
slower 0 \
supper 1 \
sinc 0.1 \
covary F \
clump.snps F \
plink /usr/local/bin/plink \
figname bmi-3 \
binary.target F 

#compare in R
read.table("~/igss/bmi-3_SCORES_AT_ALL_THRESHOLDS.txt",header=TRUE)->sc
read.table("~/igss/bmi.phen")->phen
phen[,c(2,3)]->phen
merge(phen,sc,by=1)->tmp
cor(tmp[,-1])[,1]





