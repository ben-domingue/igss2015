

awk '{if ($7<1e-5) print $0}' < ~/igss/SNP_gwas_mc_merge_nogc.tbl.uniq > ~/igss/top_hits.txt #bash, gets only hits less than 1e-5 for a top-hits score

#get rid of ambig snps
read.table("~/igss/top_hits.txt")->x
dim(x)
x[!(x[,2]=="A" & x[,3]=="T"),]->x
x[!(x[,2]=="T" & x[,3]=="A"),]->x
x[!(x[,2]=="C" & x[,3]=="G"),]->x
x[!(x[,2]=="G" & x[,3]=="C"),]->x
dim(x)
#so we're ok on ambig snps.
x[,c(1,2,5)]->x
write.table(x,file="top_hits_plink.txt",quote=FALSE,row.names=FALSE,col.names=FALSE) #necessary file for plink construction of scores

plink --bfile ~/igss/ceu-qc-small --score ~/igss/top_hits_plink.txt --out bmi-1 #bash

#now check the correlation between score and phenotype
read.table("~/igss/bmi-1.profile")->prof
read.table("~/igss/bmi.phen")->phen
prof[,c(2,6)]->prof
phen[,c(2,3)]->phen
merge(prof,phen,by=1)->tmp
tmp[,-1]->tmp
apply(tmp,2,as.numeric)->tmp
cor(tmp)

