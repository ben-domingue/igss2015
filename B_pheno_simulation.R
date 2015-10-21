## mkdir ~/igss/varying/

pheno_sim<-function(##this is a function to simulate phenotypes based on gwas results
                    ##(where noise is added to the gwas effects to create "true" effects)
                    plink.file,
                    gwas.file,
                    snp.list, #these are just the snps that you can authoritatively match between plink files & gwas results
                    S, #controls measurement error for "true" gwas effects
                    h2, #controls h2 of trait
                    out.name) {
    paste(plink.file,".bim",sep="")->bim.file
    scan(snp.list,what="character",quiet=TRUE)->snp.list
    ##
    read.table(gwas.file,header=TRUE)->gwas
    gwas[gwas[,1] %in% snp.list,]->gwas
    read.table(bim.file)->ceu
    intersect(gwas[,1],ceu[,2])->snps #we're going to work with just those snps in common between the gwas results & the ceu genetic data.
    gwas[gwas[,1] %in% snps,]->gwas
    gwas[,c(1,5)]->gwas
    sd(gwas[,2])->s #the SD of the gwas effects will provide the basis for our 
    gwas[,2]+rnorm(nrow(gwas),mean=0,sd=S*s)->gwas[,3] #the new gwas[,3] column will be the "true" (and normally unobserved) true snp effects
    write.table(gwas[,c(1,3)],file=paste("/tmp/ceu_gwas_true_effects.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
    system(paste("gcta64 --bfile ",plink.file," --simu-qt --simu-causal-loci /tmp/ceu_gwas_true_effects.txt --simu-hsq ",h2," --out ",out.name,sep=""))
    TRUE
}

#############################################################################################################################
#get snps that are common
#http://www.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz
system("awk '{print $1, $2, $3, $7, $5}' ~/igss/SNP_gwas_mc_merge_nogc.tbl.uniq > ~/igss/bmi-gwas.result")
#remove non-unique lines
read.table("~/igss/bmi-gwas.result",header=TRUE)->x
x[!duplicated(x[,1]),]->x
write.table(x,file="~/igss/bmi-gwas.result",quote=FALSE,row.names=FALSE,col.names=FALSE)
#first get unambiguous set of markers between gwas & plink files
setwd("~/igss")
plink.file<-"~/igss/ceu-qc-small"
paste(plink.file,".bim",sep="")->bim.file
read.table(bim.file)->bim 
bim[,c(2,5,6)]->bim
names(bim)<-c("MarkerName","Effect_Allele","Other_Allele")
read.table("~/igss/bmi-gwas.result" ,header=TRUE)->eff 
eff[!duplicated(eff[,1]),]->eff
names(eff)<-c("MarkerName","Effect_Allele","Other_Allele")
merge(bim,eff)->z
write.table(z$MarkerName,file="bmi-snps.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#now we're going to simulate lots of phenotypes to make sure they have reasonable properties
setwd("~/igss/varying")
for (S in c(0,1,2.5,5)) {
    for (h2 in c(0.25,.5,.75,1))
    {
        print(c(S,h2))
        pheno_sim(plink.file="~/igss/ceu-qc-small",gwas.file="~/igss/bmi-gwas.result",
                  snp.list="~/igss/bmi-snps.txt",
                  S=S,h2=h2,
                  out.name=paste("bmi",S,h2,sep="__")
          )
    }
}


setwd("~/igss/varying/")
#library(HeritHelper)
make_pgs(plink.file="~/igss/ceu-qc-small",gwas.file="~/igss/bmi-gwas.result",out.name="~/igss/varying/bmi",clump=FALSE,wd="/tmp/pgs-b2/")

options(stringsAsFactors=FALSE)
setwd("~/igss/varying")
read.table("~/igss/varying/bmi.profile",header=TRUE)->sc
sc[,c(2,6)]->sc
list.files(pattern="*.phen")->lf
coors<-list()
for (fn in lf) {
    read.table(fn,header=FALSE)->z
    z[,2:3]->z
    merge(sc,z,by=1)->tmp
    cor(tmp[,-1])[1,2]->coors[[fn]]
}
matrix(unlist(coors),ncol=1,dimnames=list(names(coors),NULL))->coors
strsplit(rownames(coors),"__")->txt
sapply(txt,"[",2)->S
sapply(txt,"[",3)->h2
gsub(".phen","",h2)->h2
data.frame(S=S,h2=h2,coors=coors[,1])->df
lm(coors~as.numeric(h2)+as.numeric(S),df)
plot(NULL,xlim=c(.25,1),ylim=c(-.1,1),xlab="trait heritability",ylab="r(pgs,trait)",main="text indicates increasing gwas measurement error")
text(df$h2,df$coors,df$S)



ln -s ~/igss/varying/bmi__2.5__0.5.phen ~/igss/bmi.phen

