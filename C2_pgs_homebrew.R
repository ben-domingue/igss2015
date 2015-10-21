##http://www.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz #don't run, this is where the GWAS results originate
awk '{print $1, $2, $3, $7, $5}' ~/igss/SNP_gwas_mc_merge_nogc.tbl.uniq > ~/igss/bmi-gwas.result #bash, get select columns of the gwas file




make_pgs<-function(plink.file=NULL,gwas.file="/tmp/GWAS.result",wd="/tmp/grs/",out.name,clump=TRUE) {
    #plink.file self-explanatory, the plink file that will be used.
    #gwas.file the set of gwas results. they need to be in the right format when you feed them to make_pgs: marker name, allele 1, allele 2, p value, effect size
    #wd the working directory. need to change if using windows.
    #clump: do you want to compute a pgs based on "clumped" snps?
    #################################
    #0. preamble
    if (is.null(plink.file)) stop("need to pass plink file")
    options(stringsAsFactors=FALSE)
    tr<-list()
    getwd()->orig.dir
    system(paste("mkdir ",wd))
    setwd(wd)
    #i'm just creating symbolic links here so that i can work in one directory. not really imperative that you do something like this, but if you don't youll need to change directory stuff downstream.
    system(paste("ln -s ",plink.file,".bed ./gen.bed",sep=""))
    system(paste("ln -s ",plink.file,".bim ./gen.bim",sep=""))
    system(paste("ln -s ",plink.file,".fam ./gen.fam",sep=""))
    #################################
    #1. remove ambiguous snps from gwas files (strand issues). remember that we're only going to use SNPs which have a quickly identifiable strand.
    system(paste("awk '{ print $1, $2 $3, $4, $5}' ",gwas.file," > temp"))
    read.table("temp")->tmp
    toupper(tmp[,2])->tmp[,2]
    write.table(tmp,file="temp",quote=FALSE,row.names=FALSE,col.names=FALSE)
    #the below "ff" strings are inserted just to make life at the bash command line easier as otherwise one has too many levels of quotes to quickly parse.
    "awk '{ if ($2 == ffACff || $2 == ffAGff || $2 == ffCAff || $2 == ffCTff || $2 == ffGAff || $2== ffGTff || $2 == ffTCff || $2 == ffTGff ) print $0}' temp > GWAS.noambig"->txt
    gsub('ff','"',txt)->txt
    system(txt)
    system("awk '{print $1}' GWAS.noambig > GWAS.snps")
    #################################
    #2. clumping!
    if (clump) { #only gets run when you are clumping
        #this just adds a header to a file.
        system('echo "SNP Allele1Allele2 P W" > head.txt')
        system("cat head.txt GWAS.noambig > GWAS2.noambig")
        #Clump data in 2 rounds using plink2
        #1st clumping & extract tops snps for 2nd round
        fun<-function(i) {
            paste("plink --bfile gen --chr ",i,"  --clump GWAS2.noambig  --clump-p1 1 --clump-p2 1 --clump-r2 .5 --clump-kb 250 --out traitX",i,".round1 --silent",sep="")->cmd
            system(cmd)
            system(paste("awk '{print $3, $5}' traitX",i,".round1.clumped > traitX",i,".round2.input",sep=""))
            system(paste("awk '{print $3}' traitX",i,".round1.clumped > traitX",i,".extract2",sep=""))
            cmd
        }
        library(parallel)
        detectCores()->nw
        makeCluster(nw)->cl 
        clusterApply(cl,1:22,fun)->garbage
        garbage[[1]]->tr$clump1
        #2nd clumping & extract tops snps for profile
        fun<-function(i) {
            paste("plink  --bfile gen --chr ",i," --extract traitX",i,".extract2 --clump traitX",i,".round2.input --clump-p1 1 --clump-p2 1 --clump-r2 .2 --clump-kb 5000 --out traitX",i,".round2 --silent",sep="")->cmd
            system(cmd)
            system(paste("awk '{print $3}' traitX",i,".round2.clumped > traitX",i,".selected",sep=""))
            cmd
        }
        clusterApply(cl,1:22,fun)->garbage
        garbage[[1]]->tr$clump2
        stopCluster(cl)
        system("cat traitX1.selected traitX2.selected traitX3.selected traitX4.selected traitX5.selected traitX6.selected traitX7.selected traitX8.selected traitX9.selected traitX10.selected traitX11.selected traitX12.selected traitX13.selected traitX14.selected traitX15.selected traitX16.selected traitX17.selected traitX18.selected traitX19.selected traitX20.selected traitX21.selected traitX22.selected > traitX.selected")
    } else { #gets run if you don't clump. just making a traitX.selected file available for later work.
        system("echo 'SNP' > /tmp/head.txt")
        system("cat /tmp/head.txt GWAS.snps > traitX.selected")
    }
    #################################
    #3. Cleaning plink files
    #this is going to prune the gwas files down to just the (1) clumped & (2) non-ambig snps (the latter are those from the GWAS.snps file)
    read.table("traitX.selected",header=TRUE)->selected
    read.table(gwas.file,header=TRUE)->effects
    #this is just in case any snps are duplicated in gwas files.
    effects[!duplicated(effects[,1]),]->effects
    #first the clumped SNPs
    effects[,1] %in% selected$SNP -> index
    effects[index,]->effects
    effects[,2]<-toupper(effects[,2])
    effects[,3]<-toupper(effects[,3])
    #now the hard part: ambig snps
    #make sure strands are aligned
    bim <- read.table(file="gen.bim",header=FALSE)
    names(bim)<-c("chr","snp","a","b","a1.bim","a2.bim")
    names(effects)<-c("snp","a1.eff","a2.eff","pv","beta")
    NULL->bim$a->bim$b
    merge(bim,effects,by="snp")->test
    #table(test$a1.bim,test$a1.eff) #useful if you want to see what is going on.
    #get rid of ambig strands from bim (already yanked from gwas data)
    test$a1.bim=="T" & test$a2.bim=="A" -> i1
    test$a1.bim=="A" & test$a2.bim=="T" -> i2
    test[!(i1 | i2),]->test
    test$a1.bim=="C" & test$a2.bim=="G" -> i1
    test$a1.bim=="G" & test$a2.bim=="C" -> i2
    test[!(i1 | i2),]->test
    #now make everything a/c. this is so that we can make sense of the risk allele
    for (nm in c("a1.bim","a2.bim","a1.eff","a2.eff")) {
        ifelse(test[[nm]]=="T","A",test[[nm]])->test[[nm]]
        ifelse(test[[nm]]=="G","C",test[[nm]])->test[[nm]]
    }
    ifelse(test$a1.bim!=test$a1.eff,-1*test$beta,test$beta)->test$beta #flipping
    test[,c("snp","a1.eff","beta")]->z
    #now read the bim file back in and overwrite the old allele names
    read.table("gen.bim")->bim 
    bim[,c(2,5)]->bim
    names(bim)<-c("snp","a1.eff")
    NULL->z$a1.eff
    merge(z,bim)->z
    z[,c("snp","a1.eff","beta")]->z
    #
    write.table(z,file="score_file.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
    nrow(z)->tr$final.n
    #################################
    #4. create score!
    setwd(orig.dir)
    system(paste("plink --bfile ",wd,"gen --score ",wd,"score_file.txt --silent --out ",out.name,sep=""))
    dump("tr",file=paste(out.name,".metadata",sep=""))
    #system(paste("rm -r ",wd))
    tr
}


setwd("~/igss/")
#first score without clumping
make_pgs(plink.file="~/igss/ceu-qc-small",gwas.file="~/igss/bmi-gwas.result",out.name="bmi-2",clump=FALSE,wd="/tmp/pgs-c2/")

#now compute a clumped score
#to do this, it will help to focus the gwas results on just those markers in our .bim file
read.table("~/igss/bmi-gwas.result",header=TRUE)->bmi.gwas
read.table("~/igss/ceu-qc-small.bim",header=FALSE)->bim
#bmi.gwas[bmi.gwas[,1] %in% bim[,2],]->bmi.gwas
write.table(bmi.gwas,file="bmi-gwas-clump.result",quote=FALSE,row.names=FALSE)

setwd("~/igss/")
make_pgs(plink.file="~/igss/ceu-qc-small",gwas.file="~/igss/bmi-gwas-clump.result",out.name="bmi-2cl",clump=TRUE,wd="/tmp/pgs-c2-cl/")

read.table("~/igss/bmi.phen")->phen
phen[,c(2,3)]->phen
read.table("~/igss/bmi-2.profile")->prof
prof[,c(2,6)]->prof1
read.table("~/igss/bmi-2cl.profile")->prof
prof[,c(2,6)]->prof2
merge(phen,prof1,by=1)->tmp
merge(tmp,prof2,by=1)->tmp
tmp[,-1]->tmp
apply(tmp,2,as.numeric)->tmp
cor(tmp)


###########################################################
#For use of make_pgs later.

## #options(repos="http://cran.stat.ucla.edu/")
## #install.packages("devtools")
## #library(devtools)
## #install_github("ben-domingue/HeritHelper")
## #library(HeritHelper)
## #make_pgs(plink.file="~/igss/ceu-qc-small",gwas.file="~/igss/bmi-gwas.result",out.name="bmi-2",clump=FALSE,wd="/tmp/pgs-c2/")

#OR
## wget https://github.com/ben-domingue/HeritHelper/archive/master.zip
## unzip master.zip
## R CMD build HeritHelper-master 
## R CMD INSTALL HeritHelper_0.0.1.tar.gz



