options(stringsAsFactors=FALSE)
read.table("~/igss/bmi.phen")->phen
phen[,c(2,3)]->phen
names(phen)[2]<-"bmi"

read.table("~/igss/bmi-1.profile")->prof
prof[,c(2,6)]->prof
names(prof)[2]<-"pgs-1"
merge(phen,prof,by=1)->df

read.table("~/igss/bmi-2.profile")->prof
prof[,c(2,6)]->prof
names(prof)[2]<-"pgs-2"
merge(df,prof,by=1)->df
read.table("~/igss/bmi-2cl.profile")->prof
prof[,c(2,6)]->prof
names(prof)[2]<-"pgs-2cl"
merge(df,prof,by=1)->df

read.table("~/igss/bmi-3_SCORES_AT_ALL_THRESHOLDS.txt",header=TRUE)->sc
merge(df,sc,by=1)->df

for (i in 2:ncol(df)) as.numeric(df[,i])->df[,i]
cor(df[,-1])->C
#correlation with trait
C[,1]
#correlation with top-hits score
C[,2]
#correlation with unclumped homebrew score
C[,3]


#zipping everything up to download
cd ~
zip -r igss.zip ~/igss/
