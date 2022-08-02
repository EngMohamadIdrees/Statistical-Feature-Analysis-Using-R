#Reading Data from csv
proteomes <- read.csv("77_cancer_proteomes_CPTAC_itraq.csv")
clinical <- read.csv("clinical_data_breast_cancer.csv")
gene_proteins <- read.csv("PAM50_proteins.csv")

#-------------------------------------Part 1 -------------------------------------------------
#Delete all row have  missing value from data
proteomes=na.omit(proteomes)
#Saving Row names of Proteomes to make it column name
n <- proteomes$RefSeq_accession_number

#Subset Column from 4 index to last index then transpose  the data frame
proteomes <- as.data.frame(t(proteomes[,4:ncol(proteomes)]))

#get back names of row  that saved in n to make it  name of column
colnames(proteomes) <- n

#bind the data frame to proteomes 
proteomes <- cbind(rownames(proteomes), data.frame(proteomes, row.names=NULL))

#name first column same as in clinical data csv
colnames(proteomes)[1] <- "Complete.TCGA.ID"

#convert the ID from AO-A12D.01TCGA to TCGA-AO-A12D as clinical data csv
#defining formula to restructure:
get.clinical.id <- function(proteome.id) {
  x = substr(proteome.id, 4, 7)
  y = substr(proteome.id, 0, 2)
  paste("TCGA",y,x,sep="-")
}

#sapply to id column in proteomes
proteomes$Complete.TCGA.ID <- sapply(proteomes$Complete.TCGA.ID, get.clinical.id)
#inner join
library(dplyr)
data <-  inner_join(clinical, proteomes, by = "Complete.TCGA.ID")

#-------------------------------------part2-----------------------------------------------------------
#2.1
bool=data$HER2.Final.Status
for(i in 1:80)
{
  
  if( bool[i] =='Negative')
  {
    bool[i]=0
    
  }
  else if(bool[i]=="Positive")
  {
    bool[i]=1
  }
  else{
    bool[i]=-1
  }
  
}
bool<-as.numeric(bool)
vector=c()
Col_name=c()
for(i in 1:100)
{
  Gene=gene_proteins[i,2]
  for(j in 31:ncol(data) )
  {
    if(Gene==colnames(data)[j])
    {
      
      final<-cor(data[,j],bool,method ="pearson")
      Col_name<-append(Col_name, Gene)
      vector<-append(vector, final)
    }
  }
}
prot=data.frame(Col_name,vector)
#2.2
prot=prot[order(abs(prot$vector), decreasing = TRUE),]

#2.3
threshold=0.05
filtered=subset(prot,prot$vector<=threshold)
#-------------------------------------part3--------------------------------------
#3
negative=subset(data,data$HER2.Final.Status=='Negative')
Positive=subset(data,data$HER2.Final.Status=='Positive')
New_data<-rbind(negative,Positive)
vector2=c()
Col_name2=c()
for(i in 1:nrow(gene_proteins))
{
  Gene=gene_proteins[i,2]
  for(j in 31:ncol(New_data))
  {
    if(Gene==colnames(New_data)[j])
    {
      mean<-mean(New_data[,j])
      final<-t.test(New_data[,j]~New_data$HER2.Final.Status,data=New_data,mean=mean)[['statistic']]
      Col_name2<-append(Col_name2, Gene)
      vector2<-append(vector2, final)
    }
  }
}
new_data_frame<-data.frame(Col_name2,vector2)
new_data_frame=new_data_frame[order(abs(new_data_frame$vector2), decreasing = TRUE),]

the_rest=subset(new_data_frame,new_data_frame$vector2<=0.05)
the_rest=the_rest[order(abs(the_rest$vector2), decreasing = TRUE),]
#NP_004439   is the most effective  protein 
#both method resulted in different feature sets



