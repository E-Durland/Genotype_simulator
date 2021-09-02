#first, the function to make random parents from an MAF pool
parent_maker<-function(male_n, female_n, MAF){
  n_par <- male_n + female_n
  t_alleles <- n_par *2
  p_MAF <- ceiling(MAF/100 * t_alleles)
  p_genpool <- c( rep( "A" , t_alleles - p_MAF) , rep( "B" , p_MAF))
  gen_indx <- c( 1:t_alleles )
  male<-c()
  female<-c()
  for ( p in 1:n_par ){
    temp_par <- sample( gen_indx, 2, replace = FALSE)
    gen_indx<-gen_indx[!(gen_indx %in% temp_par) ] #remove index numbers by identity, not position
    #assign the genotype, this time with "_" between alleles to keep them together:
    temp_par <- paste(p_genpool[ temp_par[1] ], p_genpool[ temp_par[2] ], sep = "_") 
    if(p <= male_n){
      male <- append(male, temp_par)
    }else{
      female <- append(female, temp_par)
    }
  }
  males <<- male
  females <<- female
}
#for a parental genotype pool with 2 males, 10 females and a MAF = 30% for this locus:
parent_maker(2,10,30)

#now a function to make crosses from this parental pool:
cross_sim<-function(dads,moms){
  genFreq<-c()
  crosses<-c()
  for(dad in 1:length(dads)){
    sire<-unlist(strsplit(dads[dad],split="_"))
    for(mom in 1:length(moms)){
      dam<-unlist(strsplit(moms[mom],split="_"))
      c1<-paste(sire[1],dam[1],sep="_")
      c2<-paste(sire[1],dam[2],sep="_")
      c3<-paste(sire[2],dam[1],sep="_")
      c4<-paste(sire[2],dam[2],sep="_")
      #add crosses to the list:
      crosses<-append(crosses,c1)
      crosses<-append(crosses,c2)
      crosses<-append(crosses,c3)
      crosses<-append(crosses,c4)
      #add all parental alleles to the list:
      genFreq<-append(genFreq,sire[1])
      genFreq<-append(genFreq,sire[2])
      genFreq<-append(genFreq,dam[1])
      genFreq<-append(genFreq,dam[2])
    }
  }
  #get summary stats of the pool:
  #realized MAF:
  r_MAF<-length(genFreq[genFreq=="B"])/length(genFreq)
  #each genotype proportion:
  p_AA<-length(crosses[crosses=="A_A"])/length(crosses)
  p_AB<-length(crosses[crosses=="A_B"|crosses=="B_A"])/length(crosses)
  p_BB<-length(crosses[crosses=="B_B"])/length(crosses)
  return(c("r_MAF"=round(r_MAF*100,2),"p_AA"=p_AA,"p_AB"=p_AB,"p_BB"=p_BB))
}
#to run a single cross:
sims<-cross_sim(males,females)

#to run multiple (30) simulations:
for(i in 1:30){
   parent_maker(5,20,30)
   if(i==1){
     sims<-cross_sim(males,females)
   }else{
     sims<-rbind(sims,cross_sim(males,females))
   }
}

#for multiple simulations(30) across a range of MAF values (1-50)
for(a in 1:50){
  MAF<-a
  for(i in 1:30){
    parent_maker(5,20,a)
    if(a==1){
      sims<-cross_sim(males,females)
    }else{
      sims<-rbind(sims,cross_sim(males,females))
    }
  }
}

#To plot your simulations:
library(tidyr)
library(ggplot2)
m_sims<-pivot_longer(sims,
                     cols = p_AA:p_BB,
                     names_to ="gtype",
                     values_to =  "p")%>%
  as.data.frame()
ggplot(m_sims,aes(r_MAF,p*100))+
  geom_point(aes(color=gtype),size=0.5)+
  stat_smooth(aes(color=gtype))+
  scale_x_continuous("Minor allele frequency (%)", 
                     limits = c(0,50),
                     breaks=seq(0,50,by=10))+
  scale_y_continuous("Genotypes in pool (%)")+
  scale_color_manual("Genotype", breaks=c("p_AA","p_AB","p_BB"),
                     labels=c("AA","AB","BB"),values=c("red3","green3","dodgerblue"))+
  theme_minimal()
