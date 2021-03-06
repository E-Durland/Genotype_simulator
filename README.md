# Genotype_simulator
A method to simulate genotype composition based on minor allele frequency (MAF) from poolseq data
Pooled DNA sequencing (poolseq) is a powerful tool to sequence a large
number of individuals at a fraction of the cost of traditional
approaches which use individual genotypes. This approach is also useful
when working with organisms that are too small to extract sufficient DNA
from a single individual organism.

The inherent limitation of poolseq is that the output data is
represented as the proportion of alleles at a single locus for all
individuals within the pool. For example, imagine a pool of 12
individuals with the genotypes:

> AA=6, AB=2, BB=4

The pool of alleles is calculated as:

> A = 6(AA)x2 + 2(AB) = 14

and

> B = 4(BB)x2 + 2(AB) = 10

for a final pool of A/B= 14/10 Since ‘B’ is the minor allele in this
case, the minor allele frequency (MAF) is:

> 10/(10+14)= 0.4166

In most cases with natural populations, all we would know from poolseq
is that the 12 individuals in this population have a MAF = 41.6%. In
some cases with controlled crosses, on the other hand, we know something
about how the population was created and can infer genotypes from a pool
of allele frequencies.

Take an example where a pool of offspring were generated with *MAF =
50%*. If you know that pool was created from two parents, you can
simulate what genotypes they **may** have had. The first step is to set
up the parameters:

``` r
n_par <- 2  # number of parents

n_alleles <- 2 #number of alleles (assuming diploid)

t_alleles <- n_par * n_alleles #total alleles in parental pool

MAF <- 50  #the minor allele frequency of the pool
```

Next, we create a simulated pool of alleles for the parents:

``` r
p_MAF <- MAF/100 * t_alleles #the estimated number of 'minor' alleles in the parental pool

p_genpool <- c( rep( "A" , t_alleles - p_MAF) , rep( "B" , p_MAF)) #for simplicity's sake, 'B' is always the minor allele
p_genpool # simulated 'pool' of parental alleles
```

    ## [1] "A" "A" "B" "B"

Now we randomly assign 2 alleles to each parent to simulate their
genotype:

``` r
gen_indx <- c( 1:t_alleles ) #create an 'index' of the alleles for random sampling

#for parent 1:
par1 <- sample( gen_indx, 2, replace = FALSE) #assign 2 random integers from the index
gen_indx <- gen_indx[ -par1 ]  #remove these integers from the index 
par1 <- p_genpool[ par1 ] #replace integers with alleles from the pool

#for parent 2:
par2 <- sample( gen_indx, 2, replace = FALSE) #assign 2 random integers from the index (the last two in this case)
par2 <- p_genpool[ par2 ] #replace integers with alleles from the pool
```

Our two parents have the genotypes:

``` r
par1
```

    ## [1] "B" "A"

``` r
par2
```

    ## [1] "B" "A"

If we simulate a cross from these two parents, we expect to get the
genotypes:

``` r
c1 <- paste( par1[1] , par2[1] , sep = "_" ) #allele 1 from parent 1, allele 1 from parent 2 and so on...
c2 <- paste( par1[1] , par2[2] , sep = "_" )
c3 <- paste( par1[2] , par2[1] , sep = "_" )
c4 <- paste( par1[2] , par2[2] , sep = "_" )
#list all 4 crosses:
crosses <- c(c1, c2, c3, c4)
#switch all "B_A" genotypes to "A_B" for clarity:
for (a in 1:length(crosses)){
  if(crosses[a] == "B_A"){
    crosses[a] <- "A_B"
  }
}

#the genotype composition of the offspring is:
xtabs(~crosses)
```

    ## crosses
    ## A_A A_B B_B 
    ##   1   2   1

With MAF = 50% and 2 parents there are only two parental genotypes: **AA
and BB** or **AB and AB**. The former combination will result in all
heterozygote offspring **(AB)** and the latter will create the familiar
1:2:1 ratio **(AA:AB:BB)**.

### What about more parents?

But let’s say we have a pool of 10 parents (5 male, 5 female) with equal
contribution between them (a so-called ‘factorial cross’). If a given
locus has a **MAF of 20%** we would see:

``` r
n_par <- 10  
t_alleles <- n_par * n_alleles 
MAF <- 20  
p_MAF <- MAF/100 * t_alleles
p_genpool <- c( rep( "A" , t_alleles - p_MAF) , rep( "B" , p_MAF))
p_genpool # simulated 'pool' of parental alleles
```

    ##  [1] "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "A" "B" "B" "B"
    ## [20] "B"

Again, we randomly assign 2 alleles to each parent to simulate their
genotype:

``` r
gen_indx <- c( 1:t_alleles ) 
n_males <- 5
n_females <- 5
males <- c()
females <- c()
for ( p in 1:n_par ){
  temp_par <- sample( gen_indx, 2, replace = FALSE)
  gen_indx<-gen_indx[!(gen_indx %in% temp_par) ] #remove index numbers by identity, not position
  #assign the genotype, this time with "_" between alleles to keep them together:
  temp_par <- paste(p_genpool[ temp_par[1] ], p_genpool[ temp_par[2] ], sep = "_") 
  if(p <= n_males){
    males <- append(males, temp_par)
  }else{
    females <- append(females, temp_par)
  }
}
```

Our two parental groups have the genotypes:

``` r
males
```

    ## [1] "A_A" "A_A" "B_A" "B_A" "B_A"

``` r
females
```

    ## [1] "A_A" "A_B" "A_A" "A_A" "A_A"

If we simulate a factorial cross from these parents, we expect to get
the genotypes:

``` r
crosses <- c()
for(dad in 1:length(males)){
  sire <- males[dad]%>%
    strsplit(.,split="_")%>%
    unlist()
  for(mom in 1:length(females)){
    dam <- females[mom]%>%
      strsplit(.,split="_")%>%
      unlist()
    c1<-paste(sire[1],dam[1],sep="_")
    c2<-paste(sire[1],dam[2],sep="_")
    c3<-paste(sire[2],dam[1],sep="_")
    c4<-paste(sire[2],dam[2],sep="_")
    #add crosses to the list:
    crosses<-append(crosses,c1)
    crosses<-append(crosses,c2)
    crosses<-append(crosses,c3)
    crosses<-append(crosses,c4)
  }
}

#switch all "B_A" genotypes to "A_B" for clarity:
for (a in 1:length(crosses)){
  if(crosses[a] == "B_A"){
    crosses[a] <- "A_B"
  }
}
#the genotype composition of the offspring is:
gtps<-xtabs(~crosses)
gtps
```

    ## crosses
    ## A_A A_B B_B 
    ##  63  34   3

``` r
#and, to check, the MAF of the pool is:
(gtps[[3]]*2+gtps[[2]])/(sum(gtps)*2)
```

    ## [1] 0.2

Of course, in this scenario where the parents are balanced and crossed
factorially, the offspring genotypes merely reflect what we would expect
under Hardy Weinberg Equillibrium (HWE):

> p<sup>2</sup> + 2pq + q<sup>2</sup> =1

which, when **MAF=20%**:

> p<sup>2</sup> = 0.8<sup>2</sup> =0.64

> 2pq = 2 x 0.8 x 0.2 = 0.32

> q<sup>2</sup> = 0.2<sup>2</sup> = 0.04

obviously similar to our simulated genotype composition

### What about when parents are imbalanced?

Not all crosses are balanced, however, and imbalanced cross designs can
skew genotype compositions away from those predicted by HWE. If we know
the MAF of the pool and the cross design, however, we can simulate the
range of possibilities and predict a mean probability. Predicting the
genotype composition of a pool from MAF in these scenarios requires a
five step process:

1.  Simulate parental gamete pool
2.  Randomize parental genotypes
3.  Simulate crosses
4.  Calculate genotype composition of simulated offspring
5.  Calculate MAF of simulated offspring

The first 4 steps are already familiar to us. The fifth step is
important because in imbalanced cases the journey from pool to parents
and back can actually result in an offspring pool that has a different
MAF than you started with. Of course only accurate simulations are
desired so we need this calculation as a guide to filter out the bad
ones.

Let’s simulate a pool created from pairing **5 males** to each of **20
females** (full factorial design) with a **MAF = 30%**

``` r
#let's make a function for simulating parents:
parent_maker<-function(male_n, female_n,MAF){
  n_par <- male_n + female_n
  t_alleles <- n_par *2
  p_MAF <- ceiling(MAF/100 * t_alleles) #must round up so you don't 'lose' gametes
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
```

now to simulate our 5x20 parents:

``` r
parent_maker(5,20,30)
males
```

    ## [1] "A_A" "B_B" "A_B" "A_A" "A_A"

``` r
females
```

    ##  [1] "B_A" "B_B" "B_A" "A_A" "B_A" "A_A" "B_A" "A_B" "A_A" "A_B" "B_A" "B_B"
    ## [13] "B_A" "A_A" "A_A" "A_A" "A_A" "A_A" "A_A" "A_A"

and now simulate our crosses:

``` r
#again, let's make this a function:
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
  #now, get summary stats of the pool:
  #realized MAF:
  r_MAF<-length(genFreq[genFreq=="B"])/length(genFreq)
  #each genotype proportion:
  p_AA<-length(crosses[crosses=="A_A"])/length(crosses)
  p_AB<-length(crosses[crosses=="A_B"|crosses=="B_A"])/length(crosses)
  p_BB<-length(crosses[crosses=="B_B"])/length(crosses)
  return(c("r_MAF"=round(r_MAF*100,2),"p_AA"=p_AA,"p_AB"=p_AB,"p_BB"=p_BB))
}
cross_sim(males,females)
```

    ## r_MAF  p_AA  p_AB  p_BB 
    ## 30.00  0.49  0.42  0.09

The result from the above chunk of code may or may not have given a
simulation with a *realized* MAF of the pool (‘r\_MAF’) = 30%. Since
with imbalanced crosses this can vary, we need to run the simulation
multiple times to get a ‘population’ of realistic results:

``` r
for(a in 1:50){
  parent_maker(5,20,30)
  if(a == 1){
    sims <- cross_sim(males,females)
  }else{
  sims <- rbind(sims,cross_sim(males,females))
  }
}
sims<-as.data.frame(sims)
sims[sims$r_MAF==30,]
```

    ##      r_MAF p_AA p_AB p_BB
    ## X.3     30 0.49 0.42 0.09
    ## X.6     30 0.49 0.42 0.09
    ## X.10    30 0.49 0.42 0.09
    ## X.12    30 0.49 0.42 0.09
    ## X.15    30 0.49 0.42 0.09
    ## X.16    30 0.49 0.42 0.09
    ## X.17    30 0.49 0.42 0.09
    ## X.22    30 0.49 0.42 0.09
    ## X.25    30 0.49 0.42 0.09
    ## X.26    30 0.49 0.42 0.09
    ## X.27    30 0.49 0.42 0.09
    ## X.29    30 0.49 0.42 0.09
    ## X.36    30 0.49 0.42 0.09
    ## X.39    30 0.49 0.42 0.09
    ## X.42    30 0.49 0.42 0.09
    ## X.44    30 0.49 0.42 0.09
    ## X.46    30 0.49 0.42 0.09
    ## X.48    30 0.49 0.42 0.09

In this example, most of the simulations that result in a realized MAF
of 30% have

-   AA = 48.8%
-   AB = 42.5%
-   BB = 8.8%

this is similar to HWE expectations:

-   p<sup>2</sup>= 0.7<sup>2</sup> = 0.49
-   2pq = 2 x 0.7 x 0.3 = 0.42
-   q<sup>2</sup> = 0.3<sup>2</sup> = 0.09

but there may be notable deviations from HWE in the simulation runs. To
see how this varies across a range of MAF’s, let’s use another ‘for’
loop:

``` r
for(a in 1:50){
  for(i in 1:30){
    parent_maker(5,20,a)
    if(a==1){
      sims<-cross_sim(males,females)
    }else{
      sims<-rbind(sims,cross_sim(males,females))
    }
  }
}
sims<-as.data.frame(sims)
```

Now let’s plot how this looks:

``` r
library(reshape2)
library(ggplot2)
m_sims<-melt(sims,id.vars=c("r_MAF"),variable.name = "gtype",value.name = "p")

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
```

![](https://github.com/E-Durland/Genotype_simulator/blob/main/Simulator_output.png)

We can see that even unbalanced factorial crosses *tend* to fall on HWE
expectations but there are numerous simulations that are significantly
different from the ’equillibrium model. The nature of the progeny pool
from unbalanced crosses is that heterozygotes have a tendency to be
overrepresented relative to HWE expectations for MAFs \>10 and minor
homozygotes (BB) may be nearly or completely absent when MAF\< 20.

The above functions allow for flexibility in testing different cross
designs to predict genotype pool of progeny.

### THE END
