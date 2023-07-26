

library(targets)

###############
# Run analysis
tar_make() 
###############


# Visualize analysis
tar_visnetwork() # visualizes analysis pipeline


# Examine object from analysis
tar_read(printsByPopulation, branches=1)
tar_read(sampleSize_REEX_curve)
tar_read(coverage_REEX_curve)


# Diagnostics
View(tar_meta()) # useful tool for diagnostics
View(tar_meta(targets_only = TRUE)) # simplified





# Additional checks -------------------------------------------------------

# Data summary statistics: Low-frequency dataset
l <- tar_read(data_low)
l[,.N,] # 1817 whistles
l[,length(unique(finalcategory)),] # 281 categories
l[,numWhistlesPerCat:=.N,by=finalcategory]
l[numWhistlesPerCat==1,.N,] # 155 singleton categories
l2 <- l[numWhistlesPerCat>1,,]
l2[,.N,] # 1662 whistles
l2[,length(unique(finalcategory)),] # 126 categories
range(l2[,unique(numWhistlesPerCat),by='finalcategory']$V1) # of 2plus, 2 to 95 whistles per category
median(l2[,unique(numWhistlesPerCat),by='finalcategory']$V1) # of 2plus, median 5 whistles per category
l2[,numPopsPerCat:=length(unique(population)),by='finalcategory']
l2[numPopsPerCat==1,length(unique(finalcategory)),] # 27 produced by a single population
range(l2[numPopsPerCat>1,unique(numPopsPerCat),by=finalcategory]$V1) # 2 to 11 populations
median(l2[numPopsPerCat>1,unique(numPopsPerCat),by=finalcategory]$V1) # 27 produced by a single population

# Data summary statistics: High-frequency dataset
h <- tar_read(data_high) # 1921 whistles
h[,length(unique(population)),] # 14 populations





 