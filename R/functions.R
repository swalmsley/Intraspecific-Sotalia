
# Custom functions to be called in '_targets.R'


# Read data ---------------------------------------------------------------
read_data <- function(path) {
  d <- data.table(read.csv(path))
  return(d)
}
# testing
#testData <- read_data('./input/Sotalia_covariates_final_samprates.csv')


# Process data ------------------------------------------------------------
process_data <- function(data, whistleData, subset){
  
  # clean errors
  data <- data[!is.na(finalcategory),,]
  
  # insert abbreviations for populations
  population_codes <- c('CA', 'CO', 'EC', 'JU', 'CR', 'FG', 'IG', 'PA', 'PE', 'PR', 'RN', 'SC', 'SE', 'SP', 'TO', 'VZ')
  geographical_locations <- c('Central Amazon', 'Colombian Amazon', 'Ecuador', 'Jurua', 'Cost Rica', 'French Guiana', 
                              'Ilha Grande', 'Para', 'Peru', 'Parana', 'Rio Grande do Norte', 'Santa Catarina', 
                              'Sepetiba', 'Sao Paulo', 'Tocantins', 'Venezuela')
  for (i in 1:length(population_codes)) {
    data[Geograhical.location == geographical_locations[i], pop := population_codes[i], ]
  }
  data[,population:=pop,]
  data[,pop:=NULL,]
  
  # merge in acoustic data
  data[,ctrName:=paste(Contour.name, '.ctr',sep=''),]
  #whistleData[,name:=substr(name,start=2,stop=as.integer((str_locate(name,'.ctr')))-1),]
  data <- merge(data,whistleData[,c('name','maxFreq')],by.x='ctrName',by.y='name')
  
  # subset based on low or high-frequency datasets for improved comparisons
  if (subset=='low') (data <- data[maxFreq<=22000,,])
  if (subset=='high') (data <- data[maxFreq <= 48000 & samp.rate.khz>=96,,])
  
  # calculate number of whistles by population and whistle category
  data[, num_type_loc := .N, by = c('population', 'finalcategory')]
  
  return(data)
  
}
# testing


# Build spadeR data -------------------------------------------------------
build_spadeR_data <- function(data) {
  
  pops = unique(data$population)
  num_populations = length(pops)
  categories = unique(data$finalcategory)
  num_categories = length(categories) # will be row number
  
  abundMatrix <- matrix(data=NA, nrow=num_categories, ncol=num_populations) # initialise empty matrix
  colnames(abundMatrix) <- unique(data$population)

    for(i in 1:num_categories){ # loop through each row (whistle types)
    for(j in 1:num_populations){ # loop through each column (whistle category/type)
      abundMatrix[i,j] = as.integer(as.numeric(data[population==pops[j] & finalcategory==categories[i],.N,])) # fill in the blank using the datatable
    }
  }
  return(abundMatrix)
}


# Run permutation test ----------------------------------------------------
run_permutation <- function(data) {
  
  # set up empty table to hold permutation results
  permutation_table <- as.data.table(matrix(0,nrow=1000,ncol=3))
  colnames(permutation_table) <- c('index','propPopSpecific','propSS')
  permutation_table$index <- 1:1000
  
  # add column to data
  data[,numWhistlesByCat:=.N,by='finalcategory']
  
  # loop through 1000 permutations
  for (i in 1:1000){
    
    # extract key bits
    temp <- data[numWhistlesByCat>1,c('finalcategory','population', 'Species')] # subset to whistle types with more than 1 whistle
    
    # shuffle species membership and assign membership classes
    random <- transform(temp,finalcategory=sample(finalcategory,replace=FALSE))
    random[,populationMember:=ifelse(length(unique(population))>1,"Shared","Unique"),by="finalcategory"]
    random[,speciesMember:=ifelse(length(unique(Species))>1,"Shared","Unique"),by="finalcategory"]
    
    # calculate proportions and place them in datatable
    permutation_table$propSS[i] <- (random[speciesMember=="Unique",length(unique(finalcategory)),])/(random[,length(unique(finalcategory)),])
    permutation_table$propPopSpecific[i] <- (random[populationMember=="Unique",length(unique(finalcategory)),])/(random[,length(unique(finalcategory)),])
  }
  
  # determine true proportion of whistle types that are population-specific 
  data[,populationMember:=ifelse(length(unique(population))>1,"Shared","Unique"),by="finalcategory"]
  trueProp <- (data[populationMember=="Unique" & numWhistlesByCat>1,length(unique(finalcategory))])/(data[numWhistlesByCat>1,length(unique(finalcategory)),])
  
  # calculate two-tailed p-value
  null_props <- permutation_table$propPopSpecific
  mean_prop <- mean(permutation_table$propPopSpecific)
  null_diff <- abs(null_props-mean_prop)
  true_diff <- abs(trueProp-mean_prop)
  P_val_2tailed <- (sum(null_diff>true_diff)+1) / (length(null_diff)+1)
  
  return(data.table(pValue=P_val_2tailed, trueProportion=trueProp, meanPermuted=mean_prop))

}


# Run Mantel ----------------------------------------------------------------
run_mantel <- function(geographic_distances, spade_input, spade_output, diversityMeasure) {
  
  # create geographic matrix
  geo_matrix <- as.matrix(geographic_distances[,2:ncol(geographic_distances),]) # convert into pretty matrix
  
  # extract repertoire similarity matrix
  if (diversityMeasure=="horn") (matrix <- spade_output$similarity.matrix$C12)
  if (diversityMeasure=="morisata") (matrix <- spade_output$similarity.matrix$C22)
  colnames(matrix) <- colnames(spade_input)
  rownames(matrix) <- colnames(spade_input)
  
  # convert repertoire matrix to distance matrix
  dist_matrix <- sim2dist(matrix)
  colnames(matrix) <- colnames(spade_input)
  rownames(matrix) <- colnames(spade_input)
  
  # run test
  mtest <- mantel(geo_matrix, dist_matrix, method='spearman')
  
  return(mtest)
  
}


# Run NMDS ----------------------------------------------------------------
run_NMDS <- function(spade_input, spade_output, diversityMeasure) {
  
  if (diversityMeasure=="horn") (matrix <- spade_output$similarity.matrix$C12)
  if (diversityMeasure=="morisata") (matrix <- spade_output$similarity.matrix$C22)
  
  colnames(matrix) <- colnames(spade_input)
  rownames(matrix) <- colnames(spade_input)
  
  dist_matrix <- sim2dist(matrix)
  nmds <- metaMDS(dist_matrix)
  
  return(nmds)
  
}


# Plot NMDS ---------------------------------------------------------------
plot_nmds <- function(spade_input, fitted_nmds, plot_title) {
  
  data.scores <- data.table(scores(fitted_nmds), colnames(spade_input))
  data.scores$grp <- c('SF','SF', 'SG', 'SF', 'SG', 'SG', 'SF', 'SG', 'SF', 'SG', 'SG', 'SG', 'SG', 'SG', 'UNK', 'SG') 
  
  grp.a <- data.scores[data.scores$grp == "SF", ][chull(data.scores[data.scores$grp == "SG", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
  grp.b <- data.scores[data.scores$grp == "SG", ][chull(data.scores[data.scores$grp == "SF", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
  hull.data <- rbind(grp.a, grp.b) 
  
  # Create ggplot
  g <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
    #geom_point(aes(color = grp),size=1,shape=3,alpha=0.2) +
    geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=V2,colour=grp),size=6,vjust=0) +  # add the site labels
    geom_mark_hull(expand=0,radius=0,aes(color=grp,fill=grp),alpha=0.1,linetype='dashed')+
    scale_color_manual(values=c(trait_colors[1],trait_colors[2],trait_colors[3])) +
    scale_fill_manual(name='Species',values=c(trait_colors[1],trait_colors[2],trait_colors[3]),labels=c('Tucuxi', 'Guiana dolphin', 'Unknown')) +
    labs(x='NMDS 1', y='NDMS 2') +
    ggtitle(paste(plot_title)) + 
    guides(color='none', fill=guide_legend(title='Species',override.aes = list(label=''))) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  g
  
}


# Format for iNext --------------------------------------------------------
format_for_iNEXT <- function(data) {
  sorted <- list()
  unique_populations <- unique(data$population)
  catsByPop <- data[data[ , .I[sample(.N,1)] , by = c('finalcategory','population')]$V1]
  for (pop in unique_populations) {
    sorted[[pop]] <- catsByPop[population == pop, num_type_loc][order(-catsByPop[population == pop, num_type_loc])]
  }
  return(sorted)
}


# Run iNEXT ---------------------------------------------------------------
run_iNEXT <- function(data, q_value) {
  q <- iNEXT(data, q=q_value, datatype='abundance')
  return(q)
}


# Extract cMax ------------------------------------------------------------
extract_cMax <- function(q0) {
  
  # identify coverage for each population extrapolated to 2x sample
  output <- data.table(q0$iNextEst$coverage_based$Assemblage, q0$iNextEst$coverage_based$SC)
  colnames(output) <- c('pop', 'sc')
  output[,maxCoverage_by_pop:=max(sc),by=pop]
  
  # calculate cmax as minimum of extrapolated coverages
  cMax <- min(output$maxCoverage_by_pop)
  return(cMax)
  
}
  

# Prepare modelling datasheet ---------------------------------------------
prepModelData <- function(data, group_data, estimates) {
  
  # one row per population
  pop_specific <- data[data[ , .I[sample(.N,1)] , by = 'population']$V1]
  
  # merge with estimates at cMax
  model_dt <- merge(pop_specific, estimates, by.x='population', by.y='Assemblage')
  
  # merge with median group sizes
  model_dt <- merge(model_dt, group_data[,c('population', 'median_group')], by= 'population')
  
  return(model_dt)
  
}


# Plot REEX curves --------------------------------------------
plot_REEX <- function(output, q, curve_type) {

  # curve type 1 = sample-size based, curve type 2 = coverage-based
  ylabs <- c('Richness (q=0)', 'Shannon diversity (q=1)', 'Simpson diversity (q=2)')

  g <- ggiNEXT(output, type=curve_type, se = FALSE) +
    scale_shape_manual(values=rep(20, each = 16)) +
    scale_color_manual(values=species_colors) +
    #scale_fill_manual(values=species_colors) +
    labs(x='Number of whistles analysed', y=ylabs[q+1]) +
    guides(color='none', shape='none', linetype='none') +
    scale_linewidth_manual(values = 0.25)+
    theme_minimal()

}


# Custom legend -----------------------------------------------------------
custom_legend <- function() {
  
  species_names <- c('Tucuxi', 'Guiana dolphin', 'Unknown sp.')
  species_colors <- c(trait_colors[1],trait_colors[2], trait_colors[3])
  
  # Create a dataframe for the legend
  df_legend <- data.frame(Species = species_names, Color = species_colors)
  
  # Create the legend plot
  legend_plot <- ggplot(df_legend, aes(x = 1, xend = 1.2, y = Species, yend = Species, color = Color)) +
    geom_segment(size = 3) +
    geom_text(aes(x = 1.25, y = Species, label = Species),size=10, hjust = 0) +
    scale_color_identity() +
    theme_void() +
    theme(legend.position = "none")
  
  return(legend_plot)
}



# Read shapefile ----------------------------------------------------------
read_shp <- function(path) {
  shp <- st_read(path)
  return(shp)
}



# Create map --------------------------------------------------------------
create_map <- function(countries) {
  
  # Load the rivers shapefile
  #riversCA <- st_read("./input/Map_data/CA_river/carivs.shp")
  #riversSA <- st_read("./input/Map_data/SA_river/sarivs.shp")
  
  #width_threshold <- 500
  #riversCA <- riversCA[riversCA$a_WIDTH >= width_threshold, ]
  #riversSA <- riversSA[riversSA$a_WIDTH >= width_threshold, ]
  
  # Load the location data
  locations <- read.csv('./input/Map_data/locations.csv')
  locations_df <- as.data.frame(locations)
  
  # Convert the data.frame to a spatial object
  locations <- st_as_sf(locations, coords = c("Longdec", "Latdec"), crs = 4326)
  
  # Plot the map
  g <- ggplot() +
    geom_sf(data = countries, fill = "grey97", color = "black", size = 0.1) +
    #geom_sf(data = riversCA, color='blue') +
    #geom_sf(data = riversSA, color='blue') +
    geom_text_repel(data=locations_df, aes(x=Longdec, y=Latdec, label=Pop, color=Species),size=4, show.legend = FALSE)+
    geom_sf(data = locations, aes(color = Species, shape = Species), size = 2) +
    coord_sf(xlim = c(-85, -35), ylim = c(-26, 12)) +
    scale_color_manual(values=c(trait_colors[1],trait_colors[2],trait_colors[3]),
                       labels=c('Tucuxi', 'Guiana dolphin', 'Unknown'),
                       breaks=c("S. fluviatilis", "S. guianensis", "Unknown"))+
    scale_shape_manual(values = c(16, 17, 18),
                       labels=c('Tucuxi', 'Guiana dolphin', 'Unknown'),
                       breaks=c("S. fluviatilis", "S. guianensis", "Unknown")) +
    labs(color = "Species",
         x='Longitude',
         y='Latitude',
         shape = "Species") +
    guides(size='none', text='none') +
    theme(text = element_text(size = 10),
          panel.border = element_rect(fill=NA),
          panel.background = element_blank(),
          legend.position = c(.975, .975),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          #legend.margin = margin(5, 5, 5, 5),
          legend.box.background=element_rect(),
          legend.key = element_blank())
  
  return(g)
}
  

# Plot regression results -------------------------------------------------
plot_regression <- function(fit, param_name, colnum, xlab){
  
  draws <- merge_chains(as_draws(fit, variable = param_name))
  draws <- as.data.frame(draws)
  
  if(!(param_name %in% names(draws))){
    stop("Don't be silly! param_name is not in draws")
  }  
  
  g <- ggplot(draws, aes(x = .data[[param_name]])) +
    geom_density(fill=trait_colors4[colnum],color=trait_colors4[colnum],alpha=0.25) +
    geom_vline(xintercept = 0, linetype='dashed', linewidth=0.75) +
    labs(x=xlab, y='Density')+
    theme_minimal()
  
  return(g)
  
}


# Rarefaction-Extrapolation curves ----------------------------------------
REEX_curve <- function(iNEXT_output, curveType) {
  g1 <- ggiNEXT(iNEXT_output, type=curveType) + scale_color_manual(values=pub_colors) + scale_fill_manual(values=pub_colors) 
  ggarrange(g1,nrow=1)
}



# Save plot ---------------------------------------------------------------

save_figure <- function(path, w, h, call) {
  png(path, width=w, height=h, units='in', res=800)
  print(call)
  dev.off()
}



printPop <- function(data){
  
  print(unique(data$pop))
  
}
  

# Make a basic plot -------------------------------------------------------
  
plot_xy <- function(df, xcol, ycol) {

  stopifnot("X variable not found" = xcol %in% names(df))
  stopifnot("Y variable not found" = ycol %in% names(df))

  g <- ggplot(df) + geom_point(aes(x=.data[[xcol]], y = .data[[ycol]]))

  stopifnot(class(g)[1]=="gg")

  return(g)

}



# iNEXT output table ------------------------------------------------------
iNEXT_table <- function(q0_result, result_at_cmax, data_subset) {
  
  # start with populations, alphabetical for now
  if (data_subset=="low") (d <- data.table(Assemblage=c('CA','CO','CR','EC','FG','IG','JU','PA','PE','PR','RN','SC','SE','SP','TO','VZ')))
  if (data_subset=="high") (d <- data.table(Assemblage=c('CA','CR','EC','FG','IG','JU','PA','PE','PR','RN','SC','SE','SP','TO')))
  
  # organize asymptotic estimates
  asymptotic_ests <- data.table(q0_result$AsyEst)
  a_q0 <- asymptotic_ests[Diversity=="Species richness"]
  a_q1 <- asymptotic_ests[Diversity=="Shannon diversity"]
  a_q2 <- asymptotic_ests[Diversity=="Simpson diversity"]
  
  # add asymptotic estimates
  d[,asy_0:=a_q0[,round(Estimator),]]
  d[,asy_1:=a_q1[,round(Estimator),]]
  d[,asy_2:=a_q2[,round(Estimator),]]
  
  # add repertoire sizes at cmax
  d[,cmax_q0:=result_at_cmax[Order.q==0,round(qD),]]
  d[,cmax_q1:=result_at_cmax[Order.q==1,round(qD),]]
  d[,cmax_q2:=result_at_cmax[Order.q==2,round(qD),]]
  
  # re-order based on species and alphabet for printed table
  if (data_subset=="low") (order <- data.table(Assemblage=c('CA','CO','EC','JU','PE','CR','FG','IG','PA','PR','RN','SC','SE','SP','VZ','TO')))
  if (data_subset=="high") (order <- data.table(Assemblage=c('CA','EC','JU','PE','CR','FG','IG','PA','PR','RN','SC','SE','SP','TO')))
  
  order[,position:=.I,]
  table <- merge(d, order)
  table <- table[order(position), ]

  return(nice_table(table[,1:7]))
  
}




# Coefficients table (brms) -----------------------------------------------
brms_table <- function(fit, out_path) {

  # create table
  stats.table <- data.table(as.data.frame(summary(fit)$fixed)) 
  stats.table <- cbind(row.names(summary(fit)$fixed), stats.table)
  stats.table[,c('Rhat','Bulk_ESS','Tail_ESS'):=NULL,]
  names(stats.table) <- c("Term", "Estimate", "SE", "CI-Lower", "CI-Upper")

  return(nice_table(stats.table))
  
}



print('Cleared functions')