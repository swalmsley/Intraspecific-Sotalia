
library(targets)

# Source
tar_source('R') # will do all in 'R' folder

# Options
tar_option_set()

# Seed
tar_option_set(seed = 1234)

# Variables
pop_colors <- viridis(option = "plasma",begin = 0.15, end = 0.8, n=16)   
trait_colors <- viridis(option = "plasma",begin = 0.15, end = 0.8, n=3) 
trait_colors4 <- viridis(option = "plasma",begin = 0.15, end = 0.8, n=4) 
species_colors <- c(trait_colors[1],trait_colors[1], trait_colors[2], trait_colors[1], trait_colors[2], trait_colors[2], trait_colors[1], trait_colors[2], trait_colors[1], trait_colors[2], trait_colors[2], trait_colors[2], trait_colors[2], trait_colors[2], trait_colors[3], trait_colors[2]) 

# Targets
list(
  
  # load data
  tar_target(raw, read_data('./input/Sotalia_covariates_final_samprates.csv')),
  tar_target(whistles, read_data('./input/whistleDetails_2023.csv')),
  tar_target(geo_distances, read_data('./input/distances_matrix.csv')),
  tar_target(groupSizeData, read_data('./input/groupSizes.csv')),
  
  
  # data processing
  tar_target(data_low, process_data(raw, whistles, 'low')),
  tar_target(data_high, process_data(raw, whistles, 'high')),
  tar_target(iNEXT_data, format_for_iNEXT(data_low)),
  tar_target(iNEXT_data_high, format_for_iNEXT(data_high)), # HF
  
  
  # run iNEXT
  tar_target(q0, q <- iNEXT(iNEXT_data, q=0, datatype='abundance')),
  tar_target(q1, q <- iNEXT(iNEXT_data, q=1, datatype='abundance')),
  tar_target(q2, q <- iNEXT(iNEXT_data, q=2, datatype='abundance')),
  tar_target(q0_high, q <- iNEXT(iNEXT_data_high, q=0, datatype='abundance')), # HF
  
  
  # extract cMax
  tar_target(cmax, extract_cMax(q0)),
  tar_target(cmax_high, extract_cMax(q0_high)), # HF
  
  
  
  # estimate diversity at cMax
  tar_target(estimates_raw, data.table(estimateD(iNEXT_data, 
                                             datatype = 'abundance',
                                             base='coverage',
                                             level=cmax,
                                             conf=0.95))),
  tar_target(estimates_raw_high, data.table(estimateD(iNEXT_data_high,  # HF
                                                 datatype = 'abundance',
                                                 base='coverage',
                                                 level=cmax_high,
                                                 conf=0.95))), 
  
  
  # calculate standard errors for modelling
  # note that this assumes that estimate follows a normal distribution
  tar_target(estimates, estimates_raw[,se:=((qD.UCL - qD.LCL)/(2*1.96)),]),
  tar_target(estimates_high, estimates_raw_high[,se:=((qD.UCL - qD.LCL)/(2*1.96)),]),
  
  
  # Spade R analyses
  tar_target(spadeData, build_spadeR_data(data_low)),
  tar_target(spadeQ1, SimilarityMult(data.frame(spadeData),"abundance",q=1,nboot=200,"relative")),
  tar_target(spadeQ2, SimilarityMult(data.frame(spadeData),"abundance",q=2,nboot=200,"relative")),
  
  tar_target(spadeQ1_tucuxi, SimilarityMult(data.frame(spadeData[,c('CA','CO','EC','JU','PE')]),'abundance',q=1,nboot=200,'relative')),
  tar_target(spadeQ2_tucuxi, SimilarityMult(data.frame(spadeData[,c('CA','CO','EC','JU','PE')]),'abundance',q=2,nboot=200,'relative')),
  
  tar_target(spadeQ1_guiana, SimilarityMult(data.frame(spadeData[,c('CR','FG','IG','PA','PR','RN','SC','SE','SP','VZ')]),'abundance',q=1,nboot=200,'relative')),
  tar_target(spadeQ2_guiana, SimilarityMult(data.frame(spadeData[,c('CR','FG','IG','PA','PR','RN','SC','SE','SP','VZ')]),'abundance',q=2,nboot=200,'relative')),
  
  
  # Mantel tests
  tar_target(mantel.horn, run_mantel(geo_distances, spadeData, spadeQ1, 'horn')),
  tar_target(mantel.morisata, run_mantel(geo_distances, spadeData, spadeQ2, 'morisata')),
  
  
  # Permutation test
  tar_target(permutation.test, run_permutation(data_low)),
  
  
  # NMDS analyses
  tar_target(nmds.horn, run_NMDS(spadeData, spadeQ1, 'horn')),
  tar_target(nmds.morisata, run_NMDS(spadeData, spadeQ2, 'morisata')),
  
  tar_target(plot_horn, plot_nmds(spadeData, nmds.horn, 'Horn')),
  tar_target(plot_morisata, plot_nmds(spadeData, nmds.morisata, 'Morisata')),
  tar_target(Fig2, save_figure('./figures/Figure2.png', w=12,h=6,(plot_horn | plot_morisata) + plot_layout(guides = "collect"))),
  
  
  # sample-size based REEX curves
  tar_target(g1_c, plot_REEX(q0, 0, 1)),
  tar_target(g2_c, plot_REEX(q1, 1, 1)),
  tar_target(g3_c, plot_REEX(q2, 2, 1)),
  tar_target(Fig3, save_figure('./figures/Figure3.png',w=12,h=5,((g1_c | g2_c | g3_c)+plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16))))),

  
  # prepare dataframe for modelling
  tar_target(model_dt, prepModelData(data_low, groupSizeData, estimates)), 
  
  # run regressions
  tar_target(m1, brm(qD | mi(se) ~ number.encounters.analyzed, data=model_dt[Order.q==1,,], family = 'gaussian')),
  tar_target(m2, brm(qD | mi(se) ~ number.encounters.analyzed, data=model_dt[Order.q==2,,], family = 'gaussian')),
  tar_target(m3,brm(qD | mi(se) ~ Species, data=model_dt[Order.q==1 & Species!='Sotalia sp..',,], family = 'gaussian')),
  tar_target(m4,brm(qD | mi(se) ~ Species, data=model_dt[Order.q==2 & Species!='Sotalia sp..',,], family = 'gaussian')),
  tar_target(m5, brm(qD | mi(se) ~ Habitat, data=model_dt[Order.q==1,,], family = 'gaussian')),
  tar_target(m6, brm(qD | mi(se) ~ Habitat, data=model_dt[Order.q==2,,], family = 'gaussian')),
  tar_target(m7, brm(qD | mi(se) ~ as.numeric(median_group), data=model_dt[Order.q==1 & (!is.na(median_group)),,], family = 'gaussian')),
  tar_target(m8, brm(qD | mi(se) ~ as.numeric(median_group), data=model_dt[Order.q==2 & (!is.na(median_group)),,], family = 'gaussian')),
  
  # Basic diagnostics:
  ## ppcheck(mX)
  ## plot(mX)
  
  
  # plot regression results
  tar_target(plot_m1, plot_regression(m1, 'b_number.encounters.analyzed', 1, 'Effect of number of encounters analyzed')),
  tar_target(plot_m2, plot_regression(m2, 'b_number.encounters.analyzed', 1, 'Effect of number of encounters analyzed')),
  tar_target(plot_m3, plot_regression(m3, 'b_SpeciesSotaliaguianensis', 2, 'Contrast: S Guianensis')),
  tar_target(plot_m4, plot_regression(m4, 'b_SpeciesSotaliaguianensis', 2, 'Contrast: S Guianensis')),
  tar_target(plot_m5, plot_regression(m5, 'b_HabitatRiverine', 3, 'Contrast: River')),
  tar_target(plot_m6, plot_regression(m6, 'b_HabitatRiverine', 3, 'Contrast: River')),
  tar_target(plot_m7, plot_regression(m7, 'b_as.numericmedian_group', 4, 'Effect of group size')),
  tar_target(plot_m8, plot_regression(m8, 'b_as.numericmedian_group', 4, 'Effect of group size')),

  
  # create and save map
  tar_target(countries, read_shp("./input/Map_data/countries.shp")),
  tar_target(Fig1, save_figure('./figures/Figure1.png',w=6,h=6,(create_map(countries)))),
  
  
  # save regression figure
  tar_target(Fig4, save_figure('./figures/Figure4.png',w=10,h=2.5,((plot_m2 | plot_m4 | plot_m6 | plot_m8) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16))))),
  tar_target(FigS2, save_figure('./figures/FigureS2.png',w=10,h=2.5,((plot_m1 | plot_m3 | plot_m5 | plot_m7) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16))))),
  
  
  # write tables to supplement
  tar_quarto(
    render,
    file.path('supplement', 'supplement.qmd')
  )
  
  
)


