##### -- multispeciesPP_coef_plot() -- #####
#### Function to plot standardized model coefficients from the various different models of a given species
multispeciesPP_coef_plot <- function(species_name, group = c("bird", "mamm"), mPP_out){
        group <- match.arg(group)
        species_models <- species_coefs <- mPP_out[grep(species_name, names(mPP_out))]
        for (i in seq(along = species_models)){
                species_coefs[[i]] <- data.frame(species_models[[i]][[1]], model = names(species_models)[i])
                species_coefs[[i]] <- subset(species_coefs[[i]], !(species_coefs[[i]]$variable %in% c("(Intercept)", "isPO", "ruggedness", "dist_from_urban", "dist_from_stream", "dist_from_survey")))
        }
        multi_models <- multi_coefs <- mPP_out[grep(paste(group, "multispecies", sep = "_"), names(mPP_out))]
        for (i in seq(along = multi_models)){
                multi_coefs[[i]] <- data.frame(multi_models[[i]][[1]], model = names(multi_models)[i])
                multi_coefs[[i]] <- subset(multi_coefs[[i]], multi_coefs[[i]]$species == species_name & !(multi_coefs[[i]]$variable %in% c("(Intercept)", "isPO", "ruggedness", "dist_from_urban", "dist_from_stream", "dist_from_survey")))
        }        
        species_coefs <- do.call("rbind", c(species_coefs, multi_coefs))
        species_coefs$model <- as.factor(species_coefs$model)
        species_coefs$model <- factor(species_coefs$model, 
                                      levels = levels(species_coefs$model)[c(
                                              which(!grepl("multispecies", levels(species_coefs$model)) & grepl("t1", levels(species_coefs$model)) & grepl("climate", levels(species_coefs$model))),
                                              which(!grepl("multispecies", levels(species_coefs$model)) & grepl("t1", levels(species_coefs$model)) & grepl("habitat", levels(species_coefs$model))),
                                              which(!grepl("multispecies", levels(species_coefs$model)) & grepl("t1", levels(species_coefs$model)) & grepl("full", levels(species_coefs$model))),
                                              which(grepl("multispecies", levels(species_coefs$model)) & grepl("t1", levels(species_coefs$model)) & grepl("climate", levels(species_coefs$model))),
                                              which(grepl("multispecies", levels(species_coefs$model)) & grepl("t1", levels(species_coefs$model)) & grepl("habitat", levels(species_coefs$model))),
                                              which(grepl("multispecies", levels(species_coefs$model)) & grepl("t1", levels(species_coefs$model)) & grepl("full", levels(species_coefs$model))),
                                              which(!grepl("multispecies", levels(species_coefs$model)) & grepl("t2", levels(species_coefs$model)) & grepl("climate", levels(species_coefs$model))),
                                              which(!grepl("multispecies", levels(species_coefs$model)) & grepl("t2", levels(species_coefs$model)) & grepl("habitat", levels(species_coefs$model))),
                                              which(!grepl("multispecies", levels(species_coefs$model)) & grepl("t2", levels(species_coefs$model)) & grepl("full", levels(species_coefs$model))),
                                              which(grepl("multispecies", levels(species_coefs$model)) & grepl("t2", levels(species_coefs$model)) & grepl("climate", levels(species_coefs$model))),
                                              which(grepl("multispecies", levels(species_coefs$model)) & grepl("t2", levels(species_coefs$model)) & grepl("habitat", levels(species_coefs$model))),
                                              which(grepl("multispecies", levels(species_coefs$model)) & grepl("t2", levels(species_coefs$model)) & grepl("full", levels(species_coefs$model)))
                                      )])
        levels(species_coefs$model) <- c("Historic climate-only single-species model", "Historic habitat-only single-species model",
                                         "Historic full single-species model", "Historic climate-only multi-species model",
                                         "Historic habitat-only multi-species model", "Historic full multi-species model",
                                         "Modern climate-only single-species model", "Modern habitat-only single-species model",
                                         "Modern full single-species model", "Modern climate-only multi-species model",
                                         "Modern habitat-only multi-species model", "Modern full multi-species model")
        species_coefs$variable <- as.factor(as.character(species_coefs$variable))
        species_coefs$variable <- factor(species_coefs$variable, levels = names(sort(tapply(abs(species_coefs$estimate), species_coefs$variable, mean), decreasing = TRUE)))
        species_coefs <- species_coefs[order(species_coefs$variable, species_coefs$model), ]
        species_coefs$higher <- species_coefs$estimate + (2 * species_coefs$se)
        species_coefs$lower <- species_coefs$estimate - (2 * species_coefs$se)
        coef_plot <- ggplot(species_coefs, aes(x = model, y = estimate)) +
                geom_bar(aes(fill = model), position = position_dodge(width=0.3), stat="identity", alpha=0) +
                geom_point(aes(color = model), position = position_dodge(width = .8), size = 3) +
                geom_hline(aes(yintercept = 0), linetype = 2) +
                geom_errorbar(aes(ymax = higher, ymin = lower, color = model), position = position_dodge(width = .8), size = 1, width = 0.6) +
                facet_wrap(~ variable) +
                theme_bw() + 
                ylab("Standardized regression coefficient") +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text=element_text(size=14),
                      #strip.text.x = element_blank(),
                      axis.title.y = element_text(size = 14),
                      #legend.position="none",
                      plot.margin=unit(c(0,1,1,1), "cm")
                )
        return(coef_plot)
}

##### -- multispeciesPP_dev_plot() -- #####
#### Function to produce a barplot of model deviance for each species and across all species
multispeciesPP_dev_plot <- function(mPP_out, taxon_name = NA){
        deviance_df <- data.frame(model = names(mPP_out), deviance = unlist(lapply(mPP_out, function(x) x[[2]])))
        deviance_df$predictor_set <- unlist(lapply(strsplit(as.character(deviance_df$model), "_"), function(x) x[length(x)]))
        deviance_df$time_period <- unlist(lapply(strsplit(as.character(deviance_df$model), "_"), function(x) x[length(x) - 1]))
        deviance_df$time_period <- as.factor(deviance_df$time_period)
        deviance_df$group <- unlist(lapply(strsplit(as.character(deviance_df$model), "_"), function(x) x[2]))
        deviance_df$species <- unlist(lapply(strsplit(as.character(deviance_df$model), "_"), function(x) paste(x[c(3, 4)], collapse = "_")))
        deviance_df$species[grepl("multispecies", deviance_df$species)] <- paste(deviance_df$group[grepl("multispecies", deviance_df$species)], "multispecies", sep = "_")
        deviance_df <- subset(deviance_df, species == taxon_name)
        deviance_df$model <- unlist(lapply(strsplit(as.character(deviance_df$model), "_"), function(x) paste(x[-c(length(x)-1, length(x))], collapse = "_")))
        deviance_df$model <- as.factor(deviance_df$model)
        ggplot(deviance_df, aes(x = model, y = deviance)) +
                geom_bar(aes(fill = predictor_set), position = position_dodge(width=1), stat="identity") +
                facet_wrap(~ time_period) +
                theme_bw() + 
                ylab("Unexplained deviance") +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text=element_text(size=14),
                      #strip.text.x = element_blank(),
                      axis.title.y = element_text(size = 14),
                      #legend.position="none",
                      plot.margin=unit(c(0,1,1,1), "cm")
                )
}

##### -- multispeciesPP_eval_plot() -- #####
#### Function to produce a barplot of predictive performance (meawsured using auc and cor) for each species and across all species
multispeciesPP_eval_plot <- function(mPP_eval_output, taxon_name = NA, measure = c("auc", "cor")){
        measure <- match.arg(measure)
        model_names <- names(mPP_eval_output)
        for (i in seq(along = mPP_eval_output)){
                if (grepl("t1", names(mPP_eval_output)[1])){
                        mPP_eval_output[[i]] <- mPP_eval_output[[i]][[2]]
                } else mPP_eval_output[[i]] <- mPP_eval_output[[i]][[1]]
                
                mPP_eval_output[[i]]$model <- strsplit(model_names[i], "\\.")[[1]][1]        
        }
        eval_df <- do.call("rbind", mPP_eval_output)
        eval_df$predictor_set <- unlist(lapply(strsplit(as.character(eval_df$model), "_"), function(x) x[length(x)]))
        eval_df$time_period <- unlist(lapply(strsplit(as.character(eval_df$model), "_"), function(x) x[length(x) - 1]))
        eval_df$time_period <- as.factor(eval_df$time_period)
        eval_df$group <- unlist(lapply(strsplit(as.character(eval_df$model), "_"), function(x) x[2]))
        eval_df$type <- "single species"
        eval_df$type[grepl("multispecies", eval_df$model)] <- "multispecies"
        eval_df$type <- as.factor(eval_df$type)
        eval_df <- subset(eval_df, species == taxon_name)
        eval_df$model <- unlist(lapply(strsplit(as.character(eval_df$model), "_"), function(x) paste(x[-c(length(x)-1, length(x))], collapse = "_")))
        eval_df$model <- as.factor(eval_df$model)
        ggplot(eval_df, aes_string(x = "predictor_set", y = measure)) +
                geom_bar(aes(fill = predictor_set), position = position_dodge(width=1), stat="identity") +
                facet_wrap(~ time_period + type) +
                ylim(c(0, 1)) +
                theme_bw() + 
                ylab("Unexplained deviance") +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), 
                      axis.title.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text=element_text(size=14),
                      #strip.text.x = element_blank(),
                      axis.title.y = element_text(size = 14),
                      #legend.position="none",
                      plot.margin=unit(c(0,1,1,1), "cm")
                )
}