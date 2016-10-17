######################################################
##### -- climate-vs-habitat-change-california -- #####
######################################################
##################### FUNCTIONS ######################
##### -- multispeciesPP_edit() -- #####
##### Modified multispeciesPP function: function was edited to return fit object, which includes estimates of residual deviance
##### R code from multispeciesPP by Will Fithian (https://github.com/wfithian/multispeciesPP/blob/master/R/multispeciesPP.R).
##### For more information, see 
##### Fithian et al. (2014) Bias correction in species distribution models: pooling survey and collection data for multiple species. Methods in Ecology and Evolution
multispeciesPP_edit <- 
        function (sdm.formula, bias.formula, PA, PO, BG, species = names(PO), 
                  species.PA = species, species.PO = species, quadrat.size = 1, 
                  region.size = 1, start = NULL, inverse.hessian = FALSE, penalty.l2.sdm = 0.1, 
                  penalty.l2.bias = 0.1, penalty.l2.intercept = 1e-04, weights = rep(1, 
                                                                                     n.species * nrow(x)), control = list()) 
        {
                control <- do.call("glm.control", control)
                species <- union(species.PO, species.PA)
                sdm.formula <- update(sdm.formula, ~. + 1)
                bias.formula <- update(bias.formula, ~. - 1)
                sdm.mf <- model.frame(sdm.formula, data = BG)
                bias.mf <- model.frame(bias.formula, data = BG)
                sdm.BG.model.matrix <- model.matrix(terms(sdm.mf), BG)
                sdm.means <- c(0, apply(sdm.BG.model.matrix[, -1, drop = FALSE], 
                                        2, mean))
                sdm.BG.model.matrix <- sweep(sdm.BG.model.matrix, 2, sdm.means, 
                                             "-")
                sdm.sds <- c(1, apply(sdm.BG.model.matrix[, -1, drop = FALSE], 
                                      2, sd))
                sdm.BG.model.matrix <- sweep(sdm.BG.model.matrix, 2, sdm.sds, 
                                             "/")
                sdm.standardize <- function(mat) sweep(sweep(mat, 2, sdm.means, 
                                                             "-"), 2, sdm.sds, "/")
                bias.BG.model.matrix <- model.matrix(terms(bias.mf), BG)
                bias.means <- apply(bias.BG.model.matrix, 2, mean)
                bias.BG.model.matrix <- sweep(bias.BG.model.matrix, 2, bias.means, 
                                              "-")
                bias.sds <- apply(bias.BG.model.matrix, 2, sd)
                bias.BG.model.matrix <- sweep(bias.BG.model.matrix, 2, bias.sds, 
                                              "/")
                bias.standardize <- function(mat) sweep(sweep(mat, 2, bias.means, 
                                                              "-"), 2, bias.sds, "/")
                BG.good.rows <- intersect(rownames(sdm.BG.model.matrix), 
                                          rownames(bias.BG.model.matrix))
                sdm.PA.model.matrix <- sdm.standardize(model.matrix(terms(sdm.mf), 
                                                                    PA))
                PA.good.rows <- rownames(sdm.PA.model.matrix)
                if (!is.null(species.PO)) {
                        sdm.PO.model.matrices <- lapply(as.list(species.PO), 
                                                        function(sp) sdm.standardize(model.matrix(terms(sdm.mf), 
                                                                                                  PO[[sp]])))
                        names(sdm.PO.model.matrices) <- species.PO
                        bias.PO.model.matrices <- lapply(as.list(species.PO), 
                                                         function(sp) bias.standardize(model.matrix(terms(bias.mf), 
                                                                                                    PO[[sp]])))
                        names(bias.PO.model.matrices) <- species.PO
                        PO.good.rows <- lapply(as.list(species.PO), function(sp) intersect(rownames(sdm.PO.model.matrices[[sp]]), 
                                                                                           rownames(bias.PO.model.matrices[[sp]])))
                        names(PO.good.rows) <- species.PO
                }
                n.species <- length(species)
                p.sdm <- ncol(sdm.BG.model.matrix) - 1
                p.bias <- ncol(bias.BG.model.matrix)
                sdm.margins.ab <- matrix(0, n.species, p.sdm + 1, dimnames = list(species, 
                                                                                  colnames(sdm.BG.model.matrix)))
                sdm.margins.gamma <- matrix(0, n.species, 1, dimnames = list(species, 
                                                                             "isPO"))
                bias.margins <- matrix(0, 1, p.bias, dimnames = list(NULL, 
                                                                     colnames(bias.BG.model.matrix)))
                for (sp in species.PO) {
                        k <- match(sp, species)
                        sdm.margins.ab[k, ] <- colSums(sdm.PO.model.matrices[[sp]][PO.good.rows[[sp]], 
                                                                                   , drop = FALSE])
                        sdm.margins.gamma[k, ] <- length(PO.good.rows[[sp]])
                        bias.margins <- bias.margins + colSums(bias.PO.model.matrices[[sp]][PO.good.rows[[sp]], 
                                                                                            , drop = FALSE])
                }
                abcd.from.all.coef <- function(all.coef) {
                        sdm.coef <- matrix(all.coef[1:(n.species * (p.sdm + 2))], 
                                           p.sdm + 2, n.species)
                        alpha <- sdm.coef[1, ]
                        beta <- t(sdm.coef[2:(p.sdm + 1), , drop = FALSE])
                        gamma <- sdm.coef[p.sdm + 2, ]
                        delta <- all.coef[-(1:(n.species * (p.sdm + 2)))]
                        names(alpha) <- names(gamma) <- species
                        colnames(beta) <- colnames(sdm.margins.ab)[-1]
                        rownames(beta) <- species
                        names(delta) <- colnames(bias.BG.model.matrix)
                        return(list(alpha = alpha, beta = beta, gamma = gamma, 
                                    delta = delta))
                }
                all.coef.from.abcd <- function(alpha, beta, gamma, delta) {
                        c(rbind(alpha, beta, gamma), delta)
                }
                n.PA <- length(PA.good.rows)
                n.BG <- length(BG.good.rows)
                subsamp.PA.offset <- 0
                subsamp.BG.offset <- 0
                n.sites <- n.BG + n.PA
                x <- cbind(rbind(sdm.margins.ab, 0, sdm.PA.model.matrix[PA.good.rows, 
                                                                        , drop = FALSE], sdm.BG.model.matrix[BG.good.rows, , 
                                                                                                             drop = FALSE]), c(sdm.margins.gamma, rep(0:1, c(1 + n.PA, 
                                                                                                                                                             n.BG))))
                x <- rbind(x, diag(sqrt(c(penalty.l2.intercept, rep(penalty.l2.sdm, 
                                                                    p.sdm), penalty.l2.intercept))), matrix(0, p.bias, p.sdm + 
                                                                                                                    2))
                z <- rbind(matrix(0, n.species, p.bias), bias.margins, matrix(0, 
                                                                              n.PA, p.bias), bias.BG.model.matrix[BG.good.rows, , drop = FALSE], 
                           matrix(0, p.sdm + 2, p.bias), sqrt(penalty.l2.bias/n.species) * 
                                   diag(p.bias))
                y <- rep(0, nrow(x) * n.species)
                offset <- rep(0, nrow(x) * n.species)
                for (k in 1:n.species) {
                        yk <- rep(0, nrow(x))
                        yk[1:n.species] <- 1 * (1:n.species == k)
                        yk[1 + n.species] <- 1 * (1 == k)
                        if (species[k] %in% species.PA) {
                                yk[1 + n.species + (1:n.PA)] <- PA[PA.good.rows, 
                                                                   species[k]]
                        }
                        else {
                                yk[1 + n.species + (1:n.PA)] <- NA
                        }
                        if (species[k] %in% species.PO) {
                                yk[1 + n.species + n.PA + (1:n.BG)] <- 0
                        }
                        else {
                                yk[1 + n.species + n.PA + (1:n.BG)] <- NA
                        }
                        yk[1 + n.species + n.sites + (1:(p.sdm + 2 + p.bias))] <- 0
                        y[(k - 1) * nrow(x) + 1:nrow(x)] <- yk
                        offk <- rep(0, nrow(x))
                        offk[1 + n.species + (1:n.PA)] <- log(quadrat.size)
                        offk[1 + n.species + n.PA + (1:n.BG)] <- log(region.size) - 
                                log(n.BG)
                        offset[(k - 1) * nrow(x) + 1:nrow(x)] <- offk
                }
                which.PA <- (2 + n.species):(1 + n.species + n.PA) + rep((0:(n.species - 
                                                                                     1)) * nrow(x), each = n.PA)
                which.BG <- (2 + n.species + n.PA):(1 + n.species + n.PA + 
                                                            n.BG) + rep((0:(n.species - 1)) * nrow(x), each = n.BG)
                if (is.null(start)) {
                        start.alpha <- start.gamma <- rep(0, n.species)
                        for (k in 1:n.species) {
                                if ((species[k] %in% species.PA) && sum(!is.na(PA[PA.good.rows, 
                                                                                  species[k]]) > 0)) 
                                        start.alpha[k] <- log((1 + sum(PA[PA.good.rows, 
                                                                          species[k]], na.rm = TRUE))/n.PA/quadrat.size)
                                if (species[k] %in% species.PO) 
                                        start.gamma[k] <- log1p(sdm.margins.gamma[k, 
                                                                                  ]) - start.alpha[k] - log(region.size)
                        }
                        start <- all.coef.from.abcd(start.alpha, matrix(0, p.sdm, 
                                                                        n.species), start.gamma, rep(0, p.bias))
                }
                fit <- block.glm.fit(x, z, y, weights = weights, start = start, 
                                     offset = offset, families = list(linear(), binomial(link = "cloglog"), 
                                                                      poisson(), gaussian()), row.families = rep(rep(1:4, 
                                                                                                                     c(1 + n.species, n.PA, n.BG, p.sdm + p.bias + 2)), 
                                                                                                                 n.species), control = control)
                all.coef <- fit$coefficients
                eta <- fit$linear.predictors
                mu <- fit$fitted.values
                names(all.coef)[1:(n.species * (p.sdm + 2))] <- paste(rep(species, 
                                                                          each = p.sdm + 2), c(colnames(sdm.BG.model.matrix)[1:(p.sdm + 
                                                                                                                                        1)], "isPO"), sep = ":")
                names(all.coef)[-(1:(n.species * (p.sdm + 2)))] <- paste("isPO:", 
                                                                         colnames(bias.BG.model.matrix), sep = "")
                std.errs <- fit$fit$std.errs
                names(std.errs) <- names(all.coef)
                species.coef <- matrix(all.coef[1:(n.species * (p.sdm + 2))], 
                                       p.sdm + 2, n.species, dimnames = list(c(colnames(sdm.margins.ab), 
                                                                               "isPO"), species))
                bias.coef <- all.coef[-(1:(n.species * (p.sdm + 2)))]
                names(bias.coef) <- colnames(bias.BG.model.matrix)
                fit.PA <- linear.fit.PA <- matrix(NA, nrow(PA), length(species), 
                                                  dimnames = list(dimnames(PA)[[1]], species))
                linear.fit.PA[PA.good.rows, ] <- eta[which.PA]
                fit.PA[PA.good.rows, ] <- mu[which.PA]
                fit.BG <- linear.fit.BG <- bias.fit.BG <- linear.bias.fit.BG <- matrix(NA, 
                                                                                       nrow(BG), length(species), dimnames = list(dimnames(BG)[[1]], 
                                                                                                                                  species))
                linear.fit.BG[BG.good.rows, ] <- matrix(eta[which.BG], ncol = n.species) + 
                        log(n.BG) - log(region.size)
                fit.BG[BG.good.rows, ] <- matrix(mu[which.BG], ncol = n.species) * 
                        n.BG/region.size
                linear.bias.fit.BG[BG.good.rows, ] <- c(bias.BG.model.matrix[BG.good.rows, 
                                                                             , drop = FALSE] %*% bias.coef)
                bias.fit.BG[BG.good.rows, ] <- exp(linear.bias.fit.BG[BG.good.rows, 
                                                                      ])
                fitted.sdm.margins.gamma <- colSums(fit.BG[BG.good.rows, 
                                                           , drop = FALSE]) * region.size/n.BG
                fitted.bias.margins <- colSums(t(fit.BG[BG.good.rows, species.PO, 
                                                        drop = FALSE]) %*% bias.BG.model.matrix[BG.good.rows, 
                                                                                                , drop = FALSE] * region.size/n.BG)
                score.check.gamma <- fitted.sdm.margins.gamma - sdm.margins.gamma + 
                        penalty.l2.intercept * species.coef[p.sdm + 2, ]
                score.check.gamma <- score.check.gamma[species %in% species.PO]
                score.check.bias <- fitted.bias.margins - bias.margins + 
                        penalty.l2.bias * bias.coef
                if (length(score.check.gamma) > 0) 
                        stopifnot(mean((score.check.gamma/fit$deviance)^2) < 
                                          control$epsilon)
                stopifnot(mean((score.check.bias/fit$deviance)^2) < control$epsilon)
                sd.normalizer <- c(rep(c(sdm.sds, 1), n.species), bias.sds)
                unstandardized.coef <- all.coef/sd.normalizer
                gamma.adjust <- sum(unstandardized.coef[-(1:(n.species * 
                                                                     (p.sdm + 2)))] * bias.means)
                for (k in 1:n.species) {
                        jk <- (p.sdm + 2) * (k - 1) + 1:(p.sdm + 1)
                        coef.block <- unstandardized.coef[jk]
                        unstandardized.coef[jk[1]] <- coef.block[1] - sum(coef.block[-1] * 
                                                                                  sdm.means[-1])
                        unstandardized.coef[jk[1] + p.sdm + 1] <- unstandardized.coef[jk[1] + 
                                                                                              p.sdm + 1] - gamma.adjust
                }
                unstandardized.species.coef <- matrix(unstandardized.coef[1:(n.species * 
                                                                                     (p.sdm + 2))], p.sdm + 2, n.species, dimnames = list(c(colnames(sdm.margins.ab), 
                                                                                                                                            "isPO"), species))
                unstandardized.bias.coef <- unstandardized.coef[-(1:(n.species * 
                                                                             (p.sdm + 2)))]
                names(unstandardized.bias.coef) <- colnames(bias.BG.model.matrix)
                tr <- list(sdm.formula = sdm.formula, bias.formula = bias.formula, fit = fit, 
                           normalized.species.coef = species.coef, normalized.bias.coef = bias.coef, 
                           normalized.all.coef = all.coef, normalized.std.errs = std.errs, 
                           all.coef = unstandardized.coef, std.errs = std.errs/sd.normalizer, 
                           species.coef = unstandardized.species.coef, bias.coef = unstandardized.bias.coef, 
                           linear.fit.PA = linear.fit.PA, fit.PA = fit.PA, linear.bias.fit.BG = linear.bias.fit.BG, 
                           bias.fit.BG = bias.fit.BG, linear.fit.BG = linear.fit.BG, 
                           fit.BG = fit.BG)
                class(tr) <- c("multispeciesPP", "list")
                tr
        }

##### -- multispeciesPP_wrapper() -- #####
##### Wrapper around function multispeciesPP() from library(multispeciesPP) to facilitate running of models with different types of information
##### R code from multispeciesPP by Will Fithian (https://github.com/wfithian/multispeciesPP/blob/master/R/multispeciesPP.R).
##### For more information, see 
##### Fithian et al. (2014) Bias correction in species distribution models: pooling survey and collection data for multiple species. Methods in Ecology and Evolution
multispeciesPP_wrapper <- function(pa_data = NULL, 
                                   po_data = NULL,
                                   bg = NULL,
                                   species_names = NULL,
                                   climate_predictors = paste("bio", c(1, 6, 12), sep = ""),
                                   habitat_associations = NULL,
                                   group = c("bird", "mamm", "odon"), ## Taxonomic group to model (birds/mammals/odonates)
                                   predictor_set = c("climate", "habitat", "full"), ## Use only climate, habitat, or both (full) as model predictors
                                   out_name = "out",
                                   ...){
        ### Match function arguments
        group <- match.arg(group)
        predictor_set <- match.arg(predictor_set)
        ### Create directory to save model output
        dir.create(paste(getwd(), "/output/multispeciesPP", sep = ""), showWarnings = FALSE)
        dir.create(paste(getwd(), "/output/multispeciesPP/models", sep = ""), showWarnings = FALSE)
        ### Generate useful objects
        ## Character vector of climate predictors
        climate_pred <- climate_predictors
        ## Character vector of habitat predictors
        habitat_pred <- habitat_associations
        ## Character vector of bias predictors
        bias_pred <- c("ruggedness", "dist_from_urban", "dist_from_stream", "dist_from_survey")
        ## Size of study area
        study_area <- nrow(bg)
        ### Pick appropriate variables from background object
        bg <- bg[c(intersect(names(bg), c(climate_pred, habitat_pred, bias_pred)), paste("dist_from_survey", group, sep = "_"))]
        names(bg)[grep("survey", names(bg))] <- "dist_from_survey"
        ### Select the desired species set from pa_data and po_data, if necessary
        if (!is.null(pa_data)) pa_data <- pa_data[c(species_names, climate_pred, habitat_pred)] 
        if (!is.null(po_data)) po_data <- po_data[species_names] 
        ### Standardize covariates
        if (!is.null(pa_data)) pa_data[, c(climate_pred, habitat_pred)] <- apply(pa_data[, c(climate_pred, habitat_pred)], 2, scale) %>% data.frame()
        if (!is.null(po_data)) po_data <- lapply(po_data, function(x) apply(x[c(climate_pred, habitat_pred, bias_pred)], 2, scale) %>% data.frame())
        bg[, c(climate_pred, habitat_pred, bias_pred)] <- apply(bg[c(climate_pred, habitat_pred, bias_pred)], 2, scale) %>% data.frame()
        ### Specify formulas
        climate_pred <- paste(climate_pred, collapse = " + ")
        habitat_pred <- paste(habitat_pred, collapse = " + ")
        bias_pred <- paste(bias_pred, collapse = " + ")
        ## Bias formula
        bias_formula <- as.formula(paste("~ ", bias_pred, sep = ""))
        ## SDM formula
        if (predictor_set == "full"){
                sdm_formula <- as.formula(paste("~ ", climate_pred, " + ", habitat_pred, sep = ""))
        }
        if (predictor_set == "climate"){
                sdm_formula <- as.formula(paste("~ ", climate_pred, sep = ""))
        }
        if (predictor_set == "habitat"){
                sdm_formula <- as.formula(paste("~ ", habitat_pred, sep = ""))
        }
        ### Run model
        mPP <- multispeciesPP_edit(
                sdm.formula = sdm_formula,
                bias.formula = bias_formula,
                PA = pa_data,
                PO = po_data,
                BG = bg,
                region.size = study_area,
                ...
        )
        ### Save output
        saveRDS(mPP, file = paste("output/multispeciesPP/models/mPP_", out_name, ".rds", sep = ""))
}

##### -- multispeciesPP_output() -- #####
##### Extract useful output from saved multispeciesPP models
multispeciesPP_output <- function(mPP_directory = "output/multispeciesPP/models/"){
        mPP_list <- list.files(mPP_directory)
        mPP_out <- lapply(mPP_list, function(x){
                mPP <- readRDS(paste(mPP_directory, x, sep = ""))
                # Coefficients
                coefs <- mPP$normalized.all.coef
                se <- mPP$normalized.std.errs
                summary <- data.frame(coefs, se, coefs/se, 2*pnorm(-abs(coefs/se)))                
                colnames(summary) <- c("estimate","se","z","p")
                summary$species <- factor(unlist(lapply(strsplit(row.names(summary), ':'), function(y) y[1])))
                summary$variable <- unlist(lapply(strsplit(row.names(summary), ':'), function(y) y[2]))
                summary$model <- x
                list(summary = summary,
                     deviance = mPP$fit$deviance
                )
        }
        )
        names(mPP_out) <- unlist(lapply(strsplit(mPP_list, "\\."), function(x) x[[1]]))
        return(mPP_out)
}

###############
#### roc() ####
###############
"roc" <-
        function (obsdat, preddat) 
        {
                # code adapted from Ferrier, Pearce and Watson's code, by J.Elith
                #
                # see:
                # Hanley, J.A. & McNeil, B.J. (1982) The meaning and use of the area
                # under a Receiver Operating Characteristic (ROC) curve.
                # Radiology, 143, 29-36
                #
                # Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance
                # of habitat models developed using logistic regression.
                # Ecological Modelling, 133, 225-245.
                # this is the non-parametric calculation for area under the ROC curve, 
                # using the fact that a MannWhitney U statistic is closely related to
                # the area
                #
                if (length(obsdat) != length(preddat)) 
                        stop("obs and preds must be equal lengths")
                n.x <- length(obsdat[obsdat == 0])
                n.y <- length(obsdat[obsdat == 1])
                xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
                rnk <- rank(xy)
                wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * 
                                                                                         n.y)
                return(round(wilc, 4))
        }
#########################
#### calc_deviance() ####
#########################
"calc_deviance" <-
        function(obs.values, fitted.values, weights = rep(1,length(obs.values)), family="binomial", calc.mean = TRUE)
        {
                # j. leathwick/j. elith
                #
                # version 2.1 - 5th Sept 2005
                #
                # function to calculate deviance given two vectors of raw and fitted values
                # requires a family argument which is set to binomial by default
                #
                #
                
                if (length(obs.values) != length(fitted.values)) 
                        stop("observations and predictions must be of equal length")
                
                y_i <- obs.values
                
                u_i <- fitted.values
                
                if (family == "binomial" | family == "bernoulli") {
                        
                        deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
                        deviance <- -2 * sum(deviance.contribs * weights)
                        
                }
                
                if (family == "poisson" | family == "Poisson") {
                        
                        deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
                        deviance <- 2 * sum(deviance.contribs * weights)
                        
                }
                if (family == "laplace") {
                        deviance <- sum(abs(y_i - u_i))
                }
                
                if (family == "gaussian") {
                        deviance <- sum((y_i - u_i) * (y_i - u_i))
                }
                
                
                
                if (calc.mean) deviance <- deviance/length(obs.values)
                
                return(deviance)
                
        }
#######################
##### eval_pred() #####
#######################
eval_pred <- function(obs_table = NA, pred_table = NA, species_names = NA){
        eval_table <- data.frame(species = species_names, auc = NA, cor = NA, dev = NA)
        if (nrow(eval_table) > 1){
                for (i in seq(along = species_names)){
                        obs <- obs_table[, grep(species_names[i], names(obs_table))]
                        pred <- pred_table[, grep(species_names[i], names(pred_table))]
                        eval_table$auc[i] <- roc(obs, pred)
                        eval_table$dev[i] <- calc_deviance(obs, pred)
                        eval_table$cor[i] <- cor(obs, pred, use = "complete.obs", method = "pearson")
                }                
        } else {
                obs <- obs_table[, grep(species_names, names(obs_table))]
                pred <- pred_table[, grep(species_names, names(pred_table))]
                eval_table$auc <- roc(obs, pred)
                eval_table$dev <- calc_deviance(obs, pred)
                eval_table$cor <- cor(obs, pred, use = "complete.obs", method = "pearson")                
        }
        return(eval_table)     
}

multispeciesPP_predictions <- function(mPP_directory = "output/multispeciesPP/models/"){
        # Create directory to save model predictions
        dir.create(paste(getwd(), "/output/multispeciesPP/predictions", sep = ""), showWarnings = FALSE)
        mPP_list <- list.files(mPP_directory)
        mPP_eval_output <- vector('list', length(mPP_list))
        for(i in seq(along = mPP_eval_output)){
                mPP <- readRDS(paste(mPP_directory, mPP_list[i], sep = ""))
                if (grepl("bird", mPP_list[i])){
                        t1_pa <- t1_pa_bird
                        t2_pa <- t2_pa_bird                 
                }
                if (grepl("mamm", mPP_list[i])){
                        t1_pa <- t1_pa_mamm
                        t2_pa <- t2_pa_mamm                 
                } 
                # t1_bg predictions
                predictions_t1_bg <- data.frame(t1_bg[c('longitude', 'latitude')], (1 - exp(-exp(predict.multispeciesPP(mPP, newdata = t1_bg)))))
                # t2_bg predictions
                predictions_t2_bg <- data.frame(t2_bg[c('longitude', 'latitude')], (1 - exp(-exp(predict.multispeciesPP(mPP, newdata = t2_bg)))))
                # t1_pa predictions
                predictions_t1_pa <- data.frame(t1_pa[c('longitude', 'latitude')], (1 - exp(-exp(predict.multispeciesPP(mPP, newdata = t1_pa)))))
                # t2_pa predictions
                predictions_t2_pa <- data.frame(t2_pa[c('longitude', 'latitude')], (1 - exp(-exp(predict.multispeciesPP(mPP, newdata = t2_pa)))))     
                # save predictions
                #saveRDS(predictions_t1_bg, paste("output/multispeciesPP/", strsplit(mPP_list[i], "\\.")[[1]][1], '_bg.rds', sep = ''))
                #saveRDS(predictions_t1_bg, paste('output/multispeciesPP/predictions/', strsplit(mPP_list[i], "\\.")[[1]][1], '_bg.rds', sep = ''))
                #saveRDS(predictions_t2_bg, paste('output/multispeciesPP/predictions/', strsplit(mPP_list[i], "\\.")[[1]][1], '_bg.rds', sep = ''))
                #saveRDS(predictions_t1_pa, paste('output/multispeciesPP/predictions/', strsplit(mPP_list[i], "\\.")[[1]][1], '_pa.rds', sep = ''))
                #saveRDS(predictions_t2_pa, paste('output/multispeciesPP/predictions/', strsplit(mPP_list[i], "\\.")[[1]][1], '_pa.rds', sep = ''))        
                #saveRDS(predictions_change_bg, paste('output/multispeciesPP/predictions/', strsplit(mPP_list[i], "\\.")[[1]][1], '_predictions_change_bg.rds', sep = ''))                
                eval_t1_pa <- eval_pred(obs_table = t1_pa, pred_table = predictions_t1_pa, species_names = colnames(mPP$normalized.species.coef)) 
                eval_t2_pa <- eval_pred(obs_table = t2_pa, pred_table = predictions_t2_pa, species_names = colnames(mPP$normalized.species.coef))         
                mPP_eval_output[[i]] <- list(eval_t1_pa = eval_t1_pa, eval_t2_pa = eval_t2_pa)
                names(mPP_eval_output) <- mPP_list
                rm(mPP)
        }
        saveRDS(mPP_eval_output, 'output/multispeciesPP/mPP_eval_output.rds')
        return(mPP_eval_output)
}


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
        for (i in seq(along = mPP_eval_output)){
                if (grepl("t1", names(mPP_eval_output)[1])){
                        mPP_eval_output[[i]] <- mPP_eval_output[[i]][[2]]
                } else mPP_eval_output[[i]] <- mPP_eval_output[[i]][[1]]
                
                mPP_eval_output[[i]]$model <- strsplit(names(mPP_eval_output)[i], "\\.")[[1]][1]        
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