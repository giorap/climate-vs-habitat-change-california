######################################################
##### -- climate-vs-habitat-change-california -- #####
######################################################
###################### MODELS ########################
##### -- Birds -- #####
##### Multispecies models
#### Historic models (t1)
### Full predictor set
multispeciesPP_wrapper(pa_data = t1_pa_bird,
                       po_data = t1_po_bird,
                       bg = t1_bg,
                       species_names = final_species_set_bird,
                       habitat_associations = habitat_associations_bird,
                       group = "bird", predictor_set = "full", out_name = "bird_multispecies_t1_full")
### Climate only
multispeciesPP_wrapper(pa_data = t1_pa_bird,
                       po_data = t1_po_bird,
                       bg = t1_bg,
                       species_names = final_species_set_bird,
                       habitat_associations = habitat_associations_bird,
                       group = "bird", predictor_set = "climate", out_name = "bird_multispecies_t1_climate")
### Habitat only
multispeciesPP_wrapper(pa_data = t1_pa_bird,
                       po_data = t1_po_bird,
                       bg = t1_bg,
                       species_names = final_species_set_bird,
                       habitat_associations = habitat_associations_bird,
                       group = "bird", predictor_set = "habitat", out_name = "bird_multispecies_t1_habitat")
#### Modern models (t2)
### Full predictor set
multispeciesPP_wrapper(pa_data = t2_pa_bird,
                       po_data = t2_po_bird,
                       bg = t2_bg,
                       species_names = final_species_set_bird,
                       habitat_associations = habitat_associations_bird,
                       group = "bird", predictor_set = "full", out_name = "bird_multispecies_t2_full")
### Climate only
multispeciesPP_wrapper(pa_data = t2_pa_bird,
                       po_data = t2_po_bird,
                       bg = t2_bg,
                       species_names = final_species_set_bird,
                       habitat_associations = habitat_associations_bird,
                       group = "bird", predictor_set = "climate", out_name = "bird_multispecies_t2_climate")
### Habitat only
multispeciesPP_wrapper(pa_data = t2_pa_bird,
                       po_data = t2_po_bird,
                       bg = t2_bg,
                       species_names = final_species_set_bird,
                       habitat_associations = habitat_associations_bird,
                       group = "bird", predictor_set = "habitat", out_name = "bird_multispecies_t2_habitat")
#### Single species models
#### Loop through each species
foreach (i = final_species_set_bird) %do% {
        ### Historic models (t1)
        ## Full predictor set
        multispeciesPP_wrapper(pa_data = t1_pa_bird,
                               po_data = t1_po_bird,
                               bg = t1_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_bird,
                               group = "bird", predictor_set = "full", out_name = paste("bird_", i, "_t1_full", sep = "")
                               )
        ## Climate only
        multispeciesPP_wrapper(pa_data = t1_pa_bird,
                               po_data = t1_po_bird,
                               bg = t1_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_bird,
                               group = "bird", predictor_set = "climate", out_name = paste("bird_", i, "_t1_climate", sep = "")
        )
        ## Habitat only
        multispeciesPP_wrapper(pa_data = t1_pa_bird,
                               po_data = t1_po_bird,
                               bg = t1_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_bird,
                               group = "bird", predictor_set = "habitat", out_name = paste("bird_", i, "_t1_habitat", sep = "")
        )
        ### Modern models (t2)
        ## Full predictor set
        multispeciesPP_wrapper(pa_data = t2_pa_bird,
                               po_data = t2_po_bird,
                               bg = t2_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_bird,
                               group = "bird", predictor_set = "full", out_name = paste("bird_", i, "_t2_full", sep = "")
        )
        ## Climate only
        multispeciesPP_wrapper(pa_data = t2_pa_bird,
                               po_data = t2_po_bird,
                               bg = t2_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_bird,
                               group = "bird", predictor_set = "climate", out_name = paste("bird_", i, "_t2_climate", sep = "")
        )
        ## Habitat only
        multispeciesPP_wrapper(pa_data = t2_pa_bird,
                               po_data = t2_po_bird,
                               bg = t2_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_bird,
                               group = "bird", predictor_set = "habitat", out_name = paste("bird_", i, "_t2_habitat", sep = "")
        )
}

##### -- Mammals -- #####
##### Multispecies models
#### Historic models (t1)
### Full predictor set
multispeciesPP_wrapper(pa_data = t1_pa_mamm,
                       po_data = t1_po_mamm,
                       bg = t1_bg,
                       species_names = final_species_set_mamm,
                       habitat_associations = habitat_associations_mamm,
                       group = "mamm", predictor_set = "full", out_name = "mamm_multispecies_t1_full")
### Climate only
multispeciesPP_wrapper(pa_data = t1_pa_mamm,
                       po_data = t1_po_mamm,
                       bg = t1_bg,
                       species_names = final_species_set_mamm,
                       habitat_associations = habitat_associations_mamm,
                       group = "mamm", predictor_set = "climate", out_name = "mamm_multispecies_t1_climate")
### Habitat only
multispeciesPP_wrapper(pa_data = t1_pa_mamm,
                       po_data = t1_po_mamm,
                       bg = t1_bg,
                       species_names = final_species_set_mamm,
                       habitat_associations = habitat_associations_mamm,
                       group = "mamm", predictor_set = "habitat", out_name = "mamm_multispecies_t1_habitat")
#### Modern models (t2)
### Full predictor set
multispeciesPP_wrapper(pa_data = t2_pa_mamm,
                       po_data = t2_po_mamm,
                       bg = t2_bg,
                       species_names = final_species_set_mamm,
                       habitat_associations = habitat_associations_mamm,
                       group = "mamm", predictor_set = "full", out_name = "mamm_multispecies_t2_full")
### Climate only
multispeciesPP_wrapper(pa_data = t2_pa_mamm,
                       po_data = t2_po_mamm,
                       bg = t2_bg,
                       species_names = final_species_set_mamm,
                       habitat_associations = habitat_associations_mamm,
                       group = "mamm", predictor_set = "climate", out_name = "mamm_multispecies_t2_climate")
### Habitat only
multispeciesPP_wrapper(pa_data = t2_pa_mamm,
                       po_data = t2_po_mamm,
                       bg = t2_bg,
                       species_names = final_species_set_mamm,
                       habitat_associations = habitat_associations_mamm,
                       group = "mamm", predictor_set = "habitat", out_name = "mamm_multispecies_t2_habitat")
#### Single species models
#### Loop through each species
foreach (i = final_species_set_mamm) %do% {
        ### Historic models (t1)
        ## Full predictor set
        multispeciesPP_wrapper(pa_data = t1_pa_mamm,
                               po_data = t1_po_mamm,
                               bg = t1_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_mamm,
                               group = "mamm", predictor_set = "full", out_name = paste("mamm_", i, "_t1_full", sep = "")
        )
        ## Climate only
        multispeciesPP_wrapper(pa_data = t1_pa_mamm,
                               po_data = t1_po_mamm,
                               bg = t1_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_mamm,
                               group = "mamm", predictor_set = "climate", out_name = paste("mamm_", i, "_t1_climate", sep = "")
        )
        ## Habitat only
        multispeciesPP_wrapper(pa_data = t1_pa_mamm,
                               po_data = t1_po_mamm,
                               bg = t1_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_mamm,
                               group = "mamm", predictor_set = "habitat", out_name = paste("mamm_", i, "_t1_habitat", sep = "")
        )
        ### Modern models (t2)
        ## Full predictor set
        multispeciesPP_wrapper(pa_data = t2_pa_mamm,
                               po_data = t2_po_mamm,
                               bg = t2_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_mamm,
                               group = "mamm", predictor_set = "full", out_name = paste("mamm_", i, "_t2_full", sep = "")
        )
        ## Climate only
        multispeciesPP_wrapper(pa_data = t2_pa_mamm,
                               po_data = t2_po_mamm,
                               bg = t2_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_mamm,
                               group = "mamm", predictor_set = "climate", out_name = paste("mamm_", i, "_t2_climate", sep = "")
        )
        ## Habitat only
        multispeciesPP_wrapper(pa_data = t2_pa_mamm,
                               po_data = t2_po_mamm,
                               bg = t2_bg,
                               species_names = i,
                               habitat_associations = habitat_associations_mamm,
                               group = "mamm", predictor_set = "habitat", out_name = paste("mamm_", i, "_t2_habitat", sep = "")
        )
}