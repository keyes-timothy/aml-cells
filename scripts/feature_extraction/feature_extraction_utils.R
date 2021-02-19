# feature_extraction_utils.R

# tof_tibble must be a tibble of single-cell observations
# cluster_col is the UQ name of the column defining the clusters
extract_feature_proportion <- 
  function(tof_tibble, cluster_col) { 
    
    abundances <- 
      tof_tibble %>%
      group_by({{cluster_col}}, .add = TRUE) %>% 
      summarize(abundance = n()) %>% 
      # two lines below can be compressed into a transmute if multidplyr is not being used
      mutate(
        prop = abundance / sum(abundance), 
        "{{cluster_col}}" := str_c("prop", {{cluster_col}}, sep = "@")
      ) %>% 
      select(prop, {{cluster_col}}) %>% 
      as_tibble() %>% 
      pivot_wider(
        names_from = {{cluster_col}},  
        values_from = prop,
        # note that if a cluster is not present for a given group, it will 
        # be filled in as having a relative abundance of 0. 
        values_fill = 0 
      )
    
    return(abundances)
  }

# tof_tibble must be a tibble of single-cell observations
# cluster_col is the UQ name of the column defining the clusters
#
# dtplyr version of the above function
extract_feature_proportion_dt <- 
  function(tof_tibble, cluster_col) { 
    
    existing_groups <- group_vars(tof_tibble)
    
    tof_tibble <- 
      tof_tibble %>% 
      ungroup() %>% 
      unite(col = "existing_groups", all_of(existing_groups), sep = "_@@_") %>% 
      rename("my_cluster" := {{cluster_col}})
    
    abundances <- 
      tof_tibble %>%
      lazy_dt() %>% 
      group_by(existing_groups, my_cluster) %>% 
      summarize(num_cells = n()) %>%
      mutate(abundance = num_cells / sum(num_cells)) %>% 
      select(-num_cells) %>% 
      as_tibble() %>% 
      separate(col = existing_groups, into = existing_groups, sep = "_@@_") %>%
      pivot_wider(
        names_from = my_cluster,
        values_from = abundance,
        names_prefix = "prop@",
        # note that if a cluster is not present for a given group, it will
        # be filled in as having a relative abundance of 0.
        values_fill = 0
      )
    
    return(abundances)
  }


#
#
#
extract_feature_ct <- 
  function(
    tof_tibble,
    cluster_col,
    marker_cols = where(is.numeric), 
    central_tendency = mean
  ) { 
    central_tendencies <- 
      tof_tibble %>% 
      group_by({{cluster_col}}, .add = TRUE) %>% 
      summarize(across({{marker_cols}}, central_tendency)) %>% 
      pivot_longer(
        cols = {{marker_cols}}, 
        names_to = "channel", 
        values_to = "value"
      ) %>%
      pivot_wider( 
        names_from = c(channel, {{cluster_col}}), 
        values_from = value, 
        names_sep = "@"
      )
    
    central_tendencies
  }




# 
# dtplyr version of above function
extract_feature_ct_dt <-
  function(
    tof_tibble,
    cluster_col,
    marker_cols = everything(), 
    central_tendency = mean
  ) { 
    
    # check if the incoming tibble is grouped
    existing_groups <- group_vars(tof_tibble)
    
    # only unite if there are existing groups 
    if (length(existing_groups) > 0) { 
      tof_tibble <- 
        tof_tibble %>%
        unite(col = "existing_groups", all_of(existing_groups), sep = "ϴ")
    }
    
    tof_tibble <-
      tof_tibble %>%
      rename(my_cluster = {{cluster_col}}) %>% 
      select(any_of("existing_groups"), my_cluster, {{marker_cols}})
    
    central_tendencies <-
      tof_tibble %>%
      lazy_dt() %>%
      group_by(my_cluster, existing_groups) %>%
      summarize(across(everything(), central_tendency)) %>%
      as_tibble() 
    
    if (length(existing_groups) > 0) { 
      central_tendencies <- 
        central_tendencies %>%
        separate(
          col = existing_groups,
          into = existing_groups,
          sep = "ϴ"
        )
    }
    
    central_tendencies <- 
      central_tendencies %>%
      pivot_longer(
        cols = {{marker_cols}},
        names_to = "channel",
        values_to = "value"
      ) %>%
      pivot_wider(
        names_from = c(channel, my_cluster),
        values_from = value,
        names_sep = "@"
      )
    
    central_tendencies
  }


#
#
# note: can change this to include if you want different thresholds for each channel
extract_feature_signaling_threshold <- 
  function(
    tof_tibble, 
    cluster_col, 
    stimulation_col, 
    signaling_cols = where(is.numeric), 
    threshold = asinh(10 / 5)
  ) {
    # store what the existing groups of the input tof_tibble are 
    existing_groups <- group_vars(tof_tibble)
    
    signaling_features <- 
      tof_tibble %>% 
      group_by({{cluster_col}}, {{stimulation_col}}, .add = TRUE) %>% 
      summarize(across({{signaling_cols}}, ~ mean(.x > threshold))) %>% 
      pivot_longer(
        cols = {{signaling_cols}},
        names_to = "channel", 
        values_to = "values"
      ) %>% 
      # get rid of the cluster_col grouping and replace with only the groups
      # from the input tof_tibble
      group_by(across(all_of(existing_groups))) %>% 
      transmute(
        col_names = str_c(stimulation, "_", channel, "@", {{cluster_col}}),
        values
      ) %>%
      pivot_wider(
        names_from = col_names,
        values_from = values
      )
    
    signaling_features
  } 



# first, a quick helper function to help find the kernel density estimates for each single-cell distribution

tof_find_density <- function(my_cells, ...) {
  my_cells <- 
    my_cells %>% 
    purrr::pluck("values")
  if (length(my_cells) < 2) { 
    tibble()
  } else {
    density <- density(my_cells, ...)
    tibble(value = density$x, weight = density$y)
  }
}


# Another helper functions for finding the EMD

tof_find_emd <- function(basal_densities, other_densities) {
  if (identical(basal_densities, tibble())) {
    emd <- NA_real_
  } else if (identical(other_densities, tibble())) {
    emd <- NA_real_
  } else {
    emd <- 
      quietly(emdist::emdw)(
        A = basal_densities$value, 
        wA = basal_densities$weight, 
        B = other_densities$value, 
        wB = other_densities$weight
      ) %>% 
      pluck("result")
  }
  emd
}


#
#
#
#
extract_feature_signaling_emd <- 
  function(
    tof_tibble,
    cluster_col, 
    stimulation_col,
    signaling_cols = where(is.numeric), 
    basal_level # a string indicating what the basal condition is called in stimulation_col
  ) {
    
    #existing_groups <- group_vars(tof_tibble)
    
    all_densities <- 
      tof_tibble %>%
      #mutate(across(where(is.numeric), scale)) %>% 
      pivot_longer(
        cols = {{signaling_cols}},
        names_to = "channel",
        values_to = "values"
      ) %>%
      group_by(channel, {{cluster_col}}, {{stimulation_col}}, .add = TRUE) %>%
      nest() %>%
      mutate(densities = map(.x = data, .f = tof_find_density)) %>%
      select(-data) %>%
      ungroup()
    
    basal_densities <- 
      all_densities %>% 
      filter({{stimulation_col}} == basal_level) %>% 
      select(
        -{{stimulation_col}},
        basal_densities = densities
      ) %>% 
      ungroup()
    
    stimulation_densities <- 
      all_densities %>% 
      filter({{stimulation_col}} != basal_level) %>% 
      ungroup()
    
    emd_features <- 
      stimulation_densities %>%
      quietly(left_join)(basal_densities) %>%
      pluck("result") %>%
      mutate(
        basal_densities = 
          map(basal_densities, ~ if (is.null(.x)) {tibble()} else {.x}),
        emd = map2_dbl(.x = basal_densities, .y = densities, .f = tof_find_emd), 
        col_names = str_c({{stimulation_col}}, "_", channel, "@", {{cluster_col}})
      ) %>%
      select(
        -basal_densities, 
        -densities, 
        -{{stimulation_col}}, 
        -{{cluster_col}}, 
        -channel
      ) %>% 
      pivot_wider(
        names_from = col_names, 
        values_from = emd
      )
    
    emd_features
  }



#
#
#
#
#
#first, a helper function to find the jensen-shannon index 
tof_find_js <- function(basal_densities, other_densities) {
  if (identical(basal_densities, tibble())) {
    js_index <- NA_real_
  } else if (identical(other_densities, tibble())) {
    js_index <- NA_real_
  } else {
    
    # bin the kernel density estimates by rounding to the 2nd decimal place
    basal_densities <- 
      mutate(basal_densities, value = round(value, 2))
    other_densities <- 
      mutate(other_densities, value = round(value, 2))
    
    # there might be a smarter way to compute this 
    js_index <- 
      basal_densities %>% 
      full_join(other_densities, by = "value") %>% 
      drop_na() %>% 
      select(-value) %>% 
      as.matrix() %>% 
      t() %>% 
      quietly(philentropy::JSD)(test.na = FALSE) %>% 
      pluck("result")
  }
  
  js_index
}



#
#
#
#
extract_feature_signaling_js <- 
  function(
    tof_tibble,
    cluster_col, 
    stimulation_col,
    signaling_cols = where(is.numeric), 
    basal_level # a string indicating what the basal condition is called in stimulation_col
  ) {
    
    all_densities <- 
      tof_tibble %>%
      select({{cluster_col}}, {{stimulation_col}}, {{signaling_cols}}) %>% 
      mutate(across(where(is.numeric), ~ as.numeric(scale(.x)))) %>% 
      pivot_longer(
        cols = {{signaling_cols}},
        names_to = "channel",
        values_to = "values"
      ) %>%
      group_by(channel, {{cluster_col}}, .add = TRUE) %>% 
      mutate(
        max = max(values) %>% ceiling(), 
        min = min(values) %>% floor()
      ) %>% 
      group_by({{stimulation_col}}, .add = TRUE) %>%
      nest(data = values) %>%
      mutate(
        densities =
          pmap(
            .l = list(data, min, max),
            ~ tof_find_density(my_cells = ..1, from = ..2, to = ..3)
          )
      ) %>%
      select(-data, -max, -min) %>%
      ungroup()
    
    basal_densities <- 
      all_densities %>%
      filter({{stimulation_col}} == basal_level) %>%
      select(
        -{{stimulation_col}},
        basal_densities = densities
      ) %>%
      ungroup()
    
    stimulation_densities <- 
      all_densities %>%
      filter({{stimulation_col}} != basal_level) %>%
      ungroup()
    
    js_features <- 
      stimulation_densities %>%
      quietly(left_join)(basal_densities) %>%
      pluck("result") %>%
      mutate(
        basal_densities =
          map(basal_densities, ~ if (is.null(.x)) {tibble()} else {.x}),
        js = map2_dbl(.x = basal_densities, .y = densities, .f = tof_find_js),
        col_names = str_c({{stimulation_col}}, "_", channel, "@", {{cluster_col}})
      ) %>%
      select(
        -basal_densities,
        -densities,
        -{{stimulation_col}},
        -{{cluster_col}},
        -channel
      ) %>%
      pivot_wider(
        names_from = col_names,
        values_from = js
      )
    
    js_features
  }


#
#
#
#
extract_feature_signaling <- 
  function(
    tof_tibble, 
    method = c("threshold", "js", "emd"), 
    ... # arguments to pass to extract_feature_signaling_* variant
  ) { 
    # check that the method is valid
    rlang::arg_match(method)
    
    if (method == "threshold") { 
      # compute the threshold signaling features 
      signaling_features <- 
        tof_tibble %>% 
        extract_feature_signaling_threshold(...)
      
    } else if (method == "emd") { 
      #compute the emd signaling features
      signaling_features <- 
        tof_tibble %>% 
        extract_feature_signaling_emd(...)
      
    } else { 
      # compute the js signaling features
      signaling_features <- 
        tof_tibble %>% 
        extract_feature_signaling_js(...)
    } 
  } 



#
#
#
#
tof_extract_ddpr <- 
  function(
    healthy_tibble, # tibble with all the healthy cells
    cancer_tibble, # tibble with all the cancer cells
    cluster_col, # what is the cluster column called
    stimulation_col, # what is the stimulation column called
    lineage_cols, # which are lineage markers
    signaling_cols, # which are signaling markers
    threshold = asinh(10 / 5),
    central_tendency = mean,
    only_expanded_clusters = FALSE # whether or not you only want to include the clusters expanded in the cancer samples relative to the healthy ones
  ) {
    # store existing groups to preserve them in the output
    existing_groups <- group_vars(cancer_tibble)
    
    if (only_expanded_clusters) {
      
      #find healthy abundances
      healthy_abundances <- 
        healthy_tibble %>%
        extract_feature_proportion(cluster_col = {{cluster_col}}) %>% 
        pivot_longer(
          cols = contains("prop@"), 
          names_to = "cluster", 
          values_to = "healthy_proportion", 
          names_prefix = "prop@"
        )
      # find cancer abundances
      cancer_abundances <- 
        cancer_tibble %>% 
        extract_feature_proportion(cluster_col = {{cluster_col}}) %>% 
        pivot_longer(
          cols = contains("prop@"), 
          names_to = "cluster", 
          values_to = "cancer_proportion", 
          names_prefix = "prop@"
        )
      expanded_clusters <- 
        full_join(healthy_abundances, cancer_abundances)
      # keep only the clusters for which there is at least a 
      # 5% expansion in the cancer samples overall
      filter(
        (cancer_proportion - healthy_proportion) / healthy_proportion > 0.05
      ) %>% 
        pull(cluster)
      
    } else {
      expanded_clusters <- 
        cancer_tibble %>% 
        pull({{cluster_col}}) %>% 
        unique()
    }
    
    # find abundance features
    abundance_features <- 
      cancer_tibble %>%
      filter({{cluster_col}} %in% expanded_clusters) %>% 
      extract_feature_proportion(cluster_col = {{cluster_col}})
    
    # find lineage features
    lineage_features <- 
      cancer_tibble %>% 
      filter({{cluster_col}} %in% expanded_clusters) %>% 
      select(all_of(existing_groups), {{lineage_cols}}, {{cluster_col}}) %>% 
      extract_feature_ct(
        cluster_col = {{cluster_col}}, 
        central_tendency = central_tendency
      )
    
    signaling_features <- 
      cancer_tibble %>% 
      filter({{cluster_col}} %in% expanded_clusters) %>% 
      select(
        all_of(existing_groups), 
        {{signaling_cols}}, 
        {{cluster_col}}, 
        {{stimulation_col}}
      ) %>% 
      extract_feature_signaling_threshold(
        cluster_col = {{cluster_col}},
        stimulation_col = {{stimulation_col}}, 
        signaling_cols = {{signaling_cols}}, 
        threshold = threshold
      )
    
    # put everything together 
    if (length(existing_groups) == 0) { 
      final_features <- 
        bind_cols(abundance_features, lineage_features, signaling_features)
    } else {
      final_features <- 
        abundance_features %>% 
        left_join(lineage_features) %>% 
        left_join(signaling_features)
    }
    
    final_features
  }


#
#
#
#
#
tof_extract_citrus <- 
  function(
    tof_tibble, 
    cluster_col, 
    stimulation_col, 
    marker_cols
  ) {
    citrus_tibble <- 
      tof_tibble %>% 
      group_by({{stimulation_col}})
      extract_feature_ct(
        cluster_col = {{cluster_col}}, 
        central_tendency = median
      )
    
    citrus_tibble
    
  }
    





