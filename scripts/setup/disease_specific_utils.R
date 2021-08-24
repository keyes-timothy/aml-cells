
#' Find a low-dimensional subspace representing an approximate basis for 
#' healthy cytometry data
#'
#' @param healthy_data A tibble in which rows represent cells and columns 
#' represent single-cell measurements (i.e. CyTOF antigens). Only cells taken
#' from healthy samples should be included. 
#' 
#' @param cluster_col An unquoted column name containing the cluster ids for each
#' cell in `healthy_data`, if applicable. This argument should remain NULL if the 
#' healthy basis should be computed using the healthy cells in bulk (i.e. not
#' by subtype).
#' 
#' @param pca_threshold a value between 0 and 1 representing the percentage of 
#' variance in the original data that should be preserved during PCA transformation.
#' This will determine the number of PCs in each `healthy_pcs` result matrix 
#' (see below)
#'
#' @return `healthy_pcs`, an [m x n] principal component matrix in which m
#' represents the number of channels in `healthy_data` and n represents the 
#' number of principal components ("eigencells") forming a basis for the healthy
#' subspace. If `cluster_col` is not NULL, the returned result will be a tibble 
#' with two columns (`cluster_col` and `healthy_pcs`) such that each row will 
#' represent a reduced-dimensional basis for a given cluster's healthy phenotype
#' space. 
#' 
#'
aml_find_healthy_subspace <-
  function(
    healthy_data = NULL, 
    cluster_col = NULL, 
    pca_threshold = 0.95
  ) {
    
    # group healthy_data if there is a cluster variable
    healthy_data <- 
      healthy_data %>% 
      group_by({{cluster_col}})
  
    # keep track of the antigen names in the input tibble
    col_names <- 
      healthy_data %>% 
      ungroup() %>% 
      select(-{{cluster_col}}) %>% 
      colnames()
    
    healthy_pcs <- 
      healthy_data %>%
      # better logic for the following line (to avoid the warning) is needed
      nest() %>% 
      transmute(
        data = map(.x = data, .f = t), 
        pcs = 
          map(
            .x = data, 
            .f = 
              ~recipe(.x) %>% 
              step_pca(all_numeric(), threshold = pca_threshold) %>% 
              prep() %>% 
              juice() %>% 
              as.matrix()
          )
      )
    
    if (nrow(healthy_pcs) == 1) { 
      healthy_pcs <- 
        healthy_pcs$pcs[[1]]
      row.names(healthy_pcs) <- col_names
    } else {
    healthy_pcs <- 
      healthy_pcs %>% 
      mutate(pcs = map(pcs, function(x) {row.names(x) <- col_names; return(x)}))
    }
    
    return(healthy_pcs)
  }

# will return the cell_matrix such that the column names correspond to each 
# cell's cell_id in the original input tibble.

cell_matrix_transformation <- function(cell_matrix) {
  result <- 
    cell_matrix %>% 
    t()
  
  colnames(result) <- result[1,]
  result <- result[2:nrow(result),]
  
  return(result)
}

#project_cell_matrix <- function(cell_matrix, healthy_pcs)

#' Project single-cell measurements onto a PCA-defined healthy subspace
#'
#' @param aml_data A tibble in which rows represent cells and columns represent 
#' antigens measured for each cell. Only the columns in `aml_data` corresponding 
#' to antigens used to build the healthy subspace will be used in the projection. 
#' All others will be ignored. 
#' 
#' @param healthy_pcs A matrix output of the function `aml_find_healthy_subspace`
#'
#' @return A tibble with an equal number of rows as `aml_data` and twice the number 
#' of columns. Each row will represent a cell from `aml_data`, and each column 
#' will represent either the healthy or disease-specific component of the 
#' healthy subspace projection. 
#' 
#' @export
#'
#' @examples
aml_project_onto_healthy <- 
  function(
    aml_data,
    healthy_pcs, 
    cluster_col = NULL
  ) {
    # extract antigens from healthy_pcs input object 
    if ("matrix" %in% class(healthy_pcs)) { 
      antigens <- row.names(healthy_pcs)
    } else { 
      antigens <- row.names(healthy_pcs$pcs[[1]])
    }
    
    # check that aml_data and healthy_pcs are compatible
    if (!all(antigens %in% colnames(aml_data))) { 
      stop("the columns of aml_data do not contain all the antigens used 
           to create the healthy subspace")
    }
    
    # make a cell_matrix for each cluster (if there are any) 
    # in which each row is an antigen from `antigens` and each column 
    # is a cell from `aml_data`. 
    cell_matrices <- 
      aml_data %>% 
      # give each cell an id so that the output tibble has rows in the 
      # same order as the input tibble
      mutate(cell_id = 1:nrow(aml_data)) %>% 
      select(cell_id, all_of(antigens), {{cluster_col}}) %>% 
      group_by({{cluster_col}}) %>% 
      nest() %>% 
      rename(cell_matrix = data) %>% 
      mutate(cell_matrix = map(cell_matrix, cell_matrix_transformation))
    
    # combine cell_matrices and healthy_pcs
    
    ## if there are no clusters
    if (nrow(cell_matrices) == 1L) { 
      projection_tibble <- 
        tibble(
          cell_matrix = list(cell_matrices$cell_matrix[[1]]), 
          pcs = list(healthy_pcs)
        )
    } else {
      
      ## if there are clusters
      projection_tibble <- 
        cell_matrices %>% 
        left_join(healthy_pcs)
    }
    
    # fit one linear model for each cell in the cell_matrix
    # without using an intercept! 
    projection_tibble <- 
      projection_tibble %>% 
      mutate(
        lm_result = 
          map2(.x = cell_matrix, .y = pcs, .f = function(.x, .y) lm(.x ~ .y + 0))
      ) %>% 
      transmute(
        healthy_component = 
          map(
            .x = lm_result, 
            .f = ~ 
              .x %>% 
              pluck("fitted.values") %>% 
              t() %>% 
              as_tibble() %>% 
              rename_with(.fn = ~str_c(.x, "_healthy")),
          ), 
        disease_component = 
          map(
            .x = lm_result, 
            .f = ~ 
              .x %>% 
              pluck("residuals") %>% 
              t() %>% 
              as_tibble() %>% 
              rename_with(.fn = ~str_c(.x, "_disease"))
          ), 
        cell_ids = 
          map(.x = cell_matrix, .f = ~tibble(cell_id = colnames(.x)))
      ) %>% 
      unnest(cols = c(cell_ids, healthy_component, disease_component)) %>% 
      ungroup() %>% 
      # put rows back in original order
      mutate(cell_id = as.numeric(cell_id)) %>% 
      arrange(cell_id) %>% 
      select(-cell_id, -{{cluster_col}})
    
    
    return(projection_tibble)
  }


