#marker_setup.R
#input: none
#output: none
#side-effects: Sets up global variables SURFACE_MARKERS, SIGNALING_MARKERS, TRX_FACTORS, and OTHER_MARKERS
#optimizations: set this up so that it can accept a string or an excel file with this information included

marker_setup <- function(){
  
  SURFACE_MARKERS <<- 
    c(
      "CD3_CD19",
      "CD10",
      "CD11b",
      "CD11c",
      "CD13",
      "CD14",
      "CD16",
      "CD33",
      "CD34",
      "CD38",
      "CD41",
      "CD45",
      "CD45RA",
      "CD47",
      "CD49f",
      "CD56",
      "CD61",
      "CD64",
      "CD68",
      "CD71",
      "CD90",
      "CD93",
      "CD99", 
      "CD109",
      "CD117", 
      "CD123",
      "CD135", 
      "CCR2", 
      "TIM-3", 
      "PD-L1", 
      "HLA-DR"
    )
  
  SIGNALING_MARKERS <<- c("pAkt", "pCreb", "pErk", "pS6", "pSTAT3", "pSTAT5")
  
  TRX_FACTORS <<- c("CEBPa", "GATA-1", "PU.1")
  
  OTHER_MARKERS <<- c("caspase-3", "MPO")
  
  ALL_MARKERS <<- 
    c(
      SURFACE_MARKERS, 
      SIGNALING_MARKERS, 
      TRX_FACTORS, 
      OTHER_MARKERS
    )
  
  CLASSIFIER.MARKERS <<- 
    c(
      'CD45', 
      'CD34', 
      'CD38', 
      'CD61', 
      'CD14',
      'CD135', 
      'CD45RA', 
      'CD90', 
      'HLA-DR', 
      'CD41', 
      'CD13', 
      'CD11b', 
      'CD11c'
    )
  
  CLASSIFIER.POPULATIONS <<- 
    c(
      'HSC', 
      'MPP', 
      'CMP', 
      'GMP', 
      'MEP', 
      'Monocyte', 
      'DC', 
      'Macrophage',
      'Thrombocyte'
    )
}
