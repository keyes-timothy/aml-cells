#marker_setup.R
#input: none
#output: none
#side-effects: Sets up global variables SURFACE.MARKERS, SIGNALING.MARKERS, TRX.FACTORS, and OTHER.MARKERS
#optimizations: set this up so that it can accept a string or an excel file with this information included

marker_setup <- function(){
  if(is.AML){
    SURFACE.MARKERS <<- c("CD45", "CD99", "CD45RA", "CD93", "CCR2", "CD90", "CD38", "CD34" 
                         , "CD64", "CD117", "CD123", "CD11c", "CD13", "TIM-3", "CD56"
                         , "CD10", "CD33", "PD-L1", "CD14", "CD41", "CD16", "CD68"
                         , "CD47", "CD135", "CD109", "CD49f", "HLA-DR", "CD71", "CD11b", 
                         "CD61")
    SIGNALING.MARKERS <<- c("pSTAT5", "pSTAT3", "pS6", "pCreb", "pErk", "pAkt")
    TRX.FACTORS <<- c("PU.1", "CEBPa", "GATA-1")
    OTHER.MARKERS <<- c("caspase-3", "MPO")
    
    ALL.MARKERS <<- c(SURFACE.MARKERS, SIGNALING.MARKERS, TRX.FACTORS, OTHER.MARKERS)
    
    CLASSIFIER.MARKERS <<- c('CD45', 'CD34', 'CD38', 'CD61', 'CD14',
                            'CD135', 'CD45RA', 'CD90', 'HLA-DR', 'CD41', 
                            'CD13', 'CD11b', 'CD11c')
    CLASSIFIER.POPULATIONS <<- c('HSC', 'MPP', 'CMP', 'GMP', 'MEP', 
                                'Monocyte', 'DC', 'Macrophage', 'Thrombocyte')
    
  } else {
      SURFACE.MARKERS <<- c()
      SIGNALING.MARKERS <<- c()
      TRX.FACTORS <<- c()
      OTHER.MARKERS <- c()
      
      ALL.MARKERS <<- c(SURFACE.MARKERS, SIGNALING.MARKERS, TRX.FACTORS, OTHER.MARKERS)
      
      CLASSIFIER.MARKERS <<- tolower(c('CD19', 'CD20', 
                                      'CD34', 'CD38', 
                                      'Igmi', 'Igms', 
                                      'CD179a', 'CD179b',
                                      'CD127', 'tdt',
                                       'CD45'))
      CLASSIFIER.POPULATIONS <<- c('HSC', 'Progenitor I', 'Progenitor II', 'Progenitor III', 
                                  'Pre-Pro-B', 
                                  'Pro B-I', 'Pro B-II', 
                                  'Pre B-I', 'Pre B-II', 
                                  'Immature B-I', 'Immature B-II', 
                                  'Mature B-I', 'Mature B-II', 
                                  'Early non B-I', 'Early non B-II', 
                                  'Mature non-B')
  }
}
