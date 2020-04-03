#patient_setup.R
#input: none
#output: none
#side-effects: Sets up global variables _____
#optimizations: set this up so that it can accept a string or an excel file with this information included

patient_setup <- function(){
  
  HEALTHY_CONTROLS <<- 
    c(
      "bm5721", 
      "bm5871", 
      "bm6152", 
      "bm6631"
    )
  
  PAIRED_PATIENTS <<- 
    c(
      "parbiu", 
      "parwxu", 
      "parzuu", 
      "pasbee", 
      "pasgwh", 
      "patjsp", 
      "pasrtp", 
      "pasvzc", 
      "pasxvc", 
      "patgiy", 
      "pathiw", 
      "patlhb"
    )
  
  UNPAIRED_PATIENTS <<- 
    c(
      "papwhs", 
      "parajx", 
      "parant", 
      "parbfj", 
      "parcvp", 
      "parklc", 
      "parmme", 
      "parpwl", 
      "partxh", 
      "parxmp", 
      "parzwh", 
      "paseyc", 
      "paself", 
      "pastjj", 
      "paspsv", 
      "pasptm", 
      "pasxjy", 
      "pasyyw"
    )
}
