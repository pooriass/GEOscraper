#### Example 





term <- c('(colorectal cancer)')



ids <- FindGSE(Search_Term = term, Organism = 'Homo sapiens', Study_type = 'Expression profiling by high throughput sequencing', max_result = 2000)


perioud <- seq(1, length(ids), 100)
results_list <- list()

for ( i in 1:length(perioud)){
  
  id <- perioud[i]:(perioud[i]+99)
results <- GEOscrap(GSE = ids[id], attr = c('single cell','single-cell','drug resistance', 'drug', '5FU', 'oxaliplatin', 'xenograft') ) 
results_list[[i]] <- results
}

for ( i in 7:length(perioud)){
  
  id <- perioud[i]:(perioud[i]+99)
  results <- GEOscrap(GSE = ids[id], attr = c('single cell','single-cell','drug resistance', 'drug', '5FU', 'oxaliplatin', 'xenograft') ) 
  results_list[[i]] <- results
}


for ( i in 11:length(perioud)){
  
  id <- perioud[i]:(perioud[i]+99)
  results <- GEOscrap(GSE = ids[id], attr = c('single cell','single-cell','drug resistance', 'drug', '5FU', 'oxaliplatin', 'xenograft') ) 
  results_list[[i]] <- results
}


final_result  <- rbind(data.frame(), results_list[[1]])
for (i in 1:length(results_list)){
  
final_result  <- rbind(final_result, results_list[[i]])
  
}  


final_result <- separate(data = final_result,col = 'Summary', sep = 'AND' , into = c('single cell','single-cell','drug resistance', 'drug', '5FU', 'oxaliplatin', 'xenograft'))



final_result <- separate(data = final_result,col = 'overall_design', sep = 'AND' , into = c('single cell2','single-cell2','drug resistance2', 'drug2', '5FU2', 'oxaliplatin2', 'xenograft2'))

final_result <- final_result[!duplicated(final_result$GSE_ID),]

final_result <- final_result[!is.na(final_result$GSE_ID),]


final_result <- final_result2


for ( i in 1:nrow(final_result)){
  
  # Single cell 
  if (any(c(final_result$`single cell`[i], final_result$`single-cell`[i], 
        final_result$`single cell2`[i], final_result$`single-cell2`[i]) %in% 'single cell ') ){
    final_result$`single cell`[i] <- 'Found'
  }else{final_result$`single cell`[i] <- 'Not Found'}

  
  # Drug resistance 
  if (any(c(final_result$`drug resistance`[i], final_result$`drug resistance2`[i]) %in% c('drug resistance' )  )){
    final_result$`drug resistance`[i] <- 'Found'
  }else{
    final_result$`drug resistance`[i]<- 'Not found'}

  
  # Drug 
  if (any(c(final_result$drug [i], final_result$drug2 [i]) %in% 'drug ')  ){
    final_result$drug [i] <- 'Found' 
  }else{final_result$drug [i] <- 'Not Found'}
  

  
  # 5FU
  if (any(c(final_result$`5FU`[i], final_result$`5FU2`[i]) %in% '5FU ')) {
    final_result$`5FU`[i] <- 'Found'
  } else {
    final_result$`5FU`[i] <- 'Not Found'
  }
  # Oxaliplatin
  if (any(c(final_result$oxaliplatin [i], final_result$oxaliplatin2 [i]) %in% 'oxaliplatin ' ) ){
    final_result$oxaliplatin [i] <- 'Found' 
  }else{final_result$oxaliplatin [i] <- 'Not Found'}
  
  
  
  # xenograft
  if (any(c(final_result$xenograft [i], final_result$xenograft2 [i]) %in% 'xenograft' ) ){
    final_result$xenograft [i] <- 'Found' 
  }else{final_result$xenograft [i] <- 'Not Found'}
  

  
  
  
}



final_result <- final_result[, !colnames(final_result) %in% c('single-cell', 'single-cell2', 'single cell2', 'drug resistance2', 'drug2', '5FU2', 'oxaliplatin2', 'xenograft2')]



single_title <-c( grep('single cell', final_result$title), grep('single-cell', final_result$title),
                  grep('scRNA-seq', final_result$title) , grep('scRNA seq', final_result$title),
                  grep('ScRNA-seq', final_result$title) , grep('ScRNA seq', final_result$title),
                  grep('Single cell', final_result$title), grep('Single-cell', final_result$title))




final_result$`single cell` [single_title]<- 'Found' 

write.csv(final_result,'D:/Research/Fetching GEO Data/final_result.csv')

write.csv(final_result[which(final_result$`single cell` == 'Found'),],'D:/Research/Fetching GEO Data/final_result_scrnaseq.csv')





