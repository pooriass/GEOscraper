


#install.packages('rvest')
#install.packages('httr')
#install.packages('xml2')
#install.packages('retry')
#install.packages('tidyverse')
#install.packages('stringr')



library(rvest)
library(httr)
library(xml2)
library(retry)
library(tidyverse)
library(stringr)



################################################## Find GSE #######################################


FindGSE <- function(Search_Term = Search_Term, Organism = organism, Study_type = Study_type, max_result=max_result){
  
######################################## Initial Set up ###########################################
  user_agents <- c(
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.113 Safari/537.36",
    "Mozilla/5.0 (Windows NT 6.1; WOW64; Trident/7.0; AS; rv:11.0) like Gecko",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_6) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15", 
    'Mozilla/5.0 (iPhone; CPU iPhone OS 14_5 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.0 Mobile/15E148 Safari/604.1',
    'Mozilla/5.0 (Linux; Android 10; SM-G973F) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Mobile Safari/537.36', 
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.0.3 Safari/605.1.15',
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Safari/537.36', 
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:88.0) Gecko/20100101 Firefox/88.0',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.818.42 Safari/537.36 Edg/90.0.818.42',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0', 
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Safari/537.36'
    
    
  )
  
  fetch_with_backoff <- function(url, max_attempts = 5, initial_wait = 1) {
    attempt <- 1
    success <- FALSE
    response <- NULL
    
    while (attempt <= max_attempts && !success) {
      user_agent <- sample(user_agents, 1)
      
      response <- GET(url, user_agent(user_agent))
      
      if (status_code(response) == 429) {
        wait_time <- initial_wait * (2 ^ (attempt - 1))  # Exponential backoff
        message(paste("Received 429, waiting", wait_time, "seconds before retrying..."))
        Sys.sleep(wait_time)
        attempt <- attempt + 1
      } else {
        success <- TRUE
      }
    }
    
    if (!success) {
      stop("Failed to fetch the URL after ", max_attempts, " attempts.")
    }
    
    return(content(response, as = "text"))
  }
  
  ######################### Search Process ########################################
  Search_Term <- paste0(Search_Term, collapse = '')
  Search_Term <- gsub(' ','+',Search_Term)
  
  
  
  
  ################## Temp Pause ######################
  pause_time <- sample(1:2)[1]
  
  message(paste("Pausing for", pause_time,  "seconds...") )
  Sys.sleep(pause_time)
  message("Resuming execution...")
  
  
  ## Retrivieng the main search terms 
  url1 <- paste0('https://www.ncbi.nlm.nih.gov/gds/?term=', Search_Term) 
  content <- fetch_with_backoff(url = url1)
  parsed_page <- read_html(content)
  
   main_term <- parsed_page %>% html_node("body") %>% html_node("textarea") %>% html_text()
   
   
   ######## Completing the Main_term 
   
  # Organism 
   if(length(Organism ) ==0 ) { 
     main_term <- main_term 
   }else{
       if(length(Organism) == 1){
         main_term <- paste0(main_term,  'AND \"',Organism, '\"[porgn]')
       }else{
         main_term <-  paste0(main_term,  'AND \"',Organism[1],'\"[porgn]', 'AND \"' , Organism[2], '\"[porgn]')
       }
   }
   
   
   # Study_type 
   
   
   if(length(Study_type) ==0 ) { 
     main_term <- main_term 
   }else{
     if(length(Study_type== 1) ) {
       main_term <- paste0(main_term, 'AND \"',Study_type, '\"[Filter]')
     }else{
       if(length(Study_type==3 )) {
       main_term <-  paste0(main_term,  'AND \"',Study_type[1], '\"[Filter]' ,'AND \"' , Study_type[2], '\"[Filter]')
       }else{
         main_term <-  paste0(main_term,  'AND \"',Study_type[1],'\"[Filter]' ,'AND \"' , Study_type[2],'\"[Filter]' ,'AND \"' , Study_type[3], '\"[Filter]')
         
       }
       }
   }
     
  
   
  
   
   ################## Temp Pause ######################
   pause_time <- sample(1:2)[1]
   
   message(paste("Pausing for", pause_time,  "seconds...") )
   Sys.sleep(pause_time)
   message("Resuming execution...")
   
   
  ####### Retrivieng IDs  
   
  
  max_result1 <- paste0('&retmax=',max_result )
  url1 <-paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=', main_term, max_result1) 
  url1 <- gsub(' ', '%20', url1)
  content <- fetch_with_backoff(url = url1)
  parsed_page <- read_html(content)
  


  matches <- gregexpr("<id>(\\d+)</id>", as.character(parsed_page))
  
  # Extract the matches
  matches <- regmatches(as.character(parsed_page), matches)
  
  # Flatten the list of matches
  matches <- gsub('<id>', '', unlist(matches) )  
  matches <- gsub('</id>', '', unlist(matches) )  
  
  
  ################## Temp Pause ######################
  pause_time <- sample(1:2)[1]
  
  message(paste("Pausing for", pause_time,  "seconds...") )
  Sys.sleep(pause_time)
  message("Resuming execution...")
  
  
  ##################### IDs to GSE ######################
  
  ### making a hundred list 
  
  if(length(matches) > 100 ){
  
  number_of_id <- floor(length(matches)/100)
  
  gener_number <- list()
  for ( i in 1:number_of_id){
  gener_number[i] <- 100*(i - 1 ) + 1
}
  
  gener_number <- unlist(gener_number)
  
  

  #### Untill divided numbers 
  GSE_list <- list()
  for ( i in 1:length(gener_number)) {
    number_number <- gener_number[i]
    match_numb <- number_number:(number_number+99)
    

 url2 <-  paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=', paste0(matches[match_numb], collapse = ',') )
 content <- fetch_with_backoff(url = url2)
 parsed_page <- read_html(content)
 

 # Extract the GSE number using a regular expression

 long_string <- as.character(parsed_page)
 
 matches2 <- str_extract_all(long_string, 'type="String">GSE[0-9]+</item><item')
 
 
 matches2 <- unlist( gsub('type=\"String\">', '',matches2[[1]] ) )
 GSE_list[[i]]<- gsub('</item><item', '',matches2 )
 
 
 
 pause_time <- sample(1:2)[1]
 
 message(paste("Pausing for", pause_time,  "seconds...") )
 Sys.sleep(pause_time)
 message("Resuming execution...")
 
 
  }
  
  ###### remaining numbers ####### 
  
  if( length(matches)%%100 != 0){
  base_remain <- (number_of_id*100) 
  remain <- (base_remain +1 ):(base_remain + length(matches)%%100 )
  
  }else{ remain <- (number_of_id*100) }
  
  
  url3 <-  paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=', paste0(matches[remain], collapse = ',') )
  content <- fetch_with_backoff(url = url3)
  parsed_page <- read_html(content)
  
  
  # Extract the GSE number using a regular expression
  
  long_string <- as.character(parsed_page)
  
  matches4 <- str_extract_all(long_string, 'type="String">GSE[0-9]+</item><item')
  
  
  matches4 <- unlist( gsub('type=\"String\">', '',matches4[[1]] ) )
  GSE_list2<- gsub('</item><item', '',matches4 )
  
  return( unique(c(unlist(GSE_list), GSE_list2)) )
  
  
  ###### If the legnth of result is less than 100 
  }else{
    
    
    url2 <-  paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=', paste0(matches, collapse = ',') )
    content <- fetch_with_backoff(url = url2)
    parsed_page <- read_html(content)
    
    
    # Extract the GSE number using a regular expression
    
    long_string <- as.character(parsed_page)
    
    matches <- str_extract_all(long_string, 'type="String">GSE[0-9]+</item><item')
    
    
    matches <- unlist( gsub('type=\"String\">', '',matches[[1]] ) )
    GSE_list<- gsub('</item><item', '',matches )
    
    
    return(GSE_list)

  }
  
  

}

######################################### GEOscrap #############################

GEOscrap <- function (GSE = GSE, attr = attr ){
  
  user_agents <- c(
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.113 Safari/537.36",
    "Mozilla/5.0 (Windows NT 6.1; WOW64; Trident/7.0; AS; rv:11.0) like Gecko",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_6) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15", 
    'Mozilla/5.0 (iPhone; CPU iPhone OS 14_5 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.0 Mobile/15E148 Safari/604.1',
    'Mozilla/5.0 (Linux; Android 10; SM-G973F) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Mobile Safari/537.36', 
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.0.3 Safari/605.1.15',
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Safari/537.36', 
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:88.0) Gecko/20100101 Firefox/88.0',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.818.42 Safari/537.36 Edg/90.0.818.42',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0', 
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.85 Safari/537.36'
    
    
  )
  
  fetch_with_backoff <- function(url, max_attempts = 5, initial_wait = 1) {
    attempt <- 1
    success <- FALSE
    response <- NULL
    
    while (attempt <= max_attempts && !success) {
      user_agent <- sample(user_agents, 1)
      
      response <- GET(url, user_agent(user_agent))
      
      if (status_code(response) == 429) {
        wait_time <- initial_wait * (2 ^ (attempt - 1))  # Exponential backoff
        message(paste("Received 429, waiting", wait_time, "seconds before retrying..."))
        Sys.sleep(wait_time)
        attempt <- attempt + 1
      } else {
        success <- TRUE
      }
    }
    
    if (!success) {
      stop("Failed to fetch the URL after ", max_attempts, " attempts.")
    }
    
    return(content(response, as = "text"))
  }
  
  result_df <- as.data.frame( matrix ( rep('NOT found', 15*length(GSE)),ncol = 15, nrow = length(GSE)) )
  colnames(result_df) <- c('GSE_ID','title', 'species', 'experiment_type', 'Summary','overall_design', 'PMID', 
                           'country', 'Platform', 'Bioproject', 'SRA','extracted_number_of_sample','extracted_processdata_stat', 'extracted_rawdata_stat', 'chracteristics')
  for ( i in 1:length(GSE)){
    GSE_id <- GSE[i]
    url1 <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSE_id)
    
    content <- fetch_with_backoff(url = url1)
    parsed_page <- read_html(content)
    
    
    shared_path <-  parsed_page %>% html_node("body") 
    ### Extracted data 
    
    
    # GSE ID 
    
    result_df[i,which(colnames(result_df) == 'GSE_ID')] <- GSE_id  
    
    
    #title 
    title <-shared_path %>% html_node(xpath = "//td[text()='Title']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'title')] <- ifelse(length(title == 1),  title , 'NA or more than 1 obj')  
    
    
    #species 
    species <- shared_path %>% html_node(xpath = "//td[text()='Organism']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'species')] <- ifelse(length(species == 1),  species , 'NA or more than 1 obj')  
    
    
    #experiemtn_type 
    experiment_type <- shared_path %>% html_node(xpath = "//td[text()='Experiment type']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'experiment_type')] <-ifelse(length(experiment_type == 1),  experiment_type , 'NA or more than 1 obj')  
    
    #### Summary 
    Summary <-shared_path %>% html_node(xpath = "//td[text()='Summary']/following-sibling::td") %>% html_text(trim = TRUE)
    
    # loop  for summary
    summary_attr <- list()
    for (z in 1:length(attr) ){
      summary_attr[z] <- ifelse(length(grep(attr[z], Summary) ) == 1, attr[z], paste(attr[z], 'was not found') ) 
    }
    combined_string_summary <- paste(unlist(summary_attr), collapse = " AND")
    
    result_df[i,which(colnames(result_df) == 'Summary')] <- ifelse(length(Summary == 1), combined_string_summary  , 'NA or more than 1 obj')  
    
    ####overall design 
    overall_design <-shared_path %>% html_node(xpath = "//td[text()='Overall design']/following-sibling::td") %>% html_text(trim = TRUE)
    
    # loop  for design
    design_attr <- list()
    for (z in 1:length(attr) ){
      design_attr[z] <- ifelse(length(grep(attr[z], Summary) ) == 1, attr[z], paste(attr[z], 'was not found') ) 
    }
    combined_string_design <- paste(unlist(design_attr), collapse = " AND")
    
    result_df[i,which(colnames(result_df) == 'overall_design')] <- ifelse(length(Summary == 1), combined_string_design  , 'NA or more than 1 obj')  
    
    
    
    # PMID
    PMID <- shared_path %>% html_node(xpath = "//td[text()='Citation(s)']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'PMID')] <- ifelse(length(PMID == 1),  PMID , 'NA or more than 1 obj')  
    
    #Country
    country <- shared_path %>% html_node(xpath = "//td[text()='Country']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'country')] <-ifelse(length(country == 1),  country , 'NA or more than 1 obj')  
    
    # Platform
    Platform <- shared_path %>% html_node(xpath = "//td[contains(., 'GPL')]/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'Platform')] <- ifelse(length(Platform == 1),  Platform , 'NA or more than 1 obj')  
    
    #shared_path %>% html_node(xpath = "//td[contains(text(), 'Samples')]/following-sibling::td")  %>% html_nodes(xpath = "./following-sibling::td") %>% html_text(trim = TRUE)
    
    ## Bioproject 
    
    Bioproject <- shared_path %>% html_node(xpath = "//td[text()='BioProject']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'Bioproject')] <- ifelse(length(Bioproject == 1),  Bioproject , 'NA or more than 1 obj')  
    
    ## SRA
    SRA <- shared_path %>% html_node(xpath = "//td[text()='SRA']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'SRA')] <-  ifelse(length(SRA == 1),  SRA , 'NA or more than 1 obj')  
    
    
    ## Number of samples
    number_of_sample <- shared_path  %>% 
      html_nodes('tr')  %>% html_nodes('td')   %>% html_node(xpath = "//td[contains(text(), 'Samples')]") 
    extracted_number_of_sample <- unique(gsub(".*\\((\\d+)\\).*", "\\1", number_of_sample) )
    result_df[i,which(colnames(result_df) == 'extracted_number_of_sample')] <- ifelse(length(extracted_number_of_sample == 1),  extracted_number_of_sample , 'NA or more than 1 obj')  
    
    
    ##Processed data 
    processed_data_status <- shared_path %>%  html_nodes('tr')  %>% html_nodes('td') %>%  html_node(xpath = "//td[contains(text(), 'Processed data')]") 
    extracted_processdata_stat <- unique( gsub(".*>(.*)</td>", "\\1", processed_data_status) )
    result_df[i,which(colnames(result_df) == 'extracted_processdata_stat')] <- ifelse(length(extracted_processdata_stat == 1),  extracted_processdata_stat , 'NA or more than 1 obj')  
    
    ## Raw data
    raw_data_status <-   shared_path %>%  html_nodes('tr')  %>% html_nodes('td') %>%  html_node(xpath = "//td[contains(text(), 'Raw data')]") 
    extracted_rawdata_stat <- unique( gsub(".*>(.*)</td>", "\\1", raw_data_status) )
    result_df[i,which(colnames(result_df) == 'extracted_rawdata_stat')] <- ifelse(length(extracted_rawdata_stat == 1),  extracted_rawdata_stat , 'NA or more than 1 obj')  
    
    
    ################## Temp Pause ######################
    pause_time <- sample(1:2)[1]
    
    message(paste("Pausing for", pause_time,  "seconds...") )
    Sys.sleep(pause_time)
    message("Resuming execution...")
    
    ##################  Extracting the sample MetaData
    
    sample_name <- parsed_page %>% html_node("body")  %>% html_nodes('td') %>% html_nodes('a') %>% html_text() 
    
    first_sample_name <- sample_name[grep('GSM.*.', sample_name)][1]
    
    
    url_sample <- 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
    url_sample1 <- paste0(url_sample,first_sample_name ) 
    
    
    content_sample <- fetch_with_backoff(url_sample1)
    parsed_page_sample <- read_html(content_sample)
    
    
    shared_path_sample  <-  parsed_page_sample %>% html_node("body") 
    chracteristics <- shared_path_sample  %>% html_node(xpath = "//td[text()='Characteristics']/following-sibling::td") %>% html_text(trim = TRUE)
    result_df[i,which(colnames(result_df) == 'chracteristics')] <- ifelse(length(chracteristics == 1),  chracteristics , 'NA or more than 1 obj') 
    
    ################## Temp Pause ######################
    pause_time <- sample(1:2)[2]
    
    message(paste("Pausing for", pause_time,  "seconds...") )
    Sys.sleep(pause_time)
    message("Resuming execution...")
    message(paste(length(GSE) - i, 'remaining... '))
    
    
    
  }
  
  return(result_df)
  
}







