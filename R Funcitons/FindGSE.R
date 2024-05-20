


###### FindGSE

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
  pause_time <- sample(3:5)[1]

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
  pause_time <- sample(3:5)[1]

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
  pause_time <- sample(3:5)[1]

  message(paste("Pausing for", pause_time,  "seconds...") )
  Sys.sleep(pause_time)
  message("Resuming execution...")


  ###### IDs to GSE

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




