#### Example 

term <- c('(colorectal cancer)', 'AND', '(scRNAseq)')



ids <- FindGSE(Search_Term = term, Organism = 'Homo sapiens', Study_type = 'Expression profiling by high throughput sequencing', max_result = 1000)

results <- GEOscrap(GSE = ids, attr = 'drug resistance')

