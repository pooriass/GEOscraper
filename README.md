########### GEOscraper #############
Searching Gene Expression Omnibus is one of the main challenges at the beginning of a bioinformatic project. Since I couldn't find any specific packages or functions, or any feasible facility to download and fetch the GEO metadata of different datasets, I have wrapped up some self-made functions to download and stash the metadata of GEO datasets.

Steps to accomplish your search:

1. Use the FindGSE function to find the GSE ids of the datasets based on your search term.
2. Use the GEOscrap function to get a data frame containing metadata of different datasets.
Notes:

1. Be careful about setting the arguments in terms of any possible typos.
2. Your search term should be in a vector, such as: c('(colorectal cancer)', 'AND', '(scRNAseq)'); pay attention to the parentheses.
3. There are multiple pauses throughout the functions in order to prevent too many requests (HTTP error 429) but NCBI has declared that without an API, you can send 3 requests per second. However, I preferred to seal the deal way more conservative!
4. This function has not been rigorously tested, just two or three times to get a ballpark. In terms of any further problems, bugs, etc., I would be happy if you let me know.

Any suggestions or wrapping it up into a package would make me deliriously happy.
Contact me: Salehipooria54@gmail.com

Thanks to my sweetheart, ChatGPT...
