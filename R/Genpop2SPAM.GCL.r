#' Converts Genepop files to SPAM baseline and mixture files - IN DEVELOPMENT
#'
#' Genepop files contain the information to write the baseline ".bse" and mixture ".mix" files for SPAM. Use the following considerations when using this function. 
#'
#' \itemize{
#'   \item A. This function is set up to work only with diploid alleles for SNPs and microsatellites or mtDNA haplotypes with nine or less combined SNPs.  
#'   \item B. All the files for the baseline and mixtures should be read in at the same time to create SPAM files.  This is because for each locus it will need to know how many different alleles are present across all files to properly id each one into the correct position. 
#'   \item C. Each population in the Genepop file should be separated by POP, Pop, or pop.  Any other population separator will not be recognized.  
#'   \item D. The locus names listed at the top of each Genepop file should be in the same order and have the same number of loci. 
#'    ...
#' }
#' @param baselineFile The fullpath to the baseline Genepop file for the baseline in SPAM (.mix). 
#' @param mixtureFiles A vector with the fullpath of the Genepop files to be created into seperate SPAM mixture (.mix) files.  
#' @param baselineOutPath The full path to write the baseline file to
#' @param mixtureOutPath The full path to write all of the mixture files to. 
#' @keywords  ADFG
#' @export  
Genepop2SPAM.GCL<- function(baselineFile, mixtureFiles, 
                            baselineOutPath, mixtureOutPath) {
    #load a great text manipulation package
    require(stringr)
    
    if(missing(mixtureFiles)) {mixtureFiles <- NULL}
    if(missing(baselineFile)) {baselineFile <- NULL}
    if(missing(baselineOutPath)) {baselineOutPath <- baselineFile}
    if(missing(mixtureOutPath)) {mixtureOutPath <- mixtureFiles}
    
    ###GET THE TOTAL ALLELE COUNTS WITHIN ALL FILES###
    all_files <- c(baselineFile, mixtureFiles)
    
    all_data <- c()
    file_length <- c()
    
    for(i in seq(all_files)) {
        file_data <- readLines(all_files[i])
        file_length <- c(file_length, length(file_data))
        all_data <- c(all_data, file_data)
    }
    file_start <- c(1, cumsum(file_length))[-(1 + length(file_length))]
    #find lines that seperate SILLYs in genepop file
    all_pop_lines <- which(all_data == "Pop")
    all_pop_lines <- c(all_pop_lines, which(all_data == "POP"))
    all_pop_lines <- c(all_pop_lines, which(all_data == "pop"))
    #get the assays listed in the top of the file
    all_assays <- all_data[2:(all_pop_lines[1] - 1)]
    #
    pop_start_lines <- sapply(file_start, 
                              function(x) { 
                                  all_pop_lines[which(all_pop_lines == all_pop_lines[x - all_pop_lines < 0][1])]
                              }
    )
    #lines_to_remove is all the lines that do not contain sample data
    lines_to_remove <- c()
    
    for(v in seq(pop_start_lines)) {
        lines_to_remove <- c(lines_to_remove ,file_start[v]:pop_start_lines[v])
    }
    
    #get just the data
    all_sample_data_all <- all_data[-c(lines_to_remove, all_pop_lines)]
    #write over the data to remove the sample ids
    all_sample_data <-  t(apply(cbind(all_sample_data_all),
                                1,
                                function(str){
                                    vec=strsplit(str, split=c(" ", "\t", ","))[[1]]
                                    vec[!vec%in%c(" ","\t",",","")]
                                }
    )
    )
    #get the assay factor levels (number of different alleles in the assay)
    allelevec <- apply(all_sample_data[,-1], 
                       2,
                       function(dat) {
                           
                           
                           #dat <- all_sample_data[,2]
                           if(nchar(dat[1]) == 3 | nchar(dat[1]) == 2) {
                               dat[dat == "00" | dat == "000"] <- NA   
                               length(levels(factor(dat))) 
                           }
                           else {
                               n <- nchar(dat[1])
                               first_allele  <- substr(dat, 2, n / 2)
                               first_allele[first_allele == "00" | first_allele == "000"] <- NA     
                               second_allele <- substr(dat, (n / 2) + 2, n)
                               second_allele[second_allele == "00" | second_allele == "000"] <- NA 
                               length(levels(factor(c(first_allele, second_allele))))
                           }
                       }
    )
    
    
    alleles <- list()
    for(i in seq(dim(all_sample_data[,-1])[2])) {
        
        if(nchar(all_sample_data[,i + 1]) == 3 | nchar(all_sample_data[,i + 1]) == 2) {
            alleles_per_locus <- gsub("(^|[^0-9])0+", "\\1", all_sample_data[,15], perl = TRUE)
            alleles_per_locus[alleles_per_locus == ""] <- NA
            alleles <- c(alleles, list(as.numeric(names(table(alleles_per_locus)))))
            
        }
        else {
            n <- nchar(all_sample_data[,i + 1][1])
            first_allele  <- substr(all_sample_data[,i + 1], 2, n / 2)
            first_allele[first_allele == "00" | first_allele == "000"] <- NA     
            second_allele <- substr(all_sample_data[,i + 1], (n / 2) + 2, n)
            second_allele[second_allele == "00" | second_allele == "000"] <- NA 
            alleles <- c(alleles, list(as.numeric(names(table(as.numeric(gsub("(^|[^0-9])0+", 
                                                                              "\\1", 
                                                                              c(first_allele, second_allele), perl = TRUE)))))))
        }
    }    
    allele_levels <- lapply(alleles, function(x) { levels(factor(x)) } )
    
    ### END GET THE TOTAL ALLELE COUNTS WITHIN ALL FILES (allelevec)###
    
    ####CREATE BASELINE####    
    #read the whole file in as a vector of lines
    whole_file <- readLines(baselineFile)
    #find lines that seperate SILLYs in genepop file
    pop_lines <- which(whole_file == "Pop")
    #get the assays listed in the top of the file
    assays   <- whole_file[2:(pop_lines[1] - 1)]
    #get just the data
    sample_data_all <- whole_file[-c(1:pop_lines[1], pop_lines)]
    #write over the data to remove the sample ids
    sample_data <- t(apply(cbind(sample_data_all),
                           1,
                           function(str){
                               vec=strsplit(str, split=c(" ", "\t", ","))[[1]]
                               vec[!vec%in%c(" ","\t",",","")]
                           }
    )
    )
    #get just the sillys from the sample names (sillys are the 
    #categorical code used at the AKFG Gene Conservation Lab) 
    sample_SILLYs <- unlist(lapply(sample_data[,1], 
                                   function(x) { strsplit(x, "_")[[1]][1] } 
    )
    )
    SILLYs <- levels(factor(sample_SILLYs))
    #pre-write the file to a list
    baseline <- list()
    #Iterate over each silly, find the frequency of each allele and 
    #put the right number alleles for the assay, ... and 
    #get the whole file in a list 
    for(i in seq(SILLYs)) { 
        #subset the samples by SILLY 
        samples <- sample_SILLYs == SILLYs[i]
        #remove sample names
        file_mat <-  sample_data[samples, -1]
        #
        allele_list <- apply(file_mat, 2, function(dat) {
            #dat <- file_mat[,1]
            if(nchar(dat[1]) == 3 | nchar(dat[1]) == 2) {
                as.numeric(gsub("(^|[^0-9])0+", "\\1", dat, perl = TRUE))
            }
            else {
                #take the first allele
                first_allele  <- substr(dat, 2, nchar(dat) / 2)
                #remove the leading zero if there is one
                first_allele[substring(first_allele, 1, 1) == "0"] <- substring(first_allele[substring(first_allele, 1, 1) == "0"], 2)
                #We don't want to count the zeros. 
                first_allele[first_allele == ""] <- NA
                
                #repeat for the second allele
                second_allele <- substr(dat, nchar(dat[1]) / 2 + 2, nchar(dat[1]))
                second_allele[substring(second_allele, 1, 1) == "0"] <- substring(second_allele[substring(second_allele, 1, 1) == "0"], 2)
                second_allele[second_allele == ""] <- NA
                c(first_allele, second_allele)
            }
        })
        #start out the vector of allele counts with the header line needed for each SILLY of the baseline file
        silly_chunk <-  paste("#  ", which(SILLYs == SILLYs[i]), 
                              "    ", 
                              SILLYs[i], sep = "")
        
        #sum up the different alleles for each SILLY and put them in the right position of possible alleles. 
        if(!is.list(allele_list)) { allele_list <- apply(allele_list, 2, as.list)}
        
        for(j in seq(assays)) {
            freqPop <- unname(table(factor(as.numeric(as.numeric(allele_list[[j]])), levels = alleles[[j]])))
            the_line <- paste(assays[j], paste(freqPop, collapse = " "))
            silly_chunk <- c(silly_chunk, the_line)
        }
        #write the chunk to a list element
        baseline[[i]] <- silly_chunk
    }
    
    #write the file. 
    fileConn <- file(paste(gsub(".gen",
                                "",
                                baselineOutPath),
                           sep = ""))
    writeLines(unlist(baseline), fileConn)
    close(fileConn)
    ####END OF CREATE BASELINE####
    
    ####CREATE MIXTURE####
    for(z in seq(mixtureFiles)) {
        #read the whole file in as a vector of lines
        whole_file <- readLines(mixtureFiles[z])
        #find lines that seperate SILLYs in genepop file
        pop_lines <- which(whole_file == "Pop")
        #get the assays listed in the top of the file
        assays   <- whole_file[2:(pop_lines[1] - 1)]
        #get just the data
        sample_data_all <- whole_file[-c(1:pop_lines[1], pop_lines)]
        #write over the data to remove the sample ids
        sample_data <- t(apply(cbind(sample_data_all),
                               1,
                               function(str){
                                   vec=strsplit(str, split=c(" ", "\t", ","))[[1]]
                                   vec[!vec%in%c(" ","\t",",","")]
                               }
        )
        )
        counts <- sample_data[,-1]
        alleles <- apply(counts,2, function(dat) {
            if(nchar(dat[1]) == 3 | nchar(dat[1]) == 2) {
                as.numeric(gsub("(^|[^0-9])0+", "\\1", dat, perl = TRUE))
            }
            else {
                #take the first allele
                first_allele  <- substr(dat, 2, nchar(dat) / 2)
                #remove the leading zero if there is one
                first_allele[substring(first_allele, 1, 1) == "0"] <- substring(first_allele[substring(first_allele, 1, 1) == "0"], 2)
                first_allele[first_allele == "0"] <- NA
                #We don't want to count the zeros. 
                
                #repeat for the second allele
                second_allele <- substr(dat, nchar(dat[1]) / 2 + 2, nchar(dat[1]))
                second_allele[substring(second_allele, 1, 1) == "0"] <- substring(second_allele[substring(second_allele, 1, 1) == "0"], 2)
                second_allele[second_allele == "0"] <- NA
                c(first_allele, second_allele)
            }
        })
        
        alleles[alleles == ""] <- NA
        
        
        if(!is.list(alleles)) { alleles <- apply(alleles, 2, as.list)}
        
        
        
        n <- dim(counts)[1]
        sample_alleles <- list()
        for(i in seq(dim(counts)[1])) {
            sample_alleles <- c(sample_alleles, 
                                list(lapply(alleles, 
                                            function(x) { 
                                                as.numeric(c(x[i],
                                                             x[i + n/2])
                                                )
                                            }
                                )
                                ))
        }
        
        mix_file <- c()
        for(samples in seq(length(sample_alleles))) {
            sample_line <- c()
            for(assays in seq(allelevec)) {
                sample_line <-  c(sample_line, paste(table(factor(sample_alleles[[samples]][[assays]], 
                                                                  levels = as.numeric(allele_levels[[assays]]))), 
                                                     collapse = ""))
            }
            mix_file <- c(mix_file, paste(sample_line, collapse = " ") )
        }
        mix_file <- paste(mix_file, seq(mix_file))
        
        
        #write the mixture file 
        fileConn <- mixtureOutPath[z]
        writeLines(mix_file, fileConn)
        close(fileConn)
        
        ####END OF CREATE MIXTURE####
    }
}