taxname.correct <- function(name_of_tax){
  library(dplyr)
  ###d__Bacteria
  ###phylum
  x <- name_of_tax %>% 
    sub("d__Bacteria\\.__\\.__\\.__\\.__\\.__\\.__$","d__Bacteria\\.p__\\.c__\\.o__\\.f__\\.g__\\.s__", .) %>%
    sub("d__Bacteria\\.__\\.__\\.__\\.__\\.__$","d__Bacteria\\.p__\\.c__\\.o__\\.f__\\.g__", .) %>%
    sub("d__Bacteria\\.__\\.__\\.__\\.__$","d__Bacteria\\.p__\\.c__\\.o__\\.f__", .) %>%
    sub("d__Bacteria\\.__\\.__\\.__$","d__Bacteria\\.p__\\.c__\\.o__", .) %>%
    sub("d__Bacteria\\.__\\.__$","d__Bacteria\\.p__\\.c__", .) %>%
    sub("d__Bacteria\\.__$","d__Bacteria\\.p__", .)
  ###class
  x <- x %>% 
    sub("d__Bacteria\\.p__(\\w*)\\.__\\.__\\.__\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\.o__\\.f__\\.g__\\.s__",.)%>% 
    sub("d__Bacteria\\.p__(\\w*)\\.__\\.__\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\.o__\\.f__\\.g__",.) %>%
    sub("d__Bacteria\\.p__(\\w*)\\.__\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\.o__\\.f__",.)%>%
    sub("d__Bacteria\\.p__(\\w*)\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\.o__",.)%>%
    sub("d__Bacteria\\.p__(\\w*)\\.__$","d__Bacteria.p__\\1\\.c__",.)
  ###order
  x <- x %>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.__\\.__\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\.f__\\.g__\\.s__",.)%>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.__\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\.f__\\.g__",.) %>%
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\.f__",.) %>%
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__",.) 
  ###family
  x <- x %>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.__\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\3\\.f__\\.g__\\.s__",.)%>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\3\\.f__\\.g__",.)%>%
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\3\\.f__",.)
  ###genus
  x <- x %>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.f__(\\w*)\\.__\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\3\\.f__\\4\\.g__\\.s__",.)%>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.f__(\\w*)\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\3\\.f__\\4\\.g__",.)
  ###species
  x <- x %>% 
    sub("d__Bacteria\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.f__(\\w*)\\.g__(\\w*)\\.__$","d__Bacteria.p__\\1\\.c__\\2\\.o__\\3\\.f__\\4\\.g__\\5\\.s__",.)
  #######d__Archaea
  ###phylum
  x <- x %>% 
    sub("d__Archaea\\.__\\.__\\.__\\.__\\.__\\.__$","d__Archaea\\.p__\\.c__\\.o__\\.f__\\.g__\\.s__", .) %>%
    sub("d__Archaea\\.__\\.__\\.__\\.__\\.__$","d__Archaea\\.p__\\.c__\\.o__\\.f__\\.g__", .) %>%
    sub("d__Archaea\\.__\\.__\\.__\\.__$","d__Archaea\\.p__\\.c__\\.o__\\.f__", .) %>%
    sub("d__Archaea\\.__\\.__\\.__$","d__Archaea\\.p__\\.c__\\.o__", .) %>%
    sub("d__Archaea\\.__\\.__$","d__Archaea\\.p__\\.c__", .) %>%
    sub("d__Archaea\\.__$","d__Archaea\\.p__", .)
  ###class
  x <- x %>% 
    sub("d__Archaea\\.p__(\\w*)\\.__\\.__\\.__\\.__\\.__$","d__Archaea.p__\\1\\.c__\\.o__\\.f__\\.g__\\.s__",.)%>% 
    sub("d__Archaea\\.p__(\\w*)\\.__\\.__\\.__\\.__$","d__Archaea.p__\\1\\.c__\\.o__\\.f__\\.g__",.) %>%
    sub("d__Archaea\\.p__(\\w*)\\.__\\.__\\.__$","d__Archaea.p__\\1\\.c__\\.o__\\.f__",.)%>%
    sub("d__Archaea\\.p__(\\w*)\\.__\\.__$","d__Archaea.p__\\1\\.c__\\.o__",.)%>%
    sub("d__Archaea\\.p__(\\w*)\\.__$","d__Archaea.p__\\1\\.c__",.)
  ###order
  x <- x %>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.__\\.__\\.__\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\.f__\\.g__\\.s__",.)%>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.__\\.__\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\.f__\\.g__",.) %>%
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.__\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\.f__",.) %>%
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__",.) 
  ###family
  x <- x %>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.__\\.__\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\3\\.f__\\.g__\\.s__",.)%>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.__\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\3\\.f__\\.g__",.)%>%
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\3\\.f__",.)
  ###genus
  x <- x %>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.f__(\\w*)\\.__\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\3\\.f__\\4\\.g__\\.s__",.)%>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.f__(\\w*)\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\3\\.f__\\4\\.g__",.)
  ###species
  x <- x %>% 
    sub("d__Archaea\\.p__(\\w*)\\.c__(\\w*)\\.o__(\\w*)\\.f__(\\w*)\\.g__(\\w*)\\.__$","d__Archaea.p__\\1\\.c__\\2\\.o__\\3\\.f__\\4\\.g__\\5\\.s__",.)
  
  return(x)
}


extract.genus <- function(name_of_tax){
  name_of_tax.matrix <- as.data.frame(name_of_tax)
  name_of_tax <- taxname.correct(name_of_tax)
  genus.list <- strsplit(name_of_tax,"g__")
  genus.list.length <- lengths(genus.list)
  name_of_tax.matrix$genus <- "un"
  name_of_tax.matrix$name_of_tax <- as.character(name_of_tax.matrix$name_of_tax )
  for(i in 1:length(name_of_tax)){
   if(genus.list.length[i] == 2 & !(grepl("uncultured",name_of_tax[i])) & !grepl("Incertae_Sedis",name_of_tax[i])){
     tax.name <-  unlist(lapply(genus.list[i], '[[', 2))
     tax.name <- tax.name %>% gsub("^\\.","",.) %>% gsub("\\.$","",.) 
     name_of_tax.matrix[i,"genus"] <-  tax.name
   }
    else if(grepl("uncultured",name_of_tax[i])){
      split.tax.list <- unlist(strsplit(name_of_tax.matrix$name_of_tax[i],"__"))
      tax.name <- split.tax.list[6]
      tax.name <- tax.name %>% gsub("\\.g" , "",.) %>% gsub("\\.f" , "",.) %>% gsub("\\.o" , "",.) %>% gsub("\\.c" , "",.) %>% gsub("\\.p" , "",.) %>% gsub("\\.$","",.) %>% gsub("^\\.","",.) 
      tax.name <- paste("Uc_g",tax.name,sep = "_")
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
    else if(grepl("Incertae_Sedis",name_of_tax[i])){
      split.tax.list <- unlist(strsplit(name_of_tax.matrix$name_of_tax[i],"__"))
      tax.name <- split.tax.list[6]
      tax.name <- tax.name %>% gsub("\\.g" , "",.) %>% gsub("\\.f" , "",.) %>% gsub("\\.o" , "",.) %>% gsub("\\.c" , "",.) %>% gsub("\\.p" , "",.) %>% gsub("\\.$","",.) %>% gsub("^\\.","",.) 
      tax.name <- paste(tax.name,"Incertae_Sedis",sep = "_")
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
    else{
      split.tax.list <- unlist(strsplit(name_of_tax.matrix$name_of_tax[i],"__"))
      id <- max(which(nchar(split.tax.list)>2))
      tax.name <- split.tax.list[id]
      tax.name <- tax.name %>% gsub("\\.g" , "",.) %>% gsub("\\.f" , "",.) %>% gsub("\\.o" , "",.) %>% gsub("\\.c" , "",.) %>% gsub("\\.p" , "",.) %>% gsub("\\.$","",.) %>% gsub("^\\.","",.) 
      tax.name <- paste("Un_g",tax.name,sep = "_")
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
      
  }
  return(name_of_tax.matrix)
}


extract.genus.phlym <- function(name_of_tax){
  name_of_tax.matrix <- as.data.frame(name_of_tax)
  name_of_tax <- taxname.correct(name_of_tax)
  
  phlym.list <- strsplit(name_of_tax,".c__")
  name_of_tax.matrix$phlym.name <-  unlist(lapply(phlym.list, '[[', 1)) %>% gsub("d__Bacteria.p__","",.) 
  
  genus.list <- strsplit(name_of_tax,"g__")
  genus.list.length <- lengths(genus.list)
  name_of_tax.matrix$genus <- "un"
  name_of_tax.matrix$name_of_tax <- as.character(name_of_tax.matrix$name_of_tax )
  for(i in 1:length(name_of_tax)){
    if(genus.list.length[i] == 2 & !(grepl("uncultured",name_of_tax[i])) & !grepl("Incertae_Sedis",name_of_tax[i])){
      tax.name <-  unlist(lapply(genus.list[i], '[[', 2))
      tax.name <- tax.name %>% gsub("^\\.","",.) %>% gsub("\\.$","",.) 
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
    else if(grepl("uncultured",name_of_tax[i])){
      split.tax.list <- unlist(strsplit(name_of_tax.matrix$name_of_tax[i],"__"))
      tax.name <- split.tax.list[6]
      tax.name <- tax.name %>% gsub("\\.g" , "",.) %>% gsub("\\.f" , "",.) %>% gsub("\\.o" , "",.) %>% gsub("\\.c" , "",.) %>% gsub("\\.p" , "",.) %>% gsub("\\.$","",.) %>% gsub("^\\.","",.) 
      tax.name <- paste("Uc_g",tax.name,sep = "_")
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
    else if(grepl("Incertae_Sedis",name_of_tax[i])){
      split.tax.list <- unlist(strsplit(name_of_tax.matrix$name_of_tax[i],"__"))
      tax.name <- split.tax.list[6]
      tax.name <- tax.name %>% gsub("\\.g" , "",.) %>% gsub("\\.f" , "",.) %>% gsub("\\.o" , "",.) %>% gsub("\\.c" , "",.) %>% gsub("\\.p" , "",.) %>% gsub("\\.$","",.) %>% gsub("^\\.","",.) 
      tax.name <- paste(tax.name,"Incertae_Sedis",sep = "_")
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
    else{
      split.tax.list <- unlist(strsplit(name_of_tax.matrix$name_of_tax[i],"__"))
      id <- max(which(nchar(split.tax.list)>2))
      tax.name <- split.tax.list[id]
      tax.name <- tax.name %>% gsub("\\.g" , "",.) %>% gsub("\\.f" , "",.) %>% gsub("\\.o" , "",.) %>% gsub("\\.c" , "",.) %>% gsub("\\.p" , "",.) %>% gsub("\\.$","",.) %>% gsub("^\\.","",.) 
      tax.name <- paste("Un_g",tax.name,sep = "_")
      name_of_tax.matrix[i,"genus"] <-  tax.name
    }
    
  }
  return(name_of_tax.matrix)
}
