

#' Loads a metaphlan table
#'
#' Loads a metaphlan table filtered to a certain taxonomic level
#'
#' @param taxonomic_level A taxonomic level, which loads the corresponding metaphlan table
#'
#' @export
#'
#' @example
#' metaphlan <- load_metaphlan("species")
load_metaphlan <- function(taxonomic_level = "genus") {
    metaphlan <- readr::read_delim(paste0("../data/metaphlan3/", "MetaPhlAn3_Subsample_5000000_", taxonomic_level, ".csv"), delim = "\t", show_col_types = FALSE)
    return(metaphlan)
}



#' Metaphlan taxonomy filter
#'
#' @param metaphlan 
#' @param taxonomy 
#' @param keep_unclassified 
#'
#' @return
#' @export
#'
#' @examples
filter_taxonomy <- function(metaphlan, taxonomy = "species", keep_unclassified = TRUE, remove_tax_names = TRUE){
    tax <- tibble(
        abbr_taxonomy = c("k", "p", "c", "o", "f", "g", "s", "t"), 
        full_taxonomy = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
    )
    
    tax_lvl <- filter(tax,full_taxonomy == taxonomy)$abbr_taxonomy
    next_tax_lvl <- filter(tax, abbr_taxonomy == tax$abbr_taxonomy[grep(tax_lvl, tax$abbr_taxonomy) + 1])$abbr_taxonomy
    
    
    if(keep_unclassified){
        df <- metaphlan %>% 
            filter(
                clade_name == "UNCLASSIFIED" | str_detect(clade_name, paste0(tax_lvl,"_{2}")), 
                str_detect(clade_name, paste0(next_tax_lvl,"_{2}")) == FALSE
            )
    } else {
        df <- metaphlan %>% 
            filter(
                str_detect(clade_name, paste0(tax_lvl,"_{2}")), 
                str_detect(clade_name, paste0(next_tax_lvl,"_{2}")) == FALSE
            )
    }
    
    if(taxonomy == "strain"){
        df <- metaphlan %>% 
            filter(
                clade_name == "UNCLASSIFIED" | str_detect(clade_name, paste0(tax_lvl,"_{2}"))
            )
    }
    
    if(remove_tax_names){
        df <- df %>% 
            mutate(
                clade_name = str_remove(clade_name, paste0(".*",tax_lvl,"__"))
            )
    }
    return(df)
    
    
}


#' Rename columns in metaphlan
#'
#' @param metaphlan 
#'
#' @return
#' @export
#'
#' @examples
remove_NonHuman_from_colnames <- function(metaphlan){
    colnames_old <- colnames(metaphlan)
    colnames_new <- str_remove(colnames_old, "\\_Non.*")
    
    
    return(
        metaphlan %>% 
            rename_at(colnames_old, ~ colnames_new)
    )
}



#' Calculate the species richness
#'
#' Calculate the species richness defined as number of species present in a sample. Filter_species will first filter all species with less than x relative abundance in at least one sample
#'
#' @param metaphlan A metaphlan dataframe
#' @param filter_species A cutoff value in relative abundance.
#'
#' @export
#'
#' @example
#' species_richness <- calulateSpeciesRichness(metaphlan, filter_species = 0.1)
calculateSpeciesRichness <- function(metaphlan, filter_species = 0) {
    species_count_for_each_sample <- metaphlan %>%
        filter(clade_name != "UNCLASSIFIED") %>% 
        #select(-NCBI_tax_id) %>%
        pivot_longer(-clade_name, names_to = "sample_barcode") %>%
        group_by(sample_barcode) %>%
        summarise(richness = sum(value > filter_species))

    return(species_count_for_each_sample)
}

#' transpose metaphlan table
#'
#' transpose metaphlan table
#'
#' @param metaphlan A metaphlan dataframe
#'
#' @export
#'
#' @example
#' t_metaphlan <- transposeMetaphlan(metaphlan)
transposeMetaphlan <- function(metaphlan) {
    t_metaphlan <- metaphlan %>%
        #select(-NCBI_tax_id) %>%
        pivot_longer(-clade_name, names_to = "sample_barcode") %>%
        pivot_wider(names_from = clade_name, values_from = value) %>%
        arrange(sample_barcode) %>%
        column_to_rownames(var = "sample_barcode")


    return(t_metaphlan)
}
