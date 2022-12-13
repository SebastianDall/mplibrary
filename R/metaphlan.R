

#' Loads a metaphlan table
#' #'
#' Loads a metaphlan table filtered to a certain taxonomic level
#'
#' @param taxonomic_level A taxonomic level, which loads the corresponding metaphlan table
#'
#' @export
#'
#' @example
#' metaphlan <- load_metaphlan("species")
load_metaphlan <- function(taxonomic_level = "genus") {
    metaphlan <- read_delim(paste0("../data/metaphlan3/", "MetaPhlAn3_Subsample_5000000_", taxonomic_level, ".csv"), delim = "\t")
    return(metaphlan)
}



#' Calculate the species richness
#' #'
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
        select(-NCBI_tax_id) %>%
        pivot_longer(-clade_name, names_to = "LibID") %>%
        group_by(LibID) %>%
        summarise(richness = sum(value > filter_species))

    return(species_count_for_each_sample)
}

#' transpose metaphlan table
#' #'
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
        select(-NCBI_tax_id) %>%
        pivot_longer(-clade_name, names_to = "LibID") %>%
        pivot_wider(names_from = clade_name, values_from = value) %>%
        arrange(LibID) %>%
        column_to_rownames(var = "LibID")


    return(t_metaphlan)
}
