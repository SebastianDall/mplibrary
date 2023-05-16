library(tidyverse)

################################

#' Load metadata
#'
#' Loads metadata
#'
#'
#'
#' @export
#'
#' @example metadata <- load_metadata()
#'
load_metadata <- function() {
    return(
        read_delim("../data/metadata/metadata.csv", delim = ";", show_col_types = FALSE, col_types = list(batch_3 = col_double())) %>% 
            filter(!is.na(sample_barcode)) %>% 
            distinct(sample_barcode, .keep_all = TRUE)
    )
}

################################
#' Isolates metadata for project to be analyzed
#'
#' Functionen for isolating relevant metadata. Three projects can be used: MP or NP
#'
#' @param metadata A dataframe
#' @param project_filter A project
#'
#' @export
#'
#' @example metadata_MP <- isolateProjectMetadata(metadata, "MP")
#'
isolateProjectMetadata <- function(metadata, project_filter = "MP") {
    stopifnot(project_filter %in% c("MP", "NP"))
    patientMetadataWithDonorbatches <- isolateDonorBatchesUsed(metadata, project_filter)

    DonorAndPatientMetadata <- isolateDonorAndPatientMetadata(metadata, patientMetadataWithDonorbatches, project_filter) %>%
        createXaxis()

    return(DonorAndPatientMetadata)
}

isolateDonorBatchesUsed <- function(metadata, project_filter) {
    patientMetadataWithDonorbatches <- metadata %>%
        pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>%
        filter(!is.na(batch_number), project == project_filter) %>%
        mutate(donor_batch = paste0(donor, "_", batch_number))

    return(patientMetadataWithDonorbatches)
}

isolateDonorAndPatientMetadata <- function(metadata, metadata_with_donor_batches_used, project_filter) {
    DonorAndPatientMetadata <- metadata %>%
        mutate(donor_batch = paste0(id, "_", fecal_donation_number)) %>%
        filter(project == project_filter | donor_batch %in% metadata_with_donor_batches_used$donor_batch) %>%
        mutate(group = if_else(project == "donor_batch", "FMT", group))

    return(DonorAndPatientMetadata)
}

createXaxis <- function(df) {
    dfWithXaxis <- df %>%
        mutate(x_axis = if_else(stage == "inclusion", "Pre",
            if_else(stage == "followup_30d", "Post", stage)
        )) %>%
        mutate(x_axis = if_else(project == "donor_batch", "Donor", x_axis))
}
