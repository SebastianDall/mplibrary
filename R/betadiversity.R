
beta_diversity_v1 <- function(metadata, t_taxdata, transform = F, method = "bray", abundance_threshold = 10^-34, stages = c("Pre", "Post"), seed = 1, verbose = F) {
    stopifnot("method must be either sorensen or bray (default: bray)" = method %in% c("bray", "sorensen"))

    transform <- ifelse(method == "sorensen", F, transform)

    set.seed(seed)

    t_metaphlan_tib <- t_taxdata %>%
        rownames_to_column(var = "LibID")

    metadata_patients <- metadata %>%
        filter(x_axis != "Donor") %>%
        arrange(id)

    metadata_donor <- metadata %>%
        filter(x_axis == "Donor")


    t_metaphlan_donors <- t_metaphlan_tib %>%
        left_join(select(metadata_donor, LibID, id), by = "LibID") %>%
        relocate(id) %>%
        select(-LibID) %>%
        group_by(id) %>%
        summarise(across(where(~ is.numeric(.)), ~ mean(.)))





    beta_diversity <- tibble()

    for (pt_ID in unique(metadata_patients$id)) {
        patient <- metadata_patients %>%
            filter(id == pt_ID)


        stage_filter <- all(stages == c("Pre", "Post"))

        if (stage_filter) {
            hasPrePostTaxData <- patient %>%
                filter(x_axis %in% c("Pre", "Post")) %>%
                filter(LibID %in% t_metaphlan_tib$LibID)


            if (length(hasPrePostTaxData$id) != 2) {
                if (verbose) {
                    print(paste0(pt_ID, " did not complete treatment or sample was lost"))
                }
                next
            }
        }


        patient_stages <- patient %>%
            filter(x_axis %in% stages)

        isPatientFMT <- patient_stages$group[1] == "FMT"

        if (isPatientFMT) {
            patient_donors <- patient %>%
                pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>%
                filter(!is.na(batch_number)) %>%
                mutate(donor_batch = paste0(donor, "_", batch_number))

            donors <- metadata_donor %>%
                filter(donor_batch %in% patient_donors$donor_batch)

            excluded_donors <- metadata_donor %>%
                filter(!id %in% donors$id) %>%
                filter(LibID %in% t_metaphlan_tib$LibID) %>%
                mutate(id_number = parse_number(id))
        } else {
            excluded_donors <- metadata_donor %>%
                filter(LibID %in% t_metaphlan_tib$LibID) %>%
                mutate(id_number = parse_number(id))
        }




        metaphlan_patient <- t_metaphlan_tib %>%
            filter(LibID %in% patient_stages$LibID) %>%
            column_to_rownames(var = "LibID")

        nstages <- length(rownames(metaphlan_patient))

        set.seed(1)

        isPatientMP <- patient$project[1] == "MP"
        if (isPatientMP) {
            random_donors <- excluded_donors %>%
                # filter(id_number %in% sample(unique(excluded_donors$id_number), 4)) %>%
                arrange(id) %>%
                group_by(id_number) %>%
                mutate(n_samples = row_number()) %>%
                filter(n_samples == sample(1:n(), 1)) %>%
                ungroup()


            metaphlan_random_donor <- t_metaphlan_tib %>%
                filter(LibID %in% random_donors$LibID) %>%
                column_to_rownames("LibID")


            metaphlan_random_donor_patient <- bind_rows(metaphlan_patient, metaphlan_random_donor) %>%
                abundance_filter(abundance_threshold = abundance_threshold) %>%
                hellinger_transform(transform = transform) %>%
                as.matrix()


            beta_dist_random <- beta_method(metaphlan_random_donor_patient, method)


            beta_tibble_random <- beta_dist_random %>%
                as.data.frame() %>%
                rownames_to_column("donor_comparison") %>%
                as_tibble() %>%
                filter(row_number() > nstages) %>%
                select(1:(nstages + 1)) %>%
                pivot_longer(-donor_comparison, names_to = "LibID", values_to = "beta_diversity") %>%
                mutate(comparison = "random")
        }





        if (isPatientFMT) {
            areAllDonorsInMetaphlan <- all(donors$LibID %in% t_metaphlan_tib$LibID)
            if (areAllDonorsInMetaphlan == FALSE) {
                missingDonorLib <- donors %>%
                    filter(!LibID %in% t_metaphlan_tib$LibID)

                multipleDonorBatches <- donors %>%
                    group_by(id) %>%
                    filter(n_distinct(donor_batch) > 1)

                n_missingLibs <- length(missingDonorLib$id)
                n_Batches <- length(multipleDonorBatches$id)


                # Check if donor received another batch for comparison.
                if (n_Batches - n_missingLibs == 0) {
                    missingDonorID <- donors %>%
                        filter(id %in% missingDonorLib$id)

                    metaphlan_missingDonor <- t_metaphlan_donors %>%
                        filter(id %in% missingDonorLib$id) %>%
                        column_to_rownames("LibID")

                    metaphlan_actualDonor <- t_metaphlan_tib %>%
                        filter(LibID %in% donors$LibID) %>%
                        column_to_rownames("LibID")

                    metaphlan_donor <- bind_rows(metaphlan_actualDonor, metaphlan_missingDonor)
                } else {
                    metaphlan_donor <- t_metaphlan_tib %>%
                        filter(LibID %in% donors$LibID) %>%
                        column_to_rownames("LibID")
                }
            } else {
                metaphlan_donor <- t_metaphlan_tib %>%
                    filter(LibID %in% donors$LibID) %>%
                    column_to_rownames("LibID")
            }



            metaphlan_donor_patient <- bind_rows(metaphlan_patient, metaphlan_donor) %>%
                abundance_filter(abundance_threshold = abundance_threshold) %>%
                hellinger_transform(transform = transform) %>%
                as.matrix()


            beta_dist_donor <- beta_method(metaphlan_donor_patient, method)

            beta_tibble_donor <- beta_dist_donor %>%
                as.data.frame() %>%
                rownames_to_column("donor_comparison") %>%
                as_tibble() %>%
                filter(row_number() > nstages) %>%
                select(1:(nstages + 1)) %>%
                pivot_longer(-donor_comparison, names_to = "LibID", values_to = "beta_diversity") %>%
                mutate(comparison = "real_donor")

            beta_diversity <- bind_rows(beta_diversity, beta_tibble_donor)
        }


        if (isPatientMP) {
            beta_diversity <- bind_rows(beta_diversity, beta_tibble_random)
        }
    }

    metadata_beta_diversity <- left_join(metadata_patients, beta_diversity, by = "LibID")

    return(metadata_beta_diversity)
}



################################################

#' Betadiversity dissimilarity with metadata
#' #'
#' Function creates a distance matrix, transform the matrix to a dataframe and adds projectmetadata to the dataframe.
#'
#' @param projectmetadata A metadata df that has been processed by `isolateProjectAndDonor`
#' @param t_metaphlan A transposed metahplan df.
#' @param method specifies which beta-distance should be measured. Choose between "bray" or "sorensen". Defaults to bray.
#' @param transform Boolean. Specifies if metaphlan df should be hellinger-transformed before distance calculation
#'
#' @export
#'
#' @example
#' metadata_beta <- createMetadataWBetadiversity(metadata_MP, t_metaphlan, "bray", TRUE)
#'
createMetadataWBetadiversity <- function(projectmetadata, t_metaphlan, method = "bray", transform = TRUE) {
    beta_dist_matrix <- createBetaDistMatrix(t_metaphlan, method, transform)

    beta_dist_l_with_projectmetadata <- beta_dist_matrix %>%
        rownames_to_column("LibID") %>%
        pivot_longer(-LibID, names_to = "comparison", values_to = "dissimilarity") %>%
        left_join(projectmetadata) %>%
        select(id, LibID, comparison, dissimilarity, project, group, -c(batch_1:batch_3), x_axis)

    metadata_beta <- addDonorIdtoBetaDistWMetadata(beta_dist_l_with_projectmetadata, projectmetadata)

    return(metadata_beta)
}


#' createBetaDistMatrix
#' #'
#' creates beta dissimilarity matrix
#'
#' @param t_metaphlan A metaphlan df transposed
#' @param method Either bray or sorensen
#' @param transform Boolean. Should hellinger transformation be applied or not?
#'
#' @example
#'
createBetaDistMatrix <- function(t_metaphlan, method = "bray", transform = T) {
    if (transform) {
        t_metaphlan <- vegan::decostand(t_metaphlan, "hellinger") %>% as.matrix()
        print(paste0("Performing hellinger transformation before beta-diversity calculation!"))
    } else {
        print(paste0("No transformation before beta-diversity calculation!"))
    }

    stopifnot(method %in% c("bray", "sorensen"))
    if (method == "bray") {
        beta_dist <- vegan::vegdist(t_metaphlan)
    }
    if (method == "sorensen") {
        beta_dist <- vegan::vegdist(t_metaphlan, binary = TRUE)
    }

    beta_dist <- beta_dist %>%
        as.matrix() %>%
        as.data.frame()

    return(beta_dist)
}

addDonorIdtoBetaDistWMetadata <- function(beta_distance_matrix_long, projectmetadata) {
    metadata_donors <- projectmetadata %>%
        filter(x_axis == "Donor") %>%
        select(LibID, id, fecal_donation_number) %>%
        rename(comparison = LibID, donor_id = id)

    metadata_beta <- left_join(beta_distance_matrix_long, metadata_donors)

    return(metadata_beta)
}


############################# ############################

#' Get FMT-DONOR comparisons
#' #'
#' extracts the comparison between the FMT group and actual donors received.
#'
#' @param projectmetadata A df with project related metadata created by `isolateProjectMetadata`
#' @param metadata_beta A df created by `createMetadataWBetadiversity`
#'
#' @export
#'
#' @example metadata
#'
compareFMT2ActualDonor <- function(projectmetadata, metadata_beta) {
    metadata_beta_w_donor_received <- createMetadataBetaWDonorReceived(projectmetadata, metadata_beta)

    metadata_beta_actual_donor_for_FMT <- metadata_beta_w_donor_received %>%
        filter(!is.na(donor_id), paste0(donor_id, "_", fecal_donation_number) == donor_batch) %>%
        mutate(actual_donor = "Actual Donors") %>%
        arrange(id)

    return(metadata_beta_actual_donor_for_FMT)
}

#' Get FMT-DONOR_not_received comparisons
#' #'
#' extracts the comparison between the FMT group and donors not received.
#'
#' @param projectmetadata A df with project related metadata created by `isolateProjectMetadata`
#' @param metadata_beta A df created by `createMetadataWBetadiversity`
#'
#' @export
#'
#' @example metadata
#'
compareFMT2DonorNotReceived <- function(projectmetadata, metadata_beta) {
    metadata_beta_w_donor_received <- createMetadataBetaWDonorReceived(projectmetadata, metadata_beta)

    metadata_beta_random_donor_for_FMT <- metadata_beta_w_donor_received %>%
        separate(donor_batch, into = c("donor_id_received", "donor_batch_received"), sep = "_", remove = FALSE) %>%
        filter(donor_id_received != donor_id) %>%
        distinct(LibID, comparison, .keep_all = TRUE) %>%
        select(-c(donor_batch, donor_id_received, donor_batch_received)) %>%
        mutate(actual_donor = "Donors not received") %>%
        arrange(id, donor_id, x_axis, fecal_donation_number)

    set.seed(1)
    distinct_random_donors_for_FMT <- metadata_beta_random_donor_for_FMT %>%
        distinct(id, donor_id, fecal_donation_number, .keep_all = TRUE) %>%
        group_by(id, donor_id) %>%
        filter(row_number() == sample(1:n(), 1)) %>%
        mutate(ptid_donorid_batch = paste0(id, "_", donor_id, "_", fecal_donation_number))


    metadata_beta_random_donor_filtered_for_FMT <- metadata_beta_random_donor_for_FMT %>%
        filter(paste0(id, "_", donor_id, "_", fecal_donation_number) %in% distinct_random_donors_for_FMT$ptid_donorid_batch)

    return(metadata_beta_random_donor_filtered_for_FMT)
}

#' Donor received
#' #'
#' Helper function: Adds the donorbatches the patient received during treatment to metadata_beta
#'
#' @param projectmetadata A df with project related metadata created by `isolateProjectMetadata`
#' @param metadata_beta A df created by `createMetadataWBetadiversity`
#'
#'
#'
#'
createMetadataBetaWDonorReceived <- function(projectmetadata, metadata_beta) {
    donor_batches_received_by_patient <- projectmetadata %>%
        filter(str_detect(id, "pt"), !is.na(donor)) %>%
        pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>%
        filter(!is.na(batch_number)) %>%
        mutate(
            donor_batch = paste0(donor, "_", batch_number)
        ) %>%
        select(id, donor_batch)

    metadata_beta_w_donor_received <- metadata_beta %>%
        left_join(donor_batches_received_by_patient)

    return(metadata_beta_w_donor_received)
}


########## placebo #########
#' Compare placebo group to donor group
#' #'
#' extracts the comparison between the placebo group and donors.
#'
#' @param projectmetadata A df with project related metadata created by `isolateProjectMetadata`
#' @param metadata_beta A df created by `createMetadataWBetadiversity`
#'
#' @export
#'
#' @example
#'
comparePlacebo2DonorNotReceived <- function(projectmetadata, metadata_beta) {
    metadata_donors <- projectmetadata %>%
        filter(x_axis == "Donor") %>%
        select(LibID, id, fecal_donation_number) %>%
        rename(comparison = LibID, donor_id = id)

    distinct_random_donors_for_placebo <- metadata_donors %>%
        group_by(donor_id) %>%
        filter(row_number() == sample(1:n(), 1))

    metadata_beta_placebo <- metadata_beta %>%
        filter(group == "placebo", !is.na(donor_id), comparison %in% distinct_random_donors_for_placebo$comparison) %>%
        mutate(actual_donor = "All Donors")

    return(metadata_beta_placebo)
}


#' Bind comparisons
#' #'
#' A helper function to bind_rows of the compared dataframes 
#'
#'@param projectmetadata A df from `isolateProjectMetadata`
#'@param metadata_beta A df from `createMetadataWBetadiversity`
#'
#'@export
#'
#'@example
#'
createANDcombineComparedOutputs <- function(projectmetadata, metadata_beta){
    metadata_beta_sor_FMT_and_placebo <- bind_rows(
        compareFMT2ActualDonor(projectmetadata, metadata_beta),
        compareFMT2DonorNotReceived(projectmetadata, metadata_beta),
        comparePlacebo2DonorNotReceived(projectmetadata, metadata_beta)
    )
}


beta_method <- function(df, method) {
    if (method == "bray") {
        df <- vegan::vegdist(df) %>%
            as.matrix()
    }
    if (method == "sorensen") {
        df <- vegan::vegdist(df, binary = T) %>%
            as.matrix()
    }
    return(df)
}

abundance_filter <- function(df, abundance_threshold = 0.1, verbose = F) {
    n_before <- ncol(df)

    df_filter <- df %>%
        rownames_to_column("id") %>%
        pivot_longer(-id, names_to = "clade_name", values_to = "relative_abundance") %>%
        mutate(relative_abundance = relative_abundance / 100) %>%
        group_by(clade_name) %>%
        summarise(relative_abundance = max(relative_abundance)) %>%
        filter(relative_abundance >= abundance_threshold / 100)

    df_filtered <- df %>%
        select(df_filter$clade_name)

    n_after <- ncol(df_filtered)
    if (verbose) {
        cat(paste0(
            "Species before: ", n_before, "\n",
            "Species after: ", n_after, "\n",
            "Removed species: ", n_before - n_after, "\n"
        ))
    }

    return(df_filtered)
}

hellinger_transform <- function(df, transform = transform, method = "hellinger") {
    if (transform) {
        df <- vegan::decostand(df, method = method)
    } else {
        df <- df
    }
    return(df)
}



# summariseBetaDiversityOutput <- function(beta_df, donor_id) {
#     metadata_beta_diversity_mean_donor <- beta_df %>%
#         left_join(donor_id) %>%
#         filter(!is.na(donor_comparison))

#     summarise_real_donor <- metadata_beta_diversity_mean_donor %>%
#         filter(comparison == "real_donor") %>%
#         group_by(id, group, x_axis, comparison) %>%
#         summarise(beta_diversity = median(beta_diversity))

#     summarise_random_donor <- metadata_beta_diversity_mean_donor %>%
#         filter(comparison != "real_donor") %>%
#         group_by(id, group, x_axis, comparison) %>%
#         summarise(beta_diversity = median(beta_diversity))


#     metadata_beta_diversity_summarised <- bind_rows(summarise_real_donor, summarise_random_donor) %>%
#         mutate(comparison = if_else(comparison == "random" & group == "FMT", "Donors not received",
#             if_else(comparison == "random" & group == "placebo", "All Donors", "Actual Donors")
#         )) %>%
#         arrange(id)

#     return(metadata_beta_diversity_summarised)
# }

#' Summarise betadiversity
#' #'
#' calculates the median dissimiliarity of each sample to donors received and not received.
#'
#' @param metadata_beta_compared A df created by one of the `compare...` functions
#'
#' @export
#'
#' @example
#'
summariseBetaDiversityOutput <- function(metadata_beta_compared) {
    metadata_beta_summarised <- metadata_beta_compared %>%
        group_by(id, x_axis, actual_donor, group) %>%
        summarise(median_dissimilarity = median(dissimilarity)) %>%
        arrange(id)

    return(metadata_beta_summarised)
}
