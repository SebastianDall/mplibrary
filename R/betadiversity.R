
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

#' Betadiversity function
#' #'
#' Calculates the median beta dissimilarity between donor and patient. For patients receiving FMT the beta-diversity is calculated as the median dissimilarity to the donor batch received during treatment. Furthermore, the dissimilarity to donors are also calculated.
#' For placebo patients the similarity is compared to all donors in the dataset.
#'
#' @param $$
#'
#' @export
#'
#' @example
#'
#'
createBetaDiversityDf <- function(projectmetadata, t_taxdata, transform = F, method = "bray", abundance_threshold = 10^-34, stages = c("Pre", "Post"), seed = 1, verbose = F) {
    df <- makeInputDataframes(projectmetadata, t_taxdata) %>%
        calculateBetaDiversity(stages)
}



makeInputDataframes <- function(projectmetadata, t_taxdata) {
    # Make relevant tables
    t_metaphlan_tib <- t_taxdata %>%
        rownames_to_column(var = "LibID")

    metadata_patients <- projectmetadata %>%
        filter(x_axis != "Donor") %>%
        arrange(id)

    metadata_donor <- projectmetadata %>%
        filter(x_axis == "Donor")

    # Summarise the donor relative abundances.
    t_metaphlan_donors_summarised <- t_metaphlan_tib %>%
        left_join(select(metadata_donor, LibID, id), by = "LibID") %>%
        relocate(id) %>%
        select(-LibID) %>%
        group_by(id) %>%
        summarise(across(where(~ is.numeric(.)), ~ mean(.)))

    dataframe_list <- list(
        metadata_patients = metadata_patients,
        metadata_donors = metadata_donor,
        t_metaphlan = t_metaphlan,
        t_metaphlan_donors_summarised = t_metaphlan_donors_summarised
    )
    return(dataframe_list)
}







beta_diversity_refactored <- function(metadata, t_taxdata, transform = F, method = "bray", abundance_threshold = 10^-34, stages = c("Pre", "Post"), seed = 1, verbose = F) {
    stopifnot("method must be either sorensen or bray (default: bray)" = method %in% c("bray", "sorensen"))


    transform <- ifelse(method == "sorensen", F, transform)

    set.seed(seed)

    beta_diversity <- calculateBetaDiversity(metadata_patients, stages, t_metaphlan_tib, t_metaphlan_donors)

    metadata_beta_diversity <- left_join(metadata_patients, beta_diversity, by = "LibID")

    return(metadata_beta_diversity)
}

calculateBetaDiversity <- function(df_list, stages) {
    unique_patients <- unique(df_list[["metadata_patients"]]$id)

    beta_diversity <- tibble()
    for (pt_ID in unique_patients) {
        # patient undergoing analysis
        patient <- df_list[["metadata_patients"]] %>%
            filter(id == pt_ID)

        areStagesPrePost <- all(stages == c("Pre", "Post"))
        onlyCalculateForPatientsWithPrePostData(patient, areStagesPrePost)


        patient_stages <- patient %>%
            filter(x_axis %in% stages)

        metaphlan_patient <- df_list[["t_metaphlan"]] %>%
            filter(LibID %in% patient_stages$LibID) %>%
            column_to_rownames(var = "LibID")

        nstages <- length(rownames(metaphlan_patient))

        isPatientFMT <- patient_stages$group[1] == "FMT"
        if (isPatientFMT) {
            patient_donors <- patient %>%
                pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>%
                filter(!is.na(batch_number)) %>%
                mutate(donor_batch = paste0(donor, "_", batch_number))

            donors <- df_list[["metadata_donors"]] %>%
                filter(donor_batch %in% patient_donors$donor_batch)

            excluded_donors <- df_list[["metadata_donors"]] %>%
                filter(!id %in% donors$id) %>%
                filter(LibID %in% df_list[["t_metaphlan"]]$LibID) %>%
                mutate(id_number = parse_number(id))
        } else {
            excluded_donors <- df_list[["metadata_donors"]] %>%
                filter(LibID %in% df_list[["t_metaphlan"]]$LibID) %>%
                mutate(id_number = parse_number(id))
        }






        set.seed(1)

        isPatientMP <- patient$project[1] == "MP"
        if (isPatientMP) {
            random_donors <- excluded_donors %>%
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
}


onlyCalculateForPatientsWithPrePostData <- function(patient, areStagesPrePost) {
    if (areStagesPrePost) {
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



summariseBetaDiversityOutput <- function(beta_df, donor_id) {
    metadata_beta_diversity_mean_donor <- beta_df %>%
        left_join(donor_id) %>%
        filter(!is.na(donor_comparison))

    summarise_real_donor <- metadata_beta_diversity_mean_donor %>%
        filter(comparison == "real_donor") %>%
        group_by(id, group, x_axis, comparison) %>%
        summarise(beta_diversity = median(beta_diversity))

    summarise_random_donor <- metadata_beta_diversity_mean_donor %>%
        filter(comparison != "real_donor") %>%
        group_by(id, group, x_axis, comparison) %>%
        summarise(beta_diversity = median(beta_diversity))


    metadata_beta_diversity_summarised <- bind_rows(summarise_real_donor, summarise_random_donor) %>%
        mutate(comparison = if_else(comparison == "random" & group == "FMT", "Donors not received",
            if_else(comparison == "random" & group == "placebo", "All Donors", "Actual Donors")
        )) %>%
        arrange(id)

    return(metadata_beta_diversity_summarised)
}
