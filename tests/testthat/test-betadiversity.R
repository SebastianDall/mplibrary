setwd("../../")

metadata <- load_metadata()
metaphlan <- load_metaphlan("species")

t_metaphlan <- transposeMetaphlan(metaphlan)

metadata_MP <- isolateProjectMetadata(metadata, project_filter = "MP")

test_that("Test that createDistMatrix creates the correct datamatrix", {
    beta_bray_transformed <- createBetaDistMatrix(t_metaphlan, "bray", TRUE)
    beta_bray_transformed_manuel <- vegan::decostand(t_metaphlan, "hellinger") %>%
        as.matrix() %>%
        vegan::vegdist() %>%
        as.matrix() %>%
        as.data.frame()

    beta_bray_untransformed <- createBetaDistMatrix(t_metaphlan, "bray", FALSE)
    beta_bray_untransformed_manuel <- vegan::vegdist(t_metaphlan) %>%
        as.matrix() %>%
        as.data.frame()

    beta_sor_transformed <- createBetaDistMatrix(t_metaphlan, "sorensen", TRUE)
    beta_sor_untransformed <- createBetaDistMatrix(t_metaphlan, "sorensen", FALSE)


    expect_equal(beta_bray_untransformed, beta_bray_untransformed_manuel)
    expect_equal(beta_bray_transformed, beta_bray_transformed_manuel)
    expect_equal(beta_sor_transformed, beta_sor_untransformed)
})

test_that("test that addDonorIDtoBetaDistWMetadata performs correctly", {
    beta_dist_matrix <- createBetaDistMatrix(t_metaphlan, "bray", TRUE)

    beta_dist_l_with_projectmetadata <- beta_dist_matrix %>%
        rownames_to_column("LibID") %>%
        pivot_longer(-LibID, names_to = "comparison", values_to = "dissimilarity") %>%
        left_join(metadata_MP) %>%
        select(id, LibID, comparison, dissimilarity, project, group, -c(batch_1:batch_3), x_axis)

    metadata_beta <- addDonorIdtoBetaDistWMetadata(beta_dist_l_with_projectmetadata, metadata_MP)


    expect_equal(dim(metadata_beta), c(195364, 9))
})

# old function
metadata_beta_old_real_donor <- beta_diversity_v1(metadata_MP, t_metaphlan, method = "bray", transform = TRUE) %>%
    filter(comparison == "real_donor") %>%
    arrange(LibID, donor_comparison)
##

test_that("test that createMetadataWBetadiversity creates metadata_beta", {
    metadata_beta <- createMetadataWBetadiversity(metadata_MP, t_metaphlan, method = "bray", transform = TRUE)

    metadata_beta_real_donor <- metadata_beta %>%
        filter(paste0(LibID, comparison) %in% paste0(metadata_beta_old_real_donor$LibID, metadata_beta_old_real_donor$donor_comparison)) %>%
        arrange(LibID, comparison)

    expect_equal(all(c("id", "LibID", "comparison", "dissimilarity", "donor_id") %in% colnames(metadata_beta)), T)
    expect_equal(metadata_beta_real_donor$dissimilarity, metadata_beta_old_real_donor$beta_diversity, tolerance = 10^-4)
})


metadata_beta <- createMetadataWBetadiversity(metadata_MP, t_metaphlan, "bray", TRUE)

test_that("compareFMT2ActualDonor gives the same results as old beta_diversity function", {
    metadata_beta_FMT_actual <- compareFMT2ActualDonor(metadata_MP, metadata_beta) %>%
        filter(paste0(LibID, comparison) %in% paste0(metadata_beta_old_real_donor$LibID, metadata_beta_old_real_donor$donor_comparison)) %>%
        arrange(LibID, comparison)

    expect_equal(length(metadata_beta_FMT_actual$LibID), length(metadata_beta_old_real_donor$LibID))
    expect_equal(metadata_beta_FMT_actual$dissimilarity, metadata_beta_old_real_donor$beta_diversity, tolerance = 10^-4)
})
