setwd("../..")

metadata <- load_metadata()
metaphlan <- load_metaphlan("species")

t_metaphlan <- transposeMetaphlan(metaphlan)

metadata_MP <- isolateProjectMetadata(metadata, project_filter = "MP")

metadata_beta_bray <- createMetadataWBetadiversity(metadata_MP, t_metaphlan, "bray", TRUE)


metadata_beta_bray_FMT_and_placebo <- bind_rows(
    compareFMT2ActualDonor(metadata_MP, metadata_beta_bray),
    compareFMT2DonorNotReceived(metadata_MP, metadata_beta_bray),
    comparePlacebo2DonorNotReceived(metadata_MP, metadata_beta_bray)
)

metadata_beta_bray_summarised <- metadata_beta_bray_FMT_and_placebo %>%
    group_by(id, x_axis, actual_donor, group) %>%
    summarise(median_dissimilarity = median(dissimilarity)) %>%
    arrange(id)

test_that("gg_beta fails with wrong parameters", {
    expect_error(gg_beta(metadata_beta_bray_summarised, method = "wrong method"), "method is not valid.", fixed = TRUE)
})
