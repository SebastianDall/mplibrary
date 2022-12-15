setwd("../../")

metadata <- load_metadata()
metaphlan <- load_metaphlan("species")

t_metaphlan <- transposeMetaphlan(metaphlan)

metadata_MP <- isolateProjectMetadata(metadata, project_filter = "MP")

test_that("Test that the beta-diversity function gives the correct output", {
    metadata_with_beta <- beta_diversity_v1(metadata_MP, t_metaphlan, transform = T, stages = c("Pre", "Post"), seed = 1)
    metadata_with_beta_ref <- beta_diversity_refactored(metadata_MP, t_metaphlan, transform = T, stages = c("Pre", "Post"), seed = 1)

    expect_equal(metadata_with_beta$beta_diversity[1:3], c(0.974, 0.934, 0.968), tolerance = 0.001)
    expect_equal(metadata_with_beta$beta_diversity[1:3], metadata_with_beta_ref$beta_diversity[1:3])
})


test_that("Test that makeInputDataFrames gives 4 dataframes in a list", {
    dflist <- makeInputDataframes(metadata_MP, t_metaphlan)

    expect_equal(length(dflist), 4)
    expect_equal(dim(dflist[["metadata_patients"]]), c(763, 22))
    expect_equal(dim(dflist[["metadata_donors"]]), c(55, 22))
    expect_equal(dim(dflist[["t_metaphlan"]]), c(442, 505))
    expect_equal(dim(dflist[["t_metaphlan_donors_summarised"]]), c(14, 506))
})
