setwd("../..")
library(tidyverse)

metadata <- load_metadata()

metaphlan <- load_metaphlan(taxonomic_level = "genus")

t_metaphlan <- transposeMetaphlan(metaphlan)

project_filter <- "MP"

metadata_in_metaphlan <- metadata %>%
    filter(LibID %in% colnames(metaphlan))


metaphlan_long <- metaphlan %>%
    select(-NCBI_tax_id) %>%
    pivot_longer(-clade_name, names_to = "LibID") %>%
    arrange(LibID)

stages_pre_post <- c("Pre", "Post", "Donor")

stages <- stages_pre_post


metadata_relevant <- isolateProjectMetadata(metadata_in_metaphlan, project_filter = project_filter) %>%
    mutate(
        sample_group = paste0(x_axis, " ", group),
        sample_group = if_else(str_detect(sample_group, "Donor FMT"), str_remove(sample_group, " FMT"), sample_group)
    ) %>%
    filter(x_axis %in% stages) %>%
    group_by(id) %>%
    filter(n_distinct(x_axis) == 2 | x_axis == "Donor")


metadata_relevant_with_tax <- left_join(metadata_relevant, metaphlan_long)
n_genera <- 25

df <- hellingerTransformSamples(metadata_relevant_with_tax)
df_id1 <- df %>%
    filter(id == "pt001", x_axis == "Pre")


test_that("hellingerTransformSamples: the sum of hellinger transformed relative abundance for each id is 1", {
    df <- hellingerTransformSamples(metadata_relevant_with_tax)

    df_sum <- df %>%
        group_by(interaction(id, fecal_donation_number)) %>%
        summarise(
            value = sum(value),
            hellinger = sum(hellinger)
        )

    expect_equal(df_sum$value, rep(100, length(df_sum$value)), tolerance = 10^-3)
})

metadata_tax_hellinger <- hellingerTransformSamples(metadata_relevant_with_tax)


test_that("calculateTopGenera gives the correct output", {
    df <- calculateTopGenera(metadata_tax_hellinger, n_genera = n_genera)

    df2 <- calculateTopGenera(metadata_tax_hellinger, stages = c("Post", "Donor"), n_genera = n_genera)

    expect_equal(length(df$clade_name), n_genera)
    expect_equal(df$clade_name[1:3], c("Blautia", "Escherichia", "Bifidobacterium"))
    expect_equal(df$weighted_mean_hellinger[1:3], c(0.272, 0.236, 0.225), tolerance = 10^-3)
    expect_equal(all(df$clade_name == df2$clade_name), FALSE)
})


top_genera_weighted <- calculateTopGenera(metadata_tax_hellinger, n_genera = n_genera)

test_that("arrangeTopGenera gives the correct output", {
    df <- arrangeTopGeneraInDonor(metadata_tax_hellinger, top_genera_weighted)

    expect_equal(all(df$clade_name %in% top_genera_weighted$clade_name), TRUE)
    expect_equal(all(df$clade_name == top_genera_weighted$clade_name), FALSE)
    expect_equal(length(df$clade_name) == length(top_genera_weighted$clade_name), TRUE)
})



test_that("clusterPreAndDonorSamples gives the correct output", {
    metadata_tax_clustered <- clusterPreAndDonorSamples(metadata_tax_hellinger, t_metaphlan)

    expect_equal(levels(filter(metadata_tax_clustered, group == "FMT")$id)[c(1, 2, 34)], c("pt014", "do014", "pt017"))
})
