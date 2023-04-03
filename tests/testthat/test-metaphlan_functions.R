setwd("../../src/")


test_that("Metaphlan3 loads correctly", {
    expect_equal(dim(load_metaphlan("species")), c(505, 444))
    expect_equal(dim(load_metaphlan()), c(186, 444))
    expect_equal(dim(load_metaphlan("species")), c(505, 444))
})

# metaphlan <- load_metaphlan("species")
metaphlan <- readr::read_delim(
    "../data/metaphlan3/MetaPhlAn_4.0.3_Combined_NonHuman_Subsampled_2500000_profile.txt", 
    delim = "\t", 
    show_col_types = FALSE, 
    skip = 1
)

test_that("Filter species gives expected dim", {
    expect_equal(dim(filter_taxonomy(metaphlan, "species", TRUE)), c(1306, 646))
    expect_equal(dim(filter_taxonomy(metaphlan, "species", FALSE)), c(1305, 646))
})

test_that("Filter species gives expected sum", {
    expect_equal(sum(filter_taxonomy(metaphlan, "species", TRUE)[,2]), 100, tolerance = 0.001)
    expect_equal(sum(filter_taxonomy(metaphlan, "species", FALSE)[,2]) == 100, FALSE)
})

test_that("remove_NonHuman_from_colnames - removes NonHuman", {
    expect_equal(all(is.na(str_match(colnames(remove_NonHuman_from_colnames(metaphlan)), "NonHuman"))), TRUE)
})


metaphlan <- metaphlan %>% 
    filter_taxonomy("species", TRUE) %>% 
    remove_NonHuman_from_colnames()

# Metaphlan functions

test_that("calculateSpeciesRichness has the correct output", {
    s <- calculateSpeciesRichness(metaphlan)
    expect_equal(s$richness[1], sum(metaphlan[s$sample_barcode[1]]>0) -1)
    
})


t_metaphlan <- transposeMetaphlan(metaphlan)

test_that("transpose Metaphlan transposes metaphlan", {
    expect_equal(sort(colnames(metaphlan)[-1]), rownames(t_metaphlan))
})
