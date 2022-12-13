setwd("../..")


test_that("Metaphlan3 loads correctly", {
    expect_equal(dim(load_metaphlan("species")), c(505, 444))
    expect_equal(dim(load_metaphlan()), c(186, 444))
    expect_equal(dim(load_metaphlan("species")), c(505, 444))
})

metaphlan <- load_metaphlan("species")


# Metaphlan functions

test_that("calculateSpeciesRichness has the correct output", {
    expect_equal(dim(calculateSpeciesRichness(metaphlan)), c(442, 2))
    expect_equal(colnames(calculateSpeciesRichness(metaphlan)), c("LibID", "richness"))
})


t_metaphlan <- transposeMetaphlan(metaphlan)

test_that("transpose Metaphlan transposes metaphlan", {
    expect_equal(sort(colnames(metaphlan)[-c(1, 2)]), rownames(t_metaphlan))
})
