setwd("../..")


### Loading data
test_that("Test correct metadata is loaded", {
  # Dim
  expect_equal(dim(load_metadata()), c(1277, 20))
  # colnames
  expect_setequal(
    colnames(load_metadata()),
    c("id", "sample_barcode", "LibID", "library_plate", "ext_conc", "lib_conc", "pool_ng", "stage", "project", "group", "donor", "batch_1", "batch_2", "batch_3", "fecal_donation_number", "fecal_batch_date", "pdai_score", "paired_reads", "mean_quality", "non_human_reads")
  )
})

metadata <- load_metadata()



test_that("test if isolateDonorBatchMetadata isolates donor batches", {
  patientMetadataWithDonorbatches <- isolateDonorBatchesUsed(metadata, project_filter = "MP")

  expect_equal(all(patientMetadataWithDonorbatches$project == "MP"), TRUE)
  expect_equal(all(str_detect(patientMetadataWithDonorbatches$donor_batch, "do\\d+")), TRUE) # donor_batch only contains do[number]
})


test_that("Test if isolateDonorAndPatientMetadata isolates metadata for both Patients and Donors", {
  patientMetadataWithDonorbatches <- isolateDonorBatchesUsed(metadata, project_filter = "MP")

  DonorAndPatientMetadata <- isolateDonorAndPatientMetadata(metadata, metadata_with_donor_batches_used = patientMetadataWithDonorbatches, project_filter = "MP")

  expect_equal(unique(DonorAndPatientMetadata$project), c("donor_batch", "MP"))
  expect_equal(unique(DonorAndPatientMetadata$group), c("FMT", "placebo"))
})

test_that("Test if createXaxis gives the correct Xaxis column values", {
  # setup
  patientMetadataWithDonorbatches <- isolateDonorBatchesUsed(metadata, project_filter = "MP")
  DonorAndPatientMetadata <- isolateDonorAndPatientMetadata(metadata, metadata_with_donor_batches_used = patientMetadataWithDonorbatches, project_filter = "MP")

  # dataframe to be tested
  isolatedMetadatWithXaxis <- DonorAndPatientMetadata %>%
    createXaxis()

  # test
  expect_equal(
    unique(isolatedMetadatWithXaxis$x_axis),
    c(
      "Donor", "Pre", "Post", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "followup_1m",
      "followup_6m", "followup_12m", "drop_out", "followup_3m", "treatment_1", "treatment_2", "treatment_3", "treatment_4", "treatment_6", "treatment_7", "treatment_8",
      "treatment_9", "treatment_11", "treatment_12", "treatment_13", "treatment_14", "followup_14d", "treatment_16", "treatment_17", "treatment_18", "treatment_19", "treatment_20", "end_of_study"
    )
  )
})

test_that("Test that isolateProjectMetadata works with other functions.", {
  # setup
  patientMetadataWithDonorbatches <- isolateDonorBatchesUsed(metadata, project_filter = "MP")
  DonorAndPatientMetadata <- isolateDonorAndPatientMetadata(metadata, metadata_with_donor_batches_used = patientMetadataWithDonorbatches, project_filter = "MP") %>%
    createXaxis()

  # dataframe to be tested
  isolated_metadata <- isolateProjectMetadata(metadata, project_filter = "MP")

  expect_equal(all(isolated_metadata == DonorAndPatientMetadata, na.rm = T), T)
})


## Metaphlan functions

# test_that("calculateSpeciesRichness has the correct output", {
#   expect_equivalent(dim(calculateSpeciesRichness(metaphlan)), c(442, 2))
#   expect_equivalent(colnames(calculateSpeciesRichness(metaphlan)), c("LibID", "richness"))
# })


# t_metaphlan <- transposeMetaphlan(metaphlan)

# test_that("transpose Metaphlan transposes metaphlan", {
#   expect_equivalent(sort(colnames(metaphlan)[-c(1, 2)]), rownames(t_metaphlan))
# })
