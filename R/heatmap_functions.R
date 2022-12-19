

#' Hellinger transform metadata with taxonomic data
#'
#' Helper function in generating the y-axis for the heatmap. Each sample batch is hellinger transformed
#'
#' @param projectmetadata_w_tax A df with metadata and taxonomic data.
#'
#' @export
#'
#' @example
#'
hellingerTransformSamples <- function(projectmetadata_w_tax) {
    # Defining groups used in weighted_mean - also relative abundance are hellinger-transformed.
    df <- projectmetadata_w_tax %>%
        group_by(interaction(id, fecal_donation_number)) %>%
        mutate(
            hellinger = sqrt(value / sum(value)),
            sum = sum(value)
        ) %>%
        ungroup()

    return(df)
}

#' Finding top genera
#'
#' finds the top genera based on weighted mean of hellinger transformed relative abundances.
#'
#' @param projectmetadata_w_tax_hellinger_transformed A df created by `hellingerTransformSamples`
#' @param stages A vector with stages which weighted mean will be based on. Defaults to `c("Pre", "Donor")`
#' @param n_genera Top n genera.
#'
#' @export
#'
#' @example
#'
calculateTopGenera <- function(projectmetadata_w_tax_hellinger_transformed, stages = c("Pre", "Donor"), n_genera) {
    # top genera as a function of weighted mean between pre and donor

    df <- projectmetadata_w_tax_hellinger_transformed %>%
        filter(x_axis %in% stages) %>%
        group_by(id, x_axis, clade_name) %>%
        summarise(hellinger = mean(hellinger)) %>%
        group_by(x_axis, clade_name) %>%
        summarise(
            weight = n(),
            mean_hellinger = mean(hellinger)
        ) %>%
        group_by(clade_name) %>%
        summarise(weighted_mean_hellinger = weighted.mean(mean_hellinger, weight)) %>%
        arrange(desc(weighted_mean_hellinger)) %>%
        head(n_genera)

    return(df)
}

#' Top genera ranked based on the Donor group
#'
#' Arrange the top genera based on the Donor group
#'
#' @param projectmetadata_w_tax_hellinger_transformed A df created by `hellingerTransformSamples`
#' @param top_genera A df created by `calculateTopGenera`
#'
#' @export
#'
#' @example
#'
arrangeTopGeneraInDonor <- function(projectmetadata_w_tax_hellinger_transformed, top_genera) {
    df <- projectmetadata_w_tax_hellinger_transformed %>%
        filter(clade_name %in% top_genera$clade_name, x_axis == "Donor") %>%
        group_by(id, clade_name) %>%
        summarise(hellinger = mean(hellinger)) %>%
        group_by(clade_name) %>%
        summarise(mean_hellinger = mean(hellinger)) %>%
        arrange(desc(mean_hellinger))

    return(df)
}


#' Cluster samples for heatmap
#'
#' calculates clusters based on hellinger transformed bray-curtis distance with the ward.D2 algorithm.
#'
#' @param metadata_tax_hellinger A df made by `hellingerTransformSamples`
#' @param t_metaphlan A metaphlan df transposed
#'
#' @export
#'
#' @example
#'
clusterPreAndDonorSamples <- function(metadata_tax_hellinger, t_metaphlan) {
    t_metaphlan_id <- metadata_tax_hellinger %>%
        filter(str_detect(sample_group, "Pre") | sample_group == "Donor") %>%
        arrange(LibID) %>%
        select(id, LibID, group, sample_group, x_axis) %>%
        distinct(LibID, .keep_all = T) %>%
        left_join(
            rownames_to_column(t_metaphlan, "LibID")
        )

    t_metaphlan_long_id_donor_avg <- t_metaphlan_id %>%
        pivot_longer(!id:x_axis, names_to = "clade_name", values_to = "relative_abundance") %>%
        group_by(id, x_axis, clade_name) %>%
        summarise(relative_abundance = mean(relative_abundance))

    t_metaphlan_wide_id_donor_avg <- t_metaphlan_long_id_donor_avg %>%
        ungroup() %>%
        select(!x_axis) %>%
        pivot_wider(names_from = clade_name, values_from = relative_abundance) %>%
        column_to_rownames("id")

    bray_curtis_dist <- vegan::vegdist(vegan::decostand(t_metaphlan_wide_id_donor_avg, method = "hellinger"))

    hclust_ward <- hclust(bray_curtis_dist, method = "ward.D2")

    ward_dendogram <- as.dendrogram(hclust_ward)
    ward_order <- order.dendrogram(ward_dendogram)


    id_cluster_order <- t_metaphlan_long_id_donor_avg %>%
        left_join(metadata_tax_hellinger) %>%
        mutate(
            id = factor(
                id,
                levels = hclust_ward$labels[ward_order]
            )
        ) %>%
        group_by(id, clade_name, sample_group) %>%
        distinct(id, .keep_all = T)

    return(id_cluster_order)
}
