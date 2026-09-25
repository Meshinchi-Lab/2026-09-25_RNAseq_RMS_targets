root <- rprojroot::find_root(rprojroot::has_file("config.yml"))
source(file.path(root, "R", "harmonise.R"))
source(file.path(root, "R", "rms_palette.R"))
source(file.path(root, "R", "plot_gene_distribution.R"))

fixture <- function() {
  set.seed(1)
  data.frame(
    gene_symbol = "GAPDH",
    dataset  = rep(c("PCGP", "dGTEx"), each = 20),
    sample_id = paste0("S", 1:40),
    value    = c(rlnorm(20, 5, 1), rlnorm(20, 4, 1)),
    group    = rep(c("FP", "healthy_muscle"), each = 20),
    fusion   = rep(c("FP", "normal"), each = 20),
    tissue   = rep(c(NA, "Muscle"), each = 20),
    subgroup = NA_character_,
    assay    = "pTPM",
    unit     = "pTPM",
    stringsAsFactors = FALSE)
}

test_that("returns a ggplot for each style", {
  for (s in c("violin", "box", "ridge")) {
    p <- plot_gene_distribution(fixture(), gene = "GAPDH", style = s)
    expect_s3_class(p, "ggplot")
  }
})

test_that("log2 scaling changes the y label and the data", {
  lin <- plot_gene_distribution(fixture(), "GAPDH", scale = "linear")
  lg  <- plot_gene_distribution(fixture(), "GAPDH", scale = "log2")
  expect_match(lin$labels$y, "TPM", fixed = TRUE)
  expect_match(lg$labels$y, "log2", fixed = TRUE)
  expect_false(isTRUE(all.equal(lin$data$.value, lg$data$.value)))
})

test_that("interactive = TRUE swaps in an interactive point layer", {
  st <- plot_gene_distribution(fixture(), "GAPDH", interactive = FALSE)
  iv <- plot_gene_distribution(fixture(), "GAPDH", interactive = TRUE)
  cls <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false(any(grepl("Interactive", cls(st))))
  expect_true(any(grepl("Interactive", cls(iv))))
})

test_that("static and interactive forms share title and axis labels", {
  st <- plot_gene_distribution(fixture(), "GAPDH", interactive = FALSE)
  iv <- plot_gene_distribution(fixture(), "GAPDH", interactive = TRUE)
  expect_equal(st$labels$title, iv$labels$title)
  expect_equal(st$labels$y, iv$labels$y)
})

test_that("an empty frame yields a placeholder plot rather than an error", {
  p <- plot_gene_distribution(fixture()[0, ], "NOSUCHGENE")
  expect_s3_class(p, "ggplot")
})

test_that("grouping by tissue falls back to the sequential ramp for free-text levels", {
  d <- fixture()
  d$tissue <- rep(c("Muscle", "Brain"), each = 20)
  p <- plot_gene_distribution(d, "GAPDH", group_by = "tissue")
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("every plot builds without warnings", {
  expect_silent(ggplot2::ggplot_build(plot_gene_distribution(fixture(), "GAPDH")))
})

## ------------------------------------------------- regressions from the 2026-09-25 review

mixed_fixture <- function() {
  set.seed(2)
  rbind(
    transform(fixture(), dataset = "PCGP", assay = "pTPM", unit = "pTPM",
              value = rlnorm(40, 5, 1)),
    data.frame(gene_symbol = "GAPDH", dataset = "DepMap_26Q1",
               sample_id = paste0("CL", 1:10), value = runif(10, 2, 12),
               group = "cell_line_RMS_FP", fusion = "FP", tissue = NA_character_,
               subgroup = NA_character_, assay = "log2_TPM1p",
               unit = "log2(TPM+1)", stringsAsFactors = FALSE))
}

test_that("mixing pTPM and DepMap log2 never yields a nested-log axis", {
  for (sc in c("log2", "linear")) {
    p <- plot_gene_distribution(mixed_fixture(), "GAPDH", scale = sc)
    expect_false(grepl("log2\\(log2", p$labels$y))
  }
  expect_equal(plot_gene_distribution(mixed_fixture(), "GAPDH", scale = "log2")$labels$y,
               "log2(TPM + 1)")
})

test_that("DepMap values pass through untouched on the log2 axis", {
  p <- plot_gene_distribution(mixed_fixture(), "GAPDH", scale = "log2")
  dm <- p$data[p$data$dataset == "DepMap_26Q1", ]
  expect_equal(dm$.value, dm$value)
})

test_that("group labels are humanised on the axis", {
  p <- plot_gene_distribution(fixture(), "GAPDH")
  expect_true("healthy muscle" %in% levels(p$data$.group))
  expect_false("healthy_muscle" %in% levels(p$data$.group))
})

test_that("order_levels fixes the x order and drops unlisted groups", {
  p <- plot_gene_distribution(fixture(), "GAPDH",
                              order_levels = c("FP", "healthy_muscle"))
  expect_equal(levels(p$data$.group), c("FP", "healthy muscle"))
})

test_that("tooltips carry fusion status when grouping is not fusion", {
  p <- plot_gene_distribution(fixture(), "GAPDH", group_by = "group",
                              interactive = TRUE)
  expect_true(any(grepl("fusion:", p$data$.tip)))
})

test_that("every palette entry resolves for the full group vocabulary", {
  pal <- rms_group_palette(label_groups(VALID_GROUPS))
  expect_equal(length(pal), length(VALID_GROUPS))
  expect_false(anyNA(pal))
  ## PDX is a genuinely third nominal category, so it must not share FP's or FN's hue.
  lab <- function(g) unname(GROUP_LABELS[g])
  expect_false(pal[[lab("PDX")]] == pal[[lab("FP_cellline")]])
  expect_false(pal[[lab("PDX")]] == pal[[lab("FN_cellline")]])
})

## ------------------------------------- small-n groups must stay visible (2026-09-25 rev2)

small_n_fixture <- function() {
  data.frame(
    gene_symbol = "X",
    dataset  = "Wachtel_protein",
    sample_id = c(paste0("FP", 1:4), paste0("FN", 1:6), paste0("P", 1:3),
                  "MRC-5", "Myoblasts"),
    value    = c(runif(4, 1e7, 5e8), runif(6, 1e6, 4e8), runif(3, 1e7, 5e8),
                 2.4e8, 2.5e7),
    group    = c(rep("FP_cellline", 4), rep("FN_cellline", 6),
                 rep("PDX", 3), rep("control", 2)),
    fusion   = c(rep("FP", 4), rep("FN", 6), rep("FP", 3), rep("control", 2)),
    tissue   = NA_character_, subgroup = NA_character_,
    assay    = "iTop3", unit = "iTop3 abundance",
    stringsAsFactors = FALSE)
}

test_that("every sample reaches the plot, controls included", {
  p <- plot_gene_distribution(small_n_fixture(), "X",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  expect_equal(nrow(p$data), 15L)
  expect_true(all(c("MRC-5", "Myoblasts") %in% p$data$sample_id))
  expect_equal(sum(p$data$.group == "control"), 2L)
})

test_that("a violin is not drawn for groups below MIN_N_DISTRIBUTION", {
  ## The n=2 control group must contribute points but no violin: a density over two
  ## observations renders as a hairline and makes a present group look absent.
  b <- ggplot2::ggplot_build(
    plot_gene_distribution(small_n_fixture(), "X", style = "violin",
      order_levels = c("FP_cellline", "FN_cellline", "PDX", "control")))
  violin_groups <- length(unique(b$data[[1]]$group))
  expect_equal(violin_groups, 3L)          # FP=4, FN=6, PDX=3 qualify; control=2 does not
  expect_equal(nrow(b$data[[2]]), 15L)     # all points, including both controls
})

test_that("an n=1 group still plots its point", {
  d <- small_n_fixture()
  d <- d[d$sample_id != "Myoblasts", ]     # leaves a single control
  p <- plot_gene_distribution(d, "X",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  expect_equal(sum(p$data$.group == "control"), 1L)
  expect_true("control" %in% levels(p$data$.group))
})

test_that("points are enlarged for small panels and kept small for large ones", {
  ## The threshold is 60 rows; mixed_fixture() is only 50, so build a wider frame.
  wide <- do.call(rbind, replicate(10, small_n_fixture(), simplify = FALSE))
  wide$sample_id <- paste0("s", seq_len(nrow(wide)))
  small <- plot_gene_distribution(small_n_fixture(), "X")   # 15 rows
  big   <- plot_gene_distribution(wide, "X")                # 150 rows
  pt <- function(p) p$layers[[2]]$aes_params$size
  expect_equal(pt(small), 2.1)
  expect_equal(pt(big), 0.8)
  expect_gt(pt(small), pt(big))
})

test_that("the .n_group helper counts within facet, not across the panel", {
  d <- rbind(small_n_fixture(),
             transform(small_n_fixture(), dataset = "Other",
                       sample_id = paste0("o", 1:15)))
  p <- plot_gene_distribution(d, "X")
  ctl <- unique(p$data$.n_group[p$data$.group == "control"])
  expect_equal(ctl, 2L)   # 2 per dataset, not 4 pooled
})

## ---------- not-detected convention, per 2026-08-06_RMS_protein_validation.qmd ----------

undetected_fixture <- function() {
  d <- small_n_fixture()
  ## Mirror FGFR4: RD, PDXs29 and Myoblasts below the MS floor.
  d$sample_id[c(10, 13, 15)] <- c("RD", "PDXs29", "Myoblasts")
  d$value[c(10, 13, 15)] <- NA_real_
  d
}

test_that("non-detections are kept, not dropped", {
  p <- plot_gene_distribution(undetected_fixture(), "FGFR4",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  expect_equal(nrow(p$data), 15L)
  expect_equal(sum(!p$data$.detected), 3L)
  expect_true(all(c("RD", "PDXs29", "Myoblasts") %in%
                  p$data$sample_id[!p$data$.detected]))
})

test_that("the sentinel is 30% of the minimum detected value", {
  d <- undetected_fixture()
  p <- plot_gene_distribution(d, "FGFR4", scale = "linear",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  min_det <- min(d$value, na.rm = TRUE)
  expect_equal(unique(p$data$.value[!p$data$.detected]), min_det * 0.30)
  ## The band sits strictly below every real value.
  expect_lt(max(p$data$.value[!p$data$.detected]),
            min(p$data$.value[p$data$.detected]))
})

test_that("shape encodes detection with pch 21 filled / pch 2 open", {
  ## Assert on the BUILT layer rather than the scale object: ggplot2 4.x keeps
  ## scale_shape_manual's values inside the palette closure, where they are not readable.
  p <- plot_gene_distribution(undetected_fixture(), "FGFR4",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  b <- ggplot2::ggplot_build(p)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  pts <- b$data[[which(geoms == "GeomPoint")]]
  ## pch 21 = filled circle with outline (obeys fill); pch 2 = open triangle (ignores it).
  expect_setequal(unique(pts$shape), c(21L, 2L))
  expect_equal(sum(pts$shape == 2L), 3L)    # RD, PDXs29, Myoblasts
  expect_equal(sum(pts$shape == 21L), 12L)
})

test_that("summaries exclude non-detections", {
  d <- undetected_fixture()
  p <- plot_gene_distribution(d, "FGFR4", style = "violin",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  b <- ggplot2::ggplot_build(p)
  ## control has 1 detected of 2 and PDX 2 of 3, so neither reaches
  ## MIN_N_DISTRIBUTION = 3; FP (4/4) and FN (5/6) do.
  expect_equal(length(unique(b$data[[1]]$group)), 2L)

  ## The median crossbar must sit inside the detected range, never dragged to the
  ## sentinel band by a non-detection.
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  cross <- b$data[[which(geoms == "GeomCrossbar")]]
  floor_y <- unique(p$data$.value[!p$data$.detected])
  expect_true(all(cross$y > floor_y))
  expect_true(all(cross$y >= min(p$data$.value[p$data$.detected])))
})

test_that("a separator line and one 'not detected' annotation are added, and nothing else", {
  p <- plot_gene_distribution(undetected_fixture(), "FGFR4",
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomHline" %in% geoms)
  expect_equal(sum(geoms == "GeomText"), 1L)
  ## The per-gene "k/n NA" count annotation from the notebook is deliberately NOT carried
  ## over: it belongs in the summary table, not on the figure.
  txt <- p$layers[[which(geoms == "GeomText")]]
  expect_match(as.character(txt$aes_params$label %||% ""), "not detected")
})

test_that("panels with nothing missing keep the plain dark points", {
  p <- plot_gene_distribution(small_n_fixture(), "X")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomHline" %in% geoms)
  expect_equal(sum(geoms == "GeomText"), 0L)
  expect_equal(p$layers[[2]]$aes_params$colour, "grey20")
})

test_that("tooltips mark non-detections explicitly", {
  p <- plot_gene_distribution(undetected_fixture(), "FGFR4", interactive = TRUE,
         order_levels = c("FP_cellline", "FN_cellline", "PDX", "control"))
  expect_true(all(grepl("not detected", p$data$.tip[!p$data$.detected])))
  expect_false(any(grepl("not detected", p$data$.tip[p$data$.detected])))
})
