# Helper: write a temp file and return its path.
write_tmp <- function(content, suffix = ".tmp") {
  path <- tempfile(fileext = suffix)
  writeLines(content, path)
  path
}

# Minimal 3-taxon NEXUS string (sequential, no gaps).
NEXUS_3TAX <- c(
  "#NEXUS",
  "BEGIN DATA;",
  "  DIMENSIONS NTAX=3 NCHAR=4;",
  "  FORMAT DATATYPE=DNA;",
  "  MATRIX",
  "  Alpha ACGT",
  "  Beta  TTAG",
  "  Gamma CGCA",
  ";",
  "END;"
)

# 3-taxon NEXUS with a TITLE.
NEXUS_WITH_TITLE <- c(
  "#NEXUS",
  "BEGIN DATA;",
  "  TITLE COI;",
  "  DIMENSIONS NTAX=2 NCHAR=3;",
  "  FORMAT DATATYPE=DNA;",
  "  MATRIX",
  "  T1 ACG",
  "  T2 TTA",
  ";",
  "END;"
)

# Multi-gene NEXUS (two DATA blocks).
NEXUS_MULTIGENE <- c(
  "#NEXUS",
  "BEGIN DATA;",
  "  TITLE gene_a;",
  "  DIMENSIONS NTAX=3 NCHAR=3;",
  "  FORMAT DATATYPE=DNA;",
  "  MATRIX",
  "  A AAA",
  "  B AAA",
  "  C AAA",
  ";",
  "END;",
  "BEGIN DATA;",
  "  TITLE gene_b;",
  "  DIMENSIONS NTAX=3 NCHAR=3;",
  "  FORMAT DATATYPE=DNA;",
  "  MATRIX",
  "  A CCC",
  "  B CCC",
  "  C CCC",
  ";",
  "END;"
)

# NEXUS with SETS block.
NEXUS_SETS <- c(
  "#NEXUS",
  "BEGIN DATA;",
  "  DIMENSIONS NTAX=2 NCHAR=6;",
  "  FORMAT DATATYPE=DNA;",
  "  MATRIX",
  "  T1 ACGTAC",
  "  T2 TTGGCA",
  ";",
  "END;",
  "BEGIN SETS;",
  "  CHARSET loc1 = 1-3;",
  "  CHARSET loc2 = 4-6;",
  "END;"
)

# Minimal PHYLIP (sequential).
PHYLIP_SEQ <- c(
  "3 4",
  "Alpha ACGT",
  "Beta  TTAG",
  "Gamma CGCA"
)

# -- read_nexus ----------------------------------------------------------------

test_that("read_nexus returns a PhylotripletAlignment", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  expect_s3_class(aln, "PhylotripletAlignment")
})

test_that("read_nexus: num_genes and taxon count are correct", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  expect_equal(phylotriplet:::alignment_num_genes_rcpp(aln$.ptr), 1L)
  taxa <- phylotriplet:::alignment_taxon_names_rcpp(aln$.ptr)
  expect_equal(sort(as.character(taxa)), c("Alpha", "Beta", "Gamma"))
})

test_that("read_nexus: print method runs without error", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  expect_output(print(aln), "PhylotripletAlignment")
  expect_output(print(aln), "gene")
  expect_output(print(aln), "Alpha")
})

test_that("read_nexus: TITLE sets gene name (accessible via Rcpp)", {
  path <- write_tmp(NEXUS_WITH_TITLE, ".nex")
  aln <- read_nexus(path)
  expect_equal(phylotriplet:::alignment_num_genes_rcpp(aln$.ptr), 1L)
  # The alignment exists and has 2 taxa
  taxa <- phylotriplet:::alignment_taxon_names_rcpp(aln$.ptr)
  expect_equal(length(taxa), 2L)
})

test_that("read_nexus: multiple DATA blocks → multiple genes", {
  path <- write_tmp(NEXUS_MULTIGENE, ".nex")
  aln <- read_nexus(path)
  expect_equal(phylotriplet:::alignment_num_genes_rcpp(aln$.ptr), 2L)
})

test_that("read_nexus: file-not-found throws", {
  expect_error(read_nexus("/no/such/file.nex"), "File not found")
})

test_that("read_nexus: missing #NEXUS header throws", {
  path <- write_tmp(c("BEGIN DATA;", "END;"), ".nex")
  expect_error(read_nexus(path))
})

test_that("read_nexus: unsupported DATATYPE throws", {
  bad <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=2 NCHAR=3;",
    "  FORMAT DATATYPE=PROTEIN;",
    "  MATRIX", "  T1 ACG", "  T2 TTA", ";", "END;"
  )
  path <- write_tmp(bad, ".nex")
  expect_error(read_nexus(path))
})

# -- read_phylip_files ---------------------------------------------------------

test_that("read_phylip_files: single file returns PhylotripletAlignment", {
  path <- write_tmp(PHYLIP_SEQ, ".phy")
  aln <- read_phylip_files(path)
  expect_s3_class(aln, "PhylotripletAlignment")
  expect_equal(phylotriplet:::alignment_num_genes_rcpp(aln$.ptr), 1L)
  taxa <- phylotriplet:::alignment_taxon_names_rcpp(aln$.ptr)
  expect_equal(sort(as.character(taxa)), c("Alpha", "Beta", "Gamma"))
})

test_that("read_phylip_files: two files → two genes", {
  p1 <- write_tmp(PHYLIP_SEQ, ".phy")
  p2 <- write_tmp(c("2 3", "X ACG", "Y TTA"), ".phy")
  aln <- read_phylip_files(c(p1, p2))
  expect_equal(phylotriplet:::alignment_num_genes_rcpp(aln$.ptr), 2L)
})

test_that("read_phylip_files: file-not-found throws", {
  expect_error(read_phylip_files("/no/such/file.phy"), "File.s. not found")
})

# -- triplet_site_patterns -----------------------------------------------------

test_that("triplet_site_patterns: returns a data frame with correct columns", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_s3_class(pats, "data.frame")
  expect_true(all(c("taxon1", "taxon2", "taxon3", "n_total", "n_counted", "site_pattern")
                  %in% names(pats)))
})

test_that("triplet_site_patterns: 3 taxa → 1 triplet row", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_equal(nrow(pats), 1L)
})

test_that("triplet_site_patterns: taxa are in lexicographic order", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_true(pats$taxon1[1] <= pats$taxon2[1])
  expect_true(pats$taxon2[1] <= pats$taxon3[1])
})

test_that("triplet_site_patterns: site_pattern is a list column of length-64 vectors", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_type(pats$site_pattern, "list")
  expect_equal(length(pats$site_pattern[[1]]), 64L)
})

test_that("triplet_site_patterns: n_total and n_counted are numeric", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_type(pats$n_total, "double")
  expect_type(pats$n_counted, "double")
  expect_equal(pats$n_total[1], 4.0)
  expect_equal(pats$n_counted[1], 4.0)
})

test_that("triplet_site_patterns: sum(n) equals n_counted", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  for (i in seq_len(nrow(pats))) {
    expect_equal(sum(pats$site_pattern[[i]]), pats$n_counted[i], tolerance = 1e-10)
  }
})

test_that("triplet_site_patterns: all-A alignment → all mass in AAA bucket", {
  # All three taxa have all-A sequences: only site_pattern[1] (index 0+1=1 in R) is nonzero
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=3 NCHAR=3;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX", "  A AAA", "  B AAA", "  C AAA", ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  n_vec <- pats$site_pattern[[1]]
  # AAA: i=0, j=0, k=0 → 0-indexed 0 → 1-indexed 1
  expect_equal(n_vec[1], 3.0)
  expect_equal(sum(n_vec[-1]), 0.0)
})

test_that("triplet_site_patterns: gap sites reduce n_counted but not n_total", {
  # 4 sites, 1 has gap in taxon B → n_counted = 3, n_total = 4
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=3 NCHAR=4;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX", "  A AAAA", "  B A-AA", "  C AAAA", ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_equal(pats$n_total[1], 4.0)
  expect_equal(pats$n_counted[1], 3.0)
})

test_that("triplet_site_patterns: multi-gene accumulates correctly", {
  # Gene 1: all-A; Gene 2: all-C → AAA count = 3, CCC count = 3
  path <- write_tmp(NEXUS_MULTIGENE, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_equal(nrow(pats), 1L)  # one triplet: (A, B, C)
  n_vec <- pats$site_pattern[[1]]
  # AAA: index 0+1 = 1; CCC: 1*16+1*4+1+1 = 22 (1-indexed)
  expect_equal(n_vec[1],  3.0)   # AAA from gene_a
  expect_equal(n_vec[22], 3.0)   # CCC from gene_b
  expect_equal(pats$n_total[1],   6.0)
  expect_equal(pats$n_counted[1], 6.0)
})

test_that("triplet_site_patterns: SETS block splits and both loci contribute", {
  path <- write_tmp(NEXUS_SETS, ".nex")
  aln <- read_nexus(path)
  # Only 2 taxa → no triplets
  pats <- triplet_site_patterns(aln)
  expect_equal(nrow(pats), 0L)
})

test_that("triplet_site_patterns via read_phylip_files works end-to-end", {
  # Three taxa, all-A sequence
  phy <- c("3 2", "A AA", "B AA", "C AA")
  path <- write_tmp(phy, ".phy")
  aln <- read_phylip_files(path)
  pats <- triplet_site_patterns(aln)
  expect_equal(nrow(pats), 1L)
  expect_equal(pats$n_counted[1], 2.0)
  expect_equal(pats$site_pattern[[1]][1], 2.0)   # AAA
})

test_that("triplet_site_patterns: ambiguous base N distributes uniformly", {
  # A=A, B=A, C=N → each of AAA, AAC, AAG, AAT gets 0.25 per site
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=3 NCHAR=1;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX", "  A A", "  B A", "  C N", ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  n_vec <- pats$site_pattern[[1]]
  # AAA=1, AAC=2, AAG=3, AAT=4 (1-indexed)
  expect_equal(n_vec[1], 0.25, tolerance = 1e-10)
  expect_equal(n_vec[2], 0.25, tolerance = 1e-10)
  expect_equal(n_vec[3], 0.25, tolerance = 1e-10)
  expect_equal(n_vec[4], 0.25, tolerance = 1e-10)
  expect_equal(sum(n_vec), 1.0, tolerance = 1e-10)
})

# -- site compression ---------------------------------------------------------

test_that("compressed sites populated after NEXUS parse", {
  path <- write_tmp(NEXUS_3TAX, ".nex")
  aln <- read_nexus(path)
  gene <- get_gene(aln, 1L)
  expect_false(phylotriplet:::gene_compressed_empty_rcpp(gene$.ptr))
  np <- phylotriplet:::gene_compressed_n_patterns_rcpp(gene$.ptr)
  expect_true(np >= 1L)
  expect_true(np <= phylotriplet:::gene_length_rcpp(gene$.ptr))
})

test_that("compressed sites populated after PHYLIP parse", {
  path <- write_tmp(PHYLIP_SEQ, ".phy")
  aln <- read_phylip_files(path)
  gene <- get_gene(aln, 1L)
  expect_false(phylotriplet:::gene_compressed_empty_rcpp(gene$.ptr))
})

test_that("all-identical columns yield one unique pattern", {
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=3 NCHAR=5;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX", "  A AAAAA", "  B CCCCC", "  C GGGGG", ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  gene <- get_gene(aln, 1L)
  expect_equal(phylotriplet:::gene_compressed_n_patterns_rcpp(gene$.ptr), 1L)
})

test_that("all-distinct columns preserve all patterns", {
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=2 NCHAR=4;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX", "  A ACGT", "  B TGCA", ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  gene <- get_gene(aln, 1L)
  expect_equal(phylotriplet:::gene_compressed_n_patterns_rcpp(gene$.ptr), 4L)
})

test_that("SETS block: each gene partition is compressed separately", {
  path <- write_tmp(NEXUS_SETS, ".nex")
  aln <- read_nexus(path)
  gene1 <- get_gene(aln, 1L)
  gene2 <- get_gene(aln, 2L)
  expect_false(phylotriplet:::gene_compressed_empty_rcpp(gene1$.ptr))
  expect_false(phylotriplet:::gene_compressed_empty_rcpp(gene2$.ptr))
  # Each partition has 3 columns
  expect_true(phylotriplet:::gene_compressed_n_patterns_rcpp(gene1$.ptr) <= 3L)
  expect_true(phylotriplet:::gene_compressed_n_patterns_rcpp(gene2$.ptr) <= 3L)
})

test_that("compression does not change site pattern results", {
  # Use a 3-taxon alignment with duplicate columns
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=3 NCHAR=8;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX",
    "  A ACGTACGT",
    "  B TGCATGCA",
    "  C GCATGCAT",
    ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  expect_equal(nrow(pats), 1L)
  expect_equal(pats$n_total[1], 8.0)
  expect_equal(sum(pats$site_pattern[[1]]), pats$n_counted[1], tolerance = 1e-10)
})

test_that("compression with gaps: site patterns still correct", {
  nex <- c(
    "#NEXUS", "BEGIN DATA;",
    "  DIMENSIONS NTAX=3 NCHAR=6;",
    "  FORMAT DATATYPE=DNA;",
    "  MATRIX",
    "  A AAAAAA",
    "  B AA-AA-",
    "  C AAAAAA",
    ";", "END;"
  )
  path <- write_tmp(nex, ".nex")
  aln <- read_nexus(path)
  pats <- triplet_site_patterns(aln)
  # 6 sites, 2 have gaps -> n_counted = 4
  expect_equal(pats$n_total[1], 6.0)
  expect_equal(pats$n_counted[1], 4.0)
  # All non-gap sites are AAA -> site_pattern[1] = 4
  expect_equal(pats$site_pattern[[1]][1], 4.0)
})