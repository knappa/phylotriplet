"""Tests for the Python sequence-alignment bindings."""

import os
import tempfile
import pytest
import phylotriplet as pt

# -- fixtures / helpers --------------------------------------------------------

NEXUS_3TAX = """\
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=4;
  FORMAT DATATYPE=DNA;
  MATRIX
  Alpha ACGT
  Beta  TTAG
  Gamma CGCA
;
END;
"""

NEXUS_WITH_TITLE = """\
#NEXUS
BEGIN DATA;
  TITLE COI;
  DIMENSIONS NTAX=2 NCHAR=3;
  FORMAT DATATYPE=DNA;
  MATRIX
  T1 ACG
  T2 TTA
;
END;
"""

NEXUS_MULTIGENE = """\
#NEXUS
BEGIN DATA;
  TITLE gene_a;
  DIMENSIONS NTAX=3 NCHAR=3;
  FORMAT DATATYPE=DNA;
  MATRIX
  A AAA
  B AAA
  C AAA
;
END;
BEGIN DATA;
  TITLE gene_b;
  DIMENSIONS NTAX=3 NCHAR=3;
  FORMAT DATATYPE=DNA;
  MATRIX
  A CCC
  B CCC
  C CCC
;
END;
"""

PHYLIP_SEQ = """\
3 4
Alpha ACGT
Beta  TTAG
Gamma CGCA
"""

PHYLIP_SEQ_2 = """\
2 3
X ACG
Y TTA
"""


@pytest.fixture
def nexus_3tax(tmp_path):
    p = tmp_path / "three.nex"
    p.write_text(NEXUS_3TAX)
    return str(p)


@pytest.fixture
def nexus_multigene(tmp_path):
    p = tmp_path / "multi.nex"
    p.write_text(NEXUS_MULTIGENE)
    return str(p)


@pytest.fixture
def phylip_seq(tmp_path):
    p = tmp_path / "seq.phy"
    p.write_text(PHYLIP_SEQ)
    return str(p)


# -- read_nexus_file -----------------------------------------------------------

class TestReadNexusFile:

    def test_returns_multigene_alignment(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        assert isinstance(aln, pt.MultiGeneAlignment)

    def test_num_genes(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        assert aln.num_genes == 1

    def test_taxon_names_sorted(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        names = aln.all_taxon_names()
        assert names == sorted(names)
        assert set(names) == {"Alpha", "Beta", "Gamma"}

    def test_multiple_data_blocks(self, nexus_multigene):
        aln = pt.read_nexus_file(nexus_multigene)
        assert aln.num_genes == 2

    def test_title_parsed(self, tmp_path):
        p = tmp_path / "titled.nex"
        p.write_text(NEXUS_WITH_TITLE)
        aln = pt.read_nexus_file(str(p))
        assert aln.num_genes == 1
        assert set(aln.all_taxon_names()) == {"T1", "T2"}

    def test_missing_nexus_header_raises(self, tmp_path):
        p = tmp_path / "bad.nex"
        p.write_text("BEGIN DATA;\nEND;\n")
        with pytest.raises(RuntimeError):
            pt.read_nexus_file(str(p))

    def test_unsupported_datatype_raises(self, tmp_path):
        p = tmp_path / "prot.nex"
        p.write_text(
            "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=3;\n"
            "  FORMAT DATATYPE=PROTEIN;\n  MATRIX\n  T1 ACG\n  T2 TTA\n;\nEND;\n"
        )
        with pytest.raises(RuntimeError):
            pt.read_nexus_file(str(p))

    def test_file_not_found_raises(self):
        with pytest.raises(RuntimeError):
            pt.read_nexus_file("/no/such/file.nex")

    def test_repr(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        r = repr(aln)
        assert "MultiGeneAlignment" in r
        assert "num_genes=1" in r


# -- read_phylip_file / read_phylip_files --------------------------------------

class TestReadPhylipFile:

    def test_single_file_returns_one_gene(self, phylip_seq):
        aln = pt.read_phylip_file(phylip_seq)
        assert isinstance(aln, pt.MultiGeneAlignment)
        assert aln.num_genes == 1

    def test_taxon_names(self, phylip_seq):
        aln = pt.read_phylip_file(phylip_seq)
        assert set(aln.all_taxon_names()) == {"Alpha", "Beta", "Gamma"}

    def test_file_not_found_raises(self):
        with pytest.raises(RuntimeError):
            pt.read_phylip_file("/no/such/file.phy")


class TestReadPhylipFiles:

    def test_two_files_two_genes(self, tmp_path, phylip_seq):
        p2 = tmp_path / "g2.phy"
        p2.write_text(PHYLIP_SEQ_2)
        aln = pt.read_phylip_files([phylip_seq, str(p2)])
        assert aln.num_genes == 2

    def test_gene_names_passed_through(self, tmp_path, phylip_seq):
        p2 = tmp_path / "g2.phy"
        p2.write_text(PHYLIP_SEQ_2)
        aln = pt.read_phylip_files([phylip_seq, str(p2)],
                                   gene_names=["geneA", "geneB"])
        assert aln.num_genes == 2

    def test_empty_list_gives_zero_genes(self):
        aln = pt.read_phylip_files([])
        assert aln.num_genes == 0


# -- MultiGeneAlignment.triplet_site_patterns ---------------------------------

class TestComputeSitePatterns:

    def test_returns_list_of_triplet_site_patterns(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        pats = aln.triplet_site_patterns()
        assert isinstance(pats, list)
        assert all(isinstance(p, pt.TripletSitePatterns) for p in pats)

    def test_three_taxa_one_triplet(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        assert len(aln.triplet_site_patterns()) == 1

    def test_taxa_in_lexicographic_order(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        for p in aln.triplet_site_patterns():
            assert p.taxa[0] <= p.taxa[1] <= p.taxa[2]

    def test_n_is_length_64(self, nexus_3tax):
        aln = pt.read_nexus_file(nexus_3tax)
        for p in aln.triplet_site_patterns():
            assert len(p.site_pattern) == 64

    def test_n_total_and_n_counted(self, nexus_3tax):
        p = pt.read_nexus_file(nexus_3tax).triplet_site_patterns()[0]
        assert p.site_pattern_total == pytest.approx(4.0)
        assert p.site_pattern_counted == pytest.approx(4.0)

    def test_sum_n_equals_n_counted(self, nexus_3tax):
        for p in pt.read_nexus_file(nexus_3tax).triplet_site_patterns():
            assert sum(p.site_pattern) == pytest.approx(p.site_pattern_counted)

    def test_all_A_alignment_fills_aaa_bucket(self, tmp_path):
        p = tmp_path / "aaa.nex"
        p.write_text(
            "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=3 NCHAR=3;\n"
            "  FORMAT DATATYPE=DNA;\n  MATRIX\n  A AAA\n  B AAA\n  C AAA\n;\nEND;\n"
        )
        pats = pt.read_nexus_file(str(p)).triplet_site_patterns()
        assert len(pats) == 1
        # AAA: i=0, j=0, k=0 → index 0
        assert pats[0].site_pattern[0] == pytest.approx(3.0)
        assert sum(pats[0].site_pattern[1:]) == pytest.approx(0.0)

    def test_gap_site_reduces_n_counted(self, tmp_path):
        p = tmp_path / "gap.site_patternex"
        p.write_text(
            "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=3 NCHAR=4;\n"
            "  FORMAT DATATYPE=DNA;\n  MATRIX\n"
            "  A AAAA\n  B A-AA\n  C AAAA\n;\nEND;\n"
        )
        pats = pt.read_nexus_file(str(p)).triplet_site_patterns()
        assert pats[0].site_pattern_total == pytest.approx(4.0)
        assert pats[0].site_pattern_counted == pytest.approx(3.0)

    def test_ambiguous_N_distributes_uniformly(self, tmp_path):
        # A=A, B=A, C=N → AAA, AAC, AAG, AAT each get 0.25
        p = tmp_path / "N.nex"
        p.write_text(
            "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=3 NCHAR=1;\n"
            "  FORMAT DATATYPE=DNA;\n  MATRIX\n"
            "  A A\n  B A\n  C N\n;\nEND;\n"
        )
        n = pt.read_nexus_file(str(p)).triplet_site_patterns()[0].site_pattern
        # AAA=0, AAC=1, AAG=2, AAT=3
        assert n[0] == pytest.approx(0.25)
        assert n[1] == pytest.approx(0.25)
        assert n[2] == pytest.approx(0.25)
        assert n[3] == pytest.approx(0.25)
        assert sum(n) == pytest.approx(1.0)

    def test_multigene_accumulates_patterns(self, nexus_multigene):
        # gene_a: all-A (3 sites) → AAA bucket (index 0)
        # gene_b: all-C (3 sites) → CCC bucket (index 1*16+1*4+1 = 21)
        pats = pt.read_nexus_file(nexus_multigene).triplet_site_patterns()
        assert len(pats) == 1
        p = pats[0]
        assert p.site_pattern[0]  == pytest.approx(3.0)
        assert p.site_pattern[21] == pytest.approx(3.0)
        assert p.site_pattern_total   == pytest.approx(6.0)
        assert p.site_pattern_counted == pytest.approx(6.0)

    def test_four_taxa_gives_four_triplets(self, tmp_path):
        p = tmp_path / "four.nex"
        p.write_text(
            "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=4 NCHAR=2;\n"
            "  FORMAT DATATYPE=DNA;\n  MATRIX\n"
            "  A AC\n  B TG\n  C CA\n  D GT\n;\nEND;\n"
        )
        pats = pt.read_nexus_file(str(p)).triplet_site_patterns()
        assert len(pats) == 4  # C(4,3)

    def test_repr_of_triplet_site_patterns(self, nexus_3tax):
        p = pt.read_nexus_file(nexus_3tax).triplet_site_patterns()[0]
        r = repr(p)
        assert "TripletSitePatterns" in r
        assert "n_counted" in r

    def test_end_to_end_phylip(self, phylip_seq):
        pats = pt.read_phylip_file(phylip_seq).triplet_site_patterns()
        assert len(pats) == 1
        p = pats[0]
        assert p.site_pattern_total == pytest.approx(4.0)
        assert sum(p.site_pattern) == pytest.approx(p.site_pattern_counted)