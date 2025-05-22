import pytest
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from protein_compare_toolkit.core.stats_utils import aln_dist, jensen_shannon_distance
from protein_compare_toolkit.core.stats_utils import info_content, select_diff_index, AA_COUNT

def test_aln_dist_simple_case():
    # Create a toy alignment: 3 sequences of length 4
    aln = MultipleSeqAlignment([
        SeqRecord(Seq("ACDE"), id="seq1"),
        SeqRecord(Seq("ACDE"), id="seq2"),
        SeqRecord(Seq("ACDE"), id="seq3"),
    ])

    # Run the function
    dist = aln_dist(aln)

    # Shape should be (4 positions, 20 amino acids)
    assert dist.shape == (4, AA_COUNT)

    # Each column in the 4 rows should sum to 1 (perfect agreement on A, C, D, E)
    for i in range(4):
        assert np.isclose(np.sum(dist[i]), 1.0)


def test_aln_dist_mixed_residues():
    aln = MultipleSeqAlignment([
        SeqRecord(Seq("ACDE"), id="seq1"),
        SeqRecord(Seq("ACDF"), id="seq2"),
        SeqRecord(Seq("ACDG"), id="seq3"),
    ])

    dist = aln_dist(aln)

    # First 3 positions should still be uniform
    for i in range(3):
        aa = "ACD"[i]
        idx = "ACDEFGHIKLMNPQRSTVWY".index(aa)
        assert dist[i][idx] == 1.0

    # Fourth position: 1 E, 1 F, 1 G -> each should have 1/3
    for aa in "EFG":
        idx = "ACDEFGHIKLMNPQRSTVWY".index(aa)
        assert np.isclose(dist[3][idx], 1/3)

    # All others should be near zero
    remaining = set("ACDEFGHIKLMNPQRSTVWY") - set("EFG")
    for aa in remaining:
        idx = "ACDEFGHIKLMNPQRSTVWY".index(aa)
        assert np.isclose(dist[3][idx], 0.0)


def test_aln_dist_with_gaps():
    aln = MultipleSeqAlignment([
        SeqRecord(Seq("A--E"), id="seq1"),
        SeqRecord(Seq("A--E"), id="seq2"),
        SeqRecord(Seq("A--E"), id="seq3"),
    ])

    dist = aln_dist(aln)

    # Position 0 and 3 should be valid (A and E)
    idx_a = "ACDEFGHIKLMNPQRSTVWY".index("A")
    idx_e = "ACDEFGHIKLMNPQRSTVWY".index("E")
    assert dist[0][idx_a] == 1.0
    assert dist[3][idx_e] == 1.0

    # Positions 1 and 2 are all gaps, should be all uniform distribution.
    assert np.allclose(dist[1], 1 / AA_COUNT)
    assert np.allclose(dist[2], 1 / AA_COUNT)


def test_aln_dist_weighted_counts():
    aln = MultipleSeqAlignment([
        SeqRecord(Seq("A"), id="seq1"),
        SeqRecord(Seq("A"), id="seq2"),
        SeqRecord(Seq("C"), id="seq3"),
        SeqRecord(Seq("C"), id="seq4"),
        SeqRecord(Seq("C"), id="seq5"),
    ])

    dist = aln_dist(aln)

    # Position 0: 2 A's, 3 C's => A: 2/5, C: 3/5
    idx_a = "ACDEFGHIKLMNPQRSTVWY".index("A")
    idx_c = "ACDEFGHIKLMNPQRSTVWY".index("C")
    assert np.isclose(dist[0][idx_a], 0.4)
    assert np.isclose(dist[0][idx_c], 0.6)
    # All others should be zero
    others = set(range(20)) - {idx_a, idx_c}
    assert np.allclose(dist[0][list(others)], 0.0)


def make_alignment(seq_list):
    return MultipleSeqAlignment(
        [SeqRecord(Seq(seq), id=f"seq{i}") for i, seq in enumerate(seq_list)])


def test_jsd_identical_alignments():
    aln1 = make_alignment(["ACDE", "ACDE", "ACDE"])
    aln2 = make_alignment(["ACDE", "ACDE", "ACDE"])

    jsd = jensen_shannon_distance(aln1, aln2)

    # All distances should be 0
    assert np.allclose(jsd, 0.0)


def test_jsd_completely_different():
    aln1 = make_alignment(["AAAA", "AAAA", "AAAA"])
    aln2 = make_alignment(["CCCC", "CCCC", "CCCC"])

    jsd = jensen_shannon_distance(aln1, aln2)

    # Jensen-Shannon distance between delta distributions is 1
    assert np.allclose(jsd, 1.0)


def test_jsd_half_half():
    aln1 = make_alignment(["AAAA", "AAAA", "CCCC", "CCCC"])
    aln2 = make_alignment(["AAAA", "CCCC", "CCCC", "CCCC"])

    jsd = jensen_shannon_distance(aln1, aln2)

    # Result should be between 0 and 1
    assert np.all(jsd >= 0) and np.all(jsd <= 1)
    # The function should be symmetric
    jsd_reversed = jensen_shannon_distance(aln2, aln1)
    assert np.allclose(jsd, jsd_reversed)


def test_jensen_shannon_mismatched_lengths_raises():
    aln1 = make_alignment(["ACD", "ACD"])
    aln2 = make_alignment(["ACDE", "ACDE"])

    with pytest.raises(ValueError, match="Alignment lengths differ"):
        jensen_shannon_distance(aln1, aln2)


def test_info_content_values():
    aln = make_alignment([
        "AAA",
        "AAA",
        "AAA",
        "AAA",
        "AAA",
        "AAA",
        "AAA",
        "AAA",
        "AAA",
        "AAA"
    ])
    ic = info_content(aln)

    # Should be close to log2(20) minus small correction
    expected = np.log2(20) - (1 / np.log(2)) * (19 / (2 * 10))  # r0 - e
    assert np.allclose(ic, np.full(3, expected), atol=0.01)
    assert np.all(ic >= 0)


def test_info_content_uniform():
    aln = make_alignment([
        "ACD",
        "DFG",
        "VKR",
        "LMN",
        "PQY"
    ])
    ic = info_content(aln)

    # These positions should be roughly uniformly distributed => IC near 0
    assert np.all(ic < 0.5)


def test_info_content_empty_alignment_raises():
    aln = make_alignment([])
    with pytest.raises(ValueError, match="Alignment must contain at least one"):
        info_content(aln)


def test_select_diff_index_shapes():
    aln1 = make_alignment([
        "ACD",
        "ACD",
        "ACD"
    ])
    aln2 = make_alignment([
        "AED",
        "AED",
        "AED"
    ])
    sdi1, sdi2 = select_diff_index(aln1, aln2)

    assert sdi1.shape == (3,)
    assert sdi2.shape == (3,)
    assert np.all(sdi1 >= 0)
    assert np.all(sdi2 >= 0)
