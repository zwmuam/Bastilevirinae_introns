import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from alignment_statistics import (
    minor_group, major_group, ungapped_length, trim_end_gaps,
    pairwise_identity, count_intron_groups_in_file,
    intron_aligned_lengths, group_similarity
)


def test_minor_and_major_group():
    assert minor_group("IA1__seq123") == "IA1"
    assert major_group("IA1__seq123") == "IA_total"

    assert minor_group("IA__seq123") == "IA?"
    assert minor_group("") == "?"


def test_ungapped_length():
    assert ungapped_length("A-C-G-T") == 4
    assert ungapped_length("----") == 0


def test_trim_end_gaps():
    s1, s2 = trim_end_gaps("--ACGT--", "--AC-T--")
    assert s1 == "ACGT"
    assert s2 == "AC-T"

    s1_empty, s2_empty = trim_end_gaps("----", "----")
    assert s1_empty == ""
    assert s2_empty == ""


def test_pairwise_identity():
    id1 = pairwise_identity("ACGT", "ACGT")
    assert id1 == 100.0

    id2 = pairwise_identity("ACGT", "ACTT")
    assert id2 == 75.0


def test_count_intron_groups_in_file():
    r1 = SeqRecord(Seq("ACGT"), id="IA1__1")
    r2 = SeqRecord(Seq("ACGT"), id="IA2__2")
    r3 = SeqRecord(Seq("ACGT"), id="IB1__3")

    df = count_intron_groups_in_file([r1, r2, r3])
    assert "IA1" in df.index
    assert "IA_total" in df.index
    assert df.loc["IA_total", "count"] == 2
