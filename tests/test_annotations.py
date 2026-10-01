import pytest
from pathlib import Path
from Bio import Seq
from annotations import (
    Annotation, InfernalAlignment, HmmerAlignment, DnaAlignment,
    Exon, Intron, Gene, AnnotationTrack, AnnotationBase
)


def test_annotation_basic():
    a = Annotation(seq_id="seq1", model_id="M1", model_name="Model1",
                   start=10, end=50, score=100.0, evalue=1e-5,
                   method="test", strand="+")
    assert repr(a) == "seq1__Annotation__M1__10__50__+"
    assert a.model_definition() == "M1. Model1"
    assert len(a) == 41

    ctx = a.context(5)
    assert ctx.start == 5
    assert ctx.end == 55

    rc = a.reverse_complement(seq_len=100)
    assert rc.start == 51
    assert rc.end == 91
    assert rc.strand == "-"


def test_annotation_overlaps():
    a1 = Annotation("seq1", "M1", "M1", 10, 50, 10.0, 0.01, "test", "+")
    a2 = Annotation("seq1", "M2", "M2", 30, 70, 10.0, 0.01, "test", "+")
    a3 = Annotation("seq1", "M3", "M3", 80, 100, 10.0, 0.01, "test", "+")

    assert a1.overlaps(a2) > 0
    assert a1.overlaps(a3) == 0


def test_gene_from_exons():
    exon1 = DnaAlignment("seq1", "M1", "GeneModel", 10, 40, 50.0, 1e-3, "test", "+")
    exon2 = DnaAlignment("seq1", "M1", "GeneModel", 100, 150, 60.0, 1e-4, "test", "+")

    gene = Gene.from_exons([exon1, exon2], min_intron=20)
    assert gene.start == 10
    assert gene.end == 150
    assert len(gene.children(Intron)) == 1

    intron = gene.children(Intron)[0]
    assert intron.start == 41
    assert intron.end == 99


def test_annotation_track_cull():
    a1 = Annotation("seq1", "M1", "M1", 10, 50, 100.0, 0.01, "test", "+")
    a2 = Annotation("seq1", "M2", "M2", 20, 45, 50.0, 0.01, "test", "+")  # lower score overlapping
    a3 = Annotation("seq1", "M3", "M3", 80, 120, 90.0, 0.01, "test", "+")  # non-overlapping

    track = AnnotationTrack("seq1", [a1, a2, a3])
    culled = track.cull(overlap_threshold=0.3)
    assert len(culled) == 2
    assert culled[0].model_id == "M1"
    assert culled[1].model_id == "M3"


def test_annotation_base_gff_and_sort(tmp_path):
    a1 = Annotation("seq1", "M1", "M1", 50, 100, 10.0, 0.01, "test", "+")
    a2 = Annotation("seq1", "M2", "M2", 10, 40, 20.0, 0.01, "test", "+")

    base = AnnotationBase()
    base.annotate(a1)
    base.annotate(a2)

    base.sort_annotations(by='start')
    assert base["seq1"][0].start == 10
    assert base["seq1"][1].start == 50

    gff_file = tmp_path / "out.gff"
    base.save_gff(gff_file)
    assert gff_file.exists()
    content = gff_file.read_text()
    assert "seq1" in content
    assert "M1" in content
