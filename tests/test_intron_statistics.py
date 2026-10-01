import pytest
import pandas as pd
from annotations import AnnotationBase, Gene, Intron, DnaAlignment, InfernalAlignment
from intron_statistics import (
    plot_intron_distribution, plot_intron_lengths, introns_in_genomes,
    intron_architectures, nuclease_table, taxonomic_annotations
)


def create_sample_annotations():
    base = AnnotationBase()
    gene = Gene(seq_id="phage1|contig1", strand="+", model_id="P1", model_name="Phrog1",
                start=100, end=1000, score=100.0, evalue=1e-5, method="test")

    intron = Intron(seq_id="phage1|contig1", strand="+", model_id="I1", model_name="Intron1",
                    start=200, end=500, score=50.0, evalue=1e-3, method="test")

    cds = DnaAlignment(seq_id="phage1|contig1", strand="+", model_id="phrog_123", model_name="Endonuclease",
                       start=250, end=450, score=40.0, evalue=1e-2, method="test")
    intron.host(cds)
    gene.host(intron)
    base.annotate(gene)
    return base


def test_plot_intron_distribution():
    base = create_sample_annotations()
    fig = plot_intron_distribution(base, id2name_dict={})
    assert fig is not None


def test_plot_intron_lengths():
    base = create_sample_annotations()
    fig = plot_intron_lengths(base, id2name_dict={})
    assert fig is not None


def test_introns_in_genomes():
    base = create_sample_annotations()
    df = introns_in_genomes(base, id2name_dict={})
    assert "phage1" in df.index
    assert df.loc["phage1", "n_introns"] == 1


def test_taxonomic_annotations():
    base = create_sample_annotations()
    genome_tab = introns_in_genomes(base, id2name_dict={})

    taxon_df = pd.DataFrame({
        "organism": ["Phage One"],
        "species": ["Phage species 1"],
        "genus": ["Phagegenus"]
    }, index=["phage1"])

    annotated_tab, summary = taxonomic_annotations(genome_tab, taxon_df)
    assert annotated_tab.loc["phage1", "genus"] == "Phagegenus"
    assert "Phagegenus" in summary.index
