import pytest
from pathlib import Path
from click.testing import CliRunner
from unittest.mock import patch

from sto2fasta import convert_and_save
from prune_alignment import filter_short_alignments
from alignment_statistics import alignment_stats
from find_introns import find_introns
from reannotate_introns import annotate_introns
from confirm_introns import rnaseq_analysis


def test_sto2fasta_cli(tmp_path):
    runner = CliRunner()
    sto_file = tmp_path / "align.sto"
    sto_file.write_text("""# STOCKHOLM 1.0
seq1 ACGTACGT
seq2 ACGT--GT
//
""")
    result = runner.invoke(convert_and_save, ["-s", str(sto_file)])
    assert result.exit_code == 0
    fasta_file = tmp_path / "align.fasta"
    assert fasta_file.exists()
    assert ">seq1" in fasta_file.read_text()


def test_prune_alignment_cli(tmp_path):
    runner = CliRunner()
    fasta_file = tmp_path / "align.fasta"
    fasta_file.write_text(">seq1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>seq2\nA--G\n")
    result = runner.invoke(filter_short_alignments, ["-f", str(fasta_file), "-l", "10"])
    assert result.exit_code == 0
    pruned_file = tmp_path / "align.pruned.fasta"
    assert pruned_file.exists()
    content = pruned_file.read_text()
    assert "seq1" in content
    assert "seq2" not in content


def test_alignment_statistics_cli(tmp_path):
    runner = CliRunner()
    fasta_file = tmp_path / "IA1__seq1.fasta"
    fasta_file.write_text(">IA1__seq1\nACGTACGT\n>IA2__seq2\nACGTACGT\n")
    result = runner.invoke(alignment_stats, ["-f", str(fasta_file), "-u"])
    assert result.exit_code == 0
    out_xlsx = tmp_path / "IA1__seq1_alignment_stats.xlsx"
    assert out_xlsx.exists()


@patch("find_introns.run_external")
def test_find_introns_cli(mock_run_ext, tmp_path):
    runner = CliRunner()

    fasta_file = tmp_path / "genome.fasta"
    fasta_file.write_text(">seq1\n" + "ACGT" * 2500 + "\n")

    out_dir = tmp_path / "out_find"

    # Mock infernal and hmmer tblout output
    def mock_run_side_effect(cmd, *args, **kwargs):
        if "--tblout" in cmd:
            tblout_idx = cmd.index("--tblout") + 1
            tblout = Path(cmd[tblout_idx])
            tblout.write_text("seq1 - model1 CM1 cm 1 100 3000 3500 + 1 no 1 0.0 100.0 1e-10 ! - description\n")
        elif "--domtblout" in cmd:
            domtblout_idx = cmd.index("--domtblout") + 1
            domtblout = Path(cmd[domtblout_idx])
            head = "model1 - 100 seq1__InfernalAlignment__CM1__500__6000__+___1 - 300 1e-10 100.0 0.0 1 1 1e-10 1e-10 100.0 0.0 1 100 10 40 10 40 0.95 - description\n"
            tail = "model1 - 100 seq1__InfernalAlignment__CM1__500__6000__+___1 - 300 1e-10 100.0 0.0 1 1 1e-10 1e-10 100.0 0.0 1 100 100 140 100 140 0.95 - description\n"
            domtblout.write_text(head + tail)

    mock_run_ext.side_effect = mock_run_side_effect

    cm_file = tmp_path / "dummy.cm"
    cm_file.write_text("dummy cm")
    hmm_file = tmp_path / "dummy.hmm"
    hmm_file.write_text("dummy hmm")

    phrog_table = tmp_path / "phrog.tsv"
    phrog_table.write_text("phrog\tannot\n1\tEndonuclease\n")

    result = runner.invoke(find_introns, [
        "-f", str(fasta_file),
        "-c", str(cm_file),
        "-h", str(hmm_file),
        "-o", str(out_dir),
        "-r", str(phrog_table)
    ])

    assert result.exit_code == 0
    assert out_dir.exists()


@patch("reannotate_introns.run_external")
def test_reannotate_introns_cli(mock_run_ext, tmp_path):
    runner = CliRunner()

    fasta_file = tmp_path / "introns.fasta"
    fasta_file.write_text(">intron1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")

    out_dir = tmp_path / "out_reannot"

    def mock_run_side_effect(cmd, *args, **kwargs):
        if "--tblout" in cmd:
            tblout_idx = cmd.index("--tblout") + 1
            tblout = Path(cmd[tblout_idx])
            tblout.write_text("intron1 - model1 CM1 cm 1 50 5 25 + 1 no 1 0.0 100.0 1e-10 ! - description\n")
        elif "--domtblout" in cmd:
            domtblout_idx = cmd.index("--domtblout") + 1
            domtblout = Path(cmd[domtblout_idx])
            domtblout.write_text("model1 - 100 intron1___1 - 300 1e-10 100.0 0.0 1 1 1e-10 1e-10 100.0 0.0 1 100 2 10 2 10 0.95 - description\n")

    mock_run_ext.side_effect = mock_run_side_effect

    cm_file = tmp_path / "dummy.cm"
    cm_file.write_text("dummy cm")
    hmm_file = tmp_path / "dummy.hmm"
    hmm_file.write_text("dummy hmm")

    result = runner.invoke(annotate_introns, [
        "-f", str(fasta_file),
        "-c", str(cm_file),
        "-h", str(hmm_file),
        "-o", str(out_dir)
    ])

    assert result.exit_code == 0
    assert out_dir.exists()


@patch("confirm_introns.run_external")
def test_confirm_introns_cli(mock_run_ext, tmp_path):
    runner = CliRunner()

    adapters = tmp_path / "adapters.fasta"
    adapters.write_text(">a1\nACGT\n")
    barcodes = tmp_path / "barcodes.fasta"
    barcodes.write_text(">b1\nACGT\n")

    clean_dir = tmp_path / "clean_reads"
    clean_dir.mkdir()
    read_fq = clean_dir / "ref1__bc1.fastq"
    read_fq.write_text("@read1\nACGT\n+\nIIII\n")

    ref_dir = tmp_path / "ref_dir"
    ref_dir.mkdir()
    ref_fa = ref_dir / "ref1.fasta"
    ref_fa.write_text(">ref1\nACGTACGTACGT\n")

    out_dir = tmp_path / "out_confirm"

    def mock_run_side_effect(cmd, *args, **kwargs):
        if "samtools" in cmd and "-c" in cmd:
            return b"150\n"
        elif "spliced_bam2gff" in cmd:
            return b"ref1\tspliced\tmRNA\t1\t100\t.\t+\t.\tID=m1\nref1\tspliced\texon\t1\t30\t.\t+\t.\tID=m1.e1\nref1\tspliced\texon\t70\t100\t.\t+\t.\tID=m1.e2\n"
        return b""

    mock_run_ext.side_effect = mock_run_side_effect

    result = runner.invoke(rnaseq_analysis, [
        "-ca", str(adapters),
        "-cb", str(barcodes),
        "-cr", str(clean_dir),
        "-rd", str(ref_dir),
        "-o", str(out_dir),
        "-s", "__"
    ])

    assert result.exit_code == 0
    assert out_dir.exists()
