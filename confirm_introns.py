#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses nanopore sequencing reads to confirm splicing of the predicted introns and
maps exact splicing sites for higly represented variants.
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

from pathlib import Path
from random import randint

import click
from porechop import porechop_custom

from annotations import Exon, Intron, Gene, AnnotationBase
from tweaks import run_external


@click.command(no_args_is_help=True)
@click.option("-ca", "--custom_adapters",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='fasta file with custom adapters (Forward Forward orientation >__>)')
@click.option("-cb", "--custom_barcodes",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='fasta file with custom barcodes (Forward Forward orientation >__>)')
@click.option("-fq", "--read_fastq",
              required=False,
              type=click.Path(exists=True, path_type=Path),
              help='path to fastq file')
@click.option("-cr", "--clean_reads",
              required=False,
              type=click.Path(exists=True, path_type=Path),
              help='directory with demultiplexed clean reads (skips demultiplexing and trimming steps)')
@click.option("-rd", "--reference_dir",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='directory with reference fasta files')
@click.option("-o", "--output",
              required=True,
              type=click.Path(exists=False, path_type=Path),
              help='output directory')
@click.option("-s", "--separator",
              required=True,
              type=str,
              help='separator for the primer names')
@click.option("-mr", "--max_reads",
              required=False,
              type=int,
              default=500,
              help='number of reads to subsample')
@click.option("-ml", "--min_length",
              required=False,
              type=int,
              default=500,
              help='minimum read length to keep during quality filtering')
@click.option("-mq", "--min_mean_q",
              required=False,
              type=float,
              default=90,
              help='minimum mean quality to keep during read filtering')
def rnaseq_analysis(custom_adapters: Path,
                    custom_barcodes: Path,
                    read_fastq: Path,
                    clean_reads: Path,
                    reference_dir: Path,
                    output: Path,
                    separator: str,
                    max_reads: int,
                    min_length: int,
                    min_mean_q: float):
    """
    De-multiplex reads and map them to respective references
    :param custom_adapters:
    :param custom_barcodes: reverse primers fasta
    :param read_fastq: fastq file with reads
    :param clean_reads: directory with demultiplexed clean reads
    :param reference_dir: directory with reference fasta files
    :param output: output directory
    :param separator: separator for the primer names
    :param max_reads: number of reads to subsample
    :param min_length: minimum read length to keep during quality filtering
    :param min_mean_q: minimum mean quality to keep during read filtering
    """
    assert read_fastq or clean_reads, 'Provide either fastq to trim or clean_reads'
    output.mkdir()
    if not clean_reads:
        demultiplexing_dir = output.joinpath('demultiplexing')
        demultiplexing_dir.mkdir()
        filtered_read_dir = output.joinpath(f'filtered.{min_length}bp_q{min_mean_q}')
        filtered_read_dir.mkdir()

        clean_reads = trim_and_demultiplex(read_fastq=read_fastq,
                                           custom_adapters=custom_adapters,
                                           custom_barcodes=custom_barcodes,
                                           demultiplexing_dir=demultiplexing_dir,
                                           filtered_dir=filtered_read_dir,
                                           min_length=min_length,
                                           min_mean_q=min_mean_q)
    else:
        clean_reads = [Path(f) for f in clean_reads.iterdir() if f.suffix in ('.fastq', '.fq')]

    minimap_sorted_dir = output.joinpath('minimap2_sorted')
    minimap_sorted_dir.mkdir()
    minumap_filtered_dir = output.joinpath('minimap2_filtered')
    minumap_filtered_dir.mkdir()
    raw_gff_dir = output.joinpath('raw_gff')
    raw_gff_dir.mkdir()

    filtered_introns = AnnotationBase()

    for read_fastq in clean_reads:
        reference_name = read_fastq.stem.split(separator)[0]
        reference = reference_dir.joinpath(f'{reference_name}.fasta')
        if reference.exists():
            filtered_file = map_to_reference(reference=reference,
                                             read_fastq=read_fastq,
                                             minimap_sorted_dir=minimap_sorted_dir,
                                             minumap_filtered_dir=minumap_filtered_dir)
            n_mapped_reads = count_mapped_reads(filtered_file)
            filtered_introns = splice_sites(spliced_bam=filtered_file,
                                            raw_gff_dir=raw_gff_dir,
                                            mapped_reads=n_mapped_reads,
                                            output=filtered_introns)
            if n_mapped_reads > max_reads:
                subsampled_bam = filtered_file.parent.joinpath(f'{filtered_file.stem}.{max_reads}_sample_reads.bam')
                subsample_bam(bam=filtered_file,
                              subsampled_bam=subsampled_bam,
                              max_reads=max_reads,
                              mapped_reads=n_mapped_reads)
                n_subsampled_reads = count_mapped_reads(subsampled_bam)
                print(f'sampled {n_subsampled_reads} out of {n_mapped_reads} reads')
        else:
            print(f"Reference {reference_name} not found")
    filtered_introns.save_gff(output.joinpath('Introns.gff'))

def trim_and_demultiplex(read_fastq: Path,
                         custom_adapters: Path,
                         custom_barcodes: Path,
                         demultiplexing_dir: Path,
                         filtered_dir: Path,
                         min_length: int,
                         min_mean_q: float):
    """
    Run de-multiplexing and trimming of a single fastq file using porechop
    :param read_fastq: path to fastq file to be trimmed and demultiplexed
    :param custom_adapters: path to fasta file with custom adapters
    :param custom_barcodes: path to fasta file with custom barcodes
    :param demultiplexing_dir: directory to save demultiplexed reads
    :param filtered_dir: directory to save filtered reads
    :param min_length: minimum read length to keep during quality filtering
    :param min_mean_q: minimum mean quality to keep during read filtering
    """

    porechop_custom.main(input=read_fastq,
                         barcode_dir=demultiplexing_dir,
                         custom_adapters=custom_adapters,
                         custom_barcodes=custom_barcodes,
                         require_two_barcodes=True,
                         discard_unassigned=True)
    # filtlong quality filtering
    output_files = []

    for fastq in demultiplexing_dir.glob('*.fastq'):
        filtlong_command = ['filtlong',
                            '--min_length', min_length,
                            '--min_mean_q', min_mean_q,
                            fastq]

        filtlong_output = run_external(filtlong_command, stdout='capture')
        filtered_reads = f'{fastq.stem}.fastq'
        filtered_fastq = filtered_dir.joinpath(filtered_reads)
        with filtered_fastq.open('w') as handle:
            handle.write(filtlong_output.decode())
            output_files.append(filtered_fastq)
    return output_files

def map_to_reference(reference: Path,
                     read_fastq: Path,
                     minimap_sorted_dir: Path,
                     minumap_filtered_dir: Path):
    """
    Map reads to reference gene sequence
    using split-read mapping with minimap2
    and post-process the files
    by sorting and filtering alignments
    :param reference: path to reference fasta file with gene sequence
    :param read_fastq: path to fastq file with reads to map
    :param minimap_sorted_dir: directory to save sorted bam files
    :param minumap_filtered_dir: directory to save filtered bam files
    """
    minimap_command = ['minimap2', '-ax', 'splice', '-uf', '-k14', '--splice-flank=no',  # '-p', 0.9, '--secondary=no',
                       reference, read_fastq]
    sorted_bam = minimap_sorted_dir.joinpath(read_fastq.stem + '.bam')
    sorting_command = ['samtools', 'sort', '-O', 'bam', '-o', sorted_bam]
    minimap_output = run_external(minimap_command, stdout='capture')
    run_external(sorting_command, stdin=minimap_output)

    filtered_file = minumap_filtered_dir.joinpath(read_fastq.stem + '.bam')
    filter_command = ['samtools', 'view', '-b', '-F', '4', '-q', '30', sorted_bam, '-o', filtered_file]
    run_external(filter_command)
    index_file = filtered_file.with_suffix('.bai')
    index_command = ['samtools', 'index', '-o', index_file, filtered_file]
    run_external(index_command)
    return filtered_file


def count_mapped_reads(bam: Path) -> int:
    """
    Count number of reads mapped in the bam file
    :param bam: path to bam file to count
    :return: number of mapped reads
    """
    count_command = ['samtools', 'view', '-c', '-F', '260', bam]
    mapped_reads = int(run_external(count_command, stdout='capture').decode().strip())
    return mapped_reads


def splice_sites(spliced_bam: Path,
                 raw_gff_dir: Path,
                 mapped_reads: int,
                 output: AnnotationBase,
                 minimal_fraction: float = 0.1,
                 minimal_coverage: int = 100):
    """
    Use spliced_bam2gff to extract intron borders from spliced bam file and convert them to gff format
    :param spliced_bam: path to spliced bam file
    :param raw_gff_dir: directory to save raw gff files
    :param mapped_reads: number of mapped reads in the bam file
    :param output: AnnotationBase object to store introns
    :param minimal_fraction: minimal fraction of reads supporting the intron required to report it in the final gff
    :param minimal_coverage: minimal number of reads supporting the intron required to report it in the final gff
    """
    spliced_bam2gff_command = ['spliced_bam2gff', '-M', spliced_bam]
    raw_gff = raw_gff_dir.joinpath(f'{spliced_bam.stem}.gff')
    raw_gff_lines = run_external(spliced_bam2gff_command, stdout='capture')
    with raw_gff.open('w') as handle:
        handle.write(raw_gff_lines.decode())
    annotation_lines = [l for l in raw_gff_lines.decode().splitlines() if not l.startswith('#') and l.strip()]
    genes = []
    exons = None
    while annotation_lines:
        line = annotation_lines.pop(0)
        if line.split('\t')[2] == 'mRNA':
            if exons:
                genes.append(Gene.from_exons(exons))
            exons = []
        elif line.split('\t')[2] == 'exon':
            exons.append(Exon.from_gff_line(line, in_attr_separator=' '))
    if exons:  # add the last gene
        genes.append(Gene.from_exons(exons))
    intron_versions = {}
    for gene in genes:
        possible_introns = gene.children(Intron)
        for p_intron in possible_introns:
            if (p_intron.start, p_intron.end, p_intron.strand) in intron_versions:
                intron_versions[(p_intron.start, p_intron.end, p_intron.strand)].score += 1
            else:
                intron_versions[(p_intron.start, p_intron.end, p_intron.strand)] = Intron(seq_id=p_intron.seq_id,
                                                                                          strand=p_intron.strand,
                                                                                          model_id=f'pinfish',
                                                                                          model_name=spliced_bam.stem,
                                                                                          start=p_intron.start,
                                                                                          end=p_intron.end,
                                                                                          score=1,
                                                                                          evalue=0,
                                                                                          method='averaged_pinfish')
    for intron_variant in intron_versions.values():
        supporting_read_fraction = intron_variant.score / mapped_reads
        if intron_variant.score > minimal_coverage and supporting_read_fraction > minimal_fraction:
            intron_variant.model_id = f'{intron_variant.score} / {mapped_reads}'
            intron_variant.evalue = supporting_read_fraction
            output.annotate(intron_variant)
    return output


def subsample_bam(bam: Path,
                  subsampled_bam: Path,
                  max_reads: int,
                  mapped_reads: int):
    """
    Subsample bam file to a given number of reads
    :param bam: path to bam file
    :param max_reads: number of reads to subsample
    :param mapped_reads: number of mapped reads in the bam file
    """
    assert max_reads < mapped_reads, f'Cannot subsample {max_reads} out of {mapped_reads} reads'

    fraction_to_keep = max_reads / mapped_reads
    subsample_command = ['samtools', 'view', '--subsample', fraction_to_keep, '--subsample-seed', randint(0, 100000), '-b', '-o', subsampled_bam, bam]
    run_external(subsample_command)
    index_command = ['samtools', 'index', '-o', subsampled_bam.with_suffix('.bai'), subsampled_bam]
    run_external(index_command)
    return subsampled_bam


if __name__ == '__main__':
    rnaseq_analysis()

# USAGE EXAMPLE:
# ./confirm_introns.py -ca adapter_promoter -cb barcodes_r -fq 100k_sample.fastq -rd rnaseq_references -o confirm_introns_test -s '__' -mr 1000
