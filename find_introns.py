#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "jakub.barylski@gmail.com"
__status__ = "development"

from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
from Bio import SeqIO, SeqRecord

from annotations import Annotation, Gene, HmmerAlignment, AnnotationBase, Intron
from intron_statistics import (plot_intron_distribution, plot_intron_lengths, introns_in_genomes,
                               taxonomic_annotations, intron_architectures, nuclease_table)
from tweaks import run_external, logger, log_format, default_threads


# import click
# @click.command(no_args_is_help=True)
# @click.option("-f", "--fasta",
#               required=True,
#               type=click.Path(exists=True, path_type=Path),
#               help='fasta file input genome sequences')
# @click.option("-c", "--cms",
#               required=False,
#               default=Path(__file__).parent.joinpath('databases', 'Merged.1.GISSD_IRFAM.cm'),
#               type=click.Path(path_type=Path),
#               help='InfeRNAl cm database with family models e.g. GSIID/RFAM')
# @click.option("-h", "--hmm",
#               required=False,
#               default=Path(__file__).parent.joinpath('databases', 'Phrogs4_HMMer3.hmm'),
#               type=click.Path(path_type=Path),
#               help='HMMER hmm database with family models e.g. PHROGS')
# @click.option("-p", "--mincmscore",
#               required=False,
#               default=20,
#               type=int,
#               help='minimum score for infernal hits')
# @click.option("-q", "--minhmmscore",  # https://www.biorxiv.org/content/10.1101/2021.06.24.449764v2.full
#               required=False,
#               default=20,
#               type=int,
#               help='minimum score for HMMER hits')
# # Searching sequence databases for functional homologs using profile HMMs: how to set bit-score thresholds?
# # (...). Bit scores were used as thresholds rather than E-values since they remain the same irrespective of the size of the database searched.
# @click.option("-o", "--out",
#               required=True,
#               type=click.Path(path_type=Path),
#               help='output directory for the results')
# @click.option("-n", "--context",
#               required=False,
#               default=2500,
#               type=int,
#               help='size of the neighborhood for gene structure analysis (from each side)')
# @click.option("-t", "--threads",
#               required=False,
#               default=default_threads,
#               type=int,
#               help=f'number of CPU threads to use [default: {default_threads}]')
# @click.option("-t", "--cmtblout",
#               required=False,
#               type=click.Path(exists=True, path_type=Path),
#               help='InfeRNAl tblout file with cmscan/cmsearch results (optional, skips infernal search)')
# @click.option("-x", "--cmtool",
#               required=False,
#               default='cmsearch',
#               type=click.Choice(['cmsearch', 'cmscan']),
#               help='InfeRNAl program used to search for intron RNA motifs [cmscan/cmsearch]')
# @click.option("-y", "--hmmtool",
#               required=False,
#               default='hmmscan',
#               type=click.Choice(['hmmsearch', 'hmmscan']),
#               help='HMMER program used to resolve for gene structure [hmmscan/hmmsearch]')
# @click.option("-g", "--gff",
#               required=False,
#               type=click.Path(exists=True, path_type=Path),
#               help='GFF file with annotations (optional, skips infernal search)')
# @click.option("-r", "--phrog_table",
#               required=False,
#               default=Path(__file__).parent.joinpath('databases', 'phrog_annot_v4.tsv'),
#               type=click.Path(exists=True, path_type=Path),
#               help='table with PHROG annotations')
# @click.option("-w", "--taxon_table",
#               required=False,
#               type=click.Path(exists=True, path_type=Path),
#               help='table with taxonomic annotations')
def find_introns(fasta: Path,
                 cms: Path,
                 mincmscore: int,
                 hmm: Path,
                 minhmmscore: int,
                 out: Path,
                 context: int,
                 threads: int,
                 cmtblout: Path,
                 cmtool: str,
                 hmmtool: str,
                 gff: Path,
                 phrog_table: Path,
                 taxon_table: Path):
    """
    The script searching for intron-split genes in phage genomes.
    It uses Infernal to detect sequences similar to covariance models from the custom database
    compiled from RFAM (general database of structural RNAs) and GISSD (specialized group I intron database).
    After the initial search, low-scoring hits or hits embedded in higher-scouring alignments are removed.
    Then, regions 2500bp-long regions flanking hits are extracted (merging any overlaps) and
    aligned with protein families to resolve structure of split genes.
    Additionally, the script reports the basic statistics of the search
    and plots of intron distribution and lengths.
    :param fasta: fasta file input genome sequences
    :param cms: InfeRNAl cm database with family models e.g. GSIID/RFAM
    :param mincmscore: minimum score for InfeRNAl hits
    :param hmm: HMMER hmm database with family models e.g. PHROGS
    :param minhmmscore: minimum score for HMMer3 hits
    :param out: output directory for the results
    :param context: size of the neighborhood for gene structure analysis (from each side)
    :param threads: number of CPU threads to use
    :param cmtblout: InfeRNAl tblout file with cmscan/cmsearch results (optional, skips infernal search)
    :param cmtool: InfeRNAl program used to search for intron RNA motifs [cmscan/cmsearch]
    :param hmmtool: HMMER program used to resolve for gene structure [hmmscan/hmmsearch]
    :param gff: GFF file with annotations (optional, skips infernal search)
    :param phrog_table: table with PHROG annotations
    :param taxon_table: table with taxonomic annotations
    USAGE EXAMPLE:
    ./intron_analysis.py -f /home/michalstanoch/data/Introny_bastille_paper/genomy/RefSeq_Bastillevirinae_23062023.fasta
    -c /home/michalstanoch/data/Introny_bastille_paper/intron_repo/databases/Merged.1.GISSD_IRFAM.cm
    -h /home/michalstanoch/data/Introny_bastille_paper/intron_repo/databases/Phrogs4_HMMer3.hmm
    -o /home/michalstanoch/data/Introny_bastille_paper/intron_repo/final_paper_run
    """
    # set up main log in the output directory
    logger.add(out.joinpath('intron_analysis.log').as_posix(), format=log_format)
    logger.info(f'Started with the following parameters:')
    for k, v in locals().items():
        logger.info(f'{k}: {v}')
    logger.opt(raw=True).debug("\n")

    out.mkdir(parents=True, exist_ok=True)  # todo remove exist_ok

    if cmtblout:
        infernal_tblout = cmtblout
    else:
        infernal_tblout = out.joinpath('infernal.tblout')
        infernal_command = [cmtool,
                            '-o', '/dev/null',
                            '--tblout', infernal_tblout,
                            '--noali', '--anytrunc',
                            '--cpu', threads,
                            cms, fasta]
        if not infernal_tblout.exists():  # todo remove
            run_external(infernal_command)

    infernal_alignments = AnnotationBase.from_infernal(infernal_tblout, program=cmtool)
    infernal_alignments = infernal_alignments.filter_score(threshold=mincmscore)
    culled_infernal_alignments = infernal_alignments.cull()
    culled_infernal_alignments.save_gff(out.joinpath('infernal.gff'))

    neighborhoods = culled_infernal_alignments.contexts(context)
    merged_neighborhoods = neighborhoods.merge_overlapping()
    neighborhood_fasta = out.joinpath(f'neighborhoods_{context}.fasta')
    merged_neighborhoods.annotation_sequences(input_fasta=fasta,
                                              output_fasta=neighborhood_fasta)

    forward_faa, seq_lengths = translate_fna(neighborhood_fasta,
                                             out.joinpath(f'{neighborhood_fasta.stem}.f_translation.faa'))

    if gff:
        hmm_alignments = AnnotationBase.from_gff(gff, HmmerAlignment)
    else:
        master_domtblout = out.joinpath(f'{forward_faa.stem}_X_{hmm.stem}.{hmmtool}.domtblout')
        hmmer_command = [hmmtool,
                         '-o', '/dev/null',
                         '--domtblout', master_domtblout,
                         '--noali',
                         '--cpu', threads,
                         hmm, forward_faa.as_posix()]
        if not master_domtblout.exists():  # todo remove
            run_external(hmmer_command, stdout='supress')
        hmm_alignments = AnnotationBase.from_hmmer(master_domtblout)

    hmm_alignments = hmm_alignments.filter_score(threshold=minhmmscore)
    if phrog_table:
        phrog_df = pd.read_table(phrog_table, usecols=['phrog', 'annot'])
        phrog_df['phrog'] = phrog_df['phrog'].apply(lambda x: f'phrog_{x}')
        phrog_dict = dict(zip(phrog_df['phrog'], phrog_df['annot']))
        hmm_alignments.get_model_names(phrog_dict)

    gene_structure = resolve_gene_structure(hmm_alignments,
                                            seq_lengths)

    gene_structure.save_gff(out.joinpath('gene_structure.gff'))

    if phrog_table:

        distribution_plot = plot_intron_distribution(sequence_annotations=gene_structure,
                                                     id2name_dict=phrog_dict)
        distribution_plot.write_image(out.joinpath('intron_distribution.svg'))
        length_plot = plot_intron_lengths(sequence_annotations=gene_structure,
                                          id2name_dict=phrog_dict)
        length_plot.write_image(out.joinpath('intron_lengths.svg'))
        genome2intron_table = introns_in_genomes(sequence_annotations=gene_structure,
                                                 id2name_dict=phrog_dict)
        if taxon_table:
            taxon_table = pd.read_excel(taxon_table, header=0, index_col=0)
            genome2intron_table, genus_summary = taxonomic_annotations(genome2intron_table,
                                                                       taxon_table=taxon_table)
            genus_summary.to_excel(out.joinpath('genus_summary.xlsx'))

            nuclease_tab = nuclease_table(gene_structure, taxon_table=taxon_table)
            nuclease_tab.to_excel(out.joinpath('Intron_CDS_table.xlsx'))

            architecture_table = intron_architectures(gene_annotations=gene_structure,
                                                      cm_annotations=culled_infernal_alignments,
                                                      taxon_table=taxon_table)
            architecture_table.to_excel(out.joinpath('architecture_table.xlsx'))
        genome2intron_table.to_excel(out.joinpath('intron_table.xlsx'))


def translate_fna(in_fna: Path,
                  out_faa: Path) -> Tuple[Path, Dict[str, int]]:
    """
    Translate DNA sequences to protein sequences in three forward reading frames
    :param in_fna: fasta file with input DNA sequences
    :param out_faa: output fasta file with protein sequences
    :return: output fasta file with protein sequences and dictionary with DNA sequence lengths
    """
    sequences = SeqIO.index(in_fna.as_posix(), 'fasta')
    translations = []
    lengths = {}
    for seq_id, seq in sequences.items():
        lengths[seq_id] = len(seq)
        for frame in range(1, 4):
            ft = SeqRecord.SeqRecord(seq.seq[frame - 1:]).translate(table=11, to_stop=False, stop_symbol='*')
            ft.id = f'{seq_id}___{frame}'
            translations.append(ft)
    SeqIO.write(translations, out_faa.as_posix(), 'fasta')
    return out_faa, lengths


def resolve_gene_structure(hmm_alignments: AnnotationBase,
                           dna_lengths: Dict[str, int]) -> AnnotationBase:
    """
    Parse domtblout file to find hmm alignments and use them to find exons and introns
    :param hmm_alignments: dict of HmmerAlignment objects
    :param dna_lengths: dict of DNA sequence lengths
    :return: list of annotations representing detected genes exons and introns
    """

    alignments_in_dna = AnnotationBase()
    for translation_id, alignments in hmm_alignments.items():
        dna_length = dna_lengths[translation_id.split('___')[0]]
        for alignment in alignments:
            dna_alignment = alignment.to_dna(dna_len=dna_length)
            alignments_in_dna.annotate(dna_alignment)

    resolved_genes = AnnotationBase()
    for dna_id, model_domains in alignments_in_dna.to_dict('model_id').items():
        try:
            best_split_gene = Gene.top_scoring([Gene.from_exons(ds) for ds in model_domains.values() if len(ds) > 1])
            if not best_split_gene.children(Intron):
                logger.warning(f'No introns found in {dna_id}')
            [i.dock_overlapping(alignments_in_dna[dna_id], culled=True, overlap_threshold=0.5) for i in best_split_gene.children(Intron)]
            anchor = Annotation.from_id_string(dna_id)
            aligned_gene = best_split_gene.absolute_coordinates(parent_annotation=anchor)
            resolved_genes.annotate(aligned_gene)
        except IndexError:
            logger.warning(f'Could not find split gene in {dna_id}')

    return resolved_genes


if __name__ == '__main__':
    find_introns()

# USAGE EXAMPLE:
# ./intron_analysis.py -f /home/michalstanoch/data/Introny_bastille_paper/genomy/RefSeq_Bastillevirinae_23062023.fasta
# -c /home/michalstanoch/data/Introny_bastille_paper/intron_repo/databases/Merged.1.GISSD_IRFAM.cm
# -h /home/michalstanoch/data/Introny_bastille_paper/intron_repo/databases/Phrogs4_HMMer3.hmm
# -o /home/michalstanoch/data/Introny_bastille_paper/intron_repo/final_paper_run
# -w /home/michalstanoch/data/Introny_bastille_paper/INTRON_STATS/taxonomy.xlsx
