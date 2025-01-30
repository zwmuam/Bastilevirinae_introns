#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script annotates sequences in a fasta file with protein domains and RNA families.
It was used to re-annotate Bastillevirinae and GIISSD intron sequences
with PHROG, RFAM and GSIID ribozyme and homing endonuclease families.
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

import json as js
from pathlib import Path
from typing import Dict, Tuple

import click
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqRecord
from annotations import AnnotationBase, AnnotationTrack, DnaAlignment, InfernalAlignment

from tweaks import run_external, logger, log_format, default_threads


@click.command(no_args_is_help=True)
@click.option("-f", "--fasta",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='fasta file input genome sequences')
@click.option("-c", "--cms",
              required=False,
              default=Path(__file__).parent.joinpath('databases', 'Merged.1.GISSD_IRFAM.cm'),
              type=click.Path(path_type=Path),
              help='InfeRNAl cm database with family models e.g. GSIID/RFAM')
@click.option("-h", "--hmm",
              required=False,
              default=Path(__file__).parent.joinpath('databases', 'Phrogs4_HMMer3.hmm'),
              type=click.Path(path_type=Path),
              help='HMMER hmm database with family models e.g. PHROGS')
@click.option("-p", "--mincmscore",
              required=False,
              default=20,
              type=int,
              help='minimum score for infernal hits')
@click.option("-q", "--minhmmscore",  # https://www.biorxiv.org/content/10.1101/2021.06.24.449764v2.full
              required=False,
              default=20,
              type=int,
              help='minimum score for HMMER hits')
@click.option("-o", "--out",
              required=True,
              type=click.Path(path_type=Path),
              help='output directory for the results')
@click.option("-t", "--threads",
              required=False,
              default=default_threads,
              type=int,
              help=f'number of CPU threads to use [default: {default_threads}]')
@click.option("-x", "--cmtool",
              required=False,
              default='cmsearch',
              type=click.Choice(['cmsearch', 'cmscan']),
              help='InfeRNAl program used to search for intron RNA motifs [cmscan/cmsearch]')
@click.option("-y", "--hmmtool",
              required=False,
              default='hmmscan',
              type=click.Choice(['hmmsearch', 'hmmscan']),
              help='HMMER program used to resolve gene structure [hmmscan/hmmsearch]')
@click.option("-r", "--phrog_table",
              required=False,
              default=Path(__file__).parent.joinpath('databases', 'phrog_annot_v4.tsv'),
              type=click.Path(exists=True, path_type=Path),
              help='table with PHROG annotations')
def annotate_introns(fasta: Path,
                     cms: Path,
                     mincmscore: int,
                     cmtool: str,
                     hmm: Path,
                     minhmmscore: int,
                     out: Path,
                     hmmtool: str,
                     threads: int,
                     phrog_table: Path):
    """
    This script annotates sequences in a fasta file with protein domains and RNA families.
    It uses Infernal to detect RNA families and HMMer3 to detect protein domains and
    combines the results into a single gff file and annotation files (one per sequence) compatible with r2dt.
    :param fasta: fasta file with input sequences to annotate
    :param cms: InfeRNAl cm database with family models e.g. GSIID/RFAM (default: Merged.1.GISSD_IRFAM.cm from the databases dir)
    :param mincmscore: minimum score for InfeRNAl hits (default: 20)
    :param cmtool: InfeRNAl program used to resolve gene structure [cmsearch/cmalign] (default: cmsearch)
    :param hmm: HMMER hmm database with family models e.g. PHROGS (default: phrogs4 from the databases dir)
    :param minhmmscore: minimum score for HMMer3 hits (default: 20)
    :param out: output directory for the results (included gff and ubdirectory with r2dt files)
    :param hmmtool: HMMER program used to resolve gene structure [hmmscan/hmmsearch] (default: hmmscan)
    :param threads: number of CPU threads to use (default: all available - 1)
    :param phrog_table: table with PHROG annotations (default: phrog_annot_v4.tsv from the databases dir)
    """
    # set up tmp intermediate directory and logger

    tmp_dir = out.joinpath(f'annotation_tmp')
    # tmp_dir.mkdir(parents=True)

    logger.add(tmp_dir.joinpath('annotation.log').as_posix(), format=log_format)
    logger.info(f'Started with the following parameters:')
    for k, v in locals().items():
        logger.info(f'{k}: {v}')
    logger.opt(raw=True).debug("\n")

    forward_faa, dna_lengths = translate_fna(fasta,
                                             tmp_dir.joinpath(f'{fasta.stem}.f_translation.faa'))

    hmmer_domtblout = tmp_dir.joinpath(f'{forward_faa.stem}_X_{hmm.stem}.{hmmtool}.domtblout')
    hmmer_command = [hmmtool,
                     '-o', '/dev/null',
                     '--domtblout', hmmer_domtblout,
                     '--noali',
                     '--cpu', threads,
                     hmm, forward_faa.as_posix()]
    if not hmmer_domtblout.exists():
        run_external(hmmer_command, stdout='supress')
    hmm_alignments = AnnotationBase.from_hmmer(hmmer_domtblout)
    hmm_alignments = hmm_alignments.filter_score(threshold=minhmmscore)

    if phrog_table:
        phrog_df = pd.read_table(phrog_table, usecols=['phrog', 'annot'])
        phrog_df['phrog'] = phrog_df['phrog'].apply(lambda x: f'phrog_{x}')
        phrog_dict = dict(zip(phrog_df['phrog'], phrog_df['annot']))
        hmm_alignments.get_model_names(phrog_dict)

    hammer_in_RNA = AnnotationBase()
    for translation_id, alignments in hmm_alignments.items():
        dna_length = dna_lengths[translation_id.split('___')[0]]
        for alignment in alignments:
            dna_alignment = alignment.to_dna(dna_len=dna_length)
            hammer_in_RNA.annotate(dna_alignment)

    final_annotation = hammer_in_RNA.cull()

    infernal_tblout = tmp_dir.joinpath('infernal.tblout')
    infernal_command = [cmtool,
                        '-o', '/dev/null',
                        '--tblout', infernal_tblout,
                        '--noali', '--anytrunc',
                        '--cpu', threads,
                        cms, fasta.as_posix()]
    if not infernal_tblout.exists():
        run_external(infernal_command, stdout='supress')
    infernal_alignments = AnnotationBase.from_infernal(infernal_tblout, program=cmtool)
    infernal_alignments = infernal_alignments.filter_score(threshold=mincmscore)
    culled_infernal_alignments = infernal_alignments.cull()

    for seq_id, alignments in culled_infernal_alignments.items():
        for alignment in alignments:
            final_annotation.annotate(alignment)

    r2dt_dir = out.joinpath('r2dt')
    for seq_id, seq, seq_annotations in final_annotation.with_sequences(fasta):
        seq_annotations = final_annotation[seq_id]
        r2td = r2dt_annotation(seq_annotations, seq)
        r2td_file = r2dt_dir.joinpath(f'{seq_id}.tsv')
        r2td.to_csv(r2td_file, sep='\t', index=False)

    final_annotation.save_gff(out.joinpath('annotation.gff'))


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


def r2dt_annotation(seq_annotation: AnnotationTrack,
                    seq: str):
    """
    Generate r2dt style annotation file from gff and fasta files
    :param seq_annotation: track with annotations of the single sequence (e.g. intron)
    :param seq: sequence (as a string)
    :return: r2dt style annotation
    """
    tranck_len = len(seq) + 2  # shifted to account for 5' and 3' symbols
    colour_tracks = {DnaAlignment: (np.zeros(tranck_len), 1),
                     InfernalAlignment: (np.zeros(tranck_len), 2)}
    for a in seq_annotation:
        track, colour = colour_tracks[a.__class__]
        track[a.start:a.end] = colour
    # sum the two tracks to get the final annotation
    colour_track = np.sum([track for track, _ in colour_tracks.values()], axis=0)

    annot_df = pd.DataFrame({'residue_index': np.arange(tranck_len),
                             'residue_name': ['5\''] + list(seq) + ['3\''],
                             # add the 0 to both ends to match the annotation
                             'annotation_code': colour_track.astype(int)})

    return annot_df


# colouring for the r2dt annotation
colour_dict = {'coloring': {'annotation_code': {'label': 'Intron parts',
                                                'values': {'0': 'rgb(255, 255, 255)',  # white (no annotation)
                                                           '1': 'rgb(240, 80, 80)',  # red (CDS/ORF)
                                                           '2': 'rgb(80, 80, 240)',  # blue (RFAM match)
                                                           '3': 'rgb(200, 100, 200)',  # purple (both RFAM and CDS)
                                                           '*': 'rgb(255, 255, 255)'}}}}  # white (no annotation)


def export_colour_dict(colour_dict: dict,
                       colour_dict_path: Path):
    """
    Export colour dictionary to a json file
    :param colour_dict: dictionary with colouring information
    :param colour_dict_path: output json file
    """
    with open(colour_dict_path, 'w') as f:
        js.dump(colour_dict, f, indent=4)


if __name__ == '__main__':
    annotate_introns()

"""
USAGE EXAMPLE:
reannotate_introns.py
-f GIISSD_and_Bastille.fasta 
-o Annotated_GIISSD_and_Bastille

the output includes a gff file with the annotation and a tsv file for each input sequence with the annotation in r2dt format
"""
