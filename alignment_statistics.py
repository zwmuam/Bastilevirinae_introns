#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script calculates the basic statistics for sequence alignment guided by a reference sequence/model.
It was used to analyze the alignment of group I introns in the Bastille paper and generate supplementary tables.
Most functions follow the naming convention used in the paper (intron names beggin with "IA1__", "ID2__" etc.).
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

from collections import Counter, defaultdict
from curses.ascii import isdigit
from pathlib import Path

import click
import numpy as np
import pandas as pd
from Bio import SeqIO

from tweaks import Parallel


@click.command(no_args_is_help=True)
@click.option("-f", "--fasta",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='fasta file input multiple sequence alignment')
@click.option("-u", "--unaligned",
              default=False,
              is_flag=True,
              help='Sequences are unaligned (basic fasta, skip alignment-related statistics)')
@click.option("-r", "--ref_len",
              default=251,
              required=False,
              type=int,
              help='Length of the model or sequence used as a reference for alignment.'
                   'E.g 251  fort the RF00028 group I intron model '
                   '(CLEN - length of the model in consensus residues / match states)')
@click.option("-m", "--write_matrix",
              default=False,
              is_flag=True,
              help='Write sequence identity matrix to the output file')
def alignment_stats(fasta: Path,
                    unaligned: bool,
                    ref_len: int,
                    write_matrix: bool):
    """
    Calculate basic statistics for a sequence alignment guided by a reference sequence/model (or raw fasta file)
    
    Usage example: alignment_statistics.py -f alignment.fasta -r 251 -o alignment_stats.xlsx
    """
    msa = list(SeqIO.parse(fasta, "fasta"))

    output = fasta.parent.joinpath(f'{fasta.stem}_alignment_stats.xlsx')

    with pd.ExcelWriter(output) as writer:
        intron_group_counts = count_intron_groups_in_file(msa)
        intron_group_counts.sort_index(inplace=True)
        intron_group_counts.to_excel(writer, sheet_name="group counts")

        length_stats = intron_aligned_lengths(msa, reference_len=ref_len)
        length_stats.sort_index(inplace=True)
        length_stats.to_excel(writer, sheet_name="aligned lengths")

        if not unaligned:
            matrix = identity_matrix(msa)
            if write_matrix:
                matrix.to_excel(writer, sheet_name="sequence identity matrix")

            ingroup_similarity, intergroup_similarity = group_similarity(matrix)
            ingroup_similarity.sort_index(inplace=True), intergroup_similarity.sort_index(inplace=True)
            ingroup_similarity.to_excel(writer, sheet_name="ingroup similarity")
            intergroup_similarity.to_excel(writer, sheet_name="intergroup similarity")

def minor_group(s: str) -> str:
    """
    Extracts the minor subgroup (IA1, IB2 etc.) group of the intron from its ID string
    :param s: an ID string of an intron
    :return: group of the intron
    """
    minor_t = s.strip().split("__")[0]
    if not isdigit(minor_t[-1]):
        return minor_t + '?'
    return minor_t

def major_group(s: str) -> str:
    """
    Extracts the major subgroup (IA, IB etc.) group of the intron from its ID string
    :param s: an ID string of an intron
    :return: group of the intron
    """
    mint = minor_group(s)
    # remove any numbers at the end of the string
    major_t = ''.join([i for i in mint if not i.isdigit() and i != '?'])
    return major_t + "_total"

symbols = {'nucleic': set("ACGTUN"),
           'gap': "-"}  # DNA/RNA alphabet used in the alignment


def ungapped_length(s: str) -> int:
    """
    Calculates the length of a sequence without gaps
    :param s: sequence string
    :return: length of the sequence without gaps
    """
    return sum(1 for c in s if c != symbols['gap'])


def match(c1, c2):
    """
    Checks if two characters are valid nucleic acid symbols
    :param c1: first character
    :param c2: second character
    :return: true if both characters are A/C/U/G/T, false if any of them is a gap otherwise rise an error
    """
    if c1 == symbols['gap'] and c2 == symbols['gap']:
        return None
    elif c1 == c2:
        return 1
    else:
        return 0


def trim_end_gaps(s1: str, s2: str) -> tuple[str, str]:
    """
    Trims ends of the alignment
    that are gaps in at least one of the sequences
    :param s1: first sequence string
    :param s2: second sequence string
    :return: trimmed sequences
    """
    for i in range(len(s1)):
        if s1[i] != symbols['gap'] and s2[i] != symbols['gap']:
            break
    for j in range(len(s1) - 1, 0, -1):
        if s1[j] != symbols['gap'] and s2[j] != symbols['gap']:
            break
    return s1[i:j + 1], s2[i:j + 1]


def pairwise_identity(s1: str,
                      s2: str,
                      free_end_gaps: bool = False) -> float:
    """
    Calculates the identity between two sequences without gaps
    :param s1: first sequence
    :param s2: second sequence
    :param free_end_gaps: if True end gaps are trimmed before the calculation otherwise they are treated as mismatches
                          (default: False - applies penalty for end gaps)
    :return: identity between the sequences
    """
    if free_end_gaps:
        s1, s2 = trim_end_gaps(s1, s2)
    if len(s1) != len(s2):
        raise ValueError("Sequences must be of the same length")
    identity_track = [match(c1, c2) for c1, c2 in zip(s1, s2)]
    identity_track = [i for i in identity_track if i is not None]
    if not identity_track:
        return 0
    return sum(identity_track) / len(identity_track) * 100


def identity_matrix(msa: list[SeqIO.SeqRecord]) -> pd.DataFrame:
    """
    Creates a matrix of pairwise identities (given pre-aligned sequences)
    :param msa: SeqIO-parsed alignment (list of SeqRecords)
    :return: pandas DataFrame with sequence identity values
    """

    matrix = pd.DataFrame(columns=[record.id for record in msa], index=[record.id for record in msa])

    row_jobs = Parallel(parallelized_function=_row_identity,
                        input_collection=msa,
                        kwargs={"msa": msa},
                        description="Calculating identity matrix")

    for i, row in enumerate(row_jobs.result):
        matrix.iloc[i] = row

    return matrix

def _row_identity(s1, msa: list[SeqIO.SeqRecord]):
    """
    Calculate the row of identity matrix (identities to a single sequence)
    :param s1: SeqRecord to calculate identities 
    :param msa: SeqIO-parsed alignment (list of SeqRecords)
    :return: list of identity values for the sequence
    """
    return [pairwise_identity(s1.seq, s2.seq) for s2 in msa]


def count_intron_groups_in_file(msa: list[SeqIO.SeqRecord]) -> pd.DataFrame:
    """
    Counts the number of sequences from each group in the file (e.g. group I intron subgroups)
    :param msa: SeqIO-parsed alignment (list of SeqRecords)
    :return: pandas DataFrame with counts of intron subgroups in the file
    """
    minor_grps = [minor_group(record.id) for record in msa]
    major_grps = [major_group(intron_group) for intron_group in minor_grps]
    major_grps, minor_grps = Counter(major_grps), Counter(minor_grps)
    intron_group_counts = pd.DataFrame(columns=["subgroup", "count"])
    intron_group_counts.set_index("subgroup", inplace=True)
    all_groups = major_grps + minor_grps
    for group, count in all_groups.items():
        intron_group_counts.loc[group] = count
    return intron_group_counts

def intron_aligned_lengths(msa: list[SeqIO.SeqRecord],
                           reference_len: int = None) -> pd.DataFrame:
    """
    Calculate length statistics for sequences SeqIO-parsed file
    :param msa: SeqIO-parsed alignment (list of SeqRecords)
    :param reference_len: length of the model or sequence used as an alignment reference (used to calculate coverage)
    :return: pandas DataFrame with length statistics for intron subgroups in the alignment
    """
    df = pd.DataFrame(columns=["Sequence Length"], index=[record.id for record in msa])
    for record in msa:
        df.loc[record.id] = ungapped_length(record.seq)
    if reference_len is not None:
        df["Sequence Length"] = df["Sequence Length"] / reference_len
    lengths = defaultdict(list)
    for i, row in df.iterrows():
        i_minor = minor_group(i)
        i_major = major_group(i_minor)
        if i_minor != i_major:
            lengths[i_minor].append(row["Sequence Length"])
        lengths[i_major].append(row["Sequence Length"])

    stats = pd.DataFrame(columns=["min", "max", "mean", "median", "std"])
    for intron_group, values in lengths.items():
        minim, maxim, mean = min(values), max(values), np.mean(values)
        median = np.median(values)
        std = np.std(values)
        stats.loc[intron_group] = [minim, maxim, mean, median, std]

    return stats

def group_similarity(matrix) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze similarity between and within groups in an alignment (e.g. intron subgroups)
    :param matrix: DataFrame with sequence identity values
    :return: two DataFrames, first with identity within groups, second with identity between them
    """

    major_similarities = defaultdict(list)
    minor_similarities = defaultdict(list)
    major_groups = set()
    minor_groups = set()

    for i in matrix.index:
        i_minor = minor_group(i)
        i_major = major_group(i_minor)
        major_groups.add(i_major)
        minor_groups.add(i_minor)
        for j in matrix.columns:
            if i == j:
                continue  # skip trivial self-comparisons
            j_minor = minor_group(j)
            j_major = major_group(j_minor)
            major_key = (i_major, j_major)
            identity = matrix.loc[i, j]
            major_similarities[major_key].append(identity)
            minor_key = (i_minor, j_minor)
            minor_similarities[minor_key].append(identity)

    ingroup_similarity = pd.DataFrame(columns=["min", "max", "mean", "std"])
    sorted_minor_groups = sorted(minor_groups)
    intergroup_similarity = pd.DataFrame(columns=sorted_minor_groups, index=sorted_minor_groups)

    for (i, j), values in major_similarities.items():
        if i == j:
            minim, maxim, mean = min(values), max(values), np.mean(values)
            std = np.std(values)
            ingroup_similarity.loc[i] = [minim, maxim, mean, std]

    for (i, j), values in minor_similarities.items():
        if i == j:
            minim, maxim, mean = min(values), max(values), np.mean(values)
            std = np.std(values)
            ingroup_similarity.loc[i] = [minim, maxim, mean, std]
        intergroup_similarity.loc[i, j] = np.mean(values)

    return ingroup_similarity, intergroup_similarity


if __name__ == "__main__":
    alignment_stats()
  
