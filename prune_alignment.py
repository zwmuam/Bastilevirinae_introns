#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Open fasta-formatted alignment and filter out sequences with less than N non-gap characters
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

from pathlib import Path

import click

from tweaks import parse_fasta


@click.command(no_args_is_help=True)
@click.option("-f", "--fasta",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='fasta file input alignment')
@click.option("-l", "--min_length",
              default=50,
              required=False,
              type=int,
              help='minimum length of an aligned sequence to keep')
def filter_short_alignments(fasta: Path,
                            min_length: int):
    """
    Filter out sequences with less than N aligned characters from a multiple sequence alignment
    :param fasta: path to a fasta file with MSA
    :param min_length: minimum length of an aligned sequence to keep
    """
    sequences = []
    skipped = 0
    for seq_id, seq in parse_fasta(fasta):
        non_gaps = sum(1 for c in seq if c != '-')
        if non_gaps >= min_length:
            sequences.append((seq_id, seq))
        else:
            print(f'Skipping {seq_id} with only {non_gaps} non-gap characters')
            skipped += 1
    print(f'Skipped {skipped} sequences, {len(sequences)} sequences left')
    out_file = fasta.parent.joinpath(f'{fasta.stem}.pruned.fasta')
    with out_file.open('w') as out:
        for seq_id, seq in sequences:
            out.write(f'>{seq_id}\n{seq}\n')


if __name__ == '__main__':
    filter_short_alignments()
