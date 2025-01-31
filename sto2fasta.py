#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple utility script heavily based on Biopython used to convert .sto files to .fasta alignment format.
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

from pathlib import Path

import click
from Bio import SeqIO


@click.command(no_args_is_help=True)
@click.option("-s", "--stockholm",
              required=True,
              type=click.Path(exists=True, path_type=Path),
              help='stockholm (sto) alignment file to convert')
def convert_and_save(stockholm: Path):
    """
    Convert stockholm file to fasta format and write to a new file in the same directory
    
    Usage example: sto2fasta.py -s alignment.sto
    """
    out_file = stockholm.with_suffix('.fasta')
    records = SeqIO.parse(stockholm, "stockholm")
    SeqIO.write(records, out_file, "fasta")
    print(f'Fasta file written to {out_file}')


if __name__ == '__main__':
    convert_and_save()
