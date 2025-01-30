#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
todo: add docstring
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "MIT"
__version__ = "0.1"
__email__ = "jakub.barylski@gmail.com"
__status__ = "development"

# imports
from collections import defaultdict, Counter
from copy import deepcopy
from typing import Dict, Tuple

import pandas as pd
import plotly.express as px

from annotations import AnnotationBase, Gene, Intron, DnaAlignment, InfernalAlignment
from tweaks import logger


def plot_intron_distribution(sequence_annotations: AnnotationBase,
                             id2name_dict: Dict[str, str]):
    """
    Plot how many introns are in genes aligned to each HMM
    :param sequence_annotations: AnnotationBase with genes and nested introns
    :param id2name_dict: dictionary with id as the key and name as the value
    """
    to_count = []
    for seq_id, annotations in sequence_annotations.items():
        for gene in annotations:
            assert isinstance(gene, Gene)
            outer_hmm = gene.model_definition()
            n_introns = len(gene.children(Intron))
            to_count.append((outer_hmm, n_introns))
    counted = Counter(to_count).items()
    # plot stacked bar chart of intron distribution
    out_frame = pd.DataFrame([{'HMM': h, 'n_introns': i, 'n_genes': c} for (h, i), c in counted])
    out_frame['n_introns'] = out_frame['n_introns'].astype(str)
    fig = px.bar(out_frame, x='HMM', y='n_genes', color='n_introns', title='Intron Distribution',
                 color_discrete_map={'1': 'blue', '2': 'purple', '3': 'red', '4': 'black', '5': 'pink'})
    return fig


def plot_intron_lengths(sequence_annotations: AnnotationBase,
                        id2name_dict: Dict[str, str]):
    """
    Generate a plot of intron length distributions for each homing nuclease (or no homing nuclease)
    :param sequence_annotations: AnnotationBase with genes and nested introns
    :param id2name_dict: dictionary with id as the key and name as the value
    """
    genehmm2intronlengths = []
    for seq_id, annotations in sequence_annotations.items():
        for gene in annotations:
            for intron in gene.children(Intron):
                inner_hmms = [cds.model_definition() for cds in intron.children(DnaAlignment)] if intron.children(DnaAlignment) else ['No_CDS']
                if len(inner_hmms) > 1:
                    if len(set(inner_hmms)) > 1:
                        hmm_string = ', '.join(inner_hmms)
                        logger.warning(f'Intron {intron} contain multiple inner genes: {hmm_string}')
                    else:
                        logger.warning(f'Intron {intron} contain multiple inner ORFs from the same family: {inner_hmms} '
                                       f'this may be gene duplication, missassembly, frameshift or a broken gene')
                        hmm_string, = set(inner_hmms)
                else:
                    hmm_string, = inner_hmms
                intron_length = len(intron)
                genehmm2intronlengths.append({'HMM': hmm_string, 'Intron Length': intron_length})
    plot_df = pd.DataFrame.from_records(genehmm2intronlengths)
    try:
        fig = px.strip(plot_df, x='HMM',
                       y='Intron Length',
                       color='HMM',
                       title='Intron Length Distribution',
                       width=1000,
                       height=800,
                       )
        fig.update_xaxes(tickangle=45)

    except:
        raise ValueError(plot_df)
    return fig


def introns_in_genomes(sequence_annotations: AnnotationBase,
                       id2name_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Create a table of intron counts in each taxon
    :param sequence_annotations: AnnotationBase with genes and nested introns
    :param id2name_dict: dictionary with id as the key and name as the value
    :return: pandas DataFrame with intron counts and names
    """

    phage2gene = defaultdict(list)

    for seq_id, annotations in sequence_annotations.items():
        for gene in annotations:
            phage_id = seq_id.split('|')[0]
            phage2gene[phage_id].append(gene)

    records = []
    for phage_id, genes in phage2gene.items():
        provisional_names = [g.model_definition() for g in genes]
        gene_names = []
        homing_nucleases = []
        for gene, gene_name in zip(genes, provisional_names):
            if len(gene.children(Intron)) > 1:
                gene_name = f'{gene_name} ({len(gene.children(Intron))} introns)'
            else:
                gene_name = gene_name
            for intron in gene.children(Intron):
                homing_nucleases.extend([cds.model_name for cds in intron.children(DnaAlignment)])
            gene_names.append(gene_name)
        n_introns = sum([len(g.children(Intron)) for g in genes])
        hnh_count = Counter(homing_nucleases).most_common()
        hnh_names = []
        for hnh, count in hnh_count:
            if count > 1:
                hnh_name = f'{hnh} ({count} CDSs)'
            else:
                hnh_name = hnh
            hnh_names.append(hnh_name)
        records.append({'phage_id': phage_id, 'n_genes': len(genes), 'n_introns': n_introns,
                        'gene_names': '\n'.join(gene_names), 'hnh_names': '\n'.join(hnh_names)})
    genome_intron_table = pd.DataFrame.from_records(records)
    genome_intron_table = genome_intron_table.set_index('phage_id')
    return genome_intron_table


def intron_architectures(gene_annotations: AnnotationBase,
                         cm_annotations: AnnotationBase,
                         taxon_table: pd.DataFrame) -> pd.DataFrame:
    """
    Create a table of intron architectures with internal and external gene names
    as well as the covariance models of the intron ribozyme
    :param gene_annotations: AnnotationBase with genes and nested introns
    :param cm_annotations: AnnotationBase with covariance models
    :return: pandas DataFrame with intron architectures
    """
    phage_name_dict = taxon_table['organism'].to_dict()
    phage_name_dict = {k: v.split(' ')[-1] for k, v in phage_name_dict.items()}
    # dock cm annotations to gene annotations
    merged_annotations = deepcopy(gene_annotations)
    # dock cms to overlapping introns
    introns = merged_annotations.children(Intron)
    introns.dock_overlapping(cm_annotations)
    architectures = {}
    for seq_id, annotations in merged_annotations.items():
        for gene in annotations:
            for intron in gene.children(Intron):
                intron = intron.relative_to(gene)
                intron.nested.sort_annotations()
                architecture_symbol = '__'.join([e.model_id for e in intron.nested])
                border_cm_5, border_cm_3 = [], []
                while intron.nested and isinstance(intron.nested[0], InfernalAlignment):
                    border_cm_5.append(intron.nested.pop(0))
                while intron.nested and isinstance(intron.nested[-1], InfernalAlignment):
                    border_cm_3.append(intron.nested.pop(-1))
                outer_hmm = gene.model_definition()
                if architecture_symbol not in architectures:
                    architectures[architecture_symbol] = {'5 CM': ', '.join([cm.model_name for cm in border_cm_5]),
                                                          'internal elements': ', '.join([e.model_definition() for e in intron.nested]),
                                                          '3 CM': ', '.join([cm.model_name for cm in border_cm_3]),
                                                          'genes': {outer_hmm: [seq_id]},
                                                          'count': 1}
                else:
                    if outer_hmm not in architectures[architecture_symbol]['genes']:
                        architectures[architecture_symbol]['genes'][outer_hmm] = [seq_id]
                    else:
                        architectures[architecture_symbol]['genes'][outer_hmm].append(seq_id)
                    architectures[architecture_symbol]['count'] += 1
    records = []
    for architecture in architectures.values():
        occurrence_strings = []
        for gene, genomes in architecture['genes'].items():
            genome_count = Counter(genomes).most_common()
            named_count = [f'{phage_name_dict[genome]} ({count} introns)' if count > 1 else f'{phage_name_dict[genome]}' for genome, count in genome_count]
            genome_string = '; '.join(named_count)
            occurrence_strings.append(f'{gene}: {genome_string}')
        architecture['genes'] = '\n'.join(occurrence_strings)
        records.append(architecture)
    architecture_table = pd.DataFrame.from_records(records)
    return architecture_table


def homing_nuclease_table(gene_annotations: AnnotationBase,
                          taxon_table: pd.DataFrame) -> pd.DataFrame:
    """
    Create a table of homing nuclease counts in each taxon and gene
    :param sequence_annotations: AnnotationBase with genes and nested introns
    :param taxon_table: pandas DataFrame with taxon names
    :return:
    """
    phage_name_dict = taxon_table['organism'].to_dict()
    phage_name_dict = {k: v.split(' ')[-1] for k, v in phage_name_dict.items()}
    nucleases = defaultdict(dict)
    for seq_id, annotations in gene_annotations.items():
        for gene in annotations:
            for intron in gene.children(Intron):
                intron_nucleases = [cds for cds in intron.children(DnaAlignment)]
                for cds in intron_nucleases:
                    if cds.model_id not in nucleases:
                        nucleases[cds.model_id] = {'genes': {seq_id: 1},
                                                   'count': 1}
                    else:
                        if seq_id not in nucleases[cds.model_id]['genes']:
                            nucleases[cds.model_id]['genes'][seq_id] = 1
                        else:
                            nucleases[cds.model_id]['genes'][seq_id] += 1
                        nucleases[cds.model_id]['count'] += 1




def taxonomic_annotations(genome_intron_table: pd.DataFrame,
                          taxon_table: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create a table of intron counts in each taxon
    :param genome_intron_table: pandas DataFrame with intron counts and names
    :param taxon_table: pandas DataFrame with taxon names
    :return: taxonomy-annotated intron table and summary table for genera
    """
    genome_intron_table['organism'] = taxon_table.loc[genome_intron_table.index, 'organism']
    genome_intron_table['species'] = taxon_table.loc[genome_intron_table.index, 'species']
    genome_intron_table['genus'] = taxon_table.loc[genome_intron_table.index, 'genus']
    total_genome_count = Counter(taxon_table['genus']).most_common()
    intron_genome_count = {genus: count for genus, count in Counter(genome_intron_table['genus']).items()}
    records = []
    intron_free_genera = []
    for genus, count in total_genome_count:
        if genus in intron_genome_count:
            records.append({'genus': genus, 'invaded genomes': f'{intron_genome_count[genus]} / {count}'})
        else:
            intron_free_genera.append(f'{genus} ({count} genomes)')
    summary_table = pd.DataFrame.from_records(records, index='genus')
    logger.info(f'No introns found in:\t{intron_free_genera}')
    genome_intron_table = genome_intron_table[['organism', 'species', 'genus', 'n_genes', 'n_introns', 'gene_names', 'hnh_names']]
    return genome_intron_table, summary_table
