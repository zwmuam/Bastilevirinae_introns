#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module used in "find_introns.py" to calculate and visualize the statistical analysis of the detected introns:
lengths, split genes, types of embedded nucleases, genomic and taxonomic distribution etc.
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "MIT"
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
                inner_cdss = [cds for cds in intron.children(DnaAlignment)]
                total_cds_length = sum([len(cds) for cds in inner_cdss])
                inner_hmms = [cds.model_definition() for cds in inner_cdss] if inner_cdss else ['No_CDS']
                if len(inner_hmms) > 1:
                    if len(set(inner_hmms)) > 1:
                        hmm_string = ', '.join(inner_hmms)
                        logger.warning(f'Intron {intron} contain multiple inner genes: {hmm_string}')
                    else:
                        logger.warning(
                            f'Intron {intron} contain multiple inner ORFs from the same family: {inner_hmms} '
                            f'this may be gene duplication, missassembly, frameshift or a broken gene')
                        hmm_string, = set(inner_hmms)
                else:
                    hmm_string, = inner_hmms
                intron_length = len(intron)
                genehmm2intronlengths.append(
                    {'HMM': hmm_string, 'Intron Length': intron_length, 'CDS Length': total_cds_length})
    plot_df = pd.DataFrame.from_records(genehmm2intronlengths)

    # calculate the Pearson correlation coefficient between intron and CDS lengths
    intron_cds_correlation = plot_df['Intron Length'].corr(plot_df['CDS Length'], method='pearson')
    # plot intron and CDS length distributions as a scatter plot
    fig = px.scatter(plot_df, x='Intron Length',
                     y='CDS Length',
                     color='HMM',
                     title=f'Intron and CDS Length Distribution (Pearson correlation: {intron_cds_correlation:.2f})',
                     width=1200, height=600)
    # increase marker size for better visibility
    fig.update_traces(marker=dict(size=10))
    return fig


def introns_in_genomes(sequence_annotations: AnnotationBase,
                       id2name_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Create a table of intron counts in each taxon
    :param sequence_annotations: AnnotationBase with genes and nested introns
    :param id2name_dict: dictionary with id as the key and name as the value`
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
                                                          'internal elements': ', '.join(
                                                              [e.model_definition() for e in intron.nested]),
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
            named_count = [f'{phage_name_dict[genome]} ({count} introns)' if count > 1 else f'{phage_name_dict[genome]}'
                           for genome, count in genome_count]
            genome_string = '; '.join(named_count)
            occurrence_strings.append(f'{gene}: {genome_string}')
        architecture['genes'] = '\n'.join(occurrence_strings)
        records.append(architecture)
    architecture_table = pd.DataFrame.from_records(records)
    return architecture_table


def nuclease_table(gene_annotations: AnnotationBase,
                   taxon_table: pd.DataFrame) -> pd.DataFrame:
    """
    Create a table of homing nuclease counts in each taxon and gene
    :param gene_annotations: AnnotationBase with genes and nested introns
    :param taxon_table: pandas DataFrame with taxon names
    :return:
    """
    phage_genus_dict = taxon_table['genus'].to_dict()
    nucleases = defaultdict(dict)

    for seq_id, annotations in gene_annotations.items():
        for gene in annotations:
            for intron in gene.children(Intron):
                intron_nucleases = [cds for cds in intron.children(DnaAlignment)]
                for cds in intron_nucleases:
                    if cds.model_id not in nucleases:
                        nucleases[cds.model_id] = {'PHROG_id': cds.model_id,
                                                   'name': cds.model_name,
                                                   'InterPro_domains': 'not implemented',
                                                   'genes': [gene.model_definition()],
                                                   'genera': [phage_genus_dict[seq_id]]}
                    else:
                        nucleases[cds.model_id]['genes'].append(gene.model_definition())
                        nucleases[cds.model_id]['genera'].append(phage_genus_dict[seq_id])

    records = []
    for nuclease in nucleases.values():
        gene_count = Counter(nuclease['genes']).most_common()
        genera_count = Counter(nuclease['genera']).most_common()
        nuclease['genes'] = ', '.join([f'{gene} ({count})' for gene, count in gene_count])
        nuclease['genera'] = ', '.join([f'{genus} ({count})' for genus, count in genera_count])
        records.append(nuclease)

    nuclease_df = pd.DataFrame.from_records(records)

    return nuclease_df


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
    genome_count = {genus: count for genus, count in Counter(genome_intron_table['genus']).items()}
    # sum [n_genes] column for row that have the same genus
    gene_count = genome_intron_table.groupby('genus')['n_genes'].sum()
    intron_count = genome_intron_table.groupby('genus')['n_introns'].sum()
    records = []
    intron_free_genera = []
    for genus, count in total_genome_count:
        if genus in genome_count:
            records.append({'genus': genus, 'invaded genomes': f'{genome_count[genus]} / {count}',
                            'invaded genes': gene_count[genus], 'introns': intron_count[genus]})
        else:
            intron_free_genera.append(f'{genus} ({count} genomes)')
    summary_table = pd.DataFrame.from_records(records, index='genus')
    logger.info(f'No introns found in:\t{intron_free_genera}')
    genome_intron_table = genome_intron_table[
        ['organism', 'species', 'genus', 'n_genes', 'n_introns', 'gene_names', 'hnh_names']]
    return genome_intron_table, summary_table
