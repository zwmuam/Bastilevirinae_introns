#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The module contains basic tools for handling sequence annotations
it is a part of the intron_analysis package
"""

__author__ = "Jakub Barylski"
__maintainer__ = "Jakub Barylski"
__license__ = "GNU GENERAL PUBLIC LICENSE"
__email__ = "jakub.barylski@gmail.com"

# imports
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Type, Union, Any

from Bio import SeqIO, SeqFeature as Feature, Seq

from tweaks import logger, extensions, parse_fasta, frantic_search


class Annotation:
    """
    Base class for all sequence annotations
    :ivar seq_id: sequence identifier
    :ivar model_id: an identifier of the model used for annotation (e.g. HMM, CM or reference sequence)
    :ivar model_name: name of the model used for annotation
    :ivar start: position of the first nucleotide/amino acid of the annotation
    :ivar end: position of the last nucleotide/amino acid of the annotation
    :ivar score: score of the annotation (e.g. score assigned by HMMer or Infernal)
    :ivar evalue: e-value of the annotation (e.g. e-value assigned by HMMer or Infernal)
    :ivar method: method used for annotation (e.g. cmnscan, hmmsearch etc.)
    :ivar strand: strand of the sequence (either '+' or '-')
    :ivar nested: list of nested annotations (e.g. exons in a gene)
    """

    parser_dict = None

    def __init__(self,
                 seq_id: str,
                 model_id: str,
                 model_name: str,
                 start: int,
                 end: int,
                 score: float,
                 evalue: float,
                 method: str,
                 strand: str):

        self.seq_id = seq_id
        self.model_id = model_id
        self.model_name = model_name
        self.method = method
        self.start, self.end = start, end

        self.strand = strand
        self.score = score
        self.evalue = evalue
        self.nested = AnnotationTrack(self.seq_id)

    def __repr__(self):
        """
        Annotation should be represented as 'ModelID [start - end] (score)'
        """
        return f'{self.seq_id}__{self.__class__.__name__}__{self.model_id}__{self.start}__{self.end}__{self.strand}'

    def model_definition(self):
        return f'{self.model_id}. {self.model_name}'

    def __len__(self) -> int:
        """
        Calculate length of domain in AA
        :return: length of domain in AA
        """
        return self.end - self.start + 1

    @classmethod
    def _find(cls,
              split_line: List[str],
              method: str,
              keys: List[str],
              dtype: callable = None) -> List[Any]:
        """
       Guide for a parser to find correct values in a line
       :param split_line: one line from a tabular file
       :param method: method (exact tool) used for annotation (e.g. cmnscan, hmmsearch etc.)
       :param keys: fields to be extracted from the line
       :return: list of values extracted from the line
       """
        extraction = [split_line[cls.parser_dict[method][f]] for f in keys]
        if dtype is not None:
            extraction = [dtype(e) for e in extraction]
        return extraction

    @classmethod
    def _from_tblout_line(cls,
                          program,
                          line,
                          unstranded=False):
        if program in ('mmseqs',):
            separator = '\t'
        else:
            separator = None
        line = line.rstrip('\n')
        split_line = line.split(separator)
        seq_id, model_id, model_name = cls._find(split_line, program, ['seq_id',
                                                                       'model_id',
                                                                       'model_name'])
        seq_start, seq_end = cls._find(split_line, program, ['seq_start',
                                                             'seq_end'], int)
        if unstranded:
            seq_strand = '+'
        else:
            seq_strand, = cls._find(split_line, program, ['seq_strand'])
            if seq_strand == '-':
                seq_start, seq_end = seq_end, seq_start
        score, evalue = cls._find(split_line, program, ['score', 'evalue'], float)
        return cls(seq_id=seq_id,
                   model_id=model_id,
                   model_name=model_name,
                   start=seq_start,
                   end=seq_end,
                   score=score,
                   evalue=evalue,
                   method=program,
                   strand=seq_strand)

    @classmethod
    def from_gff_line(cls,
                      gff_line,
                      between_attr_separator=';',
                      in_attr_separator='='):
        seq_id, method, _, start, end, score, strand, _, data = gff_line.split('\t')
        score = 0 if score == '.' else float(score)
        start, end = int(start), int(end)
        data = {k: v for k, v in
                [e.strip().split(in_attr_separator) for e in data.split(between_attr_separator) if e.strip()]}
        try:
            model_id = frantic_search(data, 'model', 'HMM', 'model_id', 'transcript_id')
        except KeyError:
            model_id = 'n.a.'
        try:
            model_name = frantic_search(data, 'model_name', 'HMM', 'description')
        except KeyError:
            model_name = 'n.a.'
        return cls(seq_id=seq_id,
                   model_id=model_id,
                   model_name=model_name,
                   start=start,
                   end=end,
                   strand=strand,
                   score=score,
                   evalue=10 if 'e_value' not in data else data['e_value'],
                   method=method)

    @classmethod
    def from_id_string(cls, repr_str: str) -> 'Annotation':
        """
        Create an instance from its string representation
        """
        seq_id, subtype, model_id, start, end, strand = repr_str.split('__')
        return cls(seq_id=seq_id,
                   model_id=model_id,
                   model_name=model_id,
                   start=int(start),
                   end=int(end),
                   score=0,
                   evalue=1,
                   method=f'repr_{subtype}',
                   strand=strand)

    @classmethod
    def retype(cls, parent_annotation: 'Annotation'):
        return cls(seq_id=parent_annotation.seq_id,
                   model_id=parent_annotation.model_id,
                   model_name=parent_annotation.model_name,
                   start=parent_annotation.start,
                   end=parent_annotation.end,
                   score=parent_annotation.score,
                   evalue=parent_annotation.evalue,
                   method=parent_annotation.method,
                   strand=parent_annotation.strand)

    def gff(self) -> List[str]:
        """
        Represent a domain as a GFF line
        :return: gff-formatted line (no newline (\n) at the end)
        """
        attributes = {'HMM': self.model_id, 'e_value': self.evalue, 'name': f'"{self.model_name}"'}
        attribute_string = ';'.join(f'{k}={v}' for k, v in attributes.items())
        gff_lines = ['\t'.join([self.seq_id,
                                self.method,
                                type(self).__name__,
                                str(self.start),
                                str(self.end),
                                str(self.score),
                                self.strand,
                                '.',
                                attribute_string])]
        for n in self.nested:
            gff_lines.extend(n.gff())
        return gff_lines

    def context(self, flank: int) -> 'Annotation':
        """
        Create a new domain that extends the current one
        by the specified number of residues on both sides
        :param flank: number of residues to extend the domain by
        :return: extended domain
        """
        start, end = self.start - flank, self.end + flank
        return Annotation(seq_id=self.seq_id,
                          model_id=self.model_id, model_name=self.model_name,
                          start=start, end=end, strand=self.strand,
                          score=self.score, evalue=self.evalue,
                          method=self.method)

    def reverse_complement(self, seq_len: int) -> 'Annotation':
        """
        Reverse complement the domain
        :return: reversed and complemented domain
        """
        start, end = seq_len - self.end + 1, seq_len - self.start + 1
        strand = {'+': '-', '-': '+'}[self.strand]

        return self.__class__(seq_id=self.seq_id,
                              model_id=self.model_id, model_name=self.model_name,
                              start=start, end=end, strand=strand,
                              score=self.score, evalue=self.evalue,
                              method=self.method)

    def overlaps(self,
                 other: 'Annotation') -> float:
        """
        Check if any of the two analysed domains overlap
        on more than threshold fraction of its length with the other
        :param other:  domain to compare with
        :return: are these domains overlapping?
        """
        if self.end >= other.start and self.start <= other.end:
            overlap_start = max(self.start, other.start)
            overlap_end = min(self.end, other.end)
            overlap = overlap_end - overlap_start + 1
            if overlap < 1:
                raise NotImplementedError()
            else:
                return max([overlap / len(e) for e in (self, other)])
        elif other.end >= self.start and other.start <= self.end:
            return other.overlaps(self)
        return 0

    def host(self,
             annotation: 'Annotation'):
        self.nested.append(annotation)

    def children(self,
                 subtype: Type['Annotation'],
                 recursive: bool = False) -> List['Annotation']:
        """
        Get all nested annotations of a given type
        :param subtype: the type of nested annotations
        :param recursive: include annotations of the same type nested in the deeper layers
        :return: list of nested annotations of a given type
        """
        typed_children = [n for n in self.nested if isinstance(n, subtype)]
        if recursive:
            for n in self.nested:
                typed_children.extend(n.children(subtype, recursive))
        return typed_children

    def dock_overlapping(self,
                         annotations: 'AnnotationTrack',
                         overlap_threshold: float = 0.5,
                         culled: bool = True):
        """
        Include nested annotations that overlap with the current one
        into 'self.nested' list
        """
        assert self.seq_id == annotations.seq_id, 'Cannot dock annotations from different sequences'
        overlapping = AnnotationTrack(self.seq_id, [a for a in annotations if self.overlaps(a) > overlap_threshold])
        if culled:
            overlapping = overlapping.cull()
        [self.host(a) for a in overlapping]

    def relative_to(self,
                    annotation: 'Annotation'):
        """
        Adjust coordinates of the annotation
        in relation to the other annotation
        (assume that star of the other annotation is start of the sequence)
        :return: adjusted annotation (with coordinates relative to the other annotation)
        """
        start, end = self.start - annotation.start + 1, self.end - annotation.start + 1
        anchored = self.__class__(seq_id=self.seq_id,
                                  model_id=self.model_id, model_name=self.model_name,
                                  start=start, end=end, strand=self.strand,
                                  score=self.score, evalue=self.evalue,
                                  method=self.method)

        if annotation.strand == '-':
            anchored = anchored.reverse_complement(seq_len=len(annotation))
        anchored.nested = AnnotationTrack(self.seq_id, [n.relative_to(annotation) for n in self.nested])
        return anchored

    def absolute_coordinates(self,
                             parent_annotation: 'Annotation'):
        """
        Adjust coordinates of the annotation
        that was anchored to the other
        (revert to original coordinates in relation to entire sequence)
        :param parent_annotation: annotation that was used as an anchor for relative coordinates
        :param recursive: adjust nested annotations
        :return: adjusted annotation (with absolute coordinates in the sequence)
        """
        if parent_annotation.strand == '-':
            a = self.reverse_complement(seq_len=len(parent_annotation))
        else:
            a = self
        start, end = a.start + parent_annotation.start - 1, a.end + parent_annotation.start - 1
        unbound = a.__class__(seq_id=a.seq_id,
                              model_id=a.model_id, model_name=a.model_name,
                              start=start, end=end, strand=a.strand,
                              score=a.score, evalue=a.evalue,
                              method=a.method)
        unbound.seq_id = parent_annotation.seq_id
        unbound.nested = AnnotationTrack(self.seq_id, [n.absolute_coordinates(parent_annotation) for n in self.nested])
        return unbound

    def merge(self, other: 'Annotation'):
        """
        Merge two overlapping annotations
        :param other: other annotation to merge with
        :return: merged annotation
        """
        start = min(self.start, other.start)
        end = max(self.end, other.end)
        score = self.score + other.score
        return self.__class__(seq_id=self.seq_id,
                              model_id=self.model_id, model_name=self.model_name,
                              start=start, end=end, strand=self.strand,
                              score=score, evalue=self.evalue,
                              method=self.method)

    @staticmethod
    def top_scoring(elem_list: List['Annotation']):
        """
        Given a list of annotations choose top-scoring one
        :param elem_list: list of two or more annotations
        :return: Top scoring domain from the cluster
        """
        ranking = sorted(elem_list, key=lambda e: e.score, reverse=True)
        return ranking[0]

    def get_sequence(self, seq: Seq.Seq) -> Seq.Seq:
        """
        Get the sequence of the annotation
        :param seq: sequence of the annotated contig
        :return: sequence of the annotation
        """
        annotated_seq = seq[self.start - 1:self.end]
        if self.strand == '-':
            annotated_seq = annotated_seq.reverse_complement()
        annotated_seq.id = annotated_seq.name = self.__repr__()
        annotated_seq.description = f'{self.model_name} score {self.score} evalue {self.evalue} {self.method}'
        return annotated_seq

    def get_model_name(self,
                       id2name_dict: Dict[str, str]) -> str:
        """
          Get the name of the model based on id
          :param id2name_dict: dictionary with id as key and name as value
          """
        self.model_name = id2name_dict[self.model_id]


class InfernalAlignment(Annotation):
    """
    Single alignment with InfeRNAl
    covariance model
    """

    parser_dict = {'cmsearch': {'seq_id': 0,
                                'model_name': 2,
                                'model_id': 3,
                                'seq_start': 7,
                                'seq_end': 8,
                                'seq_strand': 9,
                                'evalue': 15,
                                'score': 14},
                   'cmscan': {'model_name': 0,
                              'model_id': 1,
                              'seq_id': 2,
                              'seq_start': 7,
                              'seq_end': 8,
                              'seq_strand': 9,
                              'evalue': 15,
                              'score': 14}}

    @classmethod
    def from_cmsearch_line(cls, line):
        return cls._from_tblout_line('cmsearch', line)

    @classmethod
    def from_cmscan_line(cls, line):
        return cls._from_tblout_line('cmscan', line)


class HmmerAlignment(Annotation):
    """
    Object representing a single
    HMMer3 hit
    a domain or protein family
    aligned to an amino acid sequence
    """

    parser_dict = {'hmmsearch': {'seq_id': 0,
                                 'model_name': 3,  # todo add actual names to db (2)
                                 'model_id': 3,
                                 'seq_start': 17,
                                 'seq_end': 18,
                                 'evalue': 6,
                                 'ivalue': 12,
                                 'score': 13},
                   'hmmscan': {'model_name': 0,  # todo add actual names to db (1)
                               'model_id': 0,
                               'seq_id': 3,
                               'seq_start': 17,
                               'seq_end': 18,
                               'evalue': 6,
                               'ivalue': 12,
                               'score': 13}}

    @classmethod
    def from_hmmscan_line(cls,
                          line: str) -> 'HmmerAlignment':
        return cls._from_tblout_line('hmmscan', line, unstranded=True)

    @classmethod
    def from_hmmsearch_line(cls,
                            line: str) -> 'HmmerAlignment':
        return cls._from_tblout_line('hmmsearch', line, unstranded=True)

    def to_dna(self,
               dna_len) -> 'DnaAlignment':
        dna_id, frame = self.seq_id.split('___')
        frame = int(frame)
        if 3 >= frame >= 1:
            dna_start, dna_end = (self.start - 1) * 3 + frame, self.end * 3 + frame - 1
        elif -3 <= frame <= -1:
            dna_start, dna_end = dna_len - self.end * 3 + frame + 1, dna_len - (self.start * 3) + frame + 4
        else:
            raise ValueError(f'Invalid frame: {frame}')
        assert 1 <= dna_start < dna_end <= dna_len
        strand = '+' if frame > 0 else '-'
        return DnaAlignment(seq_id=dna_id,
                            strand=strand,
                            model_id=self.model_id,
                            model_name=self.model_name,
                            start=dna_start,
                            end=dna_end,
                            score=self.score,
                            evalue=self.evalue,
                            method=self.method)


class MmseqsAlignment(Annotation):
    """
    Object representing a single
    MMseqs2 hit an alignment to a
    domain or protein family
    """

    parser_dict = {'mmseq': {'seq_id': 1,
                             'model_id': 0,
                             'model_name': 0,
                             'start': 8,
                             'end': 9,
                             'evalue': 4,
                             'score': 2}}

    @classmethod
    def from_mmseqs_line(cls,
                         line: str):
        return cls._from_tblout_line('mmseqs', line, unstranded=True)

    def to_dna(self,
               frame,
               dna_id,
               dna_len):
        if 3 >= frame >= 1:
            dna_start, dna_end = (self.start - 1) * 3 + frame, self.end * 3 + frame - 1
        elif -3 <= frame <= -1:
            dna_start, dna_end = dna_len - self.end * 3 + frame + 1, dna_len - (self.start * 3) + frame + 4
        else:
            raise ValueError(f'Invalid frame: {frame}')
        assert 1 <= dna_start < dna_end <= dna_len
        strand = '+' if frame > 0 else '-'
        return DnaAlignment(seq_id=dna_id,
                            strand=strand,
                            model_id=self.model_id,
                            model_name=self.model_id,
                            start=dna_start,
                            end=dna_end,
                            score=self.score,
                            evalue=self.evalue,
                            method=self.method)


class DnaAlignment(Annotation): # TODO DnaAlignment
    """
    Object representing a single
    MMseqs2 hit
    a domain or protein family
    re-anchored to a DNA sequence
    """


class Intron(Annotation):
    """
    Intron - intervening sequence in the gene structure
    """


class Exon(Annotation):
    """
    Exon - coding sequence in the gene structure
    """


PutativeExon = Union[Exon, DnaAlignment]


class Gene(Annotation):
    """
    Gene - a sequence of exons and introns
    """

    def __init__(self,
                 seq_id: str,
                 strand: str,
                 model_id: str,
                 model_name: str,
                 start: int,
                 end: int,
                 score: float,
                 evalue: float,
                 method: str):
        Annotation.__init__(self,
                            seq_id=seq_id,
                            strand=strand,
                            model_id=model_id,
                            model_name=model_name,
                            start=start,
                            end=end,
                            score=score,
                            evalue=evalue,
                            method=method)

    def host_exon(self,
                  exon: Annotation):
        """
        Checks exon for compatibility,
        add it to this gene,
        update gene coordinates.
        :param exon: exon to be added
        :return: an added exon
        """
        exon = Exon.retype(exon)
        assert self.seq_id == exon.seq_id, f'Cannot add exon from different sequence {self.seq_id} {exon.seq_id}'
        assert self.strand == exon.strand, f'Cannot add exon from different strand {self.strand} {exon.strand}'
        assert self.model_id == exon.model_id, f'Cannot add exon from different model {self.model_id} {exon.model_id}'
        assert self.method == exon.method, f'Cannot add exon from different method {self.method} {exon.method}'
        overlaps = [exon.overlaps(e) for e in self.children(Exon)]
        significant_overlaps = [o for o in overlaps if o > 0.1]
        try:
            assert not any(significant_overlaps), 'Exon overlaps with existing exons'
        except AssertionError:
            significant_overlaps = [f'{o:.2f}' for o in significant_overlaps]
            logger.warning(f'Suspicious overlap in {self.seq_id} {significant_overlaps}')
        if exon.start < self.start:
            self.start = exon.start
        if exon.end > self.end:
            self.end = exon.end
        self.score += exon.score

        self.host(exon)
        return exon

    @classmethod
    def from_exons(cls,
                   domains: List[PutativeExon],
                   min_intron: int = 20): # TODO ->
        """
        TODO add docstring
        :param domains:
        :param min_intron:
        :return:
        """
        xxx = [type(d) for d in domains] # TODO WTF?
        assert all([isinstance(d, PutativeExon) for d in
                    domains]), f'method can only use Exon or DnaAlignments annotations to create a gene {xxx}'
        sorted_domains = list(sorted(domains, key=lambda d: d.start))
        first_domain, last_domain = sorted_domains[0], sorted_domains[-1]

        gene = Gene(seq_id=first_domain.seq_id,
                    strand=first_domain.strand,
                    model_id=first_domain.model_id,
                    model_name=first_domain.model_name,
                    start=first_domain.start,
                    end=last_domain.end,
                    score=0,
                    evalue=min([d.evalue for d in sorted_domains]),
                    method=first_domain.method)

        for d in sorted_domains:
            gene.host_exon(d)

        previous = gene.start
        while sorted_domains:
            current_exon = sorted_domains.pop(0)
            if current_exon.start - previous >= min_intron:
                gene.host(Intron(seq_id=gene.seq_id,
                                 strand=gene.strand,
                                 model_id=gene.model_id,
                                 model_name=gene.model_name,
                                 start=previous + 1,
                                 end=current_exon.start - 1,
                                 score=0,
                                 evalue=1,
                                 method=gene.method))
            previous = current_exon.end

        return gene


class AnnotationTrack(list):

    def __init__(self,
                 seq_id: str,
                 annotations: List[Annotation] = ()):
        list.__init__(self, annotations)
        self.seq_id = seq_id
        self.method = 'user'  # TODO why the fuck is this here?

    def filter_score(self,
                     threshold: float,
                     recursive: bool = True):
        """
        Remove all annotations with score below the threshold
        :param threshold: minimum score to keep the annotation
        :param recursive: remove annotations nested in the deeper layers (e.g. exons within genes)
        :return:
        """
        above_threshold = [a for a in self if a.score >= threshold]
        if recursive:
            for a in above_threshold:
                a.nested = a.nested.filter_score(threshold, recursive)
        return AnnotationTrack(self.seq_id, above_threshold)

    def contexts(self,
                 flank: int) -> 'AnnotationTrack':
        """
        Extend all annotations from both sides
        by a specified number of flanking residues
        :param flank: number of residues to extend the annotations by
        :return: extended annotations
        """
        return AnnotationTrack(self.seq_id, [a.context(flank) for a in self])

    def merge_overlapping(self,
                          threshold: float = 0) -> 'AnnotationTrack':
        """
        Merge overlapping annotations
        :return: merged annotations
        """
        separate = AnnotationTrack(self.seq_id)
        tmp = self.copy()
        while tmp:
            current = tmp.pop()
            overlap_found = False
            for annotation in tmp:
                if current.overlaps(annotation) > threshold:
                    tmp.remove(annotation)
                    tmp.append(current.merge(annotation))
                    overlap_found = True
                    break
            if not overlap_found:
                separate.append(current)
        return separate

    def dock_overlapping(self,
                         annotations: 'AnnotationTrack',
                         overlap_threshold: float = 0.5,
                         culled: bool = True):
        """
        Include nested annotations that overlap with any of the annotations in the track
        :param annotations:
        :param overlap_threshold:
        :param culled:
        :return:
        """
        for a in self:
            a.dock_overlapping(annotations, overlap_threshold, culled)

    def cull(self,
             overlap_threshold: float = 0.5):
        """
        Choose only (locally) the best annotations for all sequences in a single file
        :param overlap_threshold: minimum fraction of overlap to report
        """

        ordered_annotations = sorted(self, key=lambda a: a.start)
        culled_annotations = []
        while ordered_annotations:
            new_annotation = ordered_annotations.pop(0)
            overlapping_annotations = []
            non_overlapping_annotations = []
            while ordered_annotations:
                old_annotation = ordered_annotations.pop(0)
                if new_annotation.overlaps(old_annotation) > overlap_threshold:
                    overlapping_annotations.append(old_annotation)
                else:
                    non_overlapping_annotations.append(old_annotation)
            if overlapping_annotations:
                main_annotation = Annotation.top_scoring(overlapping_annotations + [new_annotation])
                ordered_annotations = non_overlapping_annotations + [main_annotation]
                ordered_annotations.sort(key=lambda a: a.start)
            else:
                culled_annotations.append(new_annotation)
                ordered_annotations = non_overlapping_annotations
        return AnnotationTrack(self.seq_id, culled_annotations)

    def sort_annotations(self,
                         by: str = 'start'):
        """
        Define sorting order for annotations on the tracks
        :param by: attribute to sort by (e.g. start, end, score)
        """
        self.sort(key=lambda a: getattr(a, by))

    def filter(self, subtype: Type[Annotation]):
        """
        Return only annotations of a given type
        :param subtype: the type of annotation to keep
        """
        return AnnotationTrack(self.seq_id, [a for a in self if isinstance(a, subtype)])

    def children(self, subtype: Type[Annotation]) -> 'AnnotationTrack':
        """
        Get all annotations of a given type nested within annotations of the track
        :param subtype: the type of nested annotations
        :return: list of nested annotations of a given type
        """
        typed_children = []
        for annotation in self:
            typed_children.extend(annotation.children(subtype))
        return AnnotationTrack(self.seq_id, typed_children)

    def to_dict(self,
                keys: str) -> Dict[str, Dict[str, Annotation]]:
        """
        Convert AnnotationBase to a dictionary
        using values of specified attribute as keys
        :param keys: attribute to use as a key
        :return: dictionary with annotations
        """
        result_dict = defaultdict(list)
        for a in self:
            result_dict[getattr(a, keys)].append(a)
        return dict(result_dict)

    def annotation_sequences(self,
                             seq: Seq.Seq) -> List[Seq.Seq]:
        """
        Get the sequence of the annotation
        :param seq: sequence of the annotated contig
        :return: sequence of the annotation
        """
        return [a.get_sequence(seq) for a in self]

    def get_model_names(self,
                        id2name_dict: Dict[str, str]):
        """
        Get the name of the model based on id
        :param id2name_dict: dictionary with id as key and name as value
        """
        [a.get_model_name(id2name_dict) for a in self]


class AnnotationBase(dict):

    def __init__(self,
                 tracks: Dict[str, AnnotationTrack] = {}):
        dict.__init__(self, tracks)

    @classmethod
    def from_infernal(cls,
                      infernal: Path,
                      program: str = 'cmscan') -> 'AnnotationBase':
        """
        Read infernal flatfile
        :param infernal: input cmscan/cmsearch file
        :param program: infernal program used to generate the file
        :return: {'contig_id': [Rna1, RNA2, (...)], (...)}
        """
        logger.info(f'loading {infernal.as_posix()}')
        with infernal.open() as handle:
            contig_to_genes = defaultdict(list)
            for line in handle:
                if line.startswith('#') or not line.strip():
                    pass
                else:
                    if program == 'cmsearch':
                        gene = InfernalAlignment.from_cmsearch_line(line)
                    elif program == 'cmscan':
                        gene = InfernalAlignment.from_cmscan_line(line)
                    contig_to_genes[gene.seq_id].append(gene)

        return cls({k: AnnotationTrack(k, v) for k, v in contig_to_genes.items()})

    @classmethod
    def from_hmmer(cls,
                   domtblout: Path,
                   program: str = 'hmmscan') -> 'AnnotationBase':
        """
        Read hmmer flatfile
        :param domtblout: input hmmscan file
        :param program: HMMer3 program used to generate the file
        :return: {'protein_id': [Aln1, Aln2, (...)], (...)}
        """
        logger.info(f'loading {domtblout.as_posix()}')
        with domtblout.open() as handle:
            contig_to_genes = defaultdict(list)
            for line in handle:
                if line.startswith('#') or not line.strip():
                    pass
                else:
                    if program == 'hmmscan':
                        gene = HmmerAlignment.from_hmmscan_line(line)
                    elif program == 'hmmsearch':
                        gene = HmmerAlignment.from_hmmsearch_line(line)
                    contig_to_genes[gene.seq_id].append(gene)
        return cls({k: AnnotationTrack(k, v) for k, v in contig_to_genes.items()})

    @classmethod
    def from_gff(cls,
                 gff: Path,
                 subtype: Type[Annotation]) -> 'AnnotationBase':
        """
        Read general feature format(GFF) file
        :param gff: input gff file
        :param subtype: type of annotation to create
        :return: {'contig_id': [Annotation1, Annotation12, (...)], (...)}
        """
        logger.info(f'loading {gff.as_posix()}')
        with gff.open() as handle:
            contig_to_genes = defaultdict(list)
            for line in handle:
                if line.startswith('#') or not line.strip():
                    pass
                else:
                    annot = subtype.from_gff_line(line)
                    contig_to_genes[annot.seq_id].append(annot)
        return cls({k: AnnotationTrack(k, v) for k, v in contig_to_genes.items()})

    def save_gff(self,
                 write_gff: Path):
        """
        Dump gff file with coordinates of annotations from contigs on the drive
        :param write_gff: path where file should be written
        """
        with write_gff.open('w') as handle:
            lines = []
            for contig, annotations in self.items():
                for a in annotations:
                    lines.extend(a.gff())
            handle.write('\n'.join(lines))

    def save_genbank(self,
                     write_gbk: Path,
                     fasta: Path):
        """
        Dump genbank file with contig (DNA) sequences and gene annotations
        :param write_gbk: path where file should be written
        :param fasta: path to the fasta file with annotated sequences
        """
        sequences = SeqIO.index(fasta.as_posix(), 'fasta')
        output = []
        for contig, genes in self.items():
            seq = sequences[contig]
            seq.annotations['molecule_type'] = 'DNA'
            seq.features = [Feature.SeqFeature(Feature.FeatureLocation(g.start,
                                                                       g.end,
                                                                       strand=int(f'{g.strand}1')),
                                               type="CDS",
                                               qualifiers={'protein_id': g.id,
                                                           'translation': g.translation}) for g in genes]
            output.append(seq)
        SeqIO.write(output, write_gbk, 'genbank')

    def annotate(self,
                 annotation: Annotation):
        """
        Add annotation to the annotation base
        :param annotation: single annotation to add
        """
        if annotation.seq_id not in self:
            self[annotation.seq_id] = AnnotationTrack(annotation.seq_id)
        self[annotation.seq_id].append(annotation)

    def sort_annotations(self,
                         by: str = 'start'):
        """
        Define sorting order for annotations on the tracks
        :param by: attribute to sort by (e.g. start, end, score)
        """
        for track in self.values():
            track.sort_annotations(key=lambda a: getattr(a, by))

    def include(self,
                other: 'AnnotationBase',
                sort: bool = True):
        """
        Merge two annotation bases together add annotations and tracks

        :param other:
        :param sort:
        """
        # add new tracks or append annotations to existing tracks
        for seq_id, track in other.items():
            if seq_id not in self:
                self[seq_id] = track
            else:
                self[seq_id].extend(track)
        if sort:
            self.sort_annotations()

    def contexts(self, flank: int) -> 'AnnotationBase':
        """
        Extend all annotations from both sides
        by a specified number of flanking residues
        :param flank: number of residues to extend the annotations by
        :return: extended annotations
        """
        return AnnotationBase({seq_id: track.contexts(flank) for seq_id, track in self.items()})

    def merge_overlapping(self, threshold: float = 0) -> 'AnnotationBase':
        """
        Merge overlapping annotations
        :return: merged annotations
        """
        return AnnotationBase({seq_id: track.merge_overlapping(threshold) for seq_id, track in self.items()})

    def dock_overlapping(self,
                         annotations: 'AnnotationBase',
                         overlap_threshold: float = 0.5,
                         culled: bool = True):
        """
        Include nested annotations that overlap with any of the annotations in the track
        :param annotations:
        :param overlap_threshold:
        :param culled:
        :return:
        """
        for seq_id, track in self.items():
            if seq_id in annotations:
                track.dock_overlapping(annotations[seq_id], overlap_threshold, culled)

    def to_dict(self,
                keys: str) -> Dict[str, Dict[str, Annotation]]:
        """
        Convert AnnotationBase to a dictionary
        using values of specified attribute as keys
        :param keys: attribute to use as a key
        :return: dictionary with annotations
        """
        return {seq_id: track.to_dict(keys) for seq_id, track in self.items()}

    def children(self, subtype: Type[Annotation]) -> 'AnnotationBase':
        """
        Get all nested annotations of a given type
        :param subtype: the type of nested annotations
        :return: list of nested annotations of a given type
        """
        return AnnotationBase({seq_id: track.children(subtype) for seq_id, track in self.items()})

    def filter_score(self,
                     threshold: float,
                     recursive: bool = True):
        """
        Remove all annotations with score below the threshold
        :param threshold: score threshold below which annotations are removed
        :param recursive: should the method be applied to nested annotations?
        :return:
        """
        return AnnotationBase(
            {seq_id: track.filter_score(threshold=threshold, recursive=recursive) for seq_id, track in self.items()})

    def cull(self):
        """
        Choose only (locally) the best annotations for all sequences in a single file
        """
        return AnnotationBase({seq_id: track.cull() for seq_id, track in self.items()})

    def annotation_sequences(self,
                             input_fasta: Path,
                             output_fasta: Path = None) -> Dict[str, List[Seq.Seq]]:
        """
        Get the sequence of the annotation
        :param input_fasta: path to the fasta file with annotated sequences
        :return: extracted sequences of all the annotations in the annotation base
        :param output_fasta: path to the fasta file with annotation sequences (optional)
        """
        input_sequences = SeqIO.index(input_fasta.as_posix(), 'fasta')
        output_sequences = []
        for contig, genes in self.items():
            seq = input_sequences[contig]
            output_sequences.extend(genes.annotation_sequences(seq))
        if output_fasta is not None:
            SeqIO.write(output_sequences, output_fasta.as_posix(), 'fasta')
        return input_sequences

    def get_model_names(self,
                        id2name_dict: Dict[str, str]):
        """
        Get the name of the model based on id
        :param id2name_dict: dictionary with id as key and name as value
        """
        [track.get_model_names(id2name_dict) for track in self.values()]

    def with_sequences(self,
                       fasta: Path):
        """
        Iterate over sequences in a fasta file
        corresponding to the annotation tracks
        to get (seq_id, sequence, annotation)
        :param fasta: path to a fasta file
        :return: generator that yields (seq_id, sequence) pairs
        """
        sequences = parse_fasta(fasta)
        used_seq_ids = set()
        for seq_id, seq in sequences:
            used_seq_ids.add(seq_id)
            yield seq_id, seq, self[seq_id]

        # check any sequence is missing
        not_in_self = used_seq_ids - set(self.keys())
        not_in_fasta = set(self.keys()) - used_seq_ids
        assert not not_in_self, f'No annotations for sequences: {not_in_self}'
        assert not not_in_fasta, f'No sequences for annotations: {not_in_fasta}'