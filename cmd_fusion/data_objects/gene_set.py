#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gene annotations.

@author: T.D. Medina
"""

from functools import total_ordering, cached_property
import gzip
import pickle
import re

from tqdm import tqdm

from cmd_fusion.utilities import _GENOME_ORDER, overlap, iterize

# %% Annotation Types
class GeneAnnotation:
    """Record of a single gene annotation."""

    # pylint: disable=too-many-instance-attributes

    def __init__(self, gene_id, seqname, start, end, strand,
                 transcript_id=None, gene_name=None, feature=None,
                 score=".", frame=".", exon_id=None, exon_number=None,
                 source=".", **kwargs):
        self.gene_id = gene_id
        self._gene_name = gene_name
        self.transcript_id = transcript_id

        self.feature = feature

        self.exon_id = exon_id
        self.exon_number = int(exon_number) if exon_number is not None else exon_number

        self.seqname = seqname
        self.start = int(start)
        self.end = int(end)

        self.strand = strand
        self.frame = frame
        self.score = score

        self.source = source
        self.other_annotations = kwargs

        self.canonical = None
        if "tag" in self.other_annotations:
            self.canonical = "Ensembl_canonical" in self.other_annotations["tag"]

        # self._hash = hash(str(self.__dict__))

    def __repr__(self):
        """Get official string representation."""
        attrs = [f"{x}={y}" for x, y in self.__dict__.items()
                 if y != "." and y is not None and not x.startswith("_")]
        string = f"{type(self).__name__}(" + ", ".join(attrs) + ")"
        return string

    def __hash__(self):
        return hash(f"{self.gene_id},{self.transcript_id},{self.exon_id},"
                    f"{self.feature},{self.range},{self.strand}")

    def _is_valid_operand(self, other):
        return isinstance(other, type(self))

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        # return self._hash == other._hash
        return str(self.__dict__) == str(other.__dict__)

    @property
    def gene_name(self):
        return self._gene_name if self._gene_name else self.gene_id

    @property
    def range(self):
        return range(self.start, self.end+1)

    def __len__(self):
        return len(self.range)

    @property
    def length(self):
        return len(self.range)

    def is_transcript(self):
        """Test if gene feature is transcript."""
        return self.feature == "transcript"


class Gene(GeneAnnotation):
    """Ensembl Gene object with unique gene ID."""

    def __init__(self, gene_id, seqname, start, end, strand,
                 transcripts=None, **kwargs):
        super().__init__(gene_id, seqname, start, end, strand, **kwargs)
        self.transcripts = transcripts if transcripts is not None else []
        # self._hash = self._make_hash()

    def __str__(self):
        return self.gene_name

    def __repr__(self):
        """Official string representation."""
        locus = f"{self.seqname}:{self.start}-{self.end}"
        string = (f"Gene(gene_id={self.gene_id}, gene_name={self._gene_name}, "
                  f"source={self.source}, locus={locus}, "
                  f"strand={self.strand}, transcripts={len(self.transcripts)})")
        return string

    def __hash__(self):
        return hash(self.gene_id)

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.gene_id == other.gene_id

    def get_exons_by_locus(self, locus, as_ids=False, as_numbers=False,
                           include_transcript_ids=True, canonical_only=True):
        if canonical_only:
            canonical = self.get_canonical_transcript()
            exons = canonical.get_exon_by_locus(locus, as_ids, as_numbers)
            if include_transcript_ids:
                exons = {canonical.transcript_id: exons}
            return exons
        exons = {transcript.transcript_id: results for transcript in self.transcripts
                 if (results := transcript.get_exon_by_locus(locus, as_ids, as_numbers))}
        if not include_transcript_ids and as_ids:
            exons = {exon for exon_list in exons.values() for exon in exon_list}
        return exons

    def get_canonical_transcript(self):
        canonical = [trans for trans in self.transcripts if trans.canonical]
        if len(canonical) == 1:
            return canonical[0]
        if not canonical:
            return None
        if len(canonical) > 1:
            raise Warning(f"Multiple transcripts marked canonical for gene {self}:\n"
                          f"{canonical}")


class Transcript(GeneAnnotation):
    """Gene sub-annotation of a single transcript.

    Designed to be nested inside a parent-level gene annotation.
    """
    def __init__(self, gene_id, transcript_id, seqname, start, end, strand,
                 exons=None, annotations=None, **kwargs):
        super().__init__(gene_id, seqname, start, end, strand,
                         transcript_id=transcript_id, **kwargs)
        self.exons = exons if exons is not None else []
        self._annotations = annotations if annotations is not None else []

    def __hash__(self):
        return hash(self.transcript_id)

    def __repr__(self):
        locus = f"{self.seqname}:{self.start}-{self.end}"
        string = (f"Transcript(gene_id={self.gene_id}, gene_name={self._gene_name}, "
                  f"source={self.source}, locus={locus}, "
                  f"strand={self.strand}, exons={len(self.exons)})")
        return string

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.transcript_id == other.transcript_id

    def get_exon_by_locus(self, locus, as_ids=False, as_numbers=False):
        exons = [exon for exon in self.exons if locus in exon.range]
        if as_ids:
            exons = [exon.exon_id for exon in exons]
        elif as_numbers:
            exons = [exon.exon_number for exon in exons]
        return exons


class Exon(GeneAnnotation):
    """Gene annotation of a single exon of a transcript."""
    def __init__(self, gene_id, transcript_id, exon_id, exon_number,
                 seqname, start, end, strand, **kwargs):
        super().__init__(gene_id, seqname, start, end, strand,
                         transcript_id=transcript_id,
                         exon_id=exon_id, exon_number=exon_number, **kwargs)

    def __str__(self):
        return str(self.exon_number)

    def __repr__(self):
        locus = f"{self.seqname}:{self.start}-{self.end}"
        string = (f"Exon(gene_id={self.gene_id}, gene_name={self._gene_name}, "
                  f"exon_id={self.exon_id}, "
                  f"source={self.source}, locus={locus}, "
                  f"strand={self.strand}, exon_number={self.exon_number})")
        return string

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.exon_id == other.exon_id

    def __hash__(self):
        return hash(f"{self.transcript_id},{self.exon_id}")


# %% SubAnnotations.
class SubAnnotation(GeneAnnotation):
    def __init__(self, gene_id, transcript_id, feature, seqname, start, end, strand,
                 **kwargs):
        super().__init__(gene_id, seqname, start, end, strand,
                         transcript_id=transcript_id, feature=feature, **kwargs)


class UTR3(SubAnnotation):
    pass


class UTR5(SubAnnotation):
    pass


class StartCodon(SubAnnotation):
    pass


class StopCodon(SubAnnotation):
    pass


class CodingSequence(SubAnnotation):
    pass


# %% GeneSet
class GeneSet:
    """Database of gene annotations."""

    def __init__(self, path=None, include=None, exclude=None, _geneset=None):
        if _geneset is not None:
            self.path = _geneset.path
            self.genes = _geneset.genes
            self.chromosomes = _geneset.chromosomes
            return
        self.path = path
        self.genes = {}
        self.chromosomes = {}

        if self.path:
            self.genes = self.make_gene_set(path, include, exclude)
            self.organize_chromosomes()

        self.size = len(self.genes)

    def __len__(self):
        return self.size

    def make_gene_set(self, path, include, exclude):
        """Construct Genes and GeneSet objects from data file."""
        print("Reading GTF file...")
        genes = self.read_annotations(path, include, exclude)
        print("Creating annotation entries...")
        genes = self.make_annotation_objects(genes)
        print("Organizing GeneSet...")
        genes = self.make_genes(genes)
        return genes

    @staticmethod
    def read_annotations(file, include, exclude):
        """Read gene annotations from file."""
        include, exclude = iterize(include), iterize(exclude)
        with gzip.open(file) as infile:
            data = []
            for line in tqdm(infile):
                line = line.decode().rstrip(";\n").split("\t")
                if line[0].replace("chr", "", 1) not in _GENOME_ORDER:
                    continue
                if include is not None and line[2] not in include:
                    continue
                if exclude is not None and line[2] in exclude:
                    continue
                line[0] = line[0].lstrip("chr")
                line[3] = int(line[3])
                line[4] = int(line[4])
                ids = dict()
                for field in line[-1].split(";"):
                    key, val = field.strip().replace('"', "").split(" ", 1)
                    if key in ids:
                        if not isinstance(ids[key], list):
                            ids[key] = [ids[key]]
                        ids[key].append(val)
                        continue
                    ids[key] = val
                line[-1] = ids
                data.append(line)
        return data


    @staticmethod
    def make_annotation_objects(data):
        """Construct gene sub-annotation objects."""
        fields = ["seqname", "source", "feature", "start", "end",
                  "score", "strand", "frame"]
        classes = {"gene": Gene, "transcript": Transcript, "exon": Exon,
                   "three_prime_utr": UTR3, "five_prime_utr": UTR5,
                   "start_codon": StartCodon, "stop_codon": StopCodon,
                   "CDS": CodingSequence, "Selenocysteine": None}
        objects = set()
        while data:
            line = data.pop()
            line_dict = dict(zip(fields, line[:-1]))
            if classes[line_dict["feature"]] is None:
                continue
            line_dict.update(line[-1])
            objects.add(classes[line_dict["feature"]](**line_dict))
        return objects

    @staticmethod
    def make_genes(data):
        """Make Gene objects from sub-annotations."""
        genes = {gene for gene in data if isinstance(gene, Gene)}
        data -= genes
        genes = {gene.gene_id: gene for gene in genes}
        transcripts = {trans for trans in data if isinstance(trans, Transcript)}
        data -= transcripts
        transcripts = {trans.transcript_id: trans for trans in transcripts}

        while data:
            entry = data.pop()
            if isinstance(entry, Exon):
                transcripts[entry.transcript_id].exons.append(entry)
            else:
                transcripts[entry.transcript_id]._annotations.append(entry)
        transcripts = list(transcripts.values())

        for transcript in transcripts:
            transcript.exons.sort(key=lambda exon: exon.exon_number)

        while transcripts:
            trans = transcripts.pop()
            genes[trans.gene_id].transcripts.append(trans)
        for gene in genes.values():
            gene.transcripts.sort(key=lambda x: x.end)
            gene.transcripts.sort(key=lambda x: x.start)
        return genes

    def organize_chromosomes(self):
        """Group genes by chromosome in a chromosome dictionary."""
        chrom_dict = {}
        for gene in self.genes.values():
            if gene.seqname not in chrom_dict:
                chrom_dict[gene.seqname] = {}
            chrom_dict[gene.seqname][gene.gene_id] = gene
        self.chromosomes = chrom_dict

    def get_locus(self, seqname, start, stop=None):
        """Return all genes that intersect a base or range."""
        if seqname not in self.chromosomes:
            return []
        if stop is None:
            stop = start
        query = range(start, stop + 1)
        results = []
        for gene in self.chromosomes[seqname].values():
            if overlap(query, range(gene.start, gene.end+1)):
                results.append(gene)
        return results

    def get_all_gene_ids(self):
        """Get gene IDs from all genes."""
        return set(self.genes.keys())

    def lookup_ensembl_id(self, gene_id):
        if gene_id in self.genes:
            return self.genes[gene_id]

    def lookup_gene_name(self, gene_name, seqname=None, exact=False):
        lookup = self.genes if seqname is None else self.chromosomes[seqname]
        if exact:
            for gene in lookup.values():
                if gene._gene_name and gene._gene_name.upper() == gene_name.upper():
                    return gene
            return
        results = [gene for gene in lookup.values()
                   if gene._gene_name and re.findall(gene_name, gene._gene_name, re.I)]
        return results

    @staticmethod
    def read_pickle(pkl):
        with open(pkl, "rb") as inpkl:
            geneset = pickle.load(inpkl)
        return geneset


@total_ordering
class Locus:
    def __init__(self, contig, position, reference):
        self.contig = contig
        self.position = int(position)
        self.reference = reference.lower()

    def _make_hash(self):
        return

    def __hash__(self):
        return self._hash

    @cached_property
    def _hash(self):
        return hash(f"{self.contig}:{self.position}:{self.reference}")

    def __repr__(self):
        string = f"Locus('{self.contig}', {self.position}, '{self.reference}')"
        return string

    def __str__(self):
        string = f"{self.contig}:{self.position}"
        return string

    def _is_valid_operand(self, other):
        if isinstance(other, type(self)):
            if self.reference == other.reference:
                return True
            raise Warning(f"Attempted to compare loci from difference references.")
        return False

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        if self.contig == other.contig:
            return self.position < other.position
        return _GENOME_ORDER[self.contig] < _GENOME_ORDER[other.contig]

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.contig == other.contig and self.position == other.position

    @staticmethod
    def convert_to_locus(locus_string, reference):
        """Convert e.g. chrX:123456 to Locus object."""
        return Locus(*locus_string.strip().replace("chr", "").split(":"),
                     reference=reference)
