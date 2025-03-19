
from collections import defaultdict
from importlib.resources import files

import pandas as pd
from pysam import VariantFile

from cmd_fusion.data_objects.gene_set import GeneSet, Locus
from cmd_fusion.utilities import FORWARD_LIFTOVER
from cmd_fusion.parsers.common import _read_keep_columns, set_standard_index, make_empty_table


class VCFFusions:
    def __init__(self, filepath, geneset: GeneSet):
        if filepath is None:
            self.table = make_empty_table()
        self.filepath = filepath
        self.table = self.read_vcf(self.filepath, geneset)

    @staticmethod
    def read_vcf(filepath, geneset):
        vcf = VariantFile(filepath)
        passing = [record for record in vcf.fetch()
                   if "PASS" in record.filter and record.info["SVTYPE"] == "Fusion"]
        fusions = defaultdict(list)
        for record in passing:
            fusions[record.id.split("_")[0]].append(record)

        table = pd.DataFrame()
        for i, j in enumerate("ab"):
            table[("General", "Breakpoints", j)] = [
                convert_hg19_to_hg38(Locus(fusion[i].contig.replace("chr", ""),
                                           fusion[i].pos,
                                           "hg19"))
                for fusion in fusions.values()
                ]
            table[("General", "Genes", j)] = [
                geneset.lookup_gene_name(fusion[i].info["GENE_NAME"], exact=True)
                for fusion in fusions.values()
                ]
            table["General", "Exons", j] = [
                exon[0]
                if (exon := gene.get_exons_by_locus(bp.position,
                                                    include_transcript_ids=False))
                else None
                for gene, bp in zip(table[("General", "Genes", j)],
                                    table[("General", "Breakpoints", j)])
                ]

        table = table[[("General", field, i) for field in ["Genes", "Exons", "Breakpoints"]
                       for i in "ab"]]
        keep_cols = _read_keep_columns("columns/vcf_columns.csv")

        for col in keep_cols:
            table[col] = [fusion[0].info[col[1]]
                          for fusion in fusions.values()]
        cosmic = ("Thermo", "Annotation", "COSMIC")
        if ("Thermo", "ANNOTATION", "") in keep_cols:
            sanger = "https://cancer.sanger.ac.uk/cosmic/fusion/summary?id="
            hyper = lambda x: f"<a target='_blank' href='{sanger}{x.replace('COSF', '')}'>{x}</a>"
            table[cosmic] = [hyper(x) if "COSF" in str(x) else x
                             for x in table[("Thermo", "ANNOTATION", "")]]

        table.columns = pd.MultiIndex.from_tuples(table.columns)
        table = set_standard_index(table)
        return table


def convert_hg19_to_hg38(locus: Locus):
    lifted = FORWARD_LIFTOVER[locus.contig][locus.position]
    if len(lifted) != 1:
        raise Warning(f"Non-singular results: {lifted}")
    lifted = Locus(lifted[0][0].replace("chr", ""), lifted[0][1], "hg38")
    return lifted