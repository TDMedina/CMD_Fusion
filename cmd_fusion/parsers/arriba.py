
import re

import pandas as pd
from pandas import MultiIndex

from cmd_fusion.parsers.common import set_standard_index, _read_keep_columns, make_empty_table
from cmd_fusion.utilities import iterize
from cmd_fusion.data_objects.gene_set import GeneSet, Locus


class ArribaFusions:
    def __init__(self, filepath, geneset: GeneSet, confidence=None):
        if filepath is None:
            self.table = make_empty_table()
            return
        self.filepath = filepath
        self.table = self.read_csv(self.filepath, geneset, confidence)

    @classmethod
    def read_csv(cls, filepath, geneset, confidence=None):
        table = pd.read_csv(filepath, sep="\t", na_values=".")
        table.columns = cls._parse_headers(table.columns)
        if confidence is not None:
            confidence = iterize(confidence)
            table = table.loc[table.confidence.isin(confidence)]
        for i in "ab":
            dex = lambda x: ("General", x, i)

            table[dex("Breakpoints")] = [Locus.convert_to_locus(locus, "hg38")
                                         for locus in table.Arriba.breakpoint[i]]

            table[dex("Genes")] = [geneset.genes[gene_id]
                                   for gene_id in table.Arriba.gene_id[i]]

            get_exon = lambda x, y: x.get_exons_by_locus(y, include_transcript_ids=False)
            genes_bps = zip(table.General.Genes[i], table.General.Breakpoints[i])
            table[dex("Exons")] = [exon[0] if (exon := get_exon(gene, bp.position)) else None
                                   for gene, bp in genes_bps]

        cols = [("General", field, i) for field in ["Genes", "Exons", "Breakpoints"]
                for i in "ab"]
        cols += _read_keep_columns("columns/arriba_columns.csv")
        table = table[cols]
        table = set_standard_index(table)
        return table

    @staticmethod
    def _parse_headers(col_headers):
        col_headers = [re.findall(r"(\S+?)(\d*)$", x.lstrip("#").split("(")[0])[0]
                       for x in col_headers]
        head_dict = dict(zip("12", "ab"))
        for i, header in enumerate(col_headers):
            if not header[1] or header[1] not in "12":
                continue
            header = list(header)
            header[1] = head_dict[header[1]]
            col_headers[i] = tuple(header)
        col_headers = MultiIndex.from_tuples([("Arriba", *headers) for headers in col_headers])
        return col_headers
