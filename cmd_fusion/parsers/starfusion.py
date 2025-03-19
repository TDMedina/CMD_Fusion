
import re

import pandas as pd
from pandas import MultiIndex

from cmd_fusion.utilities import iterize
from cmd_fusion.data_objects.gene_set import GeneSet, Locus
from cmd_fusion.parsers.common import set_standard_index, _read_keep_columns, make_empty_table


class StarFusions:
    def __init__(self, filepath, geneset: GeneSet):
        if filepath is None:
            self.table = make_empty_table()
            return
        self.filepath = filepath
        self.table = self.read_csv(self.filepath, geneset)

    @classmethod
    def read_csv(cls, filepath, geneset):
        table = pd.read_csv(filepath, sep="\t", na_values=".")
        table.columns = cls._parse_headers(table.columns)
        for i in "ab":
            dex = lambda x: ("General", x, i)

            table[dex("Breakpoints")] = [Locus.convert_to_locus(locus[:-2], "hg38")
                                         for locus in table.STAR.Breakpoint[i]]

            table[dex("Genes")] = [geneset.genes[gene_id.split("^")[1]]
                                   for gene_id in table.STAR.Gene[i]]

            get_exon = lambda x, y: x.get_exons_by_locus(y, include_transcript_ids=False)
            genes_bps = zip(table.General.Genes[i], table.General.Breakpoints[i])
            table[dex("Exons")] = [exon[0] if (exon := get_exon(gene, bp.position)) else None
                                   for gene, bp in genes_bps]
        cols = [("General", field, i) for field in ["Genes", "Exons", "Breakpoints"]
                for i in "ab"]
        cols += _read_keep_columns("columns/star_columns.csv")
        table = table[cols]
        table = set_standard_index(table)
        return table

    @staticmethod
    def _parse_headers(col_headers):
        head_dict = dict(zip(["left", "right"], "ab"))
        pattern = r"_?(left|right)_?"
        headers = []
        for col in col_headers:
            col = col.replace("#", "")
            matched = re.findall(pattern, col, re.I)
            matched = head_dict[matched[0].strip("_").lower()] if matched else ""
            replaced = re.sub(pattern, "", col, flags=re.I)
            header = (replaced, matched)
            headers.append(header)
        col_headers = MultiIndex.from_tuples([("STAR", *header) for header in headers])
        return col_headers
