
from importlib.resources import files

import pandas as pd
from pandas import MultiIndex
import yaml

from cmd_fusion.data_objects.gene_set import GeneSet, Locus
from cmd_fusion.parsers.common import set_standard_index, _read_keep_columns, make_empty_table
from cmd_fusion.utilities import FORWARD_LIFTOVER


class ThermoFusions:
    def __init__(self, filepath, geneset: GeneSet):
        if filepath is None:
            self.header = None
            self.table = make_empty_table()
            return
        self.header, self.table = self.read_csv(filepath)
        self.old_table, self.table = self.standardize_table(geneset)

    @classmethod
    def read_csv(cls, filepath):
        with open(filepath) as infile:
            header = []
            while (header_line := infile.readline()).startswith("#"):
                header.append(header_line)
        data = pd.read_csv(filepath, comment="#", sep="\t")
        header_data = cls._read_header_data()
        data = data.rename(mapper=header_data, axis=1)
        data = data[header_data.values()]
        data.columns = MultiIndex.from_tuples(data.columns)
        data[("Location", "Locus")] = [
            [convert_hg19_to_hg38(Locus.convert_to_locus(locus, "hg19"))
             for locus in loci.split("-")]
            for loci in data.Location.Locus
            ]

        cosmic = ("Annotation","COSMIC_NCBI")
        sanger = "https://cancer.sanger.ac.uk/cosmic/fusion/summary?id="
        hyper = lambda x: f"<a target='_blank' href='{sanger}{x.replace('COSF', '')}'>{x}</a>"
        data[("Annotation", "COSMIC")] = [hyper(x) if "COSF" in str(x) else x for x in data[cosmic]]

        data = data.loc[(data.Mutation.Type == "FUSION")
                        | (data.Mutation.Type == "RNAExonVariant")]
        data = data.reset_index(drop=True)
        return header, data

    @staticmethod
    def _read_header_data():
        header_file = files("cmd_fusion.parsers").joinpath("columns/thermo_header_data.yaml")
        with open(str(header_file)) as infile:
            header_data = yaml.safe_load(infile)
        header_data = {x: tuple(y) for x, y in header_data.items()}
        return header_data

    # @staticmethod
    # def _read_keep_columns():
    #     keep_col_file = files("cmd_fusion.parsers").joinpath("thermo_columns.yaml")
    #     with open(str(keep_col_file)) as infile:
    #           keep_cols = yaml.safe_load(infile)
    #     keep_cols = {tuple(k.split(",")): tuple(v.split(","))
    #                  for k, v in keep_cols.items()}
    #     return keep_cols

    def standardize_table(self, geneset: GeneSet):
        new_table = pd.DataFrame()

        breakpoints = zip(*self.table.Location.Locus)
        bp_cols = ("General", "Breakpoints", "a"), ("General", "Breakpoints", "b")
        new_table[bp_cols[0]],  new_table[bp_cols[1]] = breakpoints

        genes = [[gene.strip().split("(")[0] for gene in entry.split("-")]
                 for entry in self.table.Genes.Genes]
        gene_cols =  ("General", "Genes", "a"), ("General", "Genes", "b")
        new_table[gene_cols[0]], new_table[gene_cols[1]] = zip(*genes)
        for sub in "ab":
            new_table[("General", "Genes", sub)] = [
                geneset.lookup_gene_name(gene, bp.contig, True)
                for gene, bp in zip(new_table[("General", "Genes", sub)],
                                    new_table[("General", "Breakpoints", sub)])
                ]
            new_table[("General", "Exons", sub)] = [
                exon[0]
                if (exon := gene.get_exons_by_locus(bp.position, include_transcript_ids=False))
                else None
                for gene, bp in zip(new_table[("General", "Genes", sub)],
                                    new_table[("General", "Breakpoints", sub)])
                ]

        new_table = new_table[[("General", field, i)
                               for field in ["Genes", "Exons", "Breakpoints"]
                               for i in "ab"]]

        # keep_cols = self._read_keep_columns()
        keep_cols = _read_keep_columns("columns/thermo_columns.csv")
        for keep_col in keep_cols:
            new_table[keep_col] = self.table[keep_col[1:]]

        new_table.columns = pd.MultiIndex.from_tuples(new_table.columns)
        new_table = new_table[["General", "Thermo"]]

        new_table = set_standard_index(new_table)
        return self.table, new_table


def convert_hg19_to_hg38(locus: Locus):
    lifted = FORWARD_LIFTOVER[locus.contig][locus.position]
    if len(lifted) != 1:
        raise Warning(f"Non-singular results: {lifted}")
    lifted = Locus(lifted[0][0].replace("chr", ""), lifted[0][1], "hg38")
    return lifted
