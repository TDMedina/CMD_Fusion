
from importlib.resources import files

import pandas as pd
from pandas import DataFrame


def set_standard_index(table: DataFrame):
    dex = []
    # info = table[["Genes", "Exons", "Breakpoints"]].T.groupby(level=1).agg(tuple).T
    info = table.General.T.groupby(level=1).agg(tuple).T
    for _, row in info.iterrows():
        rowdex = [(gene.gene_id,
                   (exon is None, 1 if exon is None else exon.exon_number),
                   bp)
                  for (gene, exon, bp) in row]
        for i in range(2, -1, -1):
            rowdex.sort(key=lambda x: x[i])
        rowdex = [(gid, "None" if exon[0] else str(exon[1]), str(bp))
                  for (gid, exon, bp) in rowdex]
        rowdex = tuple("::".join(pair) for pair in zip(*rowdex))
        dex.append(rowdex)
    dex = pd.MultiIndex.from_tuples(dex, names=["genes", "exons", "breakpoints"])
    return table.set_index(dex)


def merge_results(first, second):
    merged = pd.merge(first, second, how="outer", left_index=True, right_index=True,
                      indicator=True)

    for col in ["Genes", "Exons", "Breakpoints"]:
        for i in "ab":
            dex = lambda x: (f"General{x}", col, i)
            merged[dex("")] = merged[dex("_x")].where(merged[dex("_x")].notna(),
                                                      merged[dex("_y")])
            merged.drop([dex("_x"),dex("_y")], axis=1, inplace=True)
    return merged


def merge_all_results(thermo, arriba, star):
    merge_dict = {"both": ["Thermo", "Arriba"],
                  "left_only": ["Thermo"],
                  "right_only": ["Arriba"]}
    merge_col = ("_merge", "", "")
    found_col = ("found", "", "")
    merged = merge_results(thermo, arriba)
    merged[found_col] = [merge_dict[x] for x in merged[merge_col]]
    merged.drop(merge_col, axis=1, inplace=True)

    merged = merge_results(merged, star)
    merged[found_col] = [x if isinstance(x, list) else [] for x in merged[found_col]]
    merged[merge_col] = [[] if x == "left_only" else ["STARFusion"]
                         for x in merged[merge_col]]

    merged[("General", "FoundBy", "")] = merged[found_col] + merged[merge_col]
    merged[("General", "FoundBy", "")] = [", ".join(x)
                                          for x in merged[("General", "FoundBy", "")]]
    merged.drop([found_col, merge_col], axis=1, inplace=True)

    for i in "ab":
        merged.loc[pd.isna(merged.General.Exons[i]), ("General", "Exons", i)] = "None"

    merged = merged[[col for col in ["General", "Thermo", "Arriba", "STAR"]
                     if col in merged.columns]]
    return merged


def _read_keep_columns(keep_col_file):
    keep_col_file = files("cmd_fusion.parsers").joinpath(keep_col_file)
    with open(str(keep_col_file)) as infile:
        keep_cols = infile.readlines()
    keep_cols = [tuple(line.strip().split(",")) for line in keep_cols
                 if not line.startswith("#")]
    return keep_cols


def make_empty_table():
    empty = pd.DataFrame()
    for field in ["Genes", "Exons", "Breakpoints"]:
        for i in "ab":
            empty[("General", field, i)] = []
    empty.columns = pd.MultiIndex.from_tuples(empty.columns)
    empty.index = pd.MultiIndex.from_tuples([], names=["genes", "exons", "breakpoints"])
    return empty
