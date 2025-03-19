
from glob import glob
from io import StringIO
import os
from importlib.resources import files
import yaml

import pandas as pd

from cmd_fusion.parsers import ThermoFusions, ArribaFusions, StarFusions, VCFFusions
from cmd_fusion.parsers.common import merge_all_results


def _read_all(thermo, arriba, starfusion, geneset, thermo_as_vcf=False):
    thermo = VCFFusions(thermo, geneset) if thermo_as_vcf else ThermoFusions(thermo, geneset)
    arriba = ArribaFusions(arriba, geneset)
    star = StarFusions(starfusion, geneset)
    merged = merge_all_results(thermo.table, arriba.table, star.table)
    merged.sort_index(inplace=True)
    merged["numdex"] = list(range(len(merged)))
    merged.set_index("numdex", append=True, inplace=True)
    return merged


def _read_default_cols():
    default_col_file = files("cmd_fusion.parsers").joinpath("columns/default_columns.csv")
    with open(str(default_col_file)) as infile:
        default_cols = [tuple(line.strip().split(",")) for line in infile.readlines()
                        if line.strip() and not line.startswith("#")]
    return default_cols


def _read_display_names():
    display_name_file = files("cmd_fusion.parsers").joinpath("columns/display_names.yaml")
    with open(str(display_name_file)) as infile:
        display_names = yaml.safe_load(infile)
    display_names = {tuple(x.split(",")): y.split(",") for x, y in display_names.items()}
    return display_names


def _retrieve_cases(results_dir):
    cases = sorted({case.split(".")[0] for case in os.listdir(results_dir)})
    return cases


def _read_case(case, results_dir, geneset):
    case_files = dict()
    thermo_as_vcf = False
    for tool in ["arriba", "starfusion"]:
        case_file = f"{results_dir}/{case}.{tool}.tsv"
        if os.path.isfile(case_file):
            case_files[tool] = case_file
        else:
            case_files[tool] = None
    if os.path.isfile(f"{results_dir}/{case}.thermo.tsv"):
        case_files["thermo"] = f"{results_dir}/{case}.thermo.tsv"
    elif os.path.isfile(f"{results_dir}/{case}.thermo.bcf"):
        case_files["thermo"] = f"{results_dir}/{case}.thermo.bcf"
        thermo_as_vcf = True
    else:
        case_files["thermo"] = None

    # files = [file if os.path.isfile(file := f"{results_dir}/{case}.{tool}.tsv") else None
    #          for tool in ["thermo", "arriba", "starfusion"]]
    table = _read_all(**case_files, geneset=geneset, thermo_as_vcf=thermo_as_vcf)
    for col in ["Genes", "Exons", "Breakpoints"]:
        for i in "ab":
            table[("General", col, i)] = table[("General", col, i)].astype(str)
    table = table.to_json()
    return table


def _load_json_table(table_json):
    table = pd.read_json(StringIO(table_json))
    table.columns = pd.MultiIndex.from_tuples([eval(col) for col in table.columns])
    table.index = pd.MultiIndex.from_tuples([eval(dex) for dex in table.index],
                                            names=["genes", "exons", "breakpoints", "dex"])
    return table
