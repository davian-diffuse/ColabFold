"""
Functionality for running mmseqs locally. Takes in a fasta file, outputs final.a3m
Note: Currently needs mmseqs compiled from source


Getting paired MSAs for PPIs:

1. Compute monomer MSAs for all single proteins
    * mmseqs search
    * mmseqs expandaln
    * mmseqs align
     > save the alignment output after this step - it'll be an intermediate
    ...
    * mmseqs result2msa
    get output alignment
    
2. Compute paired MSAs for all interactions
    * figure out how to concat two alignment files (check convertalis for human readable)
        "each data record consists of the alignment results for one query sequence"
    * mmseqs pairaln
    * mmseqs align
    * mmseqs pairaln
    * mmseqs result2msa
    get output paired alignment
"""

import logging
import math
import os
import shutil
import subprocess
from argparse import ArgumentParser
from pathlib import Path
from typing import List, Union

from colabfold.batch import get_queries, msa_to_str
from colabfold.utils import safe_filename

logger = logging.getLogger(__name__)


def run_mmseqs(mmseqs: Path, params: List[Union[str, Path]]):
    params_log = " ".join(str(i) for i in params)
    logger.info(f"Running {mmseqs} {params_log}")
    subprocess.check_call([mmseqs] + params)


def mmseqs_search_monomer(
    dbbase: Path,
    base: Path,
    uniref_db: Path = Path("uniref30_2302_db"),
    template_db: Path = Path(""),  # Unused by default
    metagenomic_db: Path = Path("colabfold_envdb_202108_db"),
    mmseqs: Path = Path("mmseqs"),
    # use_env: bool = True,
    # use_templates: bool = False,
    filter: bool = True,
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    s: float = 8,
    db_load_mode: int = 2,
    threads: int = 32,
):
    """Run mmseqs with a local colabfold database set

    db1: uniprot db (UniRef30)
    db2: Template (unused by default)
    db3: metagenomic db (colabfold_envdb_202108 or bfd_mgy_colabfold, the former is preferred)
    """
    if filter:
        # 0.1 was not used in benchmarks due to POSIX shell bug in line above
        #  EXPAND_EVAL=0.1
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    used_dbs = [uniref_db]
    if use_templates:
        used_dbs.append(template_db)
    if use_env:
        used_dbs.append(metagenomic_db)

    for db in used_dbs:
        if not dbbase.joinpath(f"{db}.dbtype").is_file():
            raise FileNotFoundError(f"Database {db} does not exist")
        if (
            not dbbase.joinpath(f"{db}.idx").is_file()
            and not dbbase.joinpath(f"{db}.idx.index").is_file()
        ):
            logger.info("Search does not use index")
            db_load_mode = 0
            dbSuffix1 = "_seq"
            dbSuffix2 = "_aln"
        else:
            dbSuffix1 = ".idx"
            dbSuffix2 = ".idx"

    # fmt: off
    # @formatter:off
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs", "10000"]
    if s is not None:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]
    filter_param = ["--filter-msa", str(filter), "--filter-min-enable", "1000", "--diff", str(diff), "--qid", "0.0,0.2,0.4,0.6,0.8,1.0", "--qsc", "0", "--max-seq-id", "0.95",]
    expand_param = ["--expansion-mode", "0", "-e", str(expand_eval), "--expand-filter-clusters", str(filter), "--max-seq-id", "0.95",]

    run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), dbbase.joinpath(uniref_db), base.joinpath("res"), base.joinpath("tmp"), "--threads", str(threads)] + search_param)
    run_mmseqs(mmseqs, ["mvdb", base.joinpath("tmp/latest/profile_1"), base.joinpath("prof_res")])
    run_mmseqs(mmseqs, ["lndb", base.joinpath("qdb_h"), base.joinpath("prof_res_h")])
    run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res"), dbbase.joinpath(f"{uniref_db}{dbSuffix2}"), base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + expand_param)
    run_mmseqs(mmseqs, ["align", base.joinpath("prof_res"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp"), base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads", str(threads), "--alt-ali", "10", "-a"])
    run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"),
                        base.joinpath("res_exp_realign"), base.joinpath("res_exp_realign_filter"), "--db-load-mode",
                        str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0", "--threads",
                        str(threads), "--max-seq-id", "1.0", "--filter-min-enable", "100"])
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"),
                        base.joinpath("res_exp_realign_filter"), base.joinpath("uniref.a3m"), "--msa-format-mode",
                        "6", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)
    
    # subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign")]) #! save this file as intermediate
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp")])  
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res")])
    # subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign_filter")])
    run_mmseqs(mmseqs, ["mvdb", base.joinpath("uniref.a3m"), base.joinpath("final.a3m")])
    run_mmseqs(mmseqs, ["unpackdb", base.joinpath("final.a3m"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".a3m"])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("final.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("uniref.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res")])
    # @formatter:on
    # fmt: on

    for file in base.glob("prof_res*"):
        file.unlink()
    shutil.rmtree(base.joinpath("tmp"))
    if use_templates:
        shutil.rmtree(base.joinpath("tmp2"))
    if use_env:
        shutil.rmtree(base.joinpath("tmp3"))




def mmseqs_search_pair(
    dbbase: Path,
    base: Path,
    uniref_db: Path = Path("uniref30_2302_db"),
    mmseqs: Path = Path("mmseqs"),
    s: float = 8,
    threads: int = 64,
    db_load_mode: int = 2,
    pairing_strategy: int = 0,
):
    if not dbbase.joinpath(f"{uniref_db}.dbtype").is_file():
        raise FileNotFoundError(f"Database {uniref_db} does not exist")
    if (
        not dbbase.joinpath(f"{uniref_db}.idx").is_file()
        and not dbbase.joinpath(f"{uniref_db}.idx.index").is_file()
    ):
        logger.info("Search does not use index")
        db_load_mode = 0
        dbSuffix1 = "_seq"
        dbSuffix2 = "_aln"
    else:
        dbSuffix1 = ".idx"
        dbSuffix2 = ".idx"

    # fmt: off
    # @formatter:off
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs", "10000",]
    if s is not None:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]
    expand_param = ["--expansion-mode", "0", "-e", "inf", "--expand-filter-clusters", "0", "--max-seq-id", "0.95",]
    
    #! precomputed
    # run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), dbbase.joinpath(uniref_db), base.joinpath("res"), base.joinpath("tmp"), "--threads", str(threads),] + search_param,)
    # run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res"), dbbase.joinpath(f"{uniref_db}{dbSuffix2}"), base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode), "--threads", str(threads),] + expand_param,)
    # run_mmseqs(mmseqs, ["align", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp"), base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", "0.001", "--max-accept", "1000000", "--threads", str(threads), "-c", "0.5", "--cov-mode", "1",],)
    
    #! res_exp_realign = concat(res_exp_realign1, res_exp_realign2)
    run_mmseqs(mmseqs, ["pairaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}"), base.joinpath("res_exp_realign"), base.joinpath("res_exp_realign_pair"), "--db-load-mode", str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "0", "--threads", str(threads), ],)
    run_mmseqs(mmseqs, ["align", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp_realign_pair"), base.joinpath("res_exp_realign_pair_bt"), "--db-load-mode", str(db_load_mode), "-e", "inf", "-a", "--threads", str(threads), ],)
    run_mmseqs(mmseqs, ["pairaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}"), base.joinpath("res_exp_realign_pair_bt"), base.joinpath("res_final"), "--db-load-mode", str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "1", "--threads", str(threads),],)
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_final"), base.joinpath("pair.a3m"), "--db-load-mode", str(db_load_mode), "--msa-format-mode", "5", "--threads", str(threads),],)
    run_mmseqs(mmseqs, ["unpackdb", base.joinpath("pair.a3m"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".paired.a3m",],)
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("qdb")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("qdb_h")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp_realign")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp_realign_pair")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp_realign_pair_bt")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_final")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("pair.a3m")])
    shutil.rmtree(base.joinpath("tmp"))
    # @formatter:on
    # fmt: on


