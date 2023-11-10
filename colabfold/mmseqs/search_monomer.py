"""
Functionality for running mmseqs locally. Takes in a fasta file, outputs final.a3m
Note: Currently needs mmseqs compiled from source

1. Compute monomer MSAs for all single proteins
    * mmseqs search
    * mmseqs expandaln
    * mmseqs align
     > save the alignment output after this step - it'll be an intermediate
    ...
    * mmseqs result2msa
    get output alignment
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


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "query",
        type=Path,
        help="fasta files with the queries.",
    )
    parser.add_argument(
        "dbbase",
        type=Path,
        help="The path to the database and indices you downloaded and created with setup_databases.sh",
    )
    parser.add_argument(
        "base", type=Path, help="Directory for the results (and intermediate files)"
    )
    parser.add_argument(
        "-s",
        type=float,
        default=None,
        help="MMseqs2 sensitivity. Lowering this will result in a much faster search but possibly sparser MSAs. By default, the k-mer threshold is directly set to the same one of the server, which corresponds to a sensitivity of ~8.",
    )
    # dbs are uniref, templates and environmental
    # We normally don't use templates
    parser.add_argument(
        "--db1", type=Path, default=Path("uniref30_2302_db"), help="UniRef database"
    )
    parser.add_argument("--db2", type=Path, default=Path(""), help="Templates database")
    parser.add_argument(
        "--db3",
        type=Path,
        default=Path("colabfold_envdb_202108_db"),
        help="Environmental database",
    )
    # # poor man's boolean arguments
    # parser.add_argument("--use-env", type=int, default=1, choices=[0, 1])
    # parser.add_argument("--use-templates", type=int, default=0, choices=[0, 1])
    parser.add_argument("--filter", type=int, default=1, choices=[0, 1])
    parser.add_argument(
        "--mmseqs",
        type=Path,
        default=Path("mmseqs"),
        help="Location of the mmseqs binary",
    )
    parser.add_argument("--expand-eval", type=float, default=math.inf)
    parser.add_argument("--align-eval", type=int, default=10)
    parser.add_argument("--diff", type=int, default=3000)
    parser.add_argument("--qsc", type=float, default=-20.0)
    parser.add_argument("--max-accept", type=int, default=1000000)
    parser.add_argument("--pairing_strategy", type=int, default=0)
    parser.add_argument("--db-load-mode", type=int, default=0)
    parser.add_argument("--threads", type=int, default=64)
    args = parser.parse_args()

    queries, is_complex = get_queries(args.query, None)
        """Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
        of job name, sequence and the optional a3m lines"""


    queries_unique = []
    for job_number, (raw_jobname, query_sequences, a3m_lines) in enumerate(queries):
        # remove duplicates before searching
        query_sequences = (
            [query_sequences] if isinstance(query_sequences, str) else query_sequences
        )
        query_seqs_unique = []
        for x in query_sequences:
            if x not in query_seqs_unique:
                query_seqs_unique.append(x)
        query_seqs_cardinality = [0] * len(query_seqs_unique)
        for seq in query_sequences:
            seq_idx = query_seqs_unique.index(seq)
            query_seqs_cardinality[seq_idx] += 1

        queries_unique.append([raw_jobname, query_seqs_unique, query_seqs_cardinality])

    args.base.mkdir(exist_ok=True, parents=True)
    query_file = args.base.joinpath("query.fas")
    with query_file.open("w") as f:
        for job_number, (
            raw_jobname,
            query_sequences,
            query_seqs_cardinality,
        ) in enumerate(queries_unique):
            for j, seq in enumerate(query_sequences):
                # The header of first sequence set as 101
                query_seq_headername = 101 + j
                f.write(f">{query_seq_headername}\n{seq}\n")

    run_mmseqs(
        args.mmseqs,
        ["createdb", query_file, args.base.joinpath("qdb"), "--shuffle", "0"],
    )
    with args.base.joinpath("qdb.lookup").open("w") as f:
        id = 0
        file_number = 0
        for job_number, (
            raw_jobname,
            query_sequences,
            query_seqs_cardinality,
        ) in enumerate(queries_unique):
            for seq in query_sequences:
                f.write(f"{id}\t{raw_jobname}\t{file_number}\n")
                id += 1
            file_number += 1

    mmseqs_search_monomer(
        mmseqs=args.mmseqs,
        dbbase=args.dbbase,
        base=args.base,
        uniref_db=args.db1,
        template_db=args.db2,
        metagenomic_db=args.db3,
        # use_env=args.use_env,
        # use_templates=args.use_templates,
        filter=args.filter,
        expand_eval=args.expand_eval,
        align_eval=args.align_eval,
        diff=args.diff,
        qsc=args.qsc,
        max_accept=args.max_accept,
        s=args.s,
        db_load_mode=args.db_load_mode,
        threads=args.threads,
    )
    # if is_complex is True:
    #     mmseqs_search_pair(
    #         mmseqs=args.mmseqs,
    #         dbbase=args.dbbase,
    #         base=args.base,
    #         uniref_db=args.db1,
    #         s=args.s,
    #         db_load_mode=args.db_load_mode,
    #         threads=args.threads,
    #         pairing_strategy=args.pairing_strategy,
    #     )

    #     id = 0
    #     for job_number, (
    #         raw_jobname,
    #         query_sequences,
    #         query_seqs_cardinality,
    #     ) in enumerate(queries_unique):
    #         unpaired_msa = []
    #         paired_msa = None
    #         if len(query_seqs_cardinality) > 1:
    #             paired_msa = []
    #         for seq in query_sequences:
    #             with args.base.joinpath(f"{id}.a3m").open("r") as f:
    #                 unpaired_msa.append(f.read())
    #             args.base.joinpath(f"{id}.a3m").unlink()
    #             if len(query_seqs_cardinality) > 1:
    #                 with args.base.joinpath(f"{id}.paired.a3m").open("r") as f:
    #                     paired_msa.append(f.read())
    #             args.base.joinpath(f"{id}.paired.a3m").unlink()
    #             id += 1
    #         msa = msa_to_str(
    #             unpaired_msa, paired_msa, query_sequences, query_seqs_cardinality
    #         )
    #         args.base.joinpath(f"{job_number}.a3m").write_text(msa)
    #         # add raw_jobname to the output file
    #         os.rename(
    #             args.base.joinpath(f"{job_number}.a3m"),
    #             args.base.joinpath(f"{safe_filename(raw_jobname)}.a3m"),
    #         )
    #         if args.use_templates:
    #             os.rename(
    #                 args.base.joinpath(f"{args.db2}.m8"),
    #                 args.base.joinpath(f"{safe_filename(raw_jobname)}_{args.db2}.m8"),
    #             )

    query_file.unlink()
    run_mmseqs(args.mmseqs, ["rmdb", args.base.joinpath("qdb")])
    run_mmseqs(args.mmseqs, ["rmdb", args.base.joinpath("qdb_h")])


if __name__ == "__main__":
    main()
