#!/usr/bin/env python3
"""SalmonTE - Ultra-Fast and Scalable Quantification Pipeline of Transcript Abundances from Next Generation Sequencing Data

Usage:
    SalmonTE.py index [--ref_name=ref_name] (--input_fasta=fa_file) [--te_only]
    SalmonTE.py index_loci (--ref_name=ref_name) (--genome=genome_fa) (--te_gtf=te_gtf) [--min_length=min_length] [--num_threads=numthreads]
    SalmonTE.py quant [--reference=genome] [--outpath=outpath] [--num_threads=numthreads] [--exprtype=exprtype] FILE...
    SalmonTE.py test [--inpath=inpath] [--outpath=outpath] [--tabletype=tabletype] [--figtype=figtype] [--analysis_type=analysis_type] [--conditions=conditions]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.

index_loci builds a per-locus index: every TE copy listed in a TEtranscripts-
format TE GTF (gene_id / transcript_id / family_id / class_id attributes, e.g.
hg38_rmsk_TE.gtf) is cut out of the genome FASTA and indexed as its own Salmon
target. `quant` against such an index additionally writes EXPR_LOCUS.csv
(per-locus values with genomic coordinates) and EXPR_SUBFAMILY.csv (loci summed
per subfamily). --min_length defaults to 31 (Salmon's k-mer size).
"""
from docopt import docopt
from glob import glob
from itertools import combinations
import logging
import sys
import os
import tempfile
import gzip
import shutil


def resolve_salmon_bin():
    salmon_bin = os.environ.get("SALMONTE_SALMON_BIN")
    if salmon_bin:
        return salmon_bin

    # Prefer Homebrew salmon on macOS when available to avoid picking a stale
    # conda binary from PATH.
    for candidate in ["/usr/local/bin/salmon", "/opt/homebrew/bin/salmon"]:
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return candidate

    path_salmon = shutil.which("salmon")
    if path_salmon:
        return path_salmon

    return os.path.join(os.path.dirname(__file__), "salmon/{}/bin/salmon".format(sys.platform))


def _raise_run_salmon_error(param, salmon_bin):
    index_dir = param["--reference"]
    old_index = os.path.exists(os.path.join(index_dir, "refInfo.json")) and not os.path.exists(os.path.join(index_dir, "info.json"))
    guidance = ""
    if old_index:
        guidance = (
            " Detected an old reference index format (refInfo.json). "
            "This index is not compatible with modern Salmon quant binaries. "
            "Use a Salmon version compatible with RapMap indexes or rebuild the index with the same Salmon version used for quantification."
        )

    raise RuntimeError(
        "Salmon/Snakemake quantification failed before EXPR.csv was created. "
        "Inspect Snakemake log files under .snakemake/log and rerun without stderr suppression for details. "
        "Salmon binary: {}. Reference index: {}.{}".format(salmon_bin, index_dir, guidance)
    )


def is_reference_index_dir(index_dir):
    markers = ["versionInfo.json", "refInfo.json", "info.json", "index.ssi"]
    return os.path.isdir(index_dir) and any(os.path.exists(os.path.join(index_dir, marker)) for marker in markers)

FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def is_fastq(file_name):
    # Directories are never FASTQ files, and the extension must follow a dot:
    # a folder called "fq/" or "my_fastq/" used to be mistaken for a FASTQ file.
    if os.path.isdir(file_name):
        return False
    return file_name.lower().endswith(FASTQ_EXTS)


def get_basename_noext(file_name):
    # return os.path.basename(file_name.split('.')[0])
    return os.path.basename(file_name).split('.')[0]


def get_ext(file_name):
    # return ".".join(os.path.basename(file_name).split('.')[1:])
    return ".".join(os.path.basename(file_name).split('.')[1:])

def longest_prefix(a, b):
    a, b = (b,a) if len(a) > len(b) else (a,b)
    return b[:max([i for i in range(len(a),-1,-1) if b.startswith(a[:i])]+[0])]


def get_first_readid(file_name):
    if file_name.lower().endswith(".gz"):
        with gzip.open(file_name, "rb") as inp:
            return inp.readline().decode("utf8").replace("/"," ").split()[0]
    else:
        with open(file_name, "r") as inp:
            return inp.readline().replace("/"," ").split()[0]


def correct_ext(ext):
    return "fastq.gz" if "gz" in ext else "fastq"

def collect_FASTQ_files(FILE):
    fastq_files = set()
    for file in FILE:
        collected_files = glob(file)
        if len(collected_files) == 0:
            logging.error("Failed to read : '{}' is not existed.".format(file))
            sys.exit(1)

        elif len(collected_files) == 1 and not is_fastq(collected_files[0]):
            logging.info("SalmonTE assumes that '{}' is a directory, and SalmonTE will search any FASTQ file in the directory.".format(collected_files[0]))
            collected_files = glob(collected_files[0]+"/*")

        black_list = set()
        for col_file in collected_files:
            if not is_fastq(col_file):
                # logging.error("Failed to read : '{}' seems not to be a FASTQ file.".format(file))
                logging.warn("SalmonTE found '{}' file, but this seems not to be a FASTQ file.".format(col_file))
                black_list.add(col_file)

        fastq_files |= set([os.path.abspath(file) for file in collected_files if file not in black_list])

    tmp_dir = tempfile.mkdtemp()

    paired = dict([ (file,set()) for file in fastq_files ])
    for file_a, file_b in combinations(fastq_files, 2):
        """
        if get_ext(file_a) != get_ext(file_b):
            logging.info("File extensions of all files must be same.")
            sys.exit(1)
        """

        if get_first_readid(file_a) == get_first_readid(file_b):
            paired[file_a].add(file_b)
            paired[file_b].add(file_a)

    is_paired = True
    for file in paired:
        if len(paired[file]) == 0: is_paired = False



    for file in paired:
        if is_paired and len(paired[file]) > 1:
            logging.error("One of input files can be paired to multiple files")
            sys.exit(1)
        elif not is_paired and len(paired[file]) > 0:
            logging.error("A paired-end sample and a single-end sample are placed together.")
            sys.exit(1)

    file_list = []
    if is_paired:
        for file, matched in paired.items():
            if len(matched) == 0:
                logging.error("Input dataset is supposed as a paired-ends dataset, but some files do not have thier pair.")
                sys.exit(1)

        logging.info("The input dataset is considered as a paired-ends dataset.")
        for file in sorted(paired):
            a = file
            b = list(paired[file])[0]
            if a > b: continue
            trim_a = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R1.{}".format(correct_ext(os.path.basename(a).split('.')[1:]))
            trim_b = "_".join(get_basename_noext(a).split("_")[:-1]) + "_R2.{}".format(correct_ext(os.path.basename(b).split('.')[1:]))
            os.symlink(os.path.abspath(a), os.path.join(tmp_dir, trim_a))
            os.symlink(os.path.abspath(b), os.path.join(tmp_dir, trim_b))
            file_list.append([(trim_a, trim_b)])
    else:
        logging.info("The input dataset is considered as a single-end dataset.")
        for file in sorted(fastq_files):
            file_name = os.path.join(tmp_dir, get_basename_noext(file)) + ".{}".format(correct_ext(os.path.basename(file).split('.')[1:]))
            os.symlink(os.path.abspath(file), file_name)
            file_list.append(os.path.basename(file))

    ret = dict()
    ret["paired"] = is_paired
    ret["inpath"] = tmp_dir
    ret["num_fastq"] = int(len(fastq_files) / (is_paired*2)) if is_paired else len(fastq_files)
    ret["file_list"] = file_list
    return ret

def run_salmon(param):
    os.environ.setdefault("PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION", "python")
    import snakemake
    snakefile = os.path.join(os.path.dirname(__file__), "snakemake/Snakefile.paired" if param["paired"] else "snakemake/Snakefile.single")
    salmon_bin = resolve_salmon_bin()

    ok = snakemake.snakemake(
        snakefile=snakefile,
        config={
            "input_path": param["inpath"],
            "output_path": param["--outpath"],
            "index": param["--reference"],
            "salmon": salmon_bin,
            "num_threads" : param["--num_threads"],
            "exprtype": param["--exprtype"],
        },
        quiet=True,
        lock=False
    )

    expr_file = os.path.join(param["--outpath"], "EXPR.csv")
    if not ok or not os.path.exists(expr_file):
        _raise_run_salmon_error(param, salmon_bin)

    with open(expr_file, "r") as inp:
        sample_ids = inp.readline().strip().split(',')[1:]
    with open(os.path.join(param["--outpath"], "condition.csv" ), "w") as oup:
        oup.write("SampleID,condition\n")
        oup.write(
            "\n".join([s+","+"NA" for s in sample_ids]) + "\n"
        )


def build_salmon_index(in_fasta, ref_id, te_only):
    logging.info("Building Salmon Index")
    if not os.path.exists(in_fasta):
        logging.error("Input file does not exist.")
        sys.exit(1)
    ext = get_ext(in_fasta)
    if not ext.lower() in [ "fa", "fasta" ]:
        logging.error("Input file is not a FASTA file.")
        sys.exit(1)
    ref_path = os.path.join(os.path.dirname(__file__), "reference")
    out_path = os.path.join(ref_path, ref_id)

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    clade_dict = {}
    with open(os.path.join(ref_path, "clades_extended.csv"), "r") as inp:
        for line in inp:
            anno, clade, repeat_class, category = line.strip().split(",")
            if te_only and category != "Transposable Element": continue
            clade_dict[anno] = [clade, repeat_class]


    file_fa = os.path.join(out_path, "ref.fa")
    file_clade = os.path.join(out_path, "clades.csv")

    with open(in_fasta, "r") as inp, \
         open(file_fa, "w") as oup_fa, \
         open(file_clade, "w") as oup_clade:
        oup_clade.write("name,class,clade\n")
        for line in inp:
            if line.startswith(">"):
                parts = line[1:].strip().split()
                if len(parts) == 0:
                    continue
                name = parts[0]
                anno = parts[1] if len(parts) > 1 else parts[0]
                if anno in clade_dict:
                    clade, repeat_class = clade_dict[anno]
                else: 
                    clade = "other"
                    repeat_class = "other"

                oup_fa.write(">{}\n".format(name))
                oup_clade.write("{},{},{}\n".format(name,
                                                   clade,
                                                   repeat_class))
            else:
                oup_fa.write(line.strip()+"\n")

    salmon_bin = resolve_salmon_bin()
    cmd = "{} index --transcripts {} --index {}".format(salmon_bin, file_fa, out_path)
    ret = os.system(cmd)
    if ret != 0:
        logging.error("Failed to build '{}' index using command: {}".format(ref_id, cmd))
        sys.exit(1)
    logging.info("Building '{}' index was finished!".format(ref_id))


LOCI_FILE = "loci.csv"

# class_id -> clade label, matching the labels used in the bundled clades.csv
CLASS_TO_CLADE = {
    "LINE": "Non-LTR Retrotransposon",
    "SINE": "Non-LTR Retrotransposon",
    "Retroposon": "Non-LTR Retrotransposon",
    "LTR": "LTR Retrotransposon",
    "DNA": "DNA transposon",
    "RC": "DNA transposon",
}

_REVCOMP = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def _gtf_attr(attr, key):
    tag = key + ' "'
    i = attr.find(tag)
    if i < 0:
        return None
    i += len(tag)
    return attr[i:attr.find('"', i)]


def read_te_loci(te_gtf):
    """Parse a TEtranscripts-format TE GTF into per-locus records.

    Locus IDs follow TElocal / TEfastR: transcript_id:gene_id:family_id:class_id,
    made unique with a _N suffix if the GTF repeats one (R's make.unique(sep="_")).
    """
    loci, seen = [], {}
    opener = gzip.open if te_gtf.endswith(".gz") else open
    with opener(te_gtf, "rt") as inp:
        for line in inp:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid, sub, fam, cls = (_gtf_attr(f[8], k) for k in
                                  ("transcript_id", "gene_id", "family_id", "class_id"))
            if None in (tid, sub, fam, cls):
                logging.error("TE GTF line lacks gene_id/transcript_id/family_id/class_id: %s", line.strip())
                sys.exit(1)
            name = ":".join((tid, sub, fam, cls))
            if name in seen:
                seen[name] += 1
                name = "{}_{}".format(name, seen[name])
            else:
                seen[name] = 0
            start, end = int(f[3]), int(f[4])
            loci.append({"name": name, "chrom": f[0], "start": start, "end": end,
                         "strand": f[6], "transcript_id": tid, "subfamily": sub,
                         "family": fam, "class": cls, "length": end - start + 1})
    return loci


def _iter_fasta(path):
    """Yield (chrom, sequence) one chromosome at a time."""
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as inp:
        for line in inp:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def build_locus_index(genome_fa, te_gtf, ref_id, min_length, num_threads):
    """Build a Salmon index whose targets are individual genomic TE copies."""
    for p in (genome_fa, te_gtf):
        if not os.path.exists(p):
            logging.error("Input file does not exist: %s", p)
            sys.exit(1)
    min_length = int(min_length) if min_length else 31
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference", ref_id)
    os.makedirs(out_path, exist_ok=True)

    logging.info("Reading TE loci from %s", te_gtf)
    loci = read_te_loci(te_gtf)
    by_chrom = {}
    for rec in loci:
        by_chrom.setdefault(rec["chrom"], []).append(rec)
    logging.info("%d TE loci on %d sequences", len(loci), len(by_chrom))

    file_fa = os.path.join(out_path, "ref.fa")
    written, n_short, found = [], 0, set()
    logging.info("Extracting TE sequences from %s", genome_fa)
    with open(file_fa, "w") as oup:
        for chrom, seq in _iter_fasta(genome_fa):
            if chrom not in by_chrom:
                continue
            found.add(chrom)
            for rec in by_chrom[chrom]:
                s = seq[rec["start"] - 1:rec["end"]]
                if len(s) < min_length:
                    n_short += 1
                    continue
                if rec["strand"] == "-":
                    s = s.translate(_REVCOMP)[::-1]
                oup.write(">{}\n{}\n".format(rec["name"], s))
                written.append(rec)
    missing = set(by_chrom) - found
    if missing:
        logging.warning("%d GTF sequence(s) not in the genome FASTA (chromosome naming mismatch?), "
                        "%d loci skipped, e.g. %s", len(missing),
                        sum(len(by_chrom[c]) for c in missing), sorted(missing)[:3])
    if n_short:
        logging.warning("%d loci shorter than %d bp skipped", n_short, min_length)
    if not written:
        logging.error("No TE loci were extracted; check that the GTF and genome use the same chromosome names.")
        sys.exit(1)

    cols = ["name", "chrom", "start", "end", "strand", "transcript_id", "subfamily", "family", "class", "length"]
    with open(os.path.join(out_path, LOCI_FILE), "w") as oup:
        oup.write(",".join(cols) + "\n")
        for rec in written:
            oup.write(",".join(str(rec[c]) for c in cols) + "\n")
    # clades.csv keyed by locus name, so `SalmonTE.py test` runs on locus-level EXPR.csv
    with open(os.path.join(out_path, "clades.csv"), "w") as oup:
        oup.write("name,class,clade\n")
        for rec in written:
            oup.write("{},{},{}\n".format(rec["name"], rec["class"],
                                          CLASS_TO_CLADE.get(rec["class"], rec["class"])))

    # --keepDuplicates: identical young copies must stay separate targets, otherwise
    # Salmon collapses them and those loci vanish from quant.sf.
    salmon_bin = resolve_salmon_bin()
    cmd = "{} index --transcripts {} --index {} --keepDuplicates -p {}".format(
        salmon_bin, file_fa, out_path, num_threads or 4)
    logging.info("Indexing %d TE loci: %s", len(written), cmd)
    if os.system(cmd) != 0:
        logging.error("Failed to build '%s' locus index", ref_id)
        sys.exit(1)
    logging.info("Building '%s' locus index was finished!", ref_id)


def write_locus_tables(reference, outpath):
    """For a per-locus index, add coordinates to EXPR.csv (EXPR_LOCUS.csv) and
    sum loci per subfamily (EXPR_SUBFAMILY.csv)."""
    loci_file = os.path.join(reference, LOCI_FILE)
    if not os.path.exists(loci_file):
        return
    import pandas as pd
    loci = pd.read_csv(loci_file, dtype={"chrom": str})
    expr = pd.read_csv(os.path.join(outpath, "EXPR.csv"), index_col="TE")
    samples = list(expr.columns)
    locus = loci.merge(expr, left_on="name", right_index=True, how="left")
    locus[samples] = locus[samples].fillna(0.0)
    locus = locus.rename(columns={"name": "locus_id"})
    locus.to_csv(os.path.join(outpath, "EXPR_LOCUS.csv"), index=False)
    sub = locus.groupby(["subfamily", "family", "class"], sort=True)[samples].sum().reset_index()
    sub.to_csv(os.path.join(outpath, "EXPR_SUBFAMILY.csv"), index=False)
    logging.info("Wrote per-locus table with coordinates (EXPR_LOCUS.csv, %d loci) "
                 "and per-subfamily sums (EXPR_SUBFAMILY.csv)", len(locus))


def run(args):
    if args['index_loci']:
        build_locus_index(args['--genome'], args['--te_gtf'], args['--ref_name'],
                          args['--min_length'], args['--num_threads'])

    if args['quant']:
        if args['--exprtype'] is None:
            args['--exprtype'] = "TPM"
        if args['--num_threads'] is None:
            args['--num_threads'] = 4
        if args['--outpath'] is None:
            args['--outpath'] = os.path.join(os.getcwd(), "SalmonTE_output/")
        ref_arg = args['--reference']
        if ref_arg is None or ref_arg == "hs":
            args['--reference'] = os.path.join(os.path.dirname(__file__), "reference/hs")
        elif is_reference_index_dir(ref_arg):
            args['--reference'] = os.path.abspath(ref_arg)
        else:
            bundled_ref = os.path.join(os.path.dirname(__file__), "reference", ref_arg)
            if is_reference_index_dir(bundled_ref):
                args['--reference'] = bundled_ref
            else:
                logging.error("Reference file is not found!")
                sys.exit(1)

        if args['--reference'].startswith("./"):
           args['--reference'] = args['--reference'][2:]
           
        logging.info("Starting quantification mode")
        logging.info("Collecting FASTQ files...")
        param = {**args, **collect_FASTQ_files(args['FILE'])}
        logging.info("Collected {} FASTQ files.".format(param["num_fastq"]))
        logging.info("Quantification has been finished.")
        logging.info("Running Salmon using Snakemake")

        run_salmon(param)
        os.system("cp {}/clades.csv {}".format(args['--reference'],
                                               args['--outpath']))
        write_locus_tables(args['--reference'], args['--outpath'])

    if args['index']:
        build_salmon_index(args['--input_fasta'], args['--ref_name'], args['--te_only'])

    if args['test']:
        import snakemake
        if args['--inpath'] is None:
            logging.error("Input path must be specified!")
            sys.exit(1)
        elif not os.path.exists(os.path.join(args['--inpath'], "EXPR.csv")):
            logging.error("Input path is specified incorrectly!")
            sys.exit(1)

        if args['--outpath'] is None:
            snakemake.utils.makedirs(os.path.join(os.getcwd(), "SalmonTE_output"))
        elif not os.path.exists(args['--outpath']):
            snakemake.utils.makedirs(args['--outpath'])

        if args['--tabletype'] is None:
            args['--tabletype'] = "xls"

        if args['--figtype'] is None:
            args['--figtype'] = "pdf"

        if args['--analysis_type'] is None:
            args['--analysis_type'] = "DE"
        args['--analysis_type'] = args['--analysis_type'].upper()
        if args['--analysis_type'] != "DE" and args['--analysis_type'] != "LM":
            logging.error("analysis_type must be set as 'DE' or 'LM'.")

        if args['--conditions'] is None:
            args['--conditions'] = ""
        
        #  Rscript SalmonTE_Stats.R SalmonTE_output xls PDF tmp
        ret = os.system("Rscript {} {} {} {} {} {} {}".format( os.path.join(os.path.dirname(__file__), "SalmonTE_Stats.R"),
                                        args["--inpath"],
                                        args['--tabletype'],
                                        args['--figtype'],
                                        args['--outpath'],
                                        args['--analysis_type'],
                                        args['--conditions']))
        if ret != 0:
            logging.error("The statistical analysis (SalmonTE_Stats.R) failed; see the R error above.")
            sys.exit(1)


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    args = docopt(__doc__, version='SalmonTE 0.4.1 (sulevk fork)')
    run(args)
