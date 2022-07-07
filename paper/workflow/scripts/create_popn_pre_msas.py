import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
from typing import TextIO, Set, Dict
from tempfile import TemporaryDirectory
from loguru import logger
from cyvcf2 import VCF
import subprocess
from pysam import FastxFile, FastaFile

TRANSLATE = str.maketrans("ATGC", "TACG")


def revcomp(s: str) -> str:
    return s.upper().translate(TRANSLATE)[::-1]


#############
# MAIN
############
fasta_ref = str(snakemake.input.references)
faidx = FastaFile(fasta_ref)
outdir = Path(snakemake.output[0]).absolute()
outdir.mkdir(exist_ok=True)
vcf_fname = str(snakemake.input.vcf)
vcf_rdr = VCF(vcf_fname)
samples_fname = snakemake.input.get("samples")
samples = set(vcf_rdr.samples)
vcf_rdr.close()

if samples_fname is not None:
    logger.info("Loading sample names from file...")
    with open(samples_fname) as f:
        arr = [x for x in f.readlines() if x]
        s = set()
        for sample in arr:
            if sample not in samples:
                logger.warning(f"Sample {sample} is not in the VCF [skipping]")
            else:
                s.add(sample)
        if len(s) == 0:
            logger.warning(
                "No valid samples found in {} - using all samples in VCF", samples_fname
            )
        else:
            samples = s
else:
    logger.info("Using all samples in VCF")
logger.info(f"Loaded {len(samples)} samples")

logger.info("Determining which strand each gene is on...")
strand: Dict[str, str] = dict()
for entry in FastxFile(fasta_ref):
    for field in entry.comment.rstrip().split():
        if field.startswith("strand"):
            strand[entry.name] = field[7]
    if entry.name not in strand:
        raise ValueError(f"Couldn't find strand for {entry.name}")

logger.info("Extracting consensus sequences with bcftools consensus...")
with TemporaryDirectory() as tmpdirname:
    for sample in samples:
        outname = Path(tmpdirname) / f"{sample}.fa"
        args = (
            "bcftools",
            "consensus",
            "-s",
            sample,
            "-H",
            "A",
            "-o",
            outname,
            "-f",
            fasta_ref,
            vcf_fname,
        )
        logger.debug(f"{sample} run with: {args}")

        subprocess.run(args, check=True, stderr=sys.stderr)

    logger.success("Consensus sequences created")

    logger.info("Combining consensus sequences into pre-MSA fasta files...")
    files: Dict[str, TextIO] = {}
    outdir.mkdir(exist_ok=True)

    for sample in samples:
        logger.debug("Aggregating consensus sequences for {}...", sample)
        ifname = str(Path(tmpdirname) / f"{sample}.fa")
        for entry in FastxFile(ifname):
            chrom = entry.name
            is_fwd_strand = strand[chrom] == "+"

            if chrom not in files:
                refseq = faidx.fetch(reference=chrom)
                ref_header = f">{chrom}_reference"
                seq = refseq if is_fwd_strand else revcomp(refseq)
                fh = open(outdir / f"{chrom}.fa", "w")
                files[chrom] = fh
                print(f"{ref_header}\n{seq}", file=fh)
            else:
                fh = files[chrom]

            new_header = f">{chrom}_sample={sample}"
            seq = entry.sequence if is_fwd_strand else revcomp(entry.sequence)
            print(f"{new_header}\n{seq}", file=fh)

    for f in files.values():
        f.close()

logger.success("Done")
