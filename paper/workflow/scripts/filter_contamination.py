import logging
from enum import Enum
from pathlib import Path
from typing import Set

import click
import pandas as pd
import pysam

IS_CONTAM_COL = "is_contaminant"


class RequiredIf(click.Option):
    def __init__(self, *args, **kwargs):
        self.required_if = kwargs.pop("required_if")
        assert self.required_if, "'required_if' parameter required"
        kwargs["help"] = (
            kwargs.get("help", "")
            + " NOTE: This argument is mutually inclusive with %s" % self.required_if
        ).strip()
        super(RequiredIf, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        we_are_present = self.name in opts
        other_present = self.required_if in opts

        if not other_present:
            if we_are_present:
                raise click.UsageError(
                    "Illegal usage: `%s` is mutually inclusive with `%s`"
                    % (self.name, self.required_if)
                )

        return super(RequiredIf, self).handle_parse_result(ctx, opts, args)


class Classification(Enum):
    Contaminant = "contaminant"
    Unmaped = "unmapped"
    Keep = "keep"
    Other = "other"


class Classifier:
    def __init__(self, metadata_file: str):
        self.metadata = pd.read_table(
            metadata_file,
            header=None,
            names=["organism", IS_CONTAM_COL, "accession"],
            index_col="accession",
            dtype={IS_CONTAM_COL: "bool"},
        )

    def classify(self, record: pysam.AlignedSegment) -> Classification:
        if record.is_unmapped:
            return Classification.Unmaped

        ref_id = record.reference_name
        is_contam: bool = self.metadata.at[ref_id, IS_CONTAM_COL]
        if is_contam:
            return Classification.Contaminant

        return Classification.Keep


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--samfile",
    help="{B,CR,S}AM file of reads mapped to decontamination database.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-m",
    "--metadata",
    help=(
        "TSV file containing information about each reference in the database. Column "
        "1 is the category, column 2 is whether the reference is contamination, column "
        "3 is the accession ID for the reference."
    ),
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(dir_okay=True, writable=True),
    help=(
        "Directory to write the output files to. The files written will be named "
        "unmapped.reads, contaminant.reads, and keep.reads"
    ),
    default=".",
    show_default=True,
)
@click.option(
    "--ignore-secondary/--include-secondary",
    help="Ignore organism assignments for secondary alignments?",
    default=True,
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    samfile: str, metadata: str, outdir: str, ignore_secondary: bool, verbose: bool,
):
    """This scripts classifies records in an alignment to a contamination database, with
    the help of a metadata file mapping reference names to whether they are
    contamination or not. It produces three files with a read identifier for each line:
      - unmapped reads
      - contaminated reads
      - reads to keep
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    classifier = Classifier(metadata)
    all_read_ids: Set[str] = set()
    unmapped_reads: Set[str] = set()
    contaminant_reads: Set[str] = set()
    keep_reads: Set[str] = set()

    logging.info("Classifying records in alignment file...")
    with pysam.AlignmentFile(samfile) as bam:
        for record in bam:
            read_id = record.query_name
            if read_id is None:
                logging.warning(f"Got a record with no query name\n{str(record)}")
                continue

            all_read_ids.add(read_id)
            if record.is_secondary and ignore_secondary:
                logging.debug(f"{read_id} has secondary alignment. Skipping...")
                continue

            classification = classifier.classify(record)
            if classification is Classification.Unmaped:
                unmapped_reads.add(read_id)
            elif classification is Classification.Keep:
                keep_reads.add(read_id)
            elif classification.Contaminant:
                contaminant_reads.add(read_id)
            else:
                raise NotImplementedError(
                    f"Don't know how to handle classification: {classification}"
                )

    # if any read in the pair is a "keeper" remove it from contaminants
    contaminant_reads -= keep_reads
    # reads are only unmapped if both are unmapped
    unmapped_reads -= keep_reads.union(contaminant_reads)

    assert all_read_ids == keep_reads.union(unmapped_reads, contaminant_reads)

    logging.info(f"{len(keep_reads)} reads are to be kept")
    logging.info(f"{len(contaminant_reads)} reads are contaminants")
    logging.info(f"{len(unmapped_reads)} reads are unmapped")
    logging.info("Writing output files...")

    keep_file = outdir / "keep.reads"
    keep_file.write_text("\n".join(keep_reads))
    logging.info(f"Read identifiers to keep written to {keep_file}")

    contaminant_file = outdir / "contaminant.reads"
    contaminant_file.write_text("\n".join(contaminant_reads))
    logging.info(f"Contaminant read identifiers written to {contaminant_file}")

    unmapped_file = outdir / "unmapped.reads"
    unmapped_file.write_text("\n".join(unmapped_reads))
    logging.info(f"Unmapped read identifiers written to {unmapped_file}")


if __name__ == "__main__":
    main()
