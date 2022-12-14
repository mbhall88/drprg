import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum

STOP = "*"
MISSENSE = "X"
PROT = "PROT"
DNA = "DNA"


def split_var_name(name: str) -> tuple[str, int, str]:
    items = re.match(r"([A-Z]+)([-0-9]+)([A-Z/\*]+)", name, re.I).groups()
    return items[0], int(items[1]), items[2]


@dataclass(frozen=True)
class PanelRecord:
    gene: str
    ref: str
    pos: int
    alt: str
    residue: str
    drugs: list[str]

    @staticmethod
    def from_row(row: str, delim: str = "\t") -> "PanelRecord":
        fields = row.rstrip().split(delim)
        gene = fields[0]
        drugs = fields[3].split(",")
        residue = fields[2]
        ref, pos, alt = split_var_name(fields[1])
        return PanelRecord(gene, ref, pos, alt, residue, drugs)


class RuleType(Enum):
    Missense = "missense"
    Frameshift = "frameshift"
    Nonsense = "nonsense"
    Absence = "absence"


@dataclass(frozen=True)
class Rule:
    rule_type: RuleType
    gene: str
    drugs: list[str]
    start: int | None = None  # 1-based inclusive
    end: int | None = None  # 1-based inclusive

    def __contains__(self, item: PanelRecord) -> bool:
        if item.gene != self.gene:
            return False
        start = 1 if self.start is None else self.start
        end = float("inf") if self.end is None else self.end
        if item.pos < start or item.pos > end:
            return False
        match self.rule_type:
            case RuleType.Nonsense:
                return item.alt == STOP and item.residue == PROT
            case RuleType.Missense:
                return item.residue == PROT and item.ref != item.alt
            case RuleType.Frameshift:
                len_diff = abs(len(item.ref) - len(item.alt))
                return item.residue == DNA and len_diff % 3 != 0
            case _:
                return False

    @staticmethod
    def from_row(row: str, delim: str = ",") -> "Rule":
        fields = row.rstrip().split(delim)
        rule_type = RuleType(fields[0])
        gene = fields[1]
        start = 1 if not fields[2] else int(fields[2])
        end = float("inf") if not fields[3] else int(fields[3])
        drugs = fields[4].split(";")
        return Rule(rule_type, gene, drugs, start, end)


class ExpertRules:
    def __init__(self, data: dict[str, list[Rule]]):
        self._data = data

    def __contains__(self, item: PanelRecord) -> bool:
        rules = self._data.get(item.gene, [])
        return any(item in r for r in rules)


def main():
    input_panel = sys.argv[1]
    input_rules = sys.argv[2]
    output_panel = sys.argv[3]
    output_rules = sys.argv[4]
    filter_x_amino_mutations = True

    fp_out_rules = open(output_rules, "w")
    with open(input_rules) as fp_in:
        d = defaultdict(list)
        for row in fp_in:
            rule = Rule.from_row(row)
            d[rule.gene].append(rule)
            fp_out_rules.write(row)

    expert_rules = ExpertRules(d)
    records_kept = 0
    records_filtered = 0
    extra_rules = 0
    with open(output_panel, "w") as fp_out, open(input_panel) as fp_in:
        for row in fp_in:
            record = PanelRecord.from_row(row)
            if record not in expert_rules:
                if (
                    filter_x_amino_mutations
                    and record.alt == MISSENSE
                    and record.residue == PROT
                ):
                    new_row = f"{RuleType.Missense.value},{record.gene},{record.pos},{record.pos},{';'.join(record.drugs)}"
                    print(new_row, file=fp_out_rules)
                    records_filtered += 1
                    extra_rules += 1
                else:
                    fp_out.write(row)
                    records_kept += 1
            else:
                records_filtered += 1

    fp_out_rules.close()
    print(
        f"Kept {records_kept} panel records, filtered {records_filtered}. Generated {extra_rules} extra rules",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
