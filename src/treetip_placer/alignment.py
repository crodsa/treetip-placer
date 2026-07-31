"""Fixed-column FASTA alignment adapter for SNP and concatenated-gene inputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .models import Profile
from .profiles import normalise_allele


SEQUENCE_MISSING = frozenset({"-", ".", "?", "N", "X"})


@dataclass(frozen=True)
class AlignmentTable:
    """An aligned character matrix represented through the generic profile engine."""

    profiles: tuple[Profile, ...]
    site_count: int


def read_fasta_alignment(path: str | Path) -> AlignmentTable:
    """Read a fixed-column FASTA alignment without altering any reference columns."""

    alignment_path = Path(path)
    records: list[tuple[str, str]] = []
    current_name: str | None = None
    chunks: list[str] = []
    seen_names: set[str] = set()

    def finish_record() -> None:
        if current_name is None:
            return
        sequence = "".join(chunks).replace(" ", "").upper()
        if not sequence:
            raise ValueError(f"FASTA record '{current_name}' has no sequence: {alignment_path}")
        records.append((current_name, sequence))

    with alignment_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {alignment_path}:{line_number}")
                current_name = header.split(maxsplit=1)[0]
                if current_name in seen_names:
                    raise ValueError(f"Duplicate FASTA identifier '{current_name}'")
                seen_names.add(current_name)
                chunks = []
            else:
                if current_name is None:
                    raise ValueError(
                        f"Sequence content appears before a FASTA header at {alignment_path}:{line_number}"
                    )
                chunks.append(line)
    finish_record()

    if not records:
        raise ValueError(f"FASTA alignment has no records: {alignment_path}")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise ValueError(
            f"FASTA input is not a fixed-column alignment; observed lengths: {sorted(lengths)}"
        )
    site_count = lengths.pop()
    profiles = tuple(
        Profile(
            sample_id=name,
            alleles=tuple(
                None if character in SEQUENCE_MISSING else normalise_allele(character)
                for character in sequence
            ),
        )
        for name, sequence in records
    )
    return AlignmentTable(profiles=profiles, site_count=site_count)


def require_matching_alignment_widths(
    reference: AlignmentTable, queries: AlignmentTable
) -> None:
    """Reject unaligned queries rather than silently re-aligning them."""

    if reference.site_count != queries.site_count:
        raise ValueError(
            "Reference and query alignments have different column counts "
            f"({reference.site_count} and {queries.site_count}). Align or map queries "
            "to the immutable reference coordinate system before placement."
        )
