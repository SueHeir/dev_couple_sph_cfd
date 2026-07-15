#!/usr/bin/env python3
"""Maintain identities of the public PSI context sources named by this repository.

This is citation maintenance, not model validation. It deliberately does not
download observations, calculate an error, or produce a scientific verdict. The
cited pages are context only, not a case-matched held-out dataset. It must not
be added to the executable validation manifest.
"""

from __future__ import annotations

import html
import sys
from dataclasses import dataclass
from urllib.error import URLError
from urllib.request import Request, urlopen


@dataclass(frozen=True)
class PublicRecord:
    name: str
    url: str
    required_text: tuple[str, ...]


# Identity markers only: no measured values or acceptance bands appear here.
RECORDS = (
    PublicRecord(
        "NASA NTRS 20210016650",
        "https://ntrs.nasa.gov/citations/20210016650",
        (
            "20210016650",
            "Overview of Plume-Surface Interaction Data from Subscale Inert Gas Testing",
            "Chad J. Eberhart",
        ),
    ),
    PublicRecord(
        "NASA Langley PSI announcement",
        "https://www.nasa.gov/general/what-a-blast-nasa-langley-begins-plume-surface-interaction-tests/",
        ("Plume-Surface Interaction",),
    ),
)


def fetch_text(url: str) -> str:
    request = Request(url, headers={"User-Agent": "dev_couple_sph_cfd-reference-audit/1"})
    with urlopen(request, timeout=30) as response:
        if response.status != 200:
            raise RuntimeError(f"HTTP {response.status}")
        return html.unescape(response.read().decode("utf-8", errors="replace"))


def main() -> int:
    failures = []
    for record in RECORDS:
        try:
            document = fetch_text(record.url)
        except (URLError, OSError, RuntimeError) as exc:
            failures.append(f"{record.name}: unavailable ({exc})")
            continue
        missing = [text for text in record.required_text if text.casefold() not in document.casefold()]
        if missing:
            failures.append(f"{record.name}: identity marker missing: {missing!r}")
        else:
            print(f"verified public-record identity: {record.name}")
    if failures:
        print("external-reference citation audit FAILED; re-review before relying on these records:")
        print("\n".join(f"- {failure}" for failure in failures))
        return 1
    print("citation audit completed; this is not scientific validation or plume acceptance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
