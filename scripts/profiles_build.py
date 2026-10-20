#!/usr/bin/env python3
"""Build one HMM per gene from its seed sequences (MAFFT alignment + hmmbuild via pyhmmer).

usage: profiles_build.py data/seeds/<gene>.faa data/profiles/<gene>.hmm [--threads N]
A single-sequence seed file yields a single-sequence profile (hmmbuild from one sequence), which is noted
in the HMM description so downstream reports can flag weaker profiles.
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pyhmmer

sys.path.insert(0, str(Path(__file__).parent))
from common import die, log, set_text  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("seeds")
    ap.add_argument("out")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    gene = Path(args.seeds).stem
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(args.seeds, digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    if len(seqs) == 0:
        die(f"{args.seeds}: no sequences")
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    if len(seqs) == 1:
        hmm, _, _ = builder.build(seqs[0], background)
        set_text(hmm, "name", gene)
        set_text(hmm, "description", "single-sequence profile")
    else:
        with tempfile.NamedTemporaryFile(suffix=".sto", delete=False) as tmp:
            aln_path = tmp.name
        with open(aln_path, "w") as fh:
            subprocess.run(["mafft", "--auto", "--thread", str(args.threads), "--quiet", args.seeds],
                           check=True, stdout=fh)
        with pyhmmer.easel.MSAFile(aln_path, digital=True, alphabet=alphabet) as mf:
            msa = mf.read()
        set_text(msa, "name", gene)
        hmm, _, _ = builder.build_msa(msa, background)
        set_text(hmm, "description", f"profile from {len(seqs)} seeds")
        Path(aln_path).unlink()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wb") as fh:
        hmm.write(fh)
    log(f"[{gene}] profile from {len(seqs)} sequence(s), M={hmm.M}")


if __name__ == "__main__":
    main()
