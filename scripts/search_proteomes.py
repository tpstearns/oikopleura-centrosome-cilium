#!/usr/bin/env python3
"""hmmsearch every per-gene profile against one proteome; write a hit table with profile coverage.

usage: search_proteomes.py --profiles data/profiles --proteome data/proteomes/<key>.faa \
                           --out results/hmmsearch/<key>.tsv --params config/params.yaml --threads 8
Output columns: gene, target, evalue, bits, best_domain_ievalue, profile_coverage, target_coverage, n_domains
Hits are written unfiltered below E=10 so thresholds can be changed in params without rerunning.
"""
import argparse
import sys
from pathlib import Path

import pyhmmer

sys.path.insert(0, str(Path(__file__).parent))
from common import log, text, write_tsv  # noqa: E402

FIELDS = ["gene", "target", "evalue", "bits", "best_domain_ievalue", "profile_coverage", "target_coverage",
          "n_domains", "profile_M", "profile_desc"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    alphabet = pyhmmer.easel.Alphabet.amino()
    hmms = []
    for path in sorted(Path(args.profiles).glob("*.hmm")):
        with pyhmmer.plan7.HMMFile(path) as hf:
            hmms.extend(list(hf))
    with pyhmmer.easel.SequenceFile(args.proteome, digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    lengths = {text(s.name): len(s) for s in seqs}
    log(f"{len(hmms)} profiles vs {len(seqs)} sequences")

    rows = []
    for hits in pyhmmer.hmmsearch(hmms, seqs, cpus=args.threads, E=10.0):
        gene = text(hits.query.name)
        M = hits.query.M
        for hit in hits:
            covered = set()
            tcovered = set()
            best_i = None
            for dom in hit.domains:
                if dom.i_evalue is not None:
                    best_i = dom.i_evalue if best_i is None else min(best_i, dom.i_evalue)
                a = dom.alignment
                covered.update(range(a.hmm_from, a.hmm_to + 1))
                tcovered.update(range(a.target_from, a.target_to + 1))
            tname = text(hit.name)
            rows.append({"gene": gene, "target": tname, "evalue": f"{hit.evalue:.3g}", "bits": f"{hit.score:.1f}",
                         "best_domain_ievalue": f"{best_i:.3g}" if best_i is not None else "",
                         "profile_coverage": f"{len(covered) / M:.3f}",
                         "target_coverage": f"{len(tcovered) / lengths[tname]:.3f}",
                         "n_domains": len(hit.domains), "profile_M": M,
                         "profile_desc": text(hits.query.description)})
    write_tsv(args.out, rows, FIELDS)
    log(f"wrote {len(rows)} hits to {args.out}")


if __name__ == "__main__":
    main()
