"""fix-20261030 for nm000372 events.tsv.

1. HFS run: the authors' analysis marker file (<recording>.ades.mrk) repeats the 6 acquisition markers
   (<recording>.mrk) with identical onset, duration and label. Each such pair becomes one row whose `source` is
   "acquisition_mrk;analysis_ades_mrk". No marker is lost.
2. Both runs: rows decoded from the recorded TRIG channel are added (source "TRIG"). The TRIG channel takes two levels;
   a pulse is a stretch where it leaves its most frequent (resting) level. Pulses were decoded in a Voyager Job from the
   .eeg/.vhdr fetched by git-annex key (sha256 verified) with MNE 1.x read_raw_brainvision; the decoded lists are
   *_trig_pulses.tsv next to this script (onset s, duration s, 0-based sample, levels).
   HFS: 6 pulses of 5.046 s, matching the protocol's 5-s trains (README); 5 coincide with the marked trains and the
   first (23.675 s) has no marker. trial_type `hfs_train`.
   SPES: the TRIG channel changes level only twice (0-11.570 s and 311.532-320.000 s); it does not mark the 20 single
   pulses. trial_type `trig_pulse` (meaning not documented in the release).
Usage: python events_dedup_trig.py <bids_root>
"""
import csv, math, os, sys

root = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
SF = 4096.0
for task, label in (("hfs", "hfs_train"), ("spes", "trig_pulse")):
    ev = os.path.join(root, f"sub-01/ieeg/sub-01_task-{task}_run-1_events.tsv")
    rows = list(csv.DictReader(open(ev, newline=""), delimiter="\t"))
    header = list(rows[0].keys())
    key = lambda r: (r["onset"], r["duration"], r["trial_type"], r["value"], r["sample"], r["targets"])
    acq = {key(r) for r in rows if r["source"] == "acquisition_mrk"}
    out, merged = [], 0
    for r in rows:
        if r["source"] == "analysis_ades_mrk" and key(r) in acq:
            merged += 1
            continue
        if r["source"] == "acquisition_mrk" and any(key(x) == key(r) and x["source"] == "analysis_ades_mrk" for x in rows):
            r = dict(r, source="acquisition_mrk;analysis_ades_mrk")
        out.append(r)
    for p in csv.DictReader(open(os.path.join(here, f"nm000372_{task}_trig_pulses.tsv")), delimiter="\t"):
        s = int(p["sample"])
        assert abs(float(p["onset"]) * SF - s) < 0.05  # onset printed with 6 decimals
        nd = round(float(p["duration"]) * SF)  # duration in samples (the job printed seconds to 6 decimals)
        fl = lambda x: f"{math.floor(x * 1e6) / 1e6:.6f}"  # truncate so onset + duration never exceeds the data
        out.append({"onset": fl(s / SF), "duration": fl(nd / SF), "trial_type": label,
                    "value": "n/a", "sample": str(s), "targets": "n/a", "source": "TRIG"})
    out.sort(key=lambda r: (float(r["onset"]), r["source"]))
    with open(ev, "w", newline="") as f:
        w = csv.DictWriter(f, header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(task, "rows", len(rows), "->", len(out), "merged", merged)
