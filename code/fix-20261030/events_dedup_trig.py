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
   SPES: the TRIG channel changes level only twice: it is at its higher level (251) at 0-11.570 s and 311.532-320.000 s
   and at its lower level (219) in between. Correction of 2 Nov 2026 (after J. Gadonneix's review comment): the lower
   level is the stimulator-on gate. In the HFS run the six 5.05-s trains are stretches at that same lower level; in the
   SPES run its one stretch (11.570-311.532 s) spans all 20 artefact-detected single pulses (the first 0.7 ms after
   it opens, the last 14.96 s before it closes). The two higher-level stretches written as `trig_pulse` on 30 Oct were
   the idle state. The SPES row is now that one lower-level stretch, trial_type `spes_gate` (it does not mark the
   individual pulses; see derivatives/stimartefact for those).
Usage: python events_dedup_trig.py <bids_root>
"""
import csv, math, os, sys

root = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
SF = 4096.0
for task, label in (("hfs", "hfs_train"), ("spes", "spes_gate")):
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
    pulses = []
    for p in csv.DictReader(open(os.path.join(here, f"nm000372_{task}_trig_pulses.tsv")), delimiter="\t"):
        s = int(p["sample"])
        assert abs(float(p["onset"]) * SF - s) < 0.05  # onset printed with 6 decimals
        pulses.append((s, round(float(p["duration"]) * SF), float(p["levels"])))  # duration in samples
    LOW = 0.000219  # lower TRIG level = stimulator on (both runs)
    if all(abs(lv - LOW) < 1e-9 for _, _, lv in pulses):
        gates = [(s, nd) for s, nd, _ in pulses]  # HFS: the decoded stretches are at the lower level already
    else:  # SPES: decoded stretches are at the higher (idle) level; the gate is what lies between them
        n_total = round(320.000244140625 * SF)  # samples in the SPES run (RecordingDuration x SamplingFrequency)
        idle = sorted((s, s + nd) for s, nd, _ in pulses)
        edges = [0] + [x for iv in idle for x in iv] + [n_total]
        gates = [(a, b - a) for a, b in zip(edges[0::2], edges[1::2]) if b > a]
    for s, nd in gates:
        fl = lambda x: f"{math.floor(x * 1e6) / 1e6:.6f}"  # truncate so onset + duration never exceeds the data
        out.append({"onset": fl(s / SF), "duration": fl(nd / SF), "trial_type": label,
                    "value": "n/a", "sample": str(s), "targets": "n/a", "source": "TRIG"})
    out.sort(key=lambda r: (float(r["onset"]), r["source"]))
    with open(ev, "w", newline="") as f:
        w = csv.DictWriter(f, header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(task, "rows", len(rows), "->", len(out), "merged", merged)
