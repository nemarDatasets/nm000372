"""Stimulation onsets detected from the stimulation artefact (nm000372, nm000380). Not from a source.

Input: the per-run feature files of stimfeat.py (<feat>/<stem>.npz: 'pulse' = per-sample step score, 'spec_med' =
median-over-channels log10 power spectrogram). Output: a derivative tree <out> (BIDS derivative, events tables).

Rules (fixed in this file; thresholds chosen from the score distributions of these recordings, see README):
  PULSE_THR    a stimulation pulse is a sample step whose score (90th percentile across contacts of
               |step| / robust SD of that contact's steps) exceeds 50. Typical background: median 0.7-1.4,
               99.9th percentile below 7 in recordings without stimulation pulses.
  single pulses (nm000372 task-spes, nm000380 task-sps): pulse candidates at least REFRACT s apart (first sample of
               each cluster above threshold); a candidate is kept only if it belongs to a chain of >= MIN_CHAIN pulses
               whose consecutive intervals equal the protocol interval (15 s SPES, nm000372 article; 1 s SPS,
               nm000380 article) within +-TOL s. Candidates outside such chains are listed as rejected.
  trains at 43.3 Hz (nm000372 task-hfs): pulse candidates >= 10 ms apart, grouped while gaps are <= 0.1 s; a group is
               a train if it has >= 20 pulses and a median interval of 23.1 +- 2 ms (43.3 Hz).
  100 Hz blocks (nm000380 task-stim): comb score per 0.1-s frame = mean over 100, 200 and 400 Hz of
               log10 P(f) - mean log10 P(f +- 6, 8, 10 Hz) on the median-over-contacts spectrogram (1-s windows).
               Frames with score > 0.5 (power at the stimulation lines > 3.2x the neighbouring bins) are stimulation;
               gaps <= 5 s are bridged; blocks shorter than 20 s are ignored. Block onset/offset are then refined to
               the first/last pulse step (score > 20, >= 4 ms apart) that has >= 8 further steps within the next/previous
               100 ms, searched within 1.5 s of the frame edge; if none is found the frame edge +- 0.5 s is used
               (onset_method = spectral, precision about +-0.5 s).
Usage: stimdetect.py <nm> <feat_dir> <repo> <out_derivative_dir>
"""
import csv, glob, json, os, sys
import numpy as np

NM, FEAT, REPO, OUT = sys.argv[1:5]
PULSE_THR = 50.0
STEP_THR = 20.0
COMB_THR = 0.5
summary = {}


def clusters(score, fs, thr, refract):
    idx = np.flatnonzero(score > thr)
    out = []
    last = -10 ** 12
    for i in idx:
        if i - last > refract * fs:
            out.append([i, float(score[i])])
        else:
            out[-1][1] = max(out[-1][1], float(score[i]))
        last = i
    return [(int(i) + 1, p) for i, p in out]  # +1: first sample after the step


def chains(times, period, tol, min_chain):
    keep = set()
    n = len(times)
    best = []
    i = 0
    while i < n:
        ch = [i]
        j = i
        while True:
            nxt = [k for k in range(j + 1, n) if abs(times[k] - times[j] - period) <= tol]
            if not nxt:
                break
            j = nxt[0]
            ch.append(j)
        if len(ch) >= min_chain:
            keep.update(ch)
            i = ch[-1] + 1
        else:
            i += 1
    return keep


def comb(S, F, hs=(100, 200, 400)):
    out = []
    for h in hs:
        if h + 10 > F[-1]:
            continue
        i = np.argmin(abs(F - h))
        nb = [np.argmin(abs(F - (h + d))) for d in (-10, -8, -6, 6, 8, 10)]
        out.append(S[i] - S[nb].mean(0))
    return np.mean(out, 0)


def refine(pulse, fs, t_edge, forward):
    lo = max(0, int((t_edge - 1.5) * fs)); hi = min(len(pulse), int((t_edge + 1.5) * fs))
    c = [i for i, _ in clusters(pulse[lo:hi], fs, STEP_THR, 0.004)]
    c = [lo + i for i in c]
    w = 0.1 * fs
    seq = c if forward else c[::-1]
    for i in seq:
        if forward:
            k = sum(1 for j in c if i < j <= i + w)
        else:
            k = sum(1 for j in c if i - w <= j < i)
        if k >= 8:
            return i
    return None


def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")


COLS = ["onset", "duration", "trial_type", "sample", "peak_score", "n_pulses", "method"]
for fpath in sorted(glob.glob(f"{FEAT}/*.npz")):
    stem = os.path.basename(fpath)[:-4]
    z = np.load(fpath)
    fs = float(z["fs"]); pulse = z["pulse"]
    sub = stem.split("_")[0]
    rows, rejected = [], []
    info = {"sfreq": fs, "n_contacts": int(z["n_channels"]), "duration_s": round(len(pulse) / fs, 3)}
    if "chan_types_used" in z.files:
        info["contact_types_used"] = str(z["chan_types_used"])
    if "task-spes" in stem or "task-sps" in stem:
        period, tol, refr, mc, lab = (15.0, 0.5, 1.0, 5, "spes_pulse") if "task-spes" in stem else (1.0, 0.1, 0.3, 5, "sps_pulse")
        cand = clusters(pulse, fs, PULSE_THR, refr)
        t = [s / fs for s, _ in cand]
        keep = chains(t, period, tol, mc)
        for k, (s, p) in enumerate(cand):
            r = {"onset": f"{s/fs:.6f}", "duration": 0, "trial_type": lab, "sample": s, "peak_score": f"{p:.1f}",
                 "n_pulses": 1, "method": "pulse_step"}
            (rows if k in keep else rejected).append(r)
        # trains (consecutive kept pulses at the protocol interval)
        kept = sorted(keep)
        tr = []
        for k in kept:
            if tr and abs(t[k] - t[tr[-1][-1]] - period) <= tol:
                tr[-1].append(k)
            else:
                tr.append([k])
        if lab == "sps_pulse":
            for g in tr:
                s0 = cand[g[0]][0]; s1 = cand[g[-1]][0]
                rows.append({"onset": f"{s0/fs:.6f}", "duration": f"{(s1-s0)/fs:.6f}", "trial_type": "sps_train",
                             "sample": s0, "peak_score": f"{max(cand[k][1] for k in g):.1f}", "n_pulses": len(g),
                             "method": "pulse_chain"})
        info.update(n_candidates=len(cand), n_pulses=len(keep), n_trains=len(tr), train_sizes=[len(g) for g in tr],
                    rejected_onsets_s=[float(r["onset"]) for r in rejected])
    elif "task-hfs" in stem:
        cand = clusters(pulse, fs, PULSE_THR, 0.010)
        groups = []
        for s, p in cand:
            if groups and s - groups[-1][-1][0] <= 0.1 * fs:
                groups[-1].append((s, p))
            else:
                groups.append([(s, p)])
        ntr = 0
        for g in groups:
            ipi = np.median(np.diff([s for s, _ in g])) / fs if len(g) > 1 else 0
            if len(g) >= 20 and abs(ipi - 1 / 43.3) <= 0.002:
                ntr += 1
                s0, s1 = g[0][0], g[-1][0]
                rows.append({"onset": f"{s0/fs:.6f}", "duration": f"{(s1-s0)/fs + ipi:.6f}", "trial_type": "hfs_train",
                             "sample": s0, "peak_score": f"{max(p for _, p in g):.1f}", "n_pulses": len(g),
                             "method": "pulse_train"})
            else:
                rejected.append({"onset": f"{g[0][0]/fs:.6f}", "n": len(g)})
        info.update(n_trains=ntr, train_pulses=[r["n_pulses"] for r in rows],
                    rejected_groups=len(rejected))
    elif "task-stim" in stem:
        F = z["spec_f"]; T = z["spec_t"]; sc = comb(z["spec_med"], F)
        on = sc > COMB_THR
        d = np.diff(np.r_[0, on.astype(int), 0]); st = list(np.flatnonzero(d == 1)); en = list(np.flatnonzero(d == -1) - 1)
        segs = []
        for a, b in zip(st, en):
            if segs and T[a] - T[segs[-1][1]] <= 5.0:
                segs[-1][1] = b
            else:
                segs.append([a, b])
        segs = [s for s in segs if T[s[1]] - T[s[0]] >= 20.0]
        for a, b in segs:
            t0 = float(T[a]) + 0.5; t1 = float(T[b]) - 0.5  # frame centres +- half window
            r0 = refine(pulse, fs, float(T[a]), True); r1 = refine(pulse, fs, float(T[b]), False)
            m = "pulse_refined" if (r0 is not None and r1 is not None) else "spectral"
            on_s = r0 / fs if r0 is not None else t0
            off_s = r1 / fs if r1 is not None else t1
            rows.append({"onset": f"{on_s:.6f}", "duration": f"{off_s - on_s:.6f}", "trial_type": "stim_block_100Hz",
                         "sample": int(round(on_s * fs)), "peak_score": f"{float(sc[a:b+1].max()):.2f}", "n_pulses": "n/a",
                         "method": m})
        info.update(n_blocks=len(segs), blocks=[(r["onset"], r["duration"], r["method"]) for r in rows],
                    comb_p50=float(np.median(sc)), comb_p95=float(np.percentile(sc, 95)))
    else:
        continue
    rel = glob.glob(f"{REPO}/{sub}/**/{stem}_ieeg.vhdr", recursive=True)
    assert len(rel) == 1, (stem, rel)
    d = os.path.dirname(os.path.relpath(rel[0], REPO))
    rows.sort(key=lambda r: (float(r["onset"]), r["trial_type"]))
    write_tsv(f"{OUT}/{d}/{stem}_desc-stimartefact_events.tsv", rows, COLS)
    summary[stem] = info
json.dump(summary, open(f"{OUT}/detection_summary.json", "w"), indent=1, default=str)
print(json.dumps(summary, indent=0, default=str))
