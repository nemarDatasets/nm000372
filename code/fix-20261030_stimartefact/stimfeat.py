"""Stimulation-artefact features for nm000372 / nm000380 (next-steps lane, Oct 2026). Runs in a Voyager Job.

For every *_ieeg.vhdr of a NEMAR clone: fetch the BrainVision files by annex key from the public NEMAR S3 bucket
(size + sha256 checked against the key, s3fetch.py), read them with MNE, keep the intracranial channels of
channels.tsv (type ECOG/SEEG/DBS, status not bad), and save two features to <out>/<stem>.npz:

  pulse   per sample n (native rate): 90th percentile across channels of |x[n+1]-x[n]| / robust_sd(diff of that
          channel), robust_sd = 1.4826 * median absolute deviation. Stimulation pulses are steps far above the
          channel's normal sample-to-sample change in most channels at once.
  spec    log10 power, median across channels, of a short-time Fourier transform (1 s Hann window, 0.1 s hop,
          1 Hz bins) of the signal resampled to about 1 kHz (scipy resample_poly, integer factor).

Nothing is written to the repository. Usage: stimfeat.py <repo> <out> [glob]
"""
import glob, json, os, sys, time
import numpy as np
import mne
from scipy.signal import resample_poly, stft

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s3fetch

mne.set_log_level("ERROR")
repo, out = sys.argv[1], sys.argv[2]
pat = sys.argv[3] if len(sys.argv) > 3 else "sub-*/**/*_ieeg.vhdr"
os.makedirs(out, exist_ok=True)
TMP = "/work/data"
vh = sorted(os.path.relpath(p, repo) for p in glob.glob(f"{repo}/{pat}", recursive=True))
print("runs", len(vh), flush=True)


def chan_table(rel_vhdr):
    p = os.path.join(repo, rel_vhdr.replace("_ieeg.vhdr", "_channels.tsv"))
    rows = [l.rstrip("\n").split("\t") for l in open(p, encoding="utf-8")]
    h = rows[0]
    return [dict(zip(h, r)) for r in rows[1:]]


for rel in vh:
    stem = os.path.basename(rel)[:-len("_ieeg.vhdr")]
    dst = f"{out}/{stem}.npz"
    if os.path.exists(dst):
        continue
    t0 = time.time()
    b = rel[:-len(".vhdr")]
    files = [b + ".vhdr", b + ".vmrk", b + ".eeg"]
    loc = s3fetch.materialize(repo, files, TMP, ds=os.path.basename(os.path.normpath(repo)))
    raw = mne.io.read_raw_brainvision(loc[b + ".vhdr"], preload=True)
    fs = raw.info["sfreq"]
    ct = chan_table(rel)
    keep = [r["name"] for r in ct if r.get("type", "").upper() in ("ECOG", "SEEG", "DBS")
            and r.get("status", "good") != "bad" and r["name"] in raw.ch_names]
    used = "ECOG/SEEG/DBS"
    if not keep:  # recordings whose contacts are typed MISC (type unstated in the release): use those
        keep = [r["name"] for r in ct if r.get("type", "").upper() == "MISC"
                and r.get("status", "good") != "bad" and r["name"] in raw.ch_names]
        used = "MISC"
    n = raw.n_times
    # pulse score, channel by channel (memory: one channel at a time for the robust scale; accumulate z rows)
    z = np.zeros((len(keep), n - 1), dtype=np.float32)
    good = np.zeros(len(keep), bool)
    for i, ch in enumerate(keep):
        x = raw.get_data(picks=[ch])[0]
        d = np.diff(x)
        med = np.median(d)
        s = 1.4826 * np.median(np.abs(d - med))
        if not np.isfinite(s) or s <= 0:
            continue
        z[i] = np.abs(d - med) / s
        good[i] = True
    z = z[good]
    k = int(np.floor(0.9 * (z.shape[0] - 1)))
    pulse = np.empty(n - 1, dtype=np.float32)
    step = 500000
    for a in range(0, n - 1, step):
        blk = z[:, a:a + step]
        pulse[a:a + step] = np.partition(blk, k, axis=0)[k]
    del z
    # spectral feature
    q = int(fs // 1000) if fs >= 2000 else 1
    fsd = fs / q
    nper = int(round(fsd))
    hop = max(1, int(round(fsd / 10)))
    logp = []
    for ch in np.array(keep)[good]:
        x = raw.get_data(picks=[ch])[0]
        if q > 1:
            x = resample_poly(x, 1, q)
        f, tt, Z = stft(x, fs=fsd, window="hann", nperseg=nper, noverlap=nper - hop, boundary=None, padded=False)
        logp.append(np.log10(np.abs(Z) ** 2 + 1e-30).astype(np.float32))
    logp = np.stack(logp)
    spec_med = np.median(logp, axis=0)
    spec_p90 = np.percentile(logp, 90, axis=0).astype(np.float32)
    np.savez_compressed(dst, pulse=pulse, fs=fs, n=n, spec_med=spec_med.astype(np.float32), spec_p90=spec_p90,
                        spec_f=f.astype(np.float32), spec_t=tt.astype(np.float32), fsd=fsd,
                        channels=np.array(keep)[good], n_channels=int(good.sum()), chan_types_used=used)
    for p in loc.values():
        os.remove(p)
    print(stem, f"fs={fs} n={n} ch={int(good.sum())} q={q} {time.time()-t0:.0f}s", flush=True)
open(f"{out}/DONE", "w").write("done\n")
