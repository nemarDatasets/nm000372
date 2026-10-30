"""Write the sidecars, dataset_description.json and README of derivatives/stimartefact (nm000372, nm000380) from the
output of stimdetect.py (detection_summary.json + events tables). Usage: build_derivative.py <nm> <deriv_dir>"""
import glob, json, os, sys

NM, D = sys.argv[1:3]
S = json.load(open(f"{D}/detection_summary.json"))
COLS = {
    "onset": {"Description": "Detected onset (s from the first sample of the recording): the first sample after the "
                             "first stimulation-artefact step of the pulse, train or block.", "Units": "s"},
    "duration": {"Description": "0 for single pulses; for trains and blocks the time from the first to the last "
                                "detected pulse (trains: plus one pulse interval).", "Units": "s"},
    "trial_type": {"Description": "Kind of detected stimulation.", "Levels": {}},
    "sample": {"Description": "onset x sampling frequency (0-based sample index in the data file)."},
    "peak_score": {"Description": "Largest detection score of the event: for pulses and trains the step score (90th "
                                  "percentile across contacts of |sample step| / robust SD of that contact's steps); for "
                                  "100 Hz blocks the largest comb score (log10 power ratio, see README)."},
    "n_pulses": {"Description": "Number of detected pulses in the event (n/a for 100 Hz blocks, where single pulses "
                                "are not counted)."},
    "method": {"Description": "How onset and duration were obtained.", "Levels": {
        "pulse_step": "single pulse: step score above 50",
        "pulse_chain": "train of single pulses: consecutive detected pulses at the protocol interval",
        "pulse_train": "43.3 Hz train: >= 20 pulses at 23.1 +- 2 ms intervals",
        "pulse_refined": "100 Hz block from the spectral comb score; onset/offset refined to the first/last pulse step",
        "spectral": "100 Hz block from the spectral comb score only (onset/offset precision about +-0.5 s)"}},
}
LEV = {"spes_pulse": "single-pulse electrical stimulation (SPES) pulse, detected from the artefact",
       "hfs_train": "43.3 Hz high-frequency stimulation train, detected from the artefact",
       "sps_pulse": "single pulse of single-pulse stimulation (SPS), detected from the artefact",
       "sps_train": "train of SPS pulses at 1 s intervals (first to last pulse), detected from the artefact",
       "stim_block_100Hz": "block of continuous 100 Hz stimulation, detected from the artefact"}
n_files = 0
for tsv in sorted(glob.glob(f"{D}/sub-*/**/*_desc-stimartefact_events.tsv", recursive=True)):
    types = sorted({l.split("\t")[2] for l in open(tsv).read().splitlines()[1:]})
    c = json.loads(json.dumps(COLS))
    c["trial_type"]["Levels"] = {t: LEV[t] for t in types}
    used = {r.split("\t")[-1] for r in open(tsv).read().splitlines()[1:]}
    c["method"]["Levels"] = {k: v for k, v in c["method"]["Levels"].items() if k in used}
    c["Description"] = ("Stimulation onsets DETECTED FROM THE STIMULATION ARTEFACT in the recorded signals. Not from a "
                        "source: the release does not annotate these events. Method, thresholds and counts: README of "
                        "this derivative; code: code/fix-20261030_stimartefact/.")
    json.dump(c, open(tsv[:-4] + ".json", "w"), indent=2, ensure_ascii=False)
    n_files += 1
dd = {"Name": f"{NM}: stimulation onsets detected from the stimulation artefact (not from a source)",
      "BIDSVersion": "1.9.0", "DatasetType": "derivative",
      "GeneratedBy": [{"Name": "stimfeat.py + stimdetect.py (NEMAR curation, Oct 2026)",
                       "Description": "Artefact-based detection of stimulation pulses, trains and 100 Hz blocks; "
                                      "code in code/fix-20261030_stimartefact/ of the raw dataset.",
                       "CodeURL": f"https://github.com/nemarDatasets/{NM}/tree/main/code/fix-20261030_stimartefact"}],
      "SourceDatasets": [{"DOI": f"10.82901/nemar.{NM}"}]}
json.dump(dd, open(f"{D}/dataset_description.json", "w"), indent=2, ensure_ascii=False)
L = []
A = L.append
A(f"# Stimulation onsets detected from the stimulation artefact ({NM})\n")
A("**These events are NOT from a source.** The release does not annotate them. They were detected in Oct 2026 by the "
  "NEMAR curators from the stimulation artefact visible in the recorded signals, after an external review "
  "(J. Gadonneix, Oct 2026) noted that the artefact places them. They are an analysis product: check them against "
  "your own artefact detection before relying on them. The raw dataset's own `events.tsv` files are unchanged.\n")
A("## Method\n")
A("Code: `code/fix-20261030_stimartefact/` of the raw dataset (`stimfeat.py` features, `stimdetect.py` rules, "
  "`build_derivative.py` this tree). Data were read by annex key from the NEMAR S3 bucket (size and SHA-256 checked).\n")
A("1. Contacts: channels typed ECOG, SEEG or DBS in `channels.tsv` with status not bad (recordings whose contacts "
  "are all typed MISC use the MISC channels).")
A("2. Step score per sample: 90th percentile across contacts of |x[n+1] - x[n]| divided by that contact's robust SD of "
  "steps (1.4826 x median absolute deviation). Stimulation pulses are steps far above the normal sample-to-sample "
  "change in most contacts at once.")
A("3. Single pulses: clusters of samples with step score > 50, at least 1 s (SPES) or 0.3 s (SPS) apart. A candidate "
  "is kept only if it is part of a chain of at least 5 pulses spaced at the protocol interval stated in the article "
  "(SPES: 15 s +- 0.5 s; SPS: 1 s +- 0.1 s). Other candidates are rejected and listed below.")
A("4. 43.3 Hz trains (nm000372 task-hfs): step score > 50 at least 10 ms apart, grouped while gaps are <= 0.1 s; a "
  "group with >= 20 pulses and a median interval of 23.1 +- 2 ms (43.3 Hz, article) is a train.")
A("5. 100 Hz blocks (nm000380 task-stim): signals resampled to about 1 kHz; spectrogram (1-s Hann windows, 0.1-s "
  "steps), median across contacts of log10 power; comb score = mean over 100, 200 and 400 Hz of log10 P(f) minus the "
  "mean log10 power 6-10 Hz away. Frames with comb score > 0.5 (power at the stimulation lines more than 3.2 x the "
  "neighbouring bins) count as stimulation; gaps <= 5 s are bridged and blocks shorter than 20 s ignored. Onset and "
  "offset are then refined to the first/last step (score > 20, >= 4 ms apart) followed/preceded by >= 8 steps "
  "within 100 ms, searched within 1.5 s of the spectral edge (`method` = pulse_refined); otherwise the spectral edge "
  "is used (`method` = spectral, about +-0.5 s).\n")
A("Thresholds were set once from the score distributions of these recordings (background step score: median "
  "0.7-1.4) and are the same for every run.\n")
A("## Results per recording\n")
for k, v in S.items():
    if "n_blocks" in v:
        A(f"- `{k}`: {v['n_blocks']} block(s) " + "; ".join(f"onset {float(o):.3f} s, duration {float(d):.1f} s ({m})" for o, d, m in v["blocks"])
          + f". Comb score median {v['comb_p50']:.2f}, 95th percentile {v['comb_p95']:.2f}.")
    elif "train_pulses" in v:
        A(f"- `{k}`: {v['n_trains']} trains of {', '.join(map(str, v['train_pulses']))} pulses; "
          f"{v['rejected_groups']} other artefact groups rejected (not 43.3 Hz trains).")
    else:
        A(f"- `{k}`: {v['n_pulses']} pulses kept of {v['n_candidates']} candidates, in {v['n_trains']} chain(s) of "
          f"{', '.join(map(str, v['train_sizes']))}; rejected candidate onsets (s): "
          f"{', '.join(f'{x:.3f}' for x in v['rejected_onsets_s']) or 'none'}.")
notes = sys.argv[3] if len(sys.argv) > 3 else None
if notes:
    A("\n## Notes\n")
    A(open(notes).read().strip())
open(f"{D}/README.md", "w").write("\n".join(L) + "\n")
print("sidecars", n_files)
