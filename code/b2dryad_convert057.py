"""IEEG057 (Dryad doi:10.5061/dryad.7h44j0ztk; Bratu, Oane et al. 2021 Cortex) -> iEEG-BIDS, one patient (SEEG84).

Signals: AnyWave ADES (.ades text header + .dat little-endian float32, channel-multiplexed) are byte-identical to
BrainVision IEEE_FLOAT_32 MULTIPLEXED binary, so the .dat bytes are copied unchanged as the BrainVision .eeg and a
.vhdr/.vmrk is written (resolution 1). No filtering, resampling, re-referencing or channel removal.
  HFS : Data/HFS/SEEG84.*        -> task-hfs  (43 Hz high-frequency stimulation, alternating polarity)
  SPES: Data/SPES/SEEG84 SPES.*  -> task-spes (single-pulse electrical stimulation)
Events: AnyWave marker files (.mrk = acquisition-level markers, .ades.mrk = authors' analysis markers), copied.
Channel status: AnyWave .ades.bad lists (authors' bad-channel marks).
Electrodes: Data/Contact_coordinates.xlsx columns xmri/ymri/zmri (verified = FreeSurfer tkrRAS of MRI_defaced:
x = 128 - xvox, y = zvox - 128, z = 128 - yvox).
Anatomy: Data/MRI_defaced.mgz -> anat/sub-01_T1w.nii.gz (voxel data and affine unchanged; MGH trailer tags dropped) unless --no-anat.
Privacy: the de-identified sourcedata subset excludes every zip member whose name or bytes contain the patient
name/scan-date tokens (DSI Studio Fibers/*.src.gz,*.fib.gz,*.mapping.txt; SEEG84Plan.ppr; HFS 'SEEG 84 - Z-test.xls').
Tokens are derived programmatically and never written out.
Usage: python b2dryad_convert057.py <src_dir> <bids_root> [--no-anat]
"""
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import zipfile

import numpy as np

SRC, OUT = sys.argv[1], sys.argv[2]
WITH_ANAT = "--no-anat" not in sys.argv  # anat included by default (deface check passed, receipts/IEEG057/facecheck2.png)
Z = zipfile.ZipFile(os.path.join(SRC, "Autoscopic.zip"))
NAMES = Z.namelist()
SUB = "sub-01"
D = os.path.join(OUT, SUB, "ieeg")
os.makedirs(D, exist_ok=True)


def wtsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(["n/a" if (v is None or v == "") else v for v in r])


def wjson(p, o):
    with open(p, "w") as f:
        json.dump(o, f, indent=2, ensure_ascii=False)
        f.write("\n")


def sha(b):
    return hashlib.sha256(b).hexdigest()


# --- identifying tokens (never printed)
srcname = [n for n in NAMES if n.endswith(".src.gz")][0]
b0 = os.path.basename(srcname).split(".")[0]
TOK = [t for t in b0.split("_") if re.fullmatch(r"[A-Z]{3,}", t)] + [t for t in b0.split("_") if re.fullmatch(r"\d{8}", t)]


def member_ok(n):
    if any(t in n.upper() for t in TOK):
        return False
    i = Z.getinfo(n)
    if i.is_dir() or n.endswith(".dat"):
        return True
    if i.file_size < 80e6:
        b = Z.read(n)
        if n.endswith((".gz", ".mgz")):
            try:
                b = gzip.decompress(b)
            except Exception:
                pass
        if any(t.encode() in b.upper() for t in TOK):
            return False
        if n.endswith(".xlsx"):
            zz = zipfile.ZipFile(io.BytesIO(Z.read(n)))
            if any(any(t.encode() in zz.read(m).upper() for t in TOK) for m in zz.namelist()):
                return False
    return True


def parse_ades(n):
    meta, chans = {}, []
    for l in Z.read(n).decode("latin-1").splitlines():
        l = l.strip()
        if not l or l.startswith("#") or "=" not in l:
            continue
        k, v = [x.strip() for x in l.split("=", 1)]
        if k in ("samplingRate", "numberOfSamples"):
            meta[k] = float(v) if k == "samplingRate" else int(v)
        else:
            chans.append((k, v))
    return meta, chans


def parse_mrk(n):
    rows = []
    if n not in NAMES:
        return rows
    for l in Z.read(n).decode("latin-1").splitlines():
        if not l.strip() or l.startswith("//"):
            continue
        p = [x for x in l.rstrip("\n").split("\t")]
        p = [x.strip() for x in p if x is not None]
        while p and p[-1] == "":
            p.pop()
        label, value, pos = p[0], p[1], float(p[2])
        dur = float(p[3]) if len(p) > 3 and p[3] != "" else 0.0
        targets = ",".join(p[4:]) if len(p) > 4 else "n/a"
        rows.append((pos, dur, label, value, targets))
    return rows


TYPEMAP = {"SEEG": "SEEG", "OTHER": "MISC", "TRIGGER": "TRIG"}
NONNEURAL = {"OSAT", "PR", "Pleth"}
RUNS = [("hfs", "Autoscopic/Data/HFS/SEEG84", "43 Hz high-frequency electrical stimulation (alternating polarity pulses) of contacts D08-D09 in the right posterior temporal periventricular nodular heterotopia, at 0.25-0.6 mA, with clinical testing; autoscopic hallucinations (body perception) were elicited at >=0.5 mA (marker labels)."),
        ("spes", "Autoscopic/Data/SPES/SEEG84 SPES", "Single-pulse electrical stimulation (SPES) used to map stimulation-evoked potentials / effective connectivity; marker 'SPES D08 D09'.")]
report = {"runs": []}
for task, stem_src, desc in RUNS:
    meta, chans = parse_ades(stem_src + ".ades")
    fs, ns = meta["samplingRate"], meta["numberOfSamples"]
    nch = len(chans)
    datn = stem_src + ".dat"
    size = Z.getinfo(datn).file_size
    assert size == ns * nch * 4, (task, size, ns, nch)
    stem = f"{SUB}_task-{task}_run-1"
    eeg = os.path.join(D, stem + "_ieeg.eeg")
    h = hashlib.sha256()
    with Z.open(datn) as fi, open(eeg, "wb") as fo:
        for blk in iter(lambda: fi.read(1 << 24), b""):
            h.update(blk)
            fo.write(blk)
    src_sha = h.hexdigest()
    bad = set()
    if stem_src + ".ades.bad" in NAMES:
        bad = {l.strip() for l in Z.read(stem_src + ".ades.bad").decode("latin-1").splitlines() if l.strip()}
    flt = Z.read(stem_src + ".ades.flt").decode("latin-1") if stem_src + ".ades.flt" in NAMES else ""
    # vhdr
    lines = ["Brain Vision Data Exchange Header File Version 1.0",
             f"; Written by b2dryad_convert057.py from AnyWave ADES {os.path.basename(stem_src)}.ades (binary copied byte-for-byte)", "",
             "[Common Infos]", "Codepage=UTF-8", f"DataFile={stem}_ieeg.eeg", f"MarkerFile={stem}_ieeg.vmrk",
             "DataFormat=BINARY", "DataOrientation=MULTIPLEXED", f"NumberOfChannels={nch}",
             f"SamplingInterval={1e6 / fs!r}", "", "[Binary Infos]", "BinaryFormat=IEEE_FLOAT_32", "", "[Channel Infos]"]
    for i, (nm, _) in enumerate(chans, 1):
        lines.append(f"Ch{i}={nm.replace(',', chr(92) + '1')},,1,µV")
    open(os.path.join(D, stem + "_ieeg.vhdr"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    # markers
    acq = parse_mrk(stem_src + ".mrk")
    ana = parse_mrk(stem_src + ".ades.mrk")
    vm = ["Brain Vision Data Exchange Marker File, Version 1.0", "", "[Common Infos]", "Codepage=UTF-8",
          f"DataFile={stem}_ieeg.eeg", "", "[Marker Infos]", "Mk1=New Segment,,1,1,0"]
    allm = [(p, d, l, v, t, "acquisition_mrk") for p, d, l, v, t in acq] + [(p, d, l, v, t, "analysis_ades_mrk") for p, d, l, v, t in ana]
    allm.sort(key=lambda x: (x[0], x[5]))
    for k, (p, d, l, v, t, s) in enumerate(allm, 2):
        vm.append(f"Mk{k}=Comment,{l.replace(',', chr(92) + '1')},{int(round(p * fs)) + 1},{max(1, int(round(d * fs)))},0")
    open(os.path.join(D, stem + "_ieeg.vmrk"), "w", encoding="utf-8").write("\n".join(vm) + "\n")
    wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "value", "sample", "targets", "source"],
         [["%.6f" % p, "%.6f" % d, l, ("n/a" if v == "-1" else v), int(round(p * fs)), t, s] for p, d, l, v, t, s in allm])
    # channels
    rows, cnt = [], {}
    for nm, ty in chans:
        t = "MISC" if nm in NONNEURAL else TYPEMAP.get(ty, "MISC")
        cnt[t] = cnt.get(t, 0) + 1
        grp = re.sub(r"\d+$", "", nm) if t == "SEEG" else "n/a"
        st = "bad" if nm in bad else "good"
        rows.append([nm, t, "µV" if t == "SEEG" else "n/a", "n/a", "n/a", fs, grp, st,
                     "listed in the authors' AnyWave .ades.bad file" if st == "bad" else "n/a",
                     f"ADES type {ty}" + ("; physiological monitor channel typed SEEG in the ADES header" if nm in NONNEURAL else "")])
    wtsv(os.path.join(D, stem + "_channels.tsv"),
         ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "group", "status", "status_description", "description"], rows)
    stim = sorted({l for _, _, l, _, _, _ in allm if re.search(r"mA|SPES|Hz", l)})
    wjson(os.path.join(D, stem + "_ieeg.json"), {
        "TaskName": task, "TaskDescription": desc,
        "SamplingFrequency": fs, "PowerLineFrequency": 50,
        "SoftwareFilters": {"AnyWave display filter file (.ades.flt, not applied to the stored samples)": json.loads(flt) if flt.strip().startswith("{") else "n/a"},
        "HardwareFilters": "n/a",
        "iEEGReference": "n/a (not stated in the release; signals as exported by the authors in AnyWave ADES format)",
        "RecordingDuration": ns / fs, "RecordingType": "continuous",
        "SEEGChannelCount": cnt.get("SEEG", 0), "MiscChannelCount": cnt.get("MISC", 0), "TriggerChannelCount": cnt.get("TRIG", 0),
        "ECOGChannelCount": 0,
        "ElectricalStimulation": True,
        "ElectricalStimulationParameters": "Stimulation events as labelled in the AnyWave markers: " + "; ".join(stim),
    })
    report["runs"].append({"task": task, "source": datn, "sha256_dat": src_sha, "n_channels": nch, "n_samples": ns, "sfreq": fs,
                           "n_markers_acq": len(acq), "n_markers_analysis": len(ana), "bad_channels": sorted(bad)})
    print("run", task, nch, ns, fs, len(allm), flush=True)

# electrodes (tkrRAS mm)
import openpyxl
wb = openpyxl.load_workbook(io.BytesIO(Z.read("Autoscopic/Data/Contact_coordinates.xlsx")), read_only=True, data_only=True)
ws = wb.worksheets[0]
it = ws.iter_rows(values_only=True)
hdr = [str(h) for h in next(it)]
ix = {h: i for i, h in enumerate(hdr)}
erows, bad_tkr = [], 0
for r in it:
    if r[ix["name"]] is None:
        continue
    x, y, z = (float(r[ix[k]]) for k in ("xmri", "ymri", "zmri"))
    xv, yv, zv = (float(r[ix[k]]) for k in ("xmrivox", "ymrivox", "zmrivox"))
    bad_tkr += not (abs(x - (128 - xv)) < 1e-6 and abs(y - (zv - 128)) < 1e-6 and abs(z - (128 - yv)) < 1e-6)
    erows.append([r[ix["name"]], "%.6f" % x, "%.6f" % y, "%.6f" % z, "n/a", "depth", re.sub(r"\d+$", "", str(r[ix["name"]])), r[ix["hemi"]], r[ix["exact"]]])
wtsv(os.path.join(D, f"{SUB}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "type", "group", "hemisphere", "freesurfer_label"], erows)
wjson(os.path.join(D, f"{SUB}_space-Other_electrodes.json"), {
    "hemisphere": {"Description": "Hemisphere (column 'hemi' of Contact_coordinates.xlsx)"},
    "freesurfer_label": {"Description": "Anatomical label at the contact (column 'exact' of Contact_coordinates.xlsx; FreeSurfer aseg/aparc naming)"},
    "type": {"Description": "Electrode type; the release describes stereo-EEG depth electrodes"}})
cs = {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "mm",
      "iEEGCoordinateSystemDescription": "FreeSurfer surface RAS (tkrRAS) of the patient's T1 image MRI_defaced.mgz: columns xmri/ymri/zmri of Data/Contact_coordinates.xlsx. Verified for every contact against the voxel columns: x = 128 - xmrivox, y = zmrivox - 128, z = 128 - ymrivox (the FreeSurfer conformed-space tkr transform).",
      "iEEGCoordinateProcessingDescription": "Contact coordinates as provided by the authors (Contact_coordinates.xlsx)."}
if WITH_ANAT:
    cs["IntendedFor"] = f"bids::{SUB}/anat/{SUB}_T1w.nii.gz"
wjson(os.path.join(D, f"{SUB}_space-Other_coordsystem.json"), cs)
report["electrodes"] = {"n": len(erows), "tkr_mismatch": bad_tkr}
if WITH_ANAT:
    import nibabel as nib
    a = os.path.join(OUT, SUB, "anat")
    os.makedirs(a, exist_ok=True)
    p = "/work/MRI_defaced.mgz"
    open(p, "wb").write(Z.read("Autoscopic/Data/MRI_defaced.mgz"))
    img = nib.load(p)
    nib.save(nib.Nifti1Image(np.asanyarray(img.dataobj), img.affine), os.path.join(a, f"{SUB}_T1w.nii.gz"))
    wjson(os.path.join(a, f"{SUB}_T1w.json"), {"Description": "Authors' defaced and anonymized T1-weighted MRI (Data/MRI_defaced.mgz), converted to NIfTI with voxel data and affine unchanged."})
    os.remove(p)
# participants
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "source_id", "age", "sex", "handedness"], [[SUB, "SEEG84", "n/a", "n/a", "n/a"]])
wjson(os.path.join(OUT, "participants.json"), {"source_id": {"Description": "Patient code used in the release (SEEG84)"},
                                              "age": {"Description": "not given in the Dryad release", "Units": "year"},
                                              "sex": {"Description": "not given in the Dryad release"},
                                              "handedness": {"Description": "not given in the Dryad release"}})
# de-identified sourcedata subset (member list + exclusions)
dst = os.path.join(OUT, "sourcedata", "dryad-7h44j0ztk-deidentified")
excl = []
for n in NAMES:
    if n.endswith("/"):
        continue
    if not member_ok(n):
        excl.append(n)
        continue
    if n.endswith(".dat"):
        continue  # signal bytes are already in the BIDS .eeg files (identical)
    t = os.path.join(dst, n)
    os.makedirs(os.path.dirname(t), exist_ok=True)
    with Z.open(n) as fi, open(t, "wb") as fo:
        shutil.copyfileobj(fi, fo, 1 << 22)
shutil.copy(os.path.join(SRC, "ReadMe.txt"), os.path.join(dst, "ReadMe.txt"))
report["sourcedata_excluded_count"] = len(excl)
report["sourcedata_excluded_kinds"] = sorted({("Fibers/DSI-Studio source+fib+mapping (identifying file names)" if "/Fibers/" in n and (n.endswith(".gz") and "src" in n or "mapping" in n) else os.path.basename(n) if not any(t in n.upper() for t in TOK) else "Fibers/<identifying file name>") for n in excl})
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
shutil.copy(__file__, os.path.join(OUT, "code", os.path.basename(__file__)))
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1)
wtsv(os.path.join(OUT, SUB, f"{SUB}_scans.tsv"), ["filename", "acq_time"], [[f"ieeg/{SUB}_task-{t}_run-1_ieeg.vhdr", "n/a"] for t, _, _ in RUNS] + ([[f"anat/{SUB}_T1w.nii.gz", "n/a"]] if WITH_ANAT else []))
print(json.dumps({k: v for k, v in report.items() if k != "runs"}, indent=1))
print("DONE")
