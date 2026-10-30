# Stimulation onsets detected from the stimulation artefact (nm000372)

**These events are NOT from a source.** The release does not annotate them. They were detected in Oct 2026 by the NEMAR curators from the stimulation artefact visible in the recorded signals, after an external review (J. Gadonneix, Oct 2026) noted that the artefact places them. They are an analysis product: check them against your own artefact detection before relying on them. The raw dataset's own `events.tsv` files are unchanged.

## Method

Code: `code/fix-20261030_stimartefact/` of the raw dataset (`stimfeat.py` features, `stimdetect.py` rules, `build_derivative.py` this tree). Data were read by annex key from the NEMAR S3 bucket (size and SHA-256 checked).

1. Contacts: channels typed ECOG, SEEG or DBS in `channels.tsv` with status not bad (recordings whose contacts are all typed MISC use the MISC channels).
2. Step score per sample: 90th percentile across contacts of |x[n+1] - x[n]| divided by that contact's robust SD of steps (1.4826 x median absolute deviation). Stimulation pulses are steps far above the normal sample-to-sample change in most contacts at once.
3. Single pulses: clusters of samples with step score > 50, at least 1 s (SPES) or 0.3 s (SPS) apart. A candidate is kept only if it is part of a chain of at least 5 pulses spaced at the protocol interval stated in the article (SPES: 15 s +- 0.5 s; SPS: 1 s +- 0.1 s). Other candidates are rejected and listed below.
4. 43.3 Hz trains (nm000372 task-hfs): step score > 50 at least 10 ms apart, grouped while gaps are <= 0.1 s; a group with >= 20 pulses and a median interval of 23.1 +- 2 ms (43.3 Hz, article) is a train.
5. 100 Hz blocks (nm000380 task-stim): signals resampled to about 1 kHz; spectrogram (1-s Hann windows, 0.1-s steps), median across contacts of log10 power; comb score = mean over 100, 200 and 400 Hz of log10 P(f) minus the mean log10 power 6-10 Hz away. Frames with comb score > 0.5 (power at the stimulation lines more than 3.2 x the neighbouring bins) count as stimulation; gaps <= 5 s are bridged and blocks shorter than 20 s ignored. Onset and offset are then refined to the first/last step (score > 20, >= 4 ms apart) followed/preceded by >= 8 steps within 100 ms, searched within 1.5 s of the spectral edge (`method` = pulse_refined); otherwise the spectral edge is used (`method` = spectral, about +-0.5 s).

Thresholds were set once from the score distributions of these recordings (background step score: median 0.7-1.4) and are the same for every run.

## Results per recording

- `sub-01_task-hfs_run-1`: 6 trains of 218, 218, 218, 218, 218, 218 pulses; 7 other artefact groups rejected (not 43.3 Hz trains).
- `sub-01_task-spes_run-1`: 20 pulses kept of 21 candidates, in 1 chain(s) of 20; rejected candidate onsets (s): 7.134.

## Notes

- task-hfs: the 6 detected trains start within 2 ms of the 6 trains decoded from the TRIG channel (rows with `source` = TRIG in the raw `events.tsv`), including the first train at 23.68 s that has no marker in the release; each has 218 pulses (5.05 s at 43.3 Hz). The other 7 artefact groups (e.g. at 260.7 s, next to the 'De-block start' marker) are not 43.3 Hz trains and are not listed.
- task-spes: the 20 pulses are 15.000 s apart (11.57 s to 296.57 s), matching the article (20 pulses per trial, 15 s interval). The TRIG channel does not mark them (it changes level at 0-11.57 s and 311.53-320.0 s only). One candidate at 7.13 s is not on the 15-s chain and is rejected.
