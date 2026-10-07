PPG AND PRESSURE RECORDINGS FROM A MULTILAYER VASCULAR FINGER PHANTOM
WITH PALE, MEDIUM AND DARK SKIN TONE LAYERS
======================================================================

Authors:  Laura Osorio-Sanchez, James M. May, Panicos Kyriacou
          Research Centre for Biomedical Engineering,
          City St George's, University of London, UK

Article:  "Photoplethysmography (PPG) signal quality and morphology features
          across skin pigmentation levels using a multilayer vascular finger
          phantom" (manuscript submitted)

Code:     https://github.com/Lauraosorios/ppg-pigmentation-phantom-quality-morphology

Dataset:  https://doi.org/10.5281/zenodo.23208189

Licence:  Creative Commons Attribution 4.0 International (CC BY 4.0)


1. EXPERIMENT
-------------
A multilayer vascular finger phantom with interchangeable skin tone layers
(pale, medium and dark; approximately Fitzpatrick I-II, III-IV and V-VI) was
perfused by a programmable pulsatile pump (BDC Laboratories PD-1100) with a
blood-mimicking fluid (0.1% w/v Nigrosine in deionised water).

  Skin tone layer    Pale (P), Medium (M), Dark (D)
  Pump heart rate    60, 90, 120 bpm
  Pump target flow   5, 6, 7 L/min
  Replicates         3 per condition

3 x 3 x 3 x 3 = 81 recordings.


2. FOLDER STRUCTURE
-------------------
  README.txt
  PPG.zip       (unzips to PPG/)       81 files: reflectance + transmittance PPG
  Pressure.zip  (unzips to Pressure/)  81 files: aortic + finger-phantom pressure
  metadata/
      recordings.csv     one row per recording: factors, files, sample counts
      SHA256SUMS.txt     checksums of all data files
  derived/
      recording_sqi_table.csv          per recording x channel median of
                                       9 pulse-level signal quality indices
      recording_morphology_table.csv   per recording x channel median of
                                       pyPPG morphology features


3. FILE NAMING
--------------
  <skin>R<replicate>_<heart rate>HR_<target flow>CO.csv

  Example: DR1_60HR_5CO.csv = Dark skin tone, replicate 1, 60 bpm,
  5 L/min target flow. "CO" denotes the pump target flow setting (L/min).
  The PPG and pressure files of the same recording share the same name.


4. PPG FILES (PPG/)
-------------------
Sensor:       custom finger probe; SFH7016 LEDs (green 530 nm, red 665 nm,
              infrared 940 nm) driven at 30 mA; TEMD5080X01 photodiodes.
              Reflectance detector at 4 mm source-detector separation;
              transmittance detector opposite the LEDs.
Front end:    AFE4420 (Texas Instruments), 200 Hz, all channels sampled
              simultaneously.
Length:       65,224-65,536 samples (326-328 s); see metadata/recordings.csv.
Units:        volts.

Columns:
  Reflectance RED       665 nm, reflectance detector
  Reflectance IR        940 nm, reflectance detector
  Reflectance GREEN     530 nm, reflectance detector
  Transmittance RED     665 nm, transmittance detector
  Transmittance IR      940 nm, transmittance detector
  Transmittance GREEN   530 nm, transmittance detector. Recorded but NOT
                        usable: green light does not penetrate the phantom,
                        so there is no pulsatile signal. Excluded from all
                        analysis.

Processing already applied:
  1. ADC codes converted to volts (code x 1.2 / 2^21).
  2. AFE offset-cancellation correction: the DC offset current subtracted by
     the front end was added back for each channel
     (V + 2 x I_offset x R_f, R_f = 20 kOhm), restoring the true DC level.
No filtering was applied.

NOTE: the pulsatile (AC) component is INVERTED relative to the usual PPG
convention because of the front-end sign convention. The analysis code
negates the band-passed AC component once before processing.


5. PRESSURE FILES (Pressure/)
-----------------------------
Sensors:      PendoTECH single-use pressure sensors at the ascending aorta
              segment and at the finger-phantom inlet.
Sampling:     acquired at 2 kHz, decimated to 200 Hz (anti-aliasing FIR,
              zero-phase, scipy.signal.decimate).
Units:        mmHg.

Columns:
  Aorta_Pressure     aortic segment pressure
  Finger_Pressure    finger-phantom inlet pressure


6. IMPORTANT NOTES
------------------
- PPG and pressure are NOT time-synchronised. They were logged by separate
  instruments, each started separately for every recording, and the files
  carry no timestamps. Use pressure to characterise each haemodynamic
  condition (e.g. mean arterial pressure per recording), not for
  beat-to-beat alignment with the PPG.
- Green transmittance is included for completeness but contains no usable
  signal (see section 4).


7. DERIVED TABLES (derived/)
----------------------------
Recording-level inputs to the statistical analysis in the article: one row
per recording x channel (81 x 5 = 405 rows), each value the median over all
pulses detected by pyPPG 1.0.73. They can be regenerated from PPG/ with the
code repository (scripts/01_extract_features.py).

Key columns: skin (P/M/D), replicate, hr (bpm), co_lmin (target flow,
L/min), channel.


8. CITATION
-----------
Please cite the article above and this dataset:
  https://doi.org/10.5281/zenodo.23208189
