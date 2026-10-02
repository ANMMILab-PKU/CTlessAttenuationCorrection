# CTlessAttenuationCorrection

## Lutetium-176-Aided Joint Activity and Attenuation Reconstruction in PET

This repository contains the processing and reconstruction pipeline for PET simulation data generated with **GATE**. The project investigates attenuation correction without a conventional CT scan by combining:

- **511 keV back-to-back (B2B) emission events** from the PET activity distribution;
- **Lutetium-176 (Lu-176) intrinsic radioactivity**, used to generate transmission events at approximately **202 keV** and **307 keV**;
- ground-truth labels available from simulation, including event IDs, source IDs, and phantom-Compton information;
- iterative attenuation and activity reconstruction algorithms implemented with MATLAB, C/C++, OpenMP, CUDA, and the `parallelproj` projection library.

The complete workflow converts ROOT singles data into compact binary PET data, sorts singles into coincidence events, classifies coincidences using simulation labels, generates unified histograms, and reconstructs attenuation maps using several increasingly advanced methods.

> **Research-code note:** This repository is organized around a specific TBPET/BrainPET-style simulated scanner geometry and a Linux/HPC execution environment. Many scripts contain experiment-specific paths, filenames, detector parameters, and compiled-program names. Before running the pipeline on a new dataset, update those values and verify the binary formats carefully.

---

## 1. Project Overview

The processing chain is divided into seven major stages:

```text
GATE ROOT simulation output
            |
            v
[2] ROOT Singles Extraction and Binary Conversion
            |
            v
[3] Time Sorting and Coincidence Formation
            |
            v
[4] Coincidence Classification and CDF/List-Mode Export
            |
            v
[5] Duplicate Removal and Histogram Generation
            |
            +------------------------------+
            |                              |
            v                              v
[6] MLTR transmission initialization   Unified histogram data
                                           |
                                           v
                         [7] Emission-only MLAA reconstruction
                                           |
                                           v
                         [8] MLAA-TX reconstruction
                                           |
                                           v
                         [9] TR-SCS reconstruction
```

The reconstruction stages are complementary:

- **MLTR** uses Lu-176 transmission data to obtain a coarse initial attenuation map.
- **Emission-only MLAA** estimates attenuation using 511 keV emission data and activity reconstruction.
- **MLAA-TX** jointly uses 511 keV emission data and Lu-176 transmission data, including multi-energy attenuation mapping and additive-event correction.
- **TR-SCS** uses the same physical data sources with a more conservative surrogate/cubic optimization strategy designed to improve numerical stability and handle nonconvexity or weak-curvature regions.

---

## 2. Data Sources and Physical Interpretation

### 2.1 511 keV B2B emission events

The 511 keV emission component is produced by positron annihilation in the GATE simulation. Valid emission coincidences are selected from singles that satisfy timing, energy, detector-separation, ring-difference, and line-of-response (LOR) constraints.

These events are used for:

- emission sinogram/list-mode generation;
- activity reconstruction with TOF MLEM;
- emission-only MLAA attenuation reconstruction;
- the emission term in MLAA-TX and TR-SCS.

### 2.2 Lu-176 transmission events

Lu-176 intrinsic radioactivity inside LSO/LYSO detector material provides a distributed internal source. The coincidence sorter separates transmission-related events into energy windows centered around approximately:

- **202 keV**;
- **307 keV**;
- **511 keV**, when appropriate.

Blank-scan transmission data provide the incident-count reference, while phantom-scan transmission data provide measurements attenuated by the object. The reconstruction model generally follows a Poisson form such as:

```text
p_i ~ Poisson(B_i exp(-A_i mu) + r_i)
```

where:

- `B_i` is the blank/transmission reference count;
- `A_i mu` is the LOR attenuation line integral;
- `r_i` is an additive term, such as scatter or randoms;
- `p_i` is the measured phantom transmission count.

### 2.3 Simulation labels

The ROOT extraction stage preserves event-level labels, including:

- `eventID`;
- `sourceID`;
- `comptonPhantom`.

These labels are carried through the `.reference` file and later stored in `.consist` files. They enable classification into categories such as true events and non-true/additive events.

---

## 3. Repository Organization

### `2. singlesProcessing`

Converts GATE ROOT `Singles` trees into binary singles streams.

Important programs include:

- `ROOT_singlesExtract_TBPET_22panel_1axial_6p3p4_DOI_energy_time_root_20260612_1.C`
  - reads multiple ROOT files;
  - accesses the `Singles` tree;
  - converts hierarchical detector IDs into a global crystal ID;
  - writes binary arrays for crystal IDs, time, energy, event ID, source ID, and Compton labels.
- `convert_dat4Sort_sortAfter_timeTagEnd_withEIandSIandCP_1axial_22panel.c`
  - reads the extracted arrays;
  - converts time from seconds to discrete CFD time bits;
  - sorts events by time within chunks;
  - emits sorter-compatible 64-bit singles words;
  - writes a matching `.reference` file.

### `3. sorting`

Runs the coincidence sorter. The included shell script demonstrates how the 26-bit time-tag sorter is called with a configuration file and a `.dat` singles stream.

The sorter produces several coincidence streams, including:

- `511_prompt.dat`;
- `511_delay.dat`;
- `307_prompt.dat`;
- `307_delay.dat`;
- `202_prompt.dat`;
- `202_delay.dat`;
- backscatter streams;
- matching `.consist` files.

### `4. coincidencesProcessing`

Converts sorted coincidence files into compact CDF and time-difference files while using `.consist` labels for event classification.

The transfer programs include variants for:

- all events;
- only true events with `sourceID == 0`;
- only non-true/additive events with `sourceID != 0`.

Outputs include:

- `.Cdf`: records of the form `(1, crystalID1, crystalID2)`;
- `.timeDiff`: signed time differences in picoseconds;
- `.consist`: six `uint32` values per coincidence for consistency and backscatter metadata.

### `5. transfer2histogram`

Converts CDF files into compact histogram files and aligns multiple histogram files to the same sorted crystal-ID-pair index.

Main operations:

1. normalize each detector pair so that `id1 <= id2`;
2. count duplicate occurrences of each detector pair;
3. externally sort large files in chunks;
4. merge temporary sorted files;
5. generate unified histogram files in which corresponding detector pairs occupy the same row across datasets.

The unified representation is required by the MLTR, MLAA-TX, and TR-SCS reconstruction stages.

### `6. MLTR`

Contains transmission-based maximum-likelihood attenuation reconstruction code used to generate a coarse initial attenuation map.

The main implementation:

- reads blank, transmission, and additive histograms in chunks;
- computes forward projections using `parallelproj`;
- evaluates the transmission Poisson model;
- accumulates numerator and denominator terms by backprojection;
- optionally applies a quadratic prior;
- updates the attenuation map while enforcing non-negativity.

The implementation supports different transmission energies and includes a cylindrical field-of-view mask in the reconstruction driver.

### `7. MLAA-emissionOnly`

Implements emission-only MLAA reconstruction.

The workflow alternates between:

1. attenuation-map updates using the symmetric-pixel-search (SPS) formulation;
2. activity-map updates using TOF MLEM;
3. sensitivity recomputation under the current attenuation map.

The implementation supports:

- chunked list-mode processing;
- CUDA forward and backprojection;
- TOF binning with a 7 ps nominal bin width;
- normalization factors;
- attenuation correction;
- activity and attenuation masks.

### `8. MLAA-TX`

Implements joint MLAA-TX reconstruction using:

- 511 keV emission data;
- 202 keV Lu-176 transmission data;
- 307 keV Lu-176 transmission data;
- blank-scan references;
- additive-event histograms;
- an image-domain, three-segment piecewise-linear mapping from 511 keV attenuation to lower-energy attenuation;
- optional quadratic regularization.

The main attenuation update combines emission and transmission terms using a weighted surrogate approach. The transmission terms use the chain rule through the energy-dependent attenuation mapping.

### `9. TR-SCS`

Contains a more advanced transmission-reconstruction framework based on a trust-region/sequential convex surrogate strategy.

Its main components include:

- exact or stabilized emission derivatives;
- transmission first-, second-, and third-order derivative bounds;
- Fisher-information-based data-source weight balancing;
- three-segment low-energy attenuation mapping;
- quadratic regularization;
- a box-constrained cubic surrogate solver;
- adaptive cubic-curvature scaling and trust-region acceptance logic;
- alternating attenuation and activity updates.

---

## 4. Stage 1: Extracting ROOT Singles Data

The ROOT extraction macro reads a sequence of ROOT files containing a `Singles` tree. The expected tree branches include fields such as:

```text
crystal hierarchy: level1ID ... level5ID
time:             time
energy:           energy
event labels:     eventID, sourceID, comptonPhantom
positions:        globalPosX/Y/Z, sourcePosX/Y/Z
```

The detector hierarchy is converted into a global crystal ID. For the 22-sector, one-axial-module geometry used in the supplied scripts:

```text
transaxial crystals per ring = 22 × 1 × 6 × 1 × 8 = 1056
```

The extraction macro writes separate binary files, for example:

```text
<base>.crystalID
<base>.time
<base>.energy
<base>.eventID
<base>.sourceID
<base>.comptonPhantom
```

The files must contain the same number of records and must remain aligned by index.

### Important implementation considerations

- The ROOT macro currently uses experiment-specific paths and filename templates.
- It assumes the required detector branches exist in every ROOT file.
- The generated arrays are written in binary form; downstream programs depend on the exact C/C++ data types.
- Before processing a new simulation, verify that detector ID ordering and axial/transaxial indexing match the LUT and sorter geometry.

---

## 5. Stage 2: Building Sorter-Compatible Singles Data

The converter `convert_dat4Sort_sortAfter_timeTagEnd_withEIandSIandCP_1axial_22panel.c` combines the extracted arrays and creates a sorter input stream.

### 5.1 Time conversion

The input ROOT time values are stored as `double` values in seconds. They are converted to discrete time bits using:

```text
CFD_RESOLUTION = 7e-12 seconds
 timeBits = round(time / CFD_RESOLUTION)
```

The converter sorts events by `timeBits` inside one-million-event chunks.

### 5.2 64-bit singles word layout

The current 26-bit format is:

```text
bits  0–7   : singles tag, 0x82
bits  8–18  : transaxial crystal ID, 11 bits
bits 19–28  : axial ID, 10 bits
bits 29–37  : energy bits, 9 bits
bits 38–63  : time difference/time-tag bits, 26 bits
```

Energy encoding is:

```text
energyBits = round((energy - 0.05 MeV) / 0.0025 MeV)
```

The energy value is clipped into the 9-bit range `[0, 511]`.

Time-tag words use the tag `0x86`. A `.reference` file stores three `uint32` values for every singles event:

```text
[eventID, sourceID, comptonPhantom]
```

The order of records in `.reference` must exactly match the singles words in `.dat`.

### 5.3 Output files

For a base name `<name>`, the converter creates:

```text
<name>.dat
<name>.reference
```

The `.dat` file is consumed by the sorting program; the `.reference` file is consumed by the coincidence sorter to restore event-level labels.

---

## 6. Stage 3: Coincidence Sorting

The sorter groups singles within time-tag periods and creates coincidence records using configurable timing, energy, geometry, and multiplicity rules.

The supplied sorter implementation supports:

- 26-bit time fields;
- prompt and delayed coincidence windows;
- 511 keV emission coincidences;
- 202 keV and 307 keV transmission coincidences;
- sector-difference and ring-difference constraints;
- crystal-position-based LOR/FOV tests;
- `take_all`, `take_winner`, and `kill_all` multiple-coincidence policies;
- backscatter classification;
- OpenMP processing of independent time periods;
- chunked I/O for very large `.dat` files.

The default geometry represented in the supplied configuration is approximately:

```text
crystals per ring       = 1056
crystals per rsector    = 48
number of axial rings   = 120
minimum sector difference = configurable
maximum ring difference  = configurable
```

The sorter uses a crystal LUT to validate LORs and to calculate geometric quantities. Therefore, the LUT used here must correspond to the same detector numbering convention used when generating `.crystalID` and encoding `trID`/`axID`.

### Coincidence output structure

Each coincidence `.dat` record contains six fields:

```c
uint32_t crystalID1;
uint32_t crystalID2;
float    timeDiff;
float    energy1;
float    energy2;
uint32_t type;
```

Each matching `.consist` record contains six `uint32` fields:

```c
uint32_t eid_consist;
uint32_t sid1;
uint32_t sid2;
uint32_t cp1;
uint32_t cp2;
uint32_t bs_type;
```

The `.dat` and `.consist` files must have identical record counts and identical record ordering.

---

## 7. Stage 4: Coincidence Classification and CDF Export

The programs in `4. coincidencesProcessing` read a coincidence `.dat` file together with its `.consist` file.

### 7.1 CDF format

Selected coincidences are converted into three-`uint32` records:

```text
(flag, crystalID1, crystalID2)
```

The flag is set to `1`.

### 7.2 Time-difference format

The coincidence time difference is written as a `float` in picoseconds. The supplied converters apply a sign convention of:

```text
timeDiff_ps = -timeDiff_seconds × 1e12
```

This convention must remain consistent with the TOF reconstruction code.

### 7.3 Label-based selection

The supplied variants include:

- `all`: copy all coincidence records;
- `onlyTrue_sourceID0`: select records satisfying
  `eid_consist == 1`, `sid1 == 0`, `sid2 == 0`, `cp1 == 0`, and `cp2 == 0`;
- `onlyNotTrue_sourceID0`: select the complementary records relative to the same condition.

These files are useful for constructing:

- true-event datasets;
- additive/non-true datasets;
- all-event datasets;
- validation comparisons between simulated truth and reconstructed data.

---

## 8. Stage 5: Histogram Generation and ID Alignment

The CDF files can contain a very large number of repeated crystal pairs. The histogram conversion stage reduces them to unique detector-pair records.

Each histogram record is:

```c
uint32_t counts;
uint32_t id1;
uint32_t id2;
```

Before counting, the pair is normalized:

```text
if id1 > id2: swap(id1, id2)
```

This ensures that `(id1, id2)` and `(id2, id1)` represent the same LOR.

### External merge sort

`convert_Cdf_to_histo_ExternalMergeSort_OMP.c` performs:

1. chunked input reading;
2. parallel pair normalization and local sorting;
3. local duplicate compression;
4. temporary-file generation;
5. serial multiway merging.

The implementation is designed for datasets too large to fit entirely in memory.

### Unified histogram files

`histo_unify_multiWayAlignment.c` aligns multiple histogram files using the union of all detector-ID pairs. For every unique pair, each output file receives one record. If a pair is absent from a particular input file, its count is written as zero.

This is important because MLTR and MLAA-TX compare blank, phantom, true, and additive measurements by row. The corresponding rows must refer to the same normalized LOR.

---

## 9. Stage 6: MLTR Initialization

The MLTR implementation in `6. MLTR` reconstructs a coarse attenuation map from transmission histogram data.

For each chunk of corresponding LORs, it:

1. reads blank counts, transmission counts, and additive counts;
2. scales blank counts to the transmission acquisition duration;
3. obtains crystal coordinates from the LUT;
4. forward-projects the current attenuation map;
5. evaluates the expected transmission signal;
6. adds the additive component;
7. computes the MLTR numerator and surrogate denominator;
8. backprojects the update terms;
9. optionally adds a quadratic prior;
10. updates the attenuation map with a non-negative constraint.

The principal statistical model is:

```text
p_i ~ Poisson(B_i exp(-A_i mu) + r_i)
```

The denominator uses the standard surrogate curvature choice `c_i = u_i`, where `u_i` is the expected attenuated blank signal. The implementation also supports a cylinder mask to suppress attenuation outside the physical field of view.

The MLTR output is intended primarily as an initialization for later MLAA, MLAA-TX, or TR-SCS reconstruction.

---

## 10. Stage 7: Emission-Only MLAA

The emission-only MLAA workflow in `7. MLAA-emissionOnly` alternates between attenuation and activity updates.

### 10.1 Attenuation update

`ML_attn_sps_parallelproj_lm_norm_chunk_cuda.m`:

- reads emission list-mode records in chunks;
- obtains LOR coordinates from the crystal LUT;
- accumulates measured-data backprojections;
- enumerates valid LORs in detector-pair subsets;
- computes attenuation factors using CUDA forward projection;
- applies normalization data;
- computes expected activity projections;
- backprojects expected terms;
- updates the attenuation map using the SPS ratio.

### 10.2 Activity update

`MLEM_parallelproj_TOF_bin_chunk_wAD_cuda.m` performs chunked TOF MLEM. It includes:

- attenuation correction from the current map;
- normalization-factor correction;
- TOF-bin calculation from `.timeDiff` data;
- TOF forward projection;
- measured-to-estimated ratio calculation;
- TOF backprojection;
- non-negative activity updates.

The supplied drivers use a nominal TOF bin width of 7 ps and calculate the required number of bins from the observed time-difference range.

### 10.3 Emission-only additive terms

The current `get_additive_terms.m` implementations return zero additive terms by default. They provide an interface for future scatter/random estimation or externally precomputed additive corrections.

---

## 11. Stage 8: MLAA-TX Reconstruction

The MLAA-TX implementation in `8. MLAA-TX` jointly reconstructs attenuation using emission and Lu-176 transmission information.

The main update function supports three transmission modes:

```text
1 = 202 keV only
2 = 307 keV only
3 = both 202 keV and 307 keV
```

### 11.1 Emission term

The emission component contains a dense sensitivity-like pass over valid LOR geometry and a sparse pass over measured histogram data. The measured data are combined with the current activity map, attenuation map, and normalization factors.

### 11.2 Transmission terms

For each transmission energy, the update:

1. reads phantom transmission counts;
2. reads blank counts and rescales them to the phantom acquisition duration;
3. reads additive counts;
4. maps the 511 keV attenuation map to an equivalent low-energy map;
5. computes the transmission forward projection;
6. evaluates the Poisson gradient and curvature;
7. backprojects the LOR-space quantities;
8. applies the voxel-level chain-rule derivative of the energy mapping.

### 11.3 Three-segment energy mapping

`map_mu511_to_low_energy_3segLM.m` uses three material-inspired regions:

1. air/soft tissue;
2. soft tissue/cancellous bone;
3. cancellous bone/cortical bone.

The mapping returns both:

```text
mu_low(mu_511)
eta'(mu_511)
```

The derivative is used to propagate transmission gradients and curvature back to the 511 keV attenuation map.

### 11.4 Additive correction

Transmission files are kept separate as:

- all transmission counts;
- blank-scan counts;
- additive/non-true counts.

This allows the model to use an additive term instead of treating all measured counts as purely attenuated primary transmission.

---

## 12. Stage 9: TR-SCS Reconstruction

The `9. TR-SCS` directory implements a trust-region/sequential convex surrogate approach for MLAA-TX attenuation updates.

### 12.1 Derivative-based formulation

The method accumulates image-domain terms corresponding to:

- gradient;
- Hessian or surrogate curvature;
- third-order/Lipschitz bounds.

Emission and transmission derivatives are implemented separately. Numerical safeguards include:

- clipping line integrals to a physically safe range;
- lower bounds for denominators;
- local double-precision calculations for curvature evaluation;
- non-negative curvature enforcement;
- non-negative attenuation constraints.

### 12.2 Fisher-information weight balancing

`calc_fisher_information_SCS.m` estimates the relative information content of:

- 511 keV emission;
- 202 keV transmission;
- 307 keV transmission.

It reports recommended weights that approximately balance the Fisher-information contribution of each data source. The recommendation can be multiplied by a user-defined factor to favor transmission data.

### 12.3 Cubic surrogate solver

`solve_cubic_vectorized_final_260502.m` solves a voxel-wise box-constrained cubic surrogate problem. It:

- computes candidate stationary points;
- applies physical bounds to the attenuation coefficient;
- evaluates zero, stationary, lower-bound, and upper-bound candidates;
- selects the candidate with the best surrogate value;
- reports diagnostics such as boundary hits and saddle escapes.

### 12.4 Trust-region adaptation

The SCS update maintains a cubic-curvature scaling factor. Depending on the predicted surrogate behavior, it accepts the proposed update or retains the previous map and increases the conservative scaling factor.

This design is intended to reduce overshooting in regions where the transmission likelihood has weak curvature or where the combined objective is locally nonconvex.

---

## 13. Common Binary File Formats

### 13.1 Singles `.dat`

One little-endian 64-bit word per record:

```text
bits  0–7   tag
bits  8–18  transaxial crystal ID
bits 19–28  axial ID
bits 29–37  energy bits
bits 38–63  time bits
```

Tags:

```text
0x82 = singles event
0x86 = time tag
```

### 13.2 Singles `.reference`

Three little-endian `uint32` values per singles event:

```text
[eventID, sourceID, comptonPhantom]
```

### 13.3 Coincidence `.dat`

Six-field, 24-byte records:

```text
uint32 crystalID1
uint32 crystalID2
float  timeDiff
float  energy1
float  energy2
uint32 type
```

### 13.4 Coincidence `.consist`

Six `uint32` values, 24 bytes per record:

```text
[eid_consist, sid1, sid2, cp1, cp2, bs_type]
```

### 13.5 `.Cdf`

Three `uint32` values per record:

```text
[1, crystalID1, crystalID2]
```

### 13.6 `.timeDiff`

One `float32` value per CDF record, in picoseconds.

### 13.7 `.histo`

Three `uint32` values per unique detector pair:

```text
[counts, id1, id2]
```

All histogram files that will be jointly processed must be produced with the same detector-ID convention and then aligned with the unification program.

---

## 14. Recommended End-to-End Execution Order

A typical experiment should follow this order:

### Step 1: Generate ROOT data with GATE

Create ROOT files containing the `Singles` tree and the required detector, time, energy, and label branches.

### Step 2: Extract binary singles arrays

Run the ROOT macro after updating:

- ROOT input directory;
- ROOT filename template;
- first and last file indices;
- output prefix;
- detector hierarchy parameters.

### Step 3: Encode sorter input

Compile and run the singles encoder. Confirm that all input array files contain the same number of events.

### Step 4: Run the coincidence sorter

Use a sorter configuration containing the desired:

- energy windows;
- timing windows;
- delay offset;
- detector geometry;
- multiplicity policy;
- LUT path;
- output prefix.

### Step 5: Convert coincidence files to CDF/time-difference files

Run the appropriate converter for each energy and event class, for example:

```text
511 prompt true
511 prompt non-true
202 prompt all
202 prompt non-true
307 prompt all
307 prompt non-true
```

### Step 6: Convert CDF files to histograms

Run the external merge-sort histogram converter for each CDF file.

### Step 7: Unify detector-pair IDs

Align blank, phantom, true, non-true, and additive histograms with the multi-file unification program.

### Step 8: Run MLTR

Use blank/transmission/additive histograms to create a coarse attenuation initialization.

### Step 9: Run emission-only MLAA, MLAA-TX, or TR-SCS

Load the MLTR initialization, prepare the required activity initialization, LUT, normalization factors, masks, and data paths, then run the selected iterative reconstruction driver.

---

## 15. Required Software and Hardware

The exact environment depends on the selected stage, but the project generally requires:

- GATE/Geant4 for simulation;
- ROOT for reading `TTree`/`TFile` data;
- a C/C++ compiler;
- OpenMP support;
- MATLAB for reconstruction scripts;
- CUDA-capable GPU for CUDA projection wrappers;
- the `parallelproj` projection library and compatible header/library files;
- MATLAB files containing normalization factors;
- detector crystal LUT files in the expected binary format.

For large simulations, substantial storage and memory are required. The pipeline is designed to process very large files by using:

- chunked binary I/O;
- external sorting;
- batch processing;
- OpenMP parallelization;
- CUDA forward/backprojection.

---

## 16. Configuration Checklist

Before running a new experiment, verify all of the following:

### Detector geometry

- number of rsectors;
- axial and transaxial modules;
- blocks and submodules;
- crystals per block;
- crystals per ring;
- number of axial rings;
- detector ID ordering.

### Time and energy encoding

- CFD resolution;
- time-tag interval;
- energy offset;
- energy-bin resolution;
- singles bit layout;
- sorter version and expected word layout.

### Coincidence selection

- 511, 307, and 202 keV windows;
- prompt and delayed timing windows;
- delay offset;
- sector-difference threshold;
- ring-difference threshold;
- LOR/FOV validation;
- multiplicity policy.

### Histogram alignment

- all input files must use the same acquisition geometry;
- all detector pairs must use the same zero-based or one-based convention at each interface;
- histogram files must be unified before joint reconstruction;
- blank and phantom acquisition durations must be reflected in `scale_factor`.

### Reconstruction

- image dimensions;
- voxel size;
- image origin;
- cylinder mask;
- crystal LUT;
- normalization factors;
- attenuation units, normally `mm^-1` in the supplied reconstruction scripts;
- CUDA library name and header;
- TOF bin width and time-difference sign convention;
- initial activity and attenuation maps;
- output checkpoint naming.

---

## 17. Important Consistency Requirements

The following interfaces are particularly sensitive to mismatches:

1. **ROOT crystal hierarchy to global crystal ID**  
   The crystal numbering generated during ROOT extraction must agree with the LUT and the sorter.

2. **Singles `.dat` to `.reference` alignment**  
   Exactly one reference triplet must exist for every singles word with tag `0x82`, in the same order.

3. **Coincidence `.dat` to `.consist` alignment**  
   The two files must contain the same number of records and must be written in lockstep.

4. **CDF to `.timeDiff` alignment**  
   Each CDF record must have exactly one corresponding time-difference value.

5. **Unified histogram row alignment**  
   Blank, transmission, additive, and emission histogram rows must refer to identical normalized detector pairs.

6. **TOF sign convention**  
   The sign used when writing `.timeDiff` must match the sign expected by the TOF binning and projection code.

7. **Attenuation units**  
   The projection wrappers and attenuation mapping assume a consistent unit system. Changing from `mm^-1` to another convention requires changing every related forward-projection, initialization, and energy-mapping calculation.

---

## 18. Known Limitations and Areas for Improvement

This repository is an evolving research implementation. The following limitations should be considered:

- many scripts contain hard-coded paths and experiment-specific filenames;
- some C/C++ source files require restoration of standard include directives if copied from rendered HTML/XML content;
- error handling for short or inconsistent binary files should be expanded in production use;
- some helper functions currently return zero additive terms and are placeholders for measured or modeled scatter/random estimates;
- detector geometry and file-format assumptions are not yet centralized in a single configuration system;
- several implementations exist for similar routines in different directories, so changes should be synchronized deliberately;
- the sorter and reconstruction stages require exact compatibility of crystal indexing, LUT layout, and geometry parameters;
- the supplied programs are primarily designed for Linux/HPC environments and may require adaptation on other platforms.

Recommended future improvements include:

- introducing a single machine-readable configuration file for geometry and data paths;
- adding format-validation utilities for every binary file type;
- adding unit tests for bit packing/unpacking and ID conversion;
- documenting the exact sorter configuration grammar;
- replacing experiment-specific filenames with command-line arguments or structured configuration;
- adding reproducible environment and build instructions;
- adding small synthetic regression datasets;
- recording acquisition duration and scale factors in metadata files.

---

## 19. Scientific Purpose

The overall scientific goal is to evaluate whether internal Lu-176 radioactivity can provide useful transmission information for attenuation correction in PET, reducing or eliminating dependence on an external CT scan.

The project progressively develops this idea:

1. establish a reliable simulation-to-binary-data pipeline;
2. preserve simulation truth labels for quantitative event analysis;
3. separate emission, transmission, true, non-true, and backscatter events;
4. obtain an initial attenuation estimate with MLTR;
5. recover attenuation from emission data alone with MLAA;
6. improve attenuation estimation by combining 511 keV emission and multi-energy Lu-176 transmission data;
7. improve robustness using stabilized derivatives, Fisher-information weighting, and trust-region SCS updates.

The repository therefore combines data engineering, detector-event classification, statistical reconstruction, numerical optimization, and GPU-accelerated tomography in one end-to-end PET simulation workflow.

---

## 20. License

This project is released under the **Apache License 2.0**. See the repository license file for the full license text.

---

## 21. Citation and Acknowledgment

If this repository or its processing concepts are used in academic work, please describe:

- the GATE simulation configuration;
- detector geometry and crystal numbering;
- singles and coincidence selection criteria;
- Lu-176 transmission energy windows;
- histogram scaling and unification procedures;
- reconstruction algorithm and parameter settings;
- the exact version of the code and normalization/LUT files.

Because reconstruction results depend strongly on detector geometry, normalization data, acquisition-duration scaling, and event-selection thresholds, reproducibility requires preserving those configuration files together with the generated data and output maps.
