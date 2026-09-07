# Module Reading Plan

Use this file as the checklist while reading Python asari and filling the R
skeleton.

## 1. chromatograms.py

Status: read once.

Questions to preserve:

- How are raw MS1 points converted into `(mz, scan, intensity)` tuples?
- Where does ppm tolerance enter the mass-track construction?
- Which filters happen before clustering?
- Which filters may create false negatives?
- Which parts are only RT alignment and not directly related to m/z tolerance?

## 2. mass_functions.py

Status: all top-level functions ported and checked against Python/SciPy behavior.

Questions:

- How does `nn_cluster_by_mz_seeds()` split candidate bins?
- How are seed m/z values chosen?
- What happens when no seed is found?
- Does it consider RT shape or only m/z distance?

## 3. peaks.py

Status: all 19 top-level functions ported. The SciPy-compatible peak-property
layer was checked on 200 deterministic synthetic signals, and the full peak
detection path was compared on 50 synthetic mass tracks. Integration with
`constructors.py` remains.

Questions:

- Where is a real chromatographic peak defined?
- How is noise estimated?
- How are peak boundaries chosen?
- What peak-shape checks are used?
