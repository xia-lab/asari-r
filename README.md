# asariR

`asariR` is an R implementation of [Python asari](https://github.com/shuzhao-li-lab/asari) 1.17.0.
It reads centroided mzML files, builds mass tracks, detects chromatographic peaks,
and generates feature tables for downstream statistical analysis.

Main workflow:

```text
mzML files → mass tracks → peak detection → feature tables
```

## Current status

- The R package can be installed, and its core LC-MS workflow can process real mzML data and export results.
- In one real-data file used for a strict comparison, R recovered all 2,234 preferred features reported by Python,
  but R reported four additional features. The two implementations therefore cannot yet be described as identical.
- The example file took approximately 10 seconds to process in R on the current development computer; actual runtime depends on the data and computer.
- Public R functions are available for LC, GC, DIMS, MS2, QC, and offline result visualization.

## Input Requirements

The input must be a centroided `.mzML` file.

```r
input_file <- "/absolute/path/to/sample.mzML"
file.exists(input_file)
```

The result must be `TRUE`. Do not use a nonexistent path, a directory path, or a zero-byte file.

## Installation

First install mzML reading dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("mzR")
install.packages(c("jsonlite", "remotes"))
```

Then install asariR from GitHub:

```r
remotes::install_github(
  "xia-lab/asari-r",
  upgrade = "never"
)
```

Load the package:

```r
library(asariR)
packageVersion("asariR")
```

## First Run

Only `input_file` and `output_dir` need to be changed:

```r
library(asariR)

input_file <- "/absolute/path/to/sample.mzML"
output_dir <- "/absolute/path/to/asariR-results"

stopifnot(file.exists(input_file))

result <- asari_process(
  input = input_file,
  output = output_dir,
  project_name = "my_first_project",
  mode = "pos",
  ppm = 5,
  multicores = 4,
  rt_align = FALSE,
  database_mode = "memory"
)

tables <- asari_read_results(result, table = "both")
View(tables$preferred)
```

For a single file, use `rt_align = FALSE`. For multiple samples, set `input` to the directory containing the mzML files
and use `rt_align = TRUE`.

## Where Results Are Saved

After processing, R prints the project directory and the full paths to both feature tables.

| File | Description |
|---|---|
| `preferred_Feature_table.tsv` | Preferred features that passed quality filtering; usually the primary result table |
| `export/full_Feature_table.tsv` | Complete feature list, including additional lower-quality candidate peaks |
| `project.json` | Input files, parameters, and project information |

Previously generated projects can also be read:

```r
tables <- asari_read_results(
  "/absolute/path/to/generated-project-folder",
  table = "both"
)

View(tables$preferred)
View(tables$full)
```

Identifiers such as `F1` and `F2` are only row identifiers for the current run. When comparing Python and R,
determine whether two peaks are the same using `mz` and `rtime`, rather than comparing identifiers alone.

## Commands

A command is a word entered in the Terminal that tells the software which type of task to perform. For example:

```text
process --input sample.mzML --ppm 5
```

Here `process` is the command, `--input` and `--ppm` are parameters.

| Command | Purpose | Corresponding R function |
|---|---|---|
| `process` | Process mzML data | `asari_process()` / `asari_process_dims()` |
| `analyze` | Inspect an mzML file | `asari_analyze()` |
| `annotate` | Annotate an LC or GC feature table | `asari_annotate_lc()` / `asari_annotate_gc()` |
| `viz` | Generate offline results page | `asari_dashboard()` |
| `join` | Merge two or more feature tables | `asari_merge_feature_tables()` |
| `list_workflows` | View supported workflows | `asari_available_workflows()` |

New users will mainly need `analyze` and `process`.

### Terminal format

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- COMMAND [OPTIONS]
```

Display help:

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- --help
```

Inspect an mzML file:

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- \
  analyze --input "/absolute/path/to/sample.mzML" --mode pos --ppm 5
```

Process data:

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- \
  process \
  --input "/absolute/path/to/sample.mzML" \
  --output "/absolute/path/to/results" \
  --project "my_project" \
  --mode pos --ppm 5 --multicores 4
```

List all workflows:

```bash
Rscript -e 'asariR::asari_cli(commandArgs(TRUE))' -- list_workflows
```

Use `--help` and the help pages for the corresponding R functions below to view parameters for other commands.

## Common parameters

| Parameter | Description | Common value |
|---|---|---|
| `input` | One mzML file or a directory containing mzML files | Must be specified explicitly |
| `output` | Location for saved results | An absolute path is recommended |
| `project_name` | Project name | Any short name |
| `mode` | Ion mode | `"pos"` or `"neg"` |
| `ppm` | m/z tolerance | Commonly `5` |
| `multicores` | Number of CPU cores | Start with `4` |
| `rt_align` | Whether to perform retention-time alignment | `FALSE` for one file; `TRUE` for multiple files |
| `database_mode` | Storage mode for intermediate data | `"memory"` is recommended for most users |

If you are unsure about the parameters, first process one file using the example settings.

## Other Functions

| Purpose | R function |
|---|---|
| DIMS processing | `asari_process_dims()` |
| Export mass tracks | `asari_extract_mass_tracks()` |
| Compare feature tables | `asari_compare_features()` |
| Extract target m/z | `asari_extract_targets()` |
| LC/GC annotation | `asari_annotate_lc()` / `asari_annotate_gc()` |
| Extract, match, and search MS2 spectra | `asari_extract_ms2()` / `asari_match_ms2()` / `asari_search_ms2()` |
| QC and statistical plots | `asari_feature_qc()` / `asari_qc_report()` / `asari_pca()` |
| Offline results page | `asari_dashboard()` |
| RStudio workflow assistant | `asari_gui()` |
| Thermo RAW conversion | `asari_convert_raw()` |

View the parameters and examples for a function:

```r
?asari_process
?asari_analyze
?asari_annotate_lc
?asari_extract_ms2
?asari_dashboard
```

## Python and network requirements

| Operation | Requirement |
|---|---|
| Install the R package | Internet connection required |
| Process ordinary mzML data | No Python or runtime internet connection required |
| `database_mode = "ondisk"` | Python is required to write a compatible pickle file |
| Read a Python pickle file | Python interpreter required |
| Installing the Thermo RAW converter for the first time | Internet connection required; macOS/Linux usually also requires Mono |

When centroided mzML data are already available and `database_mode = "memory"` is used,
the core processing workflow is pure R.

## Known limitations

1. The R and Python implementations produce very similar results, but complete agreement has not yet been demonstrated across many diverse datasets.
2. Peak identifiers may differ; comparisons between versions should use m/z and retention time.
3. `asari_annotate_lc()` provides local mass/RT candidate matching and is not fully equivalent to external JMS/Khipu grouping.
4. GC and MS2 results depend on the calibration table, spectral library version, and tolerances supplied by the user.
5. Before processing a large project, use a small number of files to confirm the ion mode, ppm tolerance, and peak-detection results.

## Development Testing

`testthat` is used only to check the code; it does not participate in mzML processing.

```bash
cd "/absolute/path/to/asari-r"
Rscript scripts/run_tests.R
```

Passing the tests only means that no errors were found in the covered code. It does not by itself prove complete agreement with Python for all real datasets.

## License

MIT License.

### External resources for the full test suite

The database compatibility tests require the database resources from an existing
Python asari installation. These resources are not bundled with this repository.
When that installation is outside this checkout, set its database directory explicitly:

```sh
ASARI_R_DB_DIR="/absolute/path/to/python/site-packages/asari/db" Rscript scripts/run_tests.R
```

The optional `scripts/compare_peaks_with_python.py` comparison also requires Python
asari, NumPy and SciPy in the Python environment used to run it.
