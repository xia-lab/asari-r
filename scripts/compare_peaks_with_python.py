"""用随机合成信号比较 SciPy find_peaks 与 R 私有峰查找器。"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

# Put the Matplotlib cache in a temporary directory to avoid writing to the user's home directory when importing asari.
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "asari-matplotlib")
)

from asari.peaks import stats_detect_elution_peaks


# The random seed is fixed to ensure that each validation covers the exact same 200 signals.
random = np.random.default_rng(20260718)
signals: list[np.ndarray] = []
for _ in range(200):
    number_of_scans = int(random.integers(30, 151))
    scans = np.arange(number_of_scans, dtype=float)
    signal = random.uniform(0, 20, number_of_scans)

    # Each signal is superimposed with 1 to 5 Gaussian peaks of different heights and widths.
    for _ in range(int(random.integers(1, 6))):
        center = random.uniform(3, number_of_scans - 4)
        sigma = random.uniform(1.2, 8)
        height = random.uniform(100, 5000)
        signal += height * np.exp(-((scans - center) ** 2) / (2 * sigma**2))

    # Uses integer signals while covering the most common data types for asari mass tracks.
    signals.append(signal.astype(np.int64))


# Encode input line by line to avoid relying on jsonlite which is not installed on the R side.
encoded_signals = "\n".join(
    ",".join(str(int(value)) for value in signal) for signal in signals
) + "\n"

package_root = Path(__file__).resolve().parents[1]
r_code = f"""
source({str(package_root / 'R/chromatograms.R')!r})
source({str(package_root / 'R/peaks.R')!r})
lines <- readLines(file(\"stdin\"))
encode <- function(values) paste(format(values, digits=17, scientific=FALSE, trim=TRUE), collapse=\",")
for (line in lines) {{
  values <- as.numeric(strsplit(line, \",\", fixed=TRUE)[[1L]])
  found <- .find_peaks_scipy(values, height=50, distance=3, prominence=30, width=1.5, wlen=25)
  fields <- c(
    encode(found$peaks),
    encode(found$properties$peak_heights),
    encode(found$properties$prominences),
    encode(found$properties$left_bases),
    encode(found$properties$right_bases),
    encode(found$properties$widths)
  )
  cat(paste(fields, collapse=\"|\"), \"\\n\", sep=\"\")
}}
"""


# Start R once and pass in all signals to reduce verification overhead caused by cross-language process startup.
completed = subprocess.run(
    ["Rscript", "-e", r_code],
    input=encoded_signals,
    text=True,
    capture_output=True,
    check=True,
)
r_lines = completed.stdout.splitlines()
if len(r_lines) != len(signals):
    raise AssertionError(
        f"R 返回 {len(r_lines)} 行，预期 {len(signals)} 行。\n{completed.stderr}"
    )


def parse_field(field: str) -> np.ndarray:
    """把 R 的逗号字段还原成数值向量。"""
    if field == "":
        return np.array([], dtype=float)
    return np.asarray([float(value) for value in field.split(",")])


# Discrete position and continuous peak attributes are compared on an item-by-item basis for each signal.
for case_number, (signal, r_line) in enumerate(zip(signals, r_lines), start=1):
    peaks, properties = find_peaks(
        signal,
        height=50,
        distance=3,
        prominence=30,
        width=1.5,
        wlen=25,
    )
    r_fields = [parse_field(field) for field in r_line.split("|")]
    python_fields = [
        peaks,
        properties["peak_heights"],
        properties["prominences"],
        properties["left_bases"],
        properties["right_bases"],
        properties["widths"],
    ]

    # Discrete peak and base positions must be identical; continuous attributes allow floating point rounding errors.
    for field_number, (python_values, r_values) in enumerate(
        zip(python_fields, r_fields), start=1
    ):
        if field_number in (1, 4, 5):
            matched = np.array_equal(python_values, r_values)
        else:
            matched = np.allclose(python_values, r_values, rtol=1e-12, atol=1e-12)
        if not matched:
            raise AssertionError(
                f"第 {case_number} 条信号的第 {field_number} 个字段不一致："
                f"Python={python_values!r}, R={r_values!r}"
            )

print(f"Python/R find_peaks 对照通过：{len(signals)} 条随机合成信号。")


# Then generate 50 signals close to the actual mass track magnitude to verify the complete statistical detection main process.
full_signals: list[np.ndarray] = []
for _ in range(50):
    number_of_scans = int(random.integers(101, 201))
    scans = np.arange(number_of_scans, dtype=float)
    signal = np.full(number_of_scans, 1000.0)

    # Peak centers were placed away from the boundaries, focusing on checking the consistency of detection, fit, area, and SNR.
    center = random.uniform(30, number_of_scans - 31)
    sigma = random.uniform(3, 8)
    height = random.uniform(150000, 500000)
    signal += height * np.exp(-((scans - center) ** 2) / (2 * sigma**2))
    full_signals.append(signal.astype(np.int64))

full_input = "\n".join(
    ",".join(str(int(value)) for value in signal) for signal in full_signals
) + "\n"
full_r_code = f"""
source({str(package_root / 'R/chromatograms.R')!r})
source({str(package_root / 'R/peaks.R')!r})
parameters <- list(
  min_peak_height=100000,
  min_peak_ratio=0.001,
  min_timepoints=6,
  min_intensity_threshold=1000,
  wlen=25,
  signal_noise_ratio=2,
  gaussian_shape=0.5,
  min_prominence_threshold=33000
)
lines <- readLines(file(\"stdin\"))
for (case_number in seq_along(lines)) {{
  intensity <- as.integer(strsplit(lines[[case_number]], \",\", fixed=TRUE)[[1L]])
  track <- list(
    id_number=case_number-1L,
    mz=100+case_number-1L,
    rt_scan_numbers=seq_along(intensity)-1L,
    intensity=intensity
  )
  peaks <- stats_detect_elution_peaks(list(track, length(intensity), parameters))
  if (length(peaks) == 0L) {{
    cat(\"NONE\\n\")
  }} else {{
    records <- vapply(peaks, function(peak) {{
      paste(format(c(
        peak$apex, peak$peak_area, peak$height, peak$left_base,
        peak$right_base, peak$goodness_fitting, peak$cSelectivity,
        peak$parent_masstrack_id, peak$mz, peak$snr
      ), digits=17, scientific=FALSE, trim=TRUE), collapse=\",\")
    }}, character(1))
    cat(paste(records, collapse=\";\"), \"\\n\", sep=\"\")
  }}
}}
"""

# The complete process only starts R once to avoid mistaking the process startup time for the algorithm's time consumption.
completed = subprocess.run(
    ["Rscript", "-e", full_r_code],
    input=full_input,
    text=True,
    capture_output=True,
    check=True,
)
full_r_lines = completed.stdout.splitlines()
if len(full_r_lines) != len(full_signals):
    raise AssertionError(
        f"完整流程的 R 返回 {len(full_r_lines)} 行，预期 {len(full_signals)} 行。"
    )

parameters = {
    "min_peak_height": 100000,
    "min_peak_ratio": 0.001,
    "min_timepoints": 6,
    "min_intensity_threshold": 1000,
    "wlen": 25,
    "signal_noise_ratio": 2,
    "gaussian_shape": 0.5,
    "min_prominence_threshold": 33000,
}

# Compare the ten numeric fields of the complete peak JSON.
for case_number, (signal, r_line) in enumerate(
    zip(full_signals, full_r_lines), start=1
):
    track = {
        "id_number": case_number - 1,
        "mz": 100 + case_number - 1,
        "rt_scan_numbers": np.arange(signal.size),
        "intensity": signal,
    }
    python_peaks = stats_detect_elution_peaks((track, signal.size, parameters))
    r_peaks = [] if r_line == "NONE" else [
        parse_field(record) for record in r_line.split(";")
    ]
    if len(python_peaks) != len(r_peaks):
        raise AssertionError(
            f"完整流程第 {case_number} 条信号峰数不一致："
            f"Python={len(python_peaks)}, R={len(r_peaks)}"
        )

    for peak_number, (python_peak, r_peak) in enumerate(
        zip(python_peaks, r_peaks), start=1
    ):
        python_values = np.asarray([
            python_peak["apex"],
            python_peak["peak_area"],
            python_peak["height"],
            python_peak["left_base"],
            python_peak["right_base"],
            python_peak["goodness_fitting"],
            python_peak["cSelectivity"],
            python_peak["parent_masstrack_id"],
            python_peak["mz"],
            python_peak["snr"],
        ])

        # Discrete results are required to be completely consistent, and small numerical errors across optimizers are allowed for fit scores.
        discrete_positions = [0, 1, 2, 3, 4, 7, 8, 9]
        if not np.array_equal(
            python_values[discrete_positions], r_peak[discrete_positions]
        ):
            raise AssertionError(
                f"完整流程第 {case_number} 条信号第 {peak_number} 个峰离散字段不一致："
                f"Python={python_values!r}, R={r_peak!r}"
            )
        if not np.allclose(
            python_values[[5, 6]], r_peak[[5, 6]], rtol=1e-8, atol=1e-10
        ):
            raise AssertionError(
                f"完整流程第 {case_number} 条信号第 {peak_number} 个峰质量字段不一致："
                f"Python={python_values!r}, R={r_peak!r}"
            )

print(f"Python/R 完整峰检测对照通过：{len(full_signals)} 条 mass tracks。")
