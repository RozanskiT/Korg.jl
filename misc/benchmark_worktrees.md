# Reproducing the synthesis benchmarks

[benchmark_worktrees.jl](benchmark_worktrees.jl) compares complete `Korg.synth` calls using
stock Korg, the line-absorption optimizations (PR1), shared continuum interpolation (PR2),
and both changes together. This guide records the exact source revisions, dependency versions,
settings, and latest reported measurements.

## Source versions

| Version | Branch | Commit | Directory |
| --- | --- | --- | --- |
| Stock | `main` | `5d5a4e8624d824f5f208bf22fd3213a474368e14` | `Korg-contributions` |
| PR1 | `perf-line-absorption` | `4f0fb62528c69c697d8274060e9225eb21e66761` | `Korg-pr1` |
| PR2 | `perf-shared-continuum` | `7ba8be34362b1718a66447be1cfa7f7366fdf07f` | `Korg-pr2` |
| PR1+PR2 | `perf-both` | `3051fb7aa403a36f21ae57e572d4b18b8c54fc8a` | `Korg-both` |

PR1 contains strength hoisting, accumulation through a view, and a species generator passed to
`unique`. PR2 shares continuum interpolation intervals and weights across atmospheric layers,
while retaining callback support. The combined version merges PR2 into PR1.

GitHub stores the commits and branches; each reader creates their own worktrees. The commands
below use detached worktrees pinned to the measured commits, so later changes to branch tips do
not change the versions being compared. The fork must publish `perf-measurements` and `perf-both`
as well as the two PR branches before the fresh-clone instructions can be used remotely.

## Set up a fresh clone

Use Julia **1.12.6** to reproduce the recorded environment. Start in a directory where none of
`Korg-benchmarks`, `Korg-contributions`, `Korg-pr1`, `Korg-pr2`, or `Korg-both` already exists.

```sh
git clone --branch perf-measurements \
  https://github.com/RozanskiT/Korg.jl.git Korg-benchmarks
cd Korg-benchmarks

git worktree add --detach ../Korg-contributions \
  5d5a4e8624d824f5f208bf22fd3213a474368e14

git worktree add --detach ../Korg-pr1 \
  4f0fb62528c69c697d8274060e9225eb21e66761

git worktree add --detach ../Korg-pr2 \
  7ba8be34362b1718a66447be1cfa7f7366fdf07f

git worktree add --detach ../Korg-both \
  3051fb7aa403a36f21ae57e572d4b18b8c54fc8a
```

The directory names matter: the script finds the four checkouts beside `Korg-benchmarks`.
The parent directory can be anywhere.
If these worktrees already exist, skip cloning and creation, enter `Korg-benchmarks`, and check
that their commits match the table before running.

## Install dependencies

Run the commands below from the **`Korg-benchmarks` project root**, the directory containing
`Project.toml`, `benchmark/`, and `misc/`.

**Optional: restore the dependency versions used for the recorded measurements.** Copy
[misc/benchmark_manifest.toml](benchmark_manifest.toml) into the project root and **name the copy
`Manifest.toml`**. The original `misc/benchmark_manifest.toml` stays unchanged. This command performs
both the copy and the renaming:

```sh
cp misc/benchmark_manifest.toml ./Manifest.toml
```

This replaces any existing root `Manifest.toml`. The snapshot's Korg entry uses `path = "."`,
which refers to the directory containing the restored manifest: the project root.
The root manifest remains ignored by Git; the saved snapshot in `misc/` is tracked.

Then install the dependencies, whether or not you chose to restore the snapshot:

```sh
julia --startup-file=no --project=benchmark \
  -e 'using Pkg; Pkg.instantiate()'
```

With the saved manifest, Julia installs its recorded dependency versions. On a fresh clone where
you skipped the copy, Julia resolves compatible versions and creates a new root `Manifest.toml`.
Those versions may differ from the ones used for the recorded measurements, but the benchmark
still uses the same dependencies for all four Korg versions within a run. Instantiation may
download packages and artifacts on first use.

The saved snapshot is optional; a root `Manifest.toml` must exist before running the benchmark,
because the script copies it into its temporary environment. Either setup above provides it.

Only the benchmark checkout needs this setup. The script creates a temporary environment from
`benchmark/Project.toml` and the root manifest, then changes only Korg's source path with
`Pkg.develop(...; preserve=Pkg.PRESERVE_ALL)` for each version. Each checkout runs in a separate
Julia process. It does not modify the other worktrees or their environments.

## Run the benchmarks

Run these commands from `Korg-benchmarks`:

```sh
julia --startup-file=no --project=benchmark --threads=1 \
  misc/benchmark_worktrees.jl

julia --startup-file=no --project=benchmark --threads=4 \
  misc/benchmark_worktrees.jl

julia --startup-file=no --project=benchmark --threads=10 \
  misc/benchmark_worktrees.jl
```

Each command runs all four versions and all four synthesis cases. Let one command finish before
starting another, and avoid other CPU-intensive work during measurements. Package precompilation
may appear between versions; it is outside the reported synthesis times.

To save a run's table, append `| tee benchmark-threads-10.txt` (or the matching thread count).
Progress messages are written to stderr, and the result table to stdout. To choose fewer versions
or different cases, edit the short lists in the script; the measurements below use all entries.

## Measurement settings

| Atmosphere | Teff (K) | logg | M_H |
| --- | ---: | ---: | ---: |
| Sun | 5777 | 4.44 | 0.0 |
| M dwarf | 3500 | 4.6 | 0.0 |

| Linelist | Loader | Wavelengths (Angstrom) | Step (Angstrom) |
| --- | --- | --- | ---: |
| GES | `Korg.get_GES_linelist(; include_molecules=true)` | 4700-9000 | 0.01 |
| APOGEE DR17 | `Korg.get_APOGEE_DR17_linelist(; include_water=true)` | 15000-17600 | 0.01 |

All remaining `Korg.synth` parameters use the defaults of the recorded source versions.

- Each case requests `samples=3`, `evals=1`, and `seconds=300` from BenchmarkTools. On a slower
  machine, the time budget can stop the trial before all three samples have completed.
- `run` performs one excluded warmup before the measured samples.
- `GC.gc()` runs in the setup before each sample, outside the timed synthesis. GC occurring
  during synthesis remains part of the measurement.
- Linelists are loaded once per version before benchmarking and reused across cases and samples.
  Loading, package setup, and warmup/compilation are excluded from the reported statistics.
- Every worker uses one BLAS thread and `--heap-size-hint=12G`. The heap hint is not a hard RAM
  limit. Julia thread counts are 1, 4, and 10 in the three runs.
- Versions run in the order Stock, PR1, PR2, PR1+PR2. Each runs Sun/GES, Sun/APOGEE, M dwarf/GES,
  and M dwarf/APOGEE in that order.
- Time is `median(trial).time / 1e9`, in seconds. Allocated MiB and allocation count are
  `BenchmarkTools.memory(trial) / 1024^2` and `BenchmarkTools.allocs(trial)`. BenchmarkTools 1.8.0
  records the minimum allocated bytes and allocation count across samples, independently of the
  median time. Allocated bytes are cumulative allocation per synthesis, not peak RAM usage.

## Recorded results

This section records the latest complete rerun of all four versions, at the printed precision.
The runs were performed in the order 10, 1, then 4 Julia threads; they are arranged below by thread
count. Individual sample times were not saved, so the tables do not provide uncertainty estimates.

- Machine: Apple M5, 10 CPU cores, 32 GiB RAM, macOS.
- Julia: 1.12.6, aarch64 Apple Darwin.
- BenchmarkTools: 1.8.0; remaining package versions are in the saved manifest.
- Benchmark script revision corresponding to the measurements: `6f75543bd26596b8eda3ee4a0eed2804e972d891` on
  `perf-measurements`. The documentation update changes comments only in the script.

Results depend on the machine and system load. Compare versions measured on the same machine;
small differences need repeated measurements before being interpreted as reliable improvements.

### 1 Julia thread

| Version | Star | List | Time (s) | Allocated MiB | Allocations |
| --- | --- | --- | ---: | ---: | ---: |
| Stock | Sun | GES | 35.9087 | 6062.17 | 9141881 |
| Stock | Sun | APOGEE | 6.3504 | 3160.08 | 6320003 |
| Stock | M dwarf | GES | 73.2032 | 39808.48 | 23583195 |
| Stock | M dwarf | APOGEE | 53.5855 | 79473.03 | 12460914 |
| PR1 | Sun | GES | 29.8860 | 4315.23 | 8761228 |
| PR1 | Sun | APOGEE | 5.6535 | 1733.50 | 6215673 |
| PR1 | M dwarf | GES | 65.0752 | 5535.79 | 8946407 |
| PR1 | M dwarf | APOGEE | 54.0525 | 2642.61 | 6573308 |
| PR2 | Sun | GES | 19.0698 | 6010.92 | 9141470 |
| PR2 | Sun | APOGEE | 4.2164 | 3136.74 | 6319592 |
| PR2 | M dwarf | GES | 48.4838 | 39793.79 | 23582609 |
| PR2 | M dwarf | APOGEE | 50.8337 | 79147.37 | 12460328 |
| PR1+PR2 | Sun | GES | 13.0323 | 4311.32 | 8760817 |
| PR1+PR2 | Sun | APOGEE | 3.3908 | 1733.46 | 6215262 |
| PR1+PR2 | M dwarf | GES | 37.3886 | 5535.75 | 8945821 |
| PR1+PR2 | M dwarf | APOGEE | 46.2706 | 2327.42 | 6572722 |

### 4 Julia threads

| Version | Star | List | Time (s) | Allocated MiB | Allocations |
| --- | --- | --- | ---: | ---: | ---: |
| Stock | Sun | GES | 12.3957 | 7178.73 | 9142025 |
| Stock | Sun | APOGEE | 2.9037 | 4183.19 | 6320147 |
| Stock | M dwarf | GES | 25.2555 | 41526.20 | 23583339 |
| Stock | M dwarf | APOGEE | 18.7305 | 80281.95 | 12461058 |
| PR1 | Sun | GES | 9.8410 | 5414.91 | 8761372 |
| PR1 | Sun | APOGEE | 2.6055 | 2617.96 | 6215817 |
| PR1 | M dwarf | GES | 20.9714 | 7130.24 | 8946551 |
| PR1 | M dwarf | APOGEE | 16.6729 | 3921.91 | 6573452 |
| PR2 | Sun | GES | 7.3781 | 7110.68 | 9141614 |
| PR2 | Sun | APOGEE | 2.2306 | 3766.09 | 6319736 |
| PR2 | M dwarf | GES | 18.7132 | 41474.92 | 23582753 |
| PR2 | M dwarf | APOGEE | 18.2057 | 80562.20 | 12460472 |
| PR1+PR2 | Sun | GES | 5.5777 | 5413.66 | 8760961 |
| PR1+PR2 | Sun | APOGEE | 1.9155 | 2400.05 | 6215406 |
| PR1+PR2 | M dwarf | GES | 15.0927 | 7130.19 | 8945965 |
| PR1+PR2 | M dwarf | APOGEE | 17.0458 | 3816.82 | 6572866 |

### 10 Julia threads

| Version | Star | List | Time (s) | Allocated MiB | Allocations |
| --- | --- | --- | ---: | ---: | ---: |
| Stock | Sun | GES | 8.7654 | 9338.35 | 9142313 |
| Stock | Sun | APOGEE | 2.3785 | 5736.67 | 6320435 |
| Stock | M dwarf | GES | 18.1133 | 44745.66 | 23583627 |
| Stock | M dwarf | APOGEE | 12.9345 | 83606.32 | 12461346 |
| PR1 | Sun | GES | 7.9249 | 7620.34 | 8761660 |
| PR1 | Sun | APOGEE | 2.2644 | 4386.89 | 6216105 |
| PR1 | M dwarf | GES | 16.0106 | 10321.03 | 8946839 |
| PR1 | M dwarf | APOGEE | 12.8862 | 6900.70 | 6573740 |
| PR2 | Sun | GES | 5.6252 | 9326.99 | 9141902 |
| PR2 | Sun | APOGEE | 1.9917 | 5319.49 | 6320024 |
| PR2 | M dwarf | GES | 13.8220 | 44728.80 | 23583041 |
| PR2 | M dwarf | APOGEE | 12.4115 | 82672.47 | 12460760 |
| PR1+PR2 | Sun | GES | 4.4316 | 7620.25 | 8761249 |
| PR1+PR2 | Sun | APOGEE | 1.7911 | 3832.36 | 6215694 |
| PR1+PR2 | M dwarf | GES | 10.9259 | 10319.07 | 8946253 |
| PR1+PR2 | M dwarf | APOGEE | 11.0495 | 5219.90 | 6573154 |
