#!/usr/bin/env julia
# From Korg-benchmarks: julia --project=benchmark --threads=10 misc/benchmark_worktrees.jl
# Each worktree runs in a fresh Julia process, with one BLAS thread.
# Each case gets an excluded warmup and one measured synthesis, as in the upstream suite.
# Line-list loading is excluded. Allocated MiB is cumulative allocation, not peak RAM.

using BenchmarkTools, LinearAlgebra, Pkg, Printf, Statistics

if isempty(ARGS)
    repo = dirname(@__DIR__)
    repos = dirname(repo)
    worktrees = [
        ("Stock", joinpath(repos, "Korg-contributions")),
        ("PR1", joinpath(repos, "Korg-pr1")),
        ("PR2", joinpath(repos, "Korg-pr2")),
        ("PR1+PR2", joinpath(repos, "Korg-both")),
    ]

    # Reuse the installed dependency versions in a temporary environment. Only Korg's
    # source path changes; the worktrees and their own environments remain untouched.
    ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
    rows = String[]
    mktempdir() do environment
        cp(joinpath(repo, "benchmark", "Project.toml"), joinpath(environment, "Project.toml"))
        cp(joinpath(repo, "Manifest.toml"), joinpath(environment, "Manifest.toml"))
        Pkg.activate(environment; io=devnull)
        for (name, path) in worktrees
            println(stderr, "Benchmarking $name: $path")
            Pkg.develop(; path, preserve=Pkg.PRESERVE_ALL, io=devnull)
            command = `$(Base.julia_cmd()) --startup-file=no --project=$environment
                       --threads=$(Threads.nthreads()) --heap-size-hint=12G
                       $(@__FILE__) --worker $name $path`
            push!(rows, read(command, String))
        end
    end

    @printf("\n%-10s %-9s %-8s %12s %14s %14s\n",
            "Version", "Star", "List", "Time (s)", "Allocated MiB", "Allocations")
    foreach(print, rows)
elseif length(ARGS) == 3 && ARGS[1] == "--worker"
    using Korg
    name, path = ARGS[2:3]
    @assert realpath(pkgdir(Korg)) == realpath(path) "Wrong Korg checkout loaded"
    BLAS.set_num_threads(1)

    benchmark_params = [
        ("sun", (Teff=5777, logg=4.44, M_H=0.0)),
        ("M dwarf", (Teff=3500, logg=4.6, M_H=0.0)),
    ]
    # Windows round the supplied line-centre extrema outward to multiples of 100 Å.
    wavelength_ranges_and_linelists = [
        ("GES", (4700.0, 9000.0, 0.01), Korg.get_GES_linelist(; include_molecules=true)),
        ("APOGEE", (15000.0, 17600.0, 0.01),
         Korg.get_APOGEE_DR17_linelist(; include_water=true)),
    ]

    suite = BenchmarkGroup()
    for (star, params) in benchmark_params
        suite[star] = BenchmarkGroup()
        for (list, wavelengths, linelist) in wavelength_ranges_and_linelists
            suite[star][list] = @benchmarkable Korg.synth(
                ; Teff=$(params.Teff), logg=$(params.logg), M_H=$(params.M_H),
                linelist=$linelist, wavelengths=$wavelengths) setup=(GC.gc()) seconds=300 samples=3 evals=1
        end
    end

    # Run in the same order in every process. Each run performs one excluded warmup.
    for (star, _) in benchmark_params
        for (list, _, _) in wavelength_ranges_and_linelists
            println(stderr, "  $name / $star / $list")
            trial = run(suite[star][list])
            @printf("%-10s %-9s %-8s %12.4f %14.2f %14d\n",
                    name, star, list, median(trial).time / 1e9,
                    BenchmarkTools.memory(trial) / 1024^2, BenchmarkTools.allocs(trial))
        end
    end
else
    error("Run without arguments; edit the worktrees and case lists above to change the benchmark.")
end
