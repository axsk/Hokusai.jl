# Hokusai.jl

A julia package for clustering spatial timeseries, like for example eye-tracking data, based on PCCA+.

## Installation
Requires Julia 0.4. Once the repository is public one can install the Package with
```julia
Pkg.clone("giturl")
```
Until then use
```julia
include("src/Hokusai.jl")
```

## Basic usage
```julia
cluster(data, n; tau=50, sigma=100, precluster=0, sort=:size, method=:scaling)
```

- `data` is a DataFrame containing the x-/y-coordinates in the first two columns, the relative time in the third and an identifier for aggregation in the fourth column, allowing clustering of multiple test executions.
- `n` denotes the desired number of clusters.
- `tau` denotes the timestep used for creating the markov model
- `sigma` denotes the variance used for the "spatial coupling".
- `precluster` specifies the number of preclusters to generate via kmeans, to improve performance. (use `0` for no preclustering)
- `sort` denotes how the final clusters will be ordered.
`:size` sorts by number of fixations in each cluster, and `:x` by the average horizontal position.
- `method` specifies the objective function used in the PCCA+ optimization.
Accepted values: `:scaling` (Weber), `:metastability` (Deuflhard), `:crispness` (Röblitz)

`cluster` returns a `HokusaiResult`, containing the resulting clustering in the field `assignments`.


## todo
- :kmeans as method
- fix interactive plotting
- parameter optimization
- cleanup history, move to Hokusai.jl repo (remove/keep partial data?)
- think about pi
- automatic tau/gamma
- automatic n (which method in nonreversible case?)
- write tests

## how to choose tau/sigma:
- we want tau to measure every fixation, -> min fixdur?
- sigma should seperate different fixations -> mean fixdist?
