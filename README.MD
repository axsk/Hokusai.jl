# Hokusai.jl

A julia package for clustering spatial timeseries, like for example eye-tracking data, based on PCCA+.

## References
- Bachelor thesis (https://github.com/axsk/bachelor/blob/master/ba.pdf)
- Spectral Clustering for Non-Reversible Markov Chains (https://doi.org/https://doi.org/10.1007/s40314-018-0697-0)

## Installation
Requires Julia 0.6. 
- Install the package:
```julia
Pkg.clone("https://github.com/axsk/hokusai")
```
- Copy the data file "sallsac_Hokusai.seq" into the packages data directory (~/.julia/v0.6/Hokusai/data).
- Test the package:
```julia
Pkg.test("Hokusai")
```

## Basic usage
See the test file


- `data` is a DataFrame containing the x-/y-coordinates in the first two columns, the relative time in the third and an identifier for aggregation in the fourth column, allowing clustering of multiple test executions.
- `n` denotes the desired number of clusters.
- `tau` denotes the timestep used for creating the markov jump process
- `sigma` denotes the gaussian mixture standard deviation used for the "spatial coupling".
- `precluster` specifies the number of preclusters to generate via kmeans, to improve performance. (use `0` for no preclustering)
- `sort` denotes how the final clusters will be ordered.
`:size` sorts by number of fixations in each cluster, and `:x` by the average horizontal position.
- `method` specifies the objective function used in the PCCA+ optimization.
Accepted values: `:scaling` (Weber), `:metastability` (Deuflhard), `:crispness` (Röblitz)

`cluster` returns a `HokusaiResult`, containing the resulting clustering in the field `assignments`.


## todo
- think about pi
- parameter optimization
-- automatic tau/gamma
-- automatic n (which method in nonreversible case?)

## how to choose tau/sigma:
- we want tau to measure every fixation, -> min fixdur?
- sigma should seperate different fixations -> mean fixdist?
