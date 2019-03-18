#=
CHOICES:
- how do we combine multiple timeseries:
  [x] countmatrix |> mean |> rownormalize (weights by size of timeseries)
  [ ] countmatrix |> rownormalize |> mean (weights every proband the same
      - TODO: wouldnt ^^ be better?
        - from a bayesian perspective weighting every proband the same corresponds to a uniform prior on the proband choice
        - question: how will this treat rare events/outliers
- when do we symmetrize:
  [x] countmatrix |> symmetrize |> rownormalize (has interpretation)
  [ ] countmatrix |> rownormalize |> symmetrize (lacks interpretation? => seems just wrong)
=#

struct TimeSeries
    times
    points
end

points(ts::TimeSeries) = ts.points
points(tss::Vector{TimeSeries}) = vcat((points(ts) for ts in tss)...)

times(ts::TimeSeries) = ts.times
times(tss::Vector{TimeSeries}) = vcat((times(ts) for ts in tss)...)

"high level api interface
supports multiple timeseries, preclustering, sorting, symmetrization"
function cluster(ts::Union{TimeSeries, Vector{TimeSeries}}, n, sigma, tau; precluster=0, sort=:size, method=:scaling, symmetrize=false)
    method == :kmeans && return Clustering.kmeans(points(ts)', n).assignments

    grid = points(ts)
    if precluster > 0
        km = Clustering.kmeans(grid'|>copy, precluster)
        grid = km.centers' |> copy
        kmass = km.assignments
    end

    P = transitionmatrix(ts, sigma, tau, grid, symmetrize)

    pccapResult = pccap(P, n, method=method)
    ass = pccapResult.assignments
    if precluster > 0
        ass = ass[kmass]
    end
    #ass = sortcluster(ts, ass, sort)
    return ass, pccapResult.chi
end

"given a countmatrix, compute the transitionmatrix"
transitionmatrix(C::Matrix, symmetrize) =
    rownormalize(symmetrize ? C + C' : C)

transitionmatrix(ts::Union{TimeSeries, Vector{TimeSeries}}, sigma, tau, grid, symmetrize) =
    transitionmatrix(countmatrix(ts, sigma, tau, grid), symmetrize)

"transition matrix for a single timeseries"
function countmatrix(ts::TimeSeries, sigma, tau, grid::Array)
    n = size(grid, 1) # number of gridpoints
    m = getGaussMembership(ts.points, grid, sigma) # TODO: dont know anymore why i computed this as batch and not just per fixation in the loop below...

    timeframes = div.(ts.times, tau)

    last = 1
    inds = [1]
    repeats = [1]
    for i = 2:length(timeframes)
        steps = timeframes[i] - timeframes[last]
        if steps == 0
            # still the same frame
            continue
        elseif steps > 1
            # frames skipped, count self-transitions
            push!(inds, last)
            push!(repeats, steps-1)
        end
        push!(inds, i)
        push!(repeats, 1)
        last = i
    end
    m[inds[1:end-1],:]' .* repeats[2:end]' * m[inds[2:end],:]
end

countmatrix(tss::Vector{TimeSeries}, sigma, tau, grid) =
    sum(countmatrix(ts, sigma, tau, grid) for ts in tss)

function getGaussMembership(fixations, centers, sigma)
    sqdist = pairwise(SqEuclidean(), fixations', centers')
    gausskernels = exp.(-1/(2*sigma^2) * sqdist)
    rownormalize(gausskernels)
end

rownormalize(M) = M ./ sum(M, dims=2)

function sortcluster(ts, ass, method)
    if method == :none
        return ass
    elseif method == :size
        criterion = [-count(ass.==c) for c=1:maximum(ass)]
    elseif method == :x
        # TODO: implement
    elseif method == :time
        # TODO: implement
    end

    p = invperm(sortperm(criterion)) # TODO: check if invperm belongs here
    ass = [p[a] for a in ass]
end
