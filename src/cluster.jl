# TODO: open questions:
# - how do we combine multiple timeseries:
# a countmatrix |> mean |> rownormalize (weights by size of timeseries)
# b countmatrix |> rownormalize |> mean (bad for outliers/rare points)
# - when do we symmetrize
# c countmatrix |> symmetrize |> rownormalize (has interpretation)
# d countmatrix |> rownormalize |> symmetrize (lacks interpretation?)

# the code below implements options a, d

struct TimeSeries
    times
    points
end

points(ts::TimeSeries) = ts.points
points(tss::Vector{TimeSeries}) = vcat((points(ts) for ts in tss)...)

"high level api interface
supports multiple timeseries, preclustering, sorting, symmetrization"
function cluster(ts::Union{TimeSeries, Vector{TimeSeries}}, n, sigma, tau; precluster=0, sort=:size, method=:scaling, symmetrize=false)
    method == :kmeans && return Clustering.kmeans(points(ts)', n).assignments

    grid = points(ts)
    if precluster > 0
        km = Clustering.kmeans(grid', precluster)
        grid = km.centers'
        kmass = km.assignments
    end

    C = countmatrix(ts, sigma, tau, grid)
    if symmetrize
        C = (C + C') / 2
    end
    P = rownormalize(C)

    ass = pccap(P, n, method=method).assignments
    if precluster > 0
        ass = ass[kmass]
    end
    ass = sortcluster(ts, ass, sort)
    return ass
end

"transition matrix for a single timeseries"
function countmatrix(ts::TimeSeries, sigma, tau, grid::Array)
    n = size(grid, 1)
    C = zeros(n, n)
    m = getGaussMembership(ts.points, grid, sigma) # TODO: dont know anymore why i computed this as batch and not just per fixation in the loop below...

    timeframes = div.(ts.times, sigma)

    last = 1
    for i = 1:length(timeframes)
        steps = timeframes[i] - timeframes[last]
        if steps == 0
            # still the same frame
            continue
        elseif steps > 1
            # frames skipped, count self-transitions
            C += m[last,:] * m[last,:]' * (steps - 1)
        end
        C += m[last,:] * m[i,:]'
        last = i
    end
    C
end

countmatrix(tss::Vector{TimeSeries}, sigma, tau, grid) =
    sum(countmatrix(ts, sigma, tau, grid) for ts in tss)

function getGaussMembership(fixations, centers, sigma)
    sqdist = pairwise(SqEuclidean(), fixations', centers')
    gausskernels = exp.(-1/(2*sigma^2) * sqdist)
    rownormalize(gausskernels)
end

rownormalize(M) = M ./ sum(M,2)

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
