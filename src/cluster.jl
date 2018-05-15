function getGaussMembership(fixations, centers, sigma)
    sqdist = pairwise(SqEuclidean(), fixations', centers')
    gausskernels = exp.(-sigma^(-2.) * sqdist)
    rownormalize(gausskernels)
end

rownormalize(M) = M ./ sum(M,2)

# create rate matrix
function makeW(data, sigma, grid)
    data[:index] = 1:size(data,1)
    fixations = Array(data[[:x,:y]])
    centers   = grid == :none ? fixations : grid
    gaussmemb = getGaussMembership(fixations, centers, sigma)
    n = size(centers,1)
    ratesum = zeros(Float64, n, n)
    fromsum = zeros(Float64, n)

    by(data, :groupby) do subjdata
        # user periodic boundary
        from = gaussmemb[subjdata[1:size(subjdata,1)|>collect, :index],:]
        to   = gaussmemb[subjdata[vcat(2:size(subjdata,1),1), :index],:]
        time = 1./       subjdata[1:size(subjdata,1)|>collect, :time]
        ratesum += from' * (to .* time)
        fromsum += sum(from', 2)
    end

    W = ratesum ./ fromsum

    for i=1:n
        W[i,i] = -(sum(W[i,:]) - W[i,i])
    end
    return W
end

global SymmetrizeP = false
setSymmetrizeP(b::Bool) = (SymmetrizeP = b)

function makeP(data, tau, sigma, grid, symmetrize)
    data[:index] = 1:size(data,1)
    from = (Int64)[]
    to   = (Int64)[]

    by(data, :groupby) do subjdata
        abstime = cumsum(subjdata[:time])
        offs = subjdata[1,:index] - 1
        time = 0
        iold = 1
        while (true)
            time += tau
            inew = findfirst(t-> (t>=time), abstime) # find nearest fixation after current
            inew == 0 && break
            append!(from, [iold+offs])
            append!(to,   [inew+offs])
            iold = inew
        end
    end

    fixations = Array(data[[:x,:y]])
    centers   = grid == :none ? fixations : grid
    gaussmemb = getGaussMembership(fixations, centers, sigma)

    C = gaussmemb[from,:]'*gaussmemb[to,:]
    if symmetrize
        C = (C + C') / 2
    end
    P = rownormalize(C)
end

function sortcluster!(data; sort = :size)
    methods = Dict(
                   :size=>(cl-> -size(cl,1)),
                   :x=>(cl->mean(cl[:x])))

    res = by(data, :cluster, methods[sort])
    ord = sort!(res, cols=[:x1])[:cluster]

    inv = Array{Int8}(maximum(ord))
    inv[ord]=1:size(ord,1)
    data[:cluster] = inv[data[:cluster]]
end

function precluster_grid(data, precl)
    (precl == 0) && return (:none, :none)
    precl = minimum([precl, size(data)[1]])
    km = Clustering.kmeans(Array(data[[:x, :y]])',precl)
    grid = km.centers'
    ass  = km.assignments
    grid, ass
end

function cluster(data, n; tau=50, sigma=100, precluster=0, sort=:size, method=:scaling, symmetrize=false)
    data = DataFrame(x=data[1], y=data[2], time=data[3], groupby=data[4])

    method == :kmeans && return kmeans(data, n)

    grid, kmass = precluster_grid(data, precluster)

    P = tau==0 ? makeW(data, sigma, grid) : makeP(data, tau, sigma, grid, symmetrize)

    ass = pccap(Array{Float64}(P), n, method=method).assignments
    if precluster > 0
        ass = ass[kmass]
    end

    data[:cluster] = ass
    sortcluster!(data, sort = sort)

    HokusaiResult(data, data[:cluster], n, tau, sigma, P)
end

function kmeans(data, n)
    data = DataFrame(x=data[1], y=data[2], time=data[3], groupby=data[4])
    data[:index] = 1:size(data,1)
    data[:cluster] = Clustering.kmeans(Array(data[[:x, :y]])',n).assignments
    HokusaiResult(data, data[:cluster], n, NaN, NaN, NaN)
end

