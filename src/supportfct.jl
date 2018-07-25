
function getTransitionmatrixAndPi(image, sigma, tau, symmetrize, ratematrix)
    # filter and prepare data
    data = filterdata(Hokusai.DATA, image)
    ts = TimeSeries(data)
    # get transition matrix
    grid = points(ts)
    P = transitionmatrix(ts, sigma, tau, grid, symmetrize)
    # get stationary distribution
    pi = Hokusai.stationaryDistr(P, ratematrix)
    P, pi
end

function runPerSubject(img, n, sigma, tau, i)
    d = filterdata(DATA, img)
    ts = TimeSeries(d)
    print(" subject",i)
    run(ts[i], n, sigma, tau; kwargs...)
    plotimg(parse(Int, img[1:1]))
    nothing
end


function getprojectionmatrix(assignments, chi)
    weight = mapslices(maximum, chi, 2) |> vec
    noClusters = maximum(assignments)
    countpmatrix = zeros(noClusters, noClusters)
    last = assignments[1]
    for step = 2:length(assignments)
        next = assignments[step]
        countpmatrix[last, next] += weight[step]
        last = next
    end
    countpmatrix ./= sum(countpmatrix,2)
    return countpmatrix
end
