# returns the transition matrix P and its "non-reversible part" distm, i.e. P - distm is reversible
function reversibilitydist(image, person, sigma, tau, mirrored, startPos)
    # prepare and filter data
    timeseries, ts, P, pi = prepare(image, mirrored, startPos, person, sigma, tau)
    # get weighted schur decomposition
    X, L = schurvectors(P, pi, 0)
    # calculate the non-reversible part via the non-diagonal entries of the schur value matrix L
    nondiag = L - diagm(diag(L))
    distm = zeros(size(X))
    nonzeroind = findn(nondiag .!= 0)
    n = length(nonzeroind[1])
    for i = 1:n
        c = nondiag[nonzeroind[1][i],nonzeroind[2][i]]
        row = X[:,nonzeroind[2][i]]
        column = X[:,nonzeroind[1][i]]
        m = row' .* column
        distm += c*m*diagm(pi) #abs?
    end
    distm, P
end

# deviation from revrsibility per sigma
function revPerSigma(image, person, tau, mirrored, startPos, sigmin, sigmax)
    dist = zeros(sigmax-sigmin+1,2)
    for sigma = sigmin:sigmax
        distm , P = reversibilitydist(image, person, sigma*10, tau, mirrored, startPos)
        dist[sigma-1,1] = vecnorm(distm,2) # frobenius norm
        dist[sigma-1,2] = sigma*10
    end
    # plot
    PyPlot.figure()
    PyPlot.scatter(dist[:,2], dist[:,1])
    PyPlot.xlabel("sigma", fontsize=15)
    PyPlot.ylabel("reversibility", fontsize=15)
    PyPlot.title("deviation from reversibility per sigma", fontsize=15)
    return dist
end

# deviation from revrsibility per tau
function revPerTau(image, person, sigma, mirrored, startPos, taumin, taumax)
    dist = zeros(taumax-taumin+1,2)
    for tau = taumin:taumax
        distm , P = reversibilitydist(image, person, sigma, tau*10, mirrored, startPos)
        dist[tau-1,1] = vecnorm(distm,2) # frobenius norm
        dist[tau-1,2] = tau*10
    end
    # plot
    PyPlot.figure()
    PyPlot.scatter(dist[:,2], dist[:,1])
    PyPlot.xlabel("tau", fontsize=15)
    PyPlot.ylabel("reversibility", fontsize=15)
    PyPlot.title("deviation from reversibility per tau", fontsize=15)
    return dist
end
