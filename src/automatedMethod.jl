# automated method
# person: 0- all corresponding subjects, else choose a subj from 1:length(corresponding subjects)
# mirrored: 0- not mirrored, 1 - mirrored
# startPos: 0-all, 52-left, 69-right
# plot: 0- just output derived number of clusters, 1- plot and print the computed steps between and the final image plus memberships
# pc: 0- output is derived number of clusters, 1- output also P, pi and chi needed for the calculation of Pc (projected transition matrix)
# cr: 1- output includes nc and the corresponding crispness; take care that id does not hold pc=1=cr
function automatedNumber(image, person, sigma, tau, mirrored, startPos, plot::Int=1, pc::Int=0, cr::Int=0)
    # get time series, transition matrix and stationary distribution of filtered data
    timeseries, ts, P, pi = prepare(image, mirrored, startPos, person, sigma, tau, plot)
    # get sorted ev gaps
    maxgaps, distances, d = getEvGaps(P)
    if plot == 1  println(string("maximal gaps: ",maxgaps)) end
    # algorithm is not working for number of clusters <= 2
    deleteat!(maxgaps, findall(maxgaps .==1))
    deleteat!(maxgaps, findall(maxgaps .==2))
    kmin = minimum(maxgaps[1:3]) # lower bound for n
    kmax = maximum(maxgaps[1:3]) # upper bound for n
    # calculate minchi values
    minchi = getMinchi(kmin, kmax, P, pi)
    maxminchi = trunc.(Int,minchi[sortperm(minchi[:,2],rev=true),1])
    if plot == 1  println(string("maximal minchis: ",maxminchi)) end
    # preselection of possible n from previous spectrral gap and minchi criteria
    possiblecluster = union(maxminchi[1:2],maxgaps[1:2])
    if plot == 1  println(string("posssible clusters: ",possiblecluster)) end
    # derive optimal n based on the crispness criterion of the preselected possible clusters
    crisp = getCrispness(possiblecluster, P, pi)
    if plot == 1 println(string("crispness criterion: ",crisp)) end
    # choose number of clusters with highest crispness
    noCluster = Int(crisp[findmax(crisp[:,2])[2],1])
    if plot == 1 println(string("final number of clusters: ",noCluster)) end
    # plot gaps of eigenvalues and respective minchis
    if plot == 1
        makeEVplot(distances, d)
        makeMCplot(minchi)
        # run pcca with derived number of clusters and print clustered image and corresponding membership vectors
        ass, chi = run(ts, noCluster, sigma, tau, person)
        plotimg(image)
        makeMembershipPlot(ass, chi)
    end
    if pc == 1 # needs plot=1
        return(noCluster, P, pi, chi)
    elseif cr == 1
        return(noCluster, crisp)
    else
        return(noCluster)
    end
end

# get time series, transition matrix and stationary distribution of filtered data
function prepare(image, mirrored, startPos, person, sigma, tau, plot = 1, symmetrize = false)
    data = filterdata(DATA, image, mirrored, startPos)
    if plot == 1
        if person != 0
            println(string("subject: ",unique(data[:subj])[person]))
        else
            println("all")
        end
    end
    # transform filtered data into a TimeSeries to corresponding person
    ts = TimeSeries(data)
    timeseries = person == 0 ? ts : ts[person]
        # get transition matrix
    grid = points(timeseries)
    P = transitionmatrix(timeseries, sigma, tau, grid, symmetrize)
    # get stationary distribution
    pi = stationaryDistr(P)

    timeseries, ts, P, pi
end

# get sorted schur values of transition matrix P
function getEvGaps(P)
    # get schur decomposition and calculate absolute values of schur values
    S = schur(P)
    real = []
    imaginary = []
    for x in S.values
        try
            append!(real, x.re)
        catch
            append!(real, x)
        end
        try
            append!(imaginary, x.im)
        catch
            append!(imaginary, 0)
        end
    end
    distances = sqrt.(real.^2 .+ imaginary.^2)
    # calculate maximal gaps of schur values
    d = [distances[i-1] - distances[i] for i = 2:length(distances)]
    maxgaps = sortperm(d,rev=true)[1:5]
    return(maxgaps, distances, d)
end

# calculate the minchi value for all number of clusters between kmin and kmax
function getMinchi(kmin, kmax, P, pi)
    minchi = zeros(kmax - kmin + 1, 2)
    for k = kmin:kmax
        # schurfactorization
        X, λ = schurvectors(P, pi, k)
        A = guessinit(X)
        chi = X*A
        minchi[k - kmin + 1, 2] = findmin(chi)[1]
        minchi[k - kmin + 1, 1] = k
    end
    return(minchi)
end

# calculate the crispness on set of ints of possible number of clusters
function getCrispness(possiblecluster, P, pi)
    crisp = zeros(length(possiblecluster), 2)
    i = 1
    for k in possiblecluster
        # schurfactorization
        X, λ = schurvectors(P, pi, k)
        A = guessinit(X)
        A = opt(A, X, A -> I3(A))
        # trace(S)/k
        crisp[i, 2] = I3(A)/k
        crisp[i, 1] = k
        i += 1
    end
    return(crisp)
end

# plot membership vectors
function makeMembershipPlot(ass, chi)
    println(string("assignements: ",ass))
    # plot memership vectors
    N,k = size(chi)
    PyPlot.figure()
    for i = 1:k
        PyPlot.plot(1:N,chi[:,i], label = i)
        PyPlot.legend()
    end
    PyPlot.xlabel("fixation points", fontsize = 15)
    PyPlot.ylabel("membership", fontsize = 15)
end

# calculate projected transition matrix
function calculatePc(chi, pi, P)
    return(inv(chi'*diagm(0 => pi)*chi)*(chi'*diagm(0 => pi)*P*chi))
end

# compare the different criteria for determining n
# kmin, kmax: interval for number of clusters
function compareCriteria(image, person, sigma, tau, mirrored, startPos, kmin, kmax, plot::Int=1, pc::Int=0)
    # get time series, transition matrix and stationary distribution of filtered data
    timeseries, ts, P, pi = prepare(image, mirrored, startPos, person, sigma, tau)
    # get sorted ev gaps
    maxgaps, distances, d = getEvGaps(P)
    # initialize criteria
    crisp = zeros(kmax - kmin + 1, 2)
    minchi = zeros(kmax - kmin + 1, 2)
    minPc = zeros(kmax - kmin + 1, 2)
    for k = kmin:kmax
        # schurfactorization
        X, λ = schurvectors(P, pi, k)
        # compute A via prechosen fct obj
        A = guessinit(X)
        chi = X*A
        #println("minChi: ",findmin(chi))
        minchi[k - kmin + 1, 2] = findmin(chi)[1]
        minchi[k - kmin + 1, 1] = k
        A = Hokusai.opt(A, X, A -> I3(A))
        # crispness = trace(S)/k
        crisp[k - kmin + 1, 2] = I3(A)/k
        crisp[k - kmin + 1, 1] = k

        Pc = inv(chi'*diagm(0=>pi)*chi)*(chi'*diagm(0=>pi)*P*chi)
        minPc[k - kmin + 1, 2] = findmin(Pc)[1]
        minPc[k - kmin + 1, 1] = k
    end
    # make plots
    makeEVplot(distances, d)
    makeMCplot(minchi)
    makeCrispPlot(crisp)
    makeMPplot(minPc)
end

# compute the average number of clusters
function meanNoCluster(image, mirrored, startPos, sigma, tau)
    # derive how many subjects there are
    data = filterdata(DATA, image, mirrored, startPos)
    nbSubjects = length(unique(data[:subj]))
    println(string("number of subjects: ", nbSubjects))
    # calculate number of clusters per subject
    nbC = DataFrame( person = String[], noCluster = Int[] )
    for person = 1:nbSubjects
        push!(nbC,[unique(data[:subj])[person], automatedNumber(image, person, sigma, tau, mirrored, startPos, 0)])
    end
    return(round(Statistics.mean(nbC[:noCluster]),digits = 2))
end

# compute number of clusters per sigma between sigmin*10 and sigmax*10 and plot in addition the corresponding crispness
function ncPerSigma(image, person, tau, mirrored, startPos, sigmin, sigmax)
    clustSig = zeros(sigmax-sigmin+1,3)
    for sigma = sigmin:sigmax
        clustSig[sigma-sigmin+1,1] = sigma * 10
        noCluster, crisp = automatedNumber(image, person, sigma*10, tau, mirrored, startPos, 0, 0, 1)
        clustSig[sigma-sigmin+1,2] = Int(crisp[findmax(crisp[:,2])[2],1])
        clustSig[sigma-sigmin+1,3] = crisp[findmax(crisp[:,2])[2],2]
    end
    PyPlot.figure()
    PyPlot.scatter(clustSig[:,1], clustSig[:,2])
    PyPlot.xlabel("sigma")
    PyPlot.ylabel("number of clusters")
    PyPlot.title("number clusters per sigma")

    PyPlot.figure()
    PyPlot.scatter(clustSig[:,1], clustSig[:,3])
    PyPlot.xlabel("sigma")
    PyPlot.ylabel("crispness")
    PyPlot.title("crispness per sigma")

    return(clustSig)
end

# compute number of clusters per tau between taumin*10 and taumax*10 and plot in addition the corresponding crispness
function ncPerTau(image, person, sigma, mirrored, startPos, taumin, taumax)
    clustTau = zeros(taumax-taumin+1,3)
    for tau = taumin:taumax
        clustTau[tau-taumin+1,1] = tau * 10
        noCluster, crisp = automatedNumber(image, person, sigma, tau*10, mirrored, startPos, 0, 0, 1)
        clustTau[tau-taumin+1,2] = Int(crisp[findmax(crisp[:,2])[2],1])
        clustTau[tau-taumin+1,3] = crisp[findmax(crisp[:,2])[2],2]
    end
    PyPlot.figure()
    PyPlot.scatter(clustTau[:,1], clustTau[:,2])
    PyPlot.xlabel("tau", fontsize=15)
    PyPlot.ylabel("number of clusters", fontsize=15)
    PyPlot.title("number clusters per tau", fontsize=15)

    PyPlot.figure()
    PyPlot.scatter(clustTau[:,1], clustTau[:,3])
    PyPlot.xlabel("tau", fontsize=15)
    PyPlot.ylabel("crispness", fontsize=15)
    PyPlot.title("crispness per tau", fontsize=15)

    return(clustTau)
end

function makeEVplot(distances, d)
    PyPlot.figure()
    PyPlot.scatter(collect(1:length(distances)), distances)
    PyPlot.title(string("largest gap after ",findmax(d)[2]," clusters"))
    PyPlot.ylabel("normed eigenvalues")
    PyPlot.xlabel("number of eigenvalue")
end

function makeMCplot(minchi)
    PyPlot.figure()
    PyPlot.scatter(minchi[:,1], minchi[:,2])
    PyPlot.xlabel("number of clusters")
    PyPlot.ylabel("minchi")
    PyPlot.xticks(minchi[:,1])
    PyPlot.title("minchi")
end

function makeMPplot(minPc)
    PyPlot.figure()
    PyPlot.scatter(minPc[:,1], minPc[:,2])
    PyPlot.xlabel("number of clusters", fontsize=15)
    PyPlot.ylabel("minPc", fontsize=15)
    PyPlot.xticks(minPc[:,1])
    PyPlot.title("minPc", fontsize=20)
end

function makeCrispPlot(crisp)
    PyPlot.figure()
    PyPlot.scatter(crisp[:,1], crisp[:,2])
    PyPlot.xlabel("number of clusters", fontsize=15)
    PyPlot.xticks(crisp[:,1])
    PyPlot.ylabel("crispness", fontsize=15)
    PyPlot.title("crispness", fontsize=20)
end
