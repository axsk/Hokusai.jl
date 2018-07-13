# TODO: take a look at the whole feasiblization/optimization routine

struct PccapResult
    assignments::Vector # discrete cluster assignment
    counts::Vector
    chi::Matrix         # fuzzy cluster assignment
end

function pccap(P::Matrix, n::Integer; pi=nothing, method=:scaling, ratematrix=false)
    if pi == nothing
        which = ratematrix ? :SM : :LM
        pi = eigs(P', nev=1, which=which)[2]
        @assert isreal(pi)
        pi = abs.(pi) |> vec
        pi = pi / sum(pi)                     # => first col of X is one
    end

    X, λ = schurvectors(P, pi, n, ratematrix)


    if     method == :scaling        obj = A -> I1(A,X)
    elseif method == :metastability  obj = A -> I2(A,λ)
    elseif method == :crispness      obj = A -> I3(A)
    else error("no valid pcca+ objective method")
    end

    # TODO: why n>2, what if n=2? is guessinit() already the optimium?
    A = guessinit(X)
    n>2 && (A=opt(A, X, obj))

    chi = X*A

    assignments = mapslices(indmax,chi,2) |> vec

    counts = zeros(Int, n)
    for a in assignments
        counts[a] += 1
    end
    sum(counts.>0) != n && warn("Not all clusters were assigned")
    return PccapResult(assignments, counts, chi)
end

function schurvectors(P, pi, n, ratematrix)
    Pw = diagm(sqrt.(pi))*P*diagm(1./sqrt.(pi)) # rescale to keep markov property
    Sw = schurfact!(Pw)                       # returns orthonormal vecs by def
    Xw, λ = selclusters!(Sw, n, ratematrix)
    X  = diagm(1./sqrt.(pi)) * Xw              # scale back
    X  = X[1,1]>0 ? X : -X
    X, λ
end

# select the schurvectors corresponding to the n abs-largest eigenvalues
# if reverse==true select highest abs value, otherwise select lowest (for rate matrices)
function selclusters!(S, n, ratematrix)
    ind = sortperm(abs.(S[:values]), rev=!ratematrix) # get indices for dominant eigenvalues
    select = zeros(Bool, size(ind))           # create selection vector
    select[ind[1:n]] = true
    S = ordschur!(S, select)                  # reorder selected vectors to the left
    if S.T[n+1, n] != 0                       # check if we are cutting along a schur block
        warn("conjugated eigenvector missing")
    end
    S[:vectors][:,1:n], S[:values][1:n]       # select first n vectors
end

# compute initial guess based on indexmap
guessinit(X) = feasiblize!(inv(X[indexmap(X), :]), X)

function indexmap(X)
    # get indices of rows of X to span the largest simplex
    rnorm(x) = sqrt.(sum(abs2.(x),2)) |> vec
    ind=zeros(Int, size(X,2))
    for j in 1:length(ind)
        rownorm=rnorm(X)
        # store largest row index
        ind[j]=indmax(rownorm)
        if j == 1
            # translate to origin
            X = X .- X[ind[1],:]'
        else
            # remove subspace
            X=X/rownorm[ind[j]]
            vt=X[ind[j],:]'
            X=X-X*vt'*vt
        end
    end
    return ind
end

function feasiblize!(A,X)
    A[:,1] = -sum(A[:,2:end], 2)
    A[1,:] = -minimum(X[:,2:end] * A[2:end,:], 1)
    A / sum(A[1,:])
end

# maximal scaling condition (source?)
function I1(A,X)
    sum(maximum(X*A,1))
end

# metastability criterion, cf. Deuflhard (2000)
# TODO: where does this come from?
function I2(A,λ)
    n = size(A,1)
    trace = 0
    for i = 1:n, j = 1:n
        trace += real(λ[i]) * A[i,j]^2 / A[1,j] # TODO: check if taking the real part is the right approach to complex λ
    end
    return trace
end

# crispness criterion, cf. Roeblitz (2013)
# only applies if X is normalized (eq. 8)
# TODO: check does this apply?
function I3(A)
    n = size(A,1)
    trace = 0
    for i = 1:n, j=1:n
        trace += A[i,j]^2 / A[1,j]
    end
    return trace
end

function opt(A0, X, objective)
    A = copy(A0)
    Av = view(A, 2:size(A,1), 2:size(A,2)) # view on the variable part

    function obj(a)
        Av[:] = a
        -objective(feasiblize!(A, X))
    end

    result = optimize(obj, Av[:], NelderMead())
    Av[:] = result.minimizer
    return feasiblize!(A, X)
end
