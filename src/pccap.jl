using Optim
using StatsBase

type PccapResult
    assignments::Vector # discrete cluster assignment
    counts::Vector
    chi::Matrix         # fuzzy cluster assignment
end

function pccap(P::Matrix, n::Integer; pi=nothing, method=:scaling)
    if pi == nothing
        which = ispmatrix(P)? :LM : :SM
        pi = abs.(vec((Array{Float64})(eigs(P';nev=1, which=which)[2])))
        pi = pi / sum(pi)                     # => first col of X is one
    end

    X, λ = schurvectors(P, pi, n)

    A=feasible(guess(X),X)

    if     method == :scaling        obj = A -> I1(A,X)
    elseif method == :metastability  obj = A -> I2(A,λ)
    elseif method == :crispness      obj = A -> I3(A)
	else error("no valid pcca+ objective method")
    end

    n>2 && (A=opt(A, X, obj))

    chi = X*A

    assignments = vec(mapslices(indmax,chi,2))

    #counts = hist(assignments)[2] # deprecated
    counts = zeros(Int, n)
    for a in assignments
        counts[a] += 1
    end
    sum(counts.>0) != n && warn("Not all clusters were assigned")
    return PccapResult(assignments, counts, chi)
end


function ispmatrix(P)
    # decide whether P is a probability or rate matrix
    # check whether first row sum is close to one
    return abs(sum(P[1,:]) - 1) < 0.01
end

function schurvectors(P, pi, n)	
    Pw = diagm(sqrt.(pi))*P*diagm(1./sqrt.(pi)) # rescale to keep markov property
    Sw = schurfact!(Pw)                       # returns orthonormal vecs by def
    Xw, λ = selclusters!(Sw, n, ispmatrix(P))
    X  = diagm(1./sqrt.(pi)) * Xw              # scale back
    X  = X[1,1]>0 ? X : -X
    X, λ
end

# select the schurvectors corresponding to the n abs-largest eigenvalues
# if reverse==true select highest abs value, otherwise select lowest (for rate matrices)
function selclusters!(S, n, reverse)
    ind = sortperm(abs.(S[:values]), rev=reverse) # get indices for largest eigenvalues
    select = zeros(Bool, size(ind))            # create selection vector
    select[ind[1:n]] = true
    S = ordschur!(S, select)                  # reorder selected vectors to the left
    S[:vectors][:,1:n], S[:values][1:n]       # select first n vectors
end

# compute initial guess based on indexmap
guess(X) = inv(X[indexmap(X), :])

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

function feasible(A,X)
    A[:,1] = -sum(A[:,2:end], 2)
    A[1,:] = -minimum(X[:,2:end] * A[2:end,:], 1)
    A / sum(A[1,:])
end

function I1(A,X)
    sum(maximum(X*A,1))
end

function I2(A,λ)
    n = size(A,1)
    trace = 0
    for i = 1:n, j = 1:n
        trace += real(λ[i]) * A[i,j]^2 / A[1,j] # TODO: check if takin real is the right approach to complex λ
    end
    return trace
end

function I3(A)
    n = size(A,1)
    trace = 0
    for i = 1:n, j=1:n
        trace += A[i,j]^2 / A[1,j]
    end
    return trace
end

function opt(A0, X, objective)
    function transform(A)
      # cut out the fixed part
      cA = A[2:end, 2:end]
      # flatten matrix to vector for use in optimize
      return reshape(cA,prod(size(cA)))
    end

    function transformback(tA)
      # reshape back to matrix
      cA = reshape(tA, size(A0,1)-1, size(A0,2)-1)
      # unite with the fixed part
      A = A0
      A[2:end,2:end] = cA
      return A
    end

    obj(tA) = -objective(feasible(transformback(tA), X))
    result = optimize(obj, transform(A0), NelderMead())
    return feasible(transformback(result.minimizer), X)
end
