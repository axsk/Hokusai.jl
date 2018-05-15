# test Hokusai
using Hokusai, DataFrames

data=DataFrame(x=[1, 2, 1, 2], y=[0, 0, 0, 0], time=[20, 20, 10, 10], groupby=[0, 0, 0, 0])

W=Hokusai.makeW(data, 1/100, [1 0;2 0])
#assert(W==[-0.1 0.1; ])

data = loaddata()
pccacluster!(data)
plot(data)

# test PCCAP
#using PCCAP

#P = readdlm("C:/dev/src/hokusai/ZUSE/P.txt", ',')
#pi = vec(readdlm("C:/dev/src/hokusai/ZUSE/pi.txt", ','))
#mlchi = readdlm("C:/dev/src/hokusai/ZUSE/chi.txt", ',')
#chi = genpcca(P, pi, 8)
#assert(norm(chi-mlchi,Inf)<1.e-9)
#(chi,mlchi) 

#@test_approx_eq X'*diagm(pi)*X eye(n)
