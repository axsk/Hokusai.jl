using Hokusai

## PCCA+ tests

# TODO: write better tests

P=[.5 .5 0 0; .5 .5 0 0; 0 0 1 0; 0 0 0 1]
for m in [:scaling, :metastability, :crispness]
    res = Hokusai.pccap(P, 3, method = m)
    @show res
    @assert res.assignments[1] == res.assignments[2]
    @assert res.assignments[1] != res.assignments[3]
    @assert res.assignments[1] != res.assignments[4]
    @assert res.assignments[3] != res.assignments[4]
end

Hokusai.run(1,5, 50, 50, precluster=100)
