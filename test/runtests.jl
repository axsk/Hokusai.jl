using Hokusai

## PCCA+ tests

# TODO: write better tests

P=[.5 .5 0 0; .5 .5 0 0; 0 0 1 0; 0 0 0 1]
for m in [:scaling, :metastability, :crispness]
    res = Hokusai.pccap(P, 3, method = m)
    @assert res.assignments[1] == res.assignments[2]
    @assert res.assignments[1] != res.assignments[3]
    @assert res.assignments[1] != res.assignments[4]
    @assert res.assignments[3] != res.assignments[4]
end

## Hokusai tests

img = 1
n   = 5

alldata = Hokusai.readdata!()
imgdata = Hokusai.filterdata(alldata, img)
Hokusai.cluster(imgdata[[:fposx, :fposy, :fixdur, :subj]], n, precluster=100)
Hokusai.savecl(img, n, 50, 50, precl=100, folder=tempdir())
