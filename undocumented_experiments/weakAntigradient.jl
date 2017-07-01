include("../ZernikeSuite.jl")
using ZernikeSuite
using PyPlot

maxdeg = 10
a = -0.559
fullBasis = [ZernikePoly(a, m, deg-m) for deg in 0:maxdeg for m in 0:deg]
dxFullBasis = [dx(f) for f in fullBasis]
dyFullBasis = [dy(f) for f in fullBasis]
KMat = zeros(Complex128, length(fullBasis)-1, length(fullBasis)-1)
for i in 1:length(fullBasis)-1
    for j in 1:length(fullBasis)-1
	KMat[i,j] = wip(dxFullBasis[j+1],dxFullBasis[i+1]) + wip(dyFullBasis[j+1],dyFullBasis[i+1])
    end
end
LMat = zeros(Complex128, length(fullBasis)-1, 2*(length(fullBasis)-(2*maxdeg+1)))
for i in 1:length(fullBasis)-1
    for j = 1:length(fullBasis)-(2*maxdeg+1)
	LMat[i,2*j-1] = wip(fullBasis[j],dxFullBasis[i+1])
	LMat[i,2*j+0] = wip(fullBasis[j],dyFullBasis[i+1])
    end
end

fig = figure()
GalerkinTickos = [(deg+1)*(deg+2)÷2-3/2 for deg in 1:maxdeg-1]
FunctionalTickos = [2*(deg+1)*(deg+2)÷2-1/2 for deg in 0:(maxdeg-2)-1]
ax1 = subplot(1, 4, 1)
spy(KMat)
grid("on")
axis("tight")
ax2 = subplot(1, 4, 2, sharex=ax1, sharey=ax1)
spy(inv(KMat), precision=1e-14)
grid("on")
axis("tight")
ax3 = subplot(1, 2, 2, sharey=ax1)
spy(LMat, marker="o", markersize=3)
grid("on")
axis("tight")
ax1[:yaxis][:set_ticks](GalerkinTickos)
ax1[:yaxis][:set_ticklabels]([])
ax1[:xaxis][:set_ticks](GalerkinTickos)
ax1[:xaxis][:set_ticklabels]([])
ax3[:xaxis][:set_ticks](FunctionalTickos)
ax3[:xaxis][:set_ticklabels]([])
