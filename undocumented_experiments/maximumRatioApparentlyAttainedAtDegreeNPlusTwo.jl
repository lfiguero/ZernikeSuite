# For fixed N and deltaN this code seeks the maximum of
#
# ‖u - proj_N(u)‖_{H^1_w}^2
# ─────────────────────────
#       ‖u‖_{H^1_w}^2
#
# among polynomials of degree N + k, 0 ≤ k ≤ deltaN, solving the corresponding
# eigenvalue problem. Here w is the weight that corresponds to the projection
# operator. As expected, the maximum is 0 for k = 0. Interestingly, the results
# suggest that the maximum of that ratio among all of H^1_w is essentially (or
# actually) attained with k = 2 already.
include("ZernikeSuite.jl")
using ZernikeSuite

α = 4.5
N = 7
deltaN = 8

# I want to find the maximal eigenvalue of
# <p-Proj_N(p), q-Proj_N(q)>_1 = lambda * <p, q>_1
# in Π^2_{N+d}, for d in {0, ..., deltaN}
basis = [[ZernikePoly(α, i, deg-i) for i in 0:deg] for deg in 0:N+deltaN]
basis = [basis...]
A = zeros(Complex128, (length(basis), length(basis)))
B = zeros(Complex128, (length(basis), length(basis)))
for j in 1:length(basis)
    p = basis[j]
    resp = p - proj(p, N)
    dpdx = dx(p)
    dpdy = dy(p)
    drespdx = dx(resp)
    drespdy = dy(resp)
    for i in 1:length(basis)
	# I know the Gram matrices will be Hermitian so I recycle what I can
	if i < j
	    A[i,j] = conj(A[j,i])
	    B[i,j] = conj(B[j,i])
	else
	    q = basis[i]
	    resq = q - proj(q, N)
	    dqdx = dx(q)
	    dqdy = dy(q)
	    dresqdx = dx(resq)
	    dresqdy = dy(resq)
	    A[i,j] = wip(resp, resq) + wip(drespdx, dresqdx) + wip(drespdy, dresqdy)
	    B[i,j] = wip(p, q) + wip(dpdx, dqdx) + wip(dpdy, dqdy)
	end
    end
end

for i in 0:deltaN
    Np = N + i
    spaceDim = div((Np+1)*(Np+2), 2)
    ew, ev, nconv, niter, nmult, res = eigs(A[1:spaceDim,1:spaceDim], B[1:spaceDim,1:spaceDim], nev=1, which=:LM)
    display((i, ew))
end
