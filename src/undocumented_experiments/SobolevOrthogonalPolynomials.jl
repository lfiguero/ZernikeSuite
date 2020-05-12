import LinearAlgebra: qr, svd, diag, Diagonal, dot, norm, rank, I
import Random: MersenneTwister

export sip, ssp, sp

function sip(f::ZFun, g::ZFun)
    @assert (f.α == g.α == 0) || abs(f.α-g.α)/min(abs(f.α),abs(g.α)) < 10*eps()
    2.0 * (wip(dzp(f),dzp(g)) + wip(dzs(f),dzs(g))) + wip(proj(f,0), proj(g,0))
end

function bf(f::ZFun, g::ZFun)
    @assert (f.α == g.α == 0) || abs(f.α-g.α)/min(abs(f.α),abs(g.α)) < 10*eps()
    ZRaise = ZernikeSuite.raise
    dzpzpf = ZRaise(dzp(dzp(f)))
    dzpzsf = ZRaise(dzp(dzs(f)))
    dzszpf = ZRaise(dzs(dzp(f)))
    dzszsf = ZRaise(dzs(dzs(f)))
    dzpzpg = ZRaise(dzp(dzp(g)))
    dzpzsg = ZRaise(dzp(dzs(g)))
    dzszpg = ZRaise(dzs(dzp(g)))
    dzszsg = ZRaise(dzs(dzs(g)))
    4.0 * (wip(dzpzpf,dzpzpg) + wip(dzpzsf,dzpzsg) + wip(dzszpf,dzszpg) + wip(dzszsf,dzszsg)) + 2.0 * (wip(dθ(dzp(f)), dθ(dzp(g))) + wip(dθ(dzs(f)), dθ(dzs(g))))
end

function recombine(coll::Array{ZernikeSuite.ZFun,1})
    n = length(coll)
    X = randn(n,n) + im*randn(n,n)
    # We will recombine using matrix Q in order to preserve intra-degree
    # orthogonality
    Q, R = qr(X)
    newcoll = Array{ZernikeSuite.ZFun}(undef, n)
    for i = 1:n
	newcoll[i] = Q[i,1]*coll[1]
	for j = 2:n
	    newcoll[i] += Q[i,j]*coll[j]
	end
    end
    newcoll
end

function SOP(α::Real, maxdeg::Integer, normalizeFlag::Bool=false, recombineFlag::Bool=false)
	@assert α > -1
	@assert maxdeg ≥ 0
	# We start with basis which is weighted L²-orthogonal
	basis = [ZernikePoly(α, i, k-i) for k in 0:maxdeg for i in 0:k]
	dim = length(basis)
	# We turn the basis into a weighted Sobolev-orthogonal one via a
	# Gram–Schmidt process
	squaredNorm = Array{Float64}(undef, dim)
	for i = 1:dim
		for j = 1:i-1
			basis[i] = basis[i] - sip(basis[i], basis[j])/squaredNorm[j] * basis[j]
		end
		if normalizeFlag
			basis[i] = basis[i]/sqrt(sip(basis[i],basis[i]))
			squaredNorm[i] = 1.0
		else
			squaredNorm[i] = real(sip(basis[i], basis[i]))
		end
	end
	# Intra-degree recombinations
	if recombineFlag
		for k = 0:maxdeg
			rng = ZernikeSuite.positionRange(k)
			basis[rng] = recombine(basis[rng])
		end
	end
	return basis
end

function SturmLiouvilleTest(α::Real, maxdeg::Integer, normalizeFlag::Bool=false, recombineFlag::Bool=false)
	basis = SOP(α, maxdeg, normalizeFlag, recombineFlag)
	A = Array{ComplexF64}(undef, (length(basis), length(basis)))
	B = Array{ComplexF64}(undef, (length(basis), length(basis)))
	for i = 1:length(basis)
		for j = 1:length(basis)
			A[i,j] = bf(basis[i], basis[j])
			B[i,j] = sip(basis[i], basis[j])
		end
	end
	bf_orthogonality_test = norm(A-Diagonal(diag(A)))/norm(A)
	sip_orthogonality_test = norm(B-Diagonal(diag(B)))/norm(B)
	eigenvalues = diag(A) ./ diag(B)
	return bf_orthogonality_test, sip_orthogonality_test, eigenvalues
end

function relativeComponentWeight(f::ZFun)
    weights = ZernikeSuite.h(f.α, f.degree)
    out = zeros(Float64, f.degree+1)
    for k = 0:f.degree
	rng = ZernikeSuite.positionRange(k)
	out[k+1] = real(dot(f.coefficients[rng], weights[rng].*f.coefficients[rng]))
    end
    out = sqrt.(out/sum(out))
end

function coefficientFormulaTest(α::Real, maxdeg::Integer)
	basis = SOP(α, maxdeg)
	w = real(-[basis[ZernikeSuite.positionRange(deg)[i+1]].coefficients[ZernikeSuite.positionRange(deg-2)[i]] for deg in 3:maxdeg for i in 1:deg-1])
	d = [deg for deg in 3:maxdeg for i in 1:deg-1]
	md = [abs(2*i-deg) for deg in 3:maxdeg for i in 1:deg-1]
	modelVal = (d-md).*(d+md)./(d-md.+2*α)./(d+md.+2*α)
	relRes = norm(w-modelVal)/norm(w)
	return relRes
end

# Collections which may span the Sobolev orthogonal polynomials which are
# Sobolev orthogonal to harmonic polynomials
function firstAugustTest(α::Real, deg::Integer)
	bVαmo = [ZernikePoly(α, i, deg-1-i) for i in 0:deg-1]
	collx = [dx(p)-mzp(mzs(dx(p)))-2*α*mx(p) for p in bVαmo]
	colly = [dy(p)-mzp(mzs(dy(p)))-2*α*my(p) for p in bVαmo]
	[collx;colly]
end

mbump(f::ZFun) = f - mzp(mzs(f))

# Test of a possibly useful identity
function tenthAugustTest(α::Real, m::Integer, n::Integer)
	obj = ZernikeSuite.lower(ZernikePoly(α+1, m, n))
	xsqobj = mzp(mzs(obj))
	dzsobj = dzs(obj)
	dzpobj = dzp(obj)
	oflapobj = dzp(dzsobj)
	lhs1 = -(α+1)*((m+1)*xsqobj + obj - xsqobj) + m*mbump(mzs(dzsobj)) - (α+1)*mbump(mzp(dzpobj)) + mbump(mbump(oflapobj))
	lhs2 = -(α+1)*((n+1)*xsqobj + obj - xsqobj) + n*mbump(mzp(dzpobj)) - (α+1)*mbump(mzs(dzsobj)) + mbump(mbump(oflapobj))
	rhs1 = -(m+1)*(α+1)*ZernikePoly(α, m, n)
	rhs2 = -(n+1)*(α+1)*ZernikePoly(α, m, n)
	res1 = rhs1-lhs1
	res2 = rhs2-lhs2
	relErr1 = sqrt(real(wip(res1,res1))/real(wip(rhs1,rhs1)))
	relErr2 = sqrt(real(wip(res2,res2))/real(wip(rhs2,rhs2)))
	return relErr1, relErr2
end

# Test of another possibly useful identity
function seventeenthAugustTest(α::Real, m::Integer, n::Integer)
	obj = ZernikeSuite.lower(ZernikePoly(α+1, m, n))
	lhs1 = -m*(α+1)*mzp(mzs(obj)) - (α+2)*mbump(mzp(dzp(obj))) + m*mbump(obj + mzs(dzs(obj))) + mbump(mbump(dzs(dzp(obj))))
	rhs1 = -m*(α+1)*ZernikePoly((α+1)-1, m, n)
	res1 = rhs1 - lhs1
	relErr1 = wip(res1,res1)==wip(rhs1,rhs1)==0.0 ? 0.0 : sqrt(real(wip(res1,res1))/real(wip(rhs1,rhs1)))
	lhs2 = -n*(α+1)*mzp(mzs(obj)) - (α+2)*mbump(mzs(dzs(obj))) + n*mbump(obj + mzp(dzp(obj))) + mbump(mbump(dzs(dzp(obj))))
	rhs2 = -n*(α+1)*ZernikePoly((α+1)-1, m, n)
	res2 = rhs2 - lhs2
	relErr2 = wip(res2,res2)==wip(rhs2,rhs2)==0.0 ? 0.0 : sqrt(real(wip(res2,res2))/real(wip(rhs2,rhs2)))
	return relErr1, relErr2
end

# Yet another test of a possibly useful identity
function eighteenthAugustTest(α::Real, m::Integer, n::Integer)
	obj = ZernikeSuite.lower(ZernikePoly(α+1, m, n))
	lhs1 = (α+1)*mzs(mzp(obj)) + (m+α+1)*mbump(obj) - mbump(mzp(dzp(obj)))
	rhs1 = (α+1)*ZernikePoly((α+1)-1, m, n)
	res1 = rhs1 - lhs1
	relErr1 = sqrt(real(wip(res1,res1))/real(wip(rhs1,rhs1)))
	lhs2 = (α+1)*mzs(mzp(obj)) + (n+α+1)*mbump(obj) - mbump(mzs(dzs(obj)))
	rhs2 = (α+1)*ZernikePoly((α+1)-1, m, n)
	res2 = rhs2 - lhs2
	relErr2 = sqrt(real(wip(res2,res2))/real(wip(rhs2,rhs2)))
	return relErr1, relErr2
end

# More tests
function twentysecondAugustTest(α::Real, m::Integer, n::Integer)
	obj = ZernikeSuite.lower(ZernikePoly(α+1, m-1, n-1))
	lhs = -2*α*(α+1)*mzp(mzs(obj)) + (2*m*n+α*(m+n))*mbump(obj) + α*mbump(mzs(dzs(obj))+mzp(dzp(obj)))
	Q = ZernikePoly((α+1)-1, m, n) - m*n/(m+α)/(n+α)*ZernikePoly((α+1)-1, m-1, n-1)
	rhs = -2*(α+1)*(m+α)*(n+α)/(m+n+α)*Q
	res = rhs - lhs
	relErr = sqrt(real(wip(res,res))/real(wip(rhs,rhs)))
	return lhs, rhs, relErr
end

function OPSOPTestsA(α::Real, maxdeg::Integer)
	OPBasis = [ZernikeSuite.lower(ZernikePoly(α+1, i, maxdeg-i)) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	SOPp1 = [-2*2*(α+1)*mbump(p) for p in OP]
	SOPp2 = [4*α*(α+1)*mzs(mzp(p)) for p in OP]
	SOPp3 = [-4*(α+1)*mbump(mzs(dzs(p)) + mzp(dzp(p))) for p in OP]
	SOPp4 = [4*mbump(mbump(dzs(dzp(p)))) for p in OP]
	SOP = SOPp1 + SOPp2 + SOPp3 + SOPp4
	lowerBasis = [ZernikePoly((α+1)-1, i, deg-i) for deg in 0:maxdeg+1 for i in 0:deg]
	SOTestp1 = [sip(q, Mp) for q in lowerBasis, Mp in SOPp1]
	SOTestp2 = [sip(q, Mp) for q in lowerBasis, Mp in SOPp2]
	SOTestp3 = [sip(q, Mp) for q in lowerBasis, Mp in SOPp3]
	SOTestp4 = [sip(q, Mp) for q in lowerBasis, Mp in SOPp4]
	SOTest = [sip(q, Mp) for q in lowerBasis, Mp in SOP]
	return SOTest
end

function matrixKernelStudy(mat::Matrix{ComplexF64})
	r = rank(mat)
	nrow = size(mat,1)
	ncol = size(mat,2)
	ppVt = svd(mat).Vt[r+1:ncol,:]
	makesones = ppVt*ones(ncol,1)
	aux = [makesones [zeros(1,ncol-r-1); I]]
	q = qr(aux).Q
	nonObviousKernelBasis = ppVt\q[:,2:end]
	tf = (true, false)
	allCombinations = Iterators.product(Iterators.repeated(tf, ncol)...)
	coll = []
	for combination in allCombinations
		isThisCero = zeros(nrow)
		for j in 1:ncol
			if combination[j]
				isThisCero += mat[:,j]
			end
		end
		if norm(isThisCero) < 1e-10
			push!(coll, combination)
		end
	end
	return nonObviousKernelBasis, coll
end

function OPSOPTestsB(α::Real, maxdeg::Integer)
	OPBasis = [ZernikeSuite.lower(ZernikePoly(α+1, i, maxdeg-i)) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	x = mx(ZFun((α+1)-1, 0, [1.0]))
	y = mx(ZFun((α+1)-1, 0, [1.0]))
	gSOPp1 = [4*2*(α+1)*[mx(p), my(p)] for p in OP]
	gSOPp2 = [8*α*(α+1)*[mx(p), my(p)] for p in OP]
	gSOPp3 = [-2*2*(α+1)*[mbump(dx(p)), mbump(dy(p))] for p in OP]
	gSOPp4 = [4*α*(α+1)*[mzs(mzp(dx(p))), mzs(mzp(dy(p)))] for p in OP]
	gSOPp5 = [8*(α+1)*[mx(mzs(dzs(p))+mzp(dzp(p))), my(mzs(dzs(p))+mzp(dzp(p)))] for p in OP]
	gSOPp6 = [-4*(α+1)*[mbump(dx(p)), mbump(dy(p))] for p in OP]
	gSOPp7 = [-4*(α+1)*[mbump(mx(dx(dx(p)))+my(dy(dx(p)))), mbump(mx(dx(dy(p)))+my(dy(dy(p))))] for p in OP]
	gSOPp8 = [-4*[mx(mbump(4*dzs(dzp(p)))), my(mbump(4*dzs(dzp(p))))] for p in OP]
	gSOPp9 = [[mbump(mbump(dx(4*dzs(dzp(p))))), mbump(mbump(dy(4*dzs(dzp(p)))))] for p in OP]
	lowerBasis = [ZernikePoly((α+1)-1, i, deg-i) for deg in 0:maxdeg for i in 0:deg]
	vLowerBasis = [[q1, q2] for q1 in lowerBasis, q2 in lowerBasis]
	SOTestp1 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp1]
	SOTestp2 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp2]
	SOTestp3 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp3]
	SOTestp4 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp4]
	SOTestp5 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp5]
	SOTestp6 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp6]
	SOTestp7 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp7]
	SOTestp8 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp8]
	SOTestp9 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp9]
	# I found that this matrix has rank 4 = 9-5
	return [SOTestp1[:] SOTestp2[:] SOTestp3[:] SOTestp4[:] SOTestp5[:] SOTestp6[:] SOTestp7[:] SOTestp8[:] SOTestp9[:]]
	# I found that this matrix has rank 4 = 7-3
	#return [SOTestp1[:]+SOTestp2[:] SOTestp3[:]+SOTestp6[:] SOTestp4[:] SOTestp5[:] SOTestp7[:] SOTestp8[:] SOTestp9[:]]
end

function OPSOPTestsC(α::Real, maxdeg::Integer)
	OPBasis = [ZernikeSuite.lower(ZernikePoly(α+1, i, maxdeg-i)) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	x = mx(ZFun((α+1)-1, 0, [1.0]))
	y = mx(ZFun((α+1)-1, 0, [1.0]))
	gSOPp1 = [-(2*2*α+4)*[mbump(dx(p)), mbump(dy(p))] for p in OP]
	gSOPp2 = [4*α*(α+1)*[mzs(mzp(dx(p))), mzs(mzp(dy(p)))] for p in OP]
	gSOPp3 = [8*(α+1)*[mx(mzs(dzs(p))+mzp(dzp(p))), my(mzs(dzs(p))+mzp(dzp(p)))] for p in OP]
	gSOPp4 = [-4*(α+1)*[mbump(mx(dx(dx(p)))+my(dy(dx(p)))), mbump(mx(dx(dy(p)))+my(dy(dy(p))))] for p in OP]
	gSOPp5 = [-4*[mx(mbump(4*dzs(dzp(p)))), my(mbump(4*dzs(dzp(p))))] for p in OP]
	gSOPp6 = [[mbump(mbump(dx(4*dzs(dzp(p))))), mbump(mbump(dy(4*dzs(dzp(p)))))] for p in OP]
	lowerBasis = [ZernikePoly((α+1)-1, i, deg-i) for deg in 0:maxdeg for i in 0:deg]
	vLowerBasis = [[q1, q2] for q1 in lowerBasis, q2 in lowerBasis]
	SOTestp1 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp1]
	SOTestp2 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp2]
	SOTestp3 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp3]
	SOTestp4 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp4]
	SOTestp5 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp5]
	SOTestp6 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp6]
	# I found that this matrix has rank 4 = 6-2
	return [SOTestp1[:] SOTestp2[:] SOTestp3[:] SOTestp4[:] SOTestp5[:] SOTestp6[:]]
end

function OPSOPTestsC4(α::Real, maxdeg::Integer)
	OPBasis = [ZernikeSuite.lower(ZernikePoly(α+1, i, maxdeg-i)) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	x = mx(ZFun((α+1)-1, 0, [1.0]))
	y = mx(ZFun((α+1)-1, 0, [1.0]))
	gSOPp1 = [-(2*2*α+4)*[mbump(dx(p)), mbump(dy(p))] for p in OP]
	gSOPp2 = [4*α*(α+1)*[mzs(mzp(dx(p))), mzs(mzp(dy(p)))] for p in OP]
	gSOPp3 = [8*(α+1)*[mx(mzs(dzs(p))+mzp(dzp(p))), my(mzs(dzs(p))+mzp(dzp(p)))] for p in OP]
	gSOPp5 = [-4*[mx(mbump(4*dzs(dzp(p)))), my(mbump(4*dzs(dzp(p))))] for p in OP]
	gSOPp6 = [[mbump(mbump(dx(4*dzs(dzp(p))))), mbump(mbump(dy(4*dzs(dzp(p)))))] for p in OP]
	lowerBasis = [ZernikePoly((α+1)-1, i, deg-i) for deg in 0:maxdeg for i in 0:deg]
	vLowerBasis = [[q1, q2] for q1 in lowerBasis, q2 in lowerBasis]
	SOTestp1 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp1]
	SOTestp2 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp2]
	SOTestp3 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp3]
	IBP4a = [-8*(α+1)^2*(wip(vq[1], mzs(mzp(dx(p)))) + wip(vq[2], mzs(mzp(dy(p))))) for vq in vLowerBasis, p in OP]
	IBP4b = [4*(α+1)*(wip(mx(dx(vq[1]))+my(dy(vq[1])), mbump(dx(p))) + wip(mx(dx(vq[2]))+my(dy(vq[2])), mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	IBP4c = [4*2*(α+1)*(wip(vq[1], mbump(dx(p))) + wip(vq[2], mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	SOTestp5 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp5]
	SOTestp6 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp6]
	return [SOTestp1[:] SOTestp2[:] SOTestp3[:] IBP4a[:] IBP4b[:] IBP4c[:] SOTestp5[:] SOTestp6[:]]
end

function OPSOPTestsC5(α::Real, maxdeg::Integer)
	OPBasis = [ZernikeSuite.lower(ZernikePoly(α+1, i, maxdeg-i)) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	x = mx(ZFun((α+1)-1, 0, [1.0]))
	y = mx(ZFun((α+1)-1, 0, [1.0]))
	gSOPp1 = [-(2*2*α+4)*[mbump(dx(p)), mbump(dy(p))] for p in OP]
	gSOPp2 = [4*α*(α+1)*[mzs(mzp(dx(p))), mzs(mzp(dy(p)))] for p in OP]
	gSOPp3 = [8*(α+1)*[mx(mzs(dzs(p))+mzp(dzp(p))), my(mzs(dzs(p))+mzp(dzp(p)))] for p in OP]
	gSOPp4 = [-4*(α+1)*[mbump(mx(dx(dx(p)))+my(dy(dx(p)))), mbump(mx(dx(dy(p)))+my(dy(dy(p))))] for p in OP]
	gSOPp6 = [[mbump(mbump(dx(4*dzs(dzp(p))))), mbump(mbump(dy(4*dzs(dzp(p)))))] for p in OP]
	lowerBasis = [ZernikePoly((α+1)-1, i, deg-i) for deg in 0:maxdeg for i in 0:deg]
	vLowerBasis = [[q1, q2] for q1 in lowerBasis, q2 in lowerBasis]
	SOTestp1 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp1]
	SOTestp2 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp2]
	SOTestp3 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp3]
	SOTestp4 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp4]
	IBP5a = [-8*(α+1)*wip(mx(vq[1])+my(vq[2]), mx(dx(p))+my(dy(p))) for vq in vLowerBasis, p in OP]
	IBP5b = [4*(wip(mx(dx(vq[1]))+my(dx(vq[2])), mbump(dx(p))) + wip(mx(dy(vq[1]))+my(dy(vq[2])), mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	IBP5c = [4*(wip(vq[1], mbump(dx(p))) + wip(vq[2], mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	SOTestp6 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp6]
	#                                 THIS ONE              THIS ONE
	return [SOTestp1[:] SOTestp2[:] SOTestp3[:] SOTestp4[:] IBP5a[:] IBP5b[:] IBP5c[:] SOTestp6[:]]
end

function OPSOPTestsC6(α::Real, maxdeg::Integer)
	OPBasis = [ZernikeSuite.lower(ZernikePoly(α+1, i, maxdeg-i)) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	x = mx(ZFun((α+1)-1, 0, [1.0]))
	y = mx(ZFun((α+1)-1, 0, [1.0]))
	gSOPp1 = [-(2*2*α+4)*[mbump(dx(p)), mbump(dy(p))] for p in OP]
	gSOPp2 = [4*α*(α+1)*[mzs(mzp(dx(p))), mzs(mzp(dy(p)))] for p in OP]
	gSOPp3 = [8*(α+1)*[mx(mzs(dzs(p))+mzp(dzp(p))), my(mzs(dzs(p))+mzp(dzp(p)))] for p in OP]
	gSOPp4 = [-4*(α+1)*[mbump(mx(dx(dx(p)))+my(dy(dx(p)))), mbump(mx(dx(dy(p)))+my(dy(dy(p))))] for p in OP]
	gSOPp5 = [-4*[mx(mbump(4*dzs(dzp(p)))), my(mbump(4*dzs(dzp(p))))] for p in OP]
	lowerBasis = [ZernikePoly((α+1)-1, i, deg-i) for deg in 0:maxdeg for i in 0:deg]
	vLowerBasis = [[q1, q2] for q1 in lowerBasis, q2 in lowerBasis]
	SOTestp1 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp1]
	SOTestp2 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp2]
	SOTestp3 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp3]
	SOTestp4 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp4]
	SOTestp5 = [wip(vq[1], vp[1])+wip(vq[2], vp[2]) for vq in vLowerBasis, vp in gSOPp5]
	IBP6a = [-2*2*(α+2)*(wip(vq[1], mbump(dx(p)))+wip(vq[2], mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	IBP6b = [4*(α+1)*(α+2)*(wip(vq[1],mzs(mzp(dx(p))))+wip(vq[2],mzs(mzp(dy(p))))) for vq in vLowerBasis, p in OP]
	IBP6c = [-4*(α+2)*(wip(mx(dx(vq[1]))+my(dy(vq[1])),mbump(dx(p)))+wip(mx(dx(vq[2]))+my(dy(vq[2])),mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	IBP6d = [wip(4*dzs(dzp(vq[1])),mbump(mbump(dx(p))))+wip(4*dzs(dzp(vq[2])),mbump(mbump(dy(p)))) for vq in vLowerBasis, p in OP]
	#                                                                                              vanishes
	return [SOTestp1[:] SOTestp2[:] SOTestp3[:] SOTestp4[:] SOTestp5[:] IBP6a[:] IBP6b[:] IBP6c[:] IBP6d[:]]
end

function fifteenthSeptemberTest(α::Real, maxdeg::Integer, j::Integer)
	OPBasis = [ZernikePoly(α, i, maxdeg-i) for i in 0:maxdeg]
	rng = MersenneTwister(0) # For reproducibility
	mat = randn(rng,maxdeg+1,maxdeg+1) + im*randn(rng,maxdeg+1,maxdeg+1)
	OP = [sum(mat[i,:].*OPBasis) for i in 1:maxdeg+1]
	@assert j==1 || j==2
	mj = (j==1) ? mx : my
	dj = (j==1) ? dx : dy
	ve = [0;0]; ve[j] = 1.0
	lowerBasis = [ZernikePoly(α, i, deg-i) for deg in 0:maxdeg-1 for i in 0:deg]
	vLowerBasis = [[q1, q2] for q1 in lowerBasis, q2 in lowerBasis]
	term1 = [-2*(wip(vq[1], mx(dj(p)))+wip(vq[2], my(dj(p)))) for vq in vLowerBasis, p in OP]
	term2 = [wip(vq[1],mbump(dx(dj(p))))+wip(vq[2],mbump(dy(dj(p)))) for vq in vLowerBasis, p in OP]
	term3 = [-2*α*(wip(vq[1], ve[1]*p)+wip(vq[2], ve[2]*p)) for vq in vLowerBasis, p in OP]
	term4 = [-2*α*(wip(vq[1], mj(dx(p)))+wip(vq[2], mj(dy(p)))) for vq in vLowerBasis, p in OP]
	mappedOP = [mbump(dj(p))-2*α*mj(p) for p in OP]
	full = [wip(vq[1], dx(P)) + wip(vq[2], dy(P)) for vq in vLowerBasis, P in mappedOP]
	return [term1[:] term2[:] term3[:] term4[:]], full
end

# First order parameter lowering operator
function foplo(f::ZFun, b::Real, j::Integer)
	@assert j==1 || j==2
	mj = (j==1) ? mx : my
	dj = (j==1) ? dx : dy
	return -mbump(dj(f)) + 2*(b+1)*mj(f)
end

function optosop(f::ZFun)
	α = f.α-1.0
	return ZernikeSuite.lower(foplo(foplo(f, α, 1), α-1, 1) + foplo(foplo(f, α, 2), α-1, 2))
end

# Second order Sobolev semi-inner product
function sossip(f::ZFun, g::ZFun)
	@assert (f.α == g.α == 0) || abs(f.α-g.α)/min(abs(f.α),abs(g.α)) < 10*eps()
	dzpzpf = dzp(dzp(f))
	dzpzsf = dzp(dzs(f))
	dzszpf = dzs(dzp(f))
	dzszsf = dzs(dzs(f))
	dzpzpg = dzp(dzp(g))
	dzpzsg = dzp(dzs(g))
	dzszpg = dzs(dzp(g))
	dzszsg = dzs(dzs(g))
	4.0 * (wip(dzpzpf,dzpzpg) + wip(dzpzsf,dzpzsg) + wip(dzszpf,dzszpg) + wip(dzszsf,dzszsg))
end

# Second order old Sobolev inner product
function soosip(f::ZFun, g::ZFun)
	semi = sossip(f,g)
	semi + wip(proj(f,1),proj(g,1))
end

# Second order new Sobolev inner product
function sonsip(f::ZFun, g::ZFun)
	semi = sossip(f,g)
	semi + 2.0 * (wip(proj(dzp(f),0),proj(dzp(g),0)) + wip(proj(dzs(f),0),proj(dzs(g),0))) + wip(proj(f,0),proj(g,0))
end

# Generation of second order old Sobolev orthogonal polynomials through the Gram–Schmidt process
function sooSOP(α::Real, maxdeg::Integer, normalizeFlag::Bool=false, recombineFlag::Bool=false)
	@assert α > -1
	@assert maxdeg ≥ 0
	# We start with basis which is weighted L²-orthogonal
	basis = [ZernikePoly(α, i, k-i) for k in 0:maxdeg for i in 0:k]
	dim = length(basis)
	# We turn the basis into a weighted Sobolev-orthogonal one via a
	# Gram–Schmidt process
	squaredNorm = Array{Float64}(undef, dim)
	for i = 1:dim
		for j = 1:i-1
			basis[i] = basis[i] - soosip(basis[i], basis[j])/squaredNorm[j] * basis[j]
		end
		if normalizeFlag
			basis[i] = basis[i]/sqrt(soosip(basis[i],basis[i]))
			squaredNorm[i] = 1.0
		else
			squaredNorm[i] = real(soosip(basis[i], basis[i]))
		end
	end
	# Intra-degree recombinations
	if recombineFlag
		for k = 0:maxdeg
			rng = ZernikeSuite.positionRange(k)
			basis[rng] = recombine(basis[rng])
		end
	end
	return basis
end


# Generation of second order new Sobolev orthogonal polynomials through the Gram–Schmidt process
function sonSOP(α::Real, maxdeg::Integer, normalizeFlag::Bool=false, recombineFlag::Bool=false)
	@assert α > -1
	@assert maxdeg ≥ 0
	# We start with basis which is weighted L²-orthogonal
	basis = [ZernikePoly(α, i, k-i) for k in 0:maxdeg for i in 0:k]
	dim = length(basis)
	# We turn the basis into a weighted Sobolev-orthogonal one via a
	# Gram–Schmidt process
	squaredNorm = Array{Float64}(undef, dim)
	for i = 1:dim
		for j = 1:i-1
			basis[i] = basis[i] - sonsip(basis[i], basis[j])/squaredNorm[j] * basis[j]
		end
		if normalizeFlag
			basis[i] = basis[i]/sqrt(sonsip(basis[i],basis[i]))
			squaredNorm[i] = 1.0
		else
			squaredNorm[i] = real(sonsip(basis[i], basis[i]))
		end
	end
	# Intra-degree recombinations
	if recombineFlag
		for k = 0:maxdeg
			rng = ZernikeSuite.positionRange(k)
			basis[rng] = recombine(basis[rng])
		end
	end
	return basis
end

# Sobolev slice projection
function ssp(f::ZFun, d)
	if d < 0 || d > f.degree
		return ZFun(f.α, 0, [0.0])
	end
	b = SOP(f.α, d)[positionRange(d)]
	mat = [sip(v,u) for u in b, v in b]
	vec = [sip(f,v) for v in b]
	c = mat\vec
	out = ZFun(f.α, d, zeros(ComplexF64, polyDim(d)))
	for i = 1:length(b)
		out = out + c[i]*b[i]
	end
	out
end

# Sobolev projection
function sp(f::ZFun, d)
	if d < 0
		return ZFun(f.α, 0, [0.0])
	end
	if d ≥ f.degree
		return f
	end
	b = SOP(f.α, d)
	mat = [sip(v,u) for u in b, v in b]
	vec = [sip(f,v) for v in b]
	c = mat\vec
	out = ZFun(f.α, d, zeros(ComplexF64, polyDim(d)))
	for i = 1:length(b)
		out = out + c[i]*b[i]
	end
	out
end

# Strong Sturm–Louville operator for SOP
function sslsop(f::ZFun)
	α = f.α
	main = -mbump(dx(dx(f))+dy(dy(f))) + 2*α*(mx(dx(f))+my(dy(f))) - dθ(dθ(f))
	K = 2*proj(mx(dx(f))+my(dy(f)),0)
	main + K
end
