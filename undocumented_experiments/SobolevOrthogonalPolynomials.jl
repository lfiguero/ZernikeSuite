include("../ZernikeSuite.jl")
using ZernikeSuite

function sip(f::ZFun, g::ZFun)
    @assert f.α == g.α
    2.0 * (wip(dzp(f),dzp(g)) + wip(dzs(f),dzs(g))) + wip(proj(f,0), proj(g,0))
end

function bf(f::ZFun, g::ZFun)
    @assert f.α == g.α
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
    newcoll = Array{ZernikeSuite.ZFun}(n)
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
	squaredNorm = Array{Float64}(dim)
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
	A = Array{Complex128}(length(basis), length(basis))
	B = Array{Complex128}(length(basis), length(basis))
	for i = 1:length(basis)
		for j = 1:length(basis)
			A[i,j] = bf(basis[i], basis[j])
			B[i,j] = sip(basis[i], basis[j])
		end
	end
	bf_orthogonality_test = norm(A-diagm(diag(A)))/norm(A)
	sip_orthogonality_test = norm(B-diagm(diag(B)))/norm(B)
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
    out = sqrt(out/sum(out))
end

function coefficientFormulaTest(α::Real, maxdeg::Integer)
	basis = SOP(α, maxdeg)
	w = real(-[basis[ZernikeSuite.positionRange(deg)[i+1]].coefficients[ZernikeSuite.positionRange(deg-2)[i]] for deg in 3:maxdeg for i in 1:deg-1])
	d = [deg for deg in 3:maxdeg for i in 1:deg-1]
	md = [abs(2*i-deg) for deg in 3:maxdeg for i in 1:deg-1]
	modelVal = (d-md).*(d+md)./(d-md+2*α)./(d+md+2*α)
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

# Test of a possibly useful identity
function tenthAugustTest(α::Real, m::Integer, n::Integer)
	obj = ZernikeSuite.lower(ZernikePoly(α+1, m, n))
	mbump(f::ZFun) = f - mzp(mzs(f))
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
	mbump(f::ZFun) = f - mzp(mzs(f))
	lhs1 = -m*(α+1)*mzp(mzs(obj)) - (α+2)*mbump(mzp(dzp(obj))) + m*mbump(obj + mzs(dzs(obj))) + mbump(mbump(dzs(dzp(obj))))
	rhs1 = -m*(α+1)*ZernikePoly((α+1)-1, m, n)
	res1 = rhs1 - lhs1
	relErr1 = wip(res1,res1)==wip(rhs1,rhs1)==0.0?0.0:sqrt(real(wip(res1,res1))/real(wip(rhs1,rhs1)))
	lhs2 = -n*(α+1)*mzp(mzs(obj)) - (α+2)*mbump(mzs(dzs(obj))) + n*mbump(obj + mzp(dzp(obj))) + mbump(mbump(dzs(dzp(obj))))
	rhs2 = -n*(α+1)*ZernikePoly((α+1)-1, m, n)
	res2 = rhs2 - lhs2
	relErr2 = wip(res2,res2)==wip(rhs2,rhs2)==0.0?0.0:sqrt(real(wip(res2,res2))/real(wip(rhs2,rhs2)))
	return relErr1, relErr2
end
