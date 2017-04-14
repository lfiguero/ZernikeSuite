include("../ZernikeSuite.jl")
using ZernikeSuite

function sip(f::ZFun, g::ZFun)
    @assert f.α == g.α
    wip(dx(f),dx(g)) + wip(dy(f),dy(g)) + wip(proj(f,0), proj(g,0))
end

function bf(f::ZFun, g::ZFun)
    @assert f.α == g.α
    ZRaise = ZernikeSuite.raise
    dxxf = ZRaise(dx(dx(f)))
    dxyf = ZRaise(dx(dy(f)))
    dyxf = ZRaise(dy(dx(f)))
    dyyf = ZRaise(dy(dy(f)))
    dxxg = ZRaise(dx(dx(g)))
    dxyg = ZRaise(dx(dy(g)))
    dyxg = ZRaise(dy(dx(g)))
    dyyg = ZRaise(dy(dy(g)))
    wip(dxxf,dxxg) + wip(dxyf,dxyg) + wip(dyxf,dyxg) + wip(dyyf,dyyg) + wip(dθ(dx(f)), dθ(dx(g))) + wip(dθ(dy(f)), dθ(dy(g)))
end

function SturmLiouvilleTest(α::Real, maxdeg::Integer)
    @assert α > -1
    @assert maxdeg ≥ 0
    # We start with basis which is weighted L²-orthogonal
    basis = [ZernikePoly(α, i, k-i) for k in 0:maxdeg for i in 0:k]
    dim = length(basis)
    # We turn the basis into a weighted Sobolev-orthogonal one via a
    # Gram–Schmidt process
    for i = 1:dim
	for j = 1:i-1
	    basis[i] = basis[i] - sip(basis[i], basis[j]) * basis[j]
	end
	basis[i] = basis[i]/sqrt(sip(basis[i], basis[i]))
    end
    A = Array{Complex128}(dim, dim)
    B = Array{Complex128}(dim, dim)
    for i = 1:dim
	for j = 1:dim
	    A[i,j] = bf(basis[i], basis[j])
	    B[i,j] = sip(basis[i], basis[j])
	end
    end
    bf_orthogonality_test = norm(A-diagm(diag(A)))/norm(A)
    sip_orthogonality_test = norm(B-diagm(diag(B)))/norm(B)
    eigenvalues = diag(A) ./ diag(B)
    return bf_orthogonality_test, sip_orthogonality_test, basis, eigenvalues
end
