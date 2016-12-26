include("../ZernikeSuite.jl")
using ZernikeSuite

function sip(f::ZFun, g::ZFun)
    @assert f.α == g.α
    wip(dx(f),dx(g)) + wip(dy(f),dy(g)) + wip(proj(f,0), proj(g,0))
end

function bf(f::ZFun, g::ZFun)
    @assert f.α == g.α
    wip(ZernikeSuite.raise(dx(dx(f))), ZernikeSuite.raise(dx(dx(g)))) + wip(ZernikeSuite.raise(dx(dy(f))), ZernikeSuite.raise(dx(dy(g)))) + wip(ZernikeSuite.raise(dy(dx(f))), ZernikeSuite.raise(dy(dx(g)))) + wip(ZernikeSuite.raise(dy(dy(f))), ZernikeSuite.raise(dy(dy(g)))) + wip(dθ(dx(f)), dθ(dx(g))) + wip(dθ(dy(f)), dθ(dy(g)))
end

function SturmLiouvilleTest(α::Real, maxdeg::Integer)
    @assert α > -1
    @assert maxdeg ≥ 0
    protobasis = [ZernikePoly(α, i, k-i) for k in 0:maxdeg for i in 0:k]
    dim = length(protobasis)
    # Gram–Schmidt process
    basis = protobasis
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
    bf_orthogonality_test = norm(A-diagm(diag(A)))
    sip_orthogonality_test = norm(B-diagm(diag(B)))
    eigenvalues = diag(A) ./ diag(B)
    return bf_orthogonality_test, sip_orthogonality_test, eigenvalues
end
