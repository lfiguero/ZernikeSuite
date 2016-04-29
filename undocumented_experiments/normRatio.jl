include("ZernikeSuite.jl")

col = (α,n) -> [ZernikeSuite.ZernikePoly(α+1, i, n-i) for i=0:n];
ratio = f -> (fshifted = ZernikeSuite.lower(f); real(ZernikeSuite.wip(fshifted,fshifted)/ZernikeSuite.wip(f,f)));

function randomguy(α::Real, n::Integer)
    bf = col(α,n);
    v = randn(n+1)+im*randn(n+1);
    reduce(+, map(*, v, bf))
end

function randomratio(α::Real, n::Integer)
    randomOP = randomguy(α,n);
    ratio(randomOP)
end

maxratio = (α,n) -> maximum(map(ratio, col(α,n)));
myformula = (α,n) -> (n+1)/(α+1) + 1;
