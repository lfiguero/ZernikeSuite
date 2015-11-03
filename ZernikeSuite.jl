module ZernikeSuite

export ZFun, dzs, dzp, dx, dy, proj, wip, w_sobolev_sq_sn, w_nc_sobolev_sq_sn, all_w_sobolev_sq_sn, all_w_nc_sobolev_sq_sn, ZernikePoly

function isPolySpaceDim(l::Integer)
    # Given l returns a Boolean and an integer; the Boolean is set to true if l
    # is of the form (n+1)*(n+2)/2 for some integer n and in that case the
    # integer returned is that n; otherwise the integer returned is the ceiling
    # of the non-integer solution of (n_1)*(n+2)/2 = l.
    n = (-3 + sqrt(1+8*l))/2 # This will be a float
    cn = convert(Int64, ceil(n))
    2*l == (cn+1)*(cn+2), cn
end

# ZFun type and constructors
type ZFun
    α::Real
    degree::Integer
    coefficients::Vector{Complex128}
    function ZFun(α, degree, coefficients)
	@assert α > -1
	2*length(coefficients)==(degree+1)*(degree+2) ? new(α, degree, coefficients) : error("length(coefficients) should be div((degree+1)*(degree+2),2)")
    end
end

function ZFun{T<:Number}(α::Real, coefficients::Vector{T})
    b, n = isPolySpaceDim(length(coefficients))
    if b
	newcoefficients = coefficients
    else
	cl = div((n+1)*(n+2),2)
	newcoefficients = zeros(Complex128, cl)
	newcoefficients[1:length(coefficients)] = coefficients
    end
    ZFun(α, n, newcoefficients)
end

# Unary operations
-(f::ZFun) = ZFun(f.α, f.degree, -f.coefficients)

# Binary operations
for op = (:+, :-)
    @eval begin
	function ($op)(f::ZFun, g::ZFun)
	    @assert f.α == g.α
	    fl = length(f.coefficients)
	    gl = length(g.coefficients)
	    retl = max(fl, gl)
	    retd = max(f.degree, g.degree)
	    retcoefficients = zeros(Complex128, retl)
	    retcoefficients[1:fl] = f.coefficients;
	    retcoefficients[1:gl] = ($op)(retcoefficients[1:gl], g.coefficients);
	    ZFun(f.α, retd, retcoefficients)
	end
    end
end

# Operations with scalars
for op = (:+, :-, :*, :/)
    @eval begin
	function ($op)(f::ZFun, a::Number)
	    ZFun(f.α, f.degree, ($op)(f.coefficients, a))
	end
    end
end
for op = (:+, :*)
    @eval begin
	($op)(a::Number, f::ZFun) = ($op)(f, a)
    end
end
-(a::Number, f::ZFun) = a + (-f)

# List of position ranges of coefficients up to a certain degree
positionRanges(maxdeg::Integer) = [div(k*(k+1),2)+1+(0:k) for k in 0:maxdeg]

# Raise parameter by one
function raise(f::ZFun)
    pr = positionRanges(f.degree)
    reta = f.α + 1.0
    retd = f.degree
    retc = zeros(Complex128, size(f.coefficients))
    for k = 0:f.degree
	krange = pr[k+1]
	i = [0:k]
	Adiag = (i+f.α+1) .* (k-i+f.α+1) ./ (k+f.α+1) / (f.α+1)
	v = f.coefficients[krange] .* Adiag
	if k < f.degree-1
	    kplustworange = pr[k+3]
	    i = [1:(k+1)]
	    Bdiag = i .* (k+2-i) ./ (k+f.α+3) / (f.α+1)
	    v = v - f.coefficients[kplustworange][2:k+2] .* Bdiag
	end
	retc[krange] = v
    end
    ZFun(reta, retd, retc)
end

# Lower parameter by one
function lower(f::ZFun)
    @assert f.α > 0
    pr = positionRanges(f.degree)
    reta = f.α - 1
    retd = f.degree
    retc = zeros(Complex128, size(f.coefficients))
    for k = f.degree:-1:0
	krange = pr[k+1]
	i = [0:k]
	Adiag = (i+f.α) .* (k-i+f.α) ./ (k+f.α) / f.α
	rhs = f.coefficients[krange]
	if k < f.degree-1
	    kplustworange = pr[k+3]
	    i = [1:(k+1)]
	    Bdiag = i .* (k+2-i) ./ (k+f.α+2) / f.α
	    rhs = rhs + retc[kplustworange][2:k+2] .* Bdiag
	end
	retc[krange] = rhs ./ Adiag
    end
    ZFun(reta, retd, retc)
end

# Differentiation with parameter shift
function dzsShift(f::ZFun)
    # Differential operator 0.5*d/dx + 0.5*im * d/dy
    pr = positionRanges(f.degree)
    reta = f.α + 1.0
    if f.degree == 0
	return ZFun(reta, 0, [0.0])
    end
    retd = f.degree - 1
    retc = zeros(Complex128, div((retd+1)*(retd+2),2))
    for k = 1:f.degree
	krange = pr[k+1]
	kminusonerange = pr[k]
	i = [0:k-1]
	v = f.coefficients[krange][1:k] .* (i+f.α+1) .* (k-i) / (f.α+1)
	retc[kminusonerange] = v
    end
    ZFun(reta, retd, retc)
end

function dzpShift(f::ZFun)
    # Differential operator 0.5*d/dx - 0.5*im * d/dy
    pr = positionRanges(f.degree)
    reta = f.α + 1.0
    if f.degree == 0
	return ZFun(reta, 0, [0.0])
    end
    retd = f.degree - 1
    retc = zeros(Complex128, div((retd+1)*(retd+2),2))
    for k = 1:f.degree
	krange = pr[k+1]
	kminusonerange = pr[k]
	i = [1:k]
	v = f.coefficients[krange][2:k+1] .* (k-i+f.α+1) .* i / (f.α+1)
	retc[kminusonerange] = v
    end
    ZFun(reta, retd, retc)
end

# Differentiation
dzs(f::ZFun) = lower(dzsShift(f))
dzp(f::ZFun) = lower(dzpShift(f))
dx(f::ZFun) = (dzs(f) + dzp(f))
dy(f::ZFun) = (dzs(f) - dzp(f)) / (im)
for op = (:dzs, :dzp, :dx, :dy)
    @eval begin
	function ($op)(f::ZFun, k::Integer)
	    @assert k >= 0
	    out = f
	    for i = 1:k
		out = ($op)(out)
	    end
	    out
	end
    end
end

# Weighted L^2 projection of a ZFun to the spaces of polynomials of degree
# lower than or equal to a fixed degree
function proj(f::ZFun, d::Integer)
    ed = min(f.degree, d)
    ZFun(f.α, ed, f.coefficients[1:div((ed+1)*(ed+2),2)])
end

# Squared weighted L^2 norms of the Zernike basis functions
function h(α::Real, maxdeg::Integer)
    # In order to avoid under/overflows we will compute the norms recursively
    # instead of using gamma function evaluations
    @assert α > -1 && maxdeg >= 0
    pr = positionRanges(maxdeg)
    out = zeros(Float64, div((maxdeg+1)*(maxdeg+2),2))
    out[1] = pi
    for k = 0:maxdeg-1
	for i = 0:div(k,2)
	    out[pr[k+2][i+1]] = (k-i+1)/(k-i+α+1) * out[pr[k+1][i+1]]
	end
	out[pr[k+2]] = out[pr[k+2]] + reverse(out[pr[k+2]])
	if mod(k,2)==1
	    out[pr[k+2][div(k+3,2)]] = (div(k-1,2)+1)/(div(k-1,2)+α+1) * out[pr[k+1][div(k-1,2)+1]]
	end
	out[pr[k+1]] = out[pr[k+1]]/(k+α+1)
    end
    out[pr[maxdeg+1]] = out[pr[maxdeg+1]]/(maxdeg+α+1)
    out
end

# Weighted L^2 inner product
function wip(f::ZFun, g::ZFun)
    @assert f.α == g.α
    md = min(f.degree, g.degree)
    l = div((md+1)*(md+2),2)
    pr = positionRanges(md)
    wgth = h(f.α, md)
    (g.coefficients[1:l]'*(wgth.*f.coefficients[1:l]))[1]
end

# Higher gradient collections
function hgcoll(f::ZFun, m::Integer)
    # Collection of canonical derivatives up to degree m
    @assert m >= 0
    dcoll = Array(ZFun, div((m+1)*(m+2),2))
    dcoll[1] = f
    pr = positionRanges(m)
    for k = 1:m
	for ind = 1:k
	    dcoll[pr[k+1][ind]] = dx(dcoll[pr[k][ind]])
	end
	dcoll[pr[k+1][k+1]] = dy(dcoll[pr[k][k]])
    end
    dcoll
end

function hgcollnc(f::ZFun, m::Integer)
    # Collection of non canonical derivatives up to degree m
    @assert m >= 0
    dcoll = Array(ZFun, div((m+1)*(m+2),2))
    dcoll[1] = f
    pr = positionRanges(m)
    for k = 1:m
	for ind = 1:k
	    dcoll[pr[k+1][ind]] = dzs(dcoll[pr[k][ind]])
	end
	dcoll[pr[k+1][k+1]] = dzp(dcoll[pr[k][k]])
    end
    dcoll
end

# Weighted W^{m,2} squared seminorms
function w_sobolev_sq_sn(f::ZFun, m::Integer)
    # Conventional weighted Sobolev squared seminorm
    derivatives_of_f = hgcoll(f, m)
    out = 0.0
    for i = 0:m
	out = out + real(wip(derivatives_of_f[end-i], derivatives_of_f[end-i]))
    end
    out
end

function w_nc_sobolev_sq_sn(f::ZFun, m::Integer)
    # Non conventional weighted Sobolev squared seminorm; it is equivalent
    # (with constant depending on m) to the conventional weighted Sobolev
    # squared seminorm
    derivatives_of_f = hgcollnc(f, m)
    out = 0.0
    for i = 0:m
	out = out + real(wip(derivatives_of_f[end-i], derivatives_of_f[end-i]))
    end
    out
end

# Array with all weighted W^{k,2} squared seminorms from k to m
function all_w_sobolev_sq_sn(f::ZFun, m::Integer)
    # Conventional weighted Sobolev squared seminorms
    derivatives_of_f = hgcoll(f, m)
    out = zeros(m+1)
    c = 0
    for k = 0:m
	for i = 0:k
	    out[k+1] = out[k+1] + real(wip(derivatives_of_f[c+i+1], derivatives_of_f[c+i+1]))
	end
	c = c+k+1
    end
    out
end

function all_w_nc_sobolev_sq_sn(f::ZFun, m::Integer)
    # Non conventional weighted Sobolev squared seminorms; each is equivalent
    # (with constant depending on seminorm order) to its corresponding
    # conventional weighted Sobolev squared seminorm
    derivatives_of_f = hgcollnc(f, m)
    out = zeros(m+1)
    c = 0
    for k = 0:m
	for i = 0:k
	    out[k+1] = out[k+1] + real(wip(derivatives_of_f[c+i+1], derivatives_of_f[c+i+1]))
	end
	c = c+k+1
    end
    out
end

# Zernike polynomials
function ZernikePoly(α::Real, m::Integer, n::Integer)
    retd = m+n;
    retc = zeros(Complex128, div((retd+1)*(retd+2),2))
    retc[div(retd*(retd+1),2)+1+m] = 1.0
    ZFun(α, retd, retc)
end

end
