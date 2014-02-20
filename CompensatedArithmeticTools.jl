function TwoSum(a::Float64,b::Float64)
	x = a+b
	z = x-a
	y = (a-(x-z))+(b-z)
	return (x,y)
end

function Split(a::Float64)
	p = 52
	s = 26
	factor = 67108865
	c = factor*a
	x = c - (c-a)
	y = a - x
    return(x,y)
end

function Split(a::Float32)
	p = 23
	s = 12
	factor = 4097
	c = factor*a
	x = c - (c-a)
	y = a - x
    return (x,y)
end

function TwoProduct(a::Float64,b::Float64)
	x = a*b
	[a1,a2] = Split(a)
	[b1,b2] = Split(b)
	y = (a2*b2−(((x−a1*b1)−a2*b1)−a1*b2))
end

function TwoSumComplex(x::Complex128,y::Complex128)
	a = real(x)
	b = imag(x)
	c = real(y)
	d = imag(y)
	(s1,e1) = TwoSum(a,c)
	(s2,e2) = TwoSum(b,d)
	s = s1 + im*s2
	e = e1 + im*e2
	return (s,e)
end

function TwoProductComplex(x::Complex128,y::Complex128)
	a = real(x)
	b = imag(x)
	c = real(y)
	d = imag(y)
	(z1,h1) = TwoProduct(a,c)
	(z2,h2) = TwoProduct(b,d)
	(z3,h3) = TwoProduct(a,d)
	(z4,h4) = TwoProduct(b,c)
	(z5,h5) = TwoSum(z1,-z2)
	(z6,h6) = TwoSum(z3,z4)
	p = z5+im*z6
	e = h1+im*h3
	f = -h2+im*h4
	g = h5+im*h6
	return (p,e,f,g)
end