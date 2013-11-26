function generateU(ul::Vector{Complex128},xn::Vector{Complex128})
	L = length(ul) #Num of basis functions
	N = length(xn) #signal length
	p = 1 #calculate U matrices up to this num (1) for U0,U1
	K = ifloor((N-1-p)/2)
    NPOW = 3

	ul_inv = [1/ul[i] for i=1:L]
	ul_invK = [Base.power_by_squaring(ul_inv[i],K) for i=1:L]

	#initialize memory
	g0 = zeros(Complex128,L)
	g0_K = zeros(Complex128,L)
    g0_NK = zeros(Complex128,L)
	D0 = zeros(Complex128,L)
	U0 = zeros(Complex128,(L,L))
	U1 = zeros(Complex128,(L,L))
    ul_invk = ones(Complex128,L)
	for i=0:(K)
		for j=0:(L-1)
			g0[j+1] += ul_invk[j+1]*xn[i+1]
			g0_K[j+1] += ul_invk[j+1]*xn[K+1+i+1]
            g0_NK[j+1] += ul_invk[j+1]*xn[N-K+i]
            D0[j+1] += (i + 1)*xn[i+1]*ul_invk[j+1] + (K - i)*xn[i + K + 1+1]*ul_invk[j+1]*ul_inv[j+1]*ul_invK[j+1]
            ul_invk[j+1] = ((i%NPOW)==(NPOW-1)) ? Base.power_by_squaring(ul_inv[j+1],i+1) : ul_invk[j+1]*ul_inv[j+1]
		end
	end

	for i=1:L
		for j=1:i
			U0[i,j] = (1/(ul[i]-ul[j]))*(ul[i]*g0[j] - ul[j]*g0[i] 
					+ ul_invK[j]*g0_K[i] - ul_invK[i]*g0_K[j])
			U0[j,i] = U0[i,j]
			U1[i,j] = (1/(ul[i]-ul[j]))*(ul[i]*ul[j]*g0[j] - ul[j]*ul[i]*g0[i] 
					+ ul_invK[j]*ul[i]*g0_K[i] - ul_invK[i]*ul[j]*g0_K[j])
			U1[j,i] = U1[i,j]
		end
		U0[i,i] = D0[i]
		U1[i,i] = D0[i]*ul[i]-ul[i]*g0[i]+ul_invK[i]*g0_K[i]
	end
	return (U0,U1,g0,g0_K, g0_NK)
end

function xft(xn::Vector{Complex128},ul::Vector{Complex128},Z::Vector{Complex128},q::Float64)
    
    N = length(xn)
    if length(Z)<N
        error("Length of Z must be at least the length of N")
    end
    M = ifloor((N-2)/2)+1
    (U0,U1,g0,g0_K,g0_NK) = generateU(ul,xn)
    genR(z,U0,U1) = U0-U1/z
    U0dagU0 = ctranspose(U0)*U0
    U1dagU1 = ctranspose(U1)*U1
    U1dagU0 = ctranspose(U1)*U0
    U0dagU1 = ctranspose(U0)*U1
    genRdagG0(R,g0) = ctranspose(R)*g0
    genRdagR(z,U0dagU0,U1dagU1,U1dagU0,U0dagU1) = U0dagU0 + U1dagU1 -z*U1dagU0 - U0dagU1/z
    Spec = zeros(Complex128,length(Z))
    for k=1:length(Z)
        R  = genR(Z[k],U0,U1)
        RC2 = genRdagG0(R,g0_NK*Z[k]^(-N+M))
        RR = genRdagR(Z[k],U0dagU0,U1dagU1,U1dagU0,U0dagU1)
        Spec[k] = (transpose(g0_K*Z[k]^(-M))*(\((RR+q^2),RC2)))[1]
    end
    return Spec+fft(vcat(xn,zeros(length(Z)-N)))
end

function rrt(xn::Vector{Complex128},ul::Vector{Complex128},Z::Vector{Complex128},q::Float64)
    N = length(xn)
    M = ifloor((N-2)/2)+1
    (U0,U1,g0,g0_K,g0_NK) = generateU(ul,xn)
    genR(z,U0,U1) = U0-U1/z
    U0dagU0 = ctranspose(U0)*U0
    U1dagU1 = ctranspose(U1)*U1
    U1dagU0 = ctranspose(U1)*U0
    U0dagU1 = ctranspose(U0)*U1
    genRdagG0(R,g0) = ctranspose(R)*g0
    genRdagR(z,U0dagU0,U1dagU1,U1dagU0,U0dagU1) = U0dagU0 + U1dagU1 -z*U1dagU0 - U0dagU1/z
    Spec = zeros(Complex128,length(Z))
    for k=1:length(Z)
        R  = genR(Z[k],U0,U1)
        RC = genRdagG0(R,g0)
        RR = genRdagR(Z[k],U0dagU0,U1dagU1,U1dagU0,U0dagU1)
        Spec[k] = (transpose(g0)*(\((RR+q^2),RC)))[1]
    end
    return Spec
end

function xczt(xn::Vector{Complex128},A::Complex128,W::Complex128,K::Int,q::Float64)
    N = length(xn)
    Q = N*sum(abs(xn))
    q = Q*q
    Nul = minimum((N,100))
    ul = [exp(im*2*pi/Nul)^(i-1) for i=1:Nul] #fourier basis for expanding Krylov operators U0,U1
    Z = [A*W^(-(i-1)) for i=1:K] #this is the spiral the CZT is to be evaluated on
    M = ifloor((N-2)/2)+1
    (U0,U1,g0,g0_K,g0_NK) = generateU(ul,xn)
    genR(z,U0,U1) = U0-U1/z
    U0dagU0 = ctranspose(U0)*U0
    U1dagU1 = ctranspose(U1)*U1
    U1dagU0 = ctranspose(U1)*U0
    U0dagU1 = ctranspose(U0)*U1
    genRdagG0(R,g0) = ctranspose(R)*g0
    genRdagR(z,U0dagU0,U1dagU1,U1dagU0,U0dagU1) = U0dagU0 + U1dagU1 -z*U1dagU0 - U0dagU1/z
    Spec = zeros(Complex128,length(Z))
    for k=1:length(Z)
        R  = genR(Z[k],U0,U1)
        RC2 = genRdagG0(R,g0_NK*Z[k]^(-N+M))
        RR = genRdagR(Z[k],U0dagU0,U1dagU1,U1dagU0,U0dagU1)
        Spec[k] = (transpose(g0_K*Z[k]^(-M))*(\((RR+q^2),RC2)))[1]
    end
    return Spec+czt(A,W,K,xn)
end

function rczt(xn::Vector{Complex128},A::Complex128,W::Complex128,K::Int,q::Float64)
    N = length(xn)
    Q = N*sum(abs(xn))
    q = Q*q
    Nul = minimum((N,100))
    ul = [exp(im*2*pi/Nul)^(i-1) for i=1:Nul] #fourier basis for expanding Krylov operators U0,U1
    Z = [A*W^(-(i-1)) for i=1:K] #this is the spiral the CZT is to be evaluated on
    M = ifloor((N-2)/2)+1
    (U0,U1,g0,g0_K,g0_NK) = generateU(ul,xn)
    genR(z,U0,U1) = U0-U1/z
    U0dagU0 = ctranspose(U0)*U0
    U1dagU1 = ctranspose(U1)*U1
    U1dagU0 = ctranspose(U1)*U0
    U0dagU1 = ctranspose(U0)*U1
    genRdagG0(R,g0) = ctranspose(R)*g0
    genRdagR(z,U0dagU0,U1dagU1,U1dagU0,U0dagU1) = U0dagU0 + U1dagU1 -z*U1dagU0 - U0dagU1/z
    Spec = zeros(Complex128,length(Z))
    for k=1:length(Z)
        R  = genR(Z[k],U0,U1)
        RC = genRdagG0(R,g0)
        RR = genRdagR(Z[k],U0dagU0,U1dagU1,U1dagU0,U0dagU1)
        Spec[k] = (transpose(g0)*(\((RR+q^2),RC)))[1]
    end
    return Spec
end