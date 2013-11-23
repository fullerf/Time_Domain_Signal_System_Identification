function generateU(ul::Vector{Complex128},xn::Vector{Complex128})
	L = length(ul) #Num of basis functions
	N = length(xn) #signal length
	p = 1 #calculate U matrices up to this num (1) for U0,U1
	K = int((N-1-p)/2)
    NPOW = 3

	ul_inv = [1/ul[i] for i=1:L]
	ul_invK = [Base.power_by_squaring(ul_inv[i],K) for i=1:L]

	#initialize memory
	g0 = zeros(Complex128,L)
	g0_K = zeros(Complex128,L)
	D0 = zeros(Complex128,L)
	U0 = zeros(Complex128,(L,L))
	U1 = zeros(Complex128,(L,L))
    ul_invk = ones(Complex128,L)
	for i=0:(K)
		for j=0:(L-1)
			g0[j+1] += ul_invk[j+1]*xn[i+1]
			g0_K[j+1] += ul_invk[j+1]*xn[K+1+i+1]
            D0[j+1] += (i + 1)*xn[i+1]*ul_invk[j+1] + (K - i)*xn[i + K + 1+1]*ul_invk[j+1]*ul_inv[j+1]*ul_invK[j+1]
            #D0[j+1] = D0partial[j+1] + D0partial2[j+1]
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
	return (U0,U1,g0,g0_K)
end