function czt(A::Complex128,W::Complex128,k::Int,xn::Vector{Complex128})
    #m = length of xn
    #k = desired output length
    #L : nextpow2(m-1+k))
    #yn: xn[n]*A^(-n)*W^(n^2/2), running from n = 0:(m-1), padded to L
    #vn: the inner filter W^(-n^2/2), running from n = (-m+1):(k-1), padded to L
    #vk: the outer product W^(m^2), running from m = 0:(k-1), padded to L
    #t : temporary variable
    m = length(xn)
    L = nextpow2(m+k-1)
    M = maximum([(m-1) (k-1)])
    yn = zeros(Complex128,L)
    vn = zeros(Complex128,L)
    vk = zeros(Complex128,L)
    W2 = sqrt(W)
    for i=1:(M+1)
        p = i-1
        t = W2^(p*p)
        yn[i] = (p<m) ? xn[i]*t*A^(-p) : zero(Complex128)
        vk[i] = (p<k) ? t : zero(Complex128)
        l1 = (i+m-1) 
        if (l1<L)
            vn[i+m-1] = (p<k) ? one(Complex128)/t : zero(Complex128)
        end
        l2 = m-i+1
        if l2>0
            vn[l2] = one(Complex128)/t
        end
    end
    gn = ifft(fft(vn).*fft(yn))
    return [(gn[m+i-1]*vk[i])::Complex128 for i=1:k]
end

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

function zTransformHorner(z::Complex128,xn::Vector{Complex128})
    #uses Horner's method
    x = one(Complex128)/z
    r = xn[end]
    for k=2:length(xn)
        r = r*x+xn[end-k+1]
    end
    return r*x
end

function zTransform(z::Complex128,xn::Vector{Complex128})
    x = one(Complex128)/z
    r = zero(Complex128)
    for k=1:length(xn)
        Z = x^k
        r += xn[k]*Z
    end
    return r
end

function alphaj(a::Complex128)
    return sqrt(1 - abs(a)^2)
end

function dictionaryElement(z::Complex128,a::Complex128)
    return alphaj(a)/(z-a)
end

function blaschkeElement(z::Complex128,a::Complex128)
    return (1-conj(a)*z)/(z-a)
end

function invBlaschkeElement(z::Complex128,a::Complex128)
    return (z-a)/(1-conj(a)*z)
end

function blaschkeProd(z::Complex128,aj::Vector{Complex128})
    r = one(z)
    for i=1:length(aj)
        r *= blaschkeElement(z,aj[i])
    end
    return r
end

function blaschkeProdExceptQ(z::Complex128,aj::Vector{Complex128},q::Int)
    r = one(z)
    for i=1:length(aj)
        if i!=q
            r *= blaschkeElement(z,aj[i])
        end
    end
    return r
end

function innerProduct(g::Complex128,a::Complex128)
    #a must be INTERIOR to the unit disk
    return alphaj(a)*(1/conj(a))*g
end

function innerProduct(G::Vector{Complex128},Z::Vector{Complex128})
    #Z is a collection poles INTERIOR to the unit disk
    L = length(Z)
    if length(G)!=L
        error("Length of G must match that of Z")
    end
    OverlapOld = zero(Float64)
    OverlapNew = zero(Float64)
    I = zero(Complex128)
    maxind = zero(Int64)
    f(On,Oo,Oin,Oio,Ivn,Ivo) = On > Oo ? (On,Oin,Ivn) : (Oo,Oio,Ivo)
    for i=1:L
        t = innerProduct(G[i],Z[i])
        OverlapNew = abs2(t)
        (OverlapOld,maxind,I) = f(OverlapNew,OverlapOld,i,maxind,t,I)
    end
    return (Z[maxind],I,maxind)
end

function OverlapVector(G::Vector{Complex128},Z::Vector{Complex128})
    #Z is a collection poles INTERIOR to the unit disk
    L = length(Z)
    if length(G)!=L
        error("Length of G must match that of Z")
    end
    O = zeros(Float64,L)
    for i=1:L
        O[i] = abs2(innerProduct(G[i],Z[i]))
    end
    return O
end

function g1!(G::Vector{Complex128},f::Function,Q::Vector{Complex128})
    #Q is a collection of poles EXTERIOR to the unit disk
    L = length(G)
    if length(Q)!=L
        error("input G vector must have same length as input Q vector")
    end
    for i=1:L
        G[i] = f(Q[i])
    end
end

function g1czt!(G::Vector{Complex128},xn::Vector{Complex128},Q::Vector{Complex128})
    zxn = vcat(zero(Complex128),xn)
    #Q is a collection of poles EXTERIOR to the unit disk
    lcol(i,l) = ((i-1)*l+i):((i-1)*l+((l+i-1)-2*(i-1)))
    rcol(i,l) = (((l-i)*l+((l+i-1)-2*(i-1)):-1:((l-i)*l+i)))
    brow(i,l) = (i*l+((l+i-1)-2*(i-1))):l:((l-i-1)*l+((l+i-1)-2*(i-1)))
    trow(i,l) = ((l-i-1)*l+i):-l:(i*l+i)
    L = length(G)
    l = int(sqrt(L)) #L should be a perfect square, and l is a power of 2 by construction
    if length(Q)!=L
        error("input G vector must have same length as input Q vector")
    end
    ulCorner(i,l) = ((i-1)*l+i) #upper left corner
    for i=1:int(l/2)
        side_len = l-2*(i-1)
        A = Q[ulCorner(i,l)]
        K = (4*(l-2*(i-1))-4)
        W = exp(-2*pi*im/K)
        temp = czt(A,W,K,zxn)
        G[lcol(i,l)]=temp[1:side_len]
        G[brow(i,l)]=temp[(side_len+1):(2*(side_len-1))]
        G[rcol(i,l)]=temp[(2*side_len-1):(3*side_len-2)]
        G[trow(i,l)]=temp[(3*side_len-1):(4*(side_len-1))]
    end
end

function dictBlaschke(z::Complex128,a::Complex128)
    return alphaj(a)/(1-conj(a)*z)
end


function gNext!(Gn::Vector{Complex128},Gnm1::Vector{Complex128},Anm1::Complex128,Onm1::Complex128,Q::Vector{Complex128})
    #this function modifies Gn, so returns nothing.
    #Q is a collection of poles EXTERIOR to the unit disk
    L = length(Q)
    if length(Gnm1)!=L
        error("input gnm1 vector must have same length as input z vector")
    end
    if length(Gn)!=L
        error("input gn vector must have same length as input z vector")
    end
    for i=1:L
        Gn[i] = (Gnm1[i]*invBlaschkeElement(Q[i],Anm1) - Onm1*dictBlaschke(Q[i],Anm1))
        #Gn[i] = ((Gnm1[i] - Onm1*dictionaryElement(Q[i],Anm1))*invBlaschkeElement(Q[i],Anm1))
    end
end

function impulseTakenakaMalmquist(n,aj::Vector{Complex128})
    if n<=zero(n)
        return zero(Complex128)
    end
    temp1 = zero(Complex128)
    temp2 = zero(Complex128)
    for i=1:(length(aj)-1)
        temp1 += aj[i]^(n-1)*dictionaryElement(aj[i],aj[end])*blaschkeProdExceptQ(aj[i],aj[1:(length(aj)-1)],i)*(1-abs(aj[i])^2)
    end
    temp2 += aj[end]^(n-1)*alphaj(aj[end])*blaschkeProd(aj[end],aj[1:(length(aj)-1)])
    return temp1 + temp2
end

function unitDiskGrid(n::Int,delta::Float64)
    N = nextpow2(n)
    M = N>>1
    v = linspace(delta,one(Float64)-delta,M)
    Zjk = zeros(Complex128,N,N)
    for j=1:(M)
             #bottom, running counter-clockwise
            Zjk[(M+j),((M-(j-1)):(M+1+(j-1)))] = v[j]*exp(im*linspace(5*pi/4,7*pi/4,2*j))
            #top, running clockwise
            Zjk[(M-j+1),((M-(j-1)):(M+1+(j-1)))] = v[j]*exp(im*linspace(3*pi/4,1*pi/4,2*j))
    end
    for j=1:(M)
            #right, running clockwise
            Zjk[((M-(j-1)):(M+1+(j-1))),M+j] = v[j]*exp(im*linspace(pi/4,-pi/4,2*j)) #right
            #left, running counter-clockwise
            Zjk[((M-(j-1)):(M+1+(j-1))),M-j+1] = v[j]*exp(im*linspace(3*pi/4,5*pi/4,2*j)) #left
    end
    Qjk = one(Complex128)./conj(Zjk)
    return (reshape(Zjk,N^2), reshape(Qjk,N^2))
end

function healG!(G::Vector{Complex128},x::Vector{Int},L::Int)
    for k=1:L       
        G[x[k]] = zero(Complex128)
    end
end

function reconstructTimeDomain(poles::Vector{Complex128},amps::Vector{Complex128},t)
    if length(poles)!=length(amps)
        error("Input pole vector must have same length as amplitude vector")
    end
    r = zero(Complex128)
    for k=1:length(amps)
        p = k > 1 ? poles[1:k] : [poles[1]]
        r += amps[k]*impulseTakenakaMalmquist(t,p)
    end
    return r
end

function residualSumOfSquares(model::Vector{Complex128},data::Vector{Complex128})
    if length(model)!=length(data)
        error("Input model vector must have same length as the data vector!")
    end
    r = zero(Float64)
    for k=1:length(data)
        r += abs2(data[k]-model[k])
    end
    return r
end

function akaikeInformationCriterion(RSS,numParams,dataLen)
    #This gives the sample size corrected Akaike Information Criteria
    #commonly known as AICc in the literature
    K = numParams+1
    return dataLen*log(RSS/dataLen)+2*K+2*K*((K+1)/(dataLen-K-1))
end

function afd(xn::Vector{Complex128},Nmax::Int,Nz::Int)
    Nz = int(floor(sqrt(Nz)))
    (Z,Q) = unitDiskGrid(Nz,0.005)
    L = length(Z)
    rNz = int(sqrt(L))
    Gnm1 = zeros(Complex128,L)
    Gn   = zeros(Complex128,L)
    polesList = zeros(Complex128,Nmax)
    polesAmps = zeros(Complex128,Nmax)
    poleIndList = zeros(Int64,Nmax)
    g1czt!(Gnm1,xn,Q)
    (poleNm1,poleAmpNm1,poleIndNm1) = innerProduct(Gnm1,Z)
    polesList[1] = poleNm1
    polesAmps[1] = poleAmpNm1
    poleIndList[1] = poleIndNm1
    for i=2:Nmax
        gNext!(Gn,Gnm1,poleNm1,poleAmpNm1,Q)
        healG!(Gn,poleIndList,i-1)
        (poleNm1,poleAmpNm1,poleIndNm1) = innerProduct(Gn,Z)
        polesList[i]=poleNm1
        polesAmps[i]=poleAmpNm1
        poleIndList[i] = poleIndNm1
        Gnm1 = Gn
    end
    return polesList, polesAmps
end

function xafd(xn::Vector{Complex128},Nmax::Int,Nz::Int,xN::Int,q::Float64)
    Nz = int(floor(sqrt(Nz)))
    (Z,Q) = unitDiskGrid(Nz,0.005)
    #extrapolate the time domain using xczt to xN points
    xhat = ifft(xczt(xn,one(Complex128),exp(-im*2*pi/xN),xN,q))
    L = length(Z)
    rNz = int(sqrt(L))
    Gnm1 = zeros(Complex128,L)
    Gn   = zeros(Complex128,L)
    polesList = zeros(Complex128,Nmax);
    polesAmps = zeros(Complex128,Nmax);
    poleIndList = zeros(Int64,Nmax);
    g1czt!(Gnm1,xhat,Q)
    (poleNm1,poleAmpNm1,poleIndNm1) = innerProduct(Gnm1,Z)
    polesList[1] = poleNm1
    polesAmps[1] = poleAmpNm1
    poleIndList[1] = poleIndNm1
    for i=2:Nmax
        gNext!(Gn,Gnm1,poleNm1,poleAmpNm1,Q)
        healG!(Gn,poleIndList,i-1)
        (poleNm1,poleAmpNm1,poleIndNm1) = innerProduct(Gn,Z)
        polesList[i]=poleNm1
        polesAmps[i]=poleAmpNm1
        poleIndList[i] = poleIndNm1
        Gnm1 = Gn
    end
    return polesList, polesAmps
end

function optimalAfdTruncation(poles::Vector{Complex128},amps::Vector{Complex128},xn::Vector{Complex128})
    L = length(poles)
    Nt = length(xn)
    if length(amps)!=L
        error("# of poles should be the same as number of amplitudes")
    end
    infVec = zeros(Float64,L)
    for k=1:L
        model = [reconstructTimeDomain(poles[1:k],amps[1:k],i)::Complex128 for i=linspace(1,Nt,Nt)]
        rssval = residualSumOfSquares(model,xn)
        # one complex value = 2 parameters, so each pole and amp is 4 parameters, each complex valued input
        # is worth 2 data points of input in this thinking, thus:
        infVec[k] = akaikeInformationCriterion(rssval,4*k,2*length(xn))
    end
    infMinInd = indmin(infVec)
    aicVec = exp((-infVec+infVec[infMinInd])/2)
    return indmax(aicVec)
end