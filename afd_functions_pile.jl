function czt(A::Complex128,W::Complex128,k::Int,xn::Vector{Complex128})
    #m = length of xn
    #k = desired output length
    #L : nextpow2(m-1+k))
    #yn: xn[n]*A^(-n)*W^(n^2/2), running from n = 0:(m-1), padded to L
    #vn: the inner filter W^(-n^2/2), running from n = (-m+1):(k-1), padded to L
    #vk: the outer product W^(m^2), running from m = 0:(k-1), padded to L
    #t1, t2: temporary variables
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
        yn[i] = (p<m) ? xn[i]*t*A^p : zero(Complex128)
        vk[i] = (p<k) ? t : zero(Complex128)
        vn[i+m-1] = (p<k) ? one(Complex128)/t : zero(Complex128)
        l = m-i+1
        if l>0
            vn[l] = one(Complex128)/t
        end
    end
    gn = ifft(fft(vn).*fft(yn))
    return [gn[m+i-1]*vk[i] for i=1:k]
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
    for i=1:length(G)
        G[i] = f(Q[i])
    end
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
        Gn[i] = ((Gnm1[i] - Onm1*dictionaryElement(Q[i],Anm1))*invBlaschkeElement(Q[i],Anm1))
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
    return dataLen*log(RSS/dataLen)+2*numParams+2*numParams*((numParams+1)/(dataLen-numParams-1))
end

function afd(xn,Nmax,Nz)
    Nz = int(floor(sqrt(Nz)))
    (Z,Q) = unitDiskGrid(Nz,0.005)
    L = length(Z)
    Gnm1 = zeros(Complex128,L)
    Gn   = zeros(Complex128,L)
    O = zeros(Complex128,L)
    polesList = zeros(Complex128,Nmax);
    polesAmps = zeros(Complex128,Nmax);
    poleIndList = zeros(Int64,Nmax);
    g1!(Gnm1,z->zTransform(z,xn),Q)
    (poleNm1,poleAmpNm1,poleIndNm1) = innerProduct(Gnm1,Z)
    polesList[1] = poleNm1
    polesAmps[1] = poleAmpNm1
    poleIndList[1] = poleIndNm1;
    for i=2:Nmax
        gNext!(Gn,Gnm1,poleNm1,poleAmpNm1,Q)
        healG!(Gn,poleIndList,i-1)
        (poleNm1,poleAmpNm1,poleIndNm1) = innerProduct(Gn,Z)
        polesList[i]=poleNm1
        polesAmps[i]=poleAmpNm1
        poleIndList[i] = poleIndNm1;
        Gnm1 = Gn
    end
    return polesList, polesAmps
end

function afd_with_initial_poles(xn,pinit,Nmax,Nz)
    Lpinit = length(pinit)
    if Lpinit > Nmax
        error("Nmax must be greater than the number of initial poles")
    end
    Nz = int(floor(sqrt(Nz)))
    (Z,Q) = unitDiskGrid(Nz,0.1)
    Ldisk = length(Z)
    Z = vcat(Z,pinit)
    Q = vcat(Q,1./conj(pinit))
    L = Ldisk+Lpinit
    Gnm1 = zeros(Complex128,L)
    Gn   = zeros(Complex128,L)
    O = zeros(Complex128,L)
    polesList = zeros(Complex128,Nmax)
    polesAmps = zeros(Complex128,Nmax)
    poleIndList = zeros(Int64,Nmax)
    g1!(Gnm1,z->zTransform(z,xn),Q)
    poleAmpNm1 = innerProduct(Gnm1[Ldisk+1],pinit[1])
    polesList[1] = pinit[1]
    polesAmps[1] = poleAmpNm1
    poleIndList[1] = Ldisk+1
    for i=2:Lpinit
        gNext!(Gn,Gnm1,pinit[i-1],poleAmpNm1,Q)
        healG!(Gn,poleIndList,i-1)
        poleAmpNm1 = innerProduct(Gnm1[Ldisk+i],pinit[1])
        polesList[i]=pinit[i]
        polesAmps[i]=poleAmpNm1
        poleIndList[i] = Ldisk+i
        Gnm1 = Gn
    end
    poleNm1 = pinit[end]
    for i=(Lpinit+1):Nmax
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