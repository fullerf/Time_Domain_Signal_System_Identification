function en(an,anm1,Ln)
	return sqrt((1-abs2(an))/((one(an)-abs2(anm1))*(one(an)-abs2(Ln))))
end

function ln(z::Vector{Complex128},an::Complex128,anm1::Complex128,phinm1::Vector{Complex128},phinm1star::Vector{Complex128})
    f1 = [((z[i]-anm1)/(one(Complex128)-conj(an)*z[i]))*conj(phinm1[i]) for i=1:length(z)]
    f2 = [((one(Complex128)-conj(anm1)*z[i])/(one(Complex128)-conj(an)*z[i]))*conj(phinm1star[i]) for i=1:length(z)]
    num = sum(f1)
    den = sum(f2)
    return -num/den
end

function phiNext(z::Vector{Complex128},an::Complex128,anm1::Complex128,phinm1::Vector{Complex128},phinm1star::Vector{Complex128})
    phin = zeros(Complex128,length(phinm1))
    phinstar = zeros(Complex128,length(phinm1))
    L = ln(z,an,anm1,phinm1,phinm1star)
    c = en(an,anm1,L)
    for k=1:length(z)
        q = (z[k]-anm1)
        p = (one(Complex128)-conj(anm1)*z[k])
        C = c/(z[k]-an)
        PHI = p*phinm1[k]
        PHISTAR = q*phinm1star[k]
        phin[k] = C*(PHI+conj(L)*PHISTAR)
        phinstar[k] = C*(L*PHI+PHISTAR)
    end
    return (phin, phinstar)
end

function genZ(N::Int)
	return exp((im*2*pi/N)*[0:N-1])
end

function genBasis(N::Int,a::Vector{Complex128})
    z = genZ(N)
    M = length(a)+1
    a = cat(1,zero(Complex128),a)
    B = ones(Complex128,(N,M))/sqrt(N)
    BstarOld = ones(Complex128,N)/sqrt(N)
    for k=2:M
        (B[:,k],BstarNew) = phiNext(z,a[k],a[k-1],B[:,k-1],BstarOld)
        BstarOld = BstarNew
    end
    return B
end

function genTBasis(N::Int,a::Vector{Complex128})
	B = genBasis(N,a)
	T = fft(B,1)/sqrt(N)
    T = flipud(T)
    return T
end

function genA(N::Int)
	return (rand(N)).*exp(im*2*pi*rand(N))
end

function compress2D(ay,ax,D)
	(y,x) = size(D)
	W = genBasis(y,ay)
	T = genTBasis(x,ax)
	return ((W'*(D*T)),W,T)
end

function uncompress2D(X,W,T)
	return W*X*T'
end