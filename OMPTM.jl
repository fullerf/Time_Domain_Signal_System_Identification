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

function ln(z::Vector{Complex128},a::Vector{Complex128},phinm1::Vector{Complex128})
    if length(a)==1
        an = copy(a[end])
        anm1 = zero(Complex128)
        phinm1star = phinm1
        f1 = [((z[i]-anm1)/(one(Complex128)-conj(an)*z[i]))*conj(phinm1[i]) for i=1:length(z)]
        f2 = [((one(Complex128)-conj(anm1)*z[i])/(one(Complex128)-conj(an)*z[i]))*phinm1star[i] for i=1:length(z)]
        num = sum(f1)
        den = sum(f2)
        r = -num/den
    else
        an = a[end]
        anm1 = a[end-1]
        bProd = length(a)>3 ? blaschkeProd(z,a[end-2]) : ones(Complex128,length(z))
        phinm1star = [phinm1[i]*bProd[i] for i=1:length(z)]
        f1 = [((z[i]-anm1)/(one(Complex128)-conj(an)*z[i]))*conj(phinm1[i]) for i=1:length(z)]
        f2 = [((z[i]-anm1)/(one(Complex128)-conj(an)*z[i]))*phinm1star for i=1:length(z)]
        num = sum(f1)
        den = sum(f2)
        r = -num/den
    end
    return r
end

function phiNext(z::Vector{Complex128},a::Vector{Complex128},phinm1::Vector{Complex128})
    phin = zeros(Complex128,length(phinm1))
    if length(a)==1
        an = copy(a[end])
        anm1 = zero(Complex128)
        L = ln(z,a,phinm1)
        c = en(an,anm1,L)
        for k=1:length(z)
            q = (z[k]-anm1)
            p = (one(Complex128)-conj(anm1)*z[k])
            C = c/(z[k]-an)
            PHI = p*phinm1[k]
            PHISTAR = q*phinm1[k]
            phin[k] = C*(PHI+conj(L)*PHISTAR)
        end
    else
        an = a[end]
        anm1 = a[end-1]
        L = ln(z,a,phinm1)
        c = en(an,anm1,L)
        bProd = length(a)>3 ? blaschkeProd(z,a[end-2]) : ones(Complex128,length(z))
        for k=1:length(z)
            q = (z[k]-an)
            p = (one(Complex128)-conj(anm1)*z[k])
            C = c*p/q
            phin[k] = C*(phinm1[k]+conj(L)*bProd[k]*conj(phinm1[k]))
        end
    end
    return (phin, L)
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
    return (phin, phinstar, L)
end

function phiNext(z::Vector{Complex128},an::Complex128,anm1::Complex128,phinm1::Vector{Complex128},phinm1star::Vector{Complex128},L::Complex128)
    phin = zeros(Complex128,length(phinm1))
    phinstar = zeros(Complex128,length(phinm1))
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
    return (phin, phinstar, L)
end

function blaschkeProd(z::Vector{Complex128},aj::Vector{Complex128})
    r = ones(Complex128,length(z))
    for i=1:length(aj)
        t = [(one(Complex128)-conj(aj[i])*z[k])./(z[k]-aj[i]) for k=1:length(z)]
        r = r.*t
    end
    return r
end

function conjBlaschkeProd(z::Vector{Complex128},aj::Vector{Complex128})
    r = ones(Complex128,length(z))
    for i=1:length(aj)
        t = [(z[k]-aj[i])/(one(Complex128)-conj(aj[i])*z[k]) for k=1:length(z)]
        r = r.*t
    end
    return r
end

function qrUpdate(z::Vector{Complex128},p::Vector{Complex128},B::Array{Complex128,2},Bstar::Array{Complex128,2})
    (L,C) = size(B)
    (Q,R) = qr(B)
    Lk = zeros(Complex128,C-1)
    Bnew = copy(Q)
    BstarNew = zeros(Complex128,(L,C))
    for k=2:C
        BstarNew[:,k] = conj(Bnew[:,k]).*blaschkeProd(z,[p[1:k]])
        Lk[k-1] = ln(z,p[k],p[k-1],Bnew[:,k-1],BstarNew[:,k-1])
    end
    return (Bnew[:,2:end],BstarNew[:,2:end],Lk)
end

function genZ(N::Int)
	return exp((im*2*pi/N)*[0:N-1])
end

function genZArc(N::Int,Div::Int)
    return exp((im*2*pi/(Div*N))*[0:N-1])
end

function genBasis(N::Int,z::Vector,a::Vector{Complex128})
    M = length(a)+1
    upfreq = 500;
    a = cat(1,zero(Complex128),a)
    B = ones(Complex128,(N,M))/sqrt(N)
    Bstar = ones(Complex128,(N,M))/sqrt(N)
    #BstarOld = ones(Complex128,N)/sqrt(N)
    Lk = zeros(Complex128,M-1)
    for k=2:M
        (B[:,k],Bstar[:,k],Lk[k-1]) = phiNext(z,a[k],a[k-1],B[:,k-1],Bstar[:,k-1])
        #(B[:,k],BstarNew,Lk[k-1]) = phiNext(z,a[k],a[k-1],B[:,k-1],BstarOld)
        #BstarOld = copy(BstarNew)
        if mod(k+1,upfreq)==0
            println("QR updating")
            (B[:,(k-upfreq+1):k],Bstar[:,(k-upfreq+1):k],Lk[(k-upfreq+1):k]) = qrUpdate(z,a,B[:,(k-upfreq):k],Bstar[:,(k-upfreq):k])
        end
    end
    return (B[:,2:end],Lk,Bstar[:,2:end])
end

function genBasis(N::Int,z::Vector,a::Vector{Complex128},Lk::Vector{Complex128},InitVal)
    M = length(a)+1
    a = cat(1,zero(Complex128),a)
    B = ones(Complex128,(N,M))*InitVal
    BstarOld = ones(Complex128,N)*InitVal
    for k=2:M
        (B[:,k],BstarNew) = phiNext(z,a[k],a[k-1],B[:,k-1],BstarOld,Lk[k-1])
        BstarOld = copy(BstarNew)
    end
    return B[:,2:end]
end

function genWBasis(N::Int,a::Vector{Complex128})
    z = genZ(N)
    B = genBasis(N,z,a)
    return fftshift(B,1)
end

function genTBasis(N::Int,a::Vector{Complex128})
    z = genZ(N)
	(B,L) = genBasis(N,z,a)
	T = fft(B,1)/sqrt(N)
    T = flipud(T)
    return T
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
    return (reshape(Zjk,N^2))
end

function ompTMFreqDomain(D::Array{Complex128,2},Z::Vector{Complex128},Npoles::Int,Npad::Int)
    D = cat(1,zeros(Npad,size(D,2)),D,zeros(Npad,size(D,2)))
    N = size(D,1)
    Q = size(D,2)
    M = length(Z)
    z = genZ(N)
    B = zeros(Complex128,(N,M))
    Bstar = zeros(Complex128,(N,M))
    Bnm1 = ones(Complex128,(N,M))/sqrt(N)
    Bstarnm1 = ones(Complex128,(N,M))/sqrt(N)
    poles = zeros(Complex128,Npoles+1)
    amps = zeros(Complex128,((Npoles+1),Q))
    jopt = one(Int);
    R = copy(D)
    for k=2:(Npoles+1)
        for j=1:M
            (B[:,j],Bstar[:,j]) = phiNext(z,Z[j],poles[k-1],Bnm1[:,jopt],Bstarnm1[:,jopt])
        end
        Bnm1 = copy(B)
        Bstarnm1 = copy(Bstar)
        B = fftshift(B,1)
        P = B'*R #projection on to column space of data
        O = sum(abs2(P),2) #this is the vector we wish to maximize
        jopt = indmax(O)
        poles[k] = Z[jopt]
        amps[k,:] = P[jopt,:]
        R -= B[:,jopt]*amps[k,:]
    end
    return (poles[2:end], amps[2:end,:])
end

function ompTMTimeDomain(D::Array{Complex128,2},Z::Vector{Complex128},Npoles::Int,Npad::Int)
    D = cat(1,zeros(Npad,size(D,2)),D,zeros(Npad,size(D,2)))
    N = size(D,1)
    Q = size(D,2)
    M = length(Z)
    z = genZ(N)
    B = zeros(Complex128,(N,M))
    Bstar = zeros(Complex128,(Npoles,M))
    Ltemp = zeros(Complex128,M)
    Lopt = zeros(Complex128,Npoles)
    Bnm1 = ones(Complex128,(N,M))/sqrt(N)
    Bstarnm1 = ones(Complex128,(N,M))/sqrt(N)
    poles = zeros(Complex128,Npoles+1)
    amps = zeros(Complex128,((Npoles+1),Q))
    jopt = ones(Int,Npoles+1);
    R = copy(D)
    for k=2:(Npoles+1)
        for j=1:M
            (B[:,j],Bstar[:,j],Ltemp[j]) = phiNext(z,Z[j],poles[k-1],Bnm1[:,jopt[k-1]],Bstarnm1[:,jopt[k-1]])
        end
        Bnm1 = copy(B)
        Bstarnm1 = copy(Bstar)
        B = flipud(fft(B,1))/sqrt(N)
        P = B'*R #projection on to column space of data
        O = sum(abs2(P),2) #this is the vector we wish to maximize
        O[jopt[1:(k-1)]] = 0 #prevent repeat pole selection
        jopt[k] = indmax(O)
        poles[k] = Z[jopt[k]]
        amps[k,:] = P[jopt[k],:]
        R -= B[:,jopt[k]]*amps[k,:]
        println(Ltemp[jopt[k]])
        Lopt[k-1] = Ltemp[jopt[k]]
    end
    return (poles[2:end], amps[2:end,:], Lopt)
end

function ompSepTMDecomp(Npoles::Int,Nsing::Int,Ngrid::Int,Ntarget::Int,D::Array{Complex128,2})
    (u,s,v) = svd(D)
    vp = v[:,1:Nsing]
    sp = s[1:Nsing]
    up  = u[:,1:Nsing]
    Z = unitDiskGrid(int(floor(sqrt(Ngrid))),0.001)
    (polesx,ampsx,Lk) = ompTMTimeDomain(vp,Z,Npoles,0)
    Bx = genTBasis(Ntarget,polesx)
    R = up*diagm(sp)*(Bx*ampsx)'
    return (polesx,ampsx,R,Lk)
end

