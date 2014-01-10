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
    return B[:,2:end]
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

function lsq2DTMobj(beta::Vector,grad::Vector,Nx::Int,Ny::Int,D::Array{Complex128,2})
    beta1 = beta[1:(2*Ny)]
    beta1c = [beta1[i]*exp(im*2*pi*beta1[i+1]) for i=1:2:length(beta1)]
    beta2 = beta[(2*Ny+1):end]
    beta2c = [beta2[i]*exp(im*2*pi*beta2[i+1]) for i=1:2:length(beta2)]
    (X,W,T) = compress2D(beta1c,beta2c,D)
    Dhat = uncompress2D(X,W,T)
    return sum(abs2(D[:]-Dhat[:]))
end

function lsq2DTMgenBounds(Nx::Int,Ny::Int,offset)
    betaLB = zeros(Float64,2*Nx+2*Ny)
    betaUB = ones(Float64,2*Nx+2*Ny)
    betaUB[1:2:end] -= offset
    return (betaLB,betaUB)
end

function lsq2DTMgenStart(Nx::Int,Ny::Int,offset::Float64)
    Nt = 2*Nx+2*Ny
    Nr = Nx+Ny
    betaGuess = zeros(Float64,Nt)
    betaGuess[1:2:end] = (one(Float64)-offset)*rand(Nr)
    betaGuess[2:2:end] = rand(Nr)
    return betaGuess
end

function lsq2DTMopt(Nx::Int,Ny::Int,D::Array{Complex128,2},gmethod,lmethod,ltol,gWallClockTime)
    #This function assumes that NLopt is in the namespace
    Nt = 2*Nx+2*Ny
    offset = 0.001
    opt = Opt(method,Nt)
    lopt = Opt(lmethod,Nt)
    ftol_rel!(lopt, ltol) #local stopping criteria
    maxtime!(opt,gWallClockTime) #global stopping critera
    (betalb,betaub) = lsq2DTMgenBounds(Nx,Ny,offset)
    beta0 = lsq2DTMgenStart(Nx,Ny,offset)
    local_optimizer!(opt,lopt)
    min_objective!(opt,(x,g)->lsq2DTMobj(x,g,Nx,Ny,D))
    lower_bounds!(opt,betalb)
    upper_bounds!(opt,betaub)
    (rss,betaopt,exitcode) = optimize(opt,beta0)
    betaopty = betaopt[1:(2*Ny)]
    betaoptx = betaopt[(2*Ny+1):end]
    betay = [betaopty[i]*exp(im*2*pi*betaopty[i+1]) for i=1:2:length(betaopty)]
    betax = [betaoptx[i]*exp(im*2*pi*betaoptx[i+1]) for i=1:2:length(betaoptx)]
    Wopt = genBasis(size(D,1),betay)
    Topt = genTBasis(size(D,2),betax)
    #Xopt = (W'*(D*T))
    return (betay,betax,Wopt,Topt,rss,exitcode)
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

function genOMPBasis(N::Int,Z::Vector{Complex128},vB::Vector,vBstar::Vector)
    z = genZ(N)
    M = length(Z)
    B = zeros(Complex128,(N,M))
    Bstar = zeros(Complex128,(N,M))
    for k=1:M
        (B[:,k],Bstar[:,k]) = phiNext(z,Z[k-1],zero(Complex128),vB,vBstar)
    end
    return (B,Bstar)
end

function genOMPTBasis(N::Int,Z::Vector{Complex128},vB::Vector,vBstar::Vector)
    (B,Bstar) = genOMPBasis(N,Z,vB,vBstar)
    T = fft(B,1)/sqrt(N)
    T = flipud(T)
    return (T,Bstar)
end

function indmaxmat(M)
    (x,y) = size(M)
    linind = indmax(M)
    j = linind%x
    k = fld(linind,y)+1
    return (j,k)
end

function lsq2DTMOMPobjc!(Z::Vector,D::Array{Complex128,2},vBY::Vector,vBstarY::Vector,vBX::Vector,vBstarX::Vector)
    N = length(Z)
    (by,vBstarYnew) = genOMPBasis(size(D,1),Z,vBY,vBstarY)
    (bx,vBstarXnew) = genOMPTBasis(size(D,2),Z,vBX,vBstarX)
    O = (by'*(D*bx))
    (j,k) = indmaxmat(abs2(O))
    vBY = by[:,j]
    vBstarY = vBstarYnew[:,j]
    vBX = bx[:,k]
    return (j,k,by[:,j],vBstarYnew[:,j],bx[:,k],vBstarXnew[:,k])
end

function lsq2DTMOMPresid!(byopt,bxopt,D::Array{Complex128,2})
    coeff = dot(squeeze(by,2),squeeze(D*bx,2))
    #return the residual of an omp step
    D -= coeff*by*bx'
    return coeff
end

function lsq2DTMOMPiterate(Niter::Int,Ngrid::Int,D::Array{Complex128,2})
    coeffs = zeros(Complex128,Niter)
    polesx = zeros(Complex128,Niter)
    polesy = zeros(Complex128,Niter)
    (x,y) = size(D)
    vBY = ones(Complex128,y)/sqrt(y)
    vBX = ones(Complex128,x)/sqrt(x)
    vBstarY = ones(Complex128,y)/sqrt(y)
    vBstarX = ones(Complex128,x)/sqrt(x)
    Z = unitDiskGrid(nextpow2(int(floor(sqrt(Ngrid)))),0.001)
    S = copy(D)
    for k=1:Niter
        (a,b,vBY,vBstarY,vBX,vBstarX) = lsq2DTMOMPobjc(Z,S,)
        coeffs[k] = lsq2DTMOMPresid!(Z[a],Z[b],S)
        polesy[k] = Z[a]
        polesx[k] = Z[b]
    end
    return (coeffs,polesx,polesy)
end