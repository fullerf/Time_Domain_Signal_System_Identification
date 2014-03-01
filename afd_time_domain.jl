function psiNext(p::Vector,psi::Vector)
    psi = cat(1,zero(Complex128),psi)
    N = length(psi)
    psik = zero(psi)
    if length(p)==1
        p = cat(1,zero(p[1]),p)
        c = sqrt(1-abs2(p[end]))
        #in this case we assume psi = u(t), rather than the previous psi
        for j=0:(N-1)
            psik[j+1] = c*sum([p[end]^(i)*psi[j-i+1] for i=0:(j-1)])
        end
    else
        c = sqrt((1-abs2(p[end]))/(1-abs2(p[end-1])))
        for j=0:(N-1)
            psik[j+1] = c*sum([p[end]^(i)*(psi[j-i]-conj(p[end-1])*psi[j-i+1]) for i=0:(j-1)])
        end
    end
    return psik[2:end]
end

function genPsiBasis(p::Vector,u::Vector)
    N = length(p)
    B = zeros(Complex128,(length(u),N+1))
    B[:,1] = u
    for k=1:N
        B[:,k+1] = psiNext([p[1:k]],B[:,k])
    end
    return B[:,2:end]
end

function pulseFun(k::Int)
    return cat(1,one(Complex128),zeros(Complex128,k-1))
end

function unitStep(k::Int)
    return k>=0 ? 1 : 0
end

function dictElem(p::Complex128,K::Int)
    x = [(sqrt(1-abs2(p))*(p^(k-1))*unitStep(k-1))::Complex128 for k=0:(K)]
    return x[2:end]
end

function filterR(p::Complex128,R::Vector)
    ai = 1/conj(p)
    R = cat(1,zero(Complex128),R)
    N = length(R)
    Rkp1 = zero(R)
    for k=0:(N-1)
        Rkp1[k+1]  = sum([ai^(i+1)*(p*R[k-i]-R[k-i+1]) for i=0:(k-1)])
    end
    return Rkp1[2:end]
end

function dictObj(x,g,target::Vector)
    p = x[1]*exp(im*2*pi*x[2])
    trial = dictElem(p,length(target))
    return abs2(dot(trial,target))
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
    return reshape(Zjk,N^2)
end

function afdT(target::Vector,Npoles::Int,Nstarts::Int)
    Z = unitDiskGrid(int(sqrt(Nstarts)),0.2)
    opt = Opt(:LN_SBPLX,2)
    ftol_rel!(opt,1e-6) #local stopping criteria
    lb = [0.01; 0]
    ub = [0.95; 1-eps(Float64)]
    lower_bounds!(opt,lb)
    upper_bounds!(opt,ub)
    avec = zeros(Complex128,(length(Z),2))
    ovec = zeros(Float64,length(Z))
    p = zeros(Complex128,Npoles)
    a = zeros(Complex128,Npoles)
    gk = copy(target)
    for j=1:Npoles
        max_objective!(opt,(x,g)->dictObj(x,g,gk))
        for k=1:length(Z)
            beta0 = [abs(Z[k]); ((angle(Z[k])+pi)/(2*pi))]
            (ovec[k],avec[k,:],e) = optimize(opt,beta0)
        end
        aopt = avec[indmax(ovec),:]
        p[j] = aopt[1]*exp(im*2*pi*aopt[2])
        ep = dictElem(p[j],length(target))
        a[j] = dot(ep,gk)
        gk = filterR(p[j],(gk-a[j]*ep))
    end
    return (p,a)
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



