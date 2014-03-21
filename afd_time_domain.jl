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

function psiNextRecursive(p::Vector,psi::Vector)
    psi = cat(1,zero(Complex128),psi)
    N = length(psi)
    psik = zero(psi)
    if length(p)==1
        p = cat(1,zero(p[1]),p)
        c = sqrt(1-abs2(p[end]))
        #in this case we assume psi = u(t), rather than the previous psi
        for j=1:(N-1)
            psik[j+1] = p[end]*psik[j]+c*psi[j+1]
        end
    else
        c = sqrt((1-abs2(p[end]))/(1-abs2(p[end-1])))
        for j=1:(N-1)
            psik[j+1] = p[end]*psik[j] + c*(psi[j] - conj(p[end-1])*psi[j+1])
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

function genPsiBasisRecursive(p::Vector,u::Vector)
    N = length(p)
    B = zeros(Complex128,(length(u),N+1))
    B[:,1] = u
    for k=1:N
        B[:,k+1] = psiNextRecursive([p[1:k]],B[:,k])
    end
    return B[:,2:end]
end

function pulseFun(k::Int)
    a = one(Complex128);
    return cat(1,a,zeros(Complex128,k-1))
end

function gaussianPulse(k::Int,fwhm::Float64)
    a = one(Complex128);
    x = [a*exp(-4*log(2)*(i^2/fwhm^2)) for i=1:k]
    return x/norm(x)
end

function unitStep(k::Int)
    return k>=0 ? 1 : 0
end

function dictElem(p::Complex128,K::Int)
    x = [(sqrt(1-abs2(p))*(p^(k-1))*unitStep(k-1))::Complex128 for k=0:(K)]
    return x[2:end]
end

function dictElemRecursive(p::Complex128,u::Vector)
    c = sqrt(1-abs2(p))
    u = cat(1,zero(Complex128),u)
    N = length(u)
    r = zero(u)
    for k=1:(N-1)
        r[k+1] = p*r[k]+c*u[k]
    end
    return r[2:end]
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

function filterRrecursive(p::Complex128,R)
    ai = 1/conj(p)
    Rp = cat(1,zeros(Complex128,(1,size(R,2))),R)
    N = size(Rp,1)
    Rkp1 = zero(R)
    for j=1:size(R,2)
        Rkp1j = zero(Rp[:,j])
        for k=1:(N-1)
            Rkp1j[k+1]  = ai*Rkp1[k]-ai*(Rp[k+1]-p*Rp[k])
        end
        Rkp1[:,j] = Rkp1j[2:end]
    end
    return Rkp1
end

function dictObj(x,g,target)
    p = x[1]*exp(im*2*pi*x[2])
    trial = dictElem(p,size(target,1))
    r = sum(abs2(trial'*target),2)
    return r[1]
end

function dictObjNbest(x,g,target)
    p = [(x[i]*exp(im*2*pi*x[i+1]))::Complex128 for i=1:2:length(x)]
    R = copy(target)
    for k=1:length(p)
        trial = dictElem(p[k],size(target,1))
        #R = filterR(p[k],R - trial*((trial'*R)[1]))
        R = filterRrecursive(p[k],R - trial*((trial'*R)[1]))
    end
    r = norm(R)
    return r
end

function nthRemainder(p,target)
    R = copy(target)
    L = size(target,2);
    for k=1:length(p)
        trial = dictElem(p[k],size(target,1))
        for j=1:L
            R[:,j] = filterRrecursive(p[k],R[:,j] - dot(trial,R[:,j])*trial)
        end
    end
    return R
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

function afdT(target::Array,Npoles::Int,Nstarts::Int)
    Z = unitDiskGrid(int(sqrt(Nstarts)),0.2)
    opt = Opt(:LN_SBPLX,2)
    ftol_rel!(opt,1e-6) #local stopping criteria
    lb = [0.01; 0]
    ub = [0.999; 1-eps(Float64)]
    lower_bounds!(opt,lb)
    upper_bounds!(opt,ub)
    avec = zeros(Complex128,(length(Z),2))
    ovec = zeros(Float64,length(Z))
    p = zeros(Complex128,Npoles)
    a = zeros(Complex128,(Npoles,size(target,2)))
    gk = copy(target)
    for j=1:Npoles
        max_objective!(opt,(x,g)->dictObj(x,g,gk)[1])
        for k=1:length(Z)
            beta0 = [abs(Z[k]); ((angle(Z[k])+pi)/(2*pi))]
            (ovec[k],avec[k,:],e) = optimize(opt,beta0)
        end
        aopt = avec[indmax(ovec),:]
        p[j] = aopt[1]*exp(im*2*pi*aopt[2])
        ep = dictElem(p[j],size(target,1))        
        a[j,:] = ep'*gk
        #for i=1:size(target,2)
        #    gk[:,i] = filterRrecursive(p[j],(gk-a[j,i]*ep))
        #end
        gk = filterRrecursive(p[j],(gk-ep*a[j,:]))
    end
    return (p,a)
end

function afdTCyclic(target::Array,polesIn::Vector,Nstarts::Int,CycleTol::Float64,CycleMax::Int)
    Z = unitDiskGrid(int(sqrt(Nstarts)),0.01)
    Npoles = length(polesIn)
    opt = Opt(:LN_SBPLX,2)
    ftol_rel!(opt,1e-12) #local stopping criteria
    lb = [0.01; 0]
    ub = [0.9999999999; 1-eps(Float64)]
    lower_bounds!(opt,lb)
    upper_bounds!(opt,ub)
    avec = zeros(Complex128,(length(Z),2))
    ovec = zeros(Float64,length(Z))
    p = zero(polesIn)
    prel = CycleTol+1;
    cycle_count = 0;
    while ((prel>CycleTol) && (cycle_count<CycleMax))
        cycle_count += 1
        for j=1:Npoles
            pold = copy(p)
            pprime = pold[(1:Npoles).!=j]
            gk = nthRemainder(pprime,target)
            max_objective!(opt,(x,g)->dictObj(x,g,gk)[1])
            for k=1:length(Z)
                beta0 = [abs(Z[k]); ((angle(Z[k])+pi)/(2*pi))]
                (ovec[k],avec[k,:],errcode) = optimize(opt,beta0)
            end
            aopt = avec[indmax(ovec),:]
            p[j] = aopt[1]*exp(im*2*pi*aopt[2])
            prel = norm(p-pold)/norm(pold)
            #println(prel)
        end
    end
    B = genPsiBasisRecursive(p,pulseFun(size(target,1)))
    a = B\target
    #println(cycle_count)
    return (p,a)
end


function poles2vec(p::Vector{Complex128})
    r = zeros(Float64,length(p)*2)
    r[1:2:end] = abs(p)
    r[2:2:end] = (angle(p)+pi)/(2*pi)
    return r
end

function afdTNbest(target::Array,Npoles::Int,gWallClockTime::Number,beta0::Vector)
    opt = Opt(:GN_MLSL_LDS,Npoles*2)
    lopt = Opt(:LN_SBPLX,Npoles*2)
    ftol_rel!(lopt,1e-6) #local stopping criteria
    maxtime!(opt,gWallClockTime) #global stopping critera
    lb = [0.01; 0]
    ub = [0.99999; 1-eps(Float64)]
    lb = repmat(lb,Npoles)
    ub = repmat(ub,Npoles)
    lower_bounds!(opt,lb)
    upper_bounds!(opt,ub)
    local_optimizer!(opt,lopt)
    p = zeros(Complex128,Npoles)
    a = zeros(Complex128,(Npoles,size(target,2)))
    gk = copy(target)
    min_objective!(opt,(x,g)->dictObjNbest(x,g,gk)[1])
    (r,popt,e) = optimize(opt,beta0)
    p = [(popt[i]*exp(im*2*pi*popt[i+1]))::Complex128 for i=1:2:(Npoles*2)]
    B = genPsiBasis(p,pulseFun(size(target,1)))
    a = B\target
    return (p,a)
end

function akaikeInformationCriterion(RSS,numParams,dataLen)
    #This gives the sample size corrected Akaike Information Criteria
    #commonly known as AICc in the literature
    K = numParams+1
    return dataLen*log(RSS/dataLen)+2*K+2*K*((K+1)/(dataLen-K-1))
end

function processDataStack(Y::Array{Complex128,2},minNpoles::Int,maxNpoles::Int,Nstarts::Int,Textrap::Int)
    pmat = {}
    amat = {}
    Bmat = {}
    for j=1:size(Y,2)
        AICc = {}
        ptemp   = {}
        atemp   = {}
        for k=minNpoles:maxNpoles
            (p,a) = afdT(Y[:,j],k,Nstarts)
            push!(ptemp,p)
            push!(atemp,a)
            B = genPsiBasis(p,pulseFun(size(Y,1)))
            Yhat = squeeze(B*a,2)
            R = (Y[:,j]-Yhat)
            RSS = abs(dot(R,R))
            aicval = akaikeInformationCriterion(RSS,4*k,2*size(Y,1))
            push!(AICc,aicval)
        end
        optInd = indmin(AICc)
        push!(pmat,ptemp[optInd])
        push!(amat,atemp[optInd])
        B = genPsiBasis(pmat[j],pulseFun(Textrap))
        push!(Bmat,B)
        println(optInd+minNpoles-1)
        println(AICc[optInd])
    end
    That = zeros(Complex128,(Textrap,size(Y,2)))
    for k=1:size(Y,2)
        That[:,k] = Bmat[k]*amat[k]
    end
    return (That,pmat,amat)
end

function processDataStackCyclic(Y::Array{Complex128,2},minNpoles::Int,maxNpoles::Int,maxNcycles::Int,Nstarts::Int,CycleTol::Float64,Textrap::Int)
    pmat = {}
    amat = {}
    Bmat = {}
    for j=1:size(Y,2)
        AICc = {}
        ptemp   = {}
        atemp   = {}
        for k=minNpoles:maxNpoles
            (pinit,ainit) = afdT(Y[:,j],k,Nstarts)
            (p,a) = afdTCyclic(Y[:,j],pinit,Nstarts,CycleTol,maxNcycles)
            push!(ptemp,p)
            push!(atemp,a)
            B = genPsiBasis(p,pulseFun(size(Y,1)))
            Yhat = B*a
            R = (Y[:,j]-Yhat)
            RR = R'*R
            RSS = abs(RR)[1]
            aicval = akaikeInformationCriterion(RSS,2*length(p)+2*length(a),2*size(Y,1))
            push!(AICc,aicval)
        end
        optInd = indmin(AICc)
        push!(pmat,ptemp[optInd])
        push!(amat,atemp[optInd])
        B = genPsiBasis(pmat[j],pulseFun(Textrap))
        push!(Bmat,B)
        println(optInd+minNpoles-1)
        println(AICc)
    end
    That = zeros(Complex128,(Textrap,size(Y,2)))
    for k=1:size(Y,2)
        That[:,k] = Bmat[k]*amat[k]
    end
    return (That,pmat,amat)
end

function processDataStackCyclicGlobal(Y::Array{Complex128,2},minNpoles::Int,maxNpoles::Int,maxNcycles::Int,Nstarts::Int,CycleTol::Float64,Textrap::Int)
    Bmat = {}
    AICc = {}
    ptemp   = {}
    atemp   = {}
    for k=minNpoles:maxNpoles
        (pinit,ainit) = afdT(Y,k,Nstarts)
        (p,a) = afdTCyclic(Y,pinit,Nstarts,CycleTol,maxNcycles)
        push!(ptemp,p)
        push!(atemp,a)
        B = genPsiBasis(p,pulseFun(size(Y,1)))
        Yhat = B*a
        R = (Y-Yhat)
        RSS = sum((abs(R'*R)))
        aicval = akaikeInformationCriterion(RSS,2*length(p)+2*length(a),2*length(Y))
        println(aicval)
        push!(AICc,aicval)
    end
    optInd = indmin(AICc)
    popt = ptemp[optInd]
    aopt = atemp[optInd]
    B = genPsiBasis(popt,pulseFun(Textrap))
    #println(optInd+minNpoles-1)
    #println(AICc[optInd])
    That = B*aopt
    return (That,popt,aopt)
end

function processDataStackMLSL(Y::Array{Complex128,2},minNpoles::Int,maxNpoles::Int,Nstarts::Int,SearchTime::Number,Textrap::Int)
    pmat = {}
    amat = {}
    Bmat = {}
    for j=1:size(Y,2)
        AICc = {}
        ptemp   = {}
        atemp   = {}
        for k=minNpoles:maxNpoles
            (pinit,ainit) = afdT(Y[:,j],k,Nstarts)
            (p,a) = afdTNbest(Y[:,j],k,SearchTime,poles2vec(pinit))
            push!(ptemp,p)
            push!(atemp,a)
            B = genPsiBasis(p,pulseFun(size(Y,1)))
            Yhat = B*a
            R = (Y[:,j]-Yhat)
            RSS = abs(dot(R,R))
            aicval = akaikeInformationCriterion(RSS,4*k,2*size(Y,1))
            push!(AICc,aicval)
        end
        optInd = indmin(AICc)
        push!(pmat,ptemp[optInd])
        push!(amat,atemp[optInd])
        B = genPsiBasis(pmat[j],pulseFun(Textrap))
        push!(Bmat,B)
        println(optInd+minNpoles-1)
        println(AICc[optInd])
    end
    That = zeros(Complex128,(Textrap,size(Y,2)))
    for k=1:size(Y,2)
        That[:,k] = Bmat[k]*amat[k]
    end
    return (That,pmat,amat)
end

function processDataStackNpole(Y::Array{Complex128,2},Nsing::Int,Npoles::Int,Nstarts::Int,Textrap::Int)
    pmat = {}
    amat = {}
    Bmat = {}
    for j=1:Nsing
        (p,a) = afdT(Y[:,j],k,Nstarts)
        B = genPsiBasis(p,pulseFun(size(Y,1)))
        push!(pmat,p)
        push!(amat,a)
        B = genPsiBasis(pmat[j],pulseFun(Textrap))
        push!(Bmat,B)
    end
    That = zeros(Complex128,(Textrap,Nsing))
    for k=1:Nsing
        That[:,k] = Bmat[k]*amat[k]
    end
    return (That,pmat,amat)
end

function refineDataStack(Y::Array{Complex128,2},Nsing::Int,pmat,SearchTime::Number,Textrap::Int)
    Bmat = {}
    pref = {}
    aref = {}
    for j=1:Nsing
        pin = poles2vec(pmat[j])
        (p,a) = afdTNbest(Y[:,j],length(pmat[j]),SearchTime,pin)
        push!(pref,p)
        push!(aref,a)
        B = genPsiBasis(pref[j],pulseFun(Textrap))
        push!(Bmat,B)
    end
    That = zeros(Complex128,(Textrap,Nsing))
    for k=1:Nsing
        That[:,k] = Bmat[k]*aref[k]
    end
    return (That,pref,aref)
end

function refineDataStackCyclic(Y::Array{Complex128,2},Nsing::Int,pmat,Nstarts::Int,CycleTol::Float64,CycleMax::Int,Textrap::Int)
    Bmat = {}
    pref = {}
    aref = {}
    for j=1:Nsing
        (p,a) = afdTCyclic(Y[:,j],pmat[j],Nstarts,CycleTol,CycleMax)
        push!(pref,p)
        push!(aref,a)
        B = genPsiBasis(pref[j],pulseFun(Textrap))
        push!(Bmat,B)
    end
    That = zeros(Complex128,(Textrap,Nsing))
    for k=1:Nsing
        That[:,k] = Bmat[k]*aref[k]
    end
    return (That,pref,aref)
end

function unfoldMatrix(M::Array,N::Int)
    uM = M[:,1:N]
    for k=2:int(size(M,2)/N)
        uM = cat(1,uM,M[:,((k-1)*N+1):(k*N)]);
    end
    return uM
end

function cubeMatrix(M::Array,N::Int)
    cM = zeros(typeof(M[1]),(N,size(M,2),int(size(M,1)/N)))
    for k=1:int(size(M,1)/N)
        cM[:,:,k] = M[((k-1)*N+1):(k*N),:]
    end
    return cM
end

function fftaxis(N::Int,TSpac=1.0)
    if abs(mod(N,2))>0
        a = [0:1:int((N-1)/2)]
        b = [-(N-1)/2:1:-1]
        r = cat(1,a,b)/(N*TSpac)
    else
        a = [0:1:int(N/2-1)]
        b = [-int(N/2):1:-1]
        r = cat(1,a,b)/(N*TSpac)
    end
    return r
end

function unfoldCube(M::Array)
    R = zeros(typeof(M[1]),(size(M,1)*size(M,2),size(M,3)))
    N = size(M,1);
    for j=1:size(M,2)
        R[((j-1)*N+1):(j*N),:] = squeeze(M[:,j,:],2);
    end
    return R
end

function optimalSVDTrunc(u,s,v)
    N = size(v,1)
    M = size(v,1)
    if isreal(u[1])
        f = 1
    else
        f = 2;
    end
    AICc = {}
    Y = u*diagm(s)*v';
    for k=1:(N-1)
        Yhat = u[:,1:k]*diagm(s[1:k])*v[:,1:k]'
        R = Y-Yhat
        RSS = abs(R'*R)[1]
        push!(AICc,akaikeInformationCriterion(RSS,k,f*N*M))
    end
    optInd = indmin(AICc)
    return (AICc,optInd)
end

function w2wn(w::Vector)
    c = 299.792458; #nm/fs
    return 1E7*w/(2*pi*c)
end
