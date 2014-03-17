function psiNextLebesgueMeasure(p::Vector,psi::Vector)
    N = length(psi)
    psik = zero(psi)
    if length(p)==1
        p = cat(1,zero(p[1]),p)
        c = sqrt(1-abs2(p[end]))
    else
        c = sqrt((1-abs2(p[end]))/(1-abs2(p[end-1])))
    end
    for j=2:(N)
        psik[j] = c*sum([p[end]^(j-i)*(psi[i-1]-conj(p[end-1])*psi[i]) for i=2:j])
    end
    return psik
end

function psiBasis(p::Vector,u::Vector)
    N = length(u)
    B = zeros(Complex128,(N,length(p)+1))
    B[:,1] = u
    for k=1:length(p)
        B[:,k+1] = psiNextLebesgueMeasure([p[1:k]],B[:,k])
    end
    return B[2:end,2:end]
end

function dictObjective(x::Vector,g::Vector,target::Vector,u::Vector)
    p = x[1]*exp(im*2*pi*x[2])
    trial = psiNextLebesgueMeasure([p],u)[2:end]
    return abs2(dot(target,trial))
end

function unitStep(k)
    return k>=0 ? 1 : 0
end

function invBlaschkeElem(p::Complex128,K::Int)
    return [((conj(p)^(-i))*unitStep(i-1)-conj(p)^(-i-1))::Complex128 for i=0:(K-1)]
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

function afdT(target::Vector,u::Vector,Npoles::Int,Nstarts::Int,lmethod,ltol)
    Zstarts = unitDiskGrid(int(sqrt(Nstarts)),0.3)
    opt = Opt(lmethod,2)
    ftol_rel!(opt, ltol) #local stopping criteria
    lower_bounds!(opt,[0.01; 0])
    upper_bounds!(opt,[0.9;  1-eps(Float64)])
    p = zeros(Complex128,Npoles)
    R = copy(target)
    for j=1:Npoles
        ovec = zeros(Float64,length(Zstarts))
        avec = zeros(Float64,(length(Zstarts),2))
        max_objective!(opt,(x,g)->dictObjective(x,g,R,u))
        for k=1:length(Zstarts)
            beta1 = abs(Zstarts[k]);
            beta2 = (angle(Zstarts[k])+pi)/(2*pi);
            beta0 = [beta1; beta2]
            (ovec[k],avec[k,:],e) = optimize(opt,beta0)
        end
        aopt = avec[indmax(ovec),:]
        p[j] = aopt[1]*exp(im*2*pi*aopt[2])
        B = psiBasis([p[1:j]],u)
        R = (eye(length(target))-B*pinv(B))*target
        Q = invBlaschkeElem(p[j],(length(target))*4)
        R = oneSidedConv(R,Q)[1:length(target)]
    end
    return p
end

function oneSidedConv(u::Vector,v::Vector)
    U = fft(u)
    V = fft(v)
    if length(U)<=length(V)
        return ifft(V.*cat(1,U,zeros(Complex128,length(V)-length(U))))
    else
        return ifft(U.*cat(1,V,zeros(Complex128,length(U)-length(V))))
    end
end