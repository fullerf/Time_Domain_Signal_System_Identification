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

function psiNextDiscreteMeasure(p::Vector,psi::Vector,psistar::Vector,Lk)
    N = length(psi)
    psik = zero(psi)
    psistark = zero(psistar)
    if length(p)==1
        p = cat(1,zero(p[1]),p)
        c = sqrt(1-abs2(p[end]))
    else
        c = sqrt((1-abs2(p[end]))/(1-abs2(p[end-1])))
    end
    for j=2:(N)
        psik[j] = c*sum([p[end]^(j-i)*(psi[i-1]-conj(p[end-1])*psi[i]) for i=2:j])
        psik[j] += conj(Lk)*c*sum([p[end]^(j-i)*(psistar[i]-p[end-1]*psistar[i-1]) for i=2:j])
        psistark[j] = c*sum([p[end]^(j-i)*(psistar[i]-p[end-1]*psistar[i-1]) for i=2:j])
        psistark[j] += Lk*c*sum([p[end]^(j-i)*(psi[i-1]-conj(p[end-1])*psi[i]) for i=2:j])
    end
    return (psik, psistark)
end

function unitStepFun(k::Int)
    if k>=0
        return 1
    else
        return 0
    end
end

function truncDot(x::Vector,y::Vector)
    N = length(x)
    M = length(y)
    if N<M
        return dot(x,y[1:N])
    elseif M<N
        return dot(x[1:M],y)
    else
        return dot(x,y)
    end
end

function xi(an,anm1,kmax::Int)
    return [(an^(k-1)*unitStepFun(k-1)-conj(anm1)*(an^k)*unitStepFun(k)) for k=0:(kmax-1)]
end

function xistar(an,anm1,kmax::Int)
    return [(an^(k)*unitStepFun(k)-anm1*(an^(k-1)*unitStepFun(k-1))) for k=0:(kmax-1)]
end

function lk(an,anm1,psinm1::Vector,psistarnm1::Vector,u::Vector)
    xiv = xi(an,anm1,length(psinm1))
    xistarv = xistar(an,anm1,length(psinm1))
    cnum = ifft(fft(xiv).*fft(psinm1))
    cden = ifft(fft(xistarv).*fft(psistarnm1))
    num = truncDot(u,cnum)
    den = truncDot(u,cden)
    return -num/den
end