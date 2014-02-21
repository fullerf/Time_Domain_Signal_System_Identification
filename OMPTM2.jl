function en(an,anm1,Ln)
	return sqrt((1-abs2(an))/((one(an)-abs2(anm1))*(one(an)-abs2(Ln))))
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
        an = copy(a[end])
        anm1 = copy(a[end-1])
        bProd = length(a)>=3 ? conjBlaschkeProd(z,[a[1:end-2]]) : ones(Complex128,length(z))
        f1 = [((z[i]-anm1)/(one(Complex128)-conj(an)*z[i]))*conj(phinm1[i]) for i=1:length(z)]
        f2 = [((z[i]-anm1)/(one(Complex128)-conj(an)*z[i]))*phinm1[i]*bProd[i] for i=1:length(z)]
        num = sum(f1)
        den = sum(f2)
        r = -num/den
    end
    return r
end

function phiNext(z::Vector{Complex128},a::Complex128,phinm1::Vector{Complex128},firstpole::Bool)
    phin = zeros(Complex128,length(phinm1))
    if firstpole
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
        bProd = length(a)>=3 ? blaschkeProd(z,[a[1:end-2]]) : ones(Complex128,length(z))
        for k=1:length(z)
            q = (z[k]-an)
            p = (one(Complex128)-conj(anm1)*z[k])
            C = c*p/q
            phin[k] = C*(phinm1[k]+conj(L)*bProd[k]*conj(phinm1[k]))
        end
    end
    return (phin, L)
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
        bProd = length(a)>=3 ? blaschkeProd(z,[a[1:end-2]]) : ones(Complex128,length(z))
        for k=1:length(z)
            q = (z[k]-an)
            p = (one(Complex128)-conj(anm1)*z[k])
            C = c*p/q
            phin[k] = C*(phinm1[k]+conj(L)*bProd[k]*conj(phinm1[k]))
        end
    end
    return (phin, L)
end

function phiNext(z::Vector{Complex128},a::Vector{Complex128},phinm1::Vector{Complex128},L::Complex128)
    phin = zeros(Complex128,length(phinm1))
    if length(a)==1
        an = copy(a[end])
        anm1 = zero(Complex128)
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
        c = en(an,anm1,L)
        bProd = length(a)>=3 ? blaschkeProd(z,[a[1:end-2]]) : ones(Complex128,length(z))
        for k=1:length(z)
            q = (z[k]-an)
            p = (one(Complex128)-conj(anm1)*z[k])
            C = c*p/q
            phin[k] = C*(phinm1[k]+conj(L)*bProd[k]*conj(phinm1[k]))
        end
    end
    return phin
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

function genZ(N::Int)
	return exp((im*2*pi/N)*[0:(N-1)])
end

function genZArc(N::Int,Div::Int)
    return exp((im*2*pi/(Div*N))*[0:(N-1)])
end

function genBasis(z::Vector,a::Vector{Complex128})
    P = zeros(Complex128,(length(z),length(a)+1))
    P[:,1] = ones(Complex128,length(z))/sqrt(length(z));
    L = zeros(Complex128,length(a))
    for k=1:length(a)
        (P[:,k+1],L[k]) = phiNext(z,[a[1:k]],P[:,k])
    end
    return (P[:,2:end],L)
end

function genBasisGivenL(z::Vector,a::Vector{Complex128},L::Vector{Complex128},N::Int)
    P = zeros(Complex128,(length(z),length(a)+1))
    P[:,1] = ones(Complex128,length(z))/sqrt(N);
    for k=1:length(a)
        P[:,k+1] = phiNext(z,[a[1:k]],P[:,k],L[k])
    end
    return P[:,2:end]
end