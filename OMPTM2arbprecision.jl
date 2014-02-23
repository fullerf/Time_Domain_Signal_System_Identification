function en(an,anm1,Ln)
	return sqrt((1-abs2(an))/((one(an)-abs2(anm1))*(one(an)-abs2(Ln))))
end

function ln(z::Vector{ComplexPair{BigFloat}},a::Vector{ComplexPair{BigFloat}},phinm1::Vector{ComplexPair{BigFloat}})
    if length(a)==1
        an = copy(a[end])
        anm1 = zero(ComplexPair{BigFloat})
        phinm1star = phinm1
        f1 = [((z[i]-anm1)/(one(ComplexPair{BigFloat})-conj(an)*z[i]))*conj(phinm1[i]) for i=1:length(z)]
        f2 = [((one(ComplexPair{BigFloat})-conj(anm1)*z[i])/(one(ComplexPair{BigFloat})-conj(an)*z[i]))*phinm1star[i] for i=1:length(z)]
        num = sum(f1)
        den = sum(f2)
        r = -num/den
    else
        an = copy(a[end])
        anm1 = copy(a[end-1])
        bProd = length(a)>=3 ? conjBlaschkeProd(z,[a[1:end-2]]) : ones(ComplexPair{BigFloat},length(z))
        f1 = [((z[i]-anm1)/(one(ComplexPair{BigFloat})-conj(an)*z[i]))*conj(phinm1[i]) for i=1:length(z)]
        f2 = [((z[i]-anm1)/(one(ComplexPair{BigFloat})-conj(an)*z[i]))*phinm1[i]*bProd[i] for i=1:length(z)]
        num = sum(f1)
        den = sum(f2)
        r = -num/den
    end
    return r
end

function phiNext(z::Vector{ComplexPair{BigFloat}},a::Vector{ComplexPair{BigFloat}},phinm1::Vector{ComplexPair{BigFloat}})
    phin = zeros(ComplexPair{BigFloat},length(phinm1))
    if length(a)==1
        an = copy(a[end])
        anm1 = zero(ComplexPair{BigFloat})
        L = ln(z,a,phinm1)
        c = en(an,anm1,L)
        for k=1:length(z)
            q = (z[k]-anm1)
            p = (one(ComplexPair{BigFloat})-conj(anm1)*z[k])
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
        bProd = length(a)>=3 ? blaschkeProd(z,[a[1:end-2]]) : ones(ComplexPair{BigFloat},length(z))
        for k=1:length(z)
            q = (z[k]-an)
            p = (one(ComplexPair{BigFloat})-conj(anm1)*z[k])
            C = c*p/q
            phin[k] = C*(phinm1[k]+conj(L)*bProd[k]*conj(phinm1[k]))
        end
    end
    return (phin, L)
end

function phiNext(z::Vector{ComplexPair{BigFloat}},a::Vector{ComplexPair{BigFloat}},phinm1::Vector{ComplexPair{BigFloat}},L::ComplexPair{BigFloat})
    phin = zeros(ComplexPair{BigFloat},length(phinm1))
    if length(a)==1
        an = copy(a[end])
        anm1 = zero(ComplexPair{BigFloat})
        c = en(an,anm1,L)
        for k=1:length(z)
            q = (z[k]-anm1)
            p = (one(ComplexPair{BigFloat})-conj(anm1)*z[k])
            C = c/(z[k]-an)
            PHI = p*phinm1[k]
            PHISTAR = q*phinm1[k]
            phin[k] = C*(PHI+conj(L)*PHISTAR)
        end
    else
        an = a[end]
        anm1 = a[end-1]
        c = en(an,anm1,L)
        bProd = length(a)>=3 ? blaschkeProd(z,[a[1:end-2]]) : ones(ComplexPair{BigFloat},length(z))
        for k=1:length(z)
            q = (z[k]-an)
            p = (one(ComplexPair{BigFloat})-conj(anm1)*z[k])
            C = c*p/q
            phin[k] = C*(phinm1[k]+conj(L)*bProd[k]*conj(phinm1[k]))
        end
    end
    return phin
end

function blaschkeProd(z::Vector{ComplexPair{BigFloat}},aj::Vector{ComplexPair{BigFloat}})
    r = ones(ComplexPair{BigFloat},length(z))
    for i=1:length(aj)
        t = [(one(ComplexPair{BigFloat})-conj(aj[i])*z[k])./(z[k]-aj[i]) for k=1:length(z)]
        r = r.*t
    end
    return r
end

function conjBlaschkeProd(z::Vector{ComplexPair{BigFloat}},aj::Vector{ComplexPair{BigFloat}})
    r = ones(ComplexPair{BigFloat},length(z))
    for i=1:length(aj)
        t = [(z[k]-aj[i])/(one(ComplexPair{BigFloat})-conj(aj[i])*z[k]) for k=1:length(z)]
        r = r.*t
    end
    return r
end

function genZArb(N::Int)
	return exp(convert(ComplexPair{BigFloat},(im*2*pi/N))*[0:(N-1)])
end

function genZArcArb(N::Int,Div::Int)
    return exp(convert(ComplexPair{BigFloat},(im*2*pi/(Div*N)))*[0:(N-1)])
end

function genBasis(z::Vector,a::Vector{ComplexPair{BigFloat}})
    P = zeros(ComplexPair{BigFloat},(length(z),length(a)+1))
    P[:,1] = ones(ComplexPair{BigFloat},length(z))/sqrt(length(z));
    L = zeros(ComplexPair{BigFloat},length(a))
    for k=1:length(a)
        (P[:,k+1],L[k]) = phiNext(z,[a[1:k]],P[:,k])
    end
    return (P[:,2:end],L)
end

function genBasisGivenL(z::Vector{ComplexPair{BigFloat}},a::Vector{ComplexPair{BigFloat}},L::Vector{ComplexPair{BigFloat}},N::Int)
    P = zeros(ComplexPair{BigFloat},(length(z),length(a)+1))
    P[:,1] = ones(ComplexPair{BigFloat},length(z))/sqrt(N);
    for k=1:length(a)
        P[:,k+1] = phiNext(z,[a[1:k]],P[:,k],L[k])
    end
    return P[:,2:end]
end


function unitDiskGrid(n::Int,delta::Float64)
    N = nextpow2(n)
    M = N>>1
    v = linspace(delta,one(Float64)-delta,M)
    Zjk = zeros(ComplexPair{BigFloat},N,N)
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