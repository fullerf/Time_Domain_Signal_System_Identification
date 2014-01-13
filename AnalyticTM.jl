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