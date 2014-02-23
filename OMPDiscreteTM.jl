function OMPTimeDomain(V::Array{Complex128,2},Precision::Int64,OrthoTol::Float64,gmethod,lmethod,ltol,gWallClockTime)
	set_bigfloat_precision(Precision)
	#Add a zero at the beginning of the trace to put make it 
	#compatible with a strictly proper transfer function
	#which is what we are using to fit the data
	V = cat(1,zeros(Complex128,(1,size(V,2))),V)
	L = size(V,1)
	N = size(V,1)
	zarb = genZArb(L)
	zdbl = convert(Vector{Complex128},zarb)
	#pole vectors one at arbitrary precsion, one at dbl
	parb = zeros(ComplexPair{BigFloat},N-1)
	Larb = zeros(ComplexPair{BigFloat},N-1)
	Ldbl = zeros(ComplexPair{Float64},N-1)
	pdbl = zeros(Complex128,N-1)
	AICc = zeros(Float64,N-1)
    RSS = zeros(Float64,N-1)
	B = ones(ComplexPair{BigFloat},(L,N))/sqrt(L)
	Bdbl = ones(Complex128,(L,N))/sqrt(L)
	#Set up some fitting parameters for NLOPT
	opt = Opt(gmethod,2)
    lopt = Opt(lmethod,2)
    ftol_rel!(lopt, ltol) #local stopping criteria
    maxtime!(opt,gWallClockTime) #global stopping critera
    betalb = [0.01; 0]
    betaub = [0.99; (2*pi-eps(Float64))]
    lower_bounds!(opt,betalb)
    upper_bounds!(opt,betaub)
	terminate_fit = 0;
	k = 0;
	Resid = copy(V)
	while terminate_fit<1
		k += 1
		if k==1
			anm1 = zeros(Complex128,0)
		else
			anm1 = [pdbl[k-1]]
		end
		Bnm1 = convert(Vector{Complex128},B[:,k])
		println(k)
		max_objective!(opt,(x,g)->ompObj(x,g,zdbl,anm1,Bnm1,Resid))
		beta0 = genStart(betalb,betaub)
		(oVal,betaOpt,exitcode) = optimize(opt,beta0)
		parb[k] = convert(ComplexPair{BigFloat},betaOpt[1]*exp(im*betaOpt[2]))
		pdbl[k] = convert(Complex128,parb[k])
		(Bnext,Ldbl[k]) = phiNext(zdbl,[pdbl[1:k]],Bnm1)
		(B[:,k+1],Larb[k]) = phiNext(zarb,[parb[1:k]],B[:,k])
		Bdbl[:,k+1] = convert(Vector{Complex128},B[:,k+1])
		OrthogonalityMetric = maximum(abs(abs(Bnext'*Bdbl[:,1:k+1]) - cat(1,zeros(Float64,k),one(Float64))'))
		if OrthogonalityMetric>OrthoTol
			#terminate the fit if we lose orthogonality when operating at double precision
			terminate_fit = 2; 
		end
		Tdbl = ifft(Bdbl[:,2:k+1],1)*sqrt(L)
		Resid = V-Tdbl*Tdbl'*V
		RSS[k] = norm(Resid)
		AICc[k] = akaikeInformationCriterion(RSS[k],2*k+1,2*L*size(V,2))
		if k==(N-1)
			terminate_fit = 1
		end
	end
    optind = indmin(AICc)
    Tdbl = ifft(Bdbl[:,2:k+1],1)*sqrt(L)
    A = Tdbl[:,1:(optind)]'*V
	return (zarb,parb,Larb,A,terminate_fit)
end

function akaikeInformationCriterion(RSS,numParams,dataLen)
    #This gives the sample size corrected Akaike Information Criteria
    #commonly known as AICc in the literature
    K = numParams+1
    return dataLen*log(RSS/dataLen)+2*K+2*K*((K+1)/(dataLen-K-1))
end

function genStart(lb,ub)
	bstart = zero(lb)
	for k=1:length(lb)
		bstart[k] = (ub[k]-lb[k])*rand(1)[1]+lb[k]
	end
	return bstart
end


function ompObj(beta::Vector,grad::Vector,z::Vector{Complex128},anm1::Vector{Complex128},Bnm1::Vector{Complex128},D::Array{Complex128,2})
    R = beta[1];
    T = beta[2];
    p = [anm1; R*exp(im*T)]
    (Bn,L) = phiNext(z,p,Bnm1)
    Tn = ifft(Bn,1)*sqrt(length(z))
    Proj = Tn'*D;
    return sum(abs(Proj),2)[1]
end