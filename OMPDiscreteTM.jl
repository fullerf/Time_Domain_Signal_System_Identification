function OMPTimeDomain(V::Array{Complex128,2},Precision::Int64,OrthoTol::Float64,gmethod,lmethod,ltol,gWallClockTime)
	set_bigfloat_precision(Precision)
	#Add a zero at the beginning of the trace to put make it 
	#compatible with a strictly proper transfer function
	#which is what we are using to fit the data
	cat(1,zeros(Complex128,(1,size(V,2))),V)
	(L,N) = size(V)
	zarb = genZarb(N)
	zdbl = convert(Vector{Complex128},zarb)
	#pole vectors one at arbitrary precsion, one at dbl
	parb = zeros(ComplexPair{BigFloat},N-1)
	pdbl = zeros(Complex128,N-1)
	AICc = zeros(Float64,N-1)
	B = ones(ComplexPair{BigFloat},(L,N))/sqrt(L)
	Bdbl = ones(Complex128,(L,N))/sqrt(L)
	#Set up some fitting parameters for NLOPT
	opt = Opt(gmethod,2)
    lopt = Opt(lmethod,2)
    ftol_rel!(lopt, ltol) #local stopping criteria
    maxtime!(opt,gWallClockTime) #global stopping critera
    betalb = [0.01 0]
    betaub = [0.99 (2*pi-eps(Float64))]
    lower_bounds!(opt,betalb)
    upper_bounds!(opt,betaub)
	terminate_fit = 0;
	k = 0;
	Resid = copy(V)
	while terminate_fit<1
		k += 1
		beta_init = genStart(betalb,betaub)
		if k==1
			firstpole = true
		else
			firstpole = false
		end
		Bnm1 = convert(Vector{Complex128},B[:,k])
		max_objective!(opt,(x,g)->ompObj(x,g,zdbl,Bnm1,Resid,firstpole))
		(oVal,betaOpt,exitcode) = optimize(opt,beta0)
		parb[k] = convert(ComplexPair{BigFloat},betaOpt[1]*exp(im*betaOpt[2]))
		pdbl[k] = convert(Complex128,parb[k])
		Bnext = phiNext(zdbl,[pdbl[1:k]],Bnm1)
		B[:,k+1] = phiNext(zarb,[parb[1:k]],B[:,k])
		Bdbl[:,k+1] = convert(Vector{Complex128},B[:,k+1])
		OrthogonalityMetric = maximum(abs(abs(Bnext'*Bdbl[:,1:k+1]) - cat(2,zeros(Float64,k),one(Float64))))
		if OrthogonalityMetric>OrthoTol
			#terminate the fit if we lose orthogonality when operating at double precision
			terminate_fit = 2; 
		end
		Resid = V-Bdbl[:,2:k+1]*Bdbl[:,2:k+1]'*V
		RSS = norm(Resid)
		AICc[k] = akaikeInformationCriterion(RSS,2*k+1,2*L)
		if k==(N-1)
			terminate_fit = 1
		end
	end
	optind = indmin(AICc)
	A = Bdbl[:,1:optind]'*V
	return (Bdbl[:,1:optind],A)
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
		bstart[k] = (ub[k]-lb[k])*rand(1)+lb[k]
	end
end


function ompObj(beta::Vector,grad::Vector,z::Vector{Complex128},Bnm1::Vector{Complex128},D::Array{Complex128,2},firstpole::Bool)
    R = beta[1];
    T = beta[2];
    p = R*exp(im*T)
    Bn = phiNext(z,p,Bnm1,firstpole)
    Proj = Bn'*D;
    return sum(abs(Proj),2)
end

function lsq2DTMopt(N::Int,Bnm1,gmethod,lmethod,ltol,gWallClockTime)
    #This function assumes that NLopt is in the namespace
    Nt = 2*Nx+2*Ny
    offset = 0.001
    opt = Opt(gmethod,Nt)
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