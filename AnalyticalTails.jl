#!/usr/bin/env julia

#AnalyticalTails.jl
using DifferentialEquations
using OrdinaryDiffEq
using Distributions
using LossFunctions
using NLopt
using SpecialFunctions
using Cubature
#=
To Dos:
1. Compare to built-in adjoint sensitivity analysis: http://docs.juliadiffeq.org/latest/analysis/sensitivity.html#Adjoint-Sensitivity-Analysis-1
2. Why CVODE vs. the Krylov methods?
=#

#Kd is not being optimized for some reason

#Value to beat for NM_023516 is 260.1704, which is the Rscript LSODE value.

#Extra arguments aren't needed for high level languages because they have closures (e.g. lambda in python or -> in Julia)
#

tick = [0.0 0.0]
gradRatio = ones(7)
function Simulation(vecPars,ts)
    #vecPars = log.([180, 1e-7, 1.0, 1.0,10.0, 250.0, 6.0])
    #=
    Math for this section:
    overall equation is: Xdot(t) = A * X(t)+B
    Solution is: X(t) = Ainv * (matExp(A * t) * B-Identity * B)
    Or: Solve(A * X(t) = matExp(A * t) * B-Identity * B)
    =#
    # ts now can be an array as of OrdinaryDiffEq 3.21.0

	st     = vecPars[1] #why is julia 1 based for arrays, maybe it's like c
    ka     = vecPars[2]
    kd     = vecPars[3]
    kb     = vecPars[4]
    sz     = vecPars[5]
    lc     = vecPars[6]
    sc     = vecPars[7]

    intMaxTail = 249

    #Generate matrix A
    vecDiag = [-kb/(1 + exp(-(x-lc)/sc))-kd for x in 0:intMaxTail]
    vecSubDiag = zeros(intMaxTail) + kd
    arrA = Bidiagonal(vecDiag,vecSubDiag,false)


    #Construct the source matrix
    vecB = Float64[ka*pdf(NegativeBinomial(sz,sz/(sz+st)),x) for x=intMaxTail:-1:0]
	# vecB = vecB[:,:] #better way to do this

	#solve via phiv and krylov
    # arrU = OrdinaryDiffEq.expv_timestep(ts,arrA,vecB;adaptive=true)	# important to set the time step tau? all time is used here. Verbose = true? 
    # #what is this actually solving? Why does it start with vecB and then decrease without a source term?
    # arrB = repmat(vecB,1,length(ts))
    # arrX = arrA\(arrU - arrB)

    # #Solve the steady-state time point and add it to the solution array
    # arrSS = arrA\-vecB
    # arrX = [arrX arrSS]

    # return collect(flatten(arrX'))
    
    #This is a basic banded solver that works pretty well
    T = ts[length(ts)]
    fX(vecX,p,t) = arrA*vecX+vecB
    vecX0 = zeros(intMaxTail+1)
    vecTspan = (0.0,T)
    objProb = ODEProblem(fX,vecX0,vecTspan)
    objX = solve(objProb,CVODE_BDF(linear_solver=:Band,jac_upper=1,jac_lower=1),dense=true,abtol=1e-16)
    # objX(t) = arrA\(expm(Matrix(arrA*t))*vecB - vecB)

    return(objX,arrA)
    #=
    I shouldn't be computing steady state like I am here. There's a much simpler way:

    arrSS = arrA\-vecB

    This comes from simply setting dx/dt to 0 in the original differential equation, much simpler computation, but it requires a bifurcation of the code.
    =#

    # plot(250:-1:0,arrX)

end

function dF(vecPars,intMaxTail)
    # Compute the derivatives of each parameter and return those to the adjoint system as a matrix when called.
    st     = vecPars[1] #starting tail length
    ka     = vecPars[2] #production rate
    kd     = vecPars[3] #deadenylation rate
    kb     = vecPars[4] #decapping rate
    sz     = vecPars[5] #spread for the production rate
    lc     = vecPars[6] #decapping rate position
    sc     = vecPars[7] #decapping rate spread

    #derivate with respect to kb
    arrAdb = Diagonal([-(1 + exp.(-(x-lc)/sc))^(-1.0) for x=0:(intMaxTail)])
    
    #derivate with respect to kd
    vecDiagdkd = -ones(intMaxTail+1)
    vecSubDiagdkd = ones(intMaxTail)
    arrAdkd = Bidiagonal(vecDiagdkd,vecSubDiagdkd,false)
    
    #derivative with respect to lc
    arrAdlc = Diagonal([kb*exp.((lc-x)/sc)/((1+exp.((lc-x)/sc))^2*sc) for x=0:(intMaxTail)])
    
    #derivative of sc
    arrAdsc = Diagonal([-kb*(exp.((lc - x)/sc)*(lc - x))/((1 + exp.((lc - x)/sc))^2*sc^2) for x=0:(intMaxTail)])
    
    #derivative of source with respect to ka
    vecBdka = Float64[pdf(NegativeBinomial(sz,sz/(sz+st)),xEval) for xEval=intMaxTail:-1:0]

    #derivative of source with respect to st
    vecBdst = Float64[(st-xEval)/(st*(st+sz))*(-sz)*ka*pdf(NegativeBinomial(sz,sz/(sz+st)),xEval) for xEval=intMaxTail:-1:0] 
    #derivative of source with respect to sz
    vecBdsz = zeros(intMaxTail+1)
    for xEval in intMaxTail:-1:0
        f64Term1                        = st-xEval+(st+sz)*log(sz/(st+sz))
        f64Term2                        = -(st+sz)*digamma(sz)+(st+sz)*digamma(sz+xEval)
        f64Bdsz                         = 1/(st+sz)*ka*pdf(NegativeBinomial(sz,sz/(sz+st)),xEval)*(f64Term1+f64Term2)
        vecBdsz[intMaxTail-xEval+1]     = f64Bdsz
    end

    return(vecBdst,vecBdka,arrAdkd,arrAdb,vecBdsz,arrAdlc,arrAdsc)
end

function Adjoint(vecPars,objX,arrA,arrResidualsdX,ts)
    #Run the adjoint system, with all of the integration for each time point stepwise. 
    intTsNo = length(ts)
    T = ts[intTsNo]
    intMaxTail = 249
    vecTspan = (0.0,T)
    fLambda(vecX,p,t) = arrA'*vecX

    objProb = ODEProblem(fLambda,ones(intMaxTail+1),vecTspan)
    objLambda = solve(objProb,CVODE_BDF(linear_solver=:Band,jac_upper=1,jac_lower=1),dense=true,abstol = 1e-16)
    (vecBdst,vecBdka,arrAdkd,arrAdb,vecBdsz,arrAdlc,arrAdsc) = dF(vecPars,intMaxTail)

    function initialVal(t)
        if t < 40.0
            return(arrResidualsdX[:,1])
        elseif 40.0 < t <= 60.0
            return(arrResidualsdX[:,2])
        elseif 60.0 < t <= 120.0
            return(arrResidualsdX[:,3])
        elseif 120.0 < t <= 240.0
            return(arrResidualsdX[:,4])
        elseif 240.0 < t <= 480.0
            return(arrResidualsdX[:,5])
        end
    end

    objFp(t) = -[vecBdst vecBdka arrAdkd*objX(t) arrAdb*objX(t) vecBdsz arrAdlc*objX(t) arrAdsc*objX(t)]


    function objIntegrand(t,v)
        v[:] = (objLambda(T-t).*initialVal(T-t))'*objFp(T-t)
    end

    #what kind of tolerance do I need in this integration?
    vecGradient  = hquadrature(length(vecPars),(t,v)-> objIntegrand(t,v),T,240.0,abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars),(t,v)-> objIntegrand(t,v),240.0,120.0,abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars),(t,v)-> objIntegrand(t,v),120.0,60.0,abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars),(t,v)-> objIntegrand(t,v),60.0,40.0,abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars),(t,v)-> objIntegrand(t,v),40.0,0.0,abstol=1e-16)[1] 

    # reshape(collect(Iterators.flatten([objIntegrand(x) for x in 1:10])),7,10)'
    #plot(1:40,reshape(collect(Iterators.flatten([objIntegrand(x) for x in 1:40])),7,40)')
    return(vecGradient)
end

function loss!(vecPars,grad,arrData,ts) # ! is used because this function modifies its argument
    global tick #global variables to monitor the loss
    global gradRatio

    #variance weighing the residuals
    # vecWeigh = 1.0./repeat([2.58757e-16,3.053586e-16,2.92142e-15,1.16722e-14,2.772866e-14,2.394499e-13],inner = 250) 
    vecWeigh = [2.58757e-16,3.053586e-16,2.92142e-15,1.16722e-14,2.772866e-14].^-1 #no steady state
    # vecWeigh = [2.58757e-16].^-1 #only 40min
    # println("pars: $vecPars")
    (objX, arrA) = Simulation(exp.(vecPars),ts) #this code now returns an object that can be evaluated at any time step. It's a time by species array
    arrPred = objX(ts)[:,:]
    # floatResidual = sum(value(L2DistLoss(), arrData, arrPred).*vecWeigh')
    floatResidual = sum((arrPred-arrData).^2.*vecWeigh') #is this right?
    tick = vcat(tick,[tick[size(tick)[1],1]+1 floatResidual])

    arrResidualsdX = (arrPred-arrData).*vecWeigh'.*2
    # println("Residual: $floatResidual")
    # if length(grad)>0 #finite differences gradient
    #     grad[:] =     Calculus.gradient(x -> loss!(x,[],arrData,ts)[1],vecPars) #need to check this against a finite differences gradient
    #     #what about the exponentiation of the parameters?
    #     # println("Gradient: $grad")
    # end
    if length(grad)>0
        println("Tick: $(tick[size(tick)[1],1]), Residual: $floatResidual")
        gradAdjoint = Adjoint(exp.(vecPars),objX,arrA,arrResidualsdX,ts).*exp.(vecPars)
        grad[:] = gradAdjoint
        # gradFiniteDiff = Calculus.gradient(x -> loss!(x,[],arrData,ts)[1],vecPars)
        # println("FiniteDiff: $gradFiniteDiff") #need to check this against a finite differences gradient
        #what about the exponentiation of the parameters?
        # println("Gradient: $grad")
        # println("GradRatio: $(gradFiniteDiff./gradAdjoint)")
        # gradRatio = hcat(gradRatio,gradFiniteDiff./gradAdjoint)
    end
    return (floatResidual,grad)
end
# function gradient(pars, data)
# end
function optim(vecPars,arrData,ts)
    objOpt = Opt(:LD_LBFGS,7) #optimization object, with parameters tuned below. 
    # objOpt = Opt(:LN_COBYLA,7)
    min_objective!(objOpt,(vecPars,grad) -> loss!(vecPars,grad,arrData,ts)[1]) #some issue here
    lower_bounds!(objOpt, log.([1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10]))
    upper_bounds!(objOpt, log.([250,1e2,1e1,1e1,100.0,500,500]))
    xtol_rel!(objOpt,1e-8)
    ftol_abs!(objOpt,1e-3) #this is the function tolerance, it needs to be set if the function returns values that fluctuate.
    xtol_rel!(objOpt,1e-8)
    maxeval!(objOpt, 200) #this can't be the final maxeval
    (minf,minx,ret) = optimize(objOpt, vecPars)
end

function io()
	#CSV.read can't read in the data, need another function
    tupData = readdlm("/lab/solexa_bartel/teisen/RNAseq/Scripts/JuliaScripts/TEST5_genes_miR-155_minus_sample_means_background_subtracted_v6_st_global.txt",'\t',Float64,'\n';header=true)	
    # arrData = collect(flatten(tupData[1][1:6,1:250])) #only first gene
    # arrData = tupData[1][1:5,1:250]' #only first gene, don't flatten, no steady state
    arrData = tupData[1][1:5,1:250]' #only first gene, don't flatten, only 40min, only 25 vals
	strAcc = tupData[2]
    println("Gene being optimized:\t",strAcc[1])
    return arrData
end

function main()
    ts = [40.0,60.0,120.0,240.0,480.0]
    # ts = [480.0]
    vecPars = log.([140, 1e-7, .1, .1,10.0, 264.0, 6.0])
    arrData = io()

    arrFiniteDiffGrad = Calculus.gradient(x -> loss!(x,[],arrData,ts)[1],vecPars)
    arrAnalyticalGrad = loss!(vecPars,[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],arrData,ts)[2]
    arrCompare = [arrFiniteDiffGrad arrAnalyticalGrad vecPars arrFiniteDiffGrad./(arrAnalyticalGrad)]
    optim(vecPars,arrData,ts) #log of pars

end

@time main()
