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
This code is designed to simulate a distribution of mRNA poly(A)-tail lengths over time.
It starts with an initial distribution of tail lengths in the form
of a source term (specified by a negative binomial distribution). Once present, these mRNA
are acted on by deadenylation, which decreases the tail length of each mRNA 
by one nt at a time. A logistic distribution then serves to destroy mRNAs once 
they reach short tail lengths.

Overall, the companion matrix is a bidiagonal matrix. This code simulates the solution to
these differential equations and uses an optimization routine (L-BFGS-B) along with an 
analytical gradient computed from an adjoint system. While the code is only set up for
a single gene, the goal is to fit 3000 genes in parallel. To do this I'm working on efficiency
and numerical stability of the derivates.
=#


#=
To Dos:
1. Compare to built-in adjoint sensitivity analysis: 
    http://docs.juliadiffeq.org/latest/analysis/sensitivity.html#Adjoint-Sensitivity-Analysis-1
2. Updated expv packages may allow me to stop using CVODE?
3. Why is this system unstable w.r.t. the parameter kd?
=#

#Value to beat for NM_023516 is 260.1704, which is the Rscript LSODE value.


#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()

tick = [0.0 0.0]
gradRatio = ones(7)

function Simulation(vecPars,ts)
    #=
    Takes as input the parameter vector and time steps, and returns an object.
    This object is a hermite polynomial of the solution to the ODE, with t as
    an argument. 

    Math for this section:
    overall equation is: Xdot(t) = A * X(t)+B
    Solution is: X(t) = Ainv * (matExp(A * t) * B-Identity * B)
    Or: Solve(A * X(t) = matExp(A * t) * B-Identity * B)
    =#
    # ts now can be an array as of OrdinaryDiffEq 3.21.0

	st     = vecPars[1]
    ka     = vecPars[2]
    kd     = vecPars[3]
    kb     = vecPars[4]
    sz     = vecPars[5]
    lc     = vecPars[6]
    sc     = vecPars[7]

    intMaxTail = 249 #size of matrix to initialize

    #Generate matrix A
    vecDiag = [-kb / (1 + exp( -(x - lc) / sc)) - kd for x in 0:intMaxTail]
    vecSubDiag = zeros(intMaxTail) + kd
    arrA = Bidiagonal(vecDiag, vecSubDiag, false)

    #Construct the source matrix
    vecB = Float64[ka * pdf(NegativeBinomial(sz, sz / (sz + st)), x) for 
        x = intMaxTail:-1:0]

    #### This is a basic banded solver that works pretty well
    T = ts[length(ts)]
    fX(vecX, p, t) = arrA * vecX + vecB
    vecX0 = zeros(intMaxTail + 1)
    vecTspan = (0.0, T)
    objProb = ODEProblem(fX, vecX0, vecTspan)
    objX = solve(objProb, CVODE_BDF(linear_solver =:Band,
        jac_upper = 1,jac_lower = 1),dense=true,abtol = 1e-16)

    #### This is an altnerative approach using the dense matrix exponential. 
    # objX(t) = arrA\(expm(Matrix(arrA*t))*vecB - vecB)

    return(objX,arrA)
    #=
    For the steady state solution:
    arrSS = arrA\-vecB
    =#
end

#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()
function dF(vecPars,intMaxTail)
    #=
    Compute the derivatives of each parameter and return those to the adjoint
     system as a matrix when called.
    =#
    
    st     = vecPars[1] #starting tail length
    ka     = vecPars[2] #production rate
    kd     = vecPars[3] #deadenylation rate
    kb     = vecPars[4] #decapping rate
    sz     = vecPars[5] #spread for the production rate
    lc     = vecPars[6] #decapping rate position
    sc     = vecPars[7] #decapping rate spread

    #derivate with respect to kb
    arrAdb = Diagonal([-(1 + exp.(-(x - lc) / sc))^(-1.0) for x=0:(intMaxTail)])
    
    #derivate with respect to kd
    vecDiagdkd = -ones(intMaxTail + 1)
    vecSubDiagdkd = ones(intMaxTail)
    arrAdkd = Bidiagonal(vecDiagdkd, vecSubDiagdkd, false)
    
    #derivative with respect to lc
    arrAdlc = Diagonal([kb * exp.((lc - x) / sc) / 
        ((1 + exp.((lc - x) / sc))^2 * sc) for x=0:(intMaxTail)])
    
    #derivative of sc
    arrAdsc = Diagonal([-kb * (exp.((lc - x) / sc) * (lc - x)) / ((1 + 
        exp.((lc - x) / sc))^2 * sc^2) for x=0:(intMaxTail)])
    
    #derivative of source with respect to ka
    vecBdka = Float64[pdf(NegativeBinomial(sz, sz / (sz + st)), xEval) for
        xEval = intMaxTail:-1:0]

    #derivative of source with respect to st
    vecBdst = Float64[(st - xEval) / (st * (st + sz)) * (-sz) * ka * pdf(
        NegativeBinomial(sz, sz / (sz + st)), xEval)
            for xEval=intMaxTail:-1:0] 
    #derivative of source with respect to sz
    vecBdsz = zeros(intMaxTail + 1)
    for xEval in intMaxTail:-1:0
        f64Term1                        = st - xEval + (st + sz) * 
            log(sz / (st + sz))
        f64Term2                        = -(st + sz) * digamma(sz) + (st + sz) * 
            digamma(sz + xEval)
        f64Bdsz                         = 1 / (st + sz) * ka * pdf(
            NegativeBinomial(sz, sz / (sz + st)), xEval) * (f64Term1 + f64Term2)
        vecBdsz[intMaxTail - xEval + 1] = f64Bdsz
    end

    return(vecBdst, vecBdka, arrAdkd, arrAdb, vecBdsz, arrAdlc, arrAdsc)
end

#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()
function Adjoint(vecPars, objX, arrA, arrResidualsdX, ts)
    #=Run the adjoint system, with all of the integration 
    for each time point stepwise.
    =# 
    intTsNo = length(ts)
    T = ts[intTsNo]
    intMaxTail = 249
    vecTspan = (0.0, T)
    fLambda(vecX, p, t) = arrA' * vecX

    #### These next two lines solve the ODE for lambda, the adjoint variable
    objProb = ODEProblem(fLambda, ones(intMaxTail + 1), vecTspan)
    objLambda = solve(objProb, CVODE_BDF(linear_solver =:Band, jac_upper = 1,
        jac_lower = 1), dense = true, abstol = 1e-16)

    #### Generates the derivatives of each of the pars
    (vecBdst, vecBdka, arrAdkd, arrAdb, vecBdsz, arrAdlc, arrAdsc) = 
        dF(vecPars,intMaxTail)

    #### Allows the integration to proceed such that each edge is multiplied by
    #### correct data
    function initialVal(t)
        if t < 40.0
            return(arrResidualsdX[:, 1])
        elseif 40.0 < t <= 60.0
            return(arrResidualsdX[:, 2])
        elseif 60.0 < t <= 120.0
            return(arrResidualsdX[:, 3])
        elseif 120.0 < t <= 240.0
            return(arrResidualsdX[:, 4])
        elseif 240.0 < t <= 480.0
            return(arrResidualsdX[:, 5])
        end
    end

    #### Multiply the correct parameters by their abundances, 
    #### as some derivates depend on the system
    objFp(t) = -[vecBdst vecBdka arrAdkd*objX(t) arrAdb*objX(t
        ) vecBdsz arrAdlc*objX(t) arrAdsc*objX(t)]


    function objIntegrand(t, v)
        v[:] = (objLambda(T - t).* initialVal(T - t))' * objFp(T - t)
    end

    #### what kind of tolerance do I need in this integration?
    #### Vectorized integration, stepwise between edges
    vecGradient   = hquadrature(length(vecPars), (t, v) -> 
        objIntegrand(t, v), T, 240.0, abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars), (t, v) -> 
        objIntegrand(t, v), 240.0, 120.0, abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars), (t, v) -> 
        objIntegrand(t, v), 120.0, 60.0, abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars), (t, v) -> 
        objIntegrand(t, v), 60.0, 40.0, abstol=1e-16)[1] 
    vecGradient  += hquadrature(length(vecPars), (t, v) -> 
        objIntegrand(t, v), 40.0, 0.0, abstol=1e-16)[1] 

    return(vecGradient)
end

#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()
function loss!(vecPars, grad, arrData, ts)
    #= This is the loss function for the optimization. 
    It computes updates the gradient and returns the residual.
    =#

    global tick #### Global variables to monitor the loss
    global gradRatio

    #### The residuals are variance weighted
    #### Excluding steady state for now.
    vecWeigh = [2.58757e-16, 3.053586e-16, 2.92142e-15, 1.16722e-14, 
        2.772866e-14].^-1

    #### Run simulation and compute the residual
    (objX, arrA) = Simulation(exp.(vecPars), ts)
    arrPred = objX(ts)[:, :]
    floatResidual = sum((arrPred - arrData).^2 .* vecWeigh')
    tick = vcat(tick, [tick[size(tick)[1], 1] + 1 floatResidual])
    arrResidualsdX = (arrPred - arrData).*vecWeigh'.*2

    #### Compute the gradient. 
    if length(grad)>0
        println("Tick: $(tick[size(tick)[1],1]), Residual: $floatResidual")
        gradAdjoint = Adjoint(exp.(vecPars), objX, arrA, arrResidualsdX, 
            ts).* exp.(vecPars)
        grad[:] = gradAdjoint
    end
    return (floatResidual,grad)
end
#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()#@$%^&&*()
function optim(vecPars, arrData, ts)
    #### Perform the optimization
    objOpt = Opt(:LD_LBFGS,7) #Optimization object, with parameters tuned below. 
    min_objective!(objOpt, (vecPars, grad) -> 
        loss!(vecPars, grad, arrData, ts)[1])
    lower_bounds!(objOpt, 
        log.([1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]))
    upper_bounds!(objOpt, 
        log.([250, 1e2, 1e1, 1e1, 100.0, 500, 500]))
    xtol_rel!(objOpt, 1e-8)
    ftol_abs!(objOpt, 1e-3) #### Function tolerance needs to be set as 
    #### the final value fluctuates

    xtol_rel!(objOpt, 1e-8)
    maxeval!(objOpt, 200)
    (minf, minx, ret) = optimize(objOpt, vecPars) #### Optimize
end

function io()
    #### Data input
	#### CSV.read can't read in the data, need another function
    
    ### This line reads in the data, need to include the path to the file. 
    tupData = readdlm("/lab/solexa_bartel/teisen/RNAseq/Scripts/TailKineticsModels/TEST_DATA.txt",
        '\t',Float64,'\n';header=true)

    #### Only first gene, don't flatten, exclude the steady state	
    arrData = tupData[1][1:5, 1:250]'
	strAcc = tupData[2]
    println("Gene being optimized:\t", strAcc[1])
    return arrData
end

function main()
    ts = [40.0, 60.0, 120.0, 240.0, 480.0] #### Time points to optimize
    vecPars = log.([140, 1e-7, .1, .1,10.0, 264.0, 6.0]) #### Test input pars
    arrData = io() #### Read in the data

    #### Finite differences gradient
    arrFiniteDiffGrad = Calculus.gradient(x -> 
        loss!(x, [], arrData, ts)[1], vecPars)

    #### Analytical gradient
    arrAnalyticalGrad = loss!(vecPars,[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        arrData, ts)[2]
    arrCompare = [arrFiniteDiffGrad arrAnalyticalGrad vecPars arrFiniteDiffGrad./
     (arrAnalyticalGrad)]
    optim(vecPars, arrData, ts) #log of pars
end

@time main()
