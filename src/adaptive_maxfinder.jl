# Here a build an algorithm for maximum/minimum finder on the go.
# So instead of finding the maximum of the CPR for a constant mesh of data points. 
# we adaptively guide the Ic calculation so the maximum up to the desired tolerance is found.
# Assumptions: There was no intent to make this algorithm general, it is a quick fix that solves 
# the problem for continous and periodic sinusoidal-triangular functions on a -pi:pi interval.


# The algorithm does the following. 
#   Compute y values for 5 equally spaced x values (pi/2 separation)
#   Fit using Fourier series and find the correct truncation to avoid overfitting (cross validation)
#   Find numerically the minima (actually it doesn't matter) of this fit
#   Compute the y value at this global minima
#   Add the new data to x and y and fit again until tolerance is met or hard threshold is met
#   returns xmin and ymin

"it computes the minima value of some computationally demanding `generating_function` in a fixed 
range (-π, π) using an adaptive algorithm based on numerical optimization of  K-cross validated fits
which guide adaptively the minima value search. `tolerance` sets when the algorithm should stop searching maximums. `max_iter` sets an upper limit for 
maximal search. It is important to use an initial discrete xmesh that uses the 2pi periodicity to guide the fit"
function adaptive_max_finder(generating_func::Function, model; tolerance = 1e-4, max_iter = 10, kw...) 
    println("adaptive algorithm")
    # xinit = collect(-1:0.1:1)
    # xinit = 2 .* rand(5) .- 1ß
    xinit = collect(-1:0.25:1) + (2 .* rand(9) .- 1)./10
    yinit, fitparams = generate_and_fit(xinit, model, generating_func; kw...)
    cpr_plot(xinit, yinit, model, fitparams)
    xmax = minima_search(model, fitparams)
    ymax = generating_func([xmax]; kw...)
    xmax0 = 10.0
    count = 0
    x = xinit
    y = yinit
    while abs(xmax0-xmax) > tolerance 
        # println("error:", abs(xmax0-xmax))
        err_old = abs(xmax0-xmax)
        xmax0 = xmax
        x, y, ymax, fitparams = generate_and_fit(xmax, x, y, model, generating_func; kw...)
        # cpr_plot(x, y, model, fitparams)
        xmax = minima_search(model, fitparams)
        err_new = abs(xmax0-xmax)
        err_new > err_old ? break : nothing
        count += 1
        if count > 20
            println("convergence was not reached")
            break 
        else nothing end   
    end
    return ymax
end
"generate a new data point using a given `model` and an `x` value"
function generate_and_fit(x::Array, model, generating_func; kw...)
    y = generating_func(x; kw...)
    fitparams = msqerror(x, y, model) 
    return y, fitparams
end

function generate_and_fit(x::Number, xlist, ylist, model, generating_func; kw...)
    y = generating_func([x]; kw...)
    xlist = vcat(xlist, x) ####
    ylist = vcat(ylist, y)
    fitparams = msqerror(xlist, ylist, model)
    return xlist, ylist, y, fitparams
end

"numerical search of the global minima"
function minima_search(model, fitparams)
    opt = optimize(x -> model(x, fitparams), -1, 1)
    return Optim.minimizer(opt)
end

"This function computes the optimal fitting parameters for a given dataset consisting in a 
truncated up to `custom_max_order` Fourier expansion, it performs a K-cross validation 
analysis as well, to prevent overfitting (ML algorithm)"
# Warning: Cross validation isn't working nicely check it later. Currently disabled
function msqerror(xdata, ydata, model; K = 1, custom_max_order = missing, train_tol = 1e-3)
    custom_max_order === missing ? max_order = Int64(floor(length(xdata)/2)) : 
        max_order = custom_max_order
    # data partition controlled by data availability 
    l = length(xdata)
    part = Int64(floor(l/K))
    println("partition size: ", part)
    fitted_params = []
    errorstest = zeros(max_order)
    errorstrain = zeros(max_order)
    perm_ind = shuffle(1:l) #reshufle data
    xshuffled = xdata[perm_ind]    
    yshuffled = ydata[perm_ind]
    l ≥ 2K ? nothing : K = 1   #criteria for cross validation(need at least the double of points)
    for i in 1:K #average over different K partitions
        errors_train = []
        errors_test = []
        xnew = vcat(xshuffled[Int64(part*(i-1)+1):end], xshuffled[1:Int64(part*(i-1))])
        ynew = vcat(yshuffled[Int64(part*(i-1)+1):end], yshuffled[1:Int64(part*(i-1))])
        xtrain = xnew[1:part]
        xtest = xnew[part:end]
        ytrain = ynew[1:part]
        ytest = ynew[part:end]
        # compute the fittings for each model and the train and test data set errors
        for order in 1:max_order 
            p0 = rand(Int64(2*order))
            fit = curve_fit(model, xtrain, ytrain, p0)
            push!(errors_train, sum(x -> x, (ytrain - model(xtrain, fit.param)).^2))
            push!(errors_test,sum(x -> x, (ytest - model(xtest, fit.param)).^2))
        end
        errorstest += errors_test 
        errorstrain += errors_train
    end
    good_training_indices =  findall(x-> x < train_tol, errorstrain)
    while good_training_indices == []
        train_tol *= 2
        @info "tolerance was increased to $(train_tol)"
        good_training_indices = findall(x-> x < train_tol, errorstrain)
    end
    optimised_order = good_training_indices[findmin(ifelse(K > 1, errorstest, errorstrain)[good_training_indices])[2]]
    # now we fit using the optimal order and all data points
    fit = curve_fit(model, xdata, ydata, rand(Int64(2*optimised_order))) 
    println("optimised order: ", optimised_order)
    return fit.param
end


fourier_model(x, p) = sum(i -> 
    sin_cos_model(x, p[Int64((i-1)*2 +1):Int64(i*2)], i-1), 1:Int64(length(p)/2))

# gen_model(x) = sin_model(x, [1.0 1. 2.0]) + 0.1*randn(length(x))
# gen_model2(x) = 0.5*triang(x, [1.0 1.]) + 0*randn(length(x)) + sin_model(x, [1, 1, .3])
@. sin_cos_model(x, p, i) = p[1]*sin(i*x) + p[2]*cos(i*x) 
@. sin_model(x, p) = p[1]*sin(p[2]x+p[3]) 
@. cos_model(x, p) = p[1]*cos(x/2+p[3])
@. triang(x, p) = p[1]*sawtoothwave(x+p[2])


## Usage
# xdata = range(-π, stop=π, length=20)
# ydata = sin_model(xdata, [1.0 1. 2.0]) + 0.1*randn(length(xdata)) 
# how to set the order of the expansion so there is no overfitting
# order_exp = 10 # do cross validation later to improve it.
# p0 = rand(Int64(2*order_exp)) #rand(6)

# fit = curve_fit(mod, xdata, ydata, p0)
# cpr_plot(xdata, ydata, mod, fit.param)

# msqerror(xdata, ydata, fourier_model)

#adaptive_max_finder(gen_model2, fourier_model)