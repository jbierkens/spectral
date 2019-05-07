# implement QZ-40 algorithm
# Dellnitz, M., Schütze, O., & Zheng, Q. (2002). Locating all the zeros of an analytic function in one complex variable. Journal of Computational and Applied Mathematics, 138(2), 325–333. https://doi.org/10.1016/S0377-0427(01)00371-5
# currently assumes all roots have multiplicity one

using QuadGK

function winding_number(dlogfdz::Function, rectangle::Array{Float64,1}, maxevals::Int = 10^4)::Int
    (x_min, x_max, y_min, y_max) = rectangle
    (bl, br, tr, tl) = (x_min + y_min*im, x_max + y_min*im, x_max + y_max*im, x_min + y_max*im)
    round(Int,real(1 ./(2*pi*im)*quadgk(dlogfdz,bl,br,tr,tl,bl,rtol=0, atol=0.45*2*pi,maxevals=maxevals)[1]))
end

function qz40(dlogfdz::Function, rectangle::Array{Float64,1}, verbose::Bool=false, maxevals::Int=10^4)

    rectangles = Array{Float64,1}[rectangle]
    zeroes = Complex[]
    wn = -1
    qz40_recursive(dlogfdz,rectangles,zeroes,wn,verbose,maxevals)
end

function locate_root(dlogfdz::Function, rectangle::Array{Float64,1}, maxevals::Int = 10^4)
    (x_min, x_max, y_min, y_max) = rectangle
    (bl, br, tr, tl) = (x_min + y_min*im, x_max + y_min*im, x_max + y_max*im, x_min + y_max*im)
    integrand(z::Complex) = z*dlogfdz(z)
    location = 1 ./(2*pi*im)*quadgk(integrand,bl,br,tr,tl,bl,maxevals=maxevals)[1]
    if (real(location) < x_min || real(location) > x_max || imag(location) < y_min || imag(location) > y_max)
        error("Root location $(location) outside of rectangle (x_min, x_max, y_min, y_max) $(x_min),$(x_max),$(y_min),$(y_max).")
    end
    return location
end

function random_subdivide(rectangle::Array{Float64,1})

    (x_min, x_max, y_min, y_max) = rectangle
    x_mid = x_min + rand()*(x_max - x_min)
    y_mid = y_min + rand()*(y_max - y_min)
    # rectangles = Array{Array{Float64,1},1}[]
    rectangles = Array{Array{Float64,1},1}(undef,0)
    push!(rectangles,[x_min, x_mid, y_min, y_mid]);
    push!(rectangles,[x_mid, x_max, y_min, y_mid]);
    push!(rectangles,[x_mid, x_max, y_mid, y_max]);
    push!(rectangles,[x_min, x_mid, y_mid, y_max]);
    return rectangles
end

function qz40_recursive(dlogfdz::Function, rectangles::Array{Array{Float64,1},1},zeroes::Vector{Complex},nroots::Int,verbose::Bool,maxevals::Int)
    # the rows of rectangles are the delimiters (x_min, x_max, y_min, y_max) of the rectangles
    new_rectangles = Array{Float64,1}[]
    n_rectangles = size(rectangles)[1]
    if (verbose)
        println("Number of rectangles: $(n_rectangles)")
    end
    nroots_new = length(zeroes)
    for i in 1:n_rectangles
        rect = rectangles[i]
        wn = winding_number(dlogfdz,rect,maxevals)
        nroots_new += wn
        if (wn == 1)
            push!(zeroes,locate_root(dlogfdz,rect,maxevals))
        elseif (wn > 1)
            append!(new_rectangles,random_subdivide(rect))
        elseif (wn < 0)
            error("Negative winding numbers")
        end
    end
    if (nroots > 0 && nroots_new != nroots)
        error("Expected number of roots: $(nroots), new number of roots: $(nroots_new).")
    end
    if (nroots == -1)
        println("Number of roots: $(nroots_new)")
    end
    if size(new_rectangles)[1] == 0
        return zeroes
    else
        return qz40_recursive(dlogfdz, new_rectangles, zeroes, nroots_new, verbose,maxevals)
    end
end

function newton(f, dfdz, z_init, n_iter)

    z_current = z_init
    for i in 1:n_iter
        z_current = z_current - f(z_current)/dfdz(z_current)
    end
    z_current
end
