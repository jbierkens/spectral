
function construct_ψ(U::Function, dUdx::Function, maxevals::Int=10^4)

    integrand(x::Real,z::Complex) = exp(-2*z*x-U(x))*dUdx(x)
    dz_integrand(x::Real,z::Complex) = -exp(-2*z*x-U(x))*dUdx(x)*2*x
    ψ(z::Complex) = quadgk(x -> integrand(x,z),0,Inf,maxevals=maxevals)[1]
    dψdz(z::Complex) = quadgk(x -> dz_integrand(x,z),0,Inf,maxevals=maxevals)[1]

    return(ψ,dψdz)

end

function construct_eigenfunction(U::Function, dUdx::Function, ψ::Function, γ::Complex, maxevals::Int = 10^4)

    f_plus(x) = (x < 0.0) ? ψ(γ) * exp(γ * x) : exp(γ * x + U(x)) * quadgk(ξ -> dUdx(ξ) * exp(-2 * γ * ξ - U(ξ)), x, Inf, maxevals=maxevals)[1];
    f_minus(x) = (x < 0.0) ? -ψ(γ) * exp(-γ * x + U(x)) * quadgk(ξ -> dUdx(ξ) * exp(2 * γ * ξ - U(ξ)),-Inf,x, maxevals=maxevals)[1] : exp(-γ * x);

    return (f_plus, f_minus);

end

function refresh_perturbation(U::Function, dUdx::Function, psi_plus::Function, γ::Complex, maxevals::Int = 10^4)

    (f_plus, f_minus) = construct_eigenfunction(U, dUdx, psi_plus, γ, maxevals);

    integrand_numerator = x -> exp(-U(x))* (f_plus(x)^2 + f_minus(x)^2);
    integrand_denominator = x -> exp(-U(x))* 2 * f_plus(x)*f_minus(x)

    numerator = quadgk(integrand_numerator, -30,30, maxevals=maxevals)[1];
    denominator = quadgk(integrand_denominator, -30, 30, maxevals=maxevals)[1];

    return numerator/denominator - 1;

end

function construct_eigenfunction_symmetric(U::Function, dUdx::Function, ψ::Function, γ::Complex, sgn::Int, maxevals::Int = 10^4)
    # construct +/- eigenfunction in the symmetric case, depending on value of sgn (-1/+1)

    f(x) = (x < 0.0) ? exp(γ * x) : sgn * exp(γ * x + U(x)) * quadgk(ξ -> dUdx(ξ) * exp(-2 * γ * ξ - U(ξ)), x, Inf, maxevals=maxevals)[1];

    return f

end

function refresh_perturbation_symmetric(U::Function, dUdx::Function, ψ::Function, γ::Complex, sgn::Int, maxevals::Int = 10^4)

    f = construct_eigenfunction_symmetric(U, dUdx, ψ, γ, sgn, maxevals);
    integrand_numerator = x -> exp(-U(x)) * f(x)^2;
    integrand_denominator = x -> exp(-U(x)) * f(x) * f(-x);
    numerator = quadgk(integrand_numerator,-10, 10, maxevals=maxevals)[1];
    denominator = 2 * quadgk(integrand_denominator, 0, 10, maxevals=maxevals)[1];
    return sgn * numerator/denominator - 1;

end
