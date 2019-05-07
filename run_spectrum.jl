include("complex_roots.jl");
include("spectrum.jl");
R = 2.5; rectangle = [-R,-0.001,-2*R,2*R];
maxevals = 10^4;

β = 2.5
U(x)=1/β*((1+x^2)^(β/2) - 1);
dUdx(x)=x*(1 + x^2)^(β/2-1);

(ψ,dψdz)=construct_ψ(U,dUdx);

#using QZ40
integrand_plus(z::Complex) = - dψdz(z)/(1-ψ(z));
integrand_minus(z::Complex) = dψdz(z)/(1+ψ(z));
roots_plus = qz40(integrand_plus,rectangle,true,maxevals)
roots_minus = qz40(integrand_minus,rectangle,true,maxevals)

using Plots;
using LaTeXStrings;

plot([0], [0], seriestype=:scatter,c=:white,markershape= :diamond,label="")
plot!(real(roots_minus),imag(roots_minus),seriestype=:scatter,c=:black,label=L"L^-")
plot!(real(roots_plus),imag(roots_plus),seriestype=:scatter,c=:white,markershape= :diamond,label=L"L^+")
epsilon = 0.1
perturbation_minus = zeros(Complex, size(roots_minus))
for i in 1:length(roots_minus)
    perturbation_minus[i] = epsilon * refresh_perturbation_symmetric(U,dUdx, ψ, roots_minus[i], -1)
end
perturbation_plus = zeros(Complex, size(roots_plus))
for i in 1:length(roots_plus)
    perturbation_plus[i] = epsilon * refresh_perturbation_symmetric(U,dUdx, ψ, roots_plus[i], +1)
end

quiver!(real(roots_minus),imag(roots_minus), quiver=(real(perturbation_minus),imag(perturbation_minus)),c=:black)

quiver!(real(roots_plus),imag(roots_plus), quiver=(real(perturbation_plus),imag(perturbation_plus)),c=:black,linestyle=:dot)
savefig("perturbation-beta-$(β).pdf")
