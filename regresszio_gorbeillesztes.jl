include("lin_reg_func.jl")
using LsqFit

N = 10;

x = (rand(N) .- 0.5) * 20
zaj = (randn(N) .- 0.5)

fx(x) = 2 .* x .- 3 .+ zaj
y = fx(x)

beta0, beta1 = lin_reg(x, y)
m(t, p) = p[1] .+ p[2] .* t
p0 = [0.5, 0.5]
fit = curve_fit(m, x, y, p0)
display((fit.param - [beta0, beta1]) ./ [beta0, beta1])

## Mátyás példa adatai legkisebb négyzetekre

P = [229, 323, 470, 655, 889, 1128, 1701, 2036, 2502, 3089, 3858, 4636, 5883, 7375, 9172, 10149, 12462, 15113, 17660, 21157]
t = range(1, length(P));


a, b = proba1(t, P)
m1(t, p) = p[1] * exp.(t) .+ p[2]
p0 = [0.5, 0.5]
fitProba1 = curve_fit(m1, t, P, p0)
display((fitProba1.param - [a, b]) ./ [a, b])

logP = log.(P)

A2, b2 = lin_reg(t, logP)

a2 = exp.(A2)

scatter(t, P)
(plot!(t, a2 * exp.(b2 .* t)))

m2(t,p) = p[1] .* exp.(p[2] .* t)
fitProba2 = curve_fit(m2, t, P, p0)

display((fitProba2.param - [a2, b2]) ./ [a2, b2])
(plot!(t, fitProba2.param[1] * exp.(fitProba2.param[2] .* t)))
m3(t,p) = p[1] .* exp.(p[2] .* t) .+ p[3]
fitProba3 = curve_fit(m3, t, P, [0.5,0.5,0.5])
display(plot!(t, fitProba3.param[1] * exp.(fitProba3.param[2] .* t).+fitProba3.param[3]))
