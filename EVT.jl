using Pkg
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("Plots")

using Distributions
using Random
using Plots

#### Create function for r.v.
strDist = "Pareto";
if strDist == "Normal"
    μ = 0;
    σ = 1
    distr = Normal( μ, σ );
    μDist = μ;
    σDist = σ
elseif strDist == "Pareto"
    ξ = 0.5; # ξ > 1 ⇒ mean non-existent ξ > 2 ⇒  Variance non-existent
    σ = 1; 
    μ = σ / ξ; # μ = σ / ξ ⇒ GPD = Pareto x_m = σ / ξ, α = 1 / ξ
    x_m = σ / ξ;
    distr = GeneralizedPareto(μ,  σ,ξ,)
    μDist = (1 / ξ) * x_m / ((1 / ξ) - 1);
    σDist = ( (σ / ξ)^2 / ξ ) / (( 1/ξ - 1)^2 * (1/ξ - 2))
end

X = rand(distr, 1000)

# Histogram and Cummean-Plot
histogram(X, title = "Histogram " * strDist * " r.v.",legend = false )
cumMeanPlot = plot( cumsum(X) ./ (1:length(X)),
    title = "Cummean of " * strDist * " r.v.", 
    label = "cummean of sample" )
plot!( μDist * ones(length(X),1), label = "theoretical mean",
    ylims=( μDist - σDist, μDist + σDist))
xlabel!("Number of samples")
ylabel!("(Cum-)Mean")
# Mean-Residual-Life-Plot
mrl = zeros(length(X));
X = sort(X)
for k in eachindex(X)
    #println(k) 
    mrl[k] = sum( (X .- X[k]) .* (X .> X[k]) ) ./ sum(X .> X[k] )
end
plot(X , mrl, seriestype = :scatter) # 
