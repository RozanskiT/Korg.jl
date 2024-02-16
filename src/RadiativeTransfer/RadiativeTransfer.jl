module RadiativeTransfer

# for generate_mu_grid
using FastGaussQuadrature: gausslegendre
"""
    generate_mu_grid(n_points)

Used by both radiative transfer schemes to compute quadature over μ. Returns `(μ_grid, μ_weights)`.
"""
function generate_mu_grid(n_points::Int)
    μ_grid, μ_weights = gausslegendre(n_points)
    μ_grid = @. μ_grid/2 + 0.5
    μ_weights ./= 2
    μ_grid, μ_weights
end

"""
    generate_mu_grid(μ_grid::Vector{Float64})

Generate weights for trapezoidal integration for a non-uniform grid of mu values in domain -1 to 1.

# Arguments
- `μ_grid::Vector{Float64}`: A vector of mu values, not necessarily equally spaced.

# Returns
- `μ_weights::Vector{Float64}`: A vector of weights for trapezoidal integration.
"""
function generate_mu_grid(μ_grid::AbstractVector{<:AbstractFloat})
    n_points = length(μ_grid)
    μ_weights = Vector{Float64}(undef, n_points)
    
    # First weight
    μ_weights[1] = (μ_grid[2] - μ_grid[1]) / 2 + (μ_grid[1] - 0.0)
    
    # Weights for points 2 through n-1
    for i in 2:(n_points - 1)
        left_distance = μ_grid[i] - μ_grid[i - 1]
        right_distance = μ_grid[i + 1] - μ_grid[i]
        μ_weights[i] = (left_distance + right_distance) / 2
    end
    
    # Last weight
    μ_weights[n_points] = (1.0 - μ_grid[n_points]) + (μ_grid[n_points] - μ_grid[n_points - 1]) / 2
    
    μ_grid, μ_weights
end




include("BezierTransfer.jl")
include("MoogStyleTransfer.jl")

end #module