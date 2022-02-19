mutable struct Model{T,Tv,Tt}
    n::Int  # Num vars
    m::Int  # Num cons
    x::Tv  # Starting and final solution
    x_L::Tv # Variables Lower Bound
    x_U::Tv # Variables Upper Bound
    g::Tv  # Final constraint values
    g_L::Tv # Constraints Lower Bound
    g_U::Tv # Constraints Upper Bound
    j_str::Tt
    h_str::Tt
    mult_g::Tv # lagrange multipliers on constraints
    mult_x_L::Tv # lagrange multipliers on lower bounds
    mult_x_U::Tv # lagrange multipliers on upper bounds
    obj_val::T  # Final objective
    status::Int  # Final status
    g_order::Vector{Int}
    iter::Int

    # Callbacks
    eval_f::Function
    eval_g::Function
    eval_grad_f::Function
    eval_jac_g::Function
    eval_h::Union{Function,Nothing}

    intermediate  # Can be nothing

    # For MathProgBase
    sense::Symbol
    convex_model

    parameters::Parameters
    statistics::Dict{String,Any}   # collects parameters of all iterations inside the algorithm if StatisticsFlag > 0

    Model(
        n::Int, 
        m::Int, 
        x_L::Tv, 
        x_U::Tv,
        g_L::Tv,
        g_U::Tv,
        j_str::Tt,
        h_str::Tt,
        eval_f::Function,
        eval_g::Function,
        eval_grad_f::Function,
        eval_jac_g::Function,
        eval_h::Union{Function,Nothing},
        parameters::Parameters,
        g_order, convex_model
    ) where {T, Tv<:AbstractArray{T}, Tt<:AbstractArray{Tuple{Int64,Int64}}} = new{T,Tv,Tt}(
        n, m,
        zeros(n), x_L, x_U,
        zeros(m), g_L, g_U,
        j_str, h_str,
        zeros(m), zeros(n), zeros(n),
        0.0,
        -5,
        g_order,
        0,
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h, 
        nothing, :Min, convex_model,
        parameters,
        Dict{String,Any}()
    )
end

function optimize!(model::Model)
    if isnothing(model.parameters.external_optimizer)
    	model.status = -12;
        @error "`external_optimizer` parameter must be set for subproblem solutions."
    else
        if model.parameters.method == "SLP"
            if model.parameters.algorithm == "Line Search"
                slp = SlpLS(model)
            elseif model.parameters.algorithm == "Trust Region"
                slp = SlpTR(model)
            end
            run!(slp)
        else
            @error "The method is not defined"
        end
    end
    return nothing
end

function add_statistic(model::Model, name::String, value)
    if model.parameters.StatisticsFlag == 0
        return
    end
    model.statistics[name] = value
end

function add_statistics(model::Model, name::String, value::T) where T
    if model.parameters.StatisticsFlag == 0
        return
    end
    if !haskey(model.statistics, name)
        model.statistics[name] = Array{T,1}()
    end
    push!(model.statistics[name], value)
end
