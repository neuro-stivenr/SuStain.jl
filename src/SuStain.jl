#!/usr/bin/env julia

module SuStain

using Random, Distributions, QuadGK
using CSV, DataFrames
using Base.Threads
using ProgressMeter
import YAML


# the piecewise linear function evaluation
function g(
	t::Float64,
	zvec::Vector{Float64},
	zmax::Float64,
	tvec::Vector{Float64}
)::Float64
	idx = 1 + sum(t .> tvec)
	if idx == 1
		z_event = zvec[begin]
		t_event = tvec[begin]
		return (z_event / t_event) * t
	end
	if idx == (length(zvec) + 1)
		z_event = zvec[end]
		t_event = tvec[end]
		Δz = (zmax - z_event)
		Δt = (1 - t_event)
		return z_event + (Δz / Δt) * (t - t_event)
	end
	z_event = zvec[idx]
	z_event_prev = zvec[idx - 1]
	t_event = tvec[idx]
	t_event_prev = tvec[idx - 1]
	Δz = (z_event - z_event_prev)
	Δt = (t_event - t_event_prev)
	return z_event_prev + (Δz / Δt) * (t - t_event_prev)
end

function P(
	x::Float64,
	t::Float64,
	zvec::Vector{Float64},
	zmax::Float64,
	tvec::Vector{Float64};
	σ::Float64=1.0
)
	gt = g(t, zvec, zmax, tvec)
	distr = Normal(gt, σ)
	return logpdf(distr, x)
end


function stage_to_t(k::Int64, N::Int64)
	return k / (N + 1)
end


function gen_rand_tvec(N::Int64)
	tvec = Vector{Float64}(undef, N)
	tvec[1] = rand(Uniform(0, 1))
	for n in 2:N
		tvec[n] = rand(Uniform(tvec[n - 1], 1))
	end
	return tvec
end


function model_likelihood(
	df_reg::DataFrame,
	N_stages::Int64,
	z_dict::Dict{String, Vector{Float64}},
	z_max_dict::Dict{String, Float64},
	tvec_dict::Dict{String, Vector{Float64}}
)
	j_prod = 0.0
	for j in eachrow(df_reg)
		biomarkers = names(j)
		k_sum = 0.0 # sum of likelihoods across stages for a given participant
		for stage in 0:(N_stages-1)
			t_i = stage_to_t(stage, N_stages)
			t_f = stage_to_t(stage + 1, N_stages)
			k_sum += quadgk( # integrate likelihood for this window and add to the total sum
                # first input is a function, which is to be run on every slice of time within
                # a stage, while the outputs are cumulatively summed together
				t -> 
                begin
					x_prod = 0 # product of biomarker likelihoods
                    # we implement it as a sum of log likelihoods,
                    # in order to avoid issues with arithmetic underflow
					for biomarker in biomarkers
						x_prod += P( # multiplication is adding because log
							j[biomarker], t,
							z_dict[biomarker],
							z_max_dict[biomarker],
							tvec_dict[biomarker]
						) # log likelihood
					end
                    x_prod = exp(x_prod) # bring it back into likelihood via exponent
					return x_prod
				end,
                # start and end of time window that we are integrating over
				t_i, t_f
			) |> first # output is a single value, which is the integral of likelihood over window
            # second output is not useful
		end
		# we take a log of the sum of stage likelihoods to again safely calculate the cumulative product
		j_prod += log(k_sum)
	end
	return j_prod # we keep the product of participant likelihoods logged, for future convenience
end


function gen_z_dict(df::DataFrame, biomarkers::Vector{String}; value=nothing)
	map(biomarkers) do biomarker
		if !isnothing(value)
			return biomarker => value
		end
		vec_biomarker = df[:,biomarker]
		vec_biomarker_pos = vec_biomarker[vec_biomarker .> 0]
		biomarker => [
			quantile(vec_biomarker_pos, 0.33),
			quantile(vec_biomarker_pos, 0.66),
			quantile(vec_biomarker_pos, 0.99)
		]
	end |> Dict
end

function gen_zmax_dict(df::DataFrame, biomarkers::Vector{String}; value=nothing)
	map(biomarkers) do biomarker
		if !isnothing(value)
			return biomarker => value
		end
		vec_biomarker = df[:,biomarker]
		biomarker => maximum(vec_biomarker)
	end |> Dict
end

function run_models(n_iter, regions, df_reg, N_stages, dict_z, dict_zmax)
	vec_ml = Vector{Float64}(undef, n_iter)
	vec_tvec = Vector{Dict{String, Vector{Float64}}}(undef, n_iter)
	progbar = Progress(n_iter)
	for i in 1:n_iter
		dict_tvec_rand = Dict([
			region => gen_rand_tvec(3)
			for region in regions
		]...)
		lik = model_likelihood(
			df_reg,
			N_stages,
			dict_z,
			dict_zmax,
			dict_tvec_rand
		)
		vec_ml[i] = lik
		vec_tvec[i] = dict_tvec_rand
		next!(progbar)
	end
	finish!(progbar)
	return vec_ml, vec_tvec
end

function run_parallel_models(n_iter, regions, df_reg, N_stages, dict_z, dict_zmax)
	vec_ml = Vector{Float64}(undef, n_iter)
	vec_tvec = Vector{Dict{String, Vector{Float64}}}(undef, n_iter)
	progbar = Progress(n_iter)
	Threads.@threads for i in 1:n_iter
		dict_tvec_rand = Dict([
			region => gen_rand_tvec(3)
			for region in regions
		]...)
		lik = model_likelihood(
			df_reg,
			N_stages,
			dict_z,
			dict_zmax,
			dict_tvec_rand
		)
		vec_ml[i] = lik
		vec_tvec[i] = dict_tvec_rand
		next!(progbar)
	end
	finish!(progbar)
	return vec_ml, vec_tvec
end

function plot_progression(
	dict_opt::Dict{String, Vector{Float64}},
	dict_z::Dict{String, Vector{Float64}},
	dict_zmax::Dict{String, Float64},
	regions::Vector{String}
)
	p = plot(
		xlabel="Time",
		ylabel="Z-score",
		xlims=(0, 1),
		ylims=(0, maximum(values(dict_zmax)))
	)
	for region in regions
		vec_t = [0.0, dict_opt[region]..., 1.0]
		vec_z = [0.0, dict_z[region]..., dict_zmax[region]]
		plot!(vec_t, vec_z, labels="", color="black")
		scatter!(
			vec_t[begin+1:end-1],
			vec_z[begin+1:end-1],
			labels=region,
			legend=:topleft
		)
	end
	return p
end

function coord_to_cutoffs(coord::Tuple{Float64, Float64, Float64})
    value_1 = coord[1]
    leftover_1 = 1 - value_1
    value_2 = value_1 + leftover_1 * coord[2]
    leftover_2 = 1 - value_2
    value_3 = value_2 + leftover_2 * coord[3]
    (value_1, value_2, value_3)
end

function main(datapath::String, variables::String, N_stages::Int64, N_runs::Int64, outpath::String)
    @info "Reading in data... Path: " datapath
    df_data = CSV.read(datapath, DataFrame);
    regions = String.(split(variables, ","))
    @info "Regions: " regions
    df_reg = df_data[:,regions];
    @info "Nulling out negative data values."
    df_reg = df_reg .* (df_reg .>= 0)
    @info "Generating z-score thresholds based on qunatiles."
    dict_z = SuStain.gen_z_dict(df_reg, regions)
    @info "Generating z-score maxima from the data."
    dict_zmax = SuStain.gen_zmax_dict(df_reg, regions)
    @info "Running $(N_runs) SuStain models."
    vec_loglik, vec_tvec = SuStain.run_parallel_models(N_runs, regions, df_reg, N_stages, dict_z, dict_zmax);
    @info "Identifying maximum likelihood model."
    ml_idx = argmax(vec_loglik)
    ml = vec_loglik[ml_idx]
    @info "Maximal log-likelihood: " ml
    opt_tvec = vec_tvec[ml_idx]
    output = Dict(
        "z_thresh" => dict_z,
        "z_max" => dict_zmax,
        "N_stages" => N_stages,
        "model_ll" => vec_loglik,
        "model_params" => vec_tvec,
        "opt_model_params" => opt_tvec,
        "opt_model_ll" => ml,
        "regions" => regions
    )
    @info "Writing output to: " outpath
    YAML.write_file(outpath, output)
end

if !isinteractive()
    main(ARGS[1], ARGS[2], parse(Int64, ARGS[3]), parse(Int64, ARGS[4]), ARGS[5])
end


end