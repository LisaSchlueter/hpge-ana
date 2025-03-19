function rms(x::Vector{T}) where {T <: Real} 
    sqrt(mean(x.^2))
end

"""
    filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, τ_pz::Quantity{T}, ft::Quantity{T}, grid_rt::StepRangeLen{Quantity{<:Float64}}, bl_window::ClosedInterval{<:Unitful.Time{T}}, flt_length_cusp::Quantity{T}; cusp_scale::T = 3750.0, τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs") where T<:Real
    filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, τ_pz::Quantity{T}, ft::Quantity{T}, grid_rt::StepRangeLen{Quantity{<:Float64}}, bl_window::ClosedInterval{<:Unitful.Time{T}}; kwargs...) where T<:Real
    filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, τ_pz::Quantity{T}, dsp_config::DSPConfig; kwargs...) where T<:Real

DSP filter optimization to find best rise-time (for a given flat-top time) to minimize ENC noise. This is an alternative way to calculate the enc noise compared to `dsp_trap_rt_optimization`.
Strategy: 
- Shift waveforms to have a common baseline, and deconvolute them with the pole-zero correction
- Filter waveforms with given rise-time and flat-top time
- Build histogram out of all samples in baseline from all waveforms. Remove bins at the beginning and end of the waveform to avoid edge effects.
- Calculate the RMS of the baseline noise --> ENC noise

Inputs:
- `filter_type::Symbol`: filter type (:trap or :cusp)
- `wvfs::ArrayOfRDWaveforms`: raw waveforms to be filtered
- `τ_pz::Quantity{T}`: pole-zero decay time 
- `ft::Quantity{T}`: fixed flat-top time for optimization
- `grid_rt::StepRangeLen{Quantity{<:Float64}}`: grid of rise-times to scan
- `bl_window::ClosedInterval{<:Unitful.Time{T}}`: baseline window to calculate noise
- `flt_length_cusp::Quantity{T}`: flat-top length; only relevant for cusp filter
- `cusp_scale::T`: scale factor for cusp filter; only relevant for cusp filter
- `τ_cusp::Quantity{<:AbstractFloat}`: cusp decay time; only relevant for cusp filter
Alternatively, a `dsp_config` can be given as input instead of `ft`, `grid_rt`, `bl_window`, `flt_length_cusp`.
- `dsp_config::DSPConfig`: DSP configuration object
"""
function filteropt_rt_optimization_blnoise(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig, τ_pz::Quantity{T}; 
             ft::Quantity{T}= 0.0u"µs", τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs", τ_zac::Quantity{<:AbstractFloat} = 10000000.0u"µs" ) where T<:Real

    # gather config parameters 
    grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
    if ft == 0.0u"µs"
        (_, ft) = get_fltpars(PropDict(), :trap, dsp_config)
    end
    bl_window = dsp_config.bl_window

    # cusp-specific filter parameters 
    flt_length_cusp = getproperty(dsp_config, Symbol("flt_length_cusp"))
    cusp_scale = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # zac-specific filter parameters 
    flt_length_zac = getproperty(dsp_config, Symbol("flt_length_zac"))
    zac_scale = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # waveforms: shift baseline and deconvolute (pole-zero)
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs_shift = shift_waveform.(wvfs, -bl_stats.mean)
    deconv_flt = InvCRFilter(τ_pz)
    wvfs_pz = deconv_flt.(wvfs_shift)

    # STEP 1: rise-time optimization 
    # calculate baseline rms after filtering 
    # filtered waveforms are already truncated to valid windows. no additional truncation needed.
    noise = zeros(length(grid_rt))
    for (i, rt) in enumerate(collect(grid_rt))
        if filter_type == :trap
            wvfs_flt =  TrapezoidalChargeFilter(rt, ft).(wvfs_pz)   # filtered waveforms 
        elseif filter_type == :cusp
            wvfs_flt = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale).(wvfs_pz) 
        elseif filter_type == :zac
            wvfs_flt = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale).(wvfs_pz)
        end
        valid_bins = findall(leftendpoint(bl_window) .<=  wvfs_flt[1].time .<=rightendpoint(bl_window)) # bins within baseline of filtered waveform 
        bl_trap = filter.(isfinite, map(x-> x.signal[valid_bins], wvfs_flt)) # baselines - cut bins in beginning and end 
        noise[i] = rms(vcat(bl_trap...))
    end

    # find optimal shaping time: 
    if  all(.!isfinite.(noise))
        @error "Error: no finite noise values found for filter type $filter_type - try adjusting grid_rt."
    end
    result, report = let enc = noise[findall(isfinite.(noise))], rt = ustrip.(collect(grid_rt))[findall(isfinite.(noise))]
        if length(enc) >= 4
            f_interp = BSplineKit.interpolate(rt, enc, BSplineOrder(4))
        else
            f_interp = LinearInterpolation(rt, enc)
        end
        result = optimize(f_interp, minimum(rt), maximum(rt)) 
        rt_opt = Optim.minimizer(result)*u"µs" # optimial rise time for trap filter
        min_enc = Optim.minimum(result)
        (rt = rt_opt, min_enc = min_enc), (rt = rt_opt, min_enc = min_enc, enc_grid_rt = grid_rt, enc = noise, f_interp = f_interp)
    end
    return result, report 
end

