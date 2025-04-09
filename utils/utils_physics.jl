

"""
Ge_Egap(T::Real)
bandgap energy (eV) of germanium as a function of temperature (K)
"""
function Ge_Egap(T::Real)
    return 0.744 - (4.774*1e-4 * T^2) / (T + 235)
end
# Ge_Egap(77)

"""
Ge_Energy_per_eholePair(T::Real)
energy (eV) required to create an electron-hole pair in germanium as a function of temperature (K)
"""
function Ge_Energy_per_eholePair(T::Real)
    return 2.2 * Ge_Egap(T) + 1.99 * Ge_Egap(T)^(3/2) * exp(4.75 * Ge_Egap(T) / T)
end
# Ge_Energy_per_eholePair(90)

"""
Ge_NumberChargeCarrier(E::Real, T::Real)
number of charge carriers created by an energy E (eV) in germanium at temperature T (K)
"""
function Ge_NumberChargeCarrier(E::Real, T::Real)
return E / Ge_Energy_per_eholePair(T)
end 
# Ge_NumberChargeCarrier(1e6, 77)/1e5


electron_charge = 1.60217662e-19 # C


# utility functions for DAQ and pulser (Voltage - Charge - Energy conversion)
function DAQ_ADC_to_V(ADC::Real, dynamicrange_V::Real, bits::Int)
    return (dynamicrange_V / 2^bits) * ADC
end

function pulser_ADC_to_electrons(ADC::Real, C_pulser::Real; bits::Int = 14,  dynamicrange_V::Real = 2.0, gain::Real = 1.0)
    V = DAQ_ADC_to_V(ADC, dynamicrange_V, bits)
    return  (V/gain * C_pulser) / electron_charge  # charge in electrons
end

function pulser_ADC_to_keV(ADC::Real, C_pulser::Real; bits::Int = 14,  dynamicrange_V::Real = 2.0, gain::Real = 1.0)
    return  pulser_ADC_to_electrons(ADC, C_pulser; bits = bits, dynamicrange_V = dynamicrange_V, gain = gain) * Ge_Energy_per_eholePair(90) / 1e3
end

function V_to_electrons(V::Real, C_pulser::Real; gain::Real = 1.0)
    return (V/gain * C_pulser) / electron_charge  # charge in electrons
end
