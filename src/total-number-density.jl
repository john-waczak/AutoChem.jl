const kb = 1.380649e−23 # J/K

function M(T,P)
    press_pa = 100.0 * P  # 100 Pa / mbar
    # 1 Pa = 1 J / m3. We want cm^3 so convert:
    press_final = press_pa * 1.0e-6 # there are (10^2)^3 = 10^6 cm³/m³
    Mout = press_final/(kb*T)  # we now have a stand in for number density in molecules/cm³jk
    return Mout
end

# O2 and N2 based on typical relative abundance
O2(T,P) = 0.2095 * M(T,P)
N2(T,P) = 0.7809 * M(T,P)
Ar(T,P) = 0.0094 * M(T,P)

