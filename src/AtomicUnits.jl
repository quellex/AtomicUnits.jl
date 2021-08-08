# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
# unit conversion
module AtomicUnits
using PhysicalConstants
export fac_i2e, fac_vnm2e, fac_wlenev, fac_au2eV
export i2e, e2i, vnm2e, e2vnm, i2vnm, vnm2i
export au2as, as2au, wlen2au, au2wlen, wlen2eV, eV2wlen, eV2au, au2eV
export A2i, i2A
"""
fac_i2e = 3.50944506e16
convert factor: intensity [W/cm²] ↔ electric field [a.u.]
1 W/cm² = 3.50944506e16 (1 a.u.)²
"""
const fac_i2e = 3.50944506e16   # intensity [W/cm^2] → amplitude [a.u.]
"""
fac_vnm2e = 514.221
convert factor: electric field [V/nm] ↔ electric field [a.u.]
1 a.u. = 514.221 V/nm
"""
const fac_vnm2e = 514.221# amplitude [V/nm] → amplitude [a.u.]
"""
fac_au2as = 24.1899
convert factor: time [a.u.] ↔ time [atto. sec]
1 a.u. = 24.1899 as
"""
const fac_au2as = 24.1899# time [a.u.] → time[atto sec.]
"""
fac_wleneV = 1239.84190
convert factor: wavelength [nm] ↔ photon energy [eV]
1 nm * 1 eV = 1239.84190 nm*eV
"""
const fac_wleneV = 1239.84190   # wlen [n.m.] * energy [eV]
"""
fac_au2eV = 27.2113845
convert factor: photon energy [a.u.] ↔ photon energy [eV]
1 a.u. = 27.2113845 eV
"""
const fac_au2eV = 27.2113845# energy [a.u.] → energy [eV]
## convert functions
# laser amplitude
i2e(fint) = sqrt(fint / fac_i2e)
e2i(famp) = famp^2 * fac_i2e
vnm2e(famp_vnm) = famp_vnm / fac_vnm2e
e2vnm(famp_au) = famp_au * fac_vnm2e
vnm2i(famp_vnm) = e2i(vnm2e(famp_vnm))
i2vnm(fint) = e2vnm(i2e(fint))
"""
Ponderomotive energy (intensity [W/cm²], wavelength [nm]) → Up [a.u.]
"""
Up(fint, wlen_nm) = i2e(fint)^2 / (4 * wlen2au(wlen_nm)^2)
"""
convert (vector potential [a.u.], wavelength [nm]) → intensity [W/cm²]
"""
A2i(A, wlen) = e2i(A * wlen2au(wlen))
"""
convert (intensity [W/cm²], wavelength [nm]) → vector potential [a.u.]
"""
i2A(fint, wlen) = i2e(fint) / wlen2au(wlen)
## time
au2as(t_au) = fac_au2as * t_au
as2au(t_as) = t_as / fac_au2as
## wavelength and photon energy
wlen2au(wlen_nm) = wlen2eV(wlen_nm) / fac_au2eV
au2wlen(freq_au) = fac_au2eV * freq_au
wlen2eV(wlen_nm) = fac_wleneV / wlen_nm
eV2wlen(ene_eV) = fac_wleneV / ene_eV
au2eV(ene_au) = fac_au2eV * ene_au
eV2au(ene_eV) = fac_au2eV / ene_eV
wlen2period(wlen_nm) = au2as(2π / wlen2au(wlen_nm))
end
