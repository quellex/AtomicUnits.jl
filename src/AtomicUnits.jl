# --------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
# unit conversion
module AtomicUnits
export fac_i2e, fac_vnm2au, fac_wlenev, fac_au2eV
export i2e, e2i, vnm2au, au2vnm, au2as, as2au, wlen2au, au2wlen, wlen2eV, eV2wlen, eV2au, au2eV
const fac_i2e  = 3.50944506e16   # intensity [W/cm^2] → amplitude [a.u.]
const fac_vnm2au = 514.2	# amplitude [V/nm] → amplitude [a.u.]
const fac_au2as  = 24.1899	# time [a.u.] → time[atto sec.]
const fac_wleneV = 1239.84190   # wlen [n.m.] * energy [eV]
const fac_au2eV = 27.2113845	# energy [a.u.] → energy [eV]
# convert functions
i2e(fint) = sqrt(fint / fac_i2e)
e2i(famp) = famp^2 * fac_i2e
vnm2au(famp_vnm) = fac_vnm2au * famp_vnm
au2vnm(famp_au) = famp_au / fac_vnm2au
au2as(t_au) = fac_au2as * t_au
as2au(t_as) = t_as / fac_au2as
wlen2au(wlen_nm) = wlen2eV(wlen_nm) / fac_au2eV
au2wlen(freq_au) = fac_au2eV * freq_au
wlen2eV(wlen_nm) = fac_wleneV / wlen_nm
eV2wlen(ene_eV) = fac_wleneV / ene_eV
au2eV(ene_au) = fac_au2eV/ene_au
eV2au(ene_eV) = fac_au2eV/ene_eV
end