#%%
from astropy import units as u

from astropy import constants as const

h = const.h
#%%
freq = 8.4 * u.GHz
energy = (freq * h) .to(u.MeV)
print(energy)

flux = 6.8 * u.mJy
flux_2 = flux.to(u.erg /(u.s * u.Hz * u.cm**2))
print(f"flux{flux_2}")

E_t  = freq * flux_2
print(f"E_t{E_t}")
#%%
flux_erg = E-23 #* u.erg /(u.s * u.Hz * u.cm**2)
print(flux_erg)
E_t  = flux_erg * 8.4e9
print(f"E_t{E_t}")


# %%
print((0.068 * u.Jy).to(u.erg /(u.s * u.Hz * u.cm**2)))
# %%

from astropy import units as u

from astropy import constants as const

h = const.h

unit = u.erg /(u.s * u.Hz)# * u.cm**2)
print(u.Hz)
               

flux_erg = ((69-16)*1e-3 * u.Jy).to(u.erg /(u.s * u.Hz * u.cm**2))

print(E_t)
E_t  = flux_erg * 14.7e9
print(f"E_t:{E_t}")
# %%

flux_erg = ((69-16)* u.Jy *  14.7 * u.GHz ).to(u.erg /(u.s * u.cm**2))

print(flux_erg)
# E_t  = flux_erg * 14.7e9
# print(f"E_t:{E_t}")
# %%
