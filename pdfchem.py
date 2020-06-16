#Based on Bisbas, Schruba & van Dishoeck, 2019, MNRAS, 485, 3097
#Original version written in FORTRAN by Thomas Bisbas
#Python version written by Thomas Bisbas with thanks to Pierre Nuernberger
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
import sys

s = float(sys.argv[1])
av_bar = float(sys.argv[2])

m = np.log(av_bar) - s ** 2 / 2.
print 'Mean =', av_bar  # np.exp(m + s ** 2 / 2.)
print 'Median =', np.exp(m)
print 'Mode =', np.exp(m - s ** 2)


def pdf(x, s, m):
	pdf_out = (1. / s / np.sqrt(2. * np.pi) / x) * np.exp(-(np.log(x) - m) ** 2 / 2. / s ** 2)
	return pdf_out


def av_bnd(l, x):
	av_bnd_out = (x / nh(x)) * (nh(x) + l * nh(x) * (1 - (np.log10(nh(x) / 0.1)) / (np.log10(2e6 / 0.1))))
	return av_bnd_out


def nh(av):
	nh_out = (np.log(av / 0.05) / 1.6) ** (1. / 0.12)
	return nh_out  # cm-3


def find_nearest(array, value):
	idx = (np.abs(array - value)).argmin()
	return array[idx], idx


pwrx = np.arange(-1, 2, 0.01)
av = 10 ** pwrx

#plotting av-pdf
plt.figure()
plt.loglog(av, pdf(av, s, m), '-')
plt.axvline(x=0.25, color='black')
plt.axvspan(1e-4, 0.25, alpha=0.3, color='black')
plt.xlim([1e-2, 1e2])
plt.ylim([1e-5, 5])
plt.xlabel(r'$A_V$ [mag]')
plt.ylabel(r'$A_Vp(A_V;\mu,\sigma)$')
plt.savefig('Av-pdf.png', bbox_inches='tight')

#plotting av-nh and bounds
plt.figure()
plt.loglog(nh(av), av, '-', color='C0')
plt.loglog(nh(av),av_bnd(-0.5, av), '-.', color='C0')
plt.loglog(nh(av),av_bnd(1, av), '--', color='C0')
plt.axvline(x=1, color='black')
plt.axvspan(1e-4, 1, alpha=0.3, color='black')
plt.xlim([1e-4,1e7])
plt.ylim([1e-2,1e3])
plt.xlabel(r'$n_{\rm H}$ [cm$^{-3}$]')
plt.ylabel(r'$A_V$ [mag]')
plt.savefig('Av-nH.png', bbox_inches='tight')

# getting the index for pdr sim
ipdr = np.log10(nh(av)) * 100
for i in range(ipdr.size):  # limits to within the available pdr sims
	if ipdr[i] < 0:
		ipdr[i] = 0
	if ipdr[i] > 600:
		ipdr[i] = 600

# opening hdf5 files
zcr = h5py.File('zcrs.hdf5', 'r')
fuv = h5py.File('fuvs.hdf5', 'r')
Z = h5py.File('Zs.hdf5', 'r')

# creating indices for plotting
zcrind = [1e-17, 1e-16, 1e-15, 1e-14]
fuvind = [1, 10, 100, 1000]
Zind = [0.1, 0.2, 0.5, 1]

chemid = [3, 4, 5, 6, 7] #CII, CI, CO, H2, HI

# same spacing for all sims, so we get the av list from one pdr sim
dav = zcr['data'][0, 0, :, 1]

# Cosmic-rays (zcr)
chemavg, chemavg_u, chemavg_b = [[], [], []]

chem = []
Ntot = []
for i in range(av.size):
	# maximum av integration for a particular sim
	avmax, pos = find_nearest(dav, av[i])
	length = zcr['data'][int(ipdr[i]), :, :pos, 0]
	density = zcr['data'][int(ipdr[i]), :, :pos, 2]
	tchem = zcr['data'][int(ipdr[i]), :, :pos, chemid]
	Ntot.append(np.trapz(density, length, axis=1))
	chem.append(np.trapz(tchem * density[:, :, None], length[:, :, None], axis=1))

chem = np.asarray(chem)
Ntot = np.asarray(Ntot)

chemavg = np.sum(chem * pdf(av, s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av, s, m)[:, None], axis=0)[:, None]
chemavg_u = np.sum(np.asarray(chem) * pdf(av_bnd(1.0, av),s,m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av_bnd(1.0, av),s,m)[:, None], axis=0)[:, None]
chemavg_b = np.sum(np.asarray(chem) * pdf(av_bnd(-0.5, av),s,m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av_bnd(-0.5, av),s,m)[:, None], axis=0)[:, None]

chemavg = chemavg.T
chemavg_b = chemavg_b.T
chemavg_u = chemavg_u.T

# plotting
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 9))
ax1.loglog(zcrind, chemavg[4], '-', color='C2', label='HI')
ax1.loglog(zcrind, chemavg[3], '-', color='C0', label='H2')
ax1.fill_between(zcrind, chemavg_b[4], chemavg_u[4], facecolor='C2', alpha=0.3)
ax1.fill_between(zcrind, chemavg_b[3], chemavg_u[3], facecolor='C0', alpha=0.3)
ax1.set_xlabel(r'$\zeta_{\rm CR}$ (s$^{-1}$)')
ax1.set_ylabel('Fractional Abundance')
ax1.set_title(r'$\chi/\chi_0=1$; $Z=1\,{\rm Z}_{\odot}$')
ax1.legend(loc='best')
ax2.loglog(zcrind, chemavg[2], '-', color='C0', label='CO')
ax2.loglog(zcrind, chemavg[1], '-', color='C2', label='CI')
ax2.loglog(zcrind, chemavg[0], '-', color='C1', label='CII')
ax2.fill_between(zcrind, chemavg_b[2], chemavg_u[2], facecolor='C0', alpha=0.3)
ax2.fill_between(zcrind, chemavg_b[1], chemavg_u[1], facecolor='C2', alpha=0.3)
ax2.fill_between(zcrind, chemavg_b[0], chemavg_u[0], facecolor='C1', alpha=0.3)
ax2.set_xlabel(r'$\zeta_{\rm CR}$ (s$^{-1}$)')
ax2.set_ylabel('Fractional Abundance')
ax2.legend(loc='best')
plt.savefig('zcr.png', bbox_inches='tight')

# Interstellar radiation (fuv)
chemavg, chemavg_u, chemavg_b = [[], [], []]

chem = []
Ntot = []
for i in range(av.size):
	# maximum av integration for a particular sim
	avmax, pos = find_nearest(dav, av[i])
	length = fuv['data'][int(ipdr[i]), :, :pos, 0]
	density = fuv['data'][int(ipdr[i]), :, :pos, 2]
	tchem = fuv['data'][int(ipdr[i]), :, :pos, chemid]
	Ntot.append(np.trapz(density, length, axis=1))
	chem.append(np.trapz(tchem * density[:, :, None], length[:, :, None], axis=1))

chem = np.asarray(chem)
Ntot = np.asarray(Ntot)

chemavg = np.sum(chem * pdf(av, s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av, s, m)[:, None], axis=0)[:, None]
chemavg_u = np.sum(np.asarray(chem) * pdf(av_bnd(1.0, av), s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av_bnd(1.0, av), s, m)[:, None], axis=0)[:, None]
chemavg_b = np.sum(np.asarray(chem) * pdf(av_bnd(-0.5, av), s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av_bnd(-0.5, av), s, m)[:, None], axis=0)[:, None]

chemavg = chemavg.T
chemavg_b = chemavg_b.T
chemavg_u = chemavg_u.T

# plotting
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 9))
ax1.loglog(fuvind, chemavg[4], '-', color='C2', label='HI')
ax1.loglog(fuvind, chemavg[3], '-', color='C0', label='H2')
ax1.fill_between(fuvind, chemavg_b[4], chemavg_u[4], facecolor='C2', alpha=0.3)
ax1.fill_between(fuvind, chemavg_b[3], chemavg_u[3], facecolor='C0', alpha=0.3)
ax1.set_xlabel(r'$\chi/\chi_0$')
ax1.set_ylabel('Fractional Abundance')
ax1.set_title(r'$\zeta_{\rm CR}=10^{-16}\,{\rm s}^{-1}$; $Z=1\,{\rm Z}_{\odot}$')
ax1.legend(loc='best')
ax2.loglog(fuvind, chemavg[2], '-', color='C0', label='CO')
ax2.loglog(fuvind, chemavg[1], '-', color='C2', label='CI')
ax2.loglog(fuvind, chemavg[0], '-', color='C1', label='CII')
ax2.fill_between(fuvind, chemavg_b[2], chemavg_u[2], facecolor='C0', alpha=0.3)
ax2.fill_between(fuvind, chemavg_b[1], chemavg_u[1], facecolor='C2', alpha=0.3)
ax2.fill_between(fuvind, chemavg_b[0], chemavg_u[0], facecolor='C1', alpha=0.3)
ax2.set_xlabel(r'$\chi/\chi_0$')
ax2.set_ylabel('Fractional Abundance')
ax2.legend(loc='best')
plt.savefig('fuv.png', bbox_inches='tight')

# Metallicity (Z)
chemavg, chemavg_u, chemavg_b = [[], [], []]

chem = []
Ntot = []
for i in range(av.size):
	# maximum av integration for a particular sim
	avmax, pos = find_nearest(dav, av[i])
	length = Z['data'][int(ipdr[i]), :, :pos, 0]
	density = Z['data'][int(ipdr[i]), :, :pos, 2]
	tchem = Z['data'][int(ipdr[i]), :, :pos, chemid]
	Ntot.append(np.trapz(density, length, axis=1))
	chem.append(np.trapz(tchem * density[:, :, None], length[:, :, None], axis=1))

chem = np.asarray(chem)
Ntot = np.asarray(Ntot)

chemavg = np.sum(chem * pdf(av, s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av, s, m)[:, None], axis=0)[:, None]
chemavg_u = np.sum(np.asarray(chem) * pdf(av_bnd(1.0, av), s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av_bnd(1.0, av), s, m)[:, None], axis=0)[:, None]
chemavg_b = np.sum(np.asarray(chem) * pdf(av_bnd(-0.5, av), s, m)[:, None, None], axis=0) / np.sum(Ntot * pdf(av_bnd(-0.5, av), s, m)[:, None], axis=0)[:, None]

chemavg = chemavg.T
chemavg_b = chemavg_b.T
chemavg_u = chemavg_u.T

# plotting
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 9))
ax1.semilogy(Zind, chemavg[4], '-', color='C2', label='HI')
ax1.semilogy(Zind, chemavg[3], '-', color='C0', label='H2')
ax1.fill_between(Zind, chemavg_b[4], chemavg_u[4], facecolor='C2', alpha=0.3)
ax1.fill_between(Zind, chemavg_b[3], chemavg_u[3], facecolor='C0', alpha=0.3)
ax1.set_xlabel(r'Metallicity ($Z_{\odot}$)')
ax1.set_ylabel('Fractional Abundance')
ax1.set_title(r'$\chi/\chi_0=1$; $\zeta_{\rm CR}=10^{-16}\,{\rm s}^{-1}$')
ax1.legend(loc='best')
ax2.semilogy(Zind, chemavg[2], '-', color='C0', label='CO')
ax2.semilogy(Zind, chemavg[1], '-', color='C2', label='CI')
ax2.semilogy(Zind, chemavg[0], '-', color='C1', label='CII')
ax2.fill_between(Zind, chemavg_b[2], chemavg_u[2], facecolor='C0', alpha=0.3)
ax2.fill_between(Zind, chemavg_b[1], chemavg_u[1], facecolor='C2', alpha=0.3)
ax2.fill_between(Zind, chemavg_b[0], chemavg_u[0], facecolor='C1', alpha=0.3)
ax2.set_xlabel(r'Metallicity ($Z_{\odot}$)')
ax2.set_ylabel('Fractional Abundance')
ax2.legend(loc='best')
plt.savefig('Z.png', bbox_inches='tight')

# plot all
plt.show()
