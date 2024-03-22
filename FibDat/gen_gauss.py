#! ./venv/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm
import warnings
import glob


plt.rcParams.update({
    'font.sans-serif': 'serif',
    'font.family': 'serif'
})


def path_finder(file_glob):
    return glob.glob(file_glob)

def adc_to_size(adc, c0, c1):
    return (adc-c0)/c1

def proc_cal_file(path):
    with open(path) as f:
        dat = f.readlines()
    dat = [i.replace('\n', '') for i in dat]
    dat = [float(i.split(" = ")[-1]) for i in dat]
    return {'c0': dat[1], 'c1':dat[0]}

def name_from_path(path):
    return path

def tnormal(myclip_a, myclip_b, loc=0.2, scale=1, size=1, npts=1000):
    a, b = (myclip_a - loc) / scale, (myclip_b - loc) / scale
    var, skew = truncnorm.stats(a, b, moments='vs')
    x = np.linspace(truncnorm.ppf(0.01, a, b, loc=loc, scale=scale),
                    truncnorm.ppf(0.99, a, b, loc=loc, scale=scale), 100)
    z = truncnorm(a, b, loc=loc, scale=scale)
    r = truncnorm.rvs(a, b, loc=loc, scale=scale, size=npts)
    return x, z, r, var, skew

def ed_korolev(counts, bins):
    gmd = [np.sqrt(abs(i)*abs(j)) for i, j in zip(bins[:-1], bins[1:])]
    x = [(i*j**3, i*j**2) for i, j in zip(counts, gmd)]
    return sum([i for i, j in x])/sum([j for i, j in x])

aero_cal = [path_finder(f'{i}/*-AD-*_CalData*')
            for i in path_finder('Proc/aero/*')]
aero_cal = [i for r in aero_cal for i in r]
fib_cal = [path_finder(f'{i}/*-AD-*_CalData*')
           for i in path_finder('Proc/fib/*')]
fib_cal = [i for r in fib_cal for i in r]
fib_cal = [proc_cal_file(j) for j in fib_cal]
aero_cal = [proc_cal_file(j) for j in aero_cal]

sigma = 500
n_plot = 5
n_bins = 15
n_points = 1000000
zmin, zmax = 0, 4095
mu_arr = np.linspace(500, 3500, n_plot)
fig, ax_arr = plt.subplots(n_plot, 3, figsize=(4.72, 4.72))
plt.subplots_adjust(wspace=0, hspace=0)
ax_arr[0, 0].set_title("Generated Data", fontname='serif', fontsize=8)
ax_arr[0, 1].set_title("Fibre Calibration", fontname='serif', fontsize=8)
ax_arr[0, 2].set_title("Aerosol Calibration", fontname='serif', fontsize=8)

for i, mu in enumerate(mu_arr):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning,
                                message='divide by zero encountered in log')
        x, z, r, v, s = tnormal(zmin, zmax, loc=mu, scale=sigma,
                                               size=10_000, npts=n_points)
    ax1 = ax_arr[i, 0]
    ax2 = ax_arr[i, 1]
    ax3 = ax_arr[i, 2]

    ax1.plot(x, z.pdf(x), 'k-', lw=0.5, label='pdf')
    ax1.hist(r, density=True, bins=n_bins, histtype='stepfilled', alpha=0.2,
            color='k', range=(zmin, zmax))
    h, e = np.histogram(r, bins=n_bins, density=True, range=(zmin, zmax))
    ed = ed_korolev(h, e)
    ax1.plot([ed, ed], [0, 0.1], 'k-', lw=0.5, label='ED')

    ax1.set_ylim([0, 0.001])
    ax1.set_xlim([zmin, zmax])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.get_yaxis().set_ticks([])
    ax1.get_xaxis().set_ticks([0, 1500, 3000])

    edf = []
    eda = []
    for fc, ac in zip(fib_cal, aero_cal):
        ef = [adc_to_size(i, fc['c0'], fc['c1'])*1e12 for i in e]
        edf.append(ed_korolev(h, ef))
        ax2.bar(ef[:-1], h, width=np.diff(ef),
               edgecolor="none", align="edge", color='k', alpha=0.05)
        ea = [adc_to_size(i, ac['c0'], ac['c1'])*1e12 for i in e]
        eda.append(ed_korolev(h, ea))
        ax3.bar(ea[:-1], h, width=np.diff(ea),
               edgecolor="none", align="edge", color='k', alpha=0.05)

    styles = {"linewidth": 0.7, "color": 'k'}
    pstyles = {"marker": "x"}
    ax2.boxplot(edf, vert=False, positions=[0.5e-3], widths=[4e-4],
               boxprops=styles, flierprops=pstyles, medianprops=styles,
               capprops=styles, whiskerprops=styles)
    ax3.boxplot(eda, vert=False, positions=[0.5e-3], widths=[4e-4],
               boxprops=styles, flierprops=pstyles, medianprops=styles,
               capprops=styles, whiskerprops=styles)

    ax2.set_xlim([0, 2e2])
    ax2.set_ylim([0, 1e-3])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.get_yaxis().set_ticks([])

    ax3.set_xlim([0, 2e2])
    ax3.set_ylim([0, 1e-3])
    ax3.spines['top'].set_visible(False)
    #ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.get_yaxis().set_ticks([])
    ax3.get_xaxis().set_ticks([100, 200])

    if i == 0:
        ax1.spines['top'].set_visible(True)
        ax2.spines['top'].set_visible(True)
        ax3.spines['top'].set_visible(True)
    if i != (n_plot-1):
        ax1.xaxis.set_tick_params(labelbottom=False)
        ax2.xaxis.set_tick_params(labelbottom=False)
        ax3.xaxis.set_tick_params(labelbottom=False)
        ax1.set_xticks([])
        ax2.set_xticks([])
        ax3.set_xticks([])
        ax1.spines['bottom'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax3.spines['bottom'].set_visible(False)

ax_arr[-1, 0].set_xlabel("12-bit ADC Value", fontname="serif", fontsize=8)
ax_arr[-1, 1].set_xlabel(r'$\sigma$ ($\mu$m$^{2}$)',
                         fontname="serif", fontsize=8)
ax_arr[-1, 2].set_xlabel(r'$\sigma$ ($\mu$m$^{2}$)',
                         fontname="serif", fontsize=8)
fig.suptitle('Bulk Response to Generated Data',
             fontname='serif', weight='bold', fontsize=12)
ax_arr[-1, 0].get_xaxis().set_label_coords(0.5,-0.42)
ax_arr[-1, 1].get_xaxis().set_label_coords(0.5,-0.36)
ax_arr[-1, 2].get_xaxis().set_label_coords(0.5,-0.36)
plt.show()

