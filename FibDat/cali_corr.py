#! ./venv/bin/python
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy import stats
from scipy.stats import truncnorm
import warnings


plt.rcParams.update({
    'font.sans-serif': 'serif',
    'font.family': 'serif'
})


def path_finder(file_glob):
    return glob.glob(file_glob)

def get_ucass_name(path):
    return path.split('/')[-2]

def match_file(path):
    if "aero" in path:
        t = "fib"
        un = f'{get_ucass_name(path)}_2'
    elif "fib" in path:
        t = "aero"
        un = get_ucass_name(path).replace('_2', '')
    else:
        raise ValueError
    file_glob = f'Proc/{t}/{un}/*_CalData*'
    return path_finder(file_glob)

def adc_to_size(adc, c0, c1):
    return (adc-c0)/c1

def proc_cal_file(path):
    with open(path) as f:
        dat = f.readlines()
    dat = [i.replace('\n', '') for i in dat]
    dat = [float(i.split(" = ")[-1]) for i in dat]
    return {'c0': dat[1], 'c1':dat[0]}

def rma_regression(x_list, y_list):
    """
    rma_regression -- http://doi.wiley.com/10.1002/9781118445112.stat07912
    """
    sx = stats.tstd(x_list)
    sy = stats.tstd(y_list)

    avg_x = stats.tmean(x_list)
    avg_y = stats.tmean(y_list)

    m = sy/sx
    c = avg_y - m * avg_x
    r = stats.pearsonr(x_list, y_list)[0]
    r2 = r ** 2

    x12 = [min(x_list), max(x_list)]
    y12 = [m * min(x_list) + c, m * max(x_list) + c]

    return x12, y12, r2, m

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

def y_comp(x, x12, y12):
    x1, x2 = x12
    y1, y2 = y12
    m = (y2-y1)/(x2-x1)
    c = y1 - m*x1
    return m*x+c

def rot_comp(line, ax):
    xdata, ydata = line.get_data()
    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]
    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))
    m = (sp2[1] - sp1[1])/(sp2[0] - sp1[0])
    r = np.arctan(m)*(180/np.pi)
    print(r)
    return r


aero_cal_files = [path_finder(f'{i}/*_CalData*')
            for i in path_finder('Proc/aero/*')]
aero_cal_files = [i for r in aero_cal_files for i in r]
fib_cal_files = [match_file(i) for i in aero_cal_files]
fib_cal_files = [i for r in fib_cal_files for i in r]

fib_cal = [proc_cal_file(j) for j in fib_cal_files]
aero_cal = [proc_cal_file(j) for j in aero_cal_files]

fig, ax_arr = plt.subplots(2, 2, figsize=(4.72, 4.72))
plt.subplots_adjust(wspace=0.4, hspace=0.4)
ax_arr[0, 0].set_title("c0 Fibre Versus Aerosol", fontsize=8, style='italic')
ax_arr[0, 1].set_title("c1/1E+13 Fibre Versus Aerosol", fontsize=8, style='italic')
ax_arr[1, 0].set_title("ED Fibre Versus Aerosol", fontsize=8, style='italic')
ax_arr[1, 1].set_title("c0 Versus c1/1E+13", fontsize=8, style='italic')

eda_store = []
edf_store = []
for a, f, af in zip(aero_cal, fib_cal, aero_cal_files):
    if '-AA-' in af:
        style = {'color': 'k', 'fillstyle':'none', 'marker':'o',
                 'markersize':'4'}
    else:
        style = {'color': 'k', 'fillstyle':'none', 'marker':'x',
                 'markersize':'4'}

    ax_arr[0, 0].plot(a['c0'], f['c0'], **style)
    ax_arr[0, 1].plot(a['c1']/1e13, f['c1']/1e13, **style)

    for mu in np.linspace(500, 3500, 5):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning,
                                    message='divide by zero encountered in log'
                                   )
            _, _, r, _, _ = tnormal(0, 4095, loc=mu,
                                    scale=500, size=10_000, npts=100000)
        h, e = np.histogram(r, bins=16, density=True, range=(0, 4095))
        ef = [adc_to_size(i, f['c0'], f['c1'])*1e12 for i in e]
        edf = ed_korolev(h, ef)
        edf_store.append(edf)
        ea = [adc_to_size(i, a['c0'], a['c1'])*1e12 for i in e]
        eda = ed_korolev(h, ea)
        eda_store.append(eda)
        ax_arr[1, 0].plot(eda, edf, **style)

    ax_arr[1, 1].plot(a['c1']/1e13, a['c0'], **style)
    ax_arr[1, 1].plot(f['c1']/1e13, f['c0'], **style)

x120, y120, r0, m0 = rma_regression([i['c0'] for i in aero_cal],
                                    [i['c0'] for i in fib_cal])
x121, y121, r1, m1 = rma_regression([i['c1']/1e13 for i in aero_cal],
                                    [i['c1']/1e13 for i in fib_cal])
x123, y123, r3, m3 = rma_regression(eda_store, edf_store)
x124, y124, r4, m4 = rma_regression([i['c1']/1e13 for i in aero_cal],
                                    [i['c0'] for i in aero_cal])
x125, y125, r5, m5 = rma_regression([i['c1']/1e13 for i in fib_cal],
                                    [i['c0'] for i in fib_cal])

ax_arr[0, 0].annotate(format('m=%.1f\nr=%.1f' % (m0, r0)), xy=(0.05, 0.80),
                      xycoords='axes fraction', fontsize='8')
ax_arr[0, 1].annotate(format('m=%.1f\nr=%.1f' % (m1, r1)), xy=(0.05, 0.80),
                      xycoords='axes fraction', fontsize='8')
ax_arr[1, 0].annotate(format('m=%.1f\nr=%.1f' % (m3, r3)), xy=(0.05, 0.80),
                      xycoords='axes fraction', fontsize='8')
ax_arr[1, 1].annotate(format('m=%.1f\nr=%.1f' % (m4, r4)), xy=(0.6, 0.05),
                      xycoords='axes fraction', fontsize='8')
ax_arr[1, 1].annotate(format('m=%.1f\nr=%.1f' % (m5, r5)), xy=(0.05, 0.80),
                      xycoords='axes fraction', fontsize='8')

style = {'color': 'k', 'marker':'none', 'linewidth':'0.7'}
ax_arr[0, 0].plot(x120, y120, **style, linestyle='-')
ax_arr[0, 1].plot(x121, y121, **style, linestyle='-')
ax_arr[1, 0].plot(x123, y123, **style, linestyle='-')
l4, = ax_arr[1, 1].plot(x124, y124, **style, linestyle='-')
l5, = ax_arr[1, 1].plot(x125, y125, **style, linestyle='-')

ax_arr[0, 0].set_xlabel("Aerosol Calibration", fontsize=8)
ax_arr[0, 0].set_ylabel("Fibre Calibration", fontsize=8)
ax_arr[0, 1].set_xlabel("Aerosol Calibration", fontsize=8)
ax_arr[0, 1].set_ylabel("Fibre Calibration", fontsize=8)
ax_arr[1, 0].set_xlabel(r'Aerosol $\sigma$ ($\mu$m$^{2}$)', fontsize=8)
ax_arr[1, 0].set_ylabel(r'Fibre $\sigma$ ($\mu$m$^{2}$)', fontsize=8)
ax_arr[1, 1].set_xlabel("c1", fontsize=8)
ax_arr[1, 1].set_ylabel("c0", fontsize=8)

ax_arr[0, 0].tick_params(axis='both', which='major', labelsize=6)
ax_arr[0, 0].tick_params(axis='both', which='minor', labelsize=6)
ax_arr[0, 1].tick_params(axis='both', which='major', labelsize=6)
ax_arr[0, 1].tick_params(axis='both', which='minor', labelsize=6)
ax_arr[1, 0].tick_params(axis='both', which='major', labelsize=6)
ax_arr[1, 0].tick_params(axis='both', which='minor', labelsize=6)
ax_arr[1, 1].tick_params(axis='both', which='major', labelsize=6)
ax_arr[1, 1].tick_params(axis='both', which='minor', labelsize=6)

ax_arr[0, 0].set_aspect('equal', adjustable='datalim')
ax_arr[0, 1].set_aspect('equal', adjustable='datalim')
ax_arr[1, 0].set_aspect('equal', adjustable='datalim')

arrowprops = dict(arrowstyle="->", connectionstyle=
                  "arc3,rad=1")
bbox = dict(boxstyle="round", fc="1.0")
ax_arr[1, 1].annotate("Fibre", xy=(3.1, y_comp(3.1, x125, y125)), xycoords='data',
                      xytext=(-27, 7), textcoords='offset points',
                      fontsize='8', bbox=bbox, arrowprops=arrowprops)
ax_arr[1, 1].annotate("Aerosol", xy=(4.5, y_comp(4.5, x124, y124)), xycoords='data',
                      xytext=(5, 25), textcoords='offset points',
                      fontsize='8', bbox=bbox, arrowprops=arrowprops)

ax_arr[0, 0].annotate('(a)', xy=(-0.3, 0.95),
                      xycoords='axes fraction', fontsize='8')
ax_arr[0, 1].annotate('(b)', xy=(-0.3, 0.95),
                      xycoords='axes fraction', fontsize='8')
ax_arr[1, 0].annotate('(c)', xy=(-0.3, 0.95),
                      xycoords='axes fraction', fontsize='8')
ax_arr[1, 1].annotate('(d)', xy=(-0.3, 0.95),
                      xycoords='axes fraction', fontsize='8')

for axis in ['top','bottom','left','right']:
    ax_arr[0, 0].spines[axis].set_linewidth(0.5)
    ax_arr[0, 1].spines[axis].set_linewidth(0.5)
    ax_arr[1, 0].spines[axis].set_linewidth(0.5)
    ax_arr[1, 1].spines[axis].set_linewidth(0.5)

fig.suptitle('Correlation Plots', weight='bold', fontsize=10)
plt.show()

