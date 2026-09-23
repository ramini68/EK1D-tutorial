"""Render article figures directly from retained, unfiltered EK1D run arrays.

Reference package:
    python plot_article_figures.py --data-root data --output figures_rerun
Reproduction matrix:
    python plot_article_figures.py --data-root rerun --prefix '' --output figures_rerun
No simulation, interpolation, spatial binning or current filtering is performed.
"""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit_metrics import read_run
from EK1D_tutorial import snapshot_prefix


COLORS = ['#777777', '#427aa1', '#b27920', '#a34545']


def finish(fig, output, stem):
    for ax in fig.axes:
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
    for extension in ('png', 'pdf', 'svg'):
        fig.savefig(Path(output) / f'{stem}.{extension}', dpi=400,
                    bbox_inches='tight', metadata={'Creator': 'EK1D article plotting script'}
                    if extension in ('pdf', 'svg') else None)
    plt.close(fig)


def main_figure(folder, output):
    r = read_run(folder)
    a = r['arrays']
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.65), layout='constrained')
    for t, color, style in zip((0., 21600., 86400., 129600.), COLORS, (':', '--', '-.', '-')):
        key = snapshot_prefix(t) + '_pH'
        axes[0].plot(a['x'] * 100, a[key], color=color, ls=style, lw=1.3, label=f'{t/3600:g} h')
    axes[0].set(xlabel='Distance from anode (cm)', ylabel='pH', xlim=(0, 10),
                title='(a) pH evolution')
    axes[0].legend(frameon=False, fontsize=7, loc='upper left')
    inv = a['phase_inventory_series']
    for k, color, style, label in zip(range(3), ['#b27920', '#438560', '#a34545'],
                                    ['-', '--', '-.'], ['Aqueous', 'Sorbed', 'Precipitated']):
        axes[1].plot(a['times_s']/3600, inv[:, k], color=color, ls=style, lw=1.3, label=label)
    axes[1].plot(a['times_s']/3600, inv.sum(axis=1), ':', color='#555555', lw=1.1, label='Total Cu')
    axes[1].set(xlabel='Time (h)', ylabel='Cu inventory (mol m$^{-2}$)', xlim=(0, 36),
                ylim=(0, .245), title='(b) Cu phase inventories')
    axes[1].legend(frameon=False, fontsize=7, loc='upper right', bbox_to_anchor=(1, .89))
    axes[2].plot(a['times_s']/3600, a['signed_charge_C_m2']/1000, color='#244d70', lw=1.3)
    axes[2].set(xlabel='Time (h)', ylabel='Cumulative charge (kC m$^{-2}$)',
                xlim=(0, 36), ylim=(0, None), title='(c) Cumulative charge')
    finish(fig, output, 'Figure7_Revised')


def grid_figure(folders, output):
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.4), layout='constrained')
    for folder, color, style in zip(folders, ['#787878', '#b27920', '#244d70'], [':', '--', '-']):
        r = read_run(folder)
        a = r['arrays']; label = f'{r["diagnostics"]["parameters"]["N"]} cells'
        axes[0, 0].plot(a['x']*100, a['t129600_pH'], color=color, ls=style, lw=1.3, label=label)
        axes[0, 1].stairs(a['t129600_Cu_ppt'], r['edges']*100,
                          color=color, ls=style, lw=1.3, label=label)
        axes[1, 0].plot(a['times_s']/3600, a['currents'], color=color, ls=style, lw=1.15, label=label)
        axes[1, 1].plot(a['times_s']/3600, a['signed_charge_C_m2']/1000,
                          color=color, ls=style, lw=1.3, label=label)
    axes[0, 0].set(xlabel='Distance from anode (cm)', ylabel='pH', title='(a) Final pH profiles')
    axes[0, 1].set(xlabel='Distance from anode (cm)', ylabel='Precipitated Cu (mol m$^{-3}$)',
                   title='(b) Raw cell means: front detail', xlim=(5.8, 6.6))
    axes[1, 0].set(xlabel='Time (h)', ylabel='Current density (A m$^{-2}$)',
                   title='(c) Full instantaneous current', xlim=(0, 36), yscale='log')
    axes[1, 1].set(xlabel='Time (h)', ylabel='Cumulative charge (kC m$^{-2}$)',
                   title='(d) Time integrated current', xlim=(0, 36))
    for ax in axes.flat:
        ax.legend(frameon=False, fontsize=7)
    finish(fig, output, 'FigureS1_GridSensitivity')


def time_figure(folder_a, folder_b, output):
    a, b = read_run(folder_a)['arrays'], read_run(folder_b)['arrays']
    if not np.array_equal(a['times_s'], b['times_s']):
        raise ValueError('Temporal curves must have identical saved output times.')
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), layout='constrained')
    axes[0].plot(a['times_s']/3600, a['currents'], '--', color='#b27920', lw=1.3, label='CFL 0.2')
    axes[0].plot(b['times_s']/3600, b['currents'], color='#244d70', lw=1.1, label='CFL 0.1')
    axes[0].set(xlabel='Time (h)', ylabel='Current density (A m$^{-2}$)', xlim=(0, 6),
                title='(a) Temporal comparison: first 6 h')
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(a['times_s']/3600, a['currents']-b['currents'], color='#244d70', lw=1.)
    axes[1].axhline(0, color='#777777', lw=.6)
    axes[1].set(xlabel='Time (h)', ylabel='Current difference (A m$^{-2}$)', xlim=(0, 36),
                title='(b) Difference: full 36 h')
    finish(fig, output, 'FigureS2_TimeSensitivity')


def render_article_figures(data_root, output, prefix='rk2_'):
    from assess_publication import validate_pair
    data_root, output = Path(data_root), Path(output)
    names = ['N100_cfl02', 'N200_cfl02', 'N400_cfl02', 'N400_cfl01']
    folders = [data_root/(prefix + name) for name in names]
    for folder, cells, cfl in zip(folders, (100, 200, 400, 400), (.2, .2, .2, .1)):
        p = read_run(folder)['diagnostics']['parameters']
        if any(p[k] != value for k, value in {'N': cells, 'cfl': cfl, 'L': .1,
                'tend': 129600., 'spatial_order': 2, 'time_order': 2}.items()):
            raise ValueError(f'{folder} is not the declared article run configuration.')
    validate_pair(folders[0], folders[1], 'spatial')
    validate_pair(folders[1], folders[2], 'spatial')
    validate_pair(folders[2], folders[3], 'temporal')
    output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9,
                         'legend.fontsize': 7, 'svg.fonttype': 'none', 'pdf.fonttype': 42}):
        main_figure(folders[-1], output)
        grid_figure(folders[:3], output)
        time_figure(folders[-2], folders[-1], output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', default='data')
    parser.add_argument('--prefix', default='rk2_')
    parser.add_argument('--output', default='figures_rerun')
    args = parser.parse_args()
    render_article_figures(args.data_root, args.output, args.prefix)
