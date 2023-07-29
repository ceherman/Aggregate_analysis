import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy
from scipy import integrate, interpolate, misc, stats
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys
import os

sys.path.append("/home/chase/codes/python_functions/")
import plotting as my_plot

sys.path.append(os.path.realpath('./pyilt/'))
import ilt


def get_A_val(q, viscosity=0.890e-3, temp=298):
    # q in nm-1, viscosity in Pa s, temp in K
    return 1.38e-23*temp*q**2*1e27/(6*np.pi*viscosity) # nm/s

def get_gamma_from_radius(radius, q, viscosity=0.890e-3, temp=298):
    return get_A_val(q)/radius # s-1

def get_radius_from_gamma(gamma, q, viscosity=0.890e-3, temp=298):
    return get_A_val(q)/gamma # nm

class data_file():
    def __init__(self, file):
        self.file = file
        self.name = self.file[self.file.rfind('/')+1:self.file.rfind('.')]
        self.read_data()
        self.get_g1()
        return

    def read_data(self):
        with open(self.file) as f:
            self.lines = f.readlines()

        self.angle = float(self.lines[3].split()[-1])
        self.duration_s = float(self.lines[4].split()[-1])
        average_count_rate_1 = float(self.lines[10].split()[-1])
        average_count_rate_2 = float(self.lines[11].split()[-1])
        self.average_count_rate_kHz = (average_count_rate_1 + average_count_rate_2)/2.0

        self.n = 1.33 # refractive index
        self.wavelength = 660 # wavelength
        self.q = 4*np.pi*self.n/self.wavelength*np.sin(self.angle*np.pi/180/2)

        g2_titles_line = self.lines.index('Lag time (s)         g2-1\n')
        cr_titles_line = self.lines.index('Count Rate History (KHz)  CR CHA / CR CHB\n')
        self.g2_minus_one = pd.DataFrame([l.split() for l in self.lines[g2_titles_line+1:cr_titles_line-2]],
                                         columns=['lag_time_s', 'g2-1'], dtype=float)
        self.count_rate = pd.DataFrame([l.split() for l in self.lines[cr_titles_line+1:]],
                                       columns=['t_s', 'a_khz', 'b_khz'], dtype=float)
        return

    def make_count_rate_plot(self, fig=False, ax=False):
        if not fig and not ax:
            fig, ax = my_plot.instantiate_fig(xlabel='Time [s]', ylabel='Count rate [kHz]')
        ax.plot(self.count_rate.t_s, self.count_rate.a_khz)
        my_plot.set_layout(fig, ax)
        return fig, ax

    def make_count_rate_histogram(self, fig=False, ax=False):
        if not fig and not ax:
            fig, ax = my_plot.instantiate_fig(xlabel='Count rate [kHz]', ylabel='Frequency')
        n = ax.hist(self.count_rate.a_khz, bins=100)
        my_plot.set_layout(fig, ax)
        return fig, ax

    def make_g2_plot(self, fig=False, ax=False):
        fig, ax = my_plot.instantiate_fig(ylabel=r'$g^{(2)} (\tau) - 1$', xlabel=r'$\tau$ [s]')
        ax.set_xscale('log')
        ax.plot(self.g2_minus_one.lag_time_s, self.g2_minus_one['g2-1'])
        my_plot.set_layout(fig, ax)
        return fig, ax

    def get_g1(self):
        self.g1 = self.g2_minus_one[(0 < self.g2_minus_one['g2-1']) &\
                                    (self.g2_minus_one['g2-1'] < 1)].copy()
#         self.g1 = self.g2_minus_one[(1e-6 < self.g2_minus_one.lag_time_s) &
#                                     (self.g2_minus_one.lag_time_s < 1e-2) &
#                                     (0 < self.g2_minus_one['g2-1'])].copy()

        self.g1.reset_index(inplace=True, drop=True)
        self.g1['g1'] = np.sqrt(self.g1['g2-1']/self.g1.at[0, 'g2-1'])
        return

    def fit_g1(self, alpha=1e-6, low_r_bound=1, up_r_bound=800):
        # gamma_vals in s-1, G_gamma in s
        self.alpha = alpha
        bound = np.array([get_gamma_from_radius(up_r_bound, self.q),
                          get_gamma_from_radius(low_r_bound, self.q)])
        self.gamma_vals, self.G_gamma, self.res_lsq, self.res_reg, self.g1_fit, self.C, self.Q, self.S, self.W = ilt.ilt(self.g1.lag_time_s, self.g1.g1, bound, len(self.g1)-1, alpha)
        return

    def plot_g1_fit(self):
        fig, ax = my_plot.instantiate_fig(ylabel=r'$g^{(1)} (\tau)$',
                                          xlabel=r'$\tau$ [s]')
        ax.plot(self.g1.lag_time_s, self.g1.g1, label='Exp.')
        ax.plot(self.g1.lag_time_s, self.g1_fit, label='Fit')
        ax.set_xscale('log')
        my_plot.set_layout(fig, ax, legend=True)
        return fig, ax

    def get_rh_dist(self):
        self.rh_vals = get_radius_from_gamma(self.gamma_vals, self.q)
        self.F_rh = get_A_val(self.q)*self.G_gamma/(self.rh_vals**2)
        return

    def make_rh_dist_plot(self, fig=None, ax=None, label=False, normalize=False,
                          trim=False):
        if trim:
            start, end = 7, -22
        else:
            start, end = 0, len(self.rh_vals)
        x = self.rh_vals[start:end]
        y = self.F_rh[start:end]

        if normalize:
            if fig == None and ax == None:
                fig, ax = my_plot.instantiate_fig(xlabel=r'$R_h$  [nm]',
                                                ylabel=r'$F(R_h)/F_{max}$')
            if label:
                ax.plot(x, y/max(y), label=label)
            else:
                ax.plot(x, y/max(y))

        else:
            if fig == None and ax == None:
                fig, ax = my_plot.instantiate_fig(xlabel=r'$R_h$  [nm]',
                                                ylabel=r'$F(R_h)$  [nm$^{-1}$]')
            if label:
                ax.plot(x, y, label=label)
            else:
                ax.plot(x, y)

        ax.set_xscale('log')
        my_plot.set_layout(fig, ax)
        return fig, ax

    def make_gamma_dist_plot(self, fig=None, ax=None, label=False):
        if fig == None and ax == None:
            fig, ax = my_plot.instantiate_fig(xlabel=r'$\Gamma$  [kHz]',
                                              ylabel=r'$G(\Gamma)$  [ms]')
        if label:
            ax.plot(self.gamma_vals*1e-3, self.G_gamma*1e3, label=label)
        else:
            ax.plot(self.gamma_vals*1e-3, self.G_gamma*1e3)

        ax.set_xscale('log')
        my_plot.set_layout(fig, ax)
        return fig, ax

    def draw_L_curve(self, loglog=False, asymptotes=False, fix_aspect=False):
        fig, ax = my_plot.instantiate_fig(ylabel=fr'$\Vert  x  \Vert ^2$',
                                         xlabel=fr'$\Vert  Ax - g \Vert ^2$')
        if loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')

        ax.plot(self.resid_norms, self.gamma_norms, 'o-', fillstyle='none')
        ax.plot(self.resid_norms[self.index_knee], self.gamma_norms[self.index_knee], 'o-')

        if asymptotes:
            line_1_pts = self.line_1.slope * np.array(self.resid_norms) + self.line_1.intercept
            line_2_pts = self.line_2.slope * np.array(self.resid_norms) + self.line_2.intercept
            ax.plot(self.resid_norms, line_1_pts, 'k--')
            ax.plot(self.resid_norms, line_2_pts, 'k--')
            ax.set_ylim(-0.1*max(self.gamma_norms), 1.1*max(self.gamma_norms))

        my_plot.set_layout(fig, ax)
        if fix_aspect:
            ax.set_xlim(min(self.resid_norms) - max(self.gamma_norms)*0.1, min(self.resid_norms) +  max(self.gamma_norms) * 1.1)
            ax.set_ylim(0, max(self.gamma_norms) * 1.1)
        return fig, ax

    def get_optimal_alpha(self, n_alpha=50, low_r_bound=1, up_r_bound=800):
        self.alpha_vals = np.logspace(-5, 2, n_alpha)
        self.gamma_norms = []
        self.resid_norms = []
        self.resid_v = []

        for alpha in self.alpha_vals:
            self.fit_g1(alpha, low_r_bound, up_r_bound)
            self.gamma_norms.append(np.linalg.norm(np.array(self.G_gamma))**2)
            self.resid_norms.append(np.linalg.norm(np.array(self.res_lsq))**2)
            self.resid_v.append(np.linalg.norm(np.array(self.res_lsq))**2 +\
                                alpha**2 * np.linalg.norm(np.array(self.res_reg))**2)

        n = 3
        line_1 = stats.linregress(x = self.resid_norms[-n:], y = self.gamma_norms[-n:])
        line_2 = stats.linregress(x = self.resid_norms[:n],  y = self.gamma_norms[:n])
        x_pt = (line_2.intercept - line_1.intercept)/(line_1.slope - line_2.slope)
        y_pt = line_1.slope * x_pt + line_1.intercept
        yscale = self.gamma_norms[0] - self.gamma_norms[-1]
        xscale = self.resid_norms[-1] - self.resid_norms[0]
        self.distances = [np.sqrt((yscale/xscale*(x_pt - x1))**2 + (y_pt - y1)**2) for (x1, y1) in zip(self.resid_norms, self.gamma_norms)]
        # self.distances = [np.sqrt(((x_pt - x1))**2 + (y_pt - y1)**2) for (x1, y1) in zip(self.resid_norms, self.gamma_norms)]
        self.index_knee = self.distances.index(min(self.distances))
        self.alpha_opt = self.alpha_vals[self.index_knee]
        self.fit_g1(self.alpha_opt, low_r_bound, up_r_bound)
        self.get_rh_dist()

        # For plotting
        self.line_1 = line_1
        self.line_2 = line_2
        self.x_pt, self.y_pt = x_pt, y_pt
        return











# Previous as of 2022-07-22
# def get_optimal_alpha(self):
#     self.alpha_vals = np.logspace(-5, 2, 50)
#     self.gamma_norms = []
#     self.resid_norms = []
#     self.resid_v = []

#     for alpha in self.alpha_vals:
#         self.fit_g1(alpha)
#         self.gamma_norms.append(np.linalg.norm(np.array(self.G_gamma)))
#         self.resid_norms.append(np.linalg.norm(np.array(self.res_lsq)))
#         self.resid_v.append(np.linalg.norm(np.array(self.res_lsq))**2 +\
#                             alpha**2 * np.linalg.norm(np.array(self.res_reg))**2)

#     spline = interpolate.interp1d(self.resid_norms, self.gamma_norms,
#                                   kind='quadratic', fill_value="extrapolate")
#     self.resid_norm_pts = np.linspace(min(self.resid_norms),
#                                  max(self.resid_norms), 1000)
#     self.gamma_norm_pts = spline(self.resid_norm_pts)

#     y_x = misc.derivative(spline, self.resid_norm_pts,
#                           dx=self.resid_norm_pts[2]-self.resid_norm_pts[0], n=1)
#     y_xx = misc.derivative(spline, self.resid_norm_pts,
#                           dx=self.resid_norm_pts[2]-self.resid_norm_pts[0], n=2)
#     self.k = y_xx/((1 + y_x**2)**(3/2))

#     self.x_knee = self.resid_norm_pts[np.argmax(self.k)]
#     self.y_knee = spline(self.x_knee)

#     def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):
#         distance = (y_array-y_point)**2 + (x_array-x_point)**2
#         index = np.where(distance==distance.min())
#         return index[0][0]

#     self.index_knee = find_index_of_nearest_xy(self.gamma_norms,
#                                                self.resid_norms,
#                                                self.y_knee, self.x_knee)

#     self.alpha_opt = self.alpha_vals[self.index_knee]
#     self.fit_g1(self.alpha_opt)
#     self.get_rh_dist()
#     return


# def draw_L_curve(self, loglog=False):
#     fig, ax = my_plot.instantiate_fig(ylabel=r'$\Vert  x  \Vert$',
#                                      xlabel=r'$\Vert  Ax - g \Vert$')
#     ax.scatter(self.resid_norms, self.gamma_norms)
#     ax.plot(self.resid_norm_pts, self.gamma_norm_pts, 'k--')
#     ax.scatter(self.x_knee, self.y_knee)
#     ax.scatter(self.resid_norms[self.index_knee],
#                self.gamma_norms[self.index_knee],
#                s=200, facecolors='none', edgecolors='r')
#     if loglog:
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#     my_plot.set_layout(fig, ax)
#     return fig, ax


# def draw_L_curve_power(self, n=2, loglog=False):
#     fig, ax = my_plot.instantiate_fig(ylabel=fr'$\Vert  x  \Vert ^{n}$',
#                                      xlabel=fr'$\Vert  Ax - g \Vert ^{n}$')
#     ax.scatter(np.array(self.resid_norms)**n, np.array(self.gamma_norms)**n)
#     ax.plot(np.array(self.resid_norm_pts)**n, np.array(self.gamma_norm_pts)**n, 'k--')
#     ax.scatter(np.array(self.x_knee)**n, np.array(self.y_knee)**n)
#     ax.scatter(np.array(self.resid_norms[self.index_knee])**n,
#                np.array(self.gamma_norms[self.index_knee])**n,
#                s=200, facecolors='none', edgecolors='r')
#     if loglog:
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#     my_plot.set_layout(fig, ax)
#     return fig, ax

# def get_optimal_alpha_curvature(self):
#     self.alpha_vals = np.logspace(-5, 2, 50)
#     self.gamma_norms = []
#     self.resid_norms = []
#     self.resid_v = []

#     for alpha in self.alpha_vals:
#         self.fit_g1(alpha)
#         self.gamma_norms.append(np.linalg.norm(np.array(self.G_gamma)))
#         self.resid_norms.append(np.linalg.norm(np.array(self.res_lsq)))
#         self.resid_v.append(np.linalg.norm(np.array(self.res_lsq))**2 +\
#                             alpha**2 * np.linalg.norm(np.array(self.res_reg))**2)

#     spline = interpolate.interp1d(self.resid_norms, self.gamma_norms,
#                                   kind='quadratic', fill_value="extrapolate")
#     self.resid_norm_pts = np.linspace(min(self.resid_norms),
#                                  max(self.resid_norms), 1000)
#     self.gamma_norm_pts = spline(self.resid_norm_pts)

#     y_x = misc.derivative(spline, self.resid_norm_pts,
#                           dx=self.resid_norm_pts[2]-self.resid_norm_pts[0], n=1)
#     y_xx = misc.derivative(spline, self.resid_norm_pts,
#                           dx=self.resid_norm_pts[2]-self.resid_norm_pts[0], n=2)
#     self.k = y_xx/((1 + y_x**2)**(3/2))

#     self.x_knee = self.resid_norm_pts[np.argmax(self.k)]
#     self.y_knee = spline(self.x_knee)

#     def find_index_of_nearest_xy(y_array, x_array, y_point, x_point):
#         distance = (y_array-y_point)**2 + (x_array-x_point)**2
#         index = np.where(distance==distance.min())
#         return index[0][0]

#     self.index_knee = find_index_of_nearest_xy(self.gamma_norms,
#                                                self.resid_norms,
#                                                self.y_knee, self.x_knee)

#     self.alpha_opt = self.alpha_vals[self.index_knee]
#     self.fit_g1(self.alpha_opt)
#     self.get_rh_dist()
#     return
