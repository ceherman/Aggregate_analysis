#!/usr/bin/env pythonw
import os
import sys
import re
import argparse
import tempfile
import numpy as np
import subprocess as sp
from scipy.optimize import leastsq, nnls

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx as Toolbar
import wx
import wx.aui
import matplotlib as mpl
import matplotlib.pyplot as plt

# 図の体裁
plt.style.use('classic')
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
plt_dic['axes.grid'] = True
plt_dic['font.size'] = 12
plt_dic['legend.fontsize'] = 12
plt_dic['axes.labelsize'] = 14
plt_dic['xtick.direction'] = 'in'
plt_dic['savefig.bbox'] = 'tight'
plt_dic['savefig.dpi'] = 150
plt.rcParams.update(plt_dic)

output_kugiri = "".join(["-"] * 72) + "\n"
description = "1D exponential decay analysis program by multi-components, \n"
description += "descrnon-negtive least square, contin, and BRD methods."


par = argparse.ArgumentParser(description=description)
par.add_argument('-f', '--files', metavar='decay.int',
                 nargs="+", required=True, help='1D decay data')
par.add_argument('-x', '--minimum', dest="min", default=1e-4, type=float,
                 help='mininum value for inverse space [1e-4]')
par.add_argument('-y', '--maximum', dest="max", default=1e+1, type=float,
                 help='miximum value for inverse space [1e+1]')
par.add_argument('-z', '--gridnum', dest="ng", default=400, type=int,
                 help='number of poins for inverse space [400]')
par.add_argument('-g', '--grid', default="log", choices=['linear', 'log'],
                 help='griding rule for inverse space [log]')
par.add_argument('-s', '--start', default=0, type=float,
                 help='starting time value for decay data [0.0]')
par.add_argument('-e', '--end', default=-1, type=float,
                 help='ending time for decay data [last time point]')
par.add_argument('-t', '--type', default="cpmg",
                 choices=['cpmg', 'ir', 't2g', 'diff'],
                 help='decay data type [cpmg]')
par.add_argument('-m', '--methods',
                 default=['brd', 'contin', 'nnls', 'discrete'],
                 nargs="+", choices=['brd', 'contin', 'nnls', 'discrete'],
                 help='analyzing mehods [discrete nnls brd contin]')
par.add_argument('-b', '--baseline', default=False,
                 action='store_true', help='adding constant base line')

par.add_argument('-n', '--ncomp', default=2, type=int,
                 help='number of components for discrete analysis [2]')
# contin specific
par.add_argument('-iwt', '--weight', default=1, type=int,
                 choices=[1, 2, 3],
                 help='error weighing routine on contin [1]')
par.add_argument('-inf', '--infinity', default=3, type=int,
                 help='number of last points for transform function form from ir to cpmg [1]')
# brd specific
par.add_argument('-atype', '--alpha_type', default=1, type=int,
                 choices=[1, 2, 3],
                 help='alpha determination algorism for BRD analysis [1]')
par.add_argument('-alpha', '--alpha_initial', default=1, type=float,
                 help='initial alpha for BRD analysis [1]')
par.add_argument('-alpha_loop', '--alpha_loop', default=5000, type=int,
                 help='maximu alpha loop [5000]')
par.add_argument('-fs', '--figsize', nargs=2, default=[800, 600], type=float,
                 help='figure size x and y [800 600]')

# boolean
par.add_argument('-np', '--noplot', default=False,
                 action='store_true', help='no-plot figure')
par.add_argument('-v', '--verbose', default=False,
                 action='store_true', help='verbose mode')

# file
par.add_argument('-o', '--outfile', nargs='?', default="",
                 help='specific output file name')

args = par.parse_args()

if args.grid == "log":
    args.w = np.logspace(np.log10(args.min), np.log10(args.max), args.ng)
if args.grid == "linear":
    args.w = np.linspace(args.min, args.max, args.ng)


class GetData():
    def __init__(self, infile, num):
        def nearest(v): return np.argmin(np.abs(self.x - v))
        data = np.loadtxt(infile)
        self.x = data[:, 0]
        self.y = data[:, 1]
        if args.end == -1:
            args.end = self.x[-1]
        s = nearest(args.start)
        e = nearest(args.end) + 1
        self.xs = self.x[s:e]
        self.ys = self.y[s:e]
        base, ext = os.path.splitext(infile)
        outfile = base + ".out"
        if args.outfile == "":
            self.out = open(outfile, "w")
        elif args.outfile == "-":
            self.out = sys.stdout
            outfile = "stdout"
        else:
            if len(args.files) == 1:
                outfile = args.outfile
            else:
                outfile = "%d" % num + "_" + args.outfile
            self.out = open(outfile, "w")
        self.infile = infile
        self.outfile = outfile


def DistrPlot(ax, data):
    # Contin Distr
    if args.grid == "log":
        ax.semilogx(args.w, data.f[:args.ng])
    else:
        ax.plot(args.w, data.f[:args.ng])
    ax.set_title(data.basename)
    ax.set_xlim(args.min, args.max)
    if args.type == "cpmg":
        xlabel = "T2/s"
    if args.type == "ir":
        xlabel = "T1/s"
    if args.type == "T2g":
        xlabel = "T2/s"
    if args.type == "diff":
        xlabel = "D/e-9m2/s"

    ax.set_xlabel("%s" % (xlabel))
    ax.set_ylabel("Intensity")
    text_str = "||Iexp - Isim||/N=%g" % (data.chi)
    ax.annotate(text_str, xy=(0.02, 0.98), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top')


def DecayPlot(ax, data):
    text_str = "||Iexp - Isim||/N=%g" % (data.chi)
    ax.set_xlabel("Time")
    ax.set_ylabel("Intensity")
    ax.set_title(data.basename)
    ax.annotate(text_str, xy=(0.02, 0.98), xycoords='axes fraction',
                horizontalalignment='left', verticalalignment='top')
    ax.plot(data.expx, data.expy, ".", label="Exp.")
    ax.plot(data.calx, data.caly, "-", label="Calc.")
    if args.baseline:
        xb = [data.calx[0], data.calx[-1]]
        yb = [data.baseline, data.baseline]
        ax.plot(xb, yb, color="gray", label="baseline")

    axe = ax.twinx()
    axe.plot(data.calx, data.dify, "r.", label="difference")
    for tl in axe.get_yticklabels():
        tl.set_color('red')
    axe.set_ylabel('Difference', color="red")
    axe.grid()
    ax.legend(loc="best")


def Output(f, data, string):
    f.write(output_kugiri)
    f.write("%s analysis\n" % (string))
    f.write(output_kugiri)
    f.write("baseline = %g\n" % (data.baseline))
    f.write("alpha: %s\n" % (data.alpha))
    f.write("||Ical-Iexp||/N: %g\n" % (data.chi))
    f.write("\n# Decay data (%s)\n" % (string))
    f.write("%15s" % ("time"))
    f.write("%15s" % ("exp."))
    f.write("%15s" % ("cal."))
    f.write("%15s" % ("diff."))
    f.write("\n")
    for i, x in enumerate(data.calx):
        f.write("%15e" % (x))
        f.write("%15e" % (data.expys[i]))
        f.write("%15e" % (data.caly[i]))
        f.write("%15e" % (data.dify[i]))
        f.write("\n")
    f.write("\n# Distribution data (%s)\n" % (string))
    f.write("%15s" % ("w"))
    f.write("%15s" % ("intensity"))
    f.write("\n")
    for i, w in enumerate(args.w):
        f.write("%15e" % (w))
        f.write("%15e" % (data.f[i]))
        f.write("\n")
    f.write("\n")


class DISCRETE():
    def __init__(self, data):
        print("Analyzing DISCRETE...")
        sys.stdout.flush()
        self.expx = data.x
        self.expy = data.y
        self.expys = data.ys
        self.baseline = 0.0
        self.alpha = 0.0
        p0 = self.SetInit(data.xs, data.ys)
        if args.baseline:
            p0 = np.append(p0, self.baseline)
        self.Fit(p0, data.xs, data.ys)
        self.chi = np.sum(self.dify**2 / len(self.dify))
        self.Output(data.out)
        if args.noplot is False:
            self.Plot(plotter)

    def SetInit(self, x, y):
        yv = np.abs(y)
        ymax = np.max(yv)
        p = np.argmin(np.abs(yv - 0.1 * np.max(yv)))
        xmax = x[p]
        t0 = np.linspace(xmax / args.ncomp, xmax, args.ncomp)
        m0 = [ymax / args.ncomp] * args.ncomp
        p0 = np.append(m0, t0)
        if args.baseline:
            p0 = np.append(p0, np.average(yv[-5:]))
        return p0

    def Fit(self, p0, x, y):
        def cpmg_func(p, ncomp):
            k = np.exp(np.outer(-x, 1.0 / p[ncomp:2 * ncomp]))
            yi = np.inner(k, p[:ncomp])
            return yi

        def cpmg_error(p, ncomp):
            delta = y - cpmg_func(p, ncomp)
            if args.baseline:
                delta = delta - p[-1]
            if np.any(p < 0):
                return np.ones_like(y) * 1e10
            return delta

        def diff_func(p, ncomp):
            k = np.exp(np.outer(-x, 1.0 * p[ncomp:2 * ncomp]))
            yi = np.inner(k, p[:ncomp])
            return yi

        def diff_error(p, ncomp):
            delta = y - diff_func(p, ncomp)
            if args.baseline:
                delta = delta - p[-1]
            if np.any(p < 0):
                return np.ones_like(y) * 1e10
            return delta

        def t2g_func(p, ncomp):
            k = np.exp(np.outer(-0.5 * x**2, (1.0 / p[ncomp:2 * ncomp])**2))
            yi = np.inner(k, p[:ncomp])
            return yi

        def t2g_error(p, ncomp):
            delta = y - t2g_func(p, ncomp)
            if args.baseline:
                delta = delta - p[-1]
            if np.any(p < 0):
                return np.ones_like(y) * 1e10
            return delta

        def ir_func(p, ncomp):
            k = 1.0 - 2.0 * np.exp(np.outer(-x, 1.0 / p[ncomp:2 * ncomp]))
            yi = np.inner(k, p[:ncomp])
            return yi

        def ir_error(p, ncomp):
            delta = y - ir_func(p, ncomp)
            if args.baseline:
                delta = delta - p[-1]
            if np.any(p < 0):
                return np.ones_like(y) * 1e10
            return delta
        if args.type == "t2g":
            result = leastsq(t2g_error, p0,
                             args=(args.ncomp), full_output=True)
            res = result[0]
            residual = t2g_error(res, args.ncomp)
            self.caly = t2g_func(res, args.ncomp)
        if args.type == "cpmg":
            result = leastsq(cpmg_error, p0,
                             args=(args.ncomp), full_output=True)
            res = result[0]
            residual = cpmg_error(res, args.ncomp)
            self.caly = cpmg_func(res, args.ncomp)
        if args.type == "diff":
            result = leastsq(diff_error, p0,
                             args=(args.ncomp), full_output=True)
            res = result[0]
            residual = diff_error(res, args.ncomp)
            self.caly = diff_func(res, args.ncomp)
        if args.type == "ir":
            result = leastsq(ir_error, p0,
                             args=(args.ncomp), full_output=True)
            res = result[0]
            residual = ir_error(res, args.ncomp)
            self.caly = ir_func(res, args.ncomp)
        rcs = np.sum(residual**2) / (len(x) - args.ncomp * 2)
        try:
            cover = np.diagonal(result[1] * rcs)  # rcsを掛けて共分散を求める
        except:
            cover = np.zeros(2 * args.ncomp)
        if args.verbose:
            print(result)
        # m1, m1+/-err, m1err/%, T1, T1+/-err, T1err/%,
        # m2, m2+/-err, m2err/%, T2, T2+/-err, T2err/%,
        # m3, m3+/-err, m3err/%, T3, T3+/-err, T3err/%,
        ans = np.zeros((args.ncomp, 6))
        ans[:, 0] = res[:args.ncomp]
        ans[:, 1] = cover[:args.ncomp]
        ans[:, 2] = cover[:args.ncomp] / ans[:, 0] * 100
        ans[:, 3] = res[args.ncomp:2 * args.ncomp]
        ans[:, 4] = cover[args.ncomp:2 * args.ncomp]
        ans[:, 5] = cover[args.ncomp:2 * args.ncomp] / ans[:, 3] * 100
        self.calx = x
        self.dify = y - self.caly
        self.ans = ans

        # 各成分の計算
        self.calyi = np.zeros((len(self.calx), args.ncomp))
        if args.baseline:
            self.baseline = res[-1]
            self.caly = self.caly + self.baseline
        for i in range(args.ncomp):
            param = np.array([ans[i, 0], ans[i, 3]])
            if args.type == "cpmg":
                self.calyi[:, i] = cpmg_func(param, 1) + self.baseline
            if args.type == "diff":
                self.calyi[:, i] = diff_func(param, 1) + self.baseline
            if args.type == "t2g":
                self.calyi[:, i] = t2g_func(param, 1) + self.baseline
            if args.type == "ir":
                self.calyi[:, i] = ir_func(param, 1) + self.baseline

    def Plot(self, plotter):
        ax = plotter.add('Discrete decay').gca()
        text_str = "||Iexp - Isim||/N=%g" % (self.chi)
        ax.annotate(text_str, xy=(0.02, 0.98), xycoords='axes fraction',
                    horizontalalignment='left', verticalalignment='top')
        ax.plot(self.expx, self.expy, ".", label="Delete data")
        ax.set_xlabel("Time")
        ax.set_ylabel("Intensity")
        Mrate = self.ans[:, 0] / np.sum(self.ans[:, 0]) * 100
        for i in range(args.ncomp):
            T = self.ans[i, 3]
            Mr = Mrate[i]
            ax.plot(self.calx, self.calyi[:, i], "-",
                    label="T=%g (%g%%)" % (T, Mr))
        ax.plot(self.calx, self.caly, label="Calc.")
        if args.baseline:
            xb = [self.calx[0], self.calx[-1]]
            yb = [self.baseline, self.baseline]
            ax.plot(xb, yb, color="gray", label="baseline")
        ax.legend(loc="best")
        axe = ax.twinx()
        axe.plot(self.calx, self.dify, "r.", label="difference")
        for tl in axe.get_yticklabels():
            tl.set_color('red')
        axe.set_ylabel('Difference', color="red")
        axe.grid()

    def Output(self, f):
        f.write(output_kugiri)
        f.write("Discrete analysis\n")
        f.write(output_kugiri)
        f.write("Fitting x range: %-g-%g\n" % (self.calx[0], self.calx[-1]))
        f.write("%d components anaysis\n" % (args.ncomp))
        for i in range(args.ncomp):
            M, Merr, Mper = self.ans[i, 0], self.ans[i, 1], self.ans[i, 2]
            T, Terr, Tper = self.ans[i, 3], self.ans[i, 4], self.ans[i, 5]
            Mp = M / np.sum(self.ans[:, 0]) * 100
            f.write("### %d component: population (%g%%)\n" % (i + 1, Mp))
            f.write("T: %g +/- %g (%g%%)\n" % (T, Terr, Tper))
            f.write("M: %g +/- %g (%g%%)\n" % (M, Merr, Mper))
        f.write("\n# Decay data (DISCRETE) %d componets\n" % (args.ncomp))
        f.write("%15s" % ("time"))
        f.write("%15s" % ("exp."))
        for j in range(args.ncomp):
            f.write("%16s" % "cal.%d" % (j + 1))
        f.write("%15s" % ("cal."))
        f.write("%15s" % ("diff."))
        f.write("\n")

        for i, x in enumerate(self.calx):
            f.write("%15e" % (x))
            f.write("%15e" % (self.expys[i]))
            for j in range(args.ncomp):
                f.write("%15e" % (self.calyi[i, j]))
            f.write("%15e" % (self.caly[i]))
            f.write("%15e" % (self.dify[i]))
            f.write("\n")
        f.write("\n")


class NNLS():
    def __init__(self, data):
        print("Analyzing NNLS...")
        sys.stdout.flush()
        self.expx = data.x
        self.expy = data.y
        self.expxs = data.xs
        self.expys = data.ys
        self.baseline = 0.0
        self.alpha = 0.0
        self.Fit(data.xs, data.ys)
        self.chi = np.sum(self.dify**2) / len(self.dify)
        Output(data.out, self, "NNLS")
        self.basename = os.path.basename(data.infile)
        if args.noplot == False:
            self.Plot(plotter)

    def Fit(self, x, y):
        if args.type == "cpmg":
            k = np.exp(np.outer(-x, 1.0 / args.w))
        if args.type == "diff":
            k = np.exp(np.outer(-x, 1.0 * args.w))
        if args.type == "t2g":
            k = np.exp(np.outer(-0.5 * x**2, (1.0 / args.w)**2))
        if args.type == "ir":
            k = 1.0 - 2.0 * np.exp(np.outer(-x, 1.0 / args.w))
        if args.baseline:
            k = np.append(k, np.ones((len(x), 1)), axis=1)
        self.calx = x
        self.f = np.zeros(len(k[0]))
        try:
            self.f, verbose = nnls(k, y)
            self.caly = np.inner(k, self.f)
            # print " NNLS (alpha=0) was successful."
        except (RuntimeError):
            self.f = np.zeros(len(k[0]))
            self.caly = np.zeros(len(x))
            print(" NNLS (alpha=0) was fault!")
            print(" problem can not be solved.")
        self.dify = y - self.caly
        if args.baseline:
            self.baseline = self.f[-1]
        if args.verbose:
            print(verbose)

    def Plot(self, plotter):
        ax1 = plotter.add('NNLS distr').gca()
        ax2 = plotter.add('NNLS decay').gca()
        DistrPlot(ax1, self)
        DecayPlot(ax2, self)


class BRD():
    def __init__(self, data):
        print("Analyzing BRD...")
        sys.stdout.flush()
        self.expx = data.x
        self.expy = data.y
        self.expys = data.ys
        self.baseline = 0.0
        self.alpha = 0.0
        self.Fit(data.xs, data.ys)
        self.chi = np.sum(self.dify**2) / len(self.dify)
        Output(data.out, self, "BRD")
        self.basename = os.path.basename(data.infile)
        if args.noplot is False:
            self.Plot(plotter)

    def Fit(self, x, y):
        fp = [tempfile.NamedTemporaryFile('w', delete=True)]
        if args.verbose:
            fp.append(sys.stdout)

        if args.type == "cpmg":
            k = np.exp(np.outer(-x, 1.0 / args.w))
            dtype = "T2"
        if args.type == "diff":
            k = np.exp(np.outer(-x, args.w))
            dtype = "D"
        if args.type == "ir":
            k = 1.0 - 2.0 * np.exp(np.outer(-x, 1.0 / args.w))
            dtype = "T1"
        if args.type == "t2g":
            k = np.exp(np.outer(-0.5 * x**2, (1.0 / args.w)**2))
            dtype = "T2g"
        for f in fp:
            for i, t in enumerate(x):
                s = "%g %g\n" % (t, y[i])
                f.write(s)
            f.flush()
        if args.baseline:
            dtype = dtype + "inf"
        prog = "osilap"
        com = [prog,
               "--dimension=1",
               "--output=-",
               "--sfactor=%g" % (1e-8),
               "--alpha=%g" % (args.alpha_initial),
               "--alpha_loop_max=%d" % (args.alpha_loop),
               "--Atype=%d" % (args.alpha_type),
               "--TXtype=%s" % (dtype),
               "--TXnum=%d" % (args.ng),
               "--TXmin=%f" % (args.min),
               "--TXmax=%f" % (args.max),
               fp[0].name]
        if args.grid == "log":
            com.append("--TXspace=log")
        if args.grid == "linear":
            com.append("--TXspace=linear")
        pp = sp.run(com, stdout=sp.PIPE, stderr=sp.STDOUT)
        output = pp.stdout.decode('utf-8')
        cstr = "#  Alpha :\s+(\S+)"
        obj = re.search(cstr, output, re.DOTALL)
        alpha = float(obj.group(1))
        self.alpha = float(obj.group(1))
        cstr = "# Distribution normalized coeff.:(.+?)\n#.+?\n(.+?)\n\n"
        obj = re.search(cstr, str(output), re.DOTALL)
        if args.verbose:
            print(output)
        norm = float(obj.group(1))
        data = np.array(obj.group(2).split(), dtype=float).reshape(-1, 2)
        self.f = data[:, 1] * norm
        self.calx = x
        self.caly = np.inner(k, self.f)
        self.dify = y - self.caly
        if args.baseline:
            obj = re.search("# Baseline list\n#(.+)\n", output)
            self.baseline = float(obj.group(1)) * norm
            self.caly = self.caly + self.baseline
        if args.verbose:
            output

    def Plot(self, plotter):
        ax1 = plotter.add('BRD distr').gca()
        ax2 = plotter.add('BRD decay').gca()
        # BRD plot
        DistrPlot(ax1, self)
        DecayPlot(ax2, self)


class CONTIN():
    def __init__(self, data):
        print("Analyzing CONTIN...")
        if args.type == "t2g":
            print("CONTIN can not analyze T2 gaussian function.")
            return
        sys.stdout.flush()
        self.expx = data.x
        self.expy = data.y
        self.expys = data.ys
        self.baseline = 0.0
        self.alpha = 0
        self.basename = os.path.basename(data.infile)
        # continの関数はexp(-t/T2)しかできないので式変形
        if args.type == "ir":
            ave = np.average(self.expy[-args.infinity])
            y_corr = ave - data.ys
        else:
            y_corr = data.ys
        self.Fit(data.xs, y_corr, data.ys)
        self.chi = np.sum(self.dify**2) / len(self.dify)
        Output(data.out, self, "CONTIN")
        if args.noplot is False:
            self.Plot(plotter)

    def Fit(self, x, y, yorig):
        pp = sp.Popen(['contin'], stdin=sp.PIPE, stdout=sp.PIPE,
                      stderr=sp.PIPE, universal_newlines=True, bufsize=1)
        # o = open("hoge.cin", "w")
        fp = [pp.stdin]
        if args.verbose:
            fp.append(sys.stdout)
        if args.grid == "linear":
            igrid = 1  # linear
        if args.grid == "log":
            igrid = 2    # log
        if args.baseline:
            nlinf = 1
        else:
            nlinf = 0
        if args.type == "diff":
            R22 = 1
        else:
            R22 = -1
        # 出力
        # ヘッダの出力
        for f in fp:
            f.write(" %-80s\n" % ("created by decay2distr.py"))
            f.write(" %-6s%5s%15e\n" % ("LAST",   "", 1))
            f.write(" %-6s%5s%15e\n" % ("NG",     "", args.ng))
            f.write(" %-6s%5d%15e\n" % ("GMNMX",  1, args.max))
            f.write(" %-6s%5d%15e\n" % ("GMNMX",  2, args.min))
            f.write(" %-6s%5s%15e\n" % ("IWT",    "", args.weight))
            f.write(" %-6s%5s%15e\n" % ("NERFIT", "", 10))
            f.write(" %-6s%5s%15e\n" % ("NINTT",  "", 0))
            f.write(" %-6s%5s%15e\n" % ("NLINF",  "", nlinf))
            f.write(" %-6s%5s%15e\n" % ("IQUAD",  "", 1.0))
            f.write(" %-6s%5s%15e\n" % ("IGRID",  "", igrid))
            f.write(" %-6s%5s%15s\n" % ("IFORMT", "", ""))
            f.write(" %-6s%5s%15s\n" % ("(1E15)", "", ""))
            f.write(" %-6s%5s%15s\n" % ("IFORMY", "", ""))
            f.write(" %-6s%5s%15s\n" % ("(1E15)", "", ""))
            f.write(" %-6s%5s%15e\n" % ("DOUSNQ", "", 1))
            f.write(" %-6s%5d%15e\n" % ("IUSER",  3, -1))
            f.write(" %-6s%5d%15e\n" % ("IUSER",  10, 4))
            f.write(" %-6s%5d%15e\n" % ("RUSER",  21, 1))
            f.write(" %-6s%5d%15e\n" % ("RUSER",  22, R22))
            f.write(" %-6s%5d%15e\n" % ("RUSER",  23, 0))
            f.write(" %-6s%5d%15e\n" % ("LUSER",  3, -1))
            f.write(" %-6s%5s%15s\n" % ("END", "", ""))
            # データの出力
            f.write(" %-6s%5d%15s\n" % ("NY", len(x), ""))
            for xi in x:
                f.write("%15e\n" % (xi))
            for yi in y:
                f.write("%15e\n" % (yi))
        for f in fp:
            f.close()
        output = pp.stdout.read()
        self.f = np.zeros(args.ng)
        self.f_error = np.zeros(args.ng)
        # 標準出力からデータを取得
        distrs = re.findall("##### DISTR START #####(.+?)##### END #####",
                            output, re.DOTALL)
        distr = distrs[-1].strip().split("\n")
        if args.verbose:
            print(output)
        for i, line in enumerate(distr):
            data = [float(d) for d in line.split()[:3]]
            self.f[i] = data[0]
            self.f_error[i] = data[1]
        if args.baseline:
            bases = re.findall("LINEAR COEFFICIENTS = (.+)", output)
            self.baseline = float(bases[-1].split()[0])
        # 減衰データの計算
        self.f = self.f[::-1]
        if args.type == "cpmg":
            k = np.exp(np.outer(-x, 1.0 / args.w))
        if args.type == "ir":
            k = 0.5 - 1.0 * np.exp(np.outer(-x, 1.0 / args.w))
        if args.type == "diff":
            k = np.exp(np.outer(-x, args.w))
        self.calx = x
        self.caly = np.inner(k, self.f) + self.baseline
        self.dify = yorig - self.caly

    def Plot(self, plotter):
        ax1 = plotter.add('Contin distr').gca()
        ax2 = plotter.add('Contin decay').gca()
        # Contin Distr
        DistrPlot(ax1, self)
        DecayPlot(ax2, self)


def HeaderOutput(data):
    f = data.out
    f.write(output_kugiri)
    f.write("General information\n")
    f.write(output_kugiri)
    f.write("Input file: %s\n" % (data.infile))
    f.write("Output file: %s\n" % (data.outfile))
    f.write("Selected analyzing methods: %s\n" % (args.methods))
    f.write("Exp. x range: %-g-%g\n" % (data.x[0], data.x[-1]))
    f.write("Minimum T  : %g\n" % (data.x[0]))
    f.write("Maximum T  : %g\n" % (data.x[-1]))
    f.write("Exp.  points: %d\n" % (len(data.x)))
    f.write("Fitting x range: %-g-%g\n" % (data.xs[0], data.xs[-1]))
    f.write("Minimum T  : %g\n" % (data.xs[0]))
    f.write("Maximum T  : %g\n" % (data.xs[-1]))
    f.write("Fitting  points: %d\n" % (len(data.xs)))
    f.write("Grid points: %d\n" % (args.ng))
    f.write("Grid spaces: %s\n" % (args.grid))

    f.write("\n# Decay data (EXPERIMENTAL)\n")
    f.write("%15s" % ("time"))
    f.write("%15s" % ("intensity"))
    f.write("\n")
    for i, x in enumerate(data.x):
        f.write("%15e" % (x))
        f.write("%15e" % (data.y[i]))
        f.write("\n")
    f.write("\n")
    print("Input file : %s" % (data.infile))
    print("Output file: %s" % (data.outfile))
    print(output_kugiri, end='')
    sys.stdout.flush()


# プロット
if args.noplot is False:
    # タブへのプロット

    class Plot(wx.Panel):
        def __init__(self, parent, status, id=-1, dpi=None, **kwargs):
            wx.Panel.__init__(self, parent, id=id, **kwargs)
            self.statusbar = status
            sizer = wx.BoxSizer(wx.VERTICAL)
            self.figure = mpl.figure.Figure(dpi=dpi)
            self.canvas = Canvas(self, -1, self.figure)
            self.canvas.mpl_connect('motion_notify_event',
                                    self.UpdateStatusBar)
            sizer.Add(self.canvas, 1, wx.EXPAND)
            self.toolbar = Toolbar(self.canvas)
            self.toolbar.Realize()
            sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
            self.SetSizer(sizer)
            self.text = ""

        def UpdateStatusBar(self, event):
            if event.inaxes:
                x, y = event.xdata, event.ydata
                self.text = "x=%g, y=%g" % (x, y)
            self.statusbar.SetStatusText(self.text, 0)

    class PlotNotebook(wx.Panel):
        def __init__(self, parent, status, id=-1):
            wx.Panel.__init__(self, parent, id=id)
            self.nb = wx.aui.AuiNotebook(self)
            self.status = status
            sizer = wx.BoxSizer(wx.VERTICAL)
            sizer.Add(self.nb, 1, wx.EXPAND)
            self.SetSizer(sizer)

        def add(self, name="plot"):
            page = Plot(self.nb, self.status)
            self.nb.AddPage(page, name)
            return page.figure
    app = wx.App(False)
    frame = wx.Frame(None, -1, 'Plotter', size=args.figsize)
    status = frame.CreateStatusBar()
    plotter = PlotNotebook(frame, status)


# データの読み込み
for i, f in enumerate(args.files):
    if os.path.exists(f) is False:
        print("##### %s is not found. %s is skipped. #####\n" % (f, f))
        sys.stdout.flush()
        continue
    data = GetData(f, i)
    HeaderOutput(data)
    for m in args.methods:
        if m == "nnls":
            nnls_ = NNLS(data)
        if m == "contin":
            contin = CONTIN(data)
        if m == "brd":
            brd = BRD(data)
        if m == "discrete":
            discrete = DISCRETE(data)
    if data.outfile != "stdout":
        data.out.close()
    print("\n", output_kugiri)
    sys.stdout.flush()
if args.noplot is False:
    frame.Show()
    app.MainLoop()
