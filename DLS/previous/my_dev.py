# Discarded stuff from Jupyter

from tabulate import tabulate

def to_fwf(df, fname):
    content = tabulate(df.values.tolist(), list(df.columns), tablefmt="plain")
    open(fname, "w").write(content)

pd.DataFrame.to_fwf = to_fwf

part[['lag_time_s', 'g1']].to_fwf('./test.int')

# ./osilap test.int --dimension 1 --TXmin 200 --TXmax 100000 --TXtype T2 --output test_T2

from mpmath import calculus
dir(calculus.inverselaplace.InverseLaplaceTransform)
calculus.inverselaplace.InverseLaplaceTransform.calc_time_domain_solution?
