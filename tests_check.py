"""Quick self-test: python tests_check.py  (expects ~18.3 m2/g for the sample file)."""
import n2o_analysis as na
raw = na.load_file(open("sample_data/80Cu20ZrO2_N2O_21.10.26.dat", "rb").read(), "x.dat")
cols = na.guess_columns(raw)
p = na.prepare(raw, cols)
res, _ = na.analyse(p, na.Sample(), na.auto_windows(p))
print(f"S Cu = {res.s_per_g_ox:.2f} m2/g cat  (expected 18.3)")
assert abs(res.s_per_g_ox - 18.33) < 0.2
