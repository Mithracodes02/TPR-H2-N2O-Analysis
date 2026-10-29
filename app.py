"""Streamlit app: metallic Cu surface area from N2O titration (TPR-H2 followed by N2O pulse/step), MS data.

Run locally:  streamlit run app.py
"""
import io

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import n2o_analysis as na

st.set_page_config(page_title="N2O Cu surface area", page_icon="⚗️", layout="wide")

st.title("⚗️ Metallic Cu surface area – N2O titration")
st.caption(
    "Upload one or more mass-spec exports (PV MassSpec .dat, or CSV/XLSX with time + signals). "
    "Windows are detected automatically; adjust them with the sliders if needed."
)

# ----------------------------------------------------------------------------- sidebar: common settings
with st.sidebar:
    st.header("Experiment (common to all files)")
    y_n2o = st.number_input("N2O in feed (%)", value=2.0, min_value=0.01, step=0.1) / 100
    flow = st.number_input("Total flow (mL/min)", value=50.0, min_value=0.1, step=1.0)
    vm = st.number_input("Molar volume (mL/mol)", value=22400.0, step=1.0,
                         help="22400 = 0 °C, 1 atm. Use 24465 if the flow is referenced to 25 °C.")
    site = st.number_input("Cu surface atoms (atoms/m²)", value=1.47e19, format="%.3e")
    stoich = st.number_input("Cu_s per N2O", value=2.0, step=1.0)
    st.header("Default sample")
    d_mass = st.number_input("Sample mass, oxidized (g)", value=0.3086, format="%.4f", step=0.001)
    d_load = st.number_input("Cu loading (%)", value=80.0, min_value=0.1, max_value=99.9, step=1.0) / 100
    d_basis = st.selectbox("Loading basis", na.LOADING_BASES)
    st.header("Signal treatment")
    normalize = st.checkbox("Divide by inert reference (drift correction)", value=True)
    thr = st.slider("Onset threshold (% of N2 peak)", 1, 30, 3) / 100
    pl_frac = st.slider("Plateau reached at (% of final N2O)", 85, 100, 97) / 100
    pl_len = st.number_input("Plateau window length (s)", value=240.0, step=10.0)

files = st.file_uploader("Mass-spec files", type=["dat", "txt", "csv", "xlsx"], accept_multiple_files=True)
if not files:
    st.info("Upload a file to start. The PV MassSpec 'Plot Data Export' (.dat) works directly.")
    st.stop()


@st.cache_data(show_spinner=False)
def load(data: bytes, name: str) -> pd.DataFrame:
    return na.load_file(data, name)


def fmt(x, f="{:.3g}"):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f.format(x)


def plot_raw(raw: pd.DataFrame, title="Raw signals – use this to pick the right columns"):
    fig = go.Figure()
    for c in raw.columns:
        if c != "time":
            fig.add_scatter(x=raw["time"], y=raw[c], name=str(c), line=dict(width=1.2))
    fig.update_layout(title=title, height=380, template="plotly_white", yaxis_type="log",
                      xaxis_title="Time (s)", yaxis_title="signal (log scale)", legend=dict(orientation="h"),
                      margin=dict(l=10, r=10, t=40, b=10))
    st.plotly_chart(fig, use_container_width=True)


summary_rows = []
tabs = st.tabs([f.name for f in files] + (["📋 Summary"] if len(files) > 1 else []))

for tab, up in zip(tabs, files):
    with tab:
        key = up.name
        try:
            raw = load(up.getvalue(), up.name)
        except Exception as e:  # noqa: BLE001
            st.error(f"Could not read {up.name}: {e}")
            continue

        # --- columns & sample
        guess = na.guess_columns(raw)
        sig_cols = [c for c in raw.columns if c != "time"]
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            n2 = st.selectbox("N2 (m/z 28)", [None] + sig_cols,
                              index=([None] + sig_cols).index(guess["n2"]), key=f"n2_{key}")
        with c2:
            n2o = st.selectbox("N2O (m/z 44)", sig_cols,
                               index=sig_cols.index(guess["n2o"]) if guess["n2o"] in sig_cols else 0, key=f"n2o_{key}")
        with c3:
            ref = st.selectbox("Inert reference", [None] + sig_cols,
                               index=([None] + sig_cols).index(guess["ref"]), key=f"ref_{key}")
        with c4:
            name = st.text_input("Sample name", value=up.name.rsplit(".", 1)[0], key=f"name_{key}")
        c1, c2, c3 = st.columns(3)
        with c1:
            mass = st.number_input("Mass oxidized (g)", value=d_mass, format="%.4f", key=f"m_{key}")
        with c2:
            load_ = st.number_input("Cu loading (%)", value=d_load * 100, key=f"x_{key}") / 100
        with c3:
            basis = st.selectbox("Loading basis", na.LOADING_BASES,
                                 index=na.LOADING_BASES.index(d_basis), key=f"b_{key}")

        if n2o in (n2, ref):
            st.error("N2O must be a different column from N2 and the reference.")
            plot_raw(raw)
            continue
        cols = {"n2": n2, "n2o": n2o, "ref": ref}
        p = na.prepare(raw, cols, normalize=normalize and ref is not None)
        t = p["time"].to_numpy()
        tmin, tmax = float(t.min()), float(t.max())

        try:
            auto = na.auto_windows(p, na.Windows(onset_threshold=thr, plateau_fraction=pl_frac,
                                                 plateau_length_s=pl_len))
        except Exception as e:  # noqa: BLE001
            st.error(f"Automatic window detection failed: {e}. Check the column choice, or set the windows manually.")
            plot_raw(raw)
            auto = na.Windows((tmin, tmin + 0.15 * (tmax - tmin)), tmin + 0.2 * (tmax - tmin),
                              (tmin + 0.5 * (tmax - tmin), tmax))

        # --- windows (auto + manual override)
        with st.expander("Time windows", expanded=False):
            manual = st.checkbox("Set windows manually", key=f"man_{key}")
            if manual:
                bl = st.slider("Baseline (s)", tmin, tmax, auto.baseline, key=f"bl_{key}")
                start = st.slider("Integration start (s)", tmin, tmax, auto.start, key=f"st_{key}")
                pl = st.slider("N2O plateau (start = integration end) (s)", tmin, tmax, auto.plateau, key=f"pl_{key}")
                win = na.Windows(tuple(bl), start, tuple(pl))
            else:
                win = auto
                st.write(f"Baseline **{win.baseline[0]:.0f}–{win.baseline[1]:.0f} s** · "
                         f"integration **{win.start:.1f}–{win.plateau[0]:.1f} s** · "
                         f"plateau **{win.plateau[0]:.0f}–{win.plateau[1]:.0f} s**")

        sample = na.Sample(name=name, mass_ox_g=mass, cu_loading=load_, loading_basis=basis, y_n2o=y_n2o,
                           flow_ml_min=flow, molar_volume_ml=vm, site_density_m2=site, stoich=stoich)
        try:
            res, calc = na.analyse(p, sample, win)
        except Exception as e:  # noqa: BLE001
            st.error(str(e))
            plot_raw(raw)
            continue

        # --- results
        m = st.columns(5)
        m[0].metric("S Cu (m²/g cat, oxidized)", fmt(res.s_per_g_ox, "{:.2f}"))
        m[1].metric("S Cu (m²/g cat, reduced)", fmt(res.s_per_g_red, "{:.2f}"))
        m[2].metric("S Cu (m²/g Cu)", fmt(res.s_per_g_cu, "{:.2f}"))
        m[3].metric("Dispersion", fmt(100 * res.dispersion, "{:.2f} %"))
        m[4].metric("Particle size", fmt(res.d_particle_nm, "{:.1f} nm"))
        agree = None if res.t_equiv_n2_s is None else res.t_equiv_n2_s / res.t_equiv_n2o_s - 1
        st.caption(
            f"N2O consumed {res.n_n2o_mol:.3e} mol · equivalent time {res.t_equiv_n2o_s:.2f} s · "
            f"N2 cross-check {fmt(res.t_equiv_n2_s, '{:.2f}')} s "
            f"({fmt(None if agree is None else 100 * agree, '{:+.1f} %')}) · "
            f"plateau noise {100 * res.plateau_noise:.1f} %"
        )
        for w_ in res.warnings:
            st.warning(w_)
        with st.expander("Show raw signals"):
            plot_raw(raw, "Raw signals (log scale)")

        # --- plots
        def shade(fig):
            fig.add_vrect(x0=win.baseline[0], x1=win.baseline[1], fillcolor="grey", opacity=0.10, line_width=0,
                          annotation_text="baseline", annotation_position="top left")
            fig.add_vrect(x0=win.start, x1=win.plateau[0], fillcolor="green", opacity=0.12, line_width=0,
                          annotation_text="integration", annotation_position="top left")
            fig.add_vrect(x0=win.plateau[0], x1=win.plateau[1], fillcolor="royalblue", opacity=0.10, line_width=0,
                          annotation_text="plateau", annotation_position="top left")
            fig.update_layout(height=380, margin=dict(l=10, r=10, t=40, b=10), legend=dict(orientation="h"),
                              xaxis_title="Time (s)", template="plotly_white")
            return fig

        g1, g2 = st.columns(2)
        with g1:
            fig = go.Figure()
            if n2:
                fig.add_scatter(x=t, y=p["n2"], name="N2 (m/z 28)", line=dict(color="#2E75B6", width=1.3))
            fig.add_scatter(x=t, y=p["n2o"], name="N2O (m/z 44)", line=dict(color="#C00000", width=1.3))
            fig.update_layout(title="Signals" + (" ÷ reference" if normalize and ref else ""),
                              yaxis_title="signal", yaxis_tickformat=".1e")
            st.plotly_chart(shade(fig), use_container_width=True)
        with g2:
            fig = go.Figure()
            fig.add_scatter(x=t, y=calc["f_n2o"], name="N2O out / in", line=dict(color="#C00000", width=1.3))
            if "f_n2" in calc:
                fig.add_scatter(x=t, y=calc["f_n2"], name="N2 formed", line=dict(color="#2E75B6", width=1.3))
                fig.add_scatter(x=t, y=calc["mass_balance"], name="sum (≈ 1)", line=dict(color="#7F7F7F", width=1))
            fig.update_layout(title="Fraction of N2O feed", yaxis_title="fraction", yaxis_range=[-0.2, 1.4])
            st.plotly_chart(shade(fig), use_container_width=True)

        # --- downloads
        report = na.excel_report(raw, cols, calc, sample, res)
        st.download_button("⬇️ Excel report", report, file_name=f"{name}_N2O_report.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                           key=f"dl_{key}")
        summary_rows.append({"file": up.name, "sample": name, **res.summary()})

if len(files) > 1:
    with tabs[-1]:
        df = pd.DataFrame(summary_rows)
        show = df[["sample", "s_per_g_ox", "s_per_g_red", "s_per_g_cu", "dispersion", "d_particle_nm",
                   "n_n2o_mol", "t_equiv_n2o_s", "t_equiv_n2_s", "integration_start_s", "integration_end_s"]]
        st.dataframe(show, use_container_width=True)
        buf = io.StringIO()
        df.to_csv(buf, index=False)
        st.download_button("⬇️ Summary CSV", buf.getvalue(), file_name="N2O_summary.csv", mime="text/csv")

with st.expander("Method"):
    st.markdown(
        r"""
**Reaction:** 2 Cu$_s$ + N$_2$O → Cu$_2$O + N$_2$

1. Signals are divided by the inert reference trace to remove MS drift (optional).
2. N2O out/in $f(t) = (S_{44} - S_{44}^{base}) / (S_{44}^{plateau} - S_{44}^{base})$.
3. Equivalent time of full uptake $t_{eq} = \int_{t_{onset}}^{t_{plateau}} (1 - f)\,dt$.
4. $n_{N_2O} = t_{eq} \cdot y_{N_2O} \cdot F / 60 / V_m$, $\;n_{Cu_s} = 2\,n_{N_2O}$, $\;S = n_{Cu_s} N_A / 1.47\times10^{19}$.
5. Cross-check: N2 produced (m/z 28 minus the N2O fragment, self-calibrated by mass balance).

Onset = first point where N2 exceeds the threshold; plateau = first point where the smoothed N2O reaches the
chosen fraction of its final level. No blank run is subtracted.
"""
    )

