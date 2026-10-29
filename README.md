# N2O Cu surface area – Streamlit app

Calculates the **metallic Cu surface area, dispersion and particle size** from N2O titration
(reactive frontal chromatography) after TPR-H2, using mass-spec data.

**Reaction:** 2 Cu_s + N2O → Cu2O + N2

## What it does
- Reads **PV MassSpec "Plot Data Export" (.dat)** files directly (also CSV/XLSX: first column time in s, then one column per signal).
- Guesses which curve is N2 (m/z 28), N2O (m/z 44) and the inert reference. You can change this.
- Optional **drift correction**: all signals are divided by the inert reference trace.
- **Automatic windows**: baseline, reaction onset (N2 rise) and N2O plateau. A manual override is available.
- Integrates the missing N2O → mol N2O → surface Cu → m²/g (oxidized and reduced catalyst, per g Cu), dispersion, particle size.
- **Cross-check** from the N2 produced (m/z 28 corrected for the N2O fragment). A warning appears if the two disagree by more than 15 %.
- Several files at once → summary table + CSV. An Excel report can be downloaded for each file.

## Run on your computer
```bash
pip install -r requirements.txt
streamlit run app.py
```
Your browser opens at http://localhost:8501.

Self-test: `python tests_check.py` → about 18.3 m²/g for the sample file.

## Put it online (free) with GitHub + Streamlit Community Cloud
1. Create a new repository on GitHub (e.g. `n2o-cu-surface`) and upload all the files in this folder.
2. Go to https://share.streamlit.io, sign in with GitHub and click **New app**.
3. Choose your repository, branch `main`, main file `app.py`, then click **Deploy**.
4. You get a link like `https://n2o-cu-surface.streamlit.app` that you can share with your group.

## Files
| File | Purpose |
|---|---|
| `app.py` | The web interface |
| `n2o_analysis.py` | All the calculations (can also be used from Python/Jupyter) |
| `requirements.txt` | Python packages |
| `sample_data/` | Example run (80Cu20ZrO2, 21.10.26) |
| `tests_check.py` | Quick check that the calculation gives the expected value |

## Assumptions
- Only the Cu surface is oxidized to Cu2O (no bulk oxidation) at the N2O temperature.
- No blank run is subtracted, so the slight smearing of the N2O step by the setup itself is not removed. Starting the integral at the N2 onset keeps this error small.
- Cu surface density 1.47×10¹⁹ atoms/m² (can be changed in the sidebar).

