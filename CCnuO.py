# %%
import numpy as np
import matplotlib.pyplot as plt 
import uproot
import awkward as ak
import collections

# %%
#load Sim_D2ODetector013.root
path = "/raid1/genli/G4d2o_simulations/batch_20260330_231702/job_0/Sim_D2ODetector1000.root"
file = uproot.open(path)
# Check the keys in the ROOT file
print(file.keys())

Sim_Tree = file["Sim_Tree"]
# Check the branches in the tree
print(Sim_Tree.keys())
event_data = Sim_Tree["eventData"]
print(event_data.keys())
# Show all branches in the eventData
event_data.show()

# %%
# Inspect pmtHits branch, apply quality cut, then create totalPE per event
if "pmtHits" not in event_data.keys():
    print("pmtHits branch not found. Available branches:")
    print(event_data.keys())
else:
    pmt_hits_branch = event_data["pmtHits"]
    pmt_hits = pmt_hits_branch.array()
    
    # Extract pmtNum jagged array (one variable-length list per event)
    pmt_num = pmt_hits["pmtHits.pmtNum"]
    n_events_all = len(pmt_num)
    
    # This file does not expose a direct per-hit PE branch, so each pmtHits entry
    # is treated as one hit-equivalent P.E. for the per-channel quality cut.
    # Quality cut: each PMT index (0-11) must contribute at least two P.E.
    pass_mask = []
    for evt in pmt_num:
        evt_np = np.asarray(ak.to_numpy(evt), dtype=np.int64)
        counts_12 = np.bincount(evt_np, minlength=12)[:12]
        pass_mask.append(np.all(counts_12 >= 2))
    pass_mask = np.asarray(pass_mask, dtype=bool)
    
    # Apply quality cut
    pmt_num_sel = pmt_num[pass_mask]
    pmt_hits_sel = pmt_hits[pass_mask]
    
    # New event-level entry after quality cut: totalPE = number of PMT hits in each selected event
    totalPE = ak.num(pmt_num_sel, axis=1)
    pmt_hits_sel = ak.with_field(pmt_hits_sel, totalPE, "totalPE")
    
    print(f"Events before cut: {n_events_all}")
    print(f"Events after cut : {len(pmt_num_sel)}")
    print(f"Cut efficiency   : {len(pmt_num_sel) / n_events_all:.4f}")
    print("Quality cut      : every PMT channel has >= 2 hit-equivalent P.E.")
    print("pmtNum (selected) awkward type:", ak.type(pmt_num_sel))
    print("totalPE awkward type:", ak.type(totalPE))
    print("First 10 totalPE values after cut:", ak.to_list(totalPE[:10]))
    
    print("\nFirst 5 selected events pmtNum:")
    for i, evt in enumerate(pmt_num_sel[:5]):
        print(f"Event {i}: {ak.to_list(evt)}")
    
    # Flatten selected hits across events to inspect global pmtNum usage
    pmt_num_flat = ak.to_numpy(ak.flatten(pmt_num_sel))
    print("\nTotal selected PMT hits:", len(pmt_num_flat))
    if len(pmt_num_flat) > 0:
        vals, counts = np.unique(pmt_num_flat, return_counts=True)
        print("Unique PMT numbers:", vals)
        print("Counts per PMT (selected events):")
        for v, c in zip(vals, counts):
            print(f"  PMT {v}: {c}")
    else:
        print("No PMT hits found after quality cut.")

# %%
# Thesis-style histogram of totalPE with Gaussian-like fit around the physical peak (fit range = FWHM)
from scipy.optimize import curve_fit
from scipy.stats import norm

def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def gauss_with_offset(x, A, mu, sigma, C):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + C

totalPE_np = ak.to_numpy(totalPE)
totalPE_pos = totalPE_np[totalPE_np > 0]

if len(totalPE_pos) < 20:
    raise RuntimeError("Too few nonzero totalPE events for a stable peak fit.")

# Thesis-like plotting style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 13,
    "axes.labelsize": 17,
    "axes.titlesize": 22,
    "legend.fontsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 1.2
})

# Fixed binning requested: 0 to 1000 in steps of 10 P.E.
bin_width_pe = 5
bin_min_pe = 0
bin_max_pe = 500
bin_edges_fixed = np.arange(bin_min_pe, bin_max_pe + bin_width_pe, bin_width_pe)

# Full histogram (fixed range 0-1000)
hist_full, edges_full = np.histogram(totalPE_np, bins=bin_edges_fixed)

# Fit histogram (nonzero events only, fixed range 0-1000)
hist_fit, edges_fit = np.histogram(totalPE_pos, bins=bin_edges_fixed)
centers_fit = 0.5 * (edges_fit[:-1] + edges_fit[1:])
bin_width = edges_fit[1] - edges_fit[0]

peak_idx = int(np.argmax(hist_fit))
peak_x = float(centers_fit[peak_idx])
peak_y = float(hist_fit[peak_idx])

if peak_y <= 0:
    raise RuntimeError("Unable to identify a nonzero peak for fitting in 0-1000 P.E. range.")

half_max = peak_y / 2.0

# Find FWHM bounds around nonzero peak
left_idx = peak_idx
while left_idx > 0 and hist_fit[left_idx] >= half_max:
    left_idx -= 1
right_idx = peak_idx
while right_idx < len(hist_fit) - 1 and hist_fit[right_idx] >= half_max:
    right_idx += 1

# Convert scan indices to bins inside the half-max region
left_in = left_idx + 1 if hist_fit[left_idx] < half_max else left_idx
right_in = right_idx - 1 if hist_fit[right_idx] < half_max else right_idx

# Safety fallback if FWHM interval collapses
if right_in < left_in:
    left_in = max(peak_idx - 1, 0)
    right_in = min(peak_idx + 1, len(hist_fit) - 1)

# Slightly widen fit window by one bin on each side for stability
left_in = max(left_in, 0)
right_in = min(right_in, len(hist_fit) - 1)

x_left = float(edges_fit[left_in])
x_right = float(edges_fit[right_in + 1])
fwhm = max(x_right - x_left, 1.0)

# Fit directly to histogram points in the FWHM window
fit_bin_mask = (centers_fit >= x_left) & (centers_fit <= x_right) & (hist_fit > 0)
x_fit = centers_fit[fit_bin_mask]
y_fit = hist_fit[fit_bin_mask].astype(float)

if len(x_fit) < 5:
    raise RuntimeError("Too few nonzero histogram bins in FWHM range for stable Gaussian fit.")

fit_method = "curve_fit (gauss+offset)"
try:
    c0 = max(0.0, float(np.percentile(y_fit, 20)))
    a0 = max(float(np.max(y_fit) - c0), 1.0)
    p0 = [a0, peak_x, max(fwhm / 2.355, 1.0), c0]

    y_err = np.sqrt(np.maximum(y_fit, 1.0))

    popt, pcov = curve_fit(
        gauss_with_offset,
        x_fit,
        y_fit,
        p0=p0,
        sigma=y_err,
        absolute_sigma=True,
        bounds=([0.0, x_left, 1e-6, 0.0], [np.inf, x_right, np.inf, np.inf]),
        maxfev=100000
    )

    A_fit, mu_fit, sigma_fit, C_fit = popt
    perr = np.sqrt(np.maximum(np.diag(pcov), 0.0))
    A_err, mu_err, sigma_err, C_err = perr
except Exception:
    fit_method = "norm.fit fallback"
    fit_data = totalPE_pos[(totalPE_pos >= x_left) & (totalPE_pos <= x_right)]
    if len(fit_data) < 10:
        raise RuntimeError("Fallback fit also failed due to too few events in fit range.")
    mu_fit, sigma_fit = norm.fit(fit_data)
    C_fit = 0.0
    A_fit = len(fit_data) * bin_width / (sigma_fit * np.sqrt(2 * np.pi))
    A_err = np.nan
    mu_err = sigma_fit / np.sqrt(len(fit_data))
    sigma_err = sigma_fit / np.sqrt(2 * len(fit_data))
    C_err = np.nan

x_smooth = np.linspace(x_left, x_right, 400)
y_smooth = gauss_with_offset(x_smooth, A_fit, mu_fit, sigma_fit, C_fit)

n_zero = np.sum(totalPE_np == 0)
n_all = len(totalPE_np)

# Muted thesis palette
hist_color = "#6F9C8D"
edge_color = "#1F2D3A"
fit_color = "#D95F5F"

# Wider rectangular figure
fig, ax = plt.subplots(figsize=(12.8, 6.2), dpi=140)

ax.hist(
    totalPE_np,
    bins=bin_edges_fixed,
    color=hist_color,
    edgecolor=edge_color,
    linewidth=0.55,
    alpha=0.82,
    label=f"totalPE (N={n_all:,})"
 )

ax.plot(
    x_smooth, y_smooth,
    color=fit_color,
    linewidth=2.4,
    label=f"Gaussian+offset fit ($\\mu={mu_fit:.2f}\\pm{mu_err:.2f}$, $\\sigma={sigma_fit:.2f}\\pm{sigma_err:.2f}$)"
 )

ax.axvline(x_left, color="#7A7A7A", linestyle="--", linewidth=1.25)
ax.axvline(x_right, color="#7A7A7A", linestyle="--", linewidth=1.25)

ax.set_title(r"CC $\nu_e$-O Electron Spectrum Simulation", pad=12)
ax.set_xlabel("total P.E.")
ax.set_ylabel(f"Counts/{bin_width_pe} P.E.")
ax.set_axisbelow(True)
ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)

# Keep legend compact and show fit parameters directly in the fit label
legend = ax.legend(
    loc="upper right",
    frameon=True,
    framealpha=0.96,
    fancybox=False,
    edgecolor="#303030"
 )
legend.get_frame().set_linewidth(0.9)

ax.set_xlim(bin_min_pe, bin_max_pe)
plt.tight_layout()
plt.show()

# %%
# Thesis-style histogram of eventData volume IDs with per-type ratios in the legend
from matplotlib.patches import Patch

volume_branch_candidates = [
    key for key in event_data.keys()
    if key in {"vol", "vol0", "eventData/vol", "eventData.vol", "eventData/vol0", "eventData.vol0"}
    or key.endswith(".vol") or key.endswith("/vol")
    or key.endswith(".vol0") or key.endswith("/vol0")
    or key.split(".")[-1] in {"vol", "vol0"}
    or key.split("/")[-1] in {"vol", "vol0"}
]

if not volume_branch_candidates:
    raise KeyError(
        f"No volume branch found in event_data. Available branches: {list(event_data.keys())}"
    )

volume_branch_name = volume_branch_candidates[0]
volume_values_raw = event_data[volume_branch_name].array(library="ak")
volume_values = ak.to_numpy(ak.flatten(volume_values_raw, axis=None))

if volume_values.size == 0:
    raise RuntimeError(f"Branch '{volume_branch_name}' is empty.")

volume_values = volume_values[np.isfinite(volume_values)].astype(int)

if volume_values.size == 0:
    raise RuntimeError(f"Branch '{volume_branch_name}' contains no finite entries.")

volume_counts = collections.Counter(volume_values)
sorted_volume_ids = sorted(volume_counts)
total_entries = sum(volume_counts.values())

volume_labels = {
    -1: "Others",
    0: "Vol 0",
    1: "Tailcatcher",
    2: "In Tank",
    3: "Acrylic Tank",
    4: "PMT 0",
}

bar_labels = [volume_labels.get(volume_id, f"Vol {volume_id}") for volume_id in sorted_volume_ids]
bar_counts = np.array([volume_counts[volume_id] for volume_id in sorted_volume_ids])
bar_ratios = bar_counts / total_entries

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 13,
    "axes.labelsize": 17,
    "axes.titlesize": 22,
    "legend.fontsize": 11,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 1.2,
})

color_map = {
    -1: "#8B98A7",
    0: "#9AA5B1",
    1: "#5B8E7D",
    2: "#C97C5D",
    3: "#7A6C8F",
    4: "#4C72B0",
}
bar_colors = [color_map.get(volume_id, "#4C72B0") for volume_id in sorted_volume_ids]

fig, ax = plt.subplots(figsize=(11.5, 6.2), dpi=140)
bars = ax.bar(
    bar_labels,
    bar_counts,
    color=bar_colors,
    edgecolor="#1F2D3A",
    linewidth=0.9,
    alpha=0.9,
)

for bar, ratio in zip(bars, bar_ratios):
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2.0,
        height + max(0.01 * bar_counts.max(), 1),
        f"{100.0 * ratio:.1f}%",
        ha="center",
        va="bottom",
        fontsize=11,
    )

legend_handles = [
    Patch(
        facecolor=color,
        edgecolor="#1F2D3A",
        label=f"{label}: {count:,} ({100.0 * ratio:.2f}%)",
    )
    for label, count, ratio, color in zip(bar_labels, bar_counts, bar_ratios, bar_colors)
]

ax.set_title(r"CC $\nu_e$-O Interaction Volume Distribution", pad=12)
ax.set_xlabel("Volume Category")
ax.set_ylabel("Counts")
ax.set_axisbelow(True)
ax.grid(True, axis="y", linestyle=":", linewidth=0.8, alpha=0.35)
ax.legend(
    handles=legend_handles,
    loc="upper right",
    frameon=True,
    framealpha=0.96,
    fancybox=False,
    edgecolor="#303030",
    title=f"Branch: {volume_branch_name}",
)
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

print(f"Using branch: {volume_branch_name}")
print(f"Total classified entries: {total_entries:,}")
for volume_id, label, count, ratio in zip(sorted_volume_ids, bar_labels, bar_counts, bar_ratios):
    print(f"{label} (vol={volume_id}): {count:,} entries, {100.0 * ratio:.2f}%")

# %% [markdown]
# Integrated Fluence per SNS year

# %%
Baseline_fluence = 0.0872*8.46*1e14 #1GeV, 1.4MW, 5000 hours, electron neutrino fluence at the detector per cm^2 at 20m distance

# empirical formula for nu per POT vs proton energy
def nu_per_pot(E_p): # E_p in GeV
    p0 = -0.68
    p1 = 1.79
    p2 = -1.12
    p3 = 0.28
    return p3*E_p**3 + p2*E_p**2 + p1*E_p + p0

E_p_now = 1.3
Power_now = 1.9
E_p_ref = 1.0
Power_ref = 1.4
distance_now = 19
distance_ref = 20

fluence_now = Baseline_fluence * (nu_per_pot(E_p_now) / nu_per_pot(E_p_ref)) * ((Power_now / E_p_now) / (Power_ref / E_p_ref)) * (distance_ref / distance_now)**2
print(f"Estimated electron neutrino fluence at the detector for {E_p_now} GeV protons and {Power_now} MW power: {fluence_now:.2e} per cm^2")

detector_face_area_cm2 = 35 * 77.68 * 2.54**2 # cm^2, tail-catcher area facing the target
average_cross_section = 8.257e-42 # cm^2, average CC cross section for 1 GeV electron neutrinos on oxygen
Oxygen_atoms_water = 3.12e28 # total number of oxygen atoms in fiducial and tail catcher volume
Oxygen_acrylic = 8.30e26 # total number of oxygen atoms in acrylic
oxygen_atoms_total = Oxygen_atoms_water + Oxygen_acrylic

# Fluence is already integrated per cm^2 at the detector, so using the total number
# of target atoms means the detector face area must not be multiplied again.
expected_events_SNSyear = fluence_now * average_cross_section * oxygen_atoms_total

print(f"Estimated number of CC electron neutrino interactions in the detector per SNS year: {expected_events_SNSyear:.3e}")

# %%
# Based on the given SNS exposure, sample events from the pre-cut pool,
# then keep only sampled events that pass the per-channel quality cut.
SNS_years = 3
selected_seed = 42
rng = np.random.default_rng(selected_seed)

expected_mean_events = expected_events_SNSyear * SNS_years
Total_events_SNS = expected_mean_events
sampled_event_count = rng.poisson(expected_mean_events)

if "pmt_num" not in globals():
    raise RuntimeError("pmt_num is not available. Run the earlier event-loading and PMT-processing cells first.")

source_event_count = len(pmt_num)
if source_event_count == 0:
    raise RuntimeError("No pre-cut events are available for sampling.")

sampled_indices = rng.choice(source_event_count, size=sampled_event_count, replace=True)
sampled_pmt_num = pmt_num[sampled_indices]

sampled_pass_mask = []
for evt in sampled_pmt_num:
    evt_np = np.asarray(ak.to_numpy(evt), dtype=np.int64)
    counts_12 = np.bincount(evt_np, minlength=12)[:12]
    sampled_pass_mask.append(np.all(counts_12 >= 2))
sampled_pass_mask = np.asarray(sampled_pass_mask, dtype=bool)

sampled_pmt_num_pass = sampled_pmt_num[sampled_pass_mask]
sampled_totalPE = ak.to_numpy(ak.num(sampled_pmt_num_pass, axis=1))
totalPE_source = ak.to_numpy(ak.num(pmt_num, axis=1))

if sampled_totalPE.size == 0:
    raise RuntimeError("No sampled events passed the per-channel quality cut.")

plot_bin_width_pe = 25
plot_bin_min_pe = 0
plot_bin_max_pe = max(500, int(np.ceil(sampled_totalPE.max() / plot_bin_width_pe) * plot_bin_width_pe))
plot_bin_edges = np.arange(plot_bin_min_pe, plot_bin_max_pe + plot_bin_width_pe, plot_bin_width_pe)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 13,
    "axes.labelsize": 17,
    "axes.titlesize": 22,
    "legend.fontsize": 11,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 1.2,
})

fig, ax = plt.subplots(figsize=(12.8, 6.2), dpi=140)
ax.hist(
    sampled_totalPE,
    bins=plot_bin_edges,
    color="#C97C5D",
    edgecolor="#1F2D3A",
    linewidth=0.6,
    alpha=0.88,
    label=f"Sampled events passing cut (N={sampled_totalPE.size:,})",
)

ax.set_title(f"Expected CC $\\nu_e$-O totalPE Spectrum for {SNS_years} SNS years", pad=12)
ax.set_xlabel("total P.E.")
ax.set_ylabel(f"Counts/{plot_bin_width_pe} P.E.")
ax.set_axisbelow(True)
ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
ax.legend(loc="upper right", frameon=True, framealpha=0.96, fancybox=False, edgecolor="#303030")
plt.tight_layout()
plt.show()

print(f"Seed used for reproducible sampling: {selected_seed}")
print(f"Mean expected interactions over {SNS_years} SNS years: {expected_mean_events:.2f}")
print(f"Poisson-sampled events before quality cut: {sampled_event_count}")
print(f"Sampled events passing >= 2 hit-equivalent P.E. per PMT: {sampled_totalPE.size}")
print(f"Pre-cut source event count: {source_event_count}")
print(f"Sampling source file: {path if 'path' in globals() else 'unknown'}")

# %%
# Follow the previous cell with a frequentist pseudo-experiment study.
# For each trial, sample the expected SNS exposure, fit the totalPE spectrum,
# and then plot the distributions of fitted Gaussian means and sigmas.


# -----------------------------
# User-set study parameters
# -----------------------------
times_to_sample = 10000
study_seed = 2026
study_rng = np.random.default_rng(study_seed)

# Histogram binning used inside each pseudo-experiment fit
resample_bin_width_pe = 25
resample_bin_min_pe = 0
resample_bin_max_pe = 500

# Fit-range control for each pseudo-experiment
# Options: "manual" uses fit_range_pe, "fwhm" finds the window from the sampled histogram
fit_range_mode = "manual"
fit_range_pe = (50, 250)

# Minimum statistics required before attempting a fit
minimum_events_for_fit = 20
minimum_bins_in_fit_window = 3

if "curve_fit" not in globals():
    from scipy.optimize import curve_fit
if "norm" not in globals():
    from scipy.stats import norm

if "gauss_with_offset" not in globals():
    def gauss_with_offset(x, A, mu, sigma, C):
        return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + C

if "expected_mean_events" not in globals():
    raise RuntimeError("expected_mean_events is not available. Run the previous SNS exposure cell first.")
if "totalPE" not in globals():
    raise RuntimeError("totalPE is not available. Run the PMT-processing cell first.")
if "pass_mask" not in globals():
    raise RuntimeError("pass_mask is not available. Run the PMT-processing cell first.")

# The quality cut is deterministic for each source event, so repeated
# sample-then-cut pseudo-experiments are exactly equivalent to Poisson-thinning
# the mean rate and sampling from the already selected event spectrum.
selected_totalPE_source = ak.to_numpy(totalPE).astype(float)
selected_totalPE_source = selected_totalPE_source[selected_totalPE_source > 0]
pass_efficiency = float(np.mean(np.asarray(pass_mask, dtype=float)))
expected_mean_events_after_cut = expected_mean_events * pass_efficiency

if selected_totalPE_source.size < minimum_events_for_fit:
    raise RuntimeError("Too few selected events are available for repeated fitting.")

resample_bin_edges = np.arange(
    resample_bin_min_pe,
    resample_bin_max_pe + resample_bin_width_pe,
    resample_bin_width_pe
 )

def fit_sampled_totalpe(totalpe_values, bin_edges, fit_range_mode, fit_range_pe):
    totalpe_values = np.asarray(totalpe_values, dtype=float)
    totalpe_values = totalpe_values[totalpe_values > 0]

    if totalpe_values.size < minimum_events_for_fit:
        raise RuntimeError("Too few nonzero totalPE events for a stable fit.")

    hist_fit, edges_fit = np.histogram(totalpe_values, bins=bin_edges)
    centers_fit = 0.5 * (edges_fit[:-1] + edges_fit[1:])

    peak_idx = int(np.argmax(hist_fit))
    peak_x = float(centers_fit[peak_idx])
    peak_y = float(hist_fit[peak_idx])

    if peak_y <= 0:
        raise RuntimeError("Unable to identify a nonzero peak for fitting.")

    if fit_range_mode == "manual":
        x_left, x_right = map(float, fit_range_pe)
        if x_right <= x_left:
            raise RuntimeError("fit_range_pe must satisfy upper > lower.")
    elif fit_range_mode == "fwhm":
        half_max = peak_y / 2.0
        left_idx = peak_idx
        while left_idx > 0 and hist_fit[left_idx] >= half_max:
            left_idx -= 1

        right_idx = peak_idx
        while right_idx < len(hist_fit) - 1 and hist_fit[right_idx] >= half_max:
            right_idx += 1

        left_in = left_idx + 1 if hist_fit[left_idx] < half_max else left_idx
        right_in = right_idx - 1 if hist_fit[right_idx] < half_max else right_idx

        if right_in < left_in:
            left_in = max(peak_idx - 1, 0)
            right_in = min(peak_idx + 1, len(hist_fit) - 1)

        x_left = float(edges_fit[left_in])
        x_right = float(edges_fit[right_in + 1])
    else:
        raise RuntimeError("fit_range_mode must be either 'manual' or 'fwhm'.")

    fit_bin_mask = (centers_fit >= x_left) & (centers_fit <= x_right) & (hist_fit > 0)
    x_fit = centers_fit[fit_bin_mask]
    y_fit = hist_fit[fit_bin_mask].astype(float)

    if len(x_fit) < minimum_bins_in_fit_window:
        raise RuntimeError("Too few populated bins in the fit window.")

    fit_width_guess = max((x_right - x_left) / 2.355, 1.0)
    fit_method = "curve_fit"
    try:
        c0 = max(0.0, float(np.percentile(y_fit, 20)))
        a0 = max(float(np.max(y_fit) - c0), 1.0)
        p0 = [a0, peak_x, fit_width_guess, c0]
        y_err = np.sqrt(np.maximum(y_fit, 1.0))

        popt, _ = curve_fit(
            gauss_with_offset,
            x_fit,
            y_fit,
            p0=p0,
            sigma=y_err,
            absolute_sigma=True,
            bounds=([0.0, x_left, 1e-6, 0.0], [np.inf, x_right, np.inf, np.inf]),
            maxfev=100000
        )
        _, mu_fit_trial, sigma_fit_trial, _ = popt
    except Exception:
        fit_method = "norm.fit fallback"
        fit_data = totalpe_values[(totalpe_values >= x_left) & (totalpe_values <= x_right)]
        if fit_data.size < 10:
            raise RuntimeError("Fallback fit also failed due to too few events in the fit window.")
        mu_fit_trial, sigma_fit_trial = norm.fit(fit_data)

    return float(mu_fit_trial), float(sigma_fit_trial), fit_method

fit_mu_values = []
fit_sigma_values = []
accepted_counts = []
fit_method_counts = collections.Counter()
failed_trials = 0

for trial_idx in range(times_to_sample):
    sampled_pass_count = study_rng.poisson(expected_mean_events_after_cut)
    if sampled_pass_count < minimum_events_for_fit:
        failed_trials += 1
        continue

    sampled_totalPE_trial = study_rng.choice(
        selected_totalPE_source,
        size=sampled_pass_count,
        replace=True
    )

    try:
        mu_trial, sigma_trial, fit_method_trial = fit_sampled_totalpe(
            sampled_totalPE_trial,
            resample_bin_edges,
            fit_range_mode,
            fit_range_pe
        )
    except RuntimeError:
        failed_trials += 1
        continue

    fit_mu_values.append(mu_trial)
    fit_sigma_values.append(sigma_trial)
    accepted_counts.append(sampled_pass_count)
    fit_method_counts[fit_method_trial] += 1

fit_mu_values = np.asarray(fit_mu_values, dtype=float)
fit_sigma_values = np.asarray(fit_sigma_values, dtype=float)
accepted_counts = np.asarray(accepted_counts, dtype=int)

if fit_mu_values.size == 0:
    raise RuntimeError("No pseudo-experiments produced a stable fit.")

mu_mean = fit_mu_values.mean()
mu_std = fit_mu_values.std(ddof=1)
sigma_mean = fit_sigma_values.mean()
sigma_std = fit_sigma_values.std(ddof=1)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 13,
    "axes.labelsize": 17,
    "axes.titlesize": 20,
    "legend.fontsize": 11,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 1.2,
})

fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6), dpi=140)

# mu_bins = min(40, max(12, int(np.sqrt(fit_mu_values.size))))
# sigma_bins = min(40, max(12, int(np.sqrt(fit_sigma_values.size))))
mu_binwidth = 5
sigma_binwidth = 5
mu_bins = np.arange(
    max(0, int(np.floor((fit_mu_values.min() - mu_binwidth) / mu_binwidth) * mu_binwidth)),
    int(np.ceil((fit_mu_values.max() + mu_binwidth) / mu_binwidth) * mu_binwidth) + mu_binwidth,
    mu_binwidth
)
sigma_bins = np.arange(
    max(0, int(np.floor((fit_sigma_values.min() - sigma_binwidth) / sigma_binwidth) * sigma_binwidth)),
    int(np.ceil((fit_sigma_values.max() + sigma_binwidth) / sigma_binwidth) * sigma_binwidth) + sigma_binwidth,
    sigma_binwidth
)

axes[0].hist(
    fit_mu_values,
    bins=mu_bins,
    color="#5B8E7D",
    edgecolor="#1F2D3A",
    linewidth=0.7,
    alpha=0.88,
 )
axes[0].axvline(
    mu_mean,
    color="#B03A2E",
    linestyle="--",
    linewidth=2.0,
    label=f"Mean = {mu_mean:.2f} P.E."
 )
axes[0].axvspan(
    mu_mean - mu_std,
    mu_mean + mu_std,
    color="#B03A2E",
    alpha=0.25,
    label=f"Std = {mu_std:.2f} P.E."
)
axes[0].set_title(r"Distribution of fitted $\mu$")
axes[0].set_xlabel("Fitted mean total P.E.")
axes[0].set_ylabel(f"Experiments/{mu_binwidth} P.E.")
axes[0].grid(True, linestyle=":", linewidth=0.8, alpha=0.35)
axes[0].legend()

axes[1].hist(
    fit_sigma_values,
    bins=sigma_bins,
    color="#C97C5D",
    edgecolor="#1F2D3A",
    linewidth=0.7,
    alpha=0.88,
 )
axes[1].axvline(
    sigma_mean,
    color="#B03A2E",
    linestyle="--",
    linewidth=2.0,
    label=f"Mean = {sigma_mean:.2f} P.E."
 )
axes[1].axvspan(
    sigma_mean - sigma_std,
    sigma_mean + sigma_std,
    color="#B03A2E",
    alpha=0.25,
    label=f"Std = {sigma_std:.2f} P.E."
)
axes[1].set_title(r"Distribution of fitted $\sigma$")
axes[1].set_xlabel("Fitted width in total P.E.")
axes[1].set_ylabel(f"Experiments/{sigma_binwidth} P.E.")
axes[1].grid(True, linestyle=":", linewidth=0.8, alpha=0.35)
axes[1].legend()

fig.suptitle(
    f"Frequentist pseudo-experiments for {SNS_years} SNS years (N={fit_mu_values.size} successful fits)",
    y=1.03
 )
plt.tight_layout()
plt.show()

print(f"Study seed: {study_seed}")
print(f"Requested pseudo-experiments: {times_to_sample}")
print(f"Successful fits: {fit_mu_values.size}")
print(f"Failed fits: {failed_trials}")
print(f"Pseudo-experiment bin width: {resample_bin_width_pe} P.E.")
print(f"Pseudo-experiment histogram range: [{resample_bin_min_pe}, {resample_bin_max_pe}] P.E.")
if fit_range_mode == "manual":
    print(f"Manual fit range: [{fit_range_pe[0]}, {fit_range_pe[1]}] P.E.")
else:
    print("Fit range mode: FWHM from each pseudo-experiment histogram")
print(f"Quality-cut efficiency from source sample: {pass_efficiency:.4f}")
print(f"Mean expected events before cut: {expected_mean_events:.2f}")
print(f"Mean expected events after cut : {expected_mean_events_after_cut:.2f}")
print(f"Average sampled events after cut: {accepted_counts.mean():.2f}")
print(f"Fitted mu    = {mu_mean:.3f} +/- {mu_std:.3f}")
print(f"Fitted sigma = {sigma_mean:.3f} +/- {sigma_std:.3f}")
print("Fit methods used:")
for method_name, count in fit_method_counts.items():
    print(f"  {method_name}: {count}")


