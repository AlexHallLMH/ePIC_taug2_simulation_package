#!/usr/bin/env python3
"""
process.py

Analyse UPCGEN ROOT 'particles' tree and produce histograms:
  1) Tau momentum (lab frame)
  2) Tau theta and phi (lab frame)
  3) Tau-daughter theta and phi measured in the tau-pair rest frame,
     where the z-axis is aligned with the tau- momentum for that event.

This script handles either:
 - jagged per-event branches (each tree entry is an event with arrays of particles),
 - or flat per-particle branches (use eventNumber to group particles by event).

Requires: uproot, awkward, numpy, matplotlib
"""

import sys
import math
import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import os

# Directory to save outputs
output_dir = "outputs"

# Create the directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ----------------------------
# I/O and configuration
# ----------------------------
# Optional argument: filename
if len(sys.argv) > 1:
    FILE_NAME = sys.argv[1]
else:
    FILE_NAME = "events.root"

TREE_NAME = "particles"   # from your error message the tree/key is 'particles'
OUT_PREFIX = "upcgen"     # prefix for saved plot files

# ----------------------------
# Helper / physics functions
# ----------------------------

def momentum_magnitude(px, py, pz):
    """Return magnitude of 3-momentum (works on scalars or numpy arrays)."""
    return np.sqrt(px*px + py*py + pz*pz)

def theta_angle(px, py, pz):
    """Polar angle theta in radians (array-safe). Uses safe clip for numerical issues."""
    p = momentum_magnitude(px, py, pz)
    # avoid dividing by zero
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_theta = np.clip(pz / p, -1.0, 1.0)
        theta = np.arccos(cos_theta)
    # Where p == 0 set theta to nan
    if np.isscalar(theta):
        return float(theta)
    return np.where(p > 0, theta, np.nan)

def phi_angle(px, py):
    """Azimuthal angle phi in radians in range (-pi, +pi]."""
    return np.arctan2(py, px)

def lorentz_boost(px, py, pz, E, beta_vec):
    """
    Boost 4-vectors (px,py,pz,E) by velocity -beta_vec (i.e. into frame moving with beta_vec).
    Inputs px,py,pz,E can be scalars or 1D numpy arrays of same shape.
    beta_vec is length-3 numpy array with components (beta_x, beta_y, beta_z).
    Returns px', py', pz', E' (same shape as inputs).
    """
    beta = np.asarray(beta_vec, dtype=float)
    beta2 = float(beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2])
    if beta2 >= 1.0:
        raise ValueError("Invalid boost: |beta| >= 1")
    gamma = 1.0 / math.sqrt(1.0 - beta2) if beta2 > 0 else 1.0
    # bp = beta dot p  (works if px etc are arrays)
    bp = beta[0]*px + beta[1]*py + beta[2]*pz
    gamma2 = (gamma - 1.0) / beta2 if beta2 > 0 else 0.0

    px_p = px + gamma2 * bp * beta[0] + gamma * beta[0] * E
    py_p = py + gamma2 * bp * beta[1] + gamma * beta[1] * E
    pz_p = pz + gamma2 * bp * beta[2] + gamma * beta[2] * E
    E_p  = gamma * (E + bp)
    return px_p, py_p, pz_p, E_p

def rotate_to_z(px, py, pz, target_dir):
    """
    Rotate vectors (px,py,pz) so that target_dir (3-vector) is aligned with +z axis.
    px,py,pz can be scalars or 1D numpy arrays (length Ndau). target_dir is length-3 array.
    Returns rotated px, py, pz (same shape as inputs).
    """
    # Convert target_dir to numpy and normalize
    t = np.asarray(target_dir, dtype=float)
    norm = np.linalg.norm(t)
    if norm == 0:
        raise ValueError("target_dir is zero vector (cannot define rotation)")
    t = t / norm

    # If target_dir already (almost) along +z, no rotation needed
    z = np.array([0.0, 0.0, 1.0], dtype=float)
    cross = np.cross(t, z)
    cross_norm = np.linalg.norm(cross)
    dot = np.dot(t, z)
    # If already aligned or anti-aligned with z, handle specially
    if cross_norm < 1e-12:
        # If dot ~ 1: already aligned; if dot ~ -1: axis = any orthogonal axis, rotate by pi
        if dot > 0:
            return px, py, pz
        else:
            # rotate by pi around x-axis (for example): (x,y,z) -> (-x,-y,z)
            return -px, -py, pz

    # Normalize rotation axis
    kx, ky, kz = cross / cross_norm
    angle = math.acos(np.clip(dot, -1.0, 1.0))

    c = math.cos(angle)
    s = math.sin(angle)
    # Rodrigues rotation matrix components
    R = np.array([
        [c + kx*kx*(1-c),     kx*ky*(1-c) - kz*s, kx*kz*(1-c) + ky*s],
        [ky*kx*(1-c) + kz*s,  c + ky*ky*(1-c),    ky*kz*(1-c) - kx*s],
        [kz*kx*(1-c) - ky*s,  kz*ky*(1-c) + kx*s, c + kz*kz*(1-c)]
    ], dtype=float)

    # Apply rotation: R @ [px;py;pz] for each vector (broadcasted)
    pstack = np.vstack((np.asarray(px, dtype=float), np.asarray(py, dtype=float), np.asarray(pz, dtype=float)))
    p_rot = R.dot(pstack)
    return p_rot[0], p_rot[1], p_rot[2]

# ----------------------------
# Open ROOT file and read branches
# ----------------------------
print(f"Opening file: {FILE_NAME}  tree: {TREE_NAME}")
with uproot.open(FILE_NAME) as f:
    if TREE_NAME not in f:
        raise KeyError(f"Tree '{TREE_NAME}' not found in file. Available keys: {f.keys()}")
    tree = f[TREE_NAME]
    # read the branches we need; if some are missing, this will raise
    branches_to_read = ["px", "py", "pz", "e", "pdgCode", "particleID", "motherID", "eventNumber", "statusID"]
    arrays = tree.arrays(branches_to_read, library="ak")

# ----------------------------
# Detect whether branches are jagged per-entry (per-event) or flat per-particle
# ----------------------------
px_branch = arrays["px"]
# Check jagged vs flat by trying to take len of first element.
# If px_branch[0] has a length, it's jagged (each entry is a list of particles).
# If px_branch[0] is a number, it's flat (each entry is a particle).
try:
    first = px_branch[0]
    _ = len(first)
    jagged = True
except Exception:
    jagged = False

print(f"Detected tree layout: {'jagged per-event' if jagged else 'flat per-particle'}")

# ----------------------------
# Prepare accumulators for histograms
# ----------------------------
tau_momenta = []       # combine tau- and tau+ lab momenta
tau_theta_lab = []     # lab theta for taus
tau_phi_lab = []       # lab phi for taus

daughter_theta_rf = []  # daughter theta in pair rest frame (z along tau-)
daughter_phi_rf = []    # daughter phi

n_events_processed = 0
n_events_skipped_no_taus = 0
n_events_skipped_boost_error = 0

# ----------------------------
# Event loop: jagged or flat
# ----------------------------
if jagged:
    n_events = len(arrays["px"])
    print(f"Number of events (entries) in tree: {n_events}")

    for i in range(n_events):
        # Convert event-level awkward arrays -> numpy arrays (particle lists for this event)
        event_px = np.array(ak.to_list(arrays["px"][i]), dtype=float)
        event_py = np.array(ak.to_list(arrays["py"][i]), dtype=float)
        event_pz = np.array(ak.to_list(arrays["pz"][i]), dtype=float)
        event_e  = np.array(ak.to_list(arrays["e"][i]),  dtype=float)
        event_pdg = np.array(ak.to_list(arrays["pdgCode"][i]), dtype=int)
        event_pid = np.array(ak.to_list(arrays["particleID"][i]), dtype=int)
        event_mid = np.array(ak.to_list(arrays["motherID"][i]), dtype=int)
        # statusID/eventNumber if needed:
        # event_status = np.array(ak.to_list(arrays["statusID"][i]), dtype=int)
        # event_num = arrays["eventNumber"][i]  # may be scalar per event

        # Find taus by PDG code (tau- = +15, tau+ = -15)
        tau_minus_idx = np.where(event_pdg == 15)[0]
        tau_plus_idx  = np.where(event_pdg == -15)[0]

        if len(tau_minus_idx) == 0 or len(tau_plus_idx) == 0:
            n_events_skipped_no_taus += 1
            continue  # skip event; no tau pair found

        # If multiple taus, take the first occurrence (warn if needed)
        if len(tau_minus_idx) > 1 or len(tau_plus_idx) > 1:
            # This is unusual; we pick the first tau- and first tau+
            pass

        im = int(tau_minus_idx[0])
        ip = int(tau_plus_idx[0])

        # Extract tau 4-vectors (lab frame)
        tau_mpx = event_px[im]; tau_mpy = event_py[im]; tau_mpz = event_pz[im]; tau_mE = event_e[im]
        tau_ppx = event_px[ip]; tau_ppy = event_py[ip]; tau_ppz = event_pz[ip]; tau_pE = event_e[ip]

        # Lab-frame tau kinematics
        pm_m = momentum_magnitude(tau_mpx, tau_mpy, tau_mpz)
        pm_p = momentum_magnitude(tau_ppx, tau_ppy, tau_ppz)
        tau_momenta.append(pm_m)
        tau_momenta.append(pm_p)

        tau_theta_lab.append(theta_angle(tau_mpx, tau_mpy, tau_mpz))
        tau_theta_lab.append(theta_angle(tau_ppx, tau_ppy, tau_ppz))

        tau_phi_lab.append(phi_angle(tau_mpx, tau_mpy))
        tau_phi_lab.append(phi_angle(tau_ppx, tau_ppy))

        # Find daughters (both tau- and tau+ daughters)
        pid_m = int(event_pid[im])
        pid_p = int(event_pid[ip])
        # mask of daughters whose motherID equals tau particleID
        dau_mask = (event_mid == pid_m) | (event_mid == pid_p)
        if not np.any(dau_mask):
            # no daughters stored for this event; skip daughter-angle fill
            n_events_processed += 1
            continue

        dau_px = event_px[dau_mask]
        dau_py = event_py[dau_mask]
        dau_pz = event_pz[dau_mask]
        dau_e  = event_e[dau_mask]

        # Build pair 4-vector and compute boost to pair rest frame
        pair_px = tau_mpx + tau_ppx
        pair_py = tau_mpy + tau_ppy
        pair_pz = tau_mpz + tau_ppz
        pair_E  = tau_mE  + tau_pE

        if pair_E == 0:
            n_events_skipped_boost_error += 1
            continue

        beta_vec = np.array([pair_px / pair_E, pair_py / pair_E, pair_pz / pair_E], dtype=float)
        beta2 = np.dot(beta_vec, beta_vec)
        if beta2 >= 1.0:
            # unphysical boost; skip
            n_events_skipped_boost_error += 1
            continue

        # Boost tau- into pair rest frame to get its direction there (for rotation axis)
        tm_px_rf, tm_py_rf, tm_pz_rf, tm_E_rf = lorentz_boost(tau_mpx, tau_mpy, tau_mpz, tau_mE, -beta_vec)
        tau_minus_dir_rf = np.array([tm_px_rf, tm_py_rf, tm_pz_rf], dtype=float)

        # Boost daughters into pair rest frame (vectorised for multiple daughters)
        dau_px_rf, dau_py_rf, dau_pz_rf, dau_e_rf = lorentz_boost(dau_px, dau_py, dau_pz, dau_e, -beta_vec)

        # Rotate daughters so tau- direction is along +z
        px_rot, py_rot, pz_rot = rotate_to_z(dau_px_rf, dau_py_rf, dau_pz_rf, tau_minus_dir_rf)

        # Store daughter angles
        theta_vals = theta_angle(px_rot, py_rot, pz_rot)
        phi_vals   = phi_angle(px_rot, py_rot)

        # theta_vals/phi_vals may be arrays; extend lists
        if np.isscalar(theta_vals):
            daughter_theta_rf.append(float(theta_vals))
            daughter_phi_rf.append(float(phi_vals))
        else:
            daughter_theta_rf.extend(theta_vals.tolist())
            daughter_phi_rf.extend(phi_vals.tolist())

        n_events_processed += 1

else:
    # flat per-particle layout: group by eventNumber
    # Convert branches to numpy arrays
    px_flat = np.asarray(arrays["px"].to_numpy(), dtype=float)
    py_flat = np.asarray(arrays["py"].to_numpy(), dtype=float)
    pz_flat = np.asarray(arrays["pz"].to_numpy(), dtype=float)
    e_flat  = np.asarray(arrays["e"].to_numpy(), dtype=float)
    pdg_flat = np.asarray(arrays["pdgCode"].to_numpy(), dtype=int)
    pid_flat = np.asarray(arrays["particleID"].to_numpy(), dtype=int)
    mid_flat = np.asarray(arrays["motherID"].to_numpy(), dtype=int)
    evt_flat = np.asarray(arrays["eventNumber"].to_numpy())

    unique_events = np.unique(evt_flat)
    print(f"Number of unique eventNumber values: {len(unique_events)}")

    for evt in unique_events:
        mask = (evt_flat == evt)
        event_px = px_flat[mask]
        event_py = py_flat[mask]
        event_pz = pz_flat[mask]
        event_e  = e_flat[mask]
        event_pdg = pdg_flat[mask]
        event_pid = pid_flat[mask]
        event_mid = mid_flat[mask]

        # Find taus
        tau_minus_idx = np.where(event_pdg == 15)[0]
        tau_plus_idx  = np.where(event_pdg == -15)[0]

        if len(tau_minus_idx) == 0 or len(tau_plus_idx) == 0:
            n_events_skipped_no_taus += 1
            continue

        im = int(tau_minus_idx[0])
        ip = int(tau_plus_idx[0])

        # Tau 4-vectors
        tau_mpx = event_px[im]; tau_mpy = event_py[im]; tau_mpz = event_pz[im]; tau_mE = event_e[im]
        tau_ppx = event_px[ip]; tau_ppy = event_py[ip]; tau_ppz = event_pz[ip]; tau_pE = event_e[ip]

        # Lab-frame tau kinematics
        pm_m = momentum_magnitude(tau_mpx, tau_mpy, tau_mpz)
        pm_p = momentum_magnitude(tau_ppx, tau_ppy, tau_ppz)
        tau_momenta.append(pm_m)
        tau_momenta.append(pm_p)

        tau_theta_lab.append(theta_angle(tau_mpx, tau_mpy, tau_mpz))
        tau_theta_lab.append(theta_angle(tau_ppx, tau_ppy, tau_ppz))

        tau_phi_lab.append(phi_angle(tau_mpx, tau_mpy))
        tau_phi_lab.append(phi_angle(tau_ppx, tau_ppy))

        # Daughters of both taus
        pid_m = int(event_pid[im])
        pid_p = int(event_pid[ip])
        dau_mask = (event_mid == pid_m) | (event_mid == pid_p)
        if not np.any(dau_mask):
            n_events_processed += 1
            continue

        dau_px = event_px[dau_mask]
        dau_py = event_py[dau_mask]
        dau_pz = event_pz[dau_mask]
        dau_e  = event_e[dau_mask]

        # Pair 4-vector
        pair_px = tau_mpx + tau_ppx
        pair_py = tau_mpy + tau_ppy
        pair_pz = tau_mpz + tau_ppz
        pair_E  = tau_mE  + tau_pE

        if pair_E == 0:
            n_events_skipped_boost_error += 1
            continue

        beta_vec = np.array([pair_px / pair_E, pair_py / pair_E, pair_pz / pair_E], dtype=float)
        beta2 = np.dot(beta_vec, beta_vec)
        if beta2 >= 1.0:
            n_events_skipped_boost_error += 1
            continue

        # Boost tau- to get direction in pair rest frame
        tm_px_rf, tm_py_rf, tm_pz_rf, tm_E_rf = lorentz_boost(tau_mpx, tau_mpy, tau_mpz, tau_mE, -beta_vec)
        tau_minus_dir_rf = np.array([tm_px_rf, tm_py_rf, tm_pz_rf], dtype=float)

        # Boost daughters to pair rest frame
        dau_px_rf, dau_py_rf, dau_pz_rf, dau_e_rf = lorentz_boost(dau_px, dau_py, dau_pz, dau_e, -beta_vec)

        # Rotate daughters so tau- is along +z
        px_rot, py_rot, pz_rot = rotate_to_z(dau_px_rf, dau_py_rf, dau_pz_rf, tau_minus_dir_rf)

        theta_vals = theta_angle(px_rot, py_rot, pz_rot)
        phi_vals   = phi_angle(px_rot, py_rot)

        if np.isscalar(theta_vals):
            daughter_theta_rf.append(float(theta_vals))
            daughter_phi_rf.append(float(phi_vals))
        else:
            daughter_theta_rf.extend(theta_vals.tolist())
            daughter_phi_rf.extend(phi_vals.tolist())

        n_events_processed += 1

# ----------------------------
# End event loop
# ----------------------------
print(f"Events processed (with taus): {n_events_processed}")
print(f"Events skipped (no tau pair): {n_events_skipped_no_taus}")
print(f"Events skipped (boost/geometry error): {n_events_skipped_boost_error}")

# Convert lists to numpy arrays for plotting
tau_momenta = np.array(tau_momenta, dtype=float) if len(tau_momenta) > 0 else np.array([])
tau_theta_lab = np.array(tau_theta_lab, dtype=float) if len(tau_theta_lab) > 0 else np.array([])
tau_phi_lab = np.array(tau_phi_lab, dtype=float) if len(tau_phi_lab) > 0 else np.array([])
daughter_theta_rf = np.array(daughter_theta_rf, dtype=float) if len(daughter_theta_rf) > 0 else np.array([])
daughter_phi_rf = np.array(daughter_phi_rf, dtype=float) if len(daughter_phi_rf) > 0 else np.array([])

# ----------------------------
# Plot and save histograms
# ----------------------------
def save_hist(data, nbins, xlabel, title, filename, range=None):
    plt.figure()
    plt.hist(data, bins=nbins, histtype="step")
    plt.xlabel(xlabel)
    plt.ylabel("Entries")
    plt.title(title)
    if range is not None:
        plt.xlim(range)
    plt.tight_layout()
    
    # Prepend output directory to filename
    full_path = os.path.join(output_dir, filename)
    plt.savefig(full_path)
    plt.close()
    print(f"Saved: {full_path}")

if tau_momenta.size > 0:
    save_hist(tau_momenta, 50, "Tau momentum (GeV/c)", "Tau momentum (lab frame)", f"{OUT_PREFIX}_tau_momentum_lab.png")
else:
    print("No tau momenta to plot.")

if tau_theta_lab.size > 0:
    save_hist(tau_theta_lab, 50, "Theta (rad)", "Tau theta (lab frame)", f"{OUT_PREFIX}_tau_theta_lab.png", range=(0, math.pi))
else:
    print("No tau theta (lab) to plot.")

if tau_phi_lab.size > 0:
    save_hist(tau_phi_lab, 50, "Phi (rad)", "Tau phi (lab frame)", f"{OUT_PREFIX}_tau_phi_lab.png", range=(-math.pi, math.pi))
else:
    print("No tau phi (lab) to plot.")

if daughter_theta_rf.size > 0:
    save_hist(daughter_theta_rf, 50, "Theta (rad)", "Tau-daughter theta in tau-pair rest frame (z = tau-)", f"{OUT_PREFIX}_daughter_theta_rf.png", range=(0, math.pi))
else:
    print("No daughter theta (rest frame) to plot.")

if daughter_phi_rf.size > 0:
    save_hist(daughter_phi_rf, 50, "Phi (rad)", "Tau-daughter phi in tau-pair rest frame (z = tau-)", f"{OUT_PREFIX}_daughter_phi_rf.png", range=(-math.pi, math.pi))
else:
    print("No daughter phi (rest frame) to plot.")

print("Done.")
