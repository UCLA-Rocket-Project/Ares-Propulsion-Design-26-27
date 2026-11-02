"""Press bay pressure ladders: COPV -> Dome -> LOX tank, and COPV -> Dome -> Fuel tank.

Each fitting and each pipe run (friction) is its own step on the x-axis.
Edit the INPUTS section and re-run. Graphs open in windows and are saved as PNGs.
"""
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams["font.family"] = "DejaVu Sans"

# ---------------- INPUTS ----------------
# Helium state per segment
RHO_UP = 15.68591672      # kg/m^3, COPV -> dome
V_UP = 60.22176718        # m/s,   COPV -> dome
RHO_DN = 10.47            # kg/m^3, dome -> tanks
V_DN = 80.71609355        # m/s,   dome -> tanks

ID_IN = 0.277             # -6 tube ID, in
F_DARCY = 0.00743         # Darcy friction factor (same as calc sheet)
PA_PER_PSI = 6894.757

P_COPV = 668           # psi
P_REG = 535             # dome regulator set pressure, psi

VDOT_TOTAL = 0.003078757333   # m^3/s
VDOT_LOX = 0.001504544747
VDOT_FUEL = 0.001574212587

# Pipe runs: length (in), fluid state, and flow (vdot used only for check valves).
# Major dP is computed with Darcy-Weisbach: f * (L/D) * ½ρv²
UP = dict(L=0.381 / 0.0254, rho=RHO_UP, v=V_UP, vdot=VDOT_TOTAL)
LOX = dict(L=0.147 / 0.0254, rho=RHO_DN, v=V_DN, vdot=VDOT_LOX)
FUEL = dict(L=1.32 / 0.0254, rho=RHO_DN, v=V_DN, vdot=VDOT_FUEL)

K = {"Tee": 2.0, "Long 90°": 0.7, "90° elbow": 1.5, "Reducer": 0.37,
     "Manifold (tee)": 2.0}

# Check valves: dP = SG * (Q_gpm / Cv)^2
SG_HE = 0.00672
GPM_PER_M3S = 15850.32
CV_LOX, CV_FUEL = 3.20, 1.47

UP_FITTINGS = ["Manifold (tee)", "90° elbow", "Long 90°", "Long 90°", "90° elbow"]
LOX_FITTINGS = ["Reducer", "Tee", "Check valve", "Manifold (tee)", "Long 90°", "Long 90°"]
FUEL_FITTINGS = ["Reducer", "Tee", "Long 90°", "Long 90°", "Check valve",
                 "Manifold (tee)", "Long 90°", "Long 90°", "90° elbow"]
# ----------------------------------------

AREA = math.pi / 4 * (ID_IN * 0.0254) ** 2
BLUE, DARK = "#1f77b4", "0.15"


def q_psi(rho, v):
    """Dynamic pressure ½ρv² in psi."""
    return 0.5 * rho * v ** 2 / PA_PER_PSI


def major_psi(seg):
    """Darcy-Weisbach pipe friction loss in psi."""
    return F_DARCY * (seg["L"] / ID_IN) * q_psi(seg["rho"], seg["v"])


def cv_dp(cv, vdot):
    return SG_HE * (vdot * GPM_PER_M3S / cv) ** 2


def segment_steps(seg, fittings, cv=None):
    q = q_psi(seg["rho"], seg["v"])
    steps = []
    for name in fittings:
        dp = cv_dp(cv, seg["vdot"]) if name == "Check valve" else K[name] * q
        steps.append((name, dp))
    steps.append((f"Pipe friction ({seg['L']:.1f} in)", major_psi(seg)))
    return steps, q


def ladder(title, down_seg, down_fit, tank_label, fname, cv):
    up_steps, q_up = segment_steps(UP, UP_FITTINGS)
    dn_steps, q_dn = segment_steps(down_seg, down_fit, cv)

    # Build the full step sequence: (#, label, dp, panel)
    seq = []
    p = P_COPV
    for name, dp in up_steps:
        seq.append([name, dp, "up"])
        p -= dp
    p_dome_in = p
    seq.append(["Dome regulator", p_dome_in - P_REG, "reg"])
    p = P_REG
    for name, dp in dn_steps:
        seq.append([name, dp, "dn"])
        p -= dp
    p_tank = p
    line_loss = (P_COPV - p_dome_in) + (P_REG - p_tank)

    # Pressure after each step; x = 0 is the COPV, x = i is step i
    pressures = [P_COPV]
    for name, dp, _ in seq:
        pressures.append(P_REG if name == "Dome regulator" else pressures[-1] - dp)
    n = len(seq)
    i_reg = next(i for i, s in enumerate(seq, 1) if s[2] == "reg")

    fig = plt.figure(figsize=(17, 8.5))
    gs = GridSpec(2, 2, width_ratios=[2.1, 1], height_ratios=[1, 1.15],
                  hspace=0.07, wspace=0.12, figure=fig,
                  left=0.06, right=0.98, top=0.93, bottom=0.2)
    ax_t = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_t)
    ax_tab = fig.add_subplot(gs[:, 1])
    ax_tab.axis("off")

    def staircase(i0, i1):
        xs, ys = [i0 - 0.5], [pressures[i0 - 1]]
        for i in range(i0, i1 + 1):
            xs += [i, i, i + 0.5 if i == i1 else i]
            ys += [pressures[i - 1], pressures[i], pressures[i]]
        return xs, ys

    # Top panel: COPV -> dome inlet
    xs, ys = staircase(1, i_reg - 1)
    ax_t.plot([0, 0.5], [P_COPV, P_COPV], color=BLUE, lw=2)
    ax_t.plot(xs, ys, color=BLUE, lw=2)
    ax_t.plot([0], [P_COPV], "o", color=BLUE)
    ax_t.plot([i_reg - 0.5, i_reg], [p_dome_in] * 2, color=BLUE, lw=2)
    ax_t.plot([i_reg, i_reg], [p_dome_in, p_dome_in - 100], color=DARK, lw=2, ls="--")
    ax_t.set_ylim(p_dome_in - 4, P_COPV + 5)

    # Bottom panel: dome outlet -> tank
    ax_b.plot([i_reg, i_reg], [P_REG + 100, P_REG], color=DARK, lw=2, ls="--")
    ax_b.plot([i_reg, i_reg + 0.5], [P_REG] * 2, color=BLUE, lw=2)
    xs, ys = staircase(i_reg + 1, n)
    ax_b.plot(xs, ys, color=BLUE, lw=2)
    ax_b.plot([n + 0.5, n + 1], [p_tank] * 2, color=BLUE, lw=2)
    ax_b.plot([n + 1], [p_tank], "o", color=BLUE)
    ax_b.set_ylim(p_tank - 2.5, P_REG + 2.5)

    # Broken y-axis
    ax_t.spines["bottom"].set_visible(False)
    ax_b.spines["top"].set_visible(False)
    ax_t.tick_params(labelbottom=False, bottom=False)
    d = 0.008
    kw = dict(transform=ax_t.transAxes, color="k", clip_on=False, lw=1)
    ax_t.plot((-d, +d), (-d * 2, +d * 2), **kw)
    ax_t.plot((1 - d, 1 + d), (-d * 2, +d * 2), **kw)
    kw.update(transform=ax_b.transAxes)
    ax_b.plot((-d, +d), (1 - d * 2, 1 + d * 2), **kw)
    ax_b.plot((1 - d, 1 + d), (1 - d * 2, 1 + d * 2), **kw)

    for a in (ax_t, ax_b):
        a.grid(axis="y", alpha=0.35)
        for i in range(1, n + 1):
            a.axvline(i, color="0.6", lw=0.6, ls=":", zorder=0)
        a.axvspan(i_reg - 0.5, i_reg + 0.5, color="0.9", zorder=0)

    # X-axis: one tick per step
    labels = ["COPV"] + [f"{i}. {s[0]}" for i, s in enumerate(seq, 1)] + [tank_label]
    ax_b.set_xticks(range(0, n + 2))
    ax_b.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
    ax_b.set_xlim(-0.5, n + 1.5)

    # Endpoint labels
    lab = dict(fontweight="bold", fontsize=9.5)
    ax_t.annotate(f"COPV\n{P_COPV:.1f} psi", (0, P_COPV), xytext=(6, 6),
                  textcoords="offset points", va="bottom", **lab)
    ax_t.annotate(f"Dome inlet\n{p_dome_in:.1f} psi", (i_reg, p_dome_in), xytext=(8, 4),
                  textcoords="offset points", ha="left", **lab)
    ax_b.annotate(f"Dome outlet\n{P_REG:.1f} psi", (i_reg, P_REG), xytext=(8, 6),
                  textcoords="offset points", ha="left", va="bottom", **lab)
    ax_b.annotate(f"{tank_label}\n{p_tank:.1f} psi", (n + 1, p_tank), xytext=(-6, -8),
                  textcoords="offset points", ha="right", va="top", **lab)

    ax_b.set_xlabel("Component", labelpad=4)
    fig.text(0.015, 0.56, "Pressure (psi)", rotation=90, va="center", ha="center")
    ax_t.set_title(f"{title} — line losses {line_loss:.2f} psi "
                   f"(COPV {P_COPV:.0f} → tank {p_tank:.1f} psi)", fontsize=12)

    # Table
    rows = [[str(i), s[0], f"{s[1]:.2f}", f"{pressures[i]:.1f}"] for i, s in enumerate(seq, 1)]
    rows.append(["", "Total line losses", f"{line_loss:.2f}", ""])
    tab = ax_tab.table(cellText=rows, colLabels=["#", "Component", "ΔP", "P after"],
                       loc="upper center", cellLoc="left", edges="horizontal",
                       colWidths=[0.07, 0.53, 0.19, 0.21])
    tab.auto_set_font_size(False)
    tab.set_fontsize(9.5)
    tab.scale(1, 1.45)
    for (r, col), cell in tab.get_celld().items():
        cell.PAD = 0.04
        if col >= 2:
            cell.set_text_props(ha="right")
        if r == 0 or r == len(rows):
            cell.set_text_props(fontweight="bold")
        if 0 < r <= n and seq[r - 1][2] == "reg":
            cell.set_facecolor("0.93")
    note = (f"Fittings: ΔP = K·½ρv²  ({q_up:.3f} psi/K upstream, {q_dn:.3f} psi/K downstream)\n"
            f"Pipe: ΔP = f·(L/D)·½ρv²,  f = {F_DARCY}\n"
            f"Up: ρ = {RHO_UP:.2f} kg/m³, v = {V_UP:.2f} m/s;  Down: ρ = {RHO_DN:.2f} kg/m³, v = {V_DN:.2f} m/s\n"
            f"Check valve: ΔP = SG·(Q/Cv)²,  Cv = {cv}\n"
            "Line losses exclude the regulator drop.")
    ax_tab.text(0.0, 0.0, note, fontsize=8.5, color="0.3", transform=ax_tab.transAxes)

    fig.savefig(fname, dpi=150, bbox_inches="tight")
    try:
        fig.canvas.manager.set_window_title(title)
    except AttributeError:
        pass
    return p_dome_in, p_tank, line_loss


if __name__ == "__main__":
    for args in (("LOX press ladder", LOX, LOX_FITTINGS, "LOX tank", "lox_press_ladder.png", CV_LOX),
                 ("Fuel press ladder", FUEL, FUEL_FITTINGS, "Fuel tank", "fuel_press_ladder.png", CV_FUEL)):
        print(args[0], [round(v, 3) for v in ladder(*args)])
    plt.show()  # opens each graph in its own window; close them to end the script