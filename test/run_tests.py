#!/usr/bin/env python3
"""
SHARK analytical validation suite.

    python3 run_tests.py                 # standard tests (~5-10 min)
    python3 run_tests.py --long          # + planet type-I torque (~15-25 min more)
    python3 run_tests.py --only si_linA fargo_bump
    python3 run_tests.py --threads 8

Each test builds SHARK with its own flags, runs it in test/validation_runs/<test>/,
compares with an analytical solution and prints PASS/FAIL. Exit code = number of failures.

Tests
  disk_ic_drift      t=0 dust radial drift vs terminal-velocity formula (DISK_DUST)
  disk_equilibrium   gas disk without planet stays in equilibrium (DISK_EQ, FARGO)
  dust_drift_nsh     dust-gas radial slip vs Nakagawa et al. (1986) (DISK_DUST, FARGO)
  fargo_bump         dust bump advected by the rotation vs exact solution, FARGO on/off
  si_linA            streaming-instability growth rate vs Youdin & Johansen (2007) linA
  (sod, soundwave, dustybox, dustywave, dustyshock and pscal run along x and along y, and check x/y symmetry)
  sod                Sod shock tube vs exact Riemann solution (SETUP=tests)
  soundwave          isothermal sound wave, order of convergence (32, 64, 128 cells)
  dustybox           uniform drag relaxation, 1 and 3 species (Laibe & Price 2011)
  dustywave          linear gas+dust wave vs exact solution (Laibe & Price 2011)
  dustyshock         stationary dusty shock vs steady ODE solution (Loren-Aguilar & Bate 2014)
  pscal              dust passive scalar: advection (exact, convergence, conservation) and uniform scalar
  kernel_constant    Smoluchowski coagulation, constant kernel vs exact solution
  kernel_additive    Smoluchowski coagulation, additive kernel vs exact solution (Golovin)
  fragmentation      coagulation + fragmentation: steady state and barrier scaling a_frag ~ v_frag^2
  disk_kernel        vertically integrated coagulation kernel in a disk cell (GEOM=2), well mixed and settled
  planet_torque      (--long) type-I torque vs Paardekooper et al. (2010)
  sph_hydrostatic    locally isothermal disk on the spherical grid stays hydrostatic, converges (DISK_SPH)
"""
import argparse, glob, os, re, shutil, subprocess, sys, time
import numpy as np

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # .../Disks
BIN   = os.path.join(ROOT, "bin")
SUITE = os.path.join(ROOT, "test", "test-suite")
WORK  = os.path.join(ROOT, "test", "validation_runs")
THREADS = 4
_built = None

# ----------------------------------------------------------------------------- helpers
def build(flags):
    global _built
    if flags == _built:
        return
    r = subprocess.run(f"make clean > /dev/null 2>&1; make {flags} > make.log 2>&1", shell=True, cwd=BIN)
    if r.returncode != 0 or not os.path.exists(os.path.join(BIN, "shark")):
        raise RuntimeError(f"build failed ({flags}), see {BIN}/make.log")
    _built = flags

def edit_nml(text, **groups):
    """groups: grid={...}, setup={...}, dust={...}. Replaces existing keys, inserts missing ones."""
    names = {"grid": "GRID_PARAMS", "setup": "SETUP_PARAMS", "dust": "DUST_PARAMS"}
    for g, params in groups.items():
        for k, v in params.items():
            pat = re.compile(rf"^\s*{k}\s*=.*$", re.I | re.M)
            if pat.search(text):
                text = pat.sub(f"{k} = {v}", text, count=1)
            else:
                text = re.sub(rf"^(&{names[g]}\s*)$", rf"\1\n{k} = {v}", text, count=1, flags=re.I | re.M)
    return text

def nml_value(text, key, default):
    m = re.search(rf"^\s*{key}\s*=\s*([^\s,/!]+)", text, re.I | re.M)
    return float(m.group(1).lower().replace("d", "e")) if m else default

def run(name, base_nml, flags, threads=None, **groups):
    build(flags)
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    txt = edit_nml(open(base_nml).read(), **groups)
    open(os.path.join(d, "run.nml"), "w").write(txt)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads or THREADS))
    with open(os.path.join(d, "log.txt"), "w") as log:
        r = subprocess.run([os.path.join(BIN, "shark"), "run.nml"], cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"run failed, see {d}/log.txt")
    return d, txt

def outputs(d):
    return sorted(glob.glob(os.path.join(d, "output_*")))

def info(o):
    out = {}
    for l in open(os.path.join(o, "info.dat")):
        k, _, v = l.partition("=")
        try: out[k.strip()] = float(v)
        except ValueError: pass
    return out

def field(o, var):
    i = info(o)
    return np.fromfile(os.path.join(o, var))[: int(i["NX"])*int(i["NY"])].reshape(int(i["NY"]), int(i["NX"]))

def ring_radii(o):
    return field(o, "x")[0]

FIGDIR = os.path.join(WORK, "figures")
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
except ImportError:
    plt = None

def savefig(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, name + ".png"), dpi=110); plt.close(fig)

def polar_map(ax, phi, r, f, title, vmin, vmax, cmap, rmax):
    ax.pcolormesh(phi, r, f.T, vmin=vmin, vmax=vmax, cmap=cmap, shading="nearest")
    ax.set_ylim(0, rmax); ax.set_title(title, fontsize=10); ax.set_xticklabels([]); ax.tick_params(labelsize=7)

DISK = "SETUP=disk_polar GEOM=2 GRIDSPACE=1 NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=1 SOLVERDUST=0"
SQ = np.sqrt(np.pi/8.0)*np.sqrt(2.0*np.pi)       # St = SQ * rho_grain * s / Sigma (Epstein, midplane, isothermal)

# ----------------------------------------------------------------------------- tests
def t_disk_ic_drift():
    base = os.path.join(SUITE, "DISK_DUST", "disk_dust.nml")
    d, txt = run("disk_ic_drift", base, DISK + " NX=64 NY=64 NDUST=1", grid={"tend": "1d-6"})
    o = outputs(d)[0]; i = info(o)
    h, rout, s = nml_value(txt, "hoverr", 0.1), nml_value(txt, "R_out_disk", 40.0), nml_value(txt, "scut", 0.1)
    r, sig, vd = field(o, "x"), field(o, "rho"), field(o, "vd")
    eta = 0.5*h**2*(2.0 + r/rout)
    vk  = np.sqrt(1.0/r)*i["unit_v"]
    St  = SQ*2.7*s/sig
    vexp = 2.0*(np.sqrt(1.0 - 2.0*eta) - 1.0)*vk/(St + 1.0/St)
    err = np.max(np.abs(vd/vexp - 1.0))
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].loglog(r[0], -vd.mean(0), "o", ms=3, label="SHARK"); ax[0].loglog(r[0], -vexp.mean(0), "k-", label="analytic")
        ax[0].set_xlabel("r [au]"); ax[0].set_ylabel("-v_r,dust [cm/s]"); ax[0].legend(); ax[0].set_title("Initial dust radial drift")
        ax[1].semilogy(r[0], np.abs(vd/vexp - 1).max(0) + 1e-17); ax[1].set_xlabel("r [au]"); ax[1].set_ylabel("|v/v_exp - 1|")
        ax[1].set_title("Relative error")
        savefig(fig, "disk_ic_drift")
    return err < 1e-6, f"max |v/v_exp - 1| = {err:.1e} (tol 1e-6)"

def t_disk_equilibrium():
    base = os.path.join(SUITE, "DISK_EQ", "disk.nml")
    med, prof = {}, {}
    for n in (64, 128):
        d, txt = run(f"disk_equilibrium_{n}", base, DISK + f" NX={n} NY={n} NDUST=0",
                     grid={"tend": "200.0", "fargo": ".true."})
        o0, o1 = outputs(d)[0], outputs(d)[-1]
        box = nml_value(txt, "box_l", 400.0)
        r = ring_radii(o1); free = (r > 2.0) & (r < 0.8*box)              # outside the damping zones
        dev = np.abs(field(o1, "rho")/field(o0, "rho") - 1.0)
        med[n] = np.median(dev[:, free]); prof[n] = (r, np.median(dev, axis=0))
    ratio = med[64]/med[128]
    ok = ratio > 3.0 and med[128] < 1e-2
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for n in prof:
            ax[0].loglog(*prof[n], label=f"{n}^2")
        ax[0].axvspan(1, 2, color="0.9"); ax[0].axvspan(0.8*box, box, color="0.9")
        ax[0].set_xlabel("r [au]"); ax[0].set_ylabel("median |dSigma/Sigma|"); ax[0].legend(); ax[0].set_title("Drift from equilibrium after 200 yr")
        ns = np.array(list(med)); ax[1].loglog(ns, [med[n] for n in ns], "o-", label="SHARK")
        ax[1].loglog(ns, med[64]*(ns/64.0)**-2, "k--", label="2nd order"); ax[1].set_xlabel("N"); ax[1].legend()
        ax[1].set_ylabel("median |dSigma/Sigma| (free disk)"); ax[1].set_title("Convergence")
        savefig(fig, "disk_equilibrium")
    return ok, (f"median |dSigma/Sigma| after 200 yr: {med[64]:.1e} (64^2) -> {med[128]:.1e} (128^2), "
                f"ratio {ratio:.1f} (2nd order: > 3, and < 1e-2 at 128^2)")

def t_dust_drift_nsh():
    base = os.path.join(SUITE, "DISK_DUST", "disk_dust.nml")
    d, txt = run("dust_drift_nsh", base, DISK + " NX=128 NY=128 NDUST=1",
                 grid={"tend": "10.0", "fargo": ".true."})
    o = outputs(d)[-1]; i = info(o)
    h, rout, s = nml_value(txt, "hoverr", 0.1), nml_value(txt, "R_out_disk", 40.0), nml_value(txt, "scut", 0.1)
    box = nml_value(txt, "box_l", 40.0)
    r = ring_radii(o)
    sig, sigd = field(o, "rho").mean(0), field(o, "rhod").mean(0)
    slip = (field(o, "vd") - field(o, "v")).mean(0)
    eps, St = sigd/sig, SQ*2.7*s/sig
    etavk = 0.5*h**2*(2.0 + r/rout)*np.sqrt(1.0/r)*i["unit_v"]
    nsh = -2.0*St*(1.0 + eps)/((1.0 + eps)**2 + St**2)*etavk             # dust - gas radial velocity
    sel = (r > 3.0) & (r < 0.8*box)
    rel = np.abs(slip[sel]/nsh[sel] - 1.0)
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].loglog(r, -slip, "o", ms=3, label="SHARK"); ax[0].loglog(r, -nsh, "k-", label="Nakagawa+86")
        ax[0].axvspan(0.8*box, box, color="0.9"); ax[0].set_xlabel("r [au]"); ax[0].set_ylabel("-(v_dust - v_gas) [cm/s]")
        ax[0].legend(); ax[0].set_title("Dust-gas radial slip")
        ax[1].semilogy(r[sel], rel); ax[1].set_xlabel("r [au]"); ax[1].set_ylabel("|slip/NSH - 1|"); ax[1].set_title("Relative error")
        savefig(fig, "dust_drift_nsh")
    return np.median(rel) < 0.1, f"median |slip/NSH - 1| = {np.median(rel):.2e}, max {rel.max():.2e}, 3 < r < {0.8*box:.0f} au (tol 0.1 on median)"

def t_fargo_bump():
    base = os.path.join(SUITE, "DISK_DUST", "disk_dust.nml")
    rb, w, h, rout = 10.0, 2.0, 0.1, 40.0
    setup = {"bump_amp": "1.0", "bump_r": rb, "bump_phi": "3.0", "bump_w": w}
    res, maps, geo = {}, {}, None
    for m in ("off", "on"):
        d, txt = run(f"fargo_bump_{m}", base, DISK + " NX=128 NY=128 NDUST=1",
                     grid={"tend": f"{rb**1.5:.4f}", "fargo": ".true." if m == "on" else ".false."}, setup=setup)
        o0, o1 = outputs(d)[0], outputs(d)[-1]; i = info(o1); NY = int(i["NY"]); t = i["time"]
        r, phi = ring_radii(o1), (np.arange(NY) + 0.5)*2*np.pi/NY
        d0, d1 = field(o0, "rhod"), field(o1, "rhod")
        Om = np.sqrt(1.0/r**3)*np.sqrt(1.0 - h**2*(2.0 + r/rout))
        ana = np.stack([np.interp(phi - Om[k]*t, phi, d0[:, k], period=2*np.pi) for k in range(len(r))], axis=1)
        ann = np.abs(r - rb) < 3*w; bg = np.median(d0[:, ann], axis=0)
        err = np.abs(d1[:, ann] - ana[:, ann]).sum()/(ana[:, ann] - bg).sum()
        ib = np.argmin(np.abs(r - rb))
        dph = ((phi[np.argmax(d1[:, ib])] - phi[np.argmax(ana[:, ib])] + np.pi) % (2*np.pi) - np.pi)*NY/(2*np.pi)
        res[m] = (err, dph)
        maps[m] = d1; maps["ana"] = ana; maps["bg"] = np.median(d0, axis=0)
        geo = (phi, r)
    ok = res["on"][0] < res["off"][0] and res["on"][0] < 0.8 and abs(res["on"][1]) <= 1.0
    if plt:
        phi, r = geo
        fig, ax = plt.subplots(1, 3, subplot_kw=dict(projection="polar"), figsize=(14, 5))
        for a, k, ttl in zip(ax, ("off", "on", "ana"), ("SHARK, FARGO off", "SHARK, FARGO on", "exact advection")):
            polar_map(a, phi, r, maps[k]/maps["bg"] - 1, ttl + (f"\nL1 error {res[k][0]:.3f}" if k in res else ""), 0, 1, "magma", 2*rb)
        fig.suptitle("Dust bump after one orbit at 10 au (rho_d/bg - 1)")
        savefig(fig, "fargo_bump")
    return ok, (f"L1 error: FARGO {res['on'][0]:.3f} vs no-FARGO {res['off'][0]:.3f} (FARGO < no-FARGO and < 0.8); "
                f"phase {res['on'][1]:+.1f} cells (|.| <= 1)")

def t_si_linA():
    base = os.path.join(SUITE, "SI_BA", "si.nml")
    L = 2.0*np.pi/30.0                                                       # one wavelength, K = k eta r = 30
    d, txt = run("si_linA", base, "SETUP=SI GEOM=0 NX=64 NY=64 NDUST=1 NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=2 SOLVERDUST=0",
                 grid={"tend": "10.0", "freq_out": "1000"},
                 setup={"box_l": f"{L:.10f}", "box_l_y": f"{L:.10f}", "dust2gas_species": "3.0",
                        "stokes_species": "0.1", "lin_test": "1", "lin_amp": "1d-6"})
    t, amp = [], []
    for o in outputs(d):
        rd = field(o, "rhod"); F = np.fft.fft2(rd/rd.mean() - 1.0)
        t.append(info(o)["time"]); amp.append(max(abs(F[1, 1]), abs(F[1, -1]), abs(F[-1, 1]), abs(F[-1, -1])))
    t, amp = np.array(t), np.array(amp)
    sel = t >= 2.0
    if sel.sum() < 3:
        return False, f"not enough outputs ({len(t)}) after t = 2"
    s, c = np.polyfit(t[sel], np.log(amp[sel]), 1)
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].semilogy(t, amp, "o", ms=3, label="SHARK (|FFT| of rho_d)")
        ax[0].semilogy(t, np.exp(c + s*t), "b-", lw=1, label=f"fit: s = {s:.4f}")
        ax[0].semilogy(t, np.exp(c + s*t[-1] - 0.419*(t[-1] - t)), "k--", label="linA: s = 0.4190")
        ax[0].set_xlabel("t [1/Omega]"); ax[0].set_ylabel("mode amplitude"); ax[0].legend(); ax[0].set_title("Streaming instability, linA")
        rd = field(outputs(d)[-1], "rhod"); im = ax[1].imshow(rd/rd.mean() - 1, origin="lower", cmap="RdBu_r")
        fig.colorbar(im, ax=ax[1]); ax[1].set_title(f"rho_d/<rho_d> - 1 at t = {t[-1]:.1f}"); ax[1].set_xlabel("x (radial)"); ax[1].set_ylabel("z (vertical)")
        savefig(fig, "si_linA")
    return abs(s/0.4190 - 1.0) < 0.1, f"growth rate {s:.4f} Omega vs 0.4190 (linA), error {100*(s/0.419 - 1):+.1f}% (tol 10%)"

def t_planet_torque():
    base = os.path.join(SUITE, "DISK_EQ", "disk.nml")
    q, rp, h, rout = 1e-4, 10.0, 0.1, 40.0
    d, txt = run("planet_torque", base, DISK + " NX=256 NY=256 NDUST=0",
                 grid={"tend": f"{10*rp**1.5:.4f}", "freq_out": "5000", "fargo": ".true."},
                 setup={"box_l": "30.0", "hoverr": h, "R_out_disk": rout, "Mplanet": q, "Rplanet": rp,
                        "nplanet": "1", "phi_planet": "0.0", "nramp_planet": "2.0", "a_smooth_pl": "0.4",
                        "freq_torque": "20"})
    o0 = outputs(d)[0]; i = info(o0)
    sigp = np.interp(rp, ring_radii(o0), field(o0, "rho").mean(0)/(i["unit_d"]*i["unit_l"]))
    Omp = np.sqrt((1.0 + q)/rp**3); G0 = (q/h)**2*sigp*rp**4*Omp**2
    expected = -0.85 - (1.0 + rp/rout) - 0.9
    tq = np.loadtxt(os.path.join(d, "planet_torques.dat"), ndmin=2)
    orb = tq[:, 0]*Omp/(2*np.pi); g = tq[:, 3]/G0; win = (orb >= 3) & (orb <= 10)
    meas, osc = g[win].mean(), g[win].std()
    ok = abs(meas/expected - 1.0) < 0.15 and osc < 0.1*abs(expected)
    if plt:
        o1 = outputs(d)[-1]; r = ring_radii(o1); NY = int(info(o1)["NY"]); phi = (np.arange(NY) + 0.5)*2*np.pi/NY
        fig = plt.figure(figsize=(12, 5))
        a1 = fig.add_subplot(121, projection="polar")
        polar_map(a1, phi, r, field(o1, "rho")/field(o0, "rho") - 1, "Sigma/Sigma_0 - 1", -0.1, 0.1, "RdBu_r", 2*rp)
        a2 = fig.add_subplot(122); a2.plot(orb, g, label="SHARK"); a2.axhline(expected, ls="--", c="k", label="Paardekooper+10")
        a2.axvspan(3, 10, color="0.93"); a2.set_xlabel("orbits"); a2.set_ylabel("Gamma/Gamma0"); a2.legend(); a2.set_title("Type-I torque, q = 1e-4")
        savefig(fig, "planet_torque")
    return ok, f"Gamma/Gamma0 = {meas:+.3f} +- {osc:.3f} vs {expected:+.3f} (tol 15%, oscillation < 10%)"


# ============================================================================= 1D standard tests (SETUP=tests)
DIRS = (1, 2)                                              # 1D problems run along x and along y
DIRNAME = {1: "x", 2: "y"}
BIG = 1e12                                                 # transverse box length (negligible effect on dt)

def oned(n, direction, ndust, extra=""):
    nx, ny = (n, 1) if direction == 1 else (1, n)
    return f"SETUP=tests GEOM=0 NX={nx} NY={ny} NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=2 SOLVERDUST=0 NDUST={ndust} {extra}"

def dirsetup(direction, L, **kw):
    d = {"direction": direction, "box_l": L if direction == 1 else BIG*L, "box_l_y": L if direction == 2 else BIG*L}
    d.update(kw); return d

_VARMAP = {2: {"x": "y", "v": "vy", "vd": "vdy"}}

def profile(o, var, direction=1, idust=0):
    """1D profile along the direction of the problem (dust species idust for dust fields)."""
    name = _VARMAP.get(direction, {}).get(var, var)
    i = info(o); n = int(i["NX"])*int(i["NY"])
    return np.fromfile(os.path.join(o, name))[idust*n:(idust + 1)*n]

def symmetry(a, b):
    """max |a - b| / max |a| between the x and y runs."""
    return np.max(np.abs(np.asarray(a) - np.asarray(b)))/max(np.max(np.abs(a)), 1e-300)

def sod_exact(x, t, x0=0.5, g=1.4, l=(1.0, 0.0, 1.0), r=(0.125, 0.0, 0.1)):
    """Exact Riemann solution (Toro 2009, ch. 4). l, r = (rho, u, p). Returns rho, u, p at x."""
    (rl, ul, pl), (rr, ur, pr) = l, r
    cl, cr = np.sqrt(g*pl/rl), np.sqrt(g*pr/rr)
    def f(p, rk, pk, ck):
        if p > pk:
            A, B = 2.0/((g + 1)*rk), (g - 1)/(g + 1)*pk
            return (p - pk)*np.sqrt(A/(p + B))
        return 2*ck/(g - 1)*((p/pk)**((g - 1)/(2*g)) - 1)
    lo, hi = 1e-10, 10*max(pl, pr)
    for _ in range(200):                                      # bisection on p*
        pm = 0.5*(lo + hi)
        if f(pm, rl, pl, cl) + f(pm, rr, pr, cr) + ur - ul > 0: hi = pm
        else: lo = pm
    ps = 0.5*(lo + hi); us = 0.5*(ul + ur) + 0.5*(f(ps, rr, pr, cr) - f(ps, rl, pl, cl))
    rho, u, p = np.empty_like(x), np.empty_like(x), np.empty_like(x)
    for i, xi in enumerate(x):
        s = (xi - x0)/t
        if s < us:                                            # left of contact
            if ps > pl:                                       # left shock
                rs = rl*((ps/pl + (g - 1)/(g + 1))/((g - 1)/(g + 1)*ps/pl + 1))
                sl = ul - cl*np.sqrt((g + 1)/(2*g)*ps/pl + (g - 1)/(2*g))
                rho[i], u[i], p[i] = (rl, ul, pl) if s < sl else (rs, us, ps)
            else:                                             # left rarefaction
                rs = rl*(ps/pl)**(1/g); cs_ = cl*(ps/pl)**((g - 1)/(2*g))
                if s < ul - cl: rho[i], u[i], p[i] = rl, ul, pl
                elif s > us - cs_: rho[i], u[i], p[i] = rs, us, ps
                else:
                    u[i] = 2/(g + 1)*(cl + (g - 1)/2*ul + s); c = 2/(g + 1)*(cl + (g - 1)/2*(ul - s))
                    rho[i] = rl*(c/cl)**(2/(g - 1)); p[i] = pl*(c/cl)**(2*g/(g - 1))
        else:                                                 # right of contact
            if ps > pr:                                       # right shock
                rs = rr*((ps/pr + (g - 1)/(g + 1))/((g - 1)/(g + 1)*ps/pr + 1))
                sr = ur + cr*np.sqrt((g + 1)/(2*g)*ps/pr + (g - 1)/(2*g))
                rho[i], u[i], p[i] = (rr, ur, pr) if s > sr else (rs, us, ps)
            else:                                             # right rarefaction
                rs = rr*(ps/pr)**(1/g); cs_ = cr*(ps/pr)**((g - 1)/(2*g))
                if s > ur + cr: rho[i], u[i], p[i] = rr, ur, pr
                elif s < us + cs_: rho[i], u[i], p[i] = rs, us, ps
                else:
                    u[i] = 2/(g + 1)*(-cr + (g - 1)/2*ur + s); c = 2/(g + 1)*(cr - (g - 1)/2*(ur - s))
                    rho[i] = rr*(c/cr)**(2/(g - 1)); p[i] = pr*(c/cr)**(2*g/(g - 1))
    return rho, u, p

def linear_evolve(M, c0, t):
    """c(t) = exp(M t) c0 for a diagonalisable matrix."""
    lam, V = np.linalg.eig(M)
    return V @ (np.exp(lam*t)*np.linalg.solve(V, c0))

def dustybox_exact(t, vg0, vd0, eps, ts):
    n = len(eps); M = np.zeros((n + 1, n + 1))
    for i in range(n):
        M[0, 0] -= eps[i]/ts[i]; M[0, i + 1] += eps[i]/ts[i]
        M[i + 1, 0] = 1/ts[i]; M[i + 1, i + 1] = -1/ts[i]
    return np.array([linear_evolve(M, np.r_[vg0, vd0].astype(complex), tt).real for tt in np.atleast_1d(t)])

def dustywave_exact(x, t, amp, cs, eps, ts, k):
    """Linear gas + N dust species wave, IC: all fields ~ sin(kx) (see setup.f90). Returns rho_g, v_g, rho_d[i], v_d[i]."""
    n = len(eps); N = 2 + 2*n; M = np.zeros((N, N), complex); I = 1j
    M[0, 1] = -I*k                                            # d rho_g = -ik v_g      (rho_g0 = 1)
    M[1, 0] = -I*k*cs**2                                      # d v_g   = -ik cs^2 rho_g + sum eps (v_d - v_g)/ts
    for i in range(n):
        M[1, 1] -= eps[i]/ts[i]; M[1, 3 + 2*i] += eps[i]/ts[i]
        M[2 + 2*i, 3 + 2*i] = -I*k*eps[i]                     # d rho_d = -ik eps v_d
        M[3 + 2*i, 1] = 1/ts[i]; M[3 + 2*i, 3 + 2*i] = -1/ts[i]
    c0 = -I*amp*np.r_[1.0, cs, np.ravel([[e, cs] for e in eps])]   # sin(kx) = Re[-i e^{ikx}]
    c = linear_evolve(M, c0, t)
    return [np.real(cj*np.exp(I*k*x)) for cj in c]

def dustyshock_exact(cs=1.0, K=1.0, vs=2.0, eps=1.0, L=40.0, dx=1e-3):
    """Steady isothermal dusty shock at x=0 (dust continuous, gas jumps); returns x, vg, vd (x >= 0), vdn."""
    mg, md = vs, eps*vs; Pi = mg*vs + cs**2 + md*vs
    vdn = (Pi - np.sqrt(Pi**2 - 4*(mg + md)*cs**2*mg))/(2*(mg + md))
    vg = lambda vd: ((Pi - md*vd) - np.sqrt((Pi - md*vd)**2 - 4*mg**2*cs**2))/(2*mg)
    f = lambda vd: K*(vg(vd) - vd)/md
    xs, vds, vd = [0.0], [vs], vs
    for _ in range(int(L/dx)):
        k1 = f(vd); k2 = f(vd + dx/2*k1); k3 = f(vd + dx/2*k2); k4 = f(vd + dx*k3)
        vd = max(vd + dx/6*(k1 + 2*k2 + 2*k3 + k4), vdn); xs.append(xs[-1] + dx); vds.append(vd)
    xs, vds = np.array(xs), np.array(vds)
    return xs, vg(vds), vds, vdn

def t_sod():
    err, prof = {}, {}
    for dr in DIRS:
        d, txt = run(f"sod_{DIRNAME[dr]}", os.path.join(SUITE, "SOD", "sod.nml"), oned(256, dr, 0), threads=1,
                     setup=dirsetup(dr, 1.0))
        o = outputs(d)[-1]; t = info(o)["time"]
        x, rho, v, p = (profile(o, k, dr) for k in ("x", "rho", "v", "P"))
        re, ue, pe = sod_exact(x, t)
        err[dr] = {k: np.mean(np.abs(a - b)) for k, a, b in (("rho", rho, re), ("v", v, ue), ("P", p, pe))}
        prof[dr] = (x, rho, v, p, re, ue, pe, t)
    sym = max(symmetry(prof[1][j], prof[2][j]) for j in (1, 2, 3))
    if plt:
        x, rho, v, p, re, ue, pe, t = prof[1]
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        for a, j, ex, lab in zip(ax, (1, 2, 3), (re, ue, pe), ("rho", "v", "P")):
            a.plot(x, prof[1][j], "o", ms=2, label="SHARK, x"); a.plot(prof[2][0], prof[2][j], "+", ms=3, label="SHARK, y")
            a.plot(x, ex, "k-", lw=1, label="exact"); a.set_xlabel("x or y")
            a.set_title(f"{lab}  (L1 = {err[1][lab]:.1e} x, {err[2][lab]:.1e} y)"); a.legend()
        fig.suptitle(f"Sod shock tube, 256 cells, t = {t:.3f}, x/y difference {sym:.1e}"); savefig(fig, "sod")
    ok = max(err[dr]["rho"] for dr in DIRS) < 1e-2 and sym < 1e-8
    return ok, (f"L1 rho: {err[1]['rho']:.1e} (x), {err[2]['rho']:.1e} (y) (tol 1e-2); x/y difference {sym:.1e} (tol 1e-8)")

def t_soundwave():
    base = os.path.join(SUITE, "SOUNDWAVE", "soundwave.nml"); amp = 1e-4; errs, curves = {}, {}
    for dr in DIRS:
        for n in (32, 64, 128):
            d, txt = run(f"soundwave_{n}_{DIRNAME[dr]}", base, oned(n, dr, 0), threads=1, setup=dirsetup(dr, 1.0))
            o = outputs(d)[-1]; t = info(o)["time"]; x = profile(o, "x", dr); rho = profile(o, "rho", dr)
            errs[(dr, n)] = np.mean(np.abs(rho - (1.0 + amp*np.sin(2*np.pi*(x - t)))))/amp
            curves[(dr, n)] = (x, (rho - 1)/amp)
    order = {dr: np.log2(errs[(dr, 64)]/errs[(dr, 128)]) for dr in DIRS}
    sym = symmetry(curves[(1, 128)][1], curves[(2, 128)][1])
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for n in (32, 64, 128): ax[0].plot(*curves[(1, n)], label=f"{n} cells (x)")
        ax[0].plot(*curves[(2, 128)], "+", ms=3, label="128 cells (y)")
        xx = np.linspace(0, 1, 400); ax[0].plot(xx, np.sin(2*np.pi*xx), "k--", label="exact"); ax[0].legend()
        ax[0].set_xlabel("x or y"); ax[0].set_ylabel("(rho - 1)/A"); ax[0].set_title("Sound wave after one period")
        ns = np.array([32, 64, 128])
        for dr, mk in ((1, "o-"), (2, "s--")): ax[1].loglog(ns, [errs[(dr, n)] for n in ns], mk, label=f"SHARK, {DIRNAME[dr]}")
        ax[1].loglog(ns, errs[(1, 32)]*(ns/32.0)**-2, "k:", label="2nd order"); ax[1].legend()
        ax[1].set_xlabel("N"); ax[1].set_ylabel("L1 error / A"); ax[1].set_title(f"Convergence, x/y difference {sym:.1e}")
        savefig(fig, "soundwave")
    ok = min(order.values()) > 1.7 and max(errs[(dr, 128)] for dr in DIRS) < 2e-2 and sym < 1e-8
    return ok, (f"orders {order[1]:.2f} (x), {order[2]:.2f} (y) (> 1.7); L1/A at 128: {errs[(1, 128)]:.1e}; "
                f"x/y difference {sym:.1e}")

def _dustybox(name, ndust, eps, ts, dr):
    base = os.path.join(SUITE, "DUSTYBOX", "dustybox.nml"); vg0 = 1.0
    arr = lambda a: ", ".join(str(v) for v in a)
    d, txt = run(f"{name}_{DIRNAME[dr]}", base, oned(32, dr, ndust), threads=1,
                 setup=dirsetup(dr, 1.0, eps_dust=arr(eps), ts_dust=arr(ts), vd0=arr([0.0]*ndust)))
    t = np.array([info(o)["time"] for o in outputs(d)])
    sim = np.array([[profile(o, "v", dr).mean()] + [profile(o, "vd", dr, i).mean() for i in range(ndust)]
                    for o in outputs(d)])
    ex = dustybox_exact(t, vg0, np.zeros(ndust), eps, ts)
    return t, sim, ex, np.abs(sim - ex).max()/vg0

def t_dustybox():
    cases = {"dustybox_1": (1, [1.0], [1.0]), "dustybox_3": (3, [0.5, 1.0, 1.5], [0.5, 1.0, 2.0])}
    res = {(name, dr): _dustybox(name, *c, dr) for name, c in cases.items() for dr in DIRS}
    sym = max(symmetry(res[(name, 1)][1], res[(name, 2)][1]) for name in cases)
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for a, name in zip(ax, cases):
            t, sim, ex, err = res[(name, 1)]; nd, eps, ts = cases[name]
            tt = np.linspace(0, t.max(), 300); ee = dustybox_exact(tt, 1.0, np.zeros(nd), eps, ts)
            for j in range(sim.shape[1]):
                a.plot(t, sim[:, j], "o", ms=3, label=("gas" if j == 0 else f"dust {j}") + " SHARK"); a.plot(tt, ee[:, j], "k-", lw=0.8)
            a.set_xlabel("t"); a.set_ylabel("v"); a.legend(fontsize=7)
            a.set_title(f"{name}: max error {err:.1e} (x), {res[(name, 2)][3]:.1e} (y)")
        fig.suptitle(f"Dustybox (Laibe & Price 2011), black: exact, x/y difference {sym:.1e}"); savefig(fig, "dustybox")
    errmax = max(r[3] for r in res.values())
    return errmax < 2e-2 and sym < 1e-8, f"max |v - v_exact|/v0 = {errmax:.1e} over 1/3 species, x and y (tol 2e-2); x/y difference {sym:.1e}"

def t_dustywave():
    amp, eps, ts = 1e-4, [1.0], [0.1]; err, sims = {}, {}
    for dr in DIRS:
        d, txt = run(f"dustywave_{DIRNAME[dr]}", os.path.join(SUITE, "DUSTYWAVE", "dustywave.nml"), oned(128, dr, 1),
                     threads=1, setup=dirsetup(dr, 1.0))
        o = outputs(d)[-1]; t = info(o)["time"]; x = profile(o, "x", dr)
        ex = dustywave_exact(x, t, amp, 1.0, eps, ts, 2*np.pi)
        sim = [profile(o, "rho", dr) - 1.0, profile(o, "v", dr), profile(o, "rhod", dr) - eps[0], profile(o, "vd", dr)]
        err[dr] = [np.mean(np.abs(a - b))/amp for a, b in zip(sim, ex)]; sims[dr] = (x, sim, ex, t)
    sym = max(symmetry(sims[1][1][j], sims[2][1][j]) for j in range(4))
    if plt:
        x, sim, ex, t = sims[1]
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for a, js, ttl in ((ax[0], (0, 2), "density perturbation / A"), (ax[1], (1, 3), "velocity / (A cs)")):
            for j, lab in zip(js, ("gas", "dust")):
                a.plot(x, sim[j]/amp, "o", ms=2, label=f"{lab} SHARK, x"); a.plot(sims[2][0], sims[2][1][j]/amp, "+", ms=3, label=f"{lab} SHARK, y")
                a.plot(x, ex[j]/amp, "k-", lw=0.8)
            a.set_title(ttl); a.set_xlabel("x or y"); a.legend(fontsize=7)
        fig.suptitle(f"Dustywave, eps = 1, ts = 0.1, t = {t:.2f} (black: exact), x/y difference {sym:.1e}"); savefig(fig, "dustywave")
    e = max(max(v) for v in err.values())
    return e < 0.1 and sym < 1e-8, f"max L1/A over the 4 fields = {e:.1e} (x and y, tol 0.1); x/y difference {sym:.1e}"

def t_dustyshock():
    xr, vgr, vdr, vdn = dustyshock_exact(); err, prof = {}, {}
    for dr in DIRS:
        d, txt = run(f"dustyshock_{DIRNAME[dr]}", os.path.join(SUITE, "DUSTYSHOCK", "dustyshock.nml"), oned(600, dr, 1),
                     threads=1, setup=dirsetup(dr, 60.0))
        o = outputs(d)[-1]; x = profile(o, "x", dr); vg, vd = profile(o, "v", dr), profile(o, "vd", dr)
        xs = x[np.argmax(vg < 0.5*(2.0 + 0.5))]
        sel = (x > xs - 5) & (x < xs + 30)
        vg_ref = np.where(x < xs, 2.0, np.interp(x - xs, xr, vgr)); vd_ref = np.where(x < xs, 2.0, np.interp(x - xs, xr, vdr))
        err[dr] = (np.mean(np.abs(vg - vg_ref)[sel])/(2.0 - vdn), np.mean(np.abs(vd - vd_ref)[sel])/(2.0 - vdn))
        prof[dr] = (x - xs, vg, vd, vg_ref, vd_ref)
    sym = max(symmetry(prof[1][1], prof[2][1]), symmetry(prof[1][2], prof[2][2]))
    if plt:
        xx, vg, vd, vg_ref, vd_ref = prof[1]
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(xx, vg, "o", ms=2, label="gas SHARK, x"); ax.plot(xx, vd, "s", ms=2, label="dust SHARK, x")
        ax.plot(prof[2][0], prof[2][1], "+", ms=3, label="gas SHARK, y"); ax.plot(prof[2][0], prof[2][2], "x", ms=3, label="dust SHARK, y")
        ax.plot(xx, vg_ref, "k-", lw=0.8, label="exact"); ax.plot(xx, vd_ref, "k--", lw=0.8)
        ax.set_xlim(-5, 30); ax.set_xlabel("distance from the shock"); ax.set_ylabel("v"); ax.legend(fontsize=7)
        ax.set_title(f"Dusty shock: L1 gas {err[1][0]:.1e}, dust {err[1][1]:.1e}; x/y difference {sym:.1e}"); savefig(fig, "dustyshock")
    e = max(max(v) for v in err.values())
    return e < 0.05 and sym < 1e-8, f"max L1/(v_up - v_down) = {e:.1e} (x and y, tol 5e-2); x/y difference {sym:.1e}"

def t_pscal():
    """Dust passive scalar: (1) advection with the dust at constant velocity (exact solution, 2nd-order convergence,
    conservation of rho_d s); (2) a uniform scalar stays exactly uniform in the dustywave (rho_d and v_d varying)."""
    base = os.path.join(SUITE, "PSCAL", "pscal.nml"); amp = 0.5; errs, cons, curves = {}, {}, {}
    for dr in DIRS:
        for n in (64, 128):
            d, txt = run(f"pscal_{n}_{DIRNAME[dr]}", base, oned(n, dr, 1, "NDUSTPSCAL=1"), threads=1, setup=dirsetup(dr, 1.0))
            o0, o1 = outputs(d)[0], outputs(d)[-1]; t = info(o1)["time"]; x = profile(o1, "x", dr)
            s1 = profile(o1, "dustpscal", dr); ex = 1.0 + amp*np.sin(2*np.pi*(x - t) + 1.0)
            errs[(dr, n)] = np.mean(np.abs(s1 - ex))/amp
            m0 = (profile(o0, "rhod", dr)*profile(o0, "dustpscal", dr)).sum(); m1 = (profile(o1, "rhod", dr)*s1).sum()
            cons[(dr, n)] = abs(m1/m0 - 1.0); curves[(dr, n)] = (x, s1, ex)
    order = {dr: np.log2(errs[(dr, 64)]/errs[(dr, 128)]) for dr in DIRS}
    sym = symmetry(curves[(1, 128)][1], curves[(2, 128)][1])
    # uniform scalar in the dustywave
    dev = {}
    for dr in DIRS:
        d, txt = run(f"pscal_dustywave_{DIRNAME[dr]}", os.path.join(SUITE, "DUSTYWAVE", "dustywave.nml"),
                     oned(128, dr, 1, "NDUSTPSCAL=1"), threads=1, setup=dirsetup(dr, 1.0))
        dev[dr] = np.max(np.abs(profile(outputs(d)[-1], "dustpscal", dr) - 1.0))
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for n in (64, 128): ax[0].plot(curves[(1, n)][0], curves[(1, n)][1], "o", ms=2, label=f"{n} cells (x)")
        ax[0].plot(curves[(2, 128)][0], curves[(2, 128)][1], "+", ms=3, label="128 cells (y)")
        ax[0].plot(curves[(1, 128)][0], curves[(1, 128)][2], "k-", lw=0.8, label="exact")
        ax[0].set_xlabel("x or y"); ax[0].set_ylabel("s"); ax[0].legend(); ax[0].set_title("Scalar after one crossing")
        ns = np.array([64, 128])
        for dr, mk in ((1, "o-"), (2, "s--")): ax[1].loglog(ns, [errs[(dr, n)] for n in ns], mk, label=f"SHARK, {DIRNAME[dr]}")
        ax[1].loglog(ns, errs[(1, 64)]*(ns/64.0)**-2, "k:", label="2nd order"); ax[1].legend()
        ax[1].set_xlabel("N"); ax[1].set_ylabel("L1 error / A")
        ax[1].set_title(f"uniform s in dustywave: max|s-1| = {max(dev.values()):.1e}")
        savefig(fig, "pscal")
    ok = (min(order.values()) > 1.7 and max(cons.values()) < 1e-12 and max(dev.values()) < 1e-12 and sym < 1e-8)
    return ok, (f"advection: orders {order[1]:.2f} (x), {order[2]:.2f} (y) (> 1.7), L1/A at 128 {errs[(1, 128)]:.1e}, "
                f"rho_d s conserved to {max(cons.values()):.0e} (< 1e-12); uniform s in dustywave: max|s-1| = "
                f"{max(dev.values()):.0e} (< 1e-12); x/y difference {sym:.1e}")

# ----------------------------------------------------------------------------- Smoluchowski kernels
_trap = getattr(np, "trapezoid", None) or np.trapz

def n_constant(m, tau):
    """Exact solution, K = 1, n(m,0) = exp(-m) (N0 = m0 = 1), tau = t."""
    return 4.0/(2.0 + tau)**2*np.exp(-2.0*m/(2.0 + tau))

def n_additive(m, tau):
    """Exact solution (Golovin), K = m1 + m2, n(m,0) = exp(-m), tau = t, T = 1 - exp(-tau)."""
    T = 1.0 - np.exp(-tau); th = np.linspace(0.0, np.pi, 2001); m = np.atleast_1d(m)
    out = np.empty_like(m, dtype=float)
    for k in range(0, len(m), 500):                       # exp(-m(1+T)) I1(2 m sqrt T), stable integral form
        mm = m[k:k + 500, None]
        out[k:k + 500] = _trap(np.exp(-mm*(1 + T) + 2*mm*np.sqrt(T)*np.cos(th))*np.cos(th), th, axis=1)/np.pi
    return (1 - T)*out/(m*np.sqrt(T))

def _smoluchowski(kernel):
    nd, mmin, mmax = 40, 1e-4, 1e6
    tend, box = (100.0, 1.0) if kernel == 1 else (2.0, 0.1)
    name = "kernel_constant" if kernel == 1 else "kernel_additive"
    d, txt = run(name, os.path.join(SUITE, "SMOLUCHOWSKI", "smoluchowski.nml"),
                 f"SETUP=tests GEOM=0 NX=1 NY=1 NDUST={nd} NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=2 SOLVERDUST=0", threads=1,
                 grid={"tend": tend}, setup={"box_l": box}, dust={"kernel_type": kernel})
    eta = (mmax/mmin)**(1.0/nd)
    m = mmin*eta**(np.arange(nd) + 0.5)                     # representative masses (geometric bin centres)
    t, rho = [], []
    for o in outputs(d):
        t.append(info(o)["time"]); rho.append(np.fromfile(os.path.join(o, "rhod"))[:nd])
    t, rho = np.array(t), np.array(rho)
    meanm = (rho*m).sum(1)/rho.sum(1)                        # mass-weighted mean mass M2/M1
    exact = (2.0 + t) if kernel == 1 else 2.0*np.exp(2.0*t)
    exact *= meanm[0]/2.0                                    # same discretised initial condition
    mass_err = abs(rho[-1].sum()/rho[0].sum() - 1.0); err = abs(meanm[-1]/exact[-1] - 1.0)
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        nfun = n_constant if kernel == 1 else n_additive
        mm = np.geomspace(mmin, mmax, 600)
        snaps = (0.0, 1.0, 10.0, 100.0) if kernel == 1 else (0.0, 0.5, 1.0, 2.0)
        for k, ts in enumerate(snaps):
            i = np.argmin(np.abs(t - ts)); c = f"C{k}"
            ax[0].loglog(m, rho[i]/np.log(eta), "o", ms=3, color=c, label=f"t = {t[i]:.1f}")
            ex = mm**2*(np.exp(-mm) if t[i] == 0 else nfun(mm, t[i])); ax[0].loglog(mm, ex, "-", color=c, lw=0.8)
        ax[0].set_ylim(1e-6, 2); ax[0].set_xlabel("m"); ax[0].set_ylabel("m^2 n(m)  (mass per ln m)"); ax[0].legend()
        ax[0].set_title(("Constant" if kernel == 1 else "Additive") + " kernel, 40 bins (lines: exact)")
        ax[1].semilogy(t, meanm, "o", ms=3, label="SHARK"); ax[1].semilogy(t, exact, "k-", label="exact")
        ax[1].set_xlabel("t"); ax[1].set_ylabel("mass-weighted mean mass"); ax[1].legend()
        ax[1].set_title(f"error {100*err:.1f}% at t = {t[-1]:.1f}, mass conservation {mass_err:.1e}")
        savefig(fig, name)
    return err < 0.1 and mass_err < 1e-6, f"<m>_M error {100*err:.1f}% at t = {t[-1]:.1f} (tol 10%), mass conservation {mass_err:.1e} (tol 1e-6)"

def t_kernel_constant():
    return _smoluchowski(1)

def t_kernel_additive():
    return _smoluchowski(2)

# ----------------------------------------------------------------------------- fragmentation barrier
def t_fragmentation():
    """Coagulation + Maxwellian fragmentation with turbulent dv ~ sqrt(St) (Epstein t_s ~ s):
    steady state, mass conservation, and the fragmentation-barrier scaling a_frag ~ v_frag^2 (Birnstiel+12)."""
    nd, smin, smax = 40, 1e-5, 10.0
    zeta = (smax/smin)**(1.0/nd); s = smin*zeta**(np.arange(nd) + 0.5)
    res = {}
    for vf in (100.0, 200.0):
        d, txt = run(f"fragmentation_v{int(vf)}", os.path.join(SUITE, "FRAGMENTATION", "fragmentation.nml"),
                     f"SETUP=tests GEOM=0 NX=1 NY=1 NDUST={nd} NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=2 SOLVERDUST=0",
                     threads=1, dust={"vfrag": vf})
        t, rho = [], []
        for o in outputs(d):
            t.append(info(o)["time"]); rho.append(np.fromfile(os.path.join(o, "rhod"))[:nd])
        t, rho = np.array(t), np.array(rho)
        smean = (rho*s).sum(1)/rho.sum(1)
        half = np.argmin(np.abs(t - 0.5*t[-1]))
        res[vf] = dict(t=t, rho=rho, smean=smean, steady=abs(smean[-1]/smean[half] - 1.0),
                       mass=abs(rho[-1].sum()/rho[0].sum() - 1.0))
    slope = np.log(res[200.0]["smean"][-1]/res[100.0]["smean"][-1])/np.log(2.0)
    steady = max(r["steady"] for r in res.values()); mass = max(r["mass"] for r in res.values())
    # analytic barrier of the dv model: dv(a,a) = sqrt(alpha f(1) St) cs = v_frag, St = t_s/t_L (Ormel & Cuzzi 07)
    rhog, cs, alpha, rgr, G = 1e-8, 1e5, 1e-2, 2.7, 6.67e-8
    tL = np.sqrt(3*np.pi/(32*G*rhog)); cts = np.sqrt(np.pi/8)*rgr/(rhog*cs); f1 = 3.2 - 2.0 + (2.0/2.6)
    afrag = {vf: (vf/cs)**2*tL/(alpha*f1*cts) for vf in res}
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        for k, vf in enumerate(res):
            r = res[vf]; c = f"C{k}"
            ax[0].loglog(s, r["rho"][-1]/r["rho"][-1].sum(), "o-", ms=3, color=c, label=f"v_frag = {vf:.0f} cm/s")
            ax[0].axvline(afrag[vf], color=c, ls="--", lw=0.8)
            ax[1].semilogx(r["t"], r["smean"], color=c, label=f"v_frag = {vf:.0f} cm/s")
        ax[0].set_ylim(1e-8, 1); ax[0].set_xlabel("s [cm]"); ax[0].set_ylabel("mass fraction per bin")
        ax[0].set_title("Steady state (dashed: a where dv(a,a) = v_frag)"); ax[0].legend()
        ax[1].set_yscale("log"); ax[1].set_xlabel("t [s]"); ax[1].set_ylabel("mass-weighted mean size [cm]")
        ax[1].set_title(f"<s> ~ v_frag^{slope:.2f} (barrier: 2)"); ax[1].legend()
        savefig(fig, "fragmentation")
    ok = 1.5 < slope < 2.5 and steady < 0.05 and mass < 1e-6
    return ok, (f"<s>_M = {res[100.0]['smean'][-1]:.2e}, {res[200.0]['smean'][-1]:.2e} cm: scaling v_frag^{slope:.2f} "
                f"(expected 2, tol [1.5, 2.5]); steady {100*steady:.1f}% (< 5%); mass {mass:.1e} (< 1e-6)")

# ----------------------------------------------------------------------------- disk (vertically integrated) kernel
def _kernel_grid(nd=40, mmin=1e-4, mmax=1e6):
    eta = (mmax/mmin)**(1.0/nd); mm = mmin*eta**np.arange(nd); mp = mm*eta
    return eta, mmin*eta**(np.arange(nd) + 0.5), mm, mp

def _constant_kernel_rates(rho, m, eta, vfac):
    """d rho_i/dt of dust_growth_shark for K = 1 (no fragmentation), same bin placement as SHARK."""
    nd = len(m); d = np.zeros(nd)
    for i in range(nd):
        for j in range(i + 1):
            dn = rho[i]*rho[j]/m[i]/m[j]*vfac[i, j]*(0.5 if i == j else 1.0)
            M = m[i] + m[j]
            k1 = min(max(int(np.floor(np.log(M/m[0])/np.log(eta))) + 1, 1), nd)
            while k1 > 1 and m[k1 - 1] > M: k1 -= 1
            while k1 < nd and m[min(k1, nd - 1)] <= M: k1 += 1
            k2 = min(k1 + 1, nd); e = 1.0 if k1 == k2 else min(max((m[k2 - 1] - M)/(m[k2 - 1] - m[k1 - 1]), 0), 1)
            d[i] -= m[i]*dn; d[j] -= m[j]*dn; d[k1 - 1] += e*M*dn; d[k2 - 1] += (1 - e)*M*dn
    return d

def t_disk_kernel():
    """Vertically integrated kernel (GEOM=2): dN_ij/dt = K Sigma_i Sigma_j/(m_i m_j) / (sqrt(2 pi) sqrt(H_i^2 + H_j^2)),
    H_i = H min(1, sqrt(alpha/(min(St,1/2)(1+St^2)))) (Birnstiel et al. 2010).
    (A) well-mixed dust (H_i = H): constant-kernel exact solution with tau = t/(2 sqrt(pi) H);
    (B) settled dust (H_i varies): rates after the first step vs the formula, bin by bin."""
    nd, H = 40, 1.0
    eta, m, mm, mp = _kernel_grid(nd)
    flags = f"SETUP=tests GEOM=2 NX=1 NY=1 NDUST={nd} NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=2 SOLVERDUST=0"
    base = os.path.join(SUITE, "SMOLUCHOWSKI", "smoluchowski.nml")
    setup = {"test_case": "'disk_kernel'", "h_test": H, "box_l_y": 2*np.pi}
    # (A) well mixed
    tau_end = 100.0; tend = tau_end*2*np.sqrt(np.pi)*H
    d, txt = run("disk_kernel_mixed", base, flags, threads=1, grid={"tend": f"{tend:.6f}"},
                 setup=dict(setup, box_l=1.0, ts_dust=1e-6), dust={"kernel_type": 1, "alpha_turb": 1.0})
    t, rho = [], []
    for o in outputs(d):
        t.append(info(o)["time"]); rho.append(np.fromfile(os.path.join(o, "rhod"))[:nd])
    t, rho = np.array(t), np.array(rho); tau = t/(2*np.sqrt(np.pi)*H)
    meanm = (rho*m).sum(1)/rho.sum(1); exact = (2.0 + tau)*meanm[0]/2.0
    errA = abs(meanm[-1]/exact[-1] - 1.0)
    # (B) settled: one small step, measured rates vs formula
    dt_box = 1e-6; ts0, alpha = 1e-4, 1e-3
    d, txt = run("disk_kernel_settled", base, flags, threads=1, grid={"tend": "1d-6"},
                 setup=dict(setup, box_l=dt_box, ts_dust=ts0), dust={"kernel_type": 1, "alpha_turb": alpha})
    o0, o1 = outputs(d)[0], outputs(d)[1]
    r0, r1 = np.fromfile(os.path.join(o0, "rhod"))[:nd], np.fromfile(os.path.join(o1, "rhod"))[:nd]
    dt = info(o1)["time"] - info(o0)["time"]
    s = m**(1/3); St = ts0*s/s[0]                          # t_L = H/cs = 1
    Hd = H*np.minimum(1.0, np.sqrt(alpha/(np.minimum(St, 0.5)*(1 + St**2))))
    vfac = 1.0/(np.sqrt(2*np.pi)*np.sqrt(Hd[:, None]**2 + Hd[None, :]**2))
    expect = _constant_kernel_rates(r0, m, eta, vfac); meas = (r1 - r0)/dt
    sig = np.abs(expect) > 1e-6*np.abs(expect).max()
    errB = np.max(np.abs(meas[sig]/expect[sig] - 1.0))
    errB_nofac = np.max(np.abs(meas[sig]/_constant_kernel_rates(r0, m, eta, np.full((nd, nd), 1/(2*np.sqrt(np.pi)*H)))[sig] - 1.0))
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        ax[0].semilogy(tau, meanm, "o", ms=3, label="SHARK"); ax[0].semilogy(tau, exact, "k-", label="exact")
        ax[0].set_xlabel("tau = t / (2 sqrt(pi) H)"); ax[0].set_ylabel("mass-weighted mean mass"); ax[0].legend()
        ax[0].set_title(f"(A) well-mixed dust: error {100*errA:.1f}%")
        ax[1].loglog(m[sig], np.abs(meas[sig]), "o", ms=4, label="SHARK, first step")
        ax[1].loglog(m[sig], np.abs(expect[sig]), "k+", ms=8, label="Birnstiel+10 H_i")
        ax[1].set_xlabel("m"); ax[1].set_ylabel("|d Sigma_i/dt|"); ax[1].legend()
        ax2 = ax[1].twinx(); ax2.semilogx(m, Hd/H, "r--", lw=0.8); ax2.set_ylabel("H_i / H", color="r")
        ax[1].set_title(f"(B) settled dust: max rate error {errB:.1e}")
        savefig(fig, "disk_kernel")
    ok = errA < 0.1 and errB < 1e-4
    return ok, (f"(A) well mixed: <m>_M error {100*errA:.1f}% at tau = {tau[-1]:.0f} (tol 10%); "
                f"(B) settled: max rate error {errB:.1e} (tol 1e-4; {errB_nofac:.0e} without the H_i)")

SPH = "SETUP=disk_sph GEOM=3 NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=1 SOLVERDUST=0 NDUST=0"

def t_sph_hydrostatic():
    """Locally isothermal disk on the spherical (r,theta) grid stays in hydrostatic equilibrium,
    with the residual converging with resolution (64x32 -> 128x64)."""
    base = os.path.join(SUITE, "DISK_SPH", "disk_sph.nml")
    res, maps = {}, {}
    for nx, ny in [(64, 32), (128, 64)]:
        d, txt = run(f"sph_hydrostatic_{nx}", base, SPH + f" NX={nx} NY={ny}")
        o0, o1 = outputs(d)[0], outputs(d)[-1]; i = info(o0)
        h, q, M, r0 = nml_value(txt, "hoverr", 0.05), nml_value(txt, "q_temp", -1.0), nml_value(txt, "mstar", 1.0), nml_value(txt, "r0", 1.0)
        r, th = field(o0, "x"), field(o0, "y")
        R = r*np.sin(th)
        cs = h*np.sqrt(M/r0)*(R/r0)**(0.5*q)*i["unit_v"]
        drho = field(o1, "rho")/field(o0, "rho") - 1.0
        vr, vt = field(o1, "v")/cs, field(o1, "vy")/cs
        mask = (r > 4.0) & (r < 8.0) & (np.abs(th - np.pi/2) < 0.1)
        res[nx] = np.array([np.abs(vt[mask]).max(), np.abs(vr[mask]).max(), np.abs(drho[mask]).max()])
        maps[nx] = (R, r*np.cos(th), vt, drho)
    ratio = res[64]/res[128]
    ok = res[128].max() < 1e-2 and ratio.min() > 2.5
    if plt:
        R, z, vt, drho = maps[128]
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        v = np.abs(vt).max()
        im = ax[0].pcolormesh(R, z, vt, cmap="RdBu_r", vmin=-v, vmax=v, shading="nearest"); fig.colorbar(im, ax=ax[0])
        ax[0].set_xlabel("R [au]"); ax[0].set_ylabel("z [au]"); ax[0].set_title("v_theta / c_s, 128x64")
        im = ax[1].pcolormesh(R, z, drho, cmap="RdBu_r", shading="nearest"); fig.colorbar(im, ax=ax[1])
        ax[1].set_xlabel("R [au]"); ax[1].set_title("rho/rho_0 - 1, 128x64")
        lab = ["|v_theta|/c_s", "|v_r|/c_s", "|drho/rho|"]
        for k in range(3):
            ax[2].loglog([64, 128], [res[64][k], res[128][k]], "o-", label=lab[k])
        ax[2].loglog([64, 128], [res[64].max(), res[64].max()/4], "k--", label="2nd order")
        ax[2].set_xlabel("NX"); ax[2].set_ylabel("interior max"); ax[2].legend(fontsize=8); ax[2].set_title("Convergence")
        savefig(fig, "sph_hydrostatic")
    return ok, (f"128x64, R in [4,9], |z|<2H: v_th {res[128][0]:.1e} c_s, v_r {res[128][1]:.1e} c_s, drho/rho {res[128][2]:.1e}; "
                f"64->128 ratios {ratio[0]:.1f}, {ratio[1]:.1f}, {ratio[2]:.1f} (need < 1e-2 and > 2.5)")

TESTS = [("disk_ic_drift", t_disk_ic_drift, False), ("disk_equilibrium", t_disk_equilibrium, False),
         ("dust_drift_nsh", t_dust_drift_nsh, False), ("fargo_bump", t_fargo_bump, False),
         ("si_linA", t_si_linA, False), ("sod", t_sod, False), ("soundwave", t_soundwave, False),
         ("dustybox", t_dustybox, False), ("dustywave", t_dustywave, False), ("dustyshock", t_dustyshock, False), ("pscal", t_pscal, False),
         ("kernel_constant", t_kernel_constant, False), ("kernel_additive", t_kernel_additive, False),
         ("fragmentation", t_fragmentation, False), ("disk_kernel", t_disk_kernel, False),
         ("planet_torque", t_planet_torque, True),("sph_hydrostatic", t_sph_hydrostatic, False)]

# ----------------------------------------------------------------------------- main
def main():
    global THREADS
    p = argparse.ArgumentParser(description="SHARK analytical validation suite")
    p.add_argument("--long", action="store_true", help="also run the long tests (planet torque)")
    p.add_argument("--only", nargs="+", help="run only these tests")
    p.add_argument("--threads", type=int, default=4)
    a = p.parse_args(); THREADS = a.threads
    os.makedirs(WORK, exist_ok=True)
    todo = [t for t in TESTS if (a.only and t[0] in a.only) or (not a.only and (a.long or not t[2]))]
    results = []
    for name, fn, _ in todo:
        print(f"[{name}] running...", flush=True); t0 = time.time()
        try:
            ok, msg = fn()
        except Exception as e:
            ok, msg = False, f"ERROR: {e}"
        results.append((name, ok, msg, time.time() - t0))
        print(f"[{name}] {'PASS' if ok else 'FAIL'} ({time.time() - t0:.0f} s): {msg}", flush=True)
    print("\n" + "="*100)
    for name, ok, msg, dt in results:
        print(f"{'PASS' if ok else 'FAIL'}  {name:18s} {dt:6.0f} s  {msg}")
    nfail = sum(not r[1] for r in results)
    print("="*100 + f"\n{len(results) - nfail}/{len(results)} passed. Runs in {WORK}" + (f", figures in {FIGDIR}" if plt else ", no figures (matplotlib missing)"))
    sys.exit(nfail)

if __name__ == "__main__":
    main()
