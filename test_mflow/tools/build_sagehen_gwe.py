#!/usr/bin/env python3
"""Build ex-gwf-sagehen-gwe: ex-gwf-sagehen-sfr with a GWE model on the same grid.

GWE carries the aquifer temperature by advection and conduction (ADV, CND, EST) and along
the SFR reaches (SFE). Recharge brings SUMMA's drainage temperature in through RCH's
auxiliary TEMPERATURE, which the coupler overwrites every step and SSM names as the
recharge source, and through UZE's INFILTRATION temperature on the UZF cells, which the
coupler overwrites likewise. The coupler returns GWE's water-table temperature to SUMMA as
scalarAquiferTemp, and loads ESL with SUMMA's conduction out the soil base, one row per
active cell. GWE runs in degrees Celsius, in the GWF model's metres and seconds.

The starting temperatures are a thermal equilibrium: a stand-alone run of this model with
steady flow, recharge at the long-term recharge temperature of each cell's elevation and a
basal geothermal flux, carried for EQ_YEARS. It needs mf6 (MF6, default the repository's).

Usage:
    build_sagehen_gwe.py [<output model dir>]      (default ../ex-gwf-sagehen-gwe)
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
SFR = os.path.join(HERE, "..", "ex-gwf-sagehen-sfr")
GWE = "SAGEHEN_GWE"
MF6 = os.environ.get("MF6", os.path.join(HERE, *[".."] * 6, "bin", "mf6"))

STRT = 6.0            # degC, placeholder for the equilibration's own start
T_REF = 5.2           # degC, 1991-2020 mean air temperature at Sagehen Creek Field Station
Z_REF = 1932.0        # m, the station's elevation
LAPSE = 6.5e-3        # degC m-1
GW_OFFSET = 1.0       # degC, recharge above mean air temperature, from snow insulating the ground
Q_GEO = 0.068         # W m-2, basal geothermal flux, Lake Tahoe mean 1.62 HFU (Henyey and Lee, 1976)
EQ_YEARS = 10000.0    # years of the equilibration run
RHO_SOLID = 2650.0    # kg m-3
CP_SOLID = 840.0      # J kg-1 degC-1
K_WATER = 0.58        # W m-1 degC-1
K_SOLID = 2.5         # W m-1 degC-1
K_BED = 2.0           # W m-1 degC-1, streambed
BED_THK = 0.5         # m, as build_sagehen_sfr.py's RBED_THK


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def write_array(path, a):
    np.savetxt(path, a, fmt="%9.4f")


def add_esl(model_dir, active, rate, comment):
    """One ESL row per active cell at rate (W), named ESL in the GWE name file."""
    rows, cols = np.nonzero(active)
    with open(os.path.join(model_dir, "sagehen_gwe.esl"), "w") as f:
        f.write(f"# {comment}\nBEGIN DIMENSIONS\n  MAXBOUND  {rows.size}\nEND DIMENSIONS\n\nBEGIN PERIOD 1\n")
        for r, c in zip(rows, cols):
            f.write(f"  1 {r + 1:3d} {c + 1:3d} {rate:.4f}\n")
        f.write("END PERIOD\n")
    nam = os.path.join(model_dir, "sagehen_gwe.nam")
    with open(nam) as f:
        text = f.read()
    write(nam, text.replace("  OC6 ", "  ESL6  sagehen_gwe.esl  esl\n  OC6 "))


def sfr_cells(model_dir):
    """(row, col), zero-based, of each SFR reach in order."""
    cells, inblock = [], False
    with open(os.path.join(model_dir, "sagehen.sfr")) as f:
        for line in f:
            if line.strip().startswith("BEGIN PACKAGEDATA"):
                inblock = True
            elif line.strip().startswith("END PACKAGEDATA"):
                break
            elif inblock and line.strip() and not line.lstrip().startswith("#"):
                w = line.split()
                cells.append((int(w[2]) - 1, int(w[3]) - 1))
    return cells


def read_arrays(path, nrow, ncol):
    """Every array in a single-layer MODFLOW 6 binary head or temperature file, after its 52-byte header."""
    raw = open(path, "rb").read()
    rec = 52 + 8 * nrow * ncol
    return [np.frombuffer(raw[i + 52:i + rec], dtype="<f8").reshape(nrow, ncol) for i in range(0, len(raw), rec)]


def equilibrate(model_dir, nrow, ncol):
    """Run model_dir to thermal equilibrium in a scratch copy and return (temperature, change in the last step)."""
    active = np.loadtxt(os.path.join(model_dir, "idomain1.txt")) > 0
    with tempfile.TemporaryDirectory() as tmp:
        run = os.path.join(tmp, "eq")
        shutil.copytree(model_dir, run)
        write(os.path.join(run, "sagehen.tdis"),
              f"BEGIN OPTIONS\n  TIME_UNITS  SECONDS\nEND OPTIONS\n\nBEGIN DIMENSIONS\n  NPER  1\nEND DIMENSIONS\n\n"
              f"BEGIN PERIODDATA\n  {EQ_YEARS * 365.25 * 86400.0:.6e}  60  1.1\nEND PERIODDATA\n")
        sto = os.path.join(run, "sagehen.sto")
        with open(sto) as f:
            text = f.read()
        write(sto, text.replace("  TRANSIENT\n", "  STEADY-STATE\n"))
        add_esl(run, active, Q_GEO * 90.0 * 90.0, "Basal geothermal flux (W per 90 m cell).")
        write(os.path.join(run, "sagehen_gwe.oc"),
              "BEGIN OPTIONS\n  TEMPERATURE  FILEOUT  sagehen_gwe.ucn\nEND OPTIONS\n\n"
              "BEGIN PERIOD 1\n  SAVE  TEMPERATURE  STEPS 59 60\nEND PERIOD\n")
        subprocess.run([MF6], cwd=run, check=True, stdout=subprocess.DEVNULL)
        prev, last = read_arrays(os.path.join(run, "sagehen_gwe.ucn"), nrow, ncol)
    return np.where(active, last, 0.0), np.abs(last - prev)[active].max()


def uzf_cells(model_dir):
    """(row, col), zero-based, of each UZF cell in order."""
    cells, inblock = [], False
    with open(os.path.join(model_dir, "sagehen.uzf")) as f:
        for line in f:
            if line.strip().startswith("BEGIN PACKAGEDATA"):
                inblock = True
            elif line.strip().startswith("END PACKAGEDATA"):
                break
            elif inblock and line.strip() and not line.lstrip().startswith("#"):
                w = line.split()
                cells.append((int(w[2]) - 1, int(w[3]) - 1))
    return cells


def write_uze(out_dir, strt, tinf):
    with open(os.path.join(out_dir, "sagehen_gwe.uze"), "w") as f:
        f.write("# Energy in the UZF unsaturated zone. The coupler overwrites INFILTRATION (TEMPINFL, degC)\n"
                "# with SUMMA's drainage temperature.\n")
        f.write("BEGIN OPTIONS\n  FLOW_PACKAGE_NAME  UZF\n  SAVE_FLOWS\nEND OPTIONS\n\nBEGIN PACKAGEDATA\n# uzfno strt\n")
        for n, t in enumerate(strt, 1):
            f.write(f"  {n:4d} {t:.4f}\n")
        f.write("END PACKAGEDATA\n\nBEGIN PERIOD 1\n")
        for n, t in enumerate(tinf, 1):
            f.write(f"  {n:4d} INFILTRATION {t:.4f}\n")
        f.write("END PERIOD\n")


def write_sfe(out_dir, strt):
    with open(os.path.join(out_dir, "sagehen_gwe.sfe"), "w") as f:
        f.write("# Streamflow energy along the SFR reaches, with conduction through the streambed.\n")
        f.write("BEGIN OPTIONS\n  FLOW_PACKAGE_NAME  SFR\n  SAVE_FLOWS\nEND OPTIONS\n\n")
        f.write("BEGIN PACKAGEDATA\n# rno strt ktf rbthcnd\n")
        for n, t in enumerate(strt, 1):
            f.write(f"  {n:4d} {t:.4f} {K_BED} {BED_THK}\n")
        f.write("END PACKAGEDATA\n")


def main(out_dir):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    shutil.copytree(SFR, out_dir, ignore=shutil.ignore_patterns(
        "*.lst", "*.hds", "*.cbc", "*.grb", "*.bin", ".DS_Store", "mfsim.stdout"))

    with open(os.path.join(SFR, "sagehen.sfr")) as f:
        nreach = int(re.search(r"NREACHES\s+(\d+)", f.read()).group(1))

    # recharge carries a temperature
    rch = os.path.join(out_dir, "sagehen.rcha")
    with open(rch) as f:
        text = f.read()
    text = text.replace("  READASARRAYS\n", "  READASARRAYS\n  AUXILIARY  TEMPERATURE\n")
    text = text.replace("    CONSTANT  1.0e-8\n", "    CONSTANT  1.0e-8\n  TEMPERATURE\n    OPEN/CLOSE  rchtemp1.txt\n")
    top = np.loadtxt(os.path.join(out_dir, "land1.txt"))
    active = np.loadtxt(os.path.join(out_dir, "idomain1.txt")) > 0
    trch = np.where(active, T_REF + GW_OFFSET - LAPSE * (top - Z_REF), 0.0)
    write_array(os.path.join(out_dir, "rchtemp1.txt"), trch)
    write(rch, "# TEMPERATURE (degC) is the long-term recharge temperature at each cell's elevation.\n"
               "# The coupler also overwrites TEMPERATURE (__INPUT__/SAGEHEN/RCHA/AUXVAR, degC) with SUMMA's\n"
               "# drainage temperature; GWE's SSM reads it as the recharge temperature.\n" + text)

    write(os.path.join(out_dir, "mfsim.nam"), f"""\
# MODFLOW 6 simulation name file - Sagehen with SFR and GWE, configured for the SUMMA coupler
#   * GWF model "SAGEHEN" and GWE model "{GWE}" on the same DIS grid
#   * TDIS in SECONDS, MODFLOW time step = SUMMA forcing data_step (3600 s here)
#   * RCH package "RCHA" with READASARRAYS and AUXILIARY TEMPERATURE, UZF package "UZF" with UZE
# Written by tools/build_sagehen_gwe.py.
BEGIN OPTIONS
END OPTIONS

BEGIN TIMING
  TDIS6  sagehen.tdis
END TIMING

BEGIN MODELS
  GWF6  sagehen.nam      SAGEHEN
  GWE6  sagehen_gwe.nam  {GWE}
END MODELS

BEGIN EXCHANGES
  GWF6-GWE6  sagehen.gwfgwe  SAGEHEN  {GWE}
END EXCHANGES

BEGIN SOLUTIONGROUP 1
  IMS6  sagehen.ims      SAGEHEN
  IMS6  sagehen_gwe.ims  {GWE}
END SOLUTIONGROUP 1
""")
    write(os.path.join(out_dir, "sagehen.gwfgwe"), "BEGIN OPTIONS\nEND OPTIONS\n")

    write(os.path.join(out_dir, "sagehen_gwe.nam"), """\
# GWE model name file: aquifer temperature on the GWF grid, degC.
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN PACKAGES
  DIS6  sagehen_gwe.dis  dis
  IC6   sagehen_gwe.ic   ic
  ADV6  sagehen_gwe.adv  adv
  CND6  sagehen_gwe.cnd  cnd
  EST6  sagehen_gwe.est  est
  SSM6  sagehen_gwe.ssm  ssm
  SFE6  sagehen_gwe.sfe  SFE
  UZE6  sagehen_gwe.uze  UZE
  OC6   sagehen_gwe.oc   oc
END PACKAGES
""")

    with open(os.path.join(out_dir, "sagehen.dis")) as f:
        dis = f.read()
    write(os.path.join(out_dir, "sagehen_gwe.dis"),
          dis.replace("  LENGTH_UNITS  METERS\n", "  LENGTH_UNITS  METERS\n  NOGRB\n"))

    write(os.path.join(out_dir, "sagehen_gwe.ic"),
          f"BEGIN GRIDDATA\n  STRT\n    CONSTANT  {STRT}\nEND GRIDDATA\n")
    write(os.path.join(out_dir, "sagehen_gwe.adv"), "BEGIN OPTIONS\n  SCHEME  TVD\nEND OPTIONS\n")
    write(os.path.join(out_dir, "sagehen_gwe.cnd"), f"""\
BEGIN OPTIONS
  XT3D_OFF
END OPTIONS

BEGIN GRIDDATA
  KTW
    CONSTANT  {K_WATER}
  KTS
    CONSTANT  {K_SOLID}
END GRIDDATA
""")
    write(os.path.join(out_dir, "sagehen_gwe.est"), f"""\
# POROSITY is the GWF model's specific yield.
BEGIN OPTIONS
  SAVE_FLOWS
  DENSITY_WATER        1000.0
  HEAT_CAPACITY_WATER  4184.0
END OPTIONS

BEGIN GRIDDATA
  POROSITY
    OPEN/CLOSE  sy1.txt
  HEAT_CAPACITY_SOLID
    CONSTANT  {CP_SOLID}
  DENSITY_SOLID
    CONSTANT  {RHO_SOLID}
END GRIDDATA
""")
    # SFR is an advanced package, so SFE carries it; DRN only ever takes water out
    write(os.path.join(out_dir, "sagehen_gwe.ssm"), """\
BEGIN OPTIONS
  SAVE_FLOWS
END OPTIONS

BEGIN SOURCES
  RCHA  AUX  TEMPERATURE
END SOURCES
""")
    write_sfe(out_dir, [STRT] * nreach)
    ucells = uzf_cells(out_dir)
    write_uze(out_dir, [STRT] * len(ucells), [trch[r, c] for r, c in ucells])

    write(os.path.join(out_dir, "sagehen_gwe.ims"), """\
BEGIN OPTIONS
  PRINT_OPTION  SUMMARY
END OPTIONS

BEGIN NONLINEAR
  OUTER_DVCLOSE  1.0e-4
  OUTER_MAXIMUM  50
END NONLINEAR

BEGIN LINEAR
  INNER_MAXIMUM        300
  INNER_DVCLOSE        1.0e-5
  INNER_RCLOSE         1.0e-3
  LINEAR_ACCELERATION  BICGSTAB
END LINEAR
""")
    write(os.path.join(out_dir, "sagehen_gwe.oc"), """\
BEGIN OPTIONS
  BUDGET       FILEOUT  sagehen_gwe.cbc
  TEMPERATURE  FILEOUT  sagehen_gwe.ucn
END OPTIONS

BEGIN PERIOD 1
  SAVE  TEMPERATURE  FREQUENCY 24
  SAVE  BUDGET       FREQUENCY 24
  SAVE  TEMPERATURE  LAST
  SAVE  BUDGET       LAST
END PERIOD
""")

    nrow, ncol = top.shape
    temp, change = equilibrate(out_dir, nrow, ncol)
    write_array(os.path.join(out_dir, "strt_gwe1.txt"), temp)
    write(os.path.join(out_dir, "sagehen_gwe.ic"),
          f"# Thermal equilibrium from tools/build_sagehen_gwe.py, degC: {EQ_YEARS:.0f} years of steady flow, recharge at\n"
          f"# {T_REF + GW_OFFSET} degC at {Z_REF:.0f} m less {LAPSE * 1000} degC km-1, and {Q_GEO * 1000:.0f} mW m-2 geothermal flux.\n"
          "BEGIN GRIDDATA\n  STRT\n    OPEN/CLOSE  strt_gwe1.txt\nEND GRIDDATA\n")
    write_sfe(out_dir, [temp[r, c] for r, c in sfr_cells(out_dir)])
    write_uze(out_dir, [temp[r, c] for r, c in ucells], [trch[r, c] for r, c in ucells])
    add_esl(out_dir, active, 0.0, "The coupler overwrites SENERRATE (__INPUT__/SAGEHEN_GWE/ESL/SENERRATE, W) with\n"
            "# SUMMA's conduction out the soil base, at each cell's water-table node.")
    print(f"equilibrium aquifer temperature {temp[active].min():.2f} to {temp[active].max():.2f} degC, "
          f"mean {temp[active].mean():.2f}; last step changed it by at most {change:.1e} degC")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "ex-gwf-sagehen-gwe")))
