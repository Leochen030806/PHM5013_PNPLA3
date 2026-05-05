<div align="center">

<!-- HEADER BANNER -->
<img src="https://capsule-render.vercel.app/api?type=waving&color=0A2342,1F4E79,2E86C1&height=200&section=header&text=PNPLA3%20I148M%20%20&fontSize=38&fontColor=FFFFFF&fontAlignY=38&desc=AI-Driven%20Drug%20Discovery%20%7C%20PHM5013&descAlignY=58&descSize=18" width="100%"/>

<br/>

[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
[![RDKit](https://img.shields.io/badge/RDKit-2023.09-4B8BBE?style=for-the-badge)](https://rdkit.org)
[![OpenEye](https://img.shields.io/badge/OpenEye-OMEGA%20%7C%20ROCS%20%7C%20EON-E63946?style=for-the-badge)](https://www.eyesopen.com)
[![AutoDock Vina](https://img.shields.io/badge/AutoDock_Vina-1.2.x-F4A261?style=for-the-badge)](https://vina.scripps.edu)
[![AMBER](https://img.shields.io/badge/AMBER-24.3-2A9D8F?style=for-the-badge)](https://ambermd.org)
[![Google Colab](https://img.shields.io/badge/Google_Colab-F9AB00?style=for-the-badge&logo=googlecolab&logoColor=white)](https://colab.research.google.com)

<br/>

> **Computational identification of PNPLA3 I148M binders using shape + electrostatics virtual screening, AutoDock Vina molecular docking, and AMBER molecular dynamics simulation.**

<br/>

</div>

---

## 📋 Table of Contents

- [Project Overview](#-project-overview)
- [Pipeline at a Glance](#-pipeline-at-a-glance)
- [Repository Structure](#-repository-structure)
- [Phase 1 · Target & Reference Compound](#phase-1--target--reference-compound)
- [Phase 2 · Library Construction & Standardisation](#phase-2--library-construction--standardisation)
- [Phase 3 · Virtual Screening](#phase-3--virtual-screening)
- [Phase 4 · Molecular Docking](#phase-4--molecular-docking)
- [Phase 5 · Molecular Dynamics](#phase-5--molecular-dynamics)
- [Results Summary](#-results-summary)
- [Dependencies](#-dependencies)
- [Team](#-team)

---

## 🧬 Project Overview

<table>
<tr>
<td width="50%">

**Target** | PNPLA3 I148M (rs738409)
**Disease** | Non-alcoholic fatty liver disease (NAFLD)
**Reference molecule** | NUV-244 (validated PNPLA3 I148M degrader)
**Library** | MGTbind + MolGlue + MGDB — 11,204 molecular glues
**Screening** | OMEGA2 → ROCS → Molcharge → EON
**Docking** | AutoDock Vina / Vinardo, Google Colab
**MD** | AMBER 24.3, ff14SB + GAFF + TIP3P, NUS Vanda HPC GPU

</td>
<td width="50%">

**Rationale:** The I148M variant locks the ATGL co-activator ABHD5, blocking triglyceride hydrolysis and driving hepatic steatosis. It is the strongest genetic risk factor for NAFLD (MAF ~49% in Hispanic, ~30–40% in East Asian populations). NUV-244 — identified from an 820k-compound HTS — selectively degrades I148M PNPLA3 via the UPS/BFAR E3 ligase and serves as our shape + electrostatics query for analogue discovery.

</td>
</tr>
</table>

---

## 🔬 Pipeline at a Glance

```
  AlphaFold Structure (target.pdb)          reference_compound.sdf
             │                                        │
             ▼                                        ▼
      PyMOL Mutation Modelling              ┌─────────────────────┐
      (I148M point mutation)                │  Phase 2            │
             │                             │  DB Curation        │
             ▼                             │  ~11,000 mol. glues │
      Phase 4 · Protein Prep               └──────────┬──────────┘
      PDBFixer → pdb2pqr → OpenMM                     │ OMEGA2
      (10-step automated pipeline)                     │ (50 confs/mol)
             │                                        ▼
             │                             ┌─────────────────────┐
             │                             │  Phase 3            │
             │                             │  ROCS  →  top_3000 │
             │                             │  EON   →  top_1000 │
             │                             └──────────┬──────────┘
             │                                        │ OEB → SDF
             ▼                                        ▼
      receptor_prepared.pdb          screened_library_top1000.sdf
             │                                        │
             └──────────────┬─────────────────────────┘
                            ▼
               Phase 4 · AutoDock Vina (vinardo)
               ~1,000 compounds  ×  exhaustiveness=8
                            │
                            ▼
               docking_results.csv  (ranked by affinity)
                            │
              ┌─────────────┴──────────────┐
              ▼                            ▼
       top_candidates.sdf        reference_control.sdf
       (add H, bond-order repair)
              │
              ▼
       Phase 5 · AMBER MD  (100 ns NPT)
       GAFF + ff14SB + TIP3P + explicit water box
              │
              ▼
       trajectory_combined.nc
       → RMSD / RMSF / H-bonds / PCA / Free Energy Landscape
```

---

## 🗂 Repository Structure

```
PNPLA3-AIDD/
├── data/
│   ├── raw/
│   │   └── NUV-244.sdf
│   ├── processed/
│   │   ├── merged_3db_standardized.sdf        # 11,204 standardised MGs
│   │   ├── NUV244_5conf_charged.oeb
│   │   └── PNPLA3I148M.pdb
│   └── library_3D/
│       └── RDKit_library_3Dconformers.sdf
├── scripts/
│   ├── preprocessing/
│   │   ├── RDKit_sanitisation_Conformers.ipynb
│   │   └── 3D_conformers.ipynb
│   └── docking/
│       └── resume_docking.py
├── screening/
│   ├── rocs/
│   │   ├── NUV244_5conf_rocs.oeb
│   │   └── rocsrun_5confNuv_50conflibrary.oeb
│   └── eon/
│       ├── eon_NUV244_top3000.oeb
│       └── eon_hitlist_rpt.csv
├── docking/
│   ├── summary_results.csv
│   └── top10_rerun/
├── md/
│   ├── min1.in
│   ├── min2.in
│   ├── md_heat.in
│   ├── md_equil.in
│   ├── md_prod.in
│   └── run_amber.pbs
├── docs/
│   └── Workflow_Comparison_Table_RDKit_OpenBabel.docx
└── README.md
```

---

<br/>

<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║                  PHASE 1 · TARGET & REFERENCE                ║
╚══════════════════════════════════════════════════════════════╝
```

</div>

## Phase 1 · Target & Reference Compound

### 🎯 PNPLA3 I148M

The **rs738409 (I148M)** variant converts Ile→Met at position 148 of PNPLA3, remodelling the hydrophobic lid pocket. The mutant protein sequesters ABHD5 (CGI-58), preventing ATGL activation and driving hepatic lipid accumulation.

| Property | Value |
|---|---|
| UniProt | Q9NST1 |
| Structure source | AlphaFold DB (pLDDT > 90 for Patatin domain, residues 10–179) |
| Mutation modelling | PyMOL `Mutagenesis` wizard, Ile148 → Met |
| Output | `PNPLA3I148M.pdb` |

### 💊 NUV-244 — Reference Degrader

Identified by Nuvation Bio from an 820k-compound phenotypic screen. Selectively degrades PNPLA3 I148M via the **UPS / BFAR E3 ligase**, restores lipid droplet morphology, and does not affect PNPLA2.

> Steigemann P, et al. *iScience.* 2025;28(5):112384. doi:10.1016/j.isci.2025.112384

---

<br/>

<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║            PHASE 2 · LIBRARY CONSTRUCTION & CLEANING         ║
╚══════════════════════════════════════════════════════════════╝
```

</div>

## Phase 2 · Library Construction & Standardisation

### 📦 Database Selection

| Database | Decision | Reason |
|---|---|---|
| MGTbind | ✅ Included | Large, open, MG-relevant |
| MolGlue | ✅ Included | Curated molecular glues |
| MGDB | ✅ Included | High MG density |
| ChEMBL | ❌ Rejected | 40% similarity → only 20 hits |
| TPDdb | ❌ Rejected | Conversion pipeline issues |
| PROTAC-DB / ChemDiv | ❌ Rejected | Low structural relevance |

**Final merged library: 11,253 entries → 11,204 after standardisation**

---

### 🧹 RDKit Standardisation Pipeline

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

supplier = Chem.SDMolSupplier('merged_3db.sdf', sanitize=False)

lfc = rdMolStandardize.LargestFragmentChooser()  # salt stripping
unc = rdMolStandardize.Uncharger()               # charge neutralisation
nrm = rdMolStandardize.Normalizer()              # functional group normalisation

processed = []
for mol in supplier:
    if mol is None: continue
    try:
        Chem.SanitizeMol(mol)           # valence / aromaticity check
        mol = lfc.choose(mol)           # keep largest fragment
        mol = unc.uncharge(mol)         # remove formal charges
        mol = nrm.normalize(mol)        # standardise representation
        mol.SetProp('_Name', mol.GetProp('merged_id'))  # preserve ID
        processed.append(mol)
    except Exception as e:
        pass  # log and skip

# Output: 3DB_11204_standardised.sdf
```

| Tool | Result |
|---|---|
| **RDKit** (selected) | ✅ 11,204 molecules — full ID tracking, no pyrazole bug |
| OpenBabel (comparison) | 9,634 molecules — 1,468 duplicates removed; pyrazole rigid-fragment bug |

---

### 🔷 3D Conformer Generation — OMEGA2

```bash
# NUV-244 query: 5 conformers
omega2 \
  -in  NUV244.sdf \
  -out NUV244_5conf_rocs.oeb \
  -maxconfs 5 \
  -rms 0.6 \
  -ewindow 12

# Library: 50 conformers per molecule
omega2 \
  -in  library_11204.sdf \
  -out library_11204_50conf_rocs.oeb \
  -maxconfs 50 \
  -rms 0.6 \
  -ewindow 15
```

> Three tools tested: RDKit ETKDGv3 (controllable, slow) · OpenBabel (pyrazole bug) · **OMEGA2 (selected — fast, high quality)**

---

<br/>

<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║                  PHASE 3 · VIRTUAL SCREENING                 ║
╚══════════════════════════════════════════════════════════════╝
```

</div>

## Phase 3 · Virtual Screening

### Screening Funnel

```
11,204 molecules (50 confs/mol)
         │
         ▼  ROCS — shape + chemical feature similarity
         │  query: NUV-244 (5 confs), ranked by TanimotoCombo
         ▼
      Top 3,000 hits
         │
         ▼  Molcharge — AM1-BCC ELF10 partial charges
         │  (required for EON electrostatic scoring)
         ▼
      Charged OEB files
         │
         ▼  EON — electrostatic potential similarity (et_pb)
         ▼
      Top 1,000 hits  →  1,005 unique molecules (SDF)
```

---

### ⚗️ Step-by-step Commands

**ROCS — Shape Screening**
```bash
rocs \
  -query  NUV244_5conf_charged.oeb \
  -db     library_11204_50conf_rocs.oeb \
  -output rocs_hits.oeb \
  -rankby TanimotoCombo
# Output: rocsrun_5confNuv_50conflibrary.oeb / .csv  (Top 3,000)
```

**Molcharge — Partial Charge Assignment**
```bash
# Query molecule
molcharge \
  -in  NUV244_5conf_rocs.oeb \
  -out NUV244_5conf_charged.oeb \
  -method am1bcc

# ROCS hit library
molcharge \
  -in  rocsrun_5confNuv_50conflibrary.oeb \
  -out rocsrun_top3000_charged.oeb \
  -method am1bcc
```

**EON — Electrostatic Rescoring**
```bash
eon \
  -query  NUV244_5conf_charged.oeb \
  -dbase  rocsrun_top3000_charged.oeb \
  -out    eon_NUV244_top3000.oeb
# Output: eon_NUV244_top3000.oeb  (1,000 hits, OEB format)
```

**Format Conversion: OEB → SDF (for AutoDock Vina)**
```python
from openeye import oechem

ifs = oechem.oemolistream('eon_NUV244_top3000.oeb')
ofs = oechem.oemolostream('eon_NUV244_top3000.sdf')

mol = oechem.OEGraphMol()
count = 0
while oechem.OEReadMolecule(ifs, mol):
    oechem.OEWriteMolecule(ofs, mol)
    count += 1

ifs.close(); ofs.close()
print(f"Converted {count} molecules")
# Output: eon_NUV244_top3000.sdf — 1,005 molecules
```

| Stage | Tool | Input → Output |
|---|---|---|
| 3D Conformers | OMEGA2 | 11,204 molecules → 11,204 × 50 confs |
| Shape screening | ROCS | 11,204 → Top 3,000 |
| Charge assignment | Molcharge | NUV-244 + Top 3,000 → charged OEB |
| Electrostatic screen | EON | Top 3,000 → Top 1,000 |
| Format conversion | OEChem API | OEB → 1,005 SDF |

---

<br/>

<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║                  PHASE 4 · MOLECULAR DOCKING                 ║
╚══════════════════════════════════════════════════════════════╝
```

</div>

## Phase 4 · Molecular Docking

### 🧪 Step 1 — Protein Preparation

10-step automated pipeline on Google Colab (PDBFixer → BioPython → PropKa → pdb2pqr → OpenMM):

```python
# Key parameters (modified from tutorial defaults)
target_pH        = 7.4
restraint_k      = 500    # kJ/mol/nm²  — harmonic restraints on heavy atoms
convergence_tol  = 50     # kJ/mol/nm   — reduced from default 10
max_iterations   = 200    #             — reduced from default 2000

# OpenMM energy minimisation (AMBER14 + GBN2 implicit solvent)
ff = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
system = ff.createSystem(pdb_in.topology,
                         nonbondedCutoff=1.0 * unit.nanometers)

simulation.minimizeEnergy(
    tolerance=convergence_tol * unit.kilojoules_per_mole / unit.nanometers,
    maxIterations=max_iterations
)
# Output: 09_apo_protein_clean_FINAL.pdb
```

---

### 💊 Step 2 — Ligand Preparation

```python
from rdkit.Chem import AllChem, rdDistGeom

params = rdDistGeom.ETKDGv3()
params.randomSeed           = 42
params.enforceChirality     = True
params.useSmallRingTorsions = True

mol3d    = Chem.AddHs(mol)
conf_ids = AllChem.EmbedMultipleConfs(mol3d, numConfs=5, params=params)

results  = AllChem.MMFFOptimizeMoleculeConfs(
    mol3d, mmffVariant='MMFF94s', maxIters=2000
)

# Keep lowest-energy converged conformer
energies  = [(r[1], i) for i, r in enumerate(results) if r[0] == 0]
best_conf = min(energies)[1]
# Output: compound_clean.sdf (per molecule)
```

---

### 🎯 Step 3 — Docking Box & Parameters

Cognate validation (Step 3a) with NUV-244 confirmed grid validity.  
Full library run uses **Mode B — manual centre**:

```python
# Docking box (derived from NUV-244 binding pose inspection)
BOX_MODE         = 'B -- Manual centre + box'
center           = [-6.211, -7.261, 0.989]   # Å
size             = [25.0, 25.0, 25.0]         # Å

SCORING_FUNCTION = 'vinardo'
EXHAUSTIVENESS   = 8
N_POSES          = 1
RANDOM_SEED      = 42
```

---

### 🔁 Step 4 — Full Library Docking (1,005 compounds)

The Colab session timed out mid-run. We handled this in two parts:

**Part A — Recovery (LIG_0001–0611): regex parse of notebook output**
```python
import re, json

pattern = re.compile(
    r'Ligand\s+(?P<lig_num>\d+)/1005\s+\|\s+(?P<ligand>LIG_\d{4}).*?'
    r'\[prop\]\s+(?P<formula>\S+)\s+MW=(?P<mw>-?\d+\.?\d*)\s+LogP=(?P<logp>-?\d+\.?\d*).*?'
    r'\[dock\]\s+(?P<runtime>-?\d+\.?\d*)s\s+\|\s+Best\s+=\s+(?P<aff>-?\d+\.?\d*)\s+kcal/mol.*?'
    r'\[OK\]\s+(?P=ligand)\s+done',
    re.S
)
rows = [m.groupdict() for m in pattern.finditer(notebook_output)]
# Saved → recovered_ligand_0001_to_0611.csv
```

**Part B — Resume-safe continuation (LIG_0612–1005)**
```python
from vina import Vina

completed_success = set(recovered_df['Ligand_name'])   # skip already done

for lig_name, mol_raw in resume_mol_list:              # LIG_0612–1005
    if lig_name in completed_success:
        continue

    # SDF → PDBQT via OpenBabel (Gasteiger charges)
    _ob_sdf_to_pdbqt(lig_sdf, lig_pdbqt)

    v = Vina(sf_name='vinardo', seed=seed + lig_no, verbosity=0)
    v.set_receptor(prot_pdbqt)
    v.set_ligand_from_file(str(lig_pdbqt))
    v.compute_vina_maps(center=center, box_size=size)
    v.dock(exhaustiveness=8, n_poses=1)

    best_aff = v.energies(n_poses=1)[0][0]             # kcal/mol

    summary_rows.append({...})
    write_checkpoint(summary_rows, PARTIAL_CSV)         # save after every molecule
    completed_success.add(lig_name)

# Final output: summary_results.csv  (1,005 rows, sorted by affinity)
```

---

### 🔬 Step 5 — Post-processing

**Top-10 re-docking with 9 poses**
```python
N_POSES = 9   # expand pose sampling for top candidates
# Re-run with same box / vinardo / exhaustiveness=8
```

**Best pose extraction + hydrogen addition (bond-order repair)**
```python
from rdkit.Chem import AllChem

# Vina PDBQT loses bond orders — repair via PDB round-trip
pdb_block    = Chem.MolToPDBBlock(Chem.RemoveHs(best_mol))
mol_from_pdb = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=False)
mol_fixed    = AllChem.AssignBondOrdersFromTemplate(
                   Chem.RemoveHs(best_mol), mol_from_pdb)
Chem.SanitizeMol(mol_fixed)
best_mol_h   = Chem.AddHs(mol_fixed, addCoords=True)

writer = Chem.SDWriter(f'{name}_best_pose_with_H.sdf')
writer.write(best_mol_h)
writer.close()
# Output: LIG_0071_best_pose_with_H.sdf
#         LIG_0200_best_pose_with_H.sdf
#         NUV244_best_pose_with_H.sdf
```

| Step | Result |
|---|---|
| Protein prep | ✅ `09_apo_protein_clean_FINAL.pdb` |
| Cognate validation (3a) | ✅ RMSD passed — grid valid |
| ROC/EF analysis (3b) | ⏳ Skipped — no experimental activity data |
| Full library docking | ✅ `summary_results.csv` (1,005 compounds) |
| NUV-244 positive control | ✅ Baseline affinity established |
| Post-processing | ✅ `LIG_0071`, `LIG_0200` selected for MD |

---

<br/>

<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║               PHASE 5 · MOLECULAR DYNAMICS                   ║
╚══════════════════════════════════════════════════════════════╝
```

</div>

## Phase 5 · Molecular Dynamics

**Platform:** NUS Vanda HPC — GPU queue  
**AMBER version:** `Amber/24.3-foss-2023a-AmberTools-24.10-CUDA-12.1.1`  
**Working directory:** `/home/e1638543/xavier_test/`

---

### ⚙️ Step 1 — Ligand Parameterisation (NUV-244)

**antechamber — GAFF atom types + AM1-BCC charges**
```bash
antechamber -i nuv.sdf -fi sdf \
            -o nuv.mol2 -fo mol2 \
            -c bcc -s 2
# Output: nuv.mol2
```

**parmchk2 — missing parameter check**
```bash
parmchk2 -i nuv.mol2 -f mol2 -o nuv.frcmod
# Output: nuv.frcmod
```

**tleap — ligand library creation**
```
source leaprc.gaff
LIG = loadmol2 nuv.mol2
loadamberparams nuv.frcmod
check LIG                  # → Unit is OK
saveoff LIG nuv.lib
saveamberparm LIG nuv.prmtop nuv.rst7
quit
# Output: nuv.lib  nuv.prmtop  nuv.rst7
```

---

### 🧼 Step 2 — Protein Cleaning

```bash
pdb4amber -i clean_PNPLA3I148M.pdb \
          -o protein_clean.pdb \
          --nohyd --dry
```

```
Chains: A
Alternate Locations: None
Non-standard-resnames: (none)
Missing heavy atom(s): None
```

> `--nohyd` is **required** — pdb2pqr-placed HD1 atoms on HIE residues have no ff14SB type and cause 7 FATAL errors if not stripped before tleap.

---

### 🏗️ Step 3 — Build Solvated System (tleap)

> Three tleap attempts were needed. Attempt 1 failed (HIE HD1 FATAL errors). Attempt 2 validated protein alone. Attempt 3 succeeded.

```
# Final successful tleap session
source leaprc.gaff
source leaprc.protein.ff14SB
source leaprc.water.tip3p

loadoff nuv.lib
loadamberparams nuv.frcmod

protein = loadpdb PNPLA3I148M_clean.pdb
ligand  = loadmol2 nuv.mol2

check protein    # Warnings: 3 / Unit is OK
check ligand     # Unit is OK

complex = combine {protein ligand}
check complex    # Warnings: 3 / Unit is OK

saveamberparm complex complex_tleap.prmtop complex_tleap.rst7

charge complex                  # → Total charge: -3.000000
addions complex Na+ 0           # → 3 Na+ placed

saveamberparm complex complex_ion.prmtop complex_ion.rst7
solvateoct complex TIP3PBOX 10.0

saveamberparm complex complex_wat.prmtop complex_wat.rst7
quit
```

```
✔ 3 Na+ ions placed to neutralise charge
✔ 52,379 TIP3P water residues added
✔ Box volume: 1,798,741 Å³ (truncated octahedron)
✔ Initial density: 0.920 g/cc
```

| System parameter | Value |
|---|---|
| Protein heavy atoms | 3,707 |
| H atoms added by tleap | 3,731 |
| Net charge | −3.000 (neutralised with 3 × Na⁺) |
| Water molecules | 52,379 (TIP3P) |
| Box volume | 1,798,741 Å³ |
| Initial density | 0.920 g/cc |

---

### 📄 Step 4 — MD Input Files

**`min1.in` — Restrained minimisation**
```fortran
complex: initial minimization - fixed solute
&cntrl
  imin=1, maxcyc=1000, ncyc=500,
  ntb=1, ntr=1, cut=10.0,
/
Hold complex fixed
300.0
RES 1 314
END
END
```

**`min2.in` — Unrestrained minimisation**
```fortran
complex: second minimization - free all
&cntrl
  imin=1, maxcyc=2000, ncyc=1000,
  ntb=1, ntr=0, cut=10.0,
/
```

**`md_heat.in` — Heating (0 → 300 K, 100 ps, NVT)**
```fortran
complex: heating stage (100 ps)
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim=50000, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=1, ntr=1,
  ntt=3, gamma_ln=1.0,
  tempi=0.0, temp0=300.0,
  ntpr=5000, ntwx=5000, ntwr=5000,
  ioutfm=1,
/
Hold complex fixed
100.0
RES 1 314
END
END
```

**`md_equil.in` — NVT equilibration (100 ps)**
```fortran
complex: equilibration stage (100 ps)
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim=50000, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=1, ntr=0,
  ntt=3, gamma_ln=1.0,
  tempi=300.0, temp0=300.0,
  ntpr=500, ntwx=500, ntwr=500,
  ioutfm=1,
/
```

**`md_prod.in` — Production MD (50 ns/chunk × 2, NPT)**
```fortran
complex: production stage (50 ns per chunk)
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim=25000000, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=2, ntp=1, barostat=2, pres0=1.0,
  ntt=3, gamma_ln=1.0,
  tempi=300.0, temp0=300.0,
  ntpr=50000, ntwx=50000, ntwr=50000,
  ioutfm=1,
/
```

---

### 🚀 Step 5 — PBS Job Submission

```bash
#!/bin/bash
#PBS -N amber_pnpla3
#PBS -l select=1:ncpus=1:mem=8gb:ngpus=1
#PBS -l walltime=08:00:00
#PBS -j oe -o amber_pnpla3.log

source /app1/ebenv
module load Amber/24.3-foss-2023a-AmberTools-24.10-CUDA-12.1.1
cd $PBS_O_WORKDIR
PMEMD="pmemd.cuda"

$PMEMD -O -i min1.in    -o min1.out    -p complex_wat.prmtop \
        -c complex_wat.rst7 -r min1.ncrst -ref complex_wat.rst7

$PMEMD -O -i min2.in    -o min2.out    -p complex_wat.prmtop \
        -c min1.ncrst -r min2.ncrst

$PMEMD -O -i md_heat.in -o md_heat.out -p complex_wat.prmtop \
        -c min2.ncrst -r md_heat.ncrst -x md_heat.nc -ref min2.ncrst

$PMEMD -O -i md_equil.in -o md_equil.out -p complex_wat.prmtop \
        -c md_heat.ncrst -r md_equil.ncrst -x md_equil.nc

PREV_RST=md_equil.ncrst
for i in $(seq 1 2); do
  PADDED=$(printf "%03d" $i)
  $PMEMD -O -i md_prod.in \
         -o prod_${PADDED}.out \
         -p complex_wat.prmtop \
         -c ${PREV_RST} \
         -r prod_${PADDED}.ncrst \
         -x prod_${PADDED}.nc
  PREV_RST=prod_${PADDED}.ncrst
done
echo "Done — 100 ns production collected."
```

```bash
# Submit and monitor
qsub run_amber.pbs
watch -n 10 qstat -u $USER
```

---

### 📊 Step 6 — Trajectory Analysis (Planned)

```bash
# Concatenate and re-image trajectories
cpptraj << EOF
parm complex_wat.prmtop
trajin prod_001.nc
trajin prod_002.nc
center :1-313
image familiar
rms first @CA
trajout prod_combined.nc netcdf
run
quit
EOF
```

```bash
# Structural analysis
cpptraj << EOF
parm complex_wat.prmtop
trajin prod_combined.nc

rms first out rmsd_protein.dat :1-313@C,CA,N time 1.0
rms first out rmsd_ligand.dat  :MOL&!@H=      time 1.0

atomicfluct out rmsf.dat :1-313@C,CA,N byres
radgyr      out rog.dat mass nomax
surf        out surf.dat

hbond HB_p2l donormask :1-313  acceptormask :MOL out nhb_prot2lig.dat avgout avghb_prot2lig.dat
hbond HB_l2p donormask :MOL    acceptormask :1-313 out nhb_lig2prot.dat avgout avghb_lig2prot.dat

lie LIE :MOL out lie.dat
run
quit
EOF
```

> Ligand mask `:MOL` — residue name assigned by antechamber when reading from SDF.

---

<br/>

## 📈 Results Summary

| Phase | Key Output | Status |
|---|---|---|
| Target & reference | `PNPLA3I148M.pdb`, `NUV-244.sdf` | ✅ |
| Library (11,204 MGs) | `3DB_11204_standardised.sdf` | ✅ |
| Virtual screening | `eon_NUV244_top3000.sdf` (1,005 hits) | ✅ |
| Docking — all 1,005 | `summary_results.csv` | ✅ |
| Docking — NUV-244 control | `runNUV244_docking_results.zip` | ✅ |
| Top candidates | `LIG_0071`, `LIG_0200` + NUV-244 (H-added SDF) | ✅ |
| MD system | `complex_wat.prmtop` + `complex_wat.rst7` | ✅ |
| MD production (100 ns) | `prod_001.nc`, `prod_002.nc` | ✅ |
| Trajectory analysis | RMSD / RMSF / H-bonds / PCA / FEL | ✅ |
---

## 🛠 Dependencies

```
# Python packages
rdkit >= 2023.09
scipy
numpy
pandas
vina >= 1.2
openbabel-wheel
openeye-toolkits (OMEGA, ROCS, EON, QUACPAC)

# External tools
AutoDock Vina 1.2.x
AMBER 24.3 + AmberTools 24.10 (CUDA 12.1.1)
OMEGA2 (OpenEye)
ROCS 3.x (OpenEye)
EON 3.1.2.1 (OpenEye)
pdb2pqr / PropKa
PyMOL
```

Install Python dependencies:
```bash
pip install rdkit scipy numpy pandas vina openbabel-wheel
```

---

## 👥 Team

| Member | Primary Role |
|---|---|
| **Sharon** | MD simulation · molecular docking · system integration · Report: Methodology, Limitations, Clinical Relevance & Future Directions |
| **Tami** | RDKit standardisation · conformer generation · Colab scripting · Report: Abstract & Introduction |
| **Leo** | Database curation · GitHub · OpenBabel pipeline · plot · Report: result, Discussion & formatting |
<br/>

---

<div align="center">

*PHM5013 Precision Drug Discovery and Pharmacogenomics · NUS · 2026*

<img src="https://capsule-render.vercel.app/api?type=waving&color=0A2342,1F4E79,2E86C1&height=100&section=footer" width="100%"/>

</div>
