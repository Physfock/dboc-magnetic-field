Potential energy curves (PECs) and Diagonal Born–Oppenheimer Corrections (DBOC) as functions of magnetic field strength and internuclear distance for the H2, LiH, and HeH+ molecules.

This repository provides a Podman container with a modified version of ChronusQ ([https://github.com/xsligroup/chronusq_public](https://github.com/xsligroup/chronusq_public)). The code writes the GIAO-based one- and two-electron matrix elements required for subsequent correlation calculations to HDF5 format. This data is then reused in PySCF to calculate the DBOC within the approximation discussed in the manuscript: [https://arxiv.org/pdf/2609.15327](https://arxiv.org/pdf/2609.15327) (see also [https://doi.org/10.1063/5.0269984](https://doi.org/10.1063/5.0269984)). The patch of ChronusQ only adds the ability to export the GIAO-based matrix elements to HDF5 format, leaving the original codebase unchanged.

The instructions for standalone patching ChronusQ and PySCF are also available at: [https://github.com/Physfock/cc_tesla](https://github.com/Physfock/cc_tesla). Sample output data obtained using this scheme is provided in DBOC_and_PEC.zip.

The patched and built binaries are already pre-installed in this container, allowing you to run test calculations immediately without building the packages from source.

## Instructions

### 1. Install Podman on Linux

Refer to the official Podman installation guide: [https://podman.io/docs/installation](https://podman.io/docs/installation)

For Debian-based systems (e.g., Ubuntu): The podman package is available in Debian repositories:

```bash
sudo apt-get -y install podman
```

For Fedora-based systems:

```bash
sudo dnf -y install podman
```

### 2. Download the container image

Download the container image into your home directory (`/home/<your_username>/`) using the following link: [https://app.filen.io/#/d/1000f5d4-52f0-433d-8fb6-3e28ec3a8dac%23uJ6zYPeFPCSkMtKWOBNBJAbJ4YTglNKO](https://app.filen.io/#/d/1000f5d4-52f0-433d-8fb6-3e28ec3a8dac%23uJ6zYPeFPCSkMtKWOBNBJAbJ4YTglNKO)

### 3. Navigate to the download directory

Open a terminal and navigate to your home directory (`/home/<your_username>/`), where you downloaded the `.tar` file in the previous step.

### 4. Load the image into Podman

```bash
podman load -i /home/<your_username>/cc_tesla.tar
```

### 5. Verify the image was loaded

```bash
podman images
```

### 6. Create a working directory

Create a directory named `Tesla` in your home directory (`/home/<your_username>/Tesla/`), where your programs will be executed:

```bash
mkdir Tesla
```

### 7. Run the container

```bash
podman run -it -v /home/<your_username>/Tesla/:/mnt/:Z cc_tesla bash
```

This command starts the container and mounts your local `/home/<your_username>/Tesla/` directory to `/mnt/` inside the container. Any changes made to `/mnt/` inside the container will immediately be reflected in your local `Tesla` directory, and vice versa.

You are now inside the container's terminal. All necessary tools and software (e.g., `nano`, `mc`, Python, modified PySCF, and ChronusQ) are pre-installed.

#### Troubleshooting: SELinux permissions

If you encounter permission issues or cannot see files in `/home/<your_username>/Tesla/` from within the container, your system may be running SELinux in enforcing mode. Exit the container and temporarily disable SELinux with:

```bash
sudo setenforce 0
```

#### Run test


```bash
cd mnt
```

Put basis set unc-cc-pvtz.gbs and dboc_uhf_cisd.py to mnt dir and run

```bash
python dboc_uhf_cisd.py
```

### Description of example code

#### DBOC in magnetic fields: UHF vs CISD

Computes the diagonal Born–Oppenheimer correction (DBOC) and potential energy
curves of H₂ and LiH in a magnetic field parallel to the molecular axis, at the
UHF and CISD levels.

Integrals (GIAO basis) are generated with ChronusQ; UHF, CISD and wavefunction
overlaps are computed with PySCF. The DBOC is obtained by finite differences of
phase-corrected overlaps between the reference and displaced geometries.

#### Requirements

ChronusQ, Python 3, PySCF, NumPy, SciPy, h5py, Matplotlib.

#### Usage

Set `MOLECULE`, `B_FIELD` and `BASIS_FILE` at the top of `dboc_uhf_cisd.py`, then:

```bash
python dboc_uhf_cisd.py
```

#### Output

- `DBOC_UHF.txt`, `DBOC_CISD.txt` — R (bohr), DBOC (hartree)
- `PEC_UHF.txt`, `PEC_CISD.txt` — R (bohr), energy (hartree)
- `DBOC_UHF_vs_CISD_<MOLECULE>_PRA.png` — plot

</content>
