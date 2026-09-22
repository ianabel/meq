"""Write the freegs4e reference as an MFEM mesh + GridFunction, to seed MEQ.

MEQ's [initialguess] Type = "gridfunction" reads an mfem::Mesh and an
mfem::GridFunction on it, and INTERPOLATES through meq::FieldTransfer when the
mesh is not the one being solved on.  So the guess does not have to live on
MEQ's mesh, or share its degree, or even its element type -- which means it can
be a plain Cartesian quad mesh carrying the reference at its vertices.

Both file formats are ASCII and small: a v1.0 mesh, and a GridFunction whose
header names its collection.  H1_2D_P1 puts one dof at each vertex, in vertex
order, which is what makes this a dozen lines rather than a dof-mapping job.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def write_guess(mesh_path, gf_path, R, Z, psi, box, n=96):
    """Sample `psi(Z,R)` onto an n x n quad mesh over `box` and write both files."""
    rs = np.linspace(box["rmin"], box["rmax"], n + 1)
    zs = np.linspace(box["zmin"], box["zmax"], n + 1)

    interp = RegularGridInterpolator((Z, R), psi, bounds_error=False,
                                     fill_value=None)   # extrapolate at the rim
    ZZ, RR = np.meshgrid(zs, rs, indexing="ij")
    vals = interp(np.stack([ZZ.ravel(), RR.ravel()], axis=-1))

    nv = (n + 1) * (n + 1)
    vid = lambda i, j: i * (n + 1) + j          # i over z, j over R

    with open(mesh_path, "w") as f:
        f.write("MFEM mesh v1.0\n\ndimension\n2\n\n")
        f.write(f"elements\n{n * n}\n")
        for i in range(n):
            for j in range(n):
                # attribute 1, geometry 3 = SQUARE, counterclockwise
                f.write(f"1 3 {vid(i,j)} {vid(i,j+1)} {vid(i+1,j+1)} {vid(i+1,j)}\n")
        f.write(f"\nboundary\n{4 * n}\n")
        for j in range(n):                       # geometry 1 = SEGMENT
            f.write(f"1 1 {vid(0,j)} {vid(0,j+1)}\n")
            f.write(f"1 1 {vid(n,j+1)} {vid(n,j)}\n")
        for i in range(n):
            f.write(f"1 1 {vid(i+1,0)} {vid(i,0)}\n")
            f.write(f"1 1 {vid(i,n)} {vid(i+1,n)}\n")
        f.write(f"\nvertices\n{nv}\n2\n")
        for i in range(n + 1):
            for j in range(n + 1):
                f.write(f"{rs[j]:.16e} {zs[i]:.16e}\n")

    with open(gf_path, "w") as f:
        f.write("FiniteElementSpace\n")
        f.write("FiniteElementCollection: H1_2D_P1\n")
        f.write("VDim: 1\nOrdering: 0\n\n")
        for v in vals:
            f.write(f"{v:.16e}\n")
    return float(np.nanmin(vals)), float(np.nanmax(vals))


if __name__ == "__main__":
    import json, os, sys
    npz, outdir = sys.argv[1], sys.argv[2]
    stem = os.path.splitext(os.path.basename(npz))[0]
    d = np.load(npz, allow_pickle=True)
    meta = json.load(open(os.path.join(outdir, f"{stem}-meta.json")))

    psi = np.asarray(d["psi"], float).T          # ( R, Z ) -> ( Z, R )
    psi = psi - meta["psi_surface"]              # MEQ's gauge
    lo, hi = write_guess(os.path.join(outdir, f"{stem}-guess.mesh"),
                         os.path.join(outdir, f"{stem}-guess.gf"),
                         np.asarray(d["R"], float), np.asarray(d["Z"], float),
                         psi, meta["box"])
    print(f"  {stem}: guess psi in [{lo:.4e}, {hi:.4e}], "
          f"reference axis {float(d['psi_axis']) - meta['psi_surface']:.4e}")
