# Static envelope code (copied from CZ)

Unchanged copies of `CZ/srcFortran/{kinds,const,eos_table,opacity_table,
thermo,cz_eos,mixing_length,structure,envelope}_mod.f90` (2026-10-20), the
static solar convection-zone model (hydrostatic + mixing length + radiative
diffusion, integrated from tau = 0.01 at R_sun). CZ_3D uses them only through
`tools/static_bg.f90` to build the background and initial state
(`data/static_init.dat`) on each machine with eos_hhe, so the repository holds
no EOS-generated data. With `--eos=hhe` the FreeEOS reader
(`eos_table_mod`) is compiled but never reads a file.
