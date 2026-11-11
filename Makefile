# CZ_3D: CZTR step 1 in 3D on ARBOR (Yin-Yang shell).
#
#   ./configure     once per machine (or after rebuilding ARBOR or eos_hhe elsewhere)
#   make            build cztr3d.x
#   make data       generate the inputs that are not kept in the repository:
#                   data/eos_table.bin (eos_hhe's eos_table_gen, data/eos_table.nml)
#                   and data/static_init.dat (the static envelope, tools/static_bg.f90)
#   make run        mpirun -np $(NP) ./cztr3d.x runs/static.nml
#   make smoke      the smoke test (runs/smoke.nml, 8 ranks)
#   make clean
#
# Run ./configure first (./configure --help): it finds ARBOR's install and
# eos_hhe's OpenMP build, checks that both use the same compiler, and writes
# config.mk. The MPI wrapper and the optimization flags are ARBOR's (ARBOR_FC,
# ARBOR_FFLAGS in its arbor.mk); this code adds only its own module directory,
# OpenMP and line-length flags. Intel builds go to build-intel/, cztr3d-intel.x.

CONFIG = config.mk
-include $(CONFIG)
ifeq ($(wildcard $(CONFIG)),)
  ifneq ($(MAKECMDGOALS),clean)
    $(error no $(CONFIG): run ./configure first (./configure --help))
  endif
else
  include $(ARBOR_MK)
endif

NP     = 8
FFLAGS = $(ARBOR_FFLAGS) $(LONGLINES) $(ARBOR_INC) -I$(HHE)/mod $(call moddir,$(BUILD)) -I$(BUILD)
LIBS   = $(ARBOR_LIB) $(HHE)/libeos_hhe.a $(OMP)

MODS = kinds_mod const_mod opacity_table_mod mlt_mod grey_rt_mod vmath_mod cz3d_params_mod cz3d_bg_mod \
       cz3d_transport_mod radiative_loss_mod cztr3d_atm_mod pfss_mod cz3d_app_mod

# array exp, expm1, log: Apple's Accelerate on macOS, plain loops elsewhere
ifeq ($(shell uname),Darwin)
  VMATH = src/vmath_accelerate.f90
  LIBS += -framework Accelerate
else
  VMATH = src/vmath_plain.f90
endif
OBJS = $(MODS:%=$(BUILD)/%.o)

.PHONY: all run clean data smoke
all: $(EXE)

$(BUILD):
	mkdir -p $@

# modules in dependency order (each object depends on the ones before it)
$(OBJS): $(CONFIG)
$(BUILD)/kinds_mod.o: src/kinds_mod.f90 | $(BUILD)
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/const_mod.o: src/const_mod.f90 $(BUILD)/kinds_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/opacity_table_mod.o: src/opacity_table_mod.f90 $(BUILD)/kinds_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/mlt_mod.o: src/mlt_mod.f90 $(BUILD)/const_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/grey_rt_mod.o: src/grey_rt_mod.f90 $(BUILD)/const_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/vmath_mod.o: $(VMATH) $(BUILD)/kinds_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/cz3d_params_mod.o: src/cz3d_params_mod.f90 $(ARBOR_LIB) | $(BUILD)
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/cz3d_bg_mod.o: src/cz3d_bg_mod.f90 $(BUILD)/kinds_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/cz3d_transport_mod.o: src/cz3d_transport_mod.f90 $(BUILD)/cz3d_bg_mod.o $(BUILD)/opacity_table_mod.o \
                               $(BUILD)/mlt_mod.o $(BUILD)/grey_rt_mod.o $(BUILD)/vmath_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/radiative_loss_mod.o: src/radiative_loss_mod.f90 $(BUILD)/kinds_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/cztr3d_atm_mod.o: src/cztr3d_atm_mod.f90 $(BUILD)/cz3d_params_mod.o $(BUILD)/cz3d_bg_mod.o $(BUILD)/radiative_loss_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/pfss_mod.o: src/pfss_mod.f90 $(BUILD)/kinds_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@
$(BUILD)/cz3d_app_mod.o: src/cz3d_app_mod.f90 $(BUILD)/cz3d_params_mod.o $(BUILD)/cz3d_transport_mod.o $(BUILD)/cztr3d_atm_mod.o $(BUILD)/pfss_mod.o
	$(ARBOR_FC) $(FFLAGS) -c $< -o $@

$(EXE): src/cz3d_main.f90 $(OBJS)
	$(ARBOR_FC) $(FFLAGS) $< $(OBJS) $(LIBS) -o $@

run: $(EXE)
	mpirun -np $(NP) ./$(EXE) runs/static.nml

# smoke test (runs/smoke.nml, 8 ranks, about a minute): where mpirun runs on
# this machine; on a cluster submit runs/smoke.nml as a short batch job
smoke: $(EXE)
	mpirun -np 8 ./$(EXE) runs/smoke.nml | grep -E "^step +[123]00|done:|checksum"
	@echo "expected (gfortran 14, macOS): state checksum: 38440707489373"

# ---- generated inputs ----
ENV     = tools/envelope
ENVOBJ  = $(BUILD)/env
ENVMODS = kinds_mod const_mod eos_table_mod opacity_table_mod thermo_mod cz_eos_mod mixing_length_mod \
          structure_mod envelope_mod
ENVOBJS = $(ENVMODS:%=$(ENVOBJ)/%.o)

# data/static_init.dat: on data/grid_full.dat (CZTR's grid, kept in the
# repository). Other radial grids: a parameter file data/grid_<name>.nml
# (namelist &grid, see tools/grid/grid_mod.f90) gives data/grid_<name>.dat and
# its init file data/static_init_<name>.dat with
#   make data/static_init_<name>.dat
data: data/eos_table.bin data/static_init.dat

data/eos_table.bin: data/eos_table.nml
	@test -x $(HHE)/eos_table_gen || { echo "need $(HHE)/eos_table_gen: run 'make omp' in eos_hhe"; exit 1; }
	$(HHE)/eos_table_gen data/eos_table.nml

data/static_profile.dat: $(BUILD)/static_bg.x data/opacity_merged_x70_z02.dat
	$(BUILD)/static_bg.x data/static_profile.dat 2.5e14 --eos=hhe

data/static_init.dat: data/static_profile.dat tools/build_static_init.py data/grid_full.dat
	$(PYTHON) tools/build_static_init.py data/static_profile.dat data/grid_full.dat data/static_init.dat

# (kept: the run namelists name both files)
.PRECIOUS: data/grid_%.dat
data/grid_%.dat: data/grid_%.nml data/static_profile.dat $(BUILD)/grid_gen.x
	$(BUILD)/grid_gen.x data/static_profile.dat $< $@

data/static_init_%.dat: data/grid_%.dat data/static_profile.dat tools/build_static_init.py
	$(PYTHON) tools/build_static_init.py data/static_profile.dat $< $@

# the radial grid generator (tools/grid/grid_mod.f90: CZTR's density map)
$(BUILD)/grid/grid_mod.o: tools/grid/grid_mod.f90 $(ENVOBJ)/kinds_mod.o
	mkdir -p $(BUILD)/grid
	$(ARBOR_FC) -O2 $(LONGLINES) $(call moddir,$(BUILD)/grid) -I$(BUILD)/grid -I$(ENVOBJ) -c $< -o $@
$(BUILD)/grid_gen.x: tools/grid_gen.f90 $(BUILD)/grid/grid_mod.o
	$(ARBOR_FC) -O2 $(LONGLINES) -I$(BUILD)/grid -I$(ENVOBJ) $< $(BUILD)/grid/grid_mod.o $(ENVOBJ)/kinds_mod.o -o $@

$(ENVOBJ):
	mkdir -p $@

# envelope modules in dependency order (their own module directory: some names repeat src/'s)
$(ENVOBJ)/%.o: $(ENV)/%.f90 | $(ENVOBJ)
	$(ARBOR_FC) -O2 $(OMP) $(LONGLINES) $(call moddir,$(ENVOBJ)) -I$(ENVOBJ) -I$(HHE)/mod -c $< -o $@
$(ENVOBJ)/const_mod.o: $(ENVOBJ)/kinds_mod.o
$(ENVOBJ)/eos_table_mod.o: $(ENVOBJ)/kinds_mod.o
$(ENVOBJ)/opacity_table_mod.o: $(ENVOBJ)/kinds_mod.o
$(ENVOBJ)/thermo_mod.o: $(ENVOBJ)/eos_table_mod.o
$(ENVOBJ)/cz_eos_mod.o: $(ENVOBJ)/const_mod.o $(ENVOBJ)/eos_table_mod.o $(ENVOBJ)/thermo_mod.o
$(ENVOBJ)/mixing_length_mod.o: $(ENVOBJ)/const_mod.o
$(ENVOBJ)/structure_mod.o: $(ENVOBJ)/cz_eos_mod.o $(ENVOBJ)/opacity_table_mod.o $(ENVOBJ)/mixing_length_mod.o
$(ENVOBJ)/envelope_mod.o: $(ENVOBJ)/structure_mod.o

$(BUILD)/static_bg.x: tools/static_bg.f90 $(ENVOBJS)
	$(ARBOR_FC) -O2 $(OMP) $(LONGLINES) -I$(ENVOBJ) -I$(HHE)/mod $< $(ENVOBJS) $(HHE)/libeos_hhe.a -o $@

clean:
	rm -rf $(BUILD) $(EXE)
