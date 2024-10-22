#
#  Makefile for the MsSpec Fortran 90 program
#
#         by S. Tricot and D. Sébilleau
#
#                                Last version: 18 Jun 2021 #
#  Compiler
#
FC=gfortran
#
#  Compile flags
#
#FFLAGS=
#FFLAGS=-g -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal
FFLAGS=-g -ffpe-trap=zero,invalid,overflow,underflow
#FFLAGS=-ffast-math -O3
#
#  Link flags
#
#  LDFLAGS : -L/path_to_libraries
#  LDLIBS  : -lname_of_lib
#

#LDFLAGS = -L/REXS_library/LAPACK_LIBRARY
#LDLIBS  = -llapack -lrefblas

LDFLAGS = -L/usr/lib/x86_64-linux-gnu/lapack/      #-L/usr/lib/x86_64-linux-gnu
LDLIBS  = -llapack -l:libblas.a



#  Executable name
#
EXE=spec
#
#  Building directory
#
BUILDDIR:=build

.PHONY: clean

cmn_DEPS:=REXS_library/ACCURACY_LIBRARY/accuracy.f90 \
		   REXS_library/DIMENSIONS_LIBRARY/dimensions.f90 \
		   REXS_library/UTILITIES_LIBRARY/simple_numbers.f90 \
		   REXS_library/UTILITIES_LIBRARY/mathematical_constants.f90 \
		   REXS_library/UTILITIES_LIBRARY/physical_constants.f90 \
		   REXS_library/UTILITIES_LIBRARY/locate.f90 \
		   REXS_library/UTILITIES_LIBRARY/sort1.f90 \
		   REXS_library/UTILITIES_LIBRARY/sort2.f90 \
		   REXS_library/UTILITIES_LIBRARY/atomic_properties.f90
cmn_OBJS:=$(patsubst %.f90,%.o, $(cmn_DEPS))

store_SRCS:=$(cmn_DEPS) \
		REXS_library/STORAGE_LIBRARY/Auger_rad_elts.f90 \
		REXS_library/STORAGE_LIBRARY/correlation_storage.f90 \
		REXS_library/STORAGE_LIBRARY/current_beam.f90 \
		REXS_library/STORAGE_LIBRARY/eels_rad_elts.f90 \
		REXS_library/STORAGE_LIBRARY/ex_beam.f90 \
		REXS_library/STORAGE_LIBRARY/functions_storage.f90 \
		REXS_library/STORAGE_LIBRARY/in_beam.f90 \
		REXS_library/STORAGE_LIBRARY/input_values.f90 \
		REXS_library/STORAGE_LIBRARY/ms_storage.f90 \
		REXS_library/STORAGE_LIBRARY/o1_beam.f90 \
		REXS_library/STORAGE_LIBRARY/o2_beam.f90 \
		REXS_library/STORAGE_LIBRARY/o3_beam.f90 \
		REXS_library/STORAGE_LIBRARY/path_storage.f90 \
		REXS_library/STORAGE_LIBRARY/ped_rad_elts.f90 \
		REXS_library/STORAGE_LIBRARY/spectro_storage.f90 \
		REXS_library/STORAGE_LIBRARY/symmetry_storage.f90 \
		REXS_library/STORAGE_LIBRARY/t_matrix.f90 \
		REXS_library/STORAGE_LIBRARY/various_storage.f90 \
		REXS_library/STORAGE_LIBRARY/wave_function.f90
tool_SRCS:=$(cmn_DEPS) \
		REXS_library/COMBINATORICS_LIBRARY/combinatorics.f90 \
		REXS_library/ELECTRON_CHOICE_LIBRARY/electron_choice.f90
util_SRCS:=$(cmn_DEPS) \
		REXS_library/ANGULAR_MOMENTUM_LIBRARY/angular_momentum.f90 \
		REXS_library/RENORMALIZATION_LIBRARY/renormalization.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/arc_sin.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/Euler.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/Hankel_polynomials.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/Legendre_functions.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/spherical_Bessel.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/spherical_harmonics.f90 \
		REXS_library/VARIOUS_FUNCTIONS_LIBRARY/Wigner_rotations.f90 \
		REXS_library/VECTOR_LIBRARY/vector.f90
init_SRCS:=$(cmn_DEPS) \
		REXS_library/VIBRATIONS_LIBRARY/vibrations.f90 \
		REXS_library/INPUT_OUTPUT_LIBRARY/initialize_calc.f90 \
		REXS_library/AUGER_LIBRARY/Auger_multiplet.f90
#
#  Read input data file:
#
read_SRCS:=$(cmn_DEPS) \
		REXS_library/INPUT_OUTPUT_LIBRARY/read_data.f90 \
		REXS_library/INPUT_OUTPUT_LIBRARY/read_ext_files.f90
read_OBJS:=$(patsubst %.f90,%.o, $(store_SRCS) $(tool_SRCS) $(util_SRCS) $(init_SRCS) $(read_SRCS))

#
#  Calculation:
#

calc_DEPS:=$(cmn_DEPS) \
		REXS_library/SCATTERING_AMPLITUDE_LIBRARY/scattering_amplitude.f90 \
		REXS_library/CHECK_LIBRARY/check.f90 \
		REXS_library/ANALYZER_LIBRARY/direction_ana.f90


calc_SRCS:=$(cmn_DEPS) \
		REXS_library/ELECTRON_PHOTON_INTERACTION_LIBRARY/electron_photon.f90 \
		REXS_library/SYMMETRIZATION_LIBRARY/symmetrization.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/mean_free_path.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/beam_amplitude.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/calc_tau_ce.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/calc_tau_mi.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/calc_tau_re.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/calc_tau_se.f90 \
		REXS_library/MULTIPLE_SCATTERING_LIBRARY/calc_tau_sm.f90 \
		REXS_library/SPECTROSCOPY_LIBRARY/ped_ms.f90 \
		REXS_library/CROSS_SECTION_LIBRARY/cross_section.f90

treat_SRCS:=$(cmn_DEPS) \
		REXS_library/TREATMENT_LIBRARY/treatment_ped.f90 \
		REXS_library/TREATMENT_LIBRARY/data_treatment.f90

calc_OBJS:=$(patsubst %.f90,%.o, $(calc_DEPS) $(treat_SRCS) $(calc_SRCS))


SRCS:= $(patsubst %.o,%.f90, $(read_OBJS) $(calc_OBJS))
OBJS:= $(addprefix $(BUILDDIR)/,$(notdir $(cmn_OBJS) $(read_OBJS) $(calc_OBJS)))


all: obj $(EXE)


obj: src $(OBJS)


$(EXE): $(OBJS) $(BUILDDIR)/spec.f90
	@echo "building main $@..."
	@$(FC) $(FFLAGS) -J $(BUILDDIR) -o $@ $^ $(LDFLAGS) $(LDLIBS)


%.o: %.f90
	@echo "Compiling $@..."
	@$(FC) $(FFLAGS) -J $(BUILDDIR) -I. -o $@ -c $^


src: $(SRCS) spec.f90
	@echo "updating source tree..."
	@mkdir -p $(BUILDDIR)
	@rsync -av $^ $(BUILDDIR)


clean:
	@echo "Cleaning..."
	@rm -rf $(BUILDDIR)
	@rm -rf $(EXE)
