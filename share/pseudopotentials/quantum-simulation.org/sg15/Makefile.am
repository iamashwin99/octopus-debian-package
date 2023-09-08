## Process this file with automake to produce Makefile.in

## Copyright (C) 2015 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 021.2-1.21, USA.
##

sg15dir = $(pkgdatadir)/pseudopotentials/quantum-simulation.org/sg15

dist_sg15_DATA =                      \
	Ag_ONCV_PBE-1.2.upf           \
	Al_ONCV_PBE-1.2.upf           \
	Ar_ONCV_PBE-1.2.upf           \
	As_ONCV_PBE-1.2.upf           \
	Au_ONCV_PBE-1.2.upf           \
	Ba_ONCV_PBE-1.2.upf           \
	Be_ONCV_PBE-1.2.upf           \
	Bi_ONCV_PBE-1.2.upf           \
	B_ONCV_PBE-1.2.upf            \
	Br_ONCV_PBE-1.2.upf           \
	Ca_ONCV_PBE-1.2.upf           \
	Cd_ONCV_PBE-1.2.upf           \
	Cl_ONCV_PBE-1.2.upf           \
	C_ONCV_PBE-1.2.upf            \
	Co_ONCV_PBE-1.2.upf           \
	Cr_ONCV_PBE-1.2.upf           \
	Cs_ONCV_PBE-1.2.upf           \
	Cu_ONCV_PBE-1.2.upf           \
	Fe_ONCV_PBE-1.2.upf           \
	F_ONCV_PBE-1.2.upf            \
	Ga_ONCV_PBE-1.2.upf           \
	Ge_ONCV_PBE-1.2.upf           \
	He_ONCV_PBE-1.2.upf           \
	Hf_ONCV_PBE-1.2.upf           \
	Hg_ONCV_PBE-1.2.upf           \
	H_ONCV_PBE-1.2.upf            \
	In_ONCV_PBE-1.2.upf           \
	I_ONCV_PBE-1.2.upf            \
	Ir_ONCV_PBE-1.2.upf           \
	K_ONCV_PBE-1.2.upf            \
	Kr_ONCV_PBE-1.2.upf           \
	La_ONCV_PBE-1.2.upf           \
	Li_ONCV_PBE-1.2.upf           \
	Mg_ONCV_PBE-1.2.upf           \
	Mn_ONCV_PBE-1.2.upf           \
	Mo_ONCV_PBE-1.2.upf           \
	Na_ONCV_PBE-1.2.upf           \
	Nb_ONCV_PBE-1.2.upf           \
	Ne_ONCV_PBE-1.2.upf           \
	Ni_ONCV_PBE-1.2.upf           \
	N_ONCV_PBE-1.2.upf            \
	O_ONCV_PBE-1.2.upf            \
	Os_ONCV_PBE-1.2.upf           \
	Pb_ONCV_PBE-1.2.upf           \
	Pd_ONCV_PBE-1.2.upf           \
	P_ONCV_PBE-1.2.upf            \
	Pt_ONCV_PBE-1.2.upf           \
	Rb_ONCV_PBE-1.2.upf           \
	Re_ONCV_PBE-1.2.upf           \
	Rh_ONCV_PBE-1.2.upf           \
	Ru_ONCV_PBE-1.2.upf           \
	Sb_ONCV_PBE-1.2.upf           \
	Sc_ONCV_PBE-1.2.upf           \
	Se_ONCV_PBE-1.2.upf           \
	Si_ONCV_PBE-1.2.upf           \
	Sn_ONCV_PBE-1.2.upf           \
	S_ONCV_PBE-1.2.upf            \
	Sr_ONCV_PBE-1.2.upf           \
	Ta_ONCV_PBE-1.2.upf           \
	Tc_ONCV_PBE-1.2.upf           \
	Te_ONCV_PBE-1.2.upf           \
	Ti_ONCV_PBE-1.2.upf           \
	Tl_ONCV_PBE-1.2.upf           \
	V_ONCV_PBE-1.2.upf            \
	W_ONCV_PBE-1.2.upf            \
	Xe_ONCV_PBE-1.2.upf           \
	Y_ONCV_PBE-1.2.upf            \
	Zn_ONCV_PBE-1.2.upf           \
	Zr_ONCV_PBE-1.2.upf           \
	COPYING                       \
	README.Octopus

# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_sg15_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_sg15_DATA); \
	fi

.PHONY: check
