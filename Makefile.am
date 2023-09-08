## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
## 02110-1301, USA.
##

ACLOCAL_AMFLAGS = -I m4

SUBDIRS =        \
	external_libs liboct_parser \
	src share doc testsuite scripts

EXTRA_DIST =     \
	ChangeLog PACKAGING \
	TODO \
	debian/changelog debian/control debian/copyright debian/dirs     \
	debian/README.Debian debian/rules

DISTCLEANFILES = \
	config.log config.status

DEB_DIR = deb-work

check-unittests:
	$(MAKE) -C external_libs check
	$(MAKE) -C liboct_parser check
	$(MAKE) -C src check
	$(MAKE) -C share check
	$(MAKE) -C doc check
	$(MAKE) -C testsuite check-unittests

check-short:
	$(MAKE) -C external_libs check
	$(MAKE) -C liboct_parser check
	$(MAKE) -C src check
	$(MAKE) -C share check
	$(MAKE) -C doc check
	$(MAKE) -C testsuite check-short

check-long:
	$(MAKE) -C external_libs check
	$(MAKE) -C liboct_parser check
	$(MAKE) -C src check
	$(MAKE) -C share check
	$(MAKE) -C doc check
	$(MAKE) -C testsuite check-long

check-performance:
	$(MAKE) -C share check
	$(MAKE) -C testsuite check-performance

deb: dist
	rm -rf $(DEB_DIR)
	mkdir $(DEB_DIR)
	mv -f $(PACKAGE)-$(VERSION).tar.gz $(DEB_DIR)
	(cd $(DEB_DIR); tar xzf *gz; cd $(PACKAGE)-$(VERSION); debuild -E -ePATH)
	rm -rf $(DEB_DIR)/$(PACKAGE)-$(VERSION)
	scp $(DEB_DIR)/octopus*.deb \
		marques@nautilus.fis.uc.pt:/server/www/tddft.org/programs/octopus/download/
