#//************************************************************************//
#// Copyright (C) 2012-2012 Katsuhiko Nishimra                             // 
#//                                                                        // 
#// This file is part of MolDS.                                            // 
#//                                                                        // 
#// MolDS is free software: you can redistribute it and/or modify          // 
#// it under the terms of the GNU General Public License as published by   // 
#// the Free Software Foundation, either version 3 of the License, or      // 
#// (at your option) any later version.                                    // 
#//                                                                        // 
#// MolDS is distributed in the hope that it will be useful,               // 
#// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
#// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
#// GNU General Public License for more details.                           // 
#//                                                                        // 
#// You should have received a copy of the GNU General Public License      // 
#// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
#//************************************************************************//

.PHONY: all clean

ifeq ($(OPENBLASDIR),)
   ifeq ($(BLASLIB),openblas)
      CONFIGURE_OPTIONS+=--with-blas=openblas
   endif
else
   CONFIGURE_OPTIONS+=--with-blas=openblas
   CFLAGS+= -I$(OPENBLASDIR)/include
   LDFLAGS+= -L$(OPENBLASDIR)/lib
endif

all: MolDS.out

build/Makefile:
	mkdir -p build
	cd build; \
	CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" ../configure $(CONFIGURE_OPTIONS)

build/MolDS: build/Makefile $(wildcard *.cpp) $(wildcard */*.cpp) $(wildcard */*.h)
	+make -C build

MolDS.out: build/MolDS
	cp $< $@

clean:
	rm -rf build
