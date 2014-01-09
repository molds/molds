#//************************************************************************//
#// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
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
CONFIGURE_OPTIONS =
CFLAGS =
CXXFLAGS =
LDFLAGS =
override CXXFLAGS += $(CFLAGS)
ifeq ($(OPENBLASDIR),)
   ifeq ($(BLASLIB),openblas)
      override CONFIGURE_OPTIONS += --with-blas=openblas
   endif
else
   override CONFIGURE_OPTIONS += --with-openblas=$(OPENBLASDIR)
endif

ifeq ($(BOOSTDIR),)
else
   override CONFIGURE_OPTIONS += --with-boost=$(BOOSTDIR)
endif

OUTPUT=molds

all: $(OUTPUT)

build/Makefile:
	mkdir -p build
	cd build; \
	../configure $(CONFIGURE_OPTIONS) CXXFLAGS='$(CXXFLAGS)' LDFLAGS='$(LDFLAGS)'

build/$(OUTPUT): build/Makefile $(wildcard *.cpp) $(wildcard */*.cpp) $(wildcard */*.h)
	+make -C build

$(OUTPUT): build/$(OUTPUT)
	install $< $@

clean:
	rm -rf build $(OUTPUT)
