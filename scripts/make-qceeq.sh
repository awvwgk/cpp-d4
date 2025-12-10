#!/bin/bash

INCLUDE="../include"
SRC="../src"

dosed(){
	sed -n '/namespace multicharge {/,/} \/\/ namespace multicharge/p' $1 | \
	 sed '/namespace multicharge {/d; /} \/\/ namespace multicharge/d' | \
	 awk '/./,EOF' > $2
}

##################################################################

dosed "${INCLUDE}/dftd_eeq.h" eeq.txt

cat > qceeq.h << 'EOT'
/*
 * Electronegativity equilibration (EEQ) model for DFT-D4.
 * Contains only the essential parts for DFT-D4.
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */

#pragma once

#include "qcinpdat.h" // TGeomInput class
#include "qcmat1.h"   // TRVector and TRMatrix class
#include "qcncoord.h"

namespace multicharge {

EOT

cat eeq.txt >> qceeq.h
rm eeq.txt

cat >> qceeq.h << 'EOT'
} // namespace multicharge
EOT

######

dosed "${SRC}/dftd_eeq.cpp" eeq.txt

cat > qceeq.cpp << 'EOT'
/*
 * Electronegativity equilibration (EEQ) model for DFT-D4.
 * This implementation contains only the essential parts for DFT-D4.
 *
 * Sebastian Ehlert <ehlert@thch.uni-bonn.de> (SAW190521)
 * Marvin Friede <friede@thch.uni-bonn.de> (MF161222)
 *
 * Extension for Fr, Ra and Actinides by:
 * Marvin Friede <friede@thch.uni-bonn.de> (MF121223)
 */

#include <cmath>

#include "qcinpdat.h" // TGeomInput class
#include "qclineq.h"  // Linear Algebra routines
#include "qcmat1.h"   // TRVector and TRMatrix class
#include "qcmat2.h"   // BLAS routines

// CN-related logic
#include "qcncoord.h"

// always include self
#include "qceeq.h"

// wrap all charge models in the multicharge namespace
namespace multicharge {
EOT

cat eeq.txt >> qceeq.cpp
rm eeq.txt

cat >> qceeq.cpp << 'EOT'
} // namespace multicharge
EOT

##########################################################################

sed -i "s/TMolecule/TGeomInput/" *.h
sed -i "s/TMolecule/TGeomInput/" *.cpp
