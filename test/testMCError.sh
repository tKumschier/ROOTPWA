#!/bin/bash

### BEGIN CONFIGURATION ###

	### BEGIN PHASE SPACE ###

		# mass bin [MeV/c^2]
		MASS=$1  #1800
		BINWIDTH=80

		# number of phase space events to generate
		NMB_PS_EVENTS=50000
		TESTDIR=/nfs/freenas/tuph/e18/project/compass/analysis/tkumschier/ROOTPWA/test/TEST_FOLDER_MASS/${MASS}/

		SEED_PS=${1}123456

	#-- END PHASE SPACE --#

	### BEGIN PARTICLE DATA ###

		# particle data table
		PARTICLE_DATA_TABLE="${ROOTPWA}/particleData/particleDataTable.txt"

	#-- END PARTICLE DATA --#

	### BEGIN PSEUDO DATA ###

		# number of weighted pseudo data events to generate
		NMB_PSEUDO_EVENTS=500

		SEED_PSEUDO=654321

	#-- END PSEUDO DATA --#

	### BEGIN DEWEIGHT ###

		SEED_DEWEIGHT=789012

	#-- END DEWEIGHT --#

	### BEGIN FIT ###

		SEED_FIT=210987

	#-- END FIT --#

#-- END CONFIGURATION --#


### SOME CONVENIENCE FUNCTIONS

THIS_SCRIPT=$(basename ${0})

function printInfo {
		echo ">>> ${THIS_SCRIPT}: info: ${1}"
}

function printSucc {
		echo "*** ${THIS_SCRIPT}: success: ${1}"
}

function printErr {
		echo "!!! ${THIS_SCRIPT}: ${BASH_LINENO[0]}: error: ${1}"
		exit 1
}

function testStep {
		printInfo "Test ${1} ..."
		echo "    executing: ${2}"
		if eval ${2}; then
				printSucc "${1} was successful."
				echo
		else
				printErr "${1} was not successful. Aborting."
		fi
}


### BEGIN PREPARATION ###

if [ -d "${ROOTPWA}" ]; then
	printInfo "Installation of ROOTPWA found in: ${ROOTPWA}"
	echo "    Starting MC test ..."
else
	printErr "No installation of ROOTPWA found. Aborting."
fi

if [ -d "${TESTDIR}" ]; then
	printErr "Test directory ${TESTDIR} already exists. Aborting."
fi

printInfo "Creating directories ..."

mkdir -v "${TESTDIR}"
mkdir -v "${TESTDIR}/data"
mkdir -v "${TESTDIR}/weighted_mc_data"
mkdir -v "${TESTDIR}/keyfiles"
mkdir -v "${TESTDIR}/reference_fit"
mkdir -v "${TESTDIR}/fits"
mkdir -v "${TESTDIR}/log"

echo "    Test directory and subfolders created."
printInfo "Changing to test directory ${TESTDIR}."

cd "${TESTDIR}"

exec > >(tee -i ./log/logfile.txt)
exec 2>&1

printInfo "Copying necessary files to test directory ..."

cp -v "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/generator_noBeamSimulation.conf" "./"
cp -v "${ROOTPWA}/rootpwa.config" "./"
cp -v "${ROOTPWA}/test/mcTest/reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root.example" "./reference_fit/bin65_c2pap_bestfits_converged_MASS_1800_1820_N45340.root"

# make sure the integral binning in the config file is set correctly
sed -i.bak 's/^integralBinning.*$/integralBinning                        = [ { "mass": (${MASS}/1000, (${MASS}+20)/1000) } ]/' rootpwa.config
rm -f rootpwa.config.bak

#-- END PREPARATION --#

### BEGIN MONTE CARLO GENERATION ###

# generate phase space data
testStep "generation of phase-space data" \
"${ROOTPWA}/build/bin/generatePhaseSpace \
-s ${SEED_PS} \
-n ${NMB_PS_EVENTS} \
-p \"${PARTICLE_DATA_TABLE}\" \
-M ${MASS} \
-B ${BINWIDTH} \
\"./generator_noBeamSimulation.conf\" \
-o \"./data/phase_space_MASS_${MASS}-$((MASS+BINWIDTH))_N_${NMB_PS_EVENTS}.root\""

# generate the keyfiles
export DESTINATION_DIR="${TESTDIR}/keyfiles"
export PARTICLE_DATA_TABLE
export PARTICLE_DECAY_TABLE="${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/ParticleDecays.key"
export TEMPLATE_KEY_FILES="${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/template.key"
export WAVESET_FILES="${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/wavelist.compass.2008.88waves"
testStep "generation of keyfiles" "${ROOTPWA}/userAnalysisWorkspace/3pi.--+/keyfiles/GenerateKeyfiles.sh"

# create first file manager
testStep "creation of first file manager" "${ROOTPWA}/build/bin/createFileManager"

# calculate amplitudes for generated phase-space data
testStep "calculation of amplitudes for phase-space data" "${ROOTPWA}/build/bin/calcAmplitudes -e generated"

# calculate integrals for generated phase-space data
testStep "calculation of integrals for phase-space data" "${ROOTPWA}/build/bin/calcIntegrals -e generated"

exit 0
