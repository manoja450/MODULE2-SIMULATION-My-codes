TIMESTAMP=$(date +%Y%m%d_%H%M%S)

OUTPUT_DIR="/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o/OUTPUT/run_${TIMESTAMP}"
SLURM_DIR="/home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o/SLURMOUT"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${SLURM_DIR}

cd /home/manoja450/G4WithoutLeadSheilding/MODULE2/CUSTOMOPTICALMODULE2/G4d2o

export LD_LIBRARY_PATH=""
export ROOTSYS="/usr"
export LD_LIBRARY_PATH="/usr/lib64:${LD_LIBRARY_PATH}"

export G4INSTALL="/home/manoja450/geant4-install"
source ${G4INSTALL}/bin/geant4.sh
export LD_LIBRARY_PATH="${G4INSTALL}/lib64:${LD_LIBRARY_PATH}"

export CRYHOME="/home/manoja450/cry_v1.7"
export LD_LIBRARY_PATH="${CRYHOME}/lib:${LD_LIBRARY_PATH}"

export G4NEUTRONHPDATA="/home/manoja450/geant4-data/G4NDL4.7"
export G4LEDATA="/home/manoja450/geant4-data/G4EMLOW8.6.1"
export G4LEVELGAMMADATA="/home/manoja450/geant4-data/PhotonEvaporation5.7"
export G4RADIOACTIVEDATA="/home/manoja450/geant4-data/RadioactiveDecay5.6"
export G4PARTICLEXSDATA="/home/manoja450/geant4-data/G4PARTICLEXS4.1"
export G4PIIDATA="/home/manoja450/geant4-data/G4PII1.3"
export G4REALSURFACEDATA="/home/manoja450/geant4-data/RealSurface2.2"
export G4SAIDXSDATA="/home/manoja450/geant4-data/G4SAIDDATA2.0"
export G4ABLADATA="/home/manoja450/geant4-data/G4ABLA3.3"
export G4INCLDATA="/home/manoja450/geant4-data/G4INCL1.2"

echo "========================================="
echo "Simulation started at $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo "========================================="

# Create proper Geant4 input file
echo "/run/beamOn 1000" > ${OUTPUT_DIR}/input.txt

# Run simulation
./build/G4d2o < ${OUTPUT_DIR}/input.txt > ${OUTPUT_DIR}/simulation.log 2>&1

# Save a copy of the slurm output file
cp ${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.out \
   ${SLURM_DIR}/slurm-${SLURM_JOB_ID}_${TIMESTAMP}.out

echo "========================================="
echo "Simulation finished at $(date)"
echo "========================================="

grep "Number of events processed" ${OUTPUT_DIR}/simulation.log
