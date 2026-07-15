# MODULE2-SIMULATION-My-codes
COMPILATION AND JOB SUBMIT IN COMPUTE NODE IN ORNL CLUSTER;
STEP1:Get a compute node
srun --partition=shortjobs --time=30:00 --mem=8G --cpus-per-task=4 --pty bash
STEP2. Rebuild there
cd ~/G4d2oM2NeutronCaptureVERBOSE

rm -rf build
mkdir build
cd build

cmake ..
make -j4

STEP3: Exit the interactive session
exit
STEP4: Submit your simulation
cd ~/G4d2oM2NeutronCaptureVERBOSE
sbatch simulation.sh
