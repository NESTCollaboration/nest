1.
Make a NEST build directory and a separate NEST install directory. From the build directory:
cmake -DG4=ON -DCMAKE_INSTALL_PREFIX=[path to install directory] ../relative/path/NobleElementSimulationTechnique
make; make install

2.
Include the following lines in the CMakeLists.txt of your GEANT4-based simulation project:
find_package(NEST REQUIRED)
include_directories(${NEST_INCLUDE_DIRS})

3.
In your simulation's build directory, use whatever cmake command you usually do to build the simulation, but include an additional flag:
-DNEST_DIR=[path to NEST install directory]

4.
You must set NESTStackingAction to be your simulation's Stacking Action. You probably already have a stacking action, which is set in a line that probably looks like:
SetUserAction(new myStackingAction());
If you already have a stacking action, you must modify it so it inherits from NESTStackingAction:
class myStackingAction : public NESTStackingAction
If you do not already have a stacking action, you still probably have class that inherits from G4VUserActionInitialization. In this class's Build() method, do:
SetUserAction(new NESTStackingAction());

5.
In your physics list, you'll want code that looks like this:
NEST::NESTProc* theNEST2ScintillationProcess = new NEST::NESTProc("S1",fElectromagnetic,[your electric field]);
if (theNEST2ScintillationProcess->IsApplicable(*particle)) {
  pmanager->AddProcess(theScintillationProcess, ordDefault+1, ordInActive, ordDefault+1);
}
