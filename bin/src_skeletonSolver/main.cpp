// Includes necessary headers from the KOMO framework and other required libraries
#include <KOMO/skeleton.h>
#include <KOMO/komo.h>
#include <Kin/viewer.h>
#include <PathAlgos/ConfigurationProblem.h>
#include <PathAlgos/RRT_PathFinder.h>
#include <KOMO/pathTools.h>
#include <Gui/opengl.h>

//===========================================================================

// Usage instructions for running the program
const char *USAGE =
    "\nUSAGE:  skeletonSolver <problem> -mode [path|waypoints|final] -samples [n] -collisions [true|false]"
    "\n        (set parameters in rai.cfg alternatively, see z.log.global for a log of all used options)"
    "\n";

//===========================================================================

// Main function
int main(int argc,char **argv){
  // Initialize command line arguments
  rai::initCmdLine(argc, argv);

  // Seed random number generator for reproducibility
  rnd.seed(0);

  // Print usage instructions
  cout <<USAGE <<endl;

  // Parse command line arguments or use default values
  rai::String problem = rai::getParameter<rai::String>("problem", STRING("none"));
  if(rai::argc>=2 && rai::argv[1][0]!='-') problem=rai::argv[1];

  // Set configuration and skeleton files based on the problem
  rai::String sktFile, confFile;
  if(problem=="none"){
    sktFile = rai::getParameter<rai::String>("sktFile", STRING("none"));
    confFile = rai::getParameter<rai::String>("confFile", STRING("none"));
  }else{
    confFile = problem+".g";
    sktFile = problem+".skt";
  }

  // Get parameters for the solver
  uint samples = rai::getParameter<uint>("samples", 10);
  bool collisions = rai::getParameter<bool>("collisions", false);
  bool ways = rai::getParameter<bool>("ways", true);
  bool rrt = rai::getParameter<bool>("rrt", false);
  bool path = rai::getParameter<bool>("path", false);
  double rrtStopEvals =  rai::getParameter<double>("rrtStopEvals", 10000);
  double rrtTolerance =  rai::getParameter<double>("rrtTolerance", .03);
  double rrtStepsize =  rai::getParameter<double>("rrtStepsize", .05);

  // Log the used parameters
  LOG(0) <<"used parameters: " <<rai::params();

  // Initialize the robot configuration
  rai::Configuration C;
  C.addFile(confFile);

  // Read the skeleton file and set collision detection if needed
  rai::Skeleton S;
  S.read(FILE(sktFile));
  S.collisions = collisions;
  cout <<S <<endl;

  // Create KOMO instances for different optimization problems
  std::shared_ptr<KOMO> komo_way = S.getKomo_waypoints(C);
  std::shared_ptr<KOMO> komo_path = S.getKomo_path(C);
  std::shared_ptr<KOMO> komo_final = S.getKomo_finalSlice(C);

  // Loop over the number of samples
  for(uint i=0;i<samples;i++){
    cout <<"=== SAMPLE " <<i <<" ===" <<endl;

    // Process waypoints for the robot
    komo_way->initRandom(0);
    komo_way->pathConfig.gl().setTitle("WAYPOINTS");
    komo_way->view(true, STRING("random init sample " <<i));
    if(!ways) continue;

    // Optimization process for waypoint-based planning
    {
      NLP_Solver sol;
      sol.setProblem(komo_way->nlp());
      auto ret = sol.solve();
      cout <<komo_way->report(false, true) <<endl;
      cout <<*ret <<endl;
      cout <<komo_way->getPath_qAll();

      // Visualization and interaction handling
      rai::wait(.1);
      komo_way->pathConfig.viewer()->raiseWindow();
      komo_way->view(true, STRING("solved sample " <<i <<"\n" <<*ret));
      if(!ret->feasible) continue;
      while(komo_way->view_play(true));
    }

    // Setup the path configuration for further planning
    komo_path->initWithWaypoints(komo_way->getPath_qAll());
    komo_path->pathConfig.gl().setTitle("PATH");
    komo_path->view(false, STRING("init path"));

    // RRT-based path planning
    if(rrt){
      arrA paths;
      for(uint t=0;t<komo_way->T;t++){
        rai::Configuration C;
        arr q0, qT;
        rai::Skeleton::getTwoWaypointProblem(t, C, q0, qT, *komo_way);

        // Configuration problem setup for RRT
        ConfigurationProblem cp(C, true, rrtTolerance);
        if(S.explicitCollisions.N) cp.setExplicitCollisionPairs(S.explicitCollisions);
        cp.computeAllCollisions = S.collisions;

        // Initialize and execute RRT
        RRT_PathFinder rrt(cp, q0, qT, rrtStepsize);
        if(S.verbose>1) rrt.verbose=S.verbose-2;
        rrt.verbose=2;
        rrt.maxIters=rrtStopEvals;

        // Solution handling
        arr sol = rrt.planConnect();
        if(sol.N){
          sol = path_resampleLinear(sol, komo_path->stepsPerPhase);
          komo_path->initPhaseWithDofsPath(t, C.getDofIDs(), sol, false);
          komo_path->view(false, STRING("rrt phase " <<t));
        }
      }
    }
    komo_path->view(true, STRING("path after rrts"));

    // KOMO optimization for path planning
    if(path){
      NLP_Solver sol;
      sol.setProblem(komo_path->nlp());
      auto ret = sol.solve();
      cout <<komo_path->report(false, true) <<endl;
      cout <<*ret <<endl;

      // Visualization and interaction
      rai::wait(.1);
      komo_path->pathConfig.viewer()->raiseWindow();
      komo_path->view(true, STRING("solved sample " <<i <<"\n" <<*ret));
      while(komo_path->view_play(true));
    }
  }

  // Wait for user input before exiting
  rai::wait();

  return 0;
}
