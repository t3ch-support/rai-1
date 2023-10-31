/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include "kin.h"
#include "cameraview.h"

namespace rai {

struct SimulationImp;
struct TeleopCallbacks;

//a non-threaded simulation with direct interface and stepping -- in constrast to BotSim, which is threaded (emulating real time) and has
//the default ctrl interface via low-level reference messages
struct Simulation {
  enum Engine { _noEngine, _physx, _bullet, _kinematic };
  enum ControlMode { _none, _position, _velocity, _acceleration, _posVel, _spline };
  enum ImpType { _closeGripper, _moveGripper, _depthNoise, _rgbNoise, _adversarialDropper, _objectImpulses, _blockJoints, _noPenetrations };

  std::unique_ptr<struct Simulation_self> self;

  Configuration& C;
  double time;
  Engine engine;
  Array<shared_ptr<SimulationImp>> imps; ///< list of (adversarial) imps doing things/perturbations/noise in addition to clean physics engine
  int verbose;
  FrameL grasps;
  std::shared_ptr<TeleopCallbacks> teleopCallbacks;

  Simulation(Configuration& _C, Engine _engine, int _verbose=2);
  ~Simulation();

  //== controller interface

  //-- step the simulation, optionally send a low-level control, or use the spline reference
  void step(const arr& u_control={}, double tau=.01, ControlMode u_mode = _spline);

  //-- adapt the spline reference to genreate motion (default way)
  void setSplineRef(const arr& _q, const arr& _times, bool append=true);

  //-- send a gripper command
  void moveGripper(const char* gripperFrameName, double width=.075, double speed=.3);
  void closeGripper(const char* gripperFrameName, double width=.05, double speed=.3, double force=20.);
  void closeGripperGrasp(const char* gripperFrameName, const char* objectName, double width=.05, double speed=.3, double force=20.);
  bool gripperIsDone(const char* gripperFrameName);

  //-- get state information
  const arr& get_q() { return C.getJointState(); }
  const arr& get_qDot();
  double getTimeToSplineEnd();
  double getGripperWidth(const char* gripperFrameName);
  bool getGripperIsGrasping(const char* gripperFrameName);
  bool getGripperIsClose(const char* gripperFrameName);
  bool getGripperIsOpen(const char* gripperFrameName);

  //-- get sensor information
  void getImageAndDepth(byteA& image, floatA& depth); ///< use this during stepping
  void getSegmentation(byteA& segmentation);
  CameraView& cameraview(); ///< use this if you want to initialize the sensor, etc
  rai::CameraView::Sensor& addSensor(const char* sensorName, const char* frameAttached=nullptr, uint width=640, uint height=360, double focalLength=-1., double orthoAbsHeight=-1., const arr& zRange= {}) {
    if(frameAttached && frameAttached[0]) {
      return cameraview().addSensor(sensorName, frameAttached, width, height, focalLength, orthoAbsHeight, zRange);
    } else {
      return cameraview().addSensor(sensorName);
    }
  }
  rai::CameraView::Sensor&  selectSensor(const char* name) { return cameraview().selectSensor(name); }
  byteA getScreenshot();


  //== ground truth interface
  rai::Frame* getGroundTruthFrame(const char* frame) { return C.getFrame("frame"); }


  //== perturbation/adversarial interface

  //what are really the fundamental perturbations? Their (reactive?) management should be realized by 'agents'. We need a method to add purturbation agents.
  void addImp(ImpType type, const StringA& frames, const arr& parameters);


  //== management interface

  //-- store and reset the state of the simulation
  void getState(arr& frameState, arr& q=NoArr, arr& frameVelocities=NoArr, arr& qDot=NoArr);
  void setState(const arr& frameState, const arr& q=NoArr, const arr& frameVelocities=NoArr, const arr& qDot=NoArr);
  void pushConfigurationToSimulator(const arr& frameVelocities=NoArr, const arr& qDot=NoArr);

  //-- post-hoc world manipulations
  void registerNewObjectWithEngine(rai::Frame* f);

  //allow writing pics for video
  uint& pngCount();

  Mutex& displayMutex();

  std::shared_ptr<struct PhysXInterface> hidden_physx();
  OpenGL& hidden_gl();
  void loadTeleopCallbacks();
};

//===========================================================================

struct TeleopCallbacks : OpenGL::GLClickCall, OpenGL::GLKeyCall, OpenGL::GLHoverCall{
  arr q_ref;
  bool stop=false;
  bool grab=false;
  arr oldx;
  double mouseDepth=0.;
  uint nMarkers=0;
  rai::Configuration& C;
  rai::Frame *marker=0;
  bool markerWasSet=false;

  TeleopCallbacks(rai::Configuration& C, rai::Frame *marker=0) : C(C), marker(marker) { q_ref = C.getJointState(); }

  bool hasNewMarker();
  virtual bool clickCallback(OpenGL& gl);
  virtual bool keyCallback(OpenGL& gl);
  virtual bool hoverCallback(OpenGL& gl);
};

}
