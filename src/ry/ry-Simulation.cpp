/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#ifdef RAI_PYBIND

#include "ry-Config.h"
#include "ry-Simulation.h"
#include "types.h"

#include "../Kin/frame.h"
#include "../Kin/simulation.h"
#include "../Geo/depth2PointCloud.h"

namespace rai {
struct SimulationState {
  arr frameState;
  arr frameVels;

  SimulationState(const arr& _frameState, const arr& _frameVels) : frameState(_frameState), frameVels(_frameVels) {}
};
}

void init_Simulation(pybind11::module& m) {
  pybind11::class_<rai::Simulation, std::shared_ptr<rai::Simulation>>(m, "Simulation", "A direct simulation interface to physics engines (Nvidia PhysX, Bullet) -- see https://marctoussaint.github.io/robotics-course/tutorials/simulation.html")

  .def(pybind11::init<rai::Configuration&, rai::Simulation::Engine, int>(), "create a Simulation that is associated/attached to the given configuration",
       pybind11::arg("C"),
       pybind11::arg("engine"),
       pybind11::arg("verbose") = 2 )

//  .def(pybind11::init([](shared_ptr<rai::Configuration>& C, rai::Simulation::Engine engine, int verbose) {
//    return make_shared<rai::Simulation>(*C, engine, verbose);
//  }))

  .def("step", &rai::Simulation::step,
       "",
       pybind11::arg("u_control"),
       pybind11::arg("tau") = .01,
       pybind11::arg("u_mode") = rai::Simulation::_velocity
      )

  .def("move", &rai::Simulation::move,
       "set the spline reference to genreate motion",

  .def("setSplineRef", &rai::Simulation::setSplineRef,
       "set the spline reference to genreate motion",
       pybind11::arg("path"),
       pybind11::arg("t"),
       pybind11::arg("append") = true
       )

  .def("getTimeToSplineEnd", &rai::Simulation::getTimeToSplineEnd)

  .def("get_q", &rai::Simulation::get_q)

  .def("get_qDot", &rai::Simulation::get_qDot)

  .def("openGripper", &rai::Simulation::moveGripper,
       "",
       pybind11::arg("gripperFrameName"),
       pybind11::arg("width") = .075,
       pybind11::arg("speed") = .3
      )

  .def("closeGripper", &rai::Simulation::closeGripper,
       "",
       pybind11::arg("gripperFrameName"),
       pybind11::arg("width") = .05,
       pybind11::arg("speed") = .3,
       pybind11::arg("force") = 20.
      )

  .def("getGripperWidth", &rai::Simulation::getGripperWidth,
       "",
       pybind11::arg("gripperFrameName")
      )

  .def("getGripperIsGrasping", &rai::Simulation::getGripperIsGrasping,
       "",
       pybind11::arg("gripperFrameName")
      )

  .def("getImageAndDepth", [](std::shared_ptr<rai::Simulation>& self) {
    byteA rgb;
    floatA depth;
    self->getImageAndDepth(rgb, depth);
    return pybind11::make_tuple(Array2numpy<byte>(rgb),
                                Array2numpy<float>(depth));
  })

//  .def("getSegmentation", [](std::shared_ptr<rai::Simulation>& self) {
//    byteA seg;
//    self->getSegmentation(seg);
//    return pybind11::array_t<byte>(seg.dim(), seg.p);
//  })

  .def("addSensor",  &rai::Simulation::addSensor,
       "",
       pybind11::arg("sensorName"),
       pybind11::arg("frameAttached") = std::string(),
       pybind11::arg("width") = 640,
       pybind11::arg("height") = 360,
       pybind11::arg("focalLength") = -1.,
       pybind11::arg("orthoAbsHeight") = -1.,
       pybind11::arg("zRange") = std::vector<double>()
      )
  .def("selectSensor",  &rai::Simulation::selectSensor,
       "",
       pybind11::arg("sensorName")
      )

  .def("getGroundTruthPosition", [](std::shared_ptr<rai::Simulation>& self, const char* frame) {
    rai::Frame* f = self->C.getFrame(frame);
    arr x = f->getPosition();
    return arr2numpy(x);
  })

  .def("getGroundTruthRotationMatrix", [](std::shared_ptr<rai::Simulation>& self, const char* frame) {
    rai::Frame* f = self->C.getFrame(frame);
    arr x = f->getRotationMatrix();
    return arr2numpy(x);
  })

  .def("getGroundTruthSize", [](std::shared_ptr<rai::Simulation>& self, const char* frame) {
    rai::Frame* f = self->C.getFrame(frame);
    arr x = f->getSize();
    return arr2numpy(x);
  })

  .def("addImp", &rai::Simulation::addImp)

  .def("getState", [](std::shared_ptr<rai::Simulation>& self) {
    arr X, V, x, v;
    self->getState(X, x, V, v);
    return pybind11::make_tuple(arr2numpy(X), arr2numpy(x), arr2numpy(V), arr2numpy(v));
  }, "returns a 4-tuple or frame state, joint state, frame velocities (linear & angular), joint velocities")

  .def("setState", &rai::Simulation::setState,
       "",
       pybind11::arg("frameState"),
       pybind11::arg("jointState") = NoArr,
       pybind11::arg("frameVelocities") = NoArr,
       pybind11::arg("jointVelocities") = NoArr
      )

  .def("pushConfigurationToSimulator", &rai::Simulation::pushConfigurationToSimulator,
       "set the simulator to the full (frame) state of the configuration",
       pybind11::arg("frameVelocities") = NoArr,
       pybind11::arg("jointVelocities") = NoArr
       )

  .def("depthData2pointCloud", [](std::shared_ptr<rai::Simulation>& self, const pybind11::array_t<float>& depth, const std::vector<double>& FxyCxy) {
    arr points;
    floatA _depth = numpy2arr<float>(depth);
    depthData2pointCloud(points, _depth, arr(FxyCxy, true));
    return arr2numpy(points);
  })

  .def("getScreenshot", &rai::Simulation::getScreenshot)

  .def("loadTeleopCallbacks", &rai::Simulation::loadTeleopCallbacks)


  ;

  pybind11::class_<rai::CameraView::Sensor, std::shared_ptr<rai::CameraView::Sensor>>(m, "CameraViewSensor");

}

#endif
