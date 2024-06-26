/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#ifdef RAI_PYBIND

#include "ry-Frame.h"
#include "types.h"

#include "../Core/thread.h"
#include "../Kin/kin.h"
#include "../Kin/frame.h"
#include "../Kin/viewer.h"
#include "../Gui/opengl.h"
#include "../KOMO/skeletonSymbol.h"
#include "../Kin/simulation.h"

//void checkView(shared_ptr<rai::Frame>& self, bool recopyMeshes=false){
//  if(self->C.hasView()){
//    if(recopyMeshes) self->C.gl()->recopyMeshes(self->C);
//    self->C.view();
//  }
//}

void init_Frame(pybind11::module& m) {
  pybind11::class_<rai::Frame, shared_ptr<rai::Frame>>(m, "Frame", "A (coordinate) of a configuration, which can have a parent, and associated shape, joint, and/or inertia")

    .def("setColor", &rai::Frame::setColor )

    .def("setPose",  [](shared_ptr<rai::Frame>& self, const char* pose){
       self->setPose(rai::Transformation(pose));
     })
    .def("setPosition", &rai::Frame::setPosition )
    .def("setQuaternion", &rai::Frame::setQuaternion )
    .def("setRelativePose", [](shared_ptr<rai::Frame>& self, const char* pose){
       self->setRelativePose(rai::Transformation(pose));
     })
    .def("setRelativePosition", &rai::Frame::setRelativePosition )
    .def("setRelativeQuaternion", &rai::Frame::setRelativeQuaternion )
    .def("setJoint", &rai::Frame::setJoint )
    .def("setJointState", &rai::Frame::setJointState )
    .def("setContact", &rai::Frame::setContact )
    .def("setMass", &rai::Frame::setMass )
    .def("setPointCloud", [](std::shared_ptr<rai::Frame>& self, const pybind11::array& points, const pybind11::array_t<byte>& colors) {
	arr _points = numpy2arr<double>(points);
	byteA _colors = numpy2arr<byte>(colors);
	if(self->C.viewer()->gl){
	  auto mux = self->C.gl().dataLock(RAI_HERE);
	  self->setPointCloud(_points, _colors);
	}else{
	  self->setPointCloud(_points, _colors);
	}
     }, "", pybind11::arg("points"), pybind11::arg("colors") = pybind11::array_t<byte>{} )

    .def("setShape", &rai::Frame::setShape, "", pybind11::arg("type"), pybind11::arg("size") )

    .def("setParent", &rai::Frame::setParent, "", pybind11::arg("parent"), pybind11::arg("keepAbsolutePose_and_adaptRelativePose") = false, pybind11::arg("checkForLoop") = false)
    .def("unLink", &rai::Frame::unLink )

    
    .def("setAttribute", &rai::Frame::setAttribute )
    .def("addAttributes",  [](shared_ptr<rai::Frame>& self, const pybind11::dict& D){
      if(!self->ats) self->ats = make_shared<rai::Graph>();
      self->ats->copy(dict2graph(D), true);
     }, "add/set attributes for the frame")

    .def("getAttributes", [](shared_ptr<rai::Frame>& self){
      if(!self->ats) self->ats = make_shared<rai::Graph>();
      return graph2dict(*self->ats);
     }, "get frame attributes")

    .def_readwrite("name", &rai::Frame::name )

    .def("getPosition", &rai::Frame::getPosition )
    .def("getQuaternion", &rai::Frame::getQuaternion )
    .def("getRotationMatrix", &rai::Frame::getRotationMatrix )
    .def("getRelativePosition", &rai::Frame::getRelativePosition )
    .def("getRelativeQuaternion", &rai::Frame::getRelativeQuaternion )
    .def("getJointState", &rai::Frame::getJointState )
    .def("getSize", &rai::Frame::getSize )
    .def("getMeshPoints", &rai::Frame::getMeshPoints )
    .def("getMeshTriangles", &rai::Frame::getMeshTriangles )

    .def("info", [](shared_ptr<rai::Frame>& self) {
	rai::Graph G;
	G.add<rai::String>("name", self->name);
	G.add<int>("ID", self->ID);
	self->write(G);
	if(!G["X"]) G.add<arr>("X", self->ensure_X().getArr7d());
	return graph2dict(G);
      })

    .def("setMeshAsLines", [](shared_ptr<rai::Frame>& self, const std::vector<double>& lines) {
	//    CHECK(self.frame, "this is not a valid frame");
	CHECK(self->shape, "this frame is not a mesh!");
	CHECK_EQ(self->shape->type(), rai::ST_mesh, "this frame is not a mesh!");
	uint n = lines.size()/3;
	self->shape->mesh().V = lines;
	self->shape->mesh().V.reshape(n, 3);
	uintA& T = self->shape->mesh().T;
	T.resize(n/2, 2);
	for(uint i=0; i<T.d0; i++) {
	  T(i, 0) = 2*i;
	  T(i, 1) = 2*i+1;
	}
      })
    ;

}

//===========================================================================

void init_enums(pybind11::module& m){

#undef ENUMVAL
//#define ENUMVAL(x) m.attr(#x) = pybind11::int_(int(rai::x));
#define ENUMVAL(x) .value(#x, rai::x)

  pybind11::enum_<rai::ArgWord>(m, "ArgWord")
  ENUMVAL(_left)
  ENUMVAL(_right)
  ENUMVAL(_sequence)
  ENUMVAL(_path)
      .export_values();

#undef ENUMVAL
#define ENUMVAL(pre, x) .value(#x, pre##_##x)

  pybind11::enum_<rai::JointType>(m, "JT")
  ENUMVAL(rai::JT,hingeX) ENUMVAL(rai::JT,hingeY) ENUMVAL(rai::JT,hingeZ) ENUMVAL(rai::JT,transX) ENUMVAL(rai::JT,transY) ENUMVAL(rai::JT,transZ) ENUMVAL(rai::JT,transXY) ENUMVAL(rai::JT,trans3) ENUMVAL(rai::JT,transXYPhi) ENUMVAL(rai::JT,transYPhi) ENUMVAL(rai::JT,universal) ENUMVAL(rai::JT,rigid) ENUMVAL(rai::JT,quatBall) ENUMVAL(rai::JT,phiTransXY) ENUMVAL(rai::JT,XBall) ENUMVAL(rai::JT,free) ENUMVAL(rai::JT,generic) ENUMVAL(rai::JT,tau)
  .export_values();

  pybind11::enum_<rai::ShapeType>(m, "ST")
  ENUMVAL(rai::ST, none)
  ENUMVAL(rai::ST, box)
  ENUMVAL(rai::ST, sphere)
  ENUMVAL(rai::ST, capsule)
  ENUMVAL(rai::ST, mesh)
  ENUMVAL(rai::ST, cylinder)
  ENUMVAL(rai::ST, marker)
  ENUMVAL(rai::ST, pointCloud)
  ENUMVAL(rai::ST, ssCvx)
  ENUMVAL(rai::ST, ssBox)
  ENUMVAL(rai::ST, ssCylinder)
  ENUMVAL(rai::ST, ssBoxElip)
  ENUMVAL(rai::ST, quad)
  ENUMVAL(rai::ST, camera)
  ENUMVAL(rai::ST, sdf)
    .export_values();

  pybind11::enum_<FeatureSymbol>(m, "FS")
  ENUMVAL(FS, position)
  ENUMVAL(FS, positionDiff)
  ENUMVAL(FS, positionRel)
  ENUMVAL(FS, quaternion)
  ENUMVAL(FS, quaternionDiff)
  ENUMVAL(FS, quaternionRel)
  ENUMVAL(FS, pose)
  ENUMVAL(FS, poseDiff)
  ENUMVAL(FS, poseRel)
  ENUMVAL(FS, vectorX)
  ENUMVAL(FS, vectorXDiff)
  ENUMVAL(FS, vectorXRel)
  ENUMVAL(FS, vectorY)
  ENUMVAL(FS, vectorYDiff)
  ENUMVAL(FS, vectorYRel)
  ENUMVAL(FS, vectorZ)
  ENUMVAL(FS, vectorZDiff)
  ENUMVAL(FS, vectorZRel)
  ENUMVAL(FS, scalarProductXX)
  ENUMVAL(FS, scalarProductXY)
  ENUMVAL(FS, scalarProductXZ)
  ENUMVAL(FS, scalarProductYX)
  ENUMVAL(FS, scalarProductYY)
  ENUMVAL(FS, scalarProductYZ)
  ENUMVAL(FS, scalarProductZZ)
  ENUMVAL(FS, gazeAt)

  ENUMVAL(FS, angularVel)

  ENUMVAL(FS, accumulatedCollisions)
  ENUMVAL(FS, jointLimits)
  ENUMVAL(FS, distance)
  ENUMVAL(FS, negDistance)
  ENUMVAL(FS, oppose)

  ENUMVAL(FS, qItself)
  ENUMVAL(FS, jointState)

  ENUMVAL(FS, aboveBox)
  ENUMVAL(FS, insideBox)

  ENUMVAL(FS, pairCollision_negScalar)
  ENUMVAL(FS, pairCollision_vector)
  ENUMVAL(FS, pairCollision_normal)
  ENUMVAL(FS, pairCollision_p1)
  ENUMVAL(FS, pairCollision_p2)

  ENUMVAL(FS, standingAbove)

  ENUMVAL(FS, physics)
  ENUMVAL(FS, contactConstraints)
  ENUMVAL(FS, energy)

  ENUMVAL(FS, transAccelerations)
  ENUMVAL(FS, transVelocities)
      .export_values();

  pybind11::enum_<rai::SkeletonSymbol>(m, "SY")
  ENUMVAL(rai::SY,touch) ENUMVAL(rai::SY,above) ENUMVAL(rai::SY,inside) ENUMVAL(rai::SY,oppose) ENUMVAL(rai::SY,restingOn)
  ENUMVAL(rai::SY,poseEq) ENUMVAL(rai::SY,positionEq) ENUMVAL(rai::SY,stableRelPose) ENUMVAL(rai::SY,stablePose)
  ENUMVAL(rai::SY,stable) ENUMVAL(rai::SY,stableOn) ENUMVAL(rai::SY,dynamic) ENUMVAL(rai::SY,dynamicOn) ENUMVAL(rai::SY,dynamicTrans) ENUMVAL(rai::SY,quasiStatic) ENUMVAL(rai::SY,quasiStaticOn) ENUMVAL(rai::SY,downUp) ENUMVAL(rai::SY,break) ENUMVAL(rai::SY,stableZero)
  ENUMVAL(rai::SY,contact) ENUMVAL(rai::SY,contactStick) ENUMVAL(rai::SY,contactComplementary) ENUMVAL(rai::SY,bounce) ENUMVAL(rai::SY,push)
  ENUMVAL(rai::SY,magic) ENUMVAL(rai::SY,magicTrans)
  ENUMVAL(rai::SY,pushAndPlace)
  ENUMVAL(rai::SY,topBoxGrasp) ENUMVAL(rai::SY,topBoxPlace)
  ENUMVAL(rai::SY,dampMotion)
  ENUMVAL(rai::SY,identical)
  ENUMVAL(rai::SY,alignByInt)
  ENUMVAL(rai::SY,makeFree) ENUMVAL(rai::SY,forceBalance)
  ENUMVAL(rai::SY,relPosY)
  ENUMVAL(rai::SY,touchBoxNormalX) ENUMVAL(rai::SY,touchBoxNormalY) ENUMVAL(rai::SY,touchBoxNormalZ)
  ENUMVAL(rai::SY,boxGraspX) ENUMVAL(rai::SY,boxGraspY) ENUMVAL(rai::SY,boxGraspZ)
  ENUMVAL(rai::SY,lift)
  ENUMVAL(rai::SY,stableYPhi)
  ENUMVAL(rai::SY,stableOnX)
  ENUMVAL(rai::SY,stableOnY)
  ENUMVAL(rai::SY,end)
  .export_values();

#undef ENUMVAL
#define ENUMVAL(x) .value(#x, rai::Simulation::_##x)

  pybind11::enum_<rai::Simulation::Engine>(m, "SimulationEngine")
  ENUMVAL(physx)
  ENUMVAL(bullet)
  ENUMVAL(kinematic)
  .export_values();

  pybind11::enum_<rai::Simulation::ControlMode>(m, "ControlMode")
  ENUMVAL(none)
  ENUMVAL(position)
  ENUMVAL(velocity)
  ENUMVAL(acceleration)
  ENUMVAL(spline)
  .export_values();

  pybind11::enum_<rai::Simulation::ImpType>(m, "ImpType")
  ENUMVAL(closeGripper)
  ENUMVAL(moveGripper)
  ENUMVAL(depthNoise)
  ENUMVAL(rgbNoise)
  ENUMVAL(adversarialDropper)
  ENUMVAL(objectImpulses)
  ENUMVAL(noPenetrations)
  .export_values();

}

#endif
