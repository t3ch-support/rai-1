// Generated by gencpp from file rai_msgs/MotionCompliance.msg
// DO NOT EDIT!


#ifndef RAI_MSGS_MESSAGE_MOTIONCOMPLIANCE_H
#define RAI_MSGS_MESSAGE_MOTIONCOMPLIANCE_H


#include <string>
#include <vector>
#include <map>

#include <ros/types.h>
#include <ros/serialization.h>
#include <ros/builtin_message_traits.h>
#include <ros/message_operations.h>

#include <rai_msgs/arr.h>

namespace rai_msgs
{
template <class ContainerAllocator>
struct MotionCompliance_
{
  typedef MotionCompliance_<ContainerAllocator> Type;

  MotionCompliance_()
    : K()
    , lambda(0.0)
    , xi(0.0)  {
    }
  MotionCompliance_(const ContainerAllocator& _alloc)
    : K(_alloc)
    , lambda(0.0)
    , xi(0.0)  {
  (void)_alloc;
    }



   typedef  ::rai_msgs::arr_<ContainerAllocator>  _K_type;
  _K_type K;

   typedef double _lambda_type;
  _lambda_type lambda;

   typedef double _xi_type;
  _xi_type xi;





  typedef boost::shared_ptr< ::rai_msgs::MotionCompliance_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::rai_msgs::MotionCompliance_<ContainerAllocator> const> ConstPtr;

}; // struct MotionCompliance_

typedef ::rai_msgs::MotionCompliance_<std::allocator<void> > MotionCompliance;

typedef boost::shared_ptr< ::rai_msgs::MotionCompliance > MotionCompliancePtr;
typedef boost::shared_ptr< ::rai_msgs::MotionCompliance const> MotionComplianceConstPtr;

// constants requiring out of line definition



template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const ::rai_msgs::MotionCompliance_<ContainerAllocator> & v)
{
ros::message_operations::Printer< ::rai_msgs::MotionCompliance_<ContainerAllocator> >::stream(s, "", v);
return s;
}

} // namespace rai_msgs

namespace ros
{
namespace message_traits
{



// BOOLTRAITS {'IsFixedSize': False, 'IsMessage': True, 'HasHeader': False}
// {'geometry_msgs': ['/opt/ros/kinetic/share/geometry_msgs/msg'], 'trajectory_msgs': ['/opt/ros/kinetic/share/trajectory_msgs/msg'], 'std_msgs': ['/opt/ros/kinetic/share/std_msgs/msg'], 'rai_msgs': ['/home/mtoussai/git/mlr/share/rai/rai/rai_msgs/msg']}

// !!!!!!!!!!! ['__class__', '__delattr__', '__dict__', '__doc__', '__eq__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_parsed_fields', 'constants', 'fields', 'full_name', 'has_header', 'header_present', 'names', 'package', 'parsed_fields', 'short_name', 'text', 'types']




template <class ContainerAllocator>
struct IsFixedSize< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
  : FalseType
  { };

template <class ContainerAllocator>
struct IsFixedSize< ::rai_msgs::MotionCompliance_<ContainerAllocator> const>
  : FalseType
  { };

template <class ContainerAllocator>
struct IsMessage< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::rai_msgs::MotionCompliance_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct HasHeader< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
  : FalseType
  { };

template <class ContainerAllocator>
struct HasHeader< ::rai_msgs::MotionCompliance_<ContainerAllocator> const>
  : FalseType
  { };


template<class ContainerAllocator>
struct MD5Sum< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
{
  static const char* value()
  {
    return "2124a4487f19e298feb6880df22639e8";
  }

  static const char* value(const ::rai_msgs::MotionCompliance_<ContainerAllocator>&) { return value(); }
  static const uint64_t static_value1 = 0x2124a4487f19e298ULL;
  static const uint64_t static_value2 = 0xfeb6880df22639e8ULL;
};

template<class ContainerAllocator>
struct DataType< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
{
  static const char* value()
  {
    return "rai_msgs/MotionCompliance";
  }

  static const char* value(const ::rai_msgs::MotionCompliance_<ContainerAllocator>&) { return value(); }
};

template<class ContainerAllocator>
struct Definition< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
{
  static const char* value()
  {
    return "arr K\n\
float64 lambda\n\
float64 xi\n\
\n\
================================================================================\n\
MSG: rai_msgs/arr\n\
uint32[] dim\n\
float64[] data\n\
";
  }

  static const char* value(const ::rai_msgs::MotionCompliance_<ContainerAllocator>&) { return value(); }
};

} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

  template<class ContainerAllocator> struct Serializer< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
  {
    template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
    {
      stream.next(m.K);
      stream.next(m.lambda);
      stream.next(m.xi);
    }

    ROS_DECLARE_ALLINONE_SERIALIZER
  }; // struct MotionCompliance_

} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::rai_msgs::MotionCompliance_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const ::rai_msgs::MotionCompliance_<ContainerAllocator>& v)
  {
    s << indent << "K: ";
    s << std::endl;
    Printer< ::rai_msgs::arr_<ContainerAllocator> >::stream(s, indent + "  ", v.K);
    s << indent << "lambda: ";
    Printer<double>::stream(s, indent + "  ", v.lambda);
    s << indent << "xi: ";
    Printer<double>::stream(s, indent + "  ", v.xi);
  }
};

} // namespace message_operations
} // namespace ros

#endif // RAI_MSGS_MESSAGE_MOTIONCOMPLIANCE_H
