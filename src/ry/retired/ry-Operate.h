/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include <memory>

struct RobotOperation;

namespace ry {

struct RyOperate { std::shared_ptr<RobotOperation> R; };

}
