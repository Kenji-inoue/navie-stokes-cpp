#pragma once
#include <vector>

enum class ObjectFlag { fluid, surface, inside };
using Value = double;
using Field1d = std::vector<Value>;
using Field2d = std::vector<Field1d>;
using FlagField = std::vector<std::vector<ObjectFlag>>;

struct Velocity2d {
    Field2d u;
    Field2d v;
};