#pragma once
#include "typedef.h"

class FieldUtil {
public:
    static void display(const Field1d& f, int time, int interval);
    static void display(const Field2d& f, int time, int interval);
    static void setSize(Field2d& f, int x, int y);
};