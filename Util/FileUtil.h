#pragma once
#include "typedef.h"
#include <string>

class FileUtil {
public:
    FileUtil(const char* fileName);
    ~FileUtil();
    void save(const Field2d& f, const char* title, int time, int interval);

private:
    FILE* m_fprot;
};