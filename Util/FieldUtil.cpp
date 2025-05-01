#include <iostream>
#include "FieldUtil.h"

void FieldUtil::display(const Field1d& f, int time, int interval) {
    if (time % interval != 0) {
        return;
    }

    printf("t:%d", time);
    for (const auto& value : f) {
        printf("%6.3f", value);
    }
    printf("\n");
}

void FieldUtil::display(const Field2d& f, int time, int interval) {
    if (time % interval != 0) {
        return;
    }

    printf("Iteration time:%d\n", time);
    for (const auto& row : f) {
        for (const auto& value : row) {
            printf("%6.3f", value);
        }
        printf("\n");
    }
    printf("\n");
}

void FieldUtil::setSize(Field2d& f, int x, int y) {
    f.resize(y);
    for (auto& row : f) {
        row.resize(x);
    }
}
