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
