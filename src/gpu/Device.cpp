#include "gpu/Device.hpp"

namespace simall::gpu {

bool       is_available()   noexcept { return false; }
int        device_count()   noexcept { return 0; }
DeviceInfo query(int)       noexcept { return {"none", 0, 0, 0}; }
void       synchronize()    noexcept {}

}  // namespace simall::gpu
