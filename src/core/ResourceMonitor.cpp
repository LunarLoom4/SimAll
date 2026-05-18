// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ResourceMonitor.cpp
//
// NVML is loaded dynamically; we never link against libnvidia-ml directly.
// The handle and entry points are cached in function-local statics that
// initialise once on first call. Failure is silent — `gpus` is left empty.
// =============================================================================
#include "core/ResourceMonitor.hpp"

#include <fstream>
#include <sstream>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <Windows.h>
#  include <Psapi.h>
#  pragma comment(lib, "Psapi.lib")
#else
#  include <sys/resource.h>
#  include <sys/time.h>
#  include <unistd.h>
#  include <dlfcn.h>
#endif

namespace simall::core {

namespace {

// ---------- NVML dynamic loader ----------
// Just enough of the NVML API to enumerate devices and read memory.
using nvmlReturn_t      = int;
constexpr nvmlReturn_t  NVML_SUCCESS = 0;
struct nvmlDevice_st;
using  nvmlDevice_t     = nvmlDevice_st*;
struct nvmlMemory_t { unsigned long long total; unsigned long long free; unsigned long long used; };

struct Nvml {
    bool   ok = false;
#if defined(_WIN32)
    HMODULE handle = nullptr;
#else
    void*   handle = nullptr;
#endif
    nvmlReturn_t (*Init)()                                                = nullptr;
    nvmlReturn_t (*Shutdown)()                                            = nullptr;
    nvmlReturn_t (*DeviceGetCount)(unsigned int*)                         = nullptr;
    nvmlReturn_t (*DeviceGetHandleByIndex)(unsigned int, nvmlDevice_t*)   = nullptr;
    nvmlReturn_t (*DeviceGetMemoryInfo)(nvmlDevice_t, nvmlMemory_t*)      = nullptr;
    nvmlReturn_t (*DeviceGetName)(nvmlDevice_t, char*, unsigned int)      = nullptr;
};

Nvml& nvml() {
    static Nvml n = []{
        Nvml r;
#if defined(_WIN32)
        r.handle = ::LoadLibraryA("nvml.dll");
        auto sym = [&](const char* s){
            return r.handle ? reinterpret_cast<void*>(::GetProcAddress(r.handle, s)) : nullptr;
        };
#else
        r.handle = ::dlopen("libnvidia-ml.so.1", RTLD_NOW);
        if (!r.handle) r.handle = ::dlopen("libnvidia-ml.so", RTLD_NOW);
        auto sym = [&](const char* s){
            return r.handle ? ::dlsym(r.handle, s) : nullptr;
        };
#endif
        if (!r.handle) return r;
        r.Init                   = reinterpret_cast<decltype(r.Init)>(sym("nvmlInit_v2"));
        if (!r.Init) r.Init      = reinterpret_cast<decltype(r.Init)>(sym("nvmlInit"));
        r.Shutdown               = reinterpret_cast<decltype(r.Shutdown)>(sym("nvmlShutdown"));
        r.DeviceGetCount         = reinterpret_cast<decltype(r.DeviceGetCount)>(sym("nvmlDeviceGetCount_v2"));
        if (!r.DeviceGetCount)
            r.DeviceGetCount     = reinterpret_cast<decltype(r.DeviceGetCount)>(sym("nvmlDeviceGetCount"));
        r.DeviceGetHandleByIndex = reinterpret_cast<decltype(r.DeviceGetHandleByIndex)>(sym("nvmlDeviceGetHandleByIndex_v2"));
        if (!r.DeviceGetHandleByIndex)
            r.DeviceGetHandleByIndex = reinterpret_cast<decltype(r.DeviceGetHandleByIndex)>(sym("nvmlDeviceGetHandleByIndex"));
        r.DeviceGetMemoryInfo    = reinterpret_cast<decltype(r.DeviceGetMemoryInfo)>(sym("nvmlDeviceGetMemoryInfo"));
        r.DeviceGetName          = reinterpret_cast<decltype(r.DeviceGetName)>(sym("nvmlDeviceGetName"));
        if (!r.Init || !r.DeviceGetCount || !r.DeviceGetHandleByIndex || !r.DeviceGetMemoryInfo) {
            return r;
        }
        if (r.Init() != NVML_SUCCESS) return r;
        r.ok = true;
        return r;
    }();
    return n;
}

void readProcSelfStatus(std::uint64_t& rss, std::uint64_t& vm) {
#if !defined(_WIN32)
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            std::istringstream is(line.substr(6));
            std::uint64_t kb = 0; std::string unit;
            is >> kb >> unit;
            rss = kb * 1024ull;
        } else if (line.rfind("VmSize:", 0) == 0) {
            std::istringstream is(line.substr(7));
            std::uint64_t kb = 0; std::string unit;
            is >> kb >> unit;
            vm = kb * 1024ull;
        }
    }
#else
    (void)rss; (void)vm;
#endif
}

}  // namespace

bool ResourceMonitor::nvmlAvailable() { return nvml().ok; }

ResourceSample ResourceMonitor::sample() {
    ResourceSample s;
    s.capturedAt = std::chrono::system_clock::now();

#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    if (::GetProcessMemoryInfo(::GetCurrentProcess(),
                               reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc),
                               sizeof(pmc))) {
        s.rssBytes     = pmc.WorkingSetSize;
        s.virtualBytes = pmc.PrivateUsage;
    }
    FILETIME ftCreate, ftExit, ftKernel, ftUser;
    if (::GetProcessTimes(::GetCurrentProcess(), &ftCreate, &ftExit, &ftKernel, &ftUser)) {
        auto toSec = [](FILETIME ft) {
            ULARGE_INTEGER u; u.LowPart = ft.dwLowDateTime; u.HighPart = ft.dwHighDateTime;
            return double(u.QuadPart) * 1e-7;     // 100-ns units → seconds
        };
        s.kernelCpuSeconds = toSec(ftKernel);
        s.userCpuSeconds   = toSec(ftUser);
    }
#else
    readProcSelfStatus(s.rssBytes, s.virtualBytes);
    struct rusage ru{};
    if (::getrusage(RUSAGE_SELF, &ru) == 0) {
        s.userCpuSeconds   = double(ru.ru_utime.tv_sec) + double(ru.ru_utime.tv_usec) * 1e-6;
        s.kernelCpuSeconds = double(ru.ru_stime.tv_sec) + double(ru.ru_stime.tv_usec) * 1e-6;
    }
#endif

    auto& n = nvml();
    if (n.ok) {
        unsigned int count = 0;
        if (n.DeviceGetCount(&count) == NVML_SUCCESS) {
            for (unsigned int i = 0; i < count; ++i) {
                nvmlDevice_t dev = nullptr;
                if (n.DeviceGetHandleByIndex(i, &dev) != NVML_SUCCESS) continue;
                nvmlMemory_t mem{};
                if (n.DeviceGetMemoryInfo(dev, &mem) != NVML_SUCCESS) continue;
                GpuMemorySample g;
                g.totalBytes = mem.total;
                g.freeBytes  = mem.free;
                g.usedBytes  = mem.used;
                char name[96]{};
                if (n.DeviceGetName && n.DeviceGetName(dev, name, sizeof(name)) == NVML_SUCCESS)
                    g.name = name;
                s.gpus.push_back(std::move(g));
            }
        }
    }
    return s;
}

}  // namespace simall::core
