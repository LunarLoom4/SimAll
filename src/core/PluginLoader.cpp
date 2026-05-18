// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/PluginLoader.cpp
// =============================================================================
#include "core/PluginLoader.hpp"
#include "core/PluginRegistry.hpp"
#include "core/Logger.hpp"

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <Windows.h>
namespace {
inline void* sim_dl_open(const wchar_t* p)               { return ::LoadLibraryW(p); }
inline void  sim_dl_close(void* h)                       { ::FreeLibrary(reinterpret_cast<HMODULE>(h)); }
inline void* sim_dl_sym(void* h, const char* n)          { return reinterpret_cast<void*>(::GetProcAddress(reinterpret_cast<HMODULE>(h), n)); }
inline std::string sim_dl_error() {
    DWORD e = ::GetLastError();
    return "LoadLibrary failed (code " + std::to_string(e) + ")";
}
}
#else
#  include <dlfcn.h>
namespace {
inline void* sim_dl_open(const char* p)                  { return ::dlopen(p, RTLD_NOW | RTLD_LOCAL); }
inline void  sim_dl_close(void* h)                       { ::dlclose(h); }
inline void* sim_dl_sym(void* h, const char* n)          { return ::dlsym(h, n); }
inline std::string sim_dl_error()                        { const char* e = ::dlerror(); return e ? e : "unknown"; }
}
#endif

namespace simall::core {

using CreateFn = simall::plugins::IPlugin* (*)();
using AbiFn    = int (*)();

PluginLoader::~PluginLoader() {
    unloadAll();
}

std::size_t PluginLoader::load(const std::filesystem::path& file) {
    void* handle = nullptr;
#if defined(_WIN32)
    handle = sim_dl_open(file.wstring().c_str());
#else
    handle = sim_dl_open(file.string().c_str());
#endif
    if (!handle) {
        SIMALL_LOG_ERROR("plugin", "Failed to load ", file.string(), ": ", sim_dl_error());
        return 0;
    }

    auto abiFn = reinterpret_cast<AbiFn>(sim_dl_sym(handle, "GetPluginAbi"));
    if (!abiFn) {
        SIMALL_LOG_ERROR("plugin", file.string(), " missing GetPluginAbi symbol");
        sim_dl_close(handle);
        return 0;
    }
    int reportedAbi = 0;
    try { reportedAbi = abiFn(); }
    catch (...) {
        SIMALL_LOG_ERROR("plugin", file.string(), ": GetPluginAbi threw");
        sim_dl_close(handle);
        return 0;
    }
    if (reportedAbi != plugins::kPluginAbiVersion) {
        SIMALL_LOG_ERROR("plugin", file.string(), ": ABI mismatch (host=",
                         plugins::kPluginAbiVersion, " plugin=", reportedAbi, ")");
        sim_dl_close(handle);
        return 0;
    }

    auto createFn = reinterpret_cast<CreateFn>(sim_dl_sym(handle, "CreatePlugin"));
    if (!createFn) {
        SIMALL_LOG_ERROR("plugin", file.string(), " missing CreatePlugin symbol");
        sim_dl_close(handle);
        return 0;
    }

    std::unique_ptr<plugins::IPlugin> instance;
    try {
        instance.reset(createFn());
    } catch (const std::exception& e) {
        SIMALL_LOG_ERROR("plugin", file.string(), ": CreatePlugin threw: ", e.what());
        sim_dl_close(handle);
        return 0;
    } catch (...) {
        SIMALL_LOG_ERROR("plugin", file.string(), ": CreatePlugin threw");
        sim_dl_close(handle);
        return 0;
    }
    if (!instance) {
        SIMALL_LOG_ERROR("plugin", file.string(), ": CreatePlugin returned null");
        sim_dl_close(handle);
        return 0;
    }

    auto entry = std::make_unique<LoadedPlugin>();
    entry->path    = file;
    entry->handle  = handle;
    entry->name    = instance->name();
    entry->version = instance->version();
    entry->abi     = reportedAbi;
    try { instance->on_load(); }
    catch (const std::exception& e) {
        SIMALL_LOG_ERROR("plugin", entry->name, ": on_load threw: ", e.what());
        sim_dl_close(handle);
        return 0;
    } catch (...) {
        SIMALL_LOG_ERROR("plugin", entry->name, ": on_load threw");
        sim_dl_close(handle);
        return 0;
    }
    entry->plugin = std::move(instance);
    PluginRegistry::instance().register_auto(entry->plugin.get());

    std::lock_guard lk(mtx_);
    plugins_.emplace_back(std::move(entry));
    std::size_t cookie = plugins_.size();
    SIMALL_LOG_INFO("plugin", "Loaded ", plugins_.back()->name,
                    " v", plugins_.back()->version, " from ", file.string());
    return cookie;
}

std::size_t PluginLoader::scanDirectory(const std::filesystem::path& dir) {
    if (!std::filesystem::is_directory(dir)) return 0;
    std::size_t count = 0;
    for (auto const& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        auto ext = entry.path().extension().string();
        bool match = (ext == ".simallplugin")
#if defined(_WIN32)
            || ext == ".dll"
#elif defined(__APPLE__)
            || ext == ".dylib"
#else
            || ext == ".so"
#endif
            ;
        if (!match) continue;
        if (load(entry.path()) != 0) ++count;
    }
    return count;
}

bool PluginLoader::unload(std::size_t cookie) {
    std::unique_ptr<LoadedPlugin> taken;
    {
        std::lock_guard lk(mtx_);
        if (cookie == 0 || cookie > plugins_.size()) return false;
        taken = std::move(plugins_[cookie - 1]);
        if (!taken) return false;
    }
    try { if (taken->plugin) taken->plugin->on_unload(); }
    catch (...) { SIMALL_LOG_ERROR("plugin", taken->name, ": on_unload threw"); }
    if (taken->plugin) PluginRegistry::instance().unregister(taken->plugin.get());
    taken->plugin.reset();
    if (taken->handle) sim_dl_close(taken->handle);
    return true;
}

void PluginLoader::unloadAll() {
    std::vector<std::unique_ptr<LoadedPlugin>> taken;
    {
        std::lock_guard lk(mtx_);
        taken.swap(plugins_);
    }
    for (auto it = taken.rbegin(); it != taken.rend(); ++it) {
        if (!*it) continue;
        try { if ((*it)->plugin) (*it)->plugin->on_unload(); }
        catch (...) { SIMALL_LOG_ERROR("plugin", (*it)->name, ": on_unload threw"); }
        (*it)->plugin.reset();
        if ((*it)->handle) sim_dl_close((*it)->handle);
    }
}

std::vector<LoadedPlugin> PluginLoader::snapshot() const {
    std::lock_guard lk(mtx_);
    std::vector<LoadedPlugin> out;
    out.reserve(plugins_.size());
    for (auto const& p : plugins_) {
        if (!p) continue;
        LoadedPlugin copy;
        copy.path    = p->path;
        copy.handle  = p->handle;
        copy.name    = p->name;
        copy.version = p->version;
        copy.abi     = p->abi;
        // intentionally do not move plugin pointer in snapshot
        out.push_back(std::move(copy));
    }
    return out;
}

}  // namespace simall::core
