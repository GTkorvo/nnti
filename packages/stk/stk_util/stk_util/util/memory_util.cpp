#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <stk_util/util/memory_util.hpp>

#if defined(__APPLE__)
#include<mach/task.h>
#include<mach/mach_init.h>
#endif
#ifdef __linux__
#define PROCFS
#endif

namespace stk
{

/* return current memory usage in megabytes */
size_t get_memory_usage_now()
{
  size_t memory = 0;

#if defined(__APPLE__)
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(),
      TASK_BASIC_INFO,
      (task_info_t)&t_info,
      &t_info_count))
  {
    return 0;
  }
  memory = t_info.resident_size;

#elif defined(__LIBCATAMOUNT__) && defined(REDSTORM_HEAP_INFO)
  size_t frags;
  unsigned long total_free, largest_free, total_used;
  heap_info( &frags, &total_free, &largest_free, &total_used );
  memory = total_used * 1024;

#elif defined(PROCFS)
  unsigned long rss_pages = 0;
  // Read memory size data from /proc/pid/stat
  // see "man proc" for details.
  std::ifstream proc_stat("/proc/self/stat");
  if (proc_stat)
  {
    std::string buf;
    proc_stat
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf
       >> rss_pages;
    proc_stat.close();
  }
  memory = rss_pages * sysconf( _SC_PAGESIZE);
#endif

  /* Success */
  return memory;
}

/* return current resident set size in megabytes */
void get_memory_usage(size_t & now, size_t & hwm)
{
#if defined (PROCFS)
  std::string vmhwm;
  std::string vmrss;
  std::string line(128,'\0');

  now = 0;
  hwm = 0;

  /* Read memory size data from /proc/pid/status
   * run "man proc" to get info on the contents of /proc/self/status
   */
  std::ifstream proc_status("/proc/self/status");
  if (!proc_status) return;

  while (vmhwm.empty() || vmrss.empty())
  {

    if(!std::getline(proc_status, line)) return;

    /* Find VmHWM */
    else if (line.substr(0, 6) == "VmHWM:")
    {
      vmhwm = line.substr(7);
      std::istringstream iss(vmhwm);
      iss >> hwm;
      hwm *= 1024;
    }
    /* Find VmRSS */
    else if (line.substr(0, 6) == "VmRSS:")
    {
      vmrss = line.substr(7);
      std::istringstream iss(vmrss);
      iss >> now;
      now *= 1024;
    }
  }
  proc_status.close();

  #else
  now = get_memory_usage_now();
  static size_t s_hwm = 0;
  s_hwm = std::max(now, s_hwm);
  hwm = s_hwm;
#endif
}

std::string human_bytes(size_t bytes)
{
  const size_t K = 1024;
  const size_t M = K*1024;
  const size_t G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << bytes << " M";
  } else {
    bytes /= G;
    out << bytes << " G";
  }
  return out.str();
}

} // namespace stk
