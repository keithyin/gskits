use std::{
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
    thread,
    time::{Duration, Instant},
};

use psutil;
use systemstat::{self, Platform};

#[derive(Debug, PartialEq, Eq)]
enum SysMonState {
    Init,
    Running,
    Stopped,
}

/// System resources monitor.
/// use tracing to log the resource info

pub struct SysMon {
    state: Arc<Mutex<SysMonState>>,
    need_stop: Arc<AtomicBool>,
    monitor_interval: Duration,
    prog_name: Arc<String>,
}

impl SysMon {
    pub fn new(monitor_interval: Duration, prog_name: String) -> Self {
        Self {
            state: Arc::new(Mutex::new(SysMonState::Init)),
            need_stop: Arc::new(AtomicBool::new(false)),
            monitor_interval: monitor_interval,
            prog_name: Arc::new(prog_name),
        }
    }

    /// pid: pid that need to be monitered. if None, use the program id as default
    /// worker_threads: used to compute cpu_load per thread. if None, use num threads of this program
    pub fn start_monitor(&self, pid: Option<u32>, worker_threads: Option<usize>) {
        let _state = self.state.clone();
        let _need_stop = self.need_stop.clone();
        let _monitor_interval = self.monitor_interval.clone();
        let _prog_name = self.prog_name.clone();
        thread::spawn(move || {
            monitor(
                pid,
                worker_threads,
                _state.clone(),
                _need_stop.clone(),
                _monitor_interval,
                _prog_name,
            );
        });
    }
}

impl Drop for SysMon {
    fn drop(&mut self) {
        self.need_stop.store(true, Ordering::Relaxed);

        let start = Instant::now();
        while *self.state.lock().unwrap() != SysMonState::Stopped {
            thread::sleep(Duration::from_millis(500));

            if start.elapsed() > Duration::from_secs(60) {
                // avoid monitor thread unexpected stop
                break;
            }
        }
    }
}

fn bytes2_gb(num_bytes: usize) -> usize {
    num_bytes / 1024 / 1024 / 1024
}

struct _SysMonStateGuard {
    state: Arc<Mutex<SysMonState>>,
}

impl From<Arc<Mutex<SysMonState>>> for _SysMonStateGuard {
    fn from(value: Arc<Mutex<SysMonState>>) -> Self {
        assert!(*value.lock().unwrap() == SysMonState::Init);
        *value.lock().unwrap() = SysMonState::Running;

        Self { state: value }
    }
}

impl Drop for _SysMonStateGuard {
    fn drop(&mut self) {
        *self.state.lock().unwrap() = SysMonState::Stopped;
    }
}
fn monitor(
    pid: Option<u32>,
    worker_threads: Option<usize>,
    state: Arc<Mutex<SysMonState>>,
    need_stop_flag: Arc<AtomicBool>,
    monitor_interval: Duration,
    prog_name: Arc<String>,
) {
    let _sys_mon_state_guard: _SysMonStateGuard = state.into();
    let pid = pid.unwrap_or(std::process::id());

    let sys_stat = systemstat::System::new();
    let mut monitored_process = psutil::process::Process::new(pid);
    if monitored_process.is_err() {
        tracing::warn!("Failed to get smcCore process info")
    }

    let mut peak_tot_used_memory = 0;
    let mut peak_smc_core_used_memory = 0;

    let mut now = Instant::now();

    loop {
        thread::sleep(Duration::from_secs(1));
        if now.elapsed() < monitor_interval {
            if need_stop_flag.load(Ordering::Relaxed) {
                break;
            }
            continue;
        }
        now = Instant::now();

        let mem = sys_stat.memory().unwrap();
        let tot = bytes2_gb(mem.total.0 as usize);
        let free = bytes2_gb(mem.free.0 as usize);
        let used = tot - free;
        peak_tot_used_memory = peak_tot_used_memory.max(used);
        tracing::info!(
            "System Memory Info ::> TotMem: {}GB, Used: {}GB, Free: {}GB. PeakUsed: {}GB",
            tot,
            used,
            free,
            peak_tot_used_memory
        );

        if let Ok(load_measure) = sys_stat.cpu_load() {
            thread::sleep(Duration::from_secs(1));
            match load_measure.done() {
                Ok(all_core_loads) => {
                    let num_cores = all_core_loads.len();
                    let all_core_usage: f32 = all_core_loads
                        .iter()
                        .map(|single_core_load| (1.0 - single_core_load.idle) * 100.)
                        .sum();
                    let per_core_usage = all_core_usage / num_cores as f32;
                    tracing::info!("System Cpu Load Info ::> num_cores: {}, AllCoreUsage: {:.1}%, PerCoreUsage: {:.1}%", 
                        num_cores, all_core_usage, per_core_usage);
                }

                Err(e) => tracing::warn!(
                    "System Cpu load Info. get cpu load measure faield. msg: {:?}",
                    e
                ),
            }
        }

        // logging SmcCore memory usage and cpu percentage infomation
        if let Ok(ref mut monitored_process_) = monitored_process {
            match monitored_process_.memory_info() {
                Ok(mem) => {
                    let used = bytes2_gb(mem.rss() as usize);
                    peak_smc_core_used_memory = peak_smc_core_used_memory.max(used);
                    tracing::info!(
                        "{} Memory Info ::> Used: {}GB, PeakUsed: {}GB",
                        prog_name,
                        used,
                        peak_smc_core_used_memory
                    );
                }
                Err(err) => tracing::warn!("Failed to get memory info:{:?}", err),
            }

            monitored_process_.num_threads();

            let _res = monitored_process_.cpu_percent();
            thread::sleep(Duration::from_secs(1));
            match monitored_process_.cpu_percent() {
                Ok(percent_) => {
                    let worker_threads_ =
                        worker_threads.unwrap_or(monitored_process_.num_threads() as usize);
                    let per_core_percent = percent_ / worker_threads_ as f32;
                    tracing::info!("{} Cpu Load Info ::> NumThreads:{} CpuPercent: Tot:{:.1}%, PerCpuCore:{:.1}%", 
                    prog_name, worker_threads_, percent_, per_core_percent)
                }
                Err(e) => tracing::warn!("Get {} CpuPercent error. msg: {:?}", prog_name, e),
            }
        }
    }
}
