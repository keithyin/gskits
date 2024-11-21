use std::env;
use std::{fmt::Debug, str::FromStr, time};
use uuid::Uuid;
use chrono::{Local, DateTime};

pub fn command_line_str() -> String {
    let args: Vec<String> = env::args().collect();

    args.join(" ")
}

/// parse the range str to Vector of tuple
/// range fmt: b1:e1,b2:e2,be3,b4:e4. corresponds to [[b1, e1], [b2, e2], [be3, be3], [b4, e4]]. left inclusive, right inclusive
#[allow(unused)]
pub fn range_parser<T>(range_str: &str) -> Vec<(T, T)>
where
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    range_str
        .trim()
        .split(",")
        .map(|single_range| {
            let items = single_range
                .trim()
                .split(":")
                .filter(|item| item.trim().len() != 0)
                .collect::<Vec<_>>();
            return if items.len() == 2 {
                (
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                    items[1].parse::<T>().expect(&format!("{}", items[1])),
                )
            } else {
                (
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                )
            };
        })
        .collect::<Vec<_>>()
}

pub struct Range<T> {
    range: Vec<(T, T)>,
}

impl<T> Range<T>
where
    T: PartialOrd,
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    pub fn new(range_str: &str) -> Self {
        Self {
            range: range_parser::<T>(range_str),
        }
    }

    pub fn within_range(&self, v: T) -> bool {
        let mut result = false;
        for (b, e) in &self.range {
            if *b <= v && v <= *e {
                result = true;
                break;
            }
        }
        result
    }
}



#[allow(unused)]
#[derive(Debug)]
pub struct ScopedTimer {
    iter: usize,
    elapsed: u128 // nano
}

#[allow(unused)]
impl ScopedTimer {
    pub fn new() -> Self{
        ScopedTimer{
            iter: 0,
            elapsed: 0
        }
    }

    pub fn perform_timing(&mut self) -> Timer<'_> {
        Timer::new(self)
    }

    pub fn reset(&mut self) {
        self.iter = 0;
        self.elapsed = 0;
    }

    ///
    /// None: iter/nano_sec
    /// 1000: iter/micro
    /// 1000_000: iter/milli
    /// 1000_000_000: iter/sec
    /// ...
    pub fn speed(&self, multiplier: Option<u128>) -> f64 {
        let multiplier = multiplier.unwrap_or(1);
        assert!(multiplier > 0);
        self.iter as f64 / (self.elapsed as f64 / multiplier as f64)
    }

}

impl ToString for ScopedTimer {
    fn to_string(&self) -> String {
        format!("{:.10} iter/nano_secs", (self.iter as f64 / self.elapsed as f64))
    }
}

pub struct Timer<'a> {
    scoped_timer: &'a mut ScopedTimer,
    instant: time::Instant,
}

#[allow(unused)]
impl<'a> Timer<'a> {
    fn new(scoped_timer: &'a mut ScopedTimer) -> Self {
        let instant = time::Instant::now();
        Timer { scoped_timer: scoped_timer, instant: instant }
    }

    pub fn done_with_cnt(mut self, cnt: usize) {
        let elapsed = self.instant.elapsed().as_nanos();
        self.scoped_timer.iter += cnt;
        self.scoped_timer.elapsed += elapsed;
    }
}






// the generated filename will be fname-${uuid5}-
pub fn generate_tmp_filename(fname: &str) -> String {
    let now: DateTime<Local> = Local::now();

    // 格式化为年-月-日 时:分:秒 的形式
    let formatted_str = now.format("%Y-%m-%d %H:%M:%S").to_string();
    let uuid_v5 = Uuid::new_v5(&Uuid::NAMESPACE_DNS, formatted_str.as_bytes());

    return if fname.contains(".") {
        let (prefix, suffix) = fname.rsplit_once(".").unwrap();
        format!("{}-{}.{}", prefix, uuid_v5, suffix)
    } else {

        format!("{}-{}", fname, uuid_v5)
    }

}
#[cfg(test)]
mod test {
    use super::range_parser;

    #[test]
    fn test_range_parser() {
        eprintln!("{:?}", range_parser::<i32>("1:3,2,4:10"));
        eprintln!("{:?}", range_parser::<f32>("1.0:3.2,2.1,4.1:10.4"));
    }
}
