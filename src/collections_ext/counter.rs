use std::{collections::HashMap, hash::Hash};

pub struct Counter<K> {
    counts: HashMap<K, usize>,
}

impl<K> Counter<K>
where
    K: Eq + Hash,
{
    /// 创建一个新的 Counter。
    pub fn new() -> Self {
        Counter {
            counts: HashMap::new(),
        }
    }

    /// 获取某个键的计数值。
    pub fn get(&self, key: &K) -> usize {
        *self.counts.get(key).unwrap_or(&0)
    }

    /// 将指定的键的计数器加 1。
    pub fn increment(&mut self, key: K) {
        *self.counts.entry(key).or_insert(0) += 1;
    }

    /// 重置某个键的计数器为 0。
    pub fn reset(&mut self, key: K) {
        self.counts.insert(key, 0);
    }

    /// 清空所有计数。
    pub fn clear(&mut self) {
        self.counts.clear();
    }

    /// 获取所有计数的总和。
    pub fn total(&self) -> usize {
        self.counts.values().sum()
    }

    /// 获取出现次数最多的 n 个键值对。
    pub fn mostcommon(&self, n: usize) -> Vec<(&K, usize)> {
        let mut counts_vec: Vec<(&K, usize)> = self.counts.iter().map(|(k, &v)| (k, v)).collect();
        counts_vec.sort_by(|a, b| b.1.cmp(&a.1));
        counts_vec.into_iter().take(n).collect()
    }
}
