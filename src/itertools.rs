pub fn sliding_window(
    length: usize,
    win_len: usize,
    win_ovlp: usize,
    drop_last: bool,
) -> impl Iterator<Item = (usize, usize)> {
    assert!(win_len > win_ovlp);
    let step = win_len - win_ovlp;

    (0..)
        .map(move |i| i * step)
        .take_while(move |&start| start < length)
        .map(move |start| {
            let end = (start + win_len).min(length);
            (start, end)
        })
        .filter(move |&(start, end)| !drop_last || (end - start == win_len))
}


