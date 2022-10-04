
pub(crate) struct Stats {
    global_prob_sum: usize,
    total_read_cnt: usize,
    total_multi_mapped_reads: usize,
    total_unmapped_reads: usize,
    conflicting: usize,
}

impl Stats {
    fn new() -> Stats {
        Stats { global_prob_sum: 0, total_read_cnt: 0, 
            total_multi_mapped_reads: 0, total_unmapped_reads: 0, 
            conflicting: 0 }
    }

    pub(crate) fn new_with_stats(total_read_cnt: usize, total_multi_mapped_reads: usize, total_unmapped_reads: usize) -> Stats {
        Stats { global_prob_sum: 0, total_read_cnt, 
            total_multi_mapped_reads, total_unmapped_reads, 
            conflicting: 0 }    
        }

    pub(crate) fn update(&mut self, s: &Stats) {
        self.global_prob_sum += s.global_prob_sum;
        self.total_read_cnt += s.total_read_cnt;
        self.total_multi_mapped_reads += s.total_multi_mapped_reads;
        self.total_unmapped_reads += s.total_unmapped_reads;
        self.conflicting += s.conflicting;
    }

    // prints some of the information stored in stats.
    pub(crate) fn print_stats(&self) {
        println!("Stats:");
        println!("# of mapped (and accepted) reads: {}", self.total_read_cnt);
        println!("global probsum: {}", self.global_prob_sum);
        println!("# of multi-mapped reads: {}", self.total_multi_mapped_reads);
        println!("# of conflicting reads: {}", self.conflicting);
        println!("# of unmapped reads: {}\n", self.total_unmapped_reads);
    }

    pub(crate) fn get_total_read_cnt(&self) -> usize{
        self.total_read_cnt
    }
    
    pub(crate) fn get_total_unmapped_reads(&self) -> usize {
        self.total_unmapped_reads
    }
}
