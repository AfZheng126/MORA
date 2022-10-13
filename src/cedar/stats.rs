
pub(crate) struct Stats {
    total_read_cnt: usize,
    total_multi_mapped_reads: usize,
    total_unmapped_reads: usize,
}

impl Stats {
    pub(crate) fn new_with_stats(total_read_cnt: usize, total_multi_mapped_reads: usize, total_unmapped_reads: usize) -> Stats {
        Stats {total_read_cnt, total_multi_mapped_reads, total_unmapped_reads}    
    }

    // prints some of the information stored in stats.
    pub(crate) fn print_stats(&self) {
        println!("Stats:");
        println!("# of mapped (and accepted) reads: {}", self.total_read_cnt);
        println!("# of multi-mapped reads: {}", self.total_multi_mapped_reads);
        println!("# of unmapped reads: {}\n", self.total_unmapped_reads);
    }

    pub(crate) fn get_total_read_cnt(&self) -> usize{
        self.total_read_cnt
    }
    
    pub(crate) fn get_total_unmapped_reads(&self) -> usize {
        self.total_unmapped_reads
    }
}
