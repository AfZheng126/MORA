

#[derive(PartialEq, Clone, Copy)]
pub(crate) enum Rank {
    STRAIN, FORMA, VARIETAS, // strain sub group
    SUBSPECIES, SPECIES, SPECIESSUBGROUP, SPECIESGROUP, // species sub group 
    SUBGENUS, GENUS, SUPERGENUS, // genus sub group
    SUBTRIBE, TRIBE, SUBFAMILY, FAMILY, SUPERFAMILY, // family sub group
    PARVORDER, INFRAORDER, SUBORDER, ORDER, SUPERORDER, // order sub group
    COHORT, INFRACLASS, SUBCLASS, CLASS, SUPERCLASS, // class sub group
    SUBPHYLUM, PHYLUM, SUPERPHYLUM, // phylum sub group
    SUBKINGDOM, KINGDOM, SUPERKINGDOM , DOMAINE, LIFE, // general sub group
    UNCLASSIFIED
}

#[derive(PartialEq)]
enum ReadEnd {
    LEFT, RIGHT
}
#[derive(Clone, Copy)]
struct Interval {
    begin: u32,
    end: u32,
}

impl Interval {
    fn new(begin: u32, end:u32) -> Interval {
        Interval{begin, end}
    }
}

impl PartialEq for Interval{
    fn eq(&self, other: &Self) -> bool {
        self.begin == other.begin && self.end == other.end  
    }
}

#[derive(Clone)]
pub struct TaxaNode {
    // fixed for all the mappings (commiing from taxonomy tree structure)
    id: usize,
    parent_id: usize, // Maybe change this into an enum so I can write NO_PARENT as well. 
    pub(crate) rank: Rank,
    // change per read mapping
    score: u32,
    not_incorporated_children_counter: u64,
    active_children: Vec<usize>,
    left_fw: bool,
    right_fw: bool,
    l_pos: u32,
    r_pos: u32,
    l_intervals: Vec<Interval>,
    r_intervals: Vec<Interval>,
    is_cleaned: bool,
    paried: bool,
}

impl TaxaNode {

    pub(crate) fn new(id: usize, mut parent_id: usize, mut rank: Rank) -> TaxaNode {
        if id == parent_id {
            rank = Rank::LIFE;
            parent_id = usize::MAX;                       // the equivalent of -1 for unsigned int in C++ code
        }
        TaxaNode{ id, parent_id, rank, score: 0, 
            not_incorporated_children_counter: 0, active_children: Vec::new(), 
            left_fw: true, right_fw: false, l_pos: 0, r_pos: 0, 
            l_intervals: Vec::new(), r_intervals: Vec::new(), is_cleaned: false, paried: false }
    }

    pub(crate) fn new_with_id(id: usize) -> TaxaNode {
        TaxaNode{ id, parent_id: usize::MAX, rank: Rank::STRAIN, score: 0, 
            not_incorporated_children_counter: 0, active_children: Vec::new(), 
            left_fw: true, right_fw: false, l_pos: 0, r_pos: 0, 
            l_intervals: Vec::new(), r_intervals: Vec::new(), is_cleaned: false, paried: false }    }

    pub(crate) fn is_root(&self) -> bool {
        self.rank == Rank::LIFE
    }

    fn is_ripe(&self) -> bool {
        self.not_incorporated_children_counter != 0
    }
    
    // pub(crate) fn get_intervals(&self, read_end: ReadEnd) -> Vec<Interval> {
    //     match read_end {
    //         ReadEnd::LEFT => self.l_intervals,
    //         ReadEnd::RIGHT => self.r_intervals
    //     }
    // }

    pub(crate) fn get_id(&self) -> usize {
        self.id
    }

    pub(crate) fn get_parent_id(&self) -> usize {
        self.parent_id
    }

    // fn update_intervals(&self, child: TaxaNode, read_end: ReadEnd) {
    //     let mut intervals = Vec::new();

    //     // merge two sorted interval lists into the parent and then update the parent score
    //     let child_intervals = child.get_intervals(read_end);
    //     if read_end == ReadEnd::LEFT {
    //         intervals = self.l_intervals;
    //     } else {
    //         intervals = self.r_intervals;
    //     }
        
    // }

    fn get_pos(&self, read_end: ReadEnd) -> u32{
        match read_end {
            ReadEnd::LEFT => self.l_pos,
            ReadEnd::RIGHT => self.r_pos
        }
    }

    fn update_score(&mut self) {
        for it in &self.l_intervals {
            self.score += it.end - it.begin;
        }
        for it in &self.r_intervals {
            self.score += it.end - it.begin;
        }
    }

    // 1) sorts the intervals
    // 2) checks if the intervals can be merged
    // 3) calculates the score
    fn clean_intervals(&mut self, read_end: ReadEnd) {
        let mut intervals;
        if read_end == ReadEnd::LEFT {
            intervals = self.l_intervals.clone();
        } else {
            intervals = self.r_intervals.clone();
        }
        // step 1): sorts the intervals 
        intervals = sort(intervals);

        // step 2): merges the intervals if they overlap
        for i in 0..(intervals.len()-1) {
            let mut merged = true;
            let mut k = i + 1;
            while k != intervals.len() && merged {
                if intervals[i].end >= intervals[k].begin {
                    if intervals[i].end < intervals[k].end {
                        intervals[i].end = intervals[k].end;
                    }
                    intervals.remove(k); // deletes what has been merged from intervals
                    k += 1;
                } else {
                    // there is nothing to merge
                    merged = false;
                }
            }
        }
        // step 3): updates the score and intervals 
        self.update_score(); // For some reason, the original code did not include this line
        if read_end == ReadEnd::LEFT {
            self.l_intervals = intervals;
        } else {
            self.r_intervals = intervals;
        }
        self.is_cleaned = true;
    }

   

    // adds the child only if the child is not already added
    fn add_child(&mut self, child: TaxaNode) -> bool {
        if !self.active_children.contains(&child.id){
            self.active_children.push(child.id);
            self.not_incorporated_children_counter += 1;
            true
        } else { false }
    }

    fn is_concordant(&self) -> bool {
        self.l_intervals.len() == self.r_intervals.len()
    }

    // compares if the cleaned left and right intervals are the same
    fn compare_intervals(&mut self, other: &mut TaxaNode) -> bool {
        if self.is_cleaned {
            self.clean_intervals(ReadEnd::LEFT);
            self.clean_intervals(ReadEnd::RIGHT);
        }
        if !other.is_cleaned {
            other.clean_intervals(ReadEnd::LEFT);
            other.clean_intervals(ReadEnd::RIGHT);
        }
        self.l_intervals == other.l_intervals && self.r_intervals == other.r_intervals
    }

    // resets all the values in the TaxaNode to the default values
    fn reset(&mut self) {
        self.l_intervals.clear();
        self.r_intervals.clear();
        self.active_children.clear();
        self.not_incorporated_children_counter = 0;
        self.score = 0;
    }

    pub(crate) fn str_to_rank(rankstr: &String) -> Rank {
        if rankstr == "no rank" { Rank::STRAIN } 
        else if rankstr == "varietas" { Rank::VARIETAS} 
        else if rankstr == "subspecies" { Rank::SUBSPECIES}
        else if rankstr == "species" { Rank::SPECIES}
        else if rankstr == "species subgroup" { Rank::SPECIESSUBGROUP}
        else if rankstr == "species group" { Rank::SPECIESGROUP}
        else if rankstr == "subgenus" {Rank::SUBGENUS}
        else if rankstr == "genus" {Rank::GENUS}
        else if rankstr == "supergenus" {Rank::SUPERGENUS}
        else if rankstr == "subfamily" {Rank::SUBFAMILY}
        else if rankstr == "family" {Rank::FAMILY}
        else if rankstr == "superfamily" {Rank::SUPERFAMILY}
        else if rankstr == "subtribe" {Rank::SUBTRIBE}
        else if rankstr == "tribe" {Rank::TRIBE}
        else if rankstr == "forma" {Rank::FORMA}
        else if rankstr == "cohort" {Rank::COHORT}
        else if rankstr == "parvorder" {Rank::PARVORDER}
        else if rankstr == "suborder" {Rank::SUBORDER}
        else if rankstr == "order" {Rank::ORDER}
        else if rankstr == "infraorder" {Rank::INFRAORDER}
        else if rankstr == "superorder" {Rank::SUPERORDER}
        else if rankstr == "subclass" {Rank::SUBCLASS}
        else if rankstr == "class" {Rank::CLASS}
        else if rankstr == "infraclass" {Rank::INFRACLASS}
        else if rankstr == "superclass" {Rank::SUPERCLASS}
        else if rankstr == "subphylum" { Rank::SUBPHYLUM}
        else if rankstr == "phylum" { Rank::PHYLUM}
        else if rankstr == "superphylum"{ Rank::SUPERPHYLUM}
        else if rankstr == "subkingdom" { Rank::SUBKINGDOM}
        else if rankstr == "kingdom" { Rank::KINGDOM}
        else if rankstr == "superkingdom" { Rank::SUPERKINGDOM}
        else if rankstr == "domain" { Rank::DOMAINE}
        else if rankstr == "life" { Rank::LIFE}
        else {
            panic!("ERROR: Not a valid rank: {}\n", rankstr);
        }
    }

    pub (crate) fn rank_to_string(r: &Rank) -> String {
        match r {
            Rank::STRAIN => "no rank".to_string(),
            Rank::VARIETAS => "varietas".to_string(),
            Rank::SUBSPECIES => "subspecies".to_string(),
            Rank::SPECIES => "species".to_string(),
            Rank::SPECIESSUBGROUP => "species subgroup".to_string(),
            Rank::SPECIESGROUP => "species group".to_string(),
            Rank::SUBGENUS => "subgenus".to_string(),
            Rank::GENUS => "genus".to_string(),
            Rank::SUPERGENUS => "supergenus".to_string(),
            Rank::SUBFAMILY => "subfamily".to_string(),
            Rank::FAMILY => "family".to_string(),
            Rank::SUPERFAMILY => "superfamily".to_string(),
            Rank::SUBTRIBE => "subtribe".to_string(),
            Rank::TRIBE => "tribe".to_string(),
            Rank::FORMA => "forma".to_string(),
            Rank::COHORT => "cohort".to_string(),
            Rank::PARVORDER => "parvorder".to_string(),
            Rank::SUBORDER => "suborder".to_string(),
            Rank::ORDER => "order".to_string(),
            Rank::INFRAORDER => "infraorder".to_string(),
            Rank::SUPERORDER => "superorder".to_string(),
            Rank::SUBCLASS => "subclass".to_string(),
            Rank::CLASS => "class".to_string(),
            Rank::INFRACLASS => "infraclass".to_string(),
            Rank::SUPERCLASS => "superclass".to_string(),
            Rank::SUBPHYLUM => "subphylum".to_string(),
            Rank::PHYLUM => "phylum".to_string(),
            Rank::SUPERPHYLUM => "superphylum".to_string(),
            Rank::SUBKINGDOM => "subkingdom".to_string(),
            Rank::KINGDOM => "kingdom".to_string(),
            Rank::SUPERKINGDOM => "superkingdom".to_string(),
            Rank::DOMAINE => "domain".to_string(),
            Rank::LIFE => "life".to_string(),
            Rank::UNCLASSIFIED => "unclassified".to_string()
        }
    }

}

// Sorts the intervals so that i[n] i[n+1] satisfy i[n].begin <= i[n+1].begin and if they equal, then i[n].end <= i[n+1].end
fn sort(intervals: Vec<Interval>) -> Vec<Interval> {
    let mut sorted_interval:Vec<Interval> = Vec::new();
    // first we find the first element
    for i in 0..intervals.len() {
        let mut element = intervals[0];
        if !sorted_interval.contains(&intervals[i]) { 
            // if sorted_intervals doesn't contain it, we check if it is the earliest element not yet contained in sorted_intervals
            element = intervals[i];
            // check if there is another interval earlier than it
            for other in 0..intervals.len() {
                if intervals[other].begin <= element.begin {
                    if intervals[other].end < element.end {
                        element = intervals[other];
                    }
                }
            }
        }
        sorted_interval.push(element);
    }
    sorted_interval
}

pub struct TaxaInfo {
    cnt: u64,
    sub_tree_cnt: u64,
    rank: Rank,
}

impl TaxaInfo {
    fn new() -> TaxaInfo {
        TaxaInfo { cnt: 0, sub_tree_cnt: 0, rank: Rank::STRAIN }
    }

    fn new_with_values(cnt: u64, rank: Rank) -> TaxaInfo{
        TaxaInfo { cnt, sub_tree_cnt: 0, rank }
    }
}