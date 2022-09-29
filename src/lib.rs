use bio::alphabets::dna;

pub struct Genome {
    pub id: String,
    pub desc: String,
    pub seq: Vec<u8>,
}

impl Genome {
    pub fn new(id:String, desc:String, seq: Vec<u8>) -> Self {
        Genome { id, desc, seq }
    }

    pub fn match_primers(&self,primers: Vec<Vec<u8>>, mismatch: i32, perfect_match: i32, lower: i64, upper: i64 ) {
        let mut f_starting_point: Vec<i64> = Vec::new();
        let mut r_starting_point: Vec<i64> = Vec::new();

        let genome = &self.seq;
        let perfect_match = perfect_match as usize;

        for primer in primers.iter() {
            let primer_size = primer.len() as usize;

            /* 
            let primer_size_float = primer.len() as f64;
            let mismatch_float = mismatch as f64;
            let similarity_limit = ((primer_size_float - mismatch_float) / primer_size_float) * 100.0 ;
            let similarity_limit = similarity_limit.trunc();
            println!("Similarity Limit {}%", similarity_limit);
            */

            let loop_length = genome.len() as usize - primer_size; // consider only linear DNA

            let rev_comp_primer = dna::revcomp(primer);
            let mut rev_primer = primer.clone();
            rev_primer.reverse();
            println!("{:?} {:?} {:?} {:?}",primer, rev_primer, rev_comp_primer, loop_length);


        }
    }
}