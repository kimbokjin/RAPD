use bio::alphabets::dna;
use std::str;

pub struct Genome {
    pub id: String,
    pub desc: String,
    pub seq: Vec<u8>,
}

impl Genome {
    pub fn new(id:String, desc:String, seq: Vec<u8>) -> Self {
        Genome { id, desc, seq }
    }

    pub fn match_primers(&self,primers: Vec<Vec<u8>>,total_mismatch:i32 ,mismatch: i32, perfect_match: i32, lower: i32, upper: i32 ) {
        let mut f_starting_point: Vec<i32> = Vec::new();
        let mut r_starting_point: Vec<i32> = Vec::new();

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
            let loop_length: usize = genome.len() as usize - primer_size; // consider only linear DNA
            let rev_comp_primer = dna::revcomp(primer);
            let mut rev_primer = primer.clone();
            let perfect_matching: usize = primer_size - perfect_match;
            rev_primer.reverse();
            //println!("{:?} {:?} {:?} {:?}",primer, rev_primer, rev_comp_primer, loop_length);

            for i in 0..loop_length {
                let partial_seqeunce = &genome[i..i+primer_size];
                if primer[perfect_matching..] == partial_seqeunce[perfect_matching..] {
                    let (forward_mismatch, forward_softmatch) =  guanine_hamming(primer, partial_seqeunce, true);
                    if forward_mismatch + forward_softmatch <= total_mismatch && forward_mismatch <= mismatch { //처음 3개만 보는 코드로 변경
                        // println!("Start Point:{}, # of Mismatch:{}, # of softmatch:{}",i, forward_mismatch, forward_softmatch );
                        // println!("{:?}  {:?}", &primer, &partial_seqeunce );
                        f_starting_point.push(i as i32 + 1)
                    }
                }
                
                if rev_comp_primer[..perfect_match] == partial_seqeunce[..perfect_match]{
                    let (reverse_mismatch, reverse_softmatch) = guanine_hamming(&rev_primer, partial_seqeunce, false);
                    if reverse_mismatch + reverse_softmatch <= total_mismatch && reverse_mismatch <= mismatch {
                        // let r_ref = print_u8_to_sequence(partial_seqeunce);
                        // let r_primer = print_u8_to_sequence(&rev_primer);
                        // println!("Start Point:{}, # of Mismatch:{}, # of softmatch:{}",i, reverse_mismatch, reverse_softmatch );
                        // println!("{}\n{}",r_primer,r_ref);
                        r_starting_point.push(i as i32 + 1)
                    }
                } 

            }
            let (forward_starting, reverse_starting, product_length) = calculate_distance(&f_starting_point, &r_starting_point, lower, upper);
            
            println!("{:?}", f_starting_point);

            let refe = print_u8_to_sequence(&genome[330402..330412]);
            let primer = print_u8_to_sequence(primer);
            println!("{}\n{}",primer,refe);
        }
    }
}


pub fn guanine_hamming(query: &[u8], reference: &[u8], isforward: bool)-> (i32,i32) {
    let mut mismatch: i32 = 0;
    let mut soft_match: i32 = 0;
    let mut boolean = true;
    
    let query_size = query.len();
    let reference_size = reference.len();

    let a = &b"A"[0]; //65
    let g = &b"G"[0]; //71
    let t = &b"T"[0]; //84
    let c = &b"C"[0]; //67

    // let second_base_q = &query[query_size - 2..query_size - 1][0];
    // let second_base_r = &reference[reference_size - 2..reference_size - 1][0];
    // let iter = query.iter().zip(reference.iter());
    if query.len() != reference.len() {
        // error
    } else {
        if isforward {
            for index in 0..query_size {
                let q = query[index];
                let r = reference[index];
                if q == r {
                    mismatch += 0;
                }else if (q == *g && r == *c) || (q == *g && r == *a) || (q == *g && r == *t)   {
                    mismatch += 0;
                    soft_match +=1;
                } 
                else {
                    mismatch += 1;
                }
            }           
        }
        if !isforward {
            for index in 0..query_size {
                let q = query[index];
                let r = reference[index];
                if (q == *a && r == *t)
                    || (q == *t && r == *a)
                    || (q == *g && r == *c)
                    || (q == *c && r == *g)

                {
                    mismatch += 0;
                }else if (q == *g && r == *t) || (q == *g && r == *g)
                /*|| (q == *t && r == *g) */ {
                    mismatch += 0;
                    soft_match += 1;

                } 
                else {
                        mismatch += 1
                }
            }        
        }
    }
    (mismatch, soft_match)
}


pub fn calculate_distance(
    foward: &Vec<i32>,
    reverse: &Vec<i32>,
    lower: i32,
    upper: i32,
) -> (Vec<i32>, Vec<i32>, Vec<i32>) {
    // upper lower type 변경 필요 -> gui할거면..
    let mut length: Vec<i32> = Vec::new();
    let mut forward_start: Vec<i32> = Vec::new();
    let mut reverse_start: Vec<i32> = Vec::new();
    for f in foward.iter() {
        for r in reverse.iter() {
            let distance = r - f;
            if distance >= lower && distance <= upper {
                length.push(distance);
                forward_start.push(*f);
                reverse_start.push(*r);
            }
        }
    }
    (forward_start, reverse_start, length)
}



pub fn print_u8_to_sequence(byte: &[u8]) -> String {
    let sequence = str::from_utf8(&byte)
        .expect("Found invaild UTF-8")
        .to_string();
    sequence
}