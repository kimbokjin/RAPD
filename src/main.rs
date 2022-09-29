use std::ptr::read;

use rapd::Genome;
use bio::io::fasta;

fn main() {
    let test_primer = vec![b"GGGATCCAAG".to_vec()];
    let genomes = read_genome_file("genome.fna").unwrap();
    println!("{:?}", genomes[0].desc);

    for genome in genomes.iter(){
        let a = genome.match_primers(test_primer.clone(), 0, 0, 0, 0);
    }

}

fn read_genome_file(filename: &str)  -> Result<Vec<Genome>, &str>  {
    let file : bool = filename.contains(".fasta") || filename.contains(".fna");
    let mut genome_file:Vec<Genome> = Vec::new();
    match file {
        true => {
            let sequence_file = fasta::Reader::from_file(filename).expect("Can't Loading Genome file."); 
            for sequence in sequence_file.records() {
                let record = sequence.expect("Error during fasta record parsing.");
                let seq = record.clone().seq().to_ascii_uppercase();
                let desc = String::from(record.clone().desc().unwrap_or(" "));
                let genome = Genome{
                    id: String::from(record.clone().id()),
                    desc,
                    seq,
                };
                genome_file.push(genome);
            }
            Ok(genome_file)
        }
        false => Err("Genome file format must be '.fasta' or '.fna'. "),
    } 
}