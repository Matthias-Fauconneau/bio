#[macro_use]
extern crate gb_io;
extern crate rust_spoa;

use std::io::Write;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for filename in std::env::args().skip(1)  {
        let name = &filename.split('-').next().ok_or("Expected -")?;
        for r in gb_io::reader::SeqReader::new(std::fs::File::open(&filename).expect(&filename)) {
            let r = r?;
            let gene_name = "16S";
            let sequence_features = r.features.iter().filter(|&f| {
                if f.kind != feature_kind!("rRNA") { return false; }
                let product = f.qualifier_values(qualifier_key!("product")).next().unwrap();
                product.contains(gene_name)
                });
            let sequences : Vec<_> = sequence_features.map(|s| {
            let mut sequence = r.extract_location(&s.location).unwrap();
            sequence.push(0);
            sequence
            }).collect();
            let consensus_max_length = sequences.iter().fold(0, |max_v, v| std::cmp::max(max_v, v.len()));
            let consensus = rust_spoa::poa_consensus(&sequences, consensus_max_length, 1, 5, -4, -3, -1);
            let output_path = format!("{}-{}-consensus.fa", name, gene_name);
            {
                let mut output = std::io::BufWriter::new(std::fs::File::create(&output_path)?);
                writeln!(&mut output, ">{}_{}-consensus", name, gene_name)?;
                for chunk in consensus.chunks(60) { writeln!(&mut output, "{}", std::str::from_utf8(chunk)?.to_uppercase())?; }
            }
            print!("{}", std::str::from_utf8( &std::fs::read(output_path)? )? );
        }
    }
    Ok(())
}
