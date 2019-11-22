#!cargo watch -s 'cargo +nightly run --color always 2>&1 | bat -p --paging always'
#![feature(try_trait,try_blocks,type_ascription)]

fn main() -> anyhow::Result<()> {
    for filename in std::env::args().skip(1)  {
        let name = &filename.split('-').next().ok_or(anyhow::anyhow!("Expected -"))?;
        for r in gb_io::reader::SeqReader::new(std::fs::File::open(&filename)?) {
            let r = r?;
            let gene_name = "16S";
            let sequences : Vec<_> = r.features.iter().filter(|&f| {
                (try { f.kind == gb_io::feature_kind!("rRNA") && f.qualifier_values(gb_io::qualifier_key!("product")).next()?.contains(gene_name) } : Option<_>) == Some(true)
            }).map(|s| {
                let mut sequence = r.extract_location(&s.location).unwrap();
                sequence.push(0);
                sequence
            }).collect();
            let consensus_max_length = sequences.iter().fold(0, |max_v, v| std::cmp::max(max_v, v.len()));
            let consensus = rust_spoa::poa_consensus(&sequences, consensus_max_length, 1, 5, -4, -3, -1);
            let output_path = format!("{}-{}-consensus.fa", name, gene_name);
            {
                let mut output = std::io::BufWriter::new(std::fs::File::create(&output_path)?);
                use std::io::Write;
                writeln!(&mut output, ">{}_{}-consensus", name, gene_name)?;
                for chunk in consensus.chunks(60) { writeln!(&mut output, "{}", std::str::from_utf8(chunk)?.to_uppercase())?; }
            }
            print!("{}", std::str::from_utf8( &std::fs::read(output_path)? )? );
        }
    }
    Ok(())
}
