#[macro_use]
extern crate lazy_static;
extern crate rayon;
use rayon::prelude::*;

use std::collections::HashMap;
use std::io::{self, BufRead};

// Create the lookup table for the complements
// See https://github.com/rust-bio/rust-bio/blob/0b8e1e9e012b0ea137e33ca1fc1e5c5cc4bd4f15/src/alphabets/dna.rs#L37
lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }

        for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32; // lowercase variants
        }
        comp
    };
}

/// Return complement of given DNA alphabet char (IUPAC alhpabet supported)
fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

/// Create a reverse complement of the input sequence
// TODO: fix this signature to take an &Sequence instead
fn revcomp(seq: &str) -> Sequence {
    String::from_utf8_lossy(
        seq.as_bytes()
            .iter()
            .rev()
            .map(|a| complement(*a))
            .collect::<Vec<u8>>()
            .as_ref(),
    )
    .into_owned()
}

/// Type alias for the counter
type Counter = HashMap<String, usize>;
/// Type alias of a sequence
type Sequence = String;

/// Trait for FASTA readers
pub trait FastaRead {
    fn read(&mut self, sequence: &mut Sequence) -> io::Result<()>;
}

/// Fasta reader
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    line: String,
}

impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }
}

/// Based on the rust bio fasta reader
/// https://docs.rs/bio/0.30.0/src/bio/io/fasta.rs.html#139
impl<R: io::Read> FastaRead for Reader<R> {
    fn read(&mut self, seq: &mut Sequence) -> io::Result<()> {
        seq.clear();
        // Check if at end of file
        if self.line.is_empty() {
            self.reader.read_line(&mut self.line)?;
            if self.line.is_empty() {
                return Ok(());
            }
        }
        // Get past the header
        if !self.line.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }

        // Now we are at the sequence
        loop {
            self.line.clear();
            self.reader.read_line(&mut self.line)?;
            if self.line.is_empty() || self.line.starts_with('>') {
                break;
            }
            seq.push_str(self.line.trim_end());
        }
        Ok(())
    }
}

impl<R: io::Read> IntoIterator for Reader<R> {
    type Item = io::Result<Sequence>;
    type IntoIter = IterReader<R>;

    fn into_iter(self) -> Self::IntoIter {
        IterReader { reader: self }
    }
}

pub struct IterReader<R: io::Read> {
    reader: Reader<R>,
}

impl<R: io::Read> Iterator for IterReader<R> {
    type Item = io::Result<Sequence>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut seq = Sequence::new();
        match self.reader.read(&mut seq) {
            Ok(_) if seq.is_empty() => None,
            Ok(_) => Some(Ok(seq)),
            Err(err) => Some(Err(err)),
        }
    }
}

fn count_kmer(k: usize, seq: Sequence) -> Option<Counter> {
    let mut counter: Counter = Counter::new();
    if seq.len() < k {
        return None;
    }
    for i in 0..seq.len() - k + 1 {
        let kmer_for = &seq[i..(i + k)];
        if kmer_for.contains('N') {
            continue;
        }
        let kmer_rev = revcomp(kmer_for);
        // TODO: This should work in newer rust?
        let kmer = if kmer_for < &kmer_rev {
            kmer_for.to_owned() // make kmer a string
        } else {
            kmer_rev
        };

        // Finally add the kmer to our dictionay
        let count = counter.entry(kmer).or_insert(0);
        *count += 1;
    }
    Some(counter)
}

fn count_stdin(k: usize) -> io::Result<Counter> {
    let stdin = io::stdin();
    let reader = Reader::new(stdin);

    let counter = reader
        .into_iter()
        .par_bridge()
        .flat_map(|result| {
            let seq = result.unwrap(); // TODO handle this potential explosion
            count_kmer(k, seq)
        })
        .reduce(
            || Counter::new(),
            |mut a: Counter, b: Counter| {
                for (k, v) in b {
                    let count = a.entry(k).or_insert(0);
                    *count += v;
                }
                a
            },
        );

    //for result in reader.into_par_iter() {
    //let seq = result?;
    //count_kmer(&mut counter, k, &seq)
    //}

    // Todo: read the actual seqs and do the things
    //let mut seq = Sequence::new();
    //loop {
    //match reader.read(&mut seq) {
    //Ok(_) if seq.is_empty() => break,
    //Ok(_) => count_kmer(&mut counter, k, &seq),
    //Err(err) => return Err(err),
    //}
    //}
    Ok(counter)
}

fn print_hist(counter: &Counter) {
    let mut hist = vec![0; 256];
    for (_kmer, count) in counter {
        let idx = if *count > 255 {
            255 as usize
        } else {
            *count as usize
        };
        hist[idx] += 1;
    }
    for i in 0..256 {
        println!("{}\t{}", i, hist[i]);
    }
}

fn main() -> io::Result<()> {
    let counter = count_stdin(31)?;
    print_hist(&counter);
    Ok(())
}
