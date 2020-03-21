// Use fastq crate to read the gzipped file diretly and offload the reading to a new thread
// Do the counting sequentially in this thread

#[macro_use]
extern crate lazy_static;
use fastq::{parse_path_fa, Record};
use fnv::FnvHashMap;
use std::collections::HashMap;
use std::env::args;
use std::io::{self};

//type Counter = HashMap<Vec<u8>, usize>;
type Counter = FnvHashMap<Vec<u8>, usize>;

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
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|a| complement(*a))
        .collect::<Vec<u8>>()
}

fn count_kmer(counter: &mut Counter, k: usize, seq: &[u8]) {
    if seq.len() < k {
        return;
    }
    for i in 0..seq.len() - k + 1 {
        let kmer_for = &seq[i..(i + k)];
        if kmer_for.contains(&('N' as u8)) {
            continue;
        }
        let kmer_rev = revcomp(kmer_for);
        let kmer = if kmer_for < &kmer_rev {
            kmer_for.to_owned() // make kmer a string
        } else {
            kmer_rev
        };

        // Finally add the kmer to our dictionay
        let count = counter.entry(kmer).or_insert(0);
        *count += 1;
    }
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
    // Get the file from
    let mut counter = Counter::default();
    let filename = args().nth(1);
    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => None,
        Some(name) => Some(name),
    };
    parse_path_fa(path, |parser| {
        parser
            .each(|rec| {
                count_kmer(&mut counter, 31, rec.seq());
                true
            })
            .unwrap();
    })?;
    print_hist(&counter);
    std::mem::forget(counter);
    Ok(())
}
