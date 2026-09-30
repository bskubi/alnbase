//! Store the reference as `Seq` bytes in a file, and read it through a memory
//! map.
//!
//! # On-disk layout
//!
//! Written back-to-front so the writer never needs to seek: the data body is
//! streamed, then the per-contig table, then a fixed-size trailer.
//!
//! ```text
//! ┌─ data body ────────────────┬─ header body ──────┬─ header start ─┐
//! │ name0 seq0 name1 seq1 ...  │ 48 B × n_contigs   │ 24 B trailer   │
//! │                            │  seq_size    u64   │ n_contigs u64  │
//! │ names: raw UTF-8, no NUL   │  seq_offset  u64   │ version   u64  │
//! │ seqs:  1 byte/base (Seq)   │  name_size   u64   │ magic     u64  │
//! │   IUPAC bases only; no     │  name_offset u64   │                │
//! │   flag bit is ever set     │  md5      16 B     │                │
//! └────────────────────────────┴────────────────────┴────────────────┘
//! ```
//!
//! All integers are little-endian. Offsets are absolute from the start of the
//! file. Nothing is alignment-padded; the reader does unaligned reads via
//! `u64::from_le_bytes`, and `Seq` is `align_of == 1` so sequence slices can be
//! cast straight out of the mapping.
//!
//! Reading is backwards from EOF: trailer at `len-24`, table at
//! `len - 24 - n_contigs*48`, everything before that is the data body.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::iter::repeat_n;
use std::path::Path;

use md5::{Digest, Md5};
use memmap2::Mmap;
use rust_htslib::bgzf;

use crate::seq::Seq;

/// `b"REFRMM\x1a\n"` — the 0x1a stops `cat` on DOS, the \n catches CRLF mangling.
const MAGIC: u64 = u64::from_le_bytes(*b"REFRMM\x1a\n");
/// Refused unless it matches exactly; rebuild an index with `alnbase index`.
const VERSION: u64 = 3;
const TRAILER: usize = 24;
const ENTRY: usize = 48;

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

// ---------------------------------------------------------------- writing --

struct Entry {
    seq_size: u64,
    seq_off: u64,
    name_size: u64,
    name_off: u64,
    md5: [u8; 16],
}

/// The SAM specification's `@SQ M5` digest, fed one FASTA line at a time: MD5
/// over the sequence with every byte outside `!`..=`~` removed and lowercase
/// letters uppercased (SAMv1 section 1.2.1).
fn update_m5(hasher: &mut Md5, line: &[u8], scratch: &mut Vec<u8>) {
    scratch.clear();
    scratch.extend(line.iter().filter(|b| (b'!'..=b'~').contains(*b)).map(u8::to_ascii_uppercase));
    hasher.update(scratch);
}

/// Strip a trailing line terminator (handles both `\n` and `\r\n`).
fn chomp(mut l: &[u8]) -> &[u8] {
    while let Some((&last, rest)) = l.split_last() {
        if last == b'\n' || last == b'\r' {
            l = rest;
        } else {
            break;
        }
    }
    l
}

/// Convert a FASTA (plain, gzip, or BGZF) into the `Aref` format.
///
/// Single pass, streaming: memory use is one line plus the contig table,
/// regardless of genome size.
pub fn write_from_fasta<P: AsRef<Path>, Q: AsRef<Path>>(fasta: P, out: Q) -> io::Result<()> {
    // Written beside the destination and renamed into place only once complete,
    // so a failed run never truncates or half-replaces an existing index.
    let out = out.as_ref();
    let mut tmp_name = out.as_os_str().to_owned();
    tmp_name.push(".partial");
    let tmp = std::path::PathBuf::from(tmp_name);
    let result = write_index(fasta.as_ref(), &tmp).and_then(|()| std::fs::rename(&tmp, out));
    if result.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    result
}

fn write_index(fasta: &Path, out: &Path) -> io::Result<()> {
    // bgzf_open sniffs the input, so this is the same call for .fa and .fa.gz.
    let hts = bgzf::Reader::from_path(fasta)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut fa = BufReader::with_capacity(1 << 20, hts);
    let mut w = BufWriter::with_capacity(1 << 20, File::create(out)?);
    let mut hasher = Md5::new();
    let mut m5_scratch: Vec<u8> = Vec::with_capacity(1 << 16);
    let mut seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();

    let mut pos: u64 = 0; // bytes written so far == current file offset
    let mut entries: Vec<Entry> = Vec::new();
    let mut line: Vec<u8> = Vec::with_capacity(1 << 16);
    let mut conv: Vec<Seq> = Vec::with_capacity(1 << 16);
    let mut line_no: u64 = 0;
    let mut contig_name: Vec<u8> = Vec::new();

    loop {
        line.clear();
        if fa.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_no += 1;
        let l = chomp(&line);
        if l.is_empty() || l[0] == b';' {
            continue;
        }

        if l[0] == b'>' {
            // close out the previous contig
            if let Some(e) = entries.last_mut() {
                e.seq_size = pos - e.seq_off;
                e.md5 = std::mem::replace(&mut hasher, Md5::new()).finalize().into();
            }
            let name = l[1..]
                .split(|c: &u8| c.is_ascii_whitespace())
                .next()
                .unwrap_or(&[]);
            if name.is_empty() {
                return Err(invalid("FASTA record with an empty name"));
            }
            // A BAM names contigs, and a name that means two sequences cannot
            // be resolved; silently keeping one would call against the wrong one.
            if !seen.insert(name.to_vec()) {
                return Err(invalid(&format!(
                    "the FASTA has more than one contig named {}",
                    String::from_utf8_lossy(name)
                )));
            }
            contig_name = name.to_vec();
            let name_off = pos;
            w.write_all(name)?;
            pos += name.len() as u64;
            entries.push(Entry {
                seq_size: 0, // patched when the record ends
                seq_off: pos,
                name_size: name.len() as u64,
                name_off,
                md5: [0; 16], // set when the record ends
            });
        } else {
            if entries.is_empty() {
                return Err(invalid("sequence data before the first '>' header"));
            }
            update_m5(&mut hasher, l, &mut m5_scratch);
            conv.clear();
            for &b in l.iter().filter(|b| !b.is_ascii_whitespace()) {
                let Some(base) = Seq::FROM_FASTA[b as usize] else {
                    let e = entries.last().expect("sequence before a header is refused above");
                    let at = pos - e.seq_off + conv.len() as u64 + 1;
                    return Err(invalid(&format!(
                        "line {line_no}: {:?} at {}:{at} (1-based) is not an IUPAC base; a reference \
                         may hold only A C G T U R Y S W K M B D H V N, in either case",
                        b as char,
                        String::from_utf8_lossy(&contig_name),
                    )));
                };
                conv.push(base);
            }
            w.write_all(bytemuck::cast_slice(&conv))?;
            pos += conv.len() as u64;
        }
    }
    if let Some(e) = entries.last_mut() {
        e.seq_size = pos - e.seq_off;
        e.md5 = hasher.finalize().into();
    }

    // header body
    for e in &entries {
        for v in [e.seq_size, e.seq_off, e.name_size, e.name_off] {
            w.write_all(&v.to_le_bytes())?;
        }
        w.write_all(&e.md5)?;
    }
    // header start
    w.write_all(&(entries.len() as u64).to_le_bytes())?;
    w.write_all(&VERSION.to_le_bytes())?;
    w.write_all(&MAGIC.to_le_bytes())?;
    w.flush()?;
    Ok(())
}

// ---------------------------------------------------------------- reading --

struct Contig {
    seq_off: usize,
    len: usize,
    name_off: usize,
    name_len: usize,
    md5: [u8; 16],
}

pub struct Aref {
    mmap: Mmap,
    contigs: Vec<Contig>,
    by_name: HashMap<Box<str>, usize>,
    /// Read from the trailer and checked on open; reported by `alnbase info`.
    version: u64,
}

#[inline]
fn rd_u64(b: &[u8], at: usize) -> io::Result<u64> {
    let c = b.get(at..at + 8).ok_or_else(|| invalid("truncated file"))?;
    Ok(u64::from_le_bytes(c.try_into().unwrap()))
}

impl Aref {
    /// Map an existing `.aref` file.
    ///
    /// # Safety
    ///
    /// The mapping is UB if the file is modified or truncated while open.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        let n = mmap.len();
        if n < TRAILER {
            return Err(invalid("file too short to hold a trailer"));
        }
        if rd_u64(&mmap, n - 8)? != MAGIC {
            return Err(invalid("bad magic"));
        }
        // A newer index is refused because this binary cannot know what
        // changed. An older one is refused because its table layout differs:
        // version 1 entries had no MD5. Rebuilding is cheap.
        let version = rd_u64(&mmap, n - 16)?;
        if version > VERSION {
            return Err(invalid(&format!(
                "index format version {version}, but this alnbase understands \
                 {VERSION}; rebuild the index or use a newer alnbase"
            )));
        }
        if version != VERSION {
            return Err(invalid(&format!(
                "index format version {version}, but this alnbase understands \
                 {VERSION}; rebuild it with `alnbase index`"
            )));
        }
        let n_contigs = rd_u64(&mmap, n - 24)? as usize;

        let table = n_contigs
            .checked_mul(ENTRY)
            .and_then(|t| n.checked_sub(TRAILER + t))
            .ok_or_else(|| invalid("contig table overruns the file"))?;
        let data_end = table; // data body occupies [0, table)

        let mut contigs = Vec::with_capacity(n_contigs);
        let mut by_name = HashMap::with_capacity(n_contigs);
        for i in 0..n_contigs {
            let at = table + i * ENTRY;
            let len = rd_u64(&mmap, at)? as usize;
            let seq_off = rd_u64(&mmap, at + 8)? as usize;
            let name_len = rd_u64(&mmap, at + 16)? as usize;
            let name_off = rd_u64(&mmap, at + 24)? as usize;
            let md5: [u8; 16] = mmap[at + 32..at + 48].try_into().unwrap();

            let ok = |off: usize, l: usize| off.checked_add(l).is_some_and(|e| e <= data_end);
            if !ok(seq_off, len) || !ok(name_off, name_len) {
                return Err(invalid("contig entry points outside the data body"));
            }
            let name = std::str::from_utf8(&mmap[name_off..name_off + name_len])
                .map_err(|_| invalid("contig name is not valid UTF-8"))?;
            by_name.insert(name.into(), i);
            contigs.push(Contig { seq_off, len, name_off, name_len, md5 });
        }
        Ok(Self { mmap, contigs, by_name, version })
    }

    #[inline]
    pub fn version(&self) -> u64 {
        self.version
    }
    /// The contig's `@SQ M5` digest as lowercase hex.
    pub fn md5_hex(&self, tid: usize) -> Option<String> {
        self.contigs.get(tid).map(|c| c.md5.iter().map(|b| format!("{b:02x}")).collect())
    }
    #[inline]
    pub fn n_contigs(&self) -> usize {
        self.contigs.len()
    }
    #[inline]
    pub fn len_of(&self, tid: usize) -> Option<usize> {
        self.contigs.get(tid).map(|c| c.len)
    }
    #[inline]
    pub fn tid(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied()
    }

    pub fn name(&self, tid: usize) -> Option<&str> {
        let c = self.contigs.get(tid)?;
        // validated as UTF-8 in `open`
        Some(unsafe {
            std::str::from_utf8_unchecked(&self.mmap[c.name_off..c.name_off + c.name_len])
        })
    }

    /// The whole contig, borrowed straight out of the mapping.
    pub fn contig(&self, tid: usize) -> Option<&[Seq]> {
        let c = self.contigs.get(tid)?;
        Some(bytemuck::cast_slice(&self.mmap[c.seq_off..c.seq_off + c.len]))
    }

    /// Borrowed slice for a fully in-bounds region — no padding, no iterator.
    /// Returns `None` if any part of `[start, end)` falls off the contig.
    pub fn slice(&self, tid: usize, start: i64, end: i64) -> Option<&[Seq]> {
        let c = self.contigs.get(tid)?;
        if start < 0 || end < start || end as u64 > c.len as u64 {
            return None;
        }
        Some(&self.contig(tid)?[start as usize..end as usize])
    }

    /// Zero-indexed, half-open `[start, end)`. Coordinates may be negative or
    /// run past the contig end; anything off-contig yields [`Seq::PAD`].
    ///
    /// Allocation-free: the in-bounds middle is borrowed from the mapping and
    /// the flanks are `repeat_n`, so nothing is copied. `None` only when the
    /// reference has no contig with this `tid`.
    pub fn fetch(&self, tid: usize, start: i64, end: i64)
    -> Option<impl DoubleEndedIterator<Item = Seq> + Clone + '_> {
        let seq = self.contig(tid)?;
        let len = seq.len() as i64;
        let end = end.max(start);

        let total = end.saturating_sub(start);
        let left = 0i64.saturating_sub(start).clamp(0, total);
        let right = end.saturating_sub(len).clamp(0, total - left);

        let lo = start.clamp(0, len) as usize;
        let hi = end.clamp(0, len).max(start.clamp(0, len)) as usize;

        Some(
            repeat_n(Seq::PAD, left as usize)
                .chain(seq[lo..hi].iter().copied())
                .chain(repeat_n(Seq::PAD, right as usize)),
        )
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(tag: &str) -> Aref {
        let dir = crate::test_support::temp_dir();
        let fa = dir.join(format!("aref_{tag}.fa"));
        let mm = dir.join(format!("aref_{tag}.aref"));
        std::fs::write(&fa, b">chr1 some description\nACGTN\nacgtu\n>chr2\nTTTT\n").unwrap();
        write_from_fasta(&fa, &mm).unwrap();
        Aref::open(&mm).unwrap()
    }

    #[test]
    fn layout() {
        let r = roundtrip("layout");
        assert_eq!(r.n_contigs(), 2);
        assert_eq!(r.name(0), Some("chr1")); // description dropped
        assert_eq!(r.len_of(0), Some(10));
        assert_eq!(r.len_of(1), Some(4));
        assert_eq!(r.tid("chr2"), Some(1));
    }

    #[test]
    fn encoding() {
        let r = roundtrip("encoding");
        let s = r.contig(0).unwrap();
        assert_eq!(s[0], Seq::A);
        assert_eq!(s[4], Seq::N);
        assert_eq!(s[5], Seq::A); // lowercase folded
        assert_eq!(s[9], Seq::T); // 'u'
        assert!(s.iter().all(|b| b.0 & Seq::FLAGS == 0), "an index holds bases only");
    }

    /// A character that is not an IUPAC base stops indexing, naming where it is,
    /// and leaves no index behind.
    #[test]
    fn a_character_that_is_not_a_base_is_refused() {
        for (bad, at) in [("*", "chr2:3"), ("-", "chr2:3"), (".", "chr2:3"), ("0", "chr2:3"), ("X", "chr2:3")] {
            let fa = crate::test_support::temp_path("aref_bad", "fa");
            let mm = crate::test_support::temp_path("aref_bad", "aref");
            let _ = std::fs::remove_file(&mm);
            std::fs::write(&fa, format!(">chr1\nACGT\n>chr2 desc\nAC{bad}T\n")).unwrap();
            let e = write_from_fasta(&fa, &mm).unwrap_err().to_string();
            assert!(e.contains(&format!("line 4: '{bad}' at {at} (1-based)")), "{e}");
            assert!(!mm.exists(), "no index is left behind for {bad}");
        }
    }

    #[test]
    fn padding() {
        let r = roundtrip("padding");
        let g = |a, b| r.fetch(1, a, b).unwrap().collect::<Vec<_>>();
        assert_eq!(g(0, 4), vec![Seq::T; 4]);
        assert_eq!(g(-2, 2), vec![Seq::PAD, Seq::PAD, Seq::T, Seq::T]);
        assert_eq!(g(2, 6), vec![Seq::T, Seq::T, Seq::PAD, Seq::PAD]);
        assert_eq!(g(-1, 5), vec![Seq::PAD, Seq::T, Seq::T, Seq::T, Seq::T, Seq::PAD]);
        assert_eq!(g(-9, -5), vec![Seq::PAD; 4]); // entirely left of the contig
        assert_eq!(g(100, 104), vec![Seq::PAD; 4]); // entirely right
        assert_eq!(g(2, 2).len(), 0); // empty
        assert_eq!(g(3, 1).len(), 0); // inverted
        assert!(r.fetch(9, 0, 4).is_none()); // no such contig
    }

    #[test]
    fn borrowed_slice() {
        let r = roundtrip("borrowed_slice");
        assert_eq!(r.slice(1, 1, 3), Some(&[Seq::T, Seq::T][..]));
        assert_eq!(r.slice(1, -1, 3), None);
        assert_eq!(r.slice(1, 1, 99), None);
    }
    /// An index from a future alnbase is refused by version, not misread.
    /// The trailer holds the version in the 8 bytes before the magic.
    #[test]
    fn an_index_version_this_binary_does_not_know_is_refused() {
        let path = crate::test_support::temp_path("aref_version", "aref");
        let src = crate::test_support::temp_path("aref_version", "fa");
        std::fs::write(&src, ">c1\nACGTACGT\n").unwrap();
        write_from_fasta(&src, &path).unwrap();
        assert!(Aref::open(&path).is_ok(), "the index this binary wrote opens");

        let mut bytes = std::fs::read(&path).unwrap();
        let n = bytes.len();
        let newer = crate::test_support::temp_path("aref_version_newer", "aref");
        bytes[n - 16..n - 8].copy_from_slice(&(VERSION + 1).to_le_bytes());
        std::fs::write(&newer, &bytes).unwrap();
        let e = match Aref::open(&newer) { Err(e) => e.to_string(), Ok(_) => panic!("opened a newer index") };
        assert!(e.contains(&format!("version {}", VERSION + 1)), "{e}");
        assert!(e.contains("newer alnbase"), "{e}");

        let older = crate::test_support::temp_path("aref_version_older", "aref");
        bytes[n - 16..n - 8].copy_from_slice(&0u64.to_le_bytes());
        std::fs::write(&older, &bytes).unwrap();
        let e = match Aref::open(&older) { Err(e) => e.to_string(), Ok(_) => panic!("opened an older index") };
        assert!(e.contains("version 0") && e.contains("alnbase index"), "{e}");
    }

}