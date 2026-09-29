//! Open a BAM, including one that arrives on a pipe.
//!
//! The only thing shared between the two subcommands besides the reference
//! store: plumbing, not logic. It lives in one place so `-` means the same
//! thing to both, and so the decision about what a pipe implies is made once.
//!
//! # What streaming costs
//!
//! Nothing this tool needs. Both subcommands make a single forward pass and
//! never seek, so a pipe is as good as a file. What a pipe does remove is the
//! chance to check anything before reading -- an index, a second pass to
//! measure something -- so any precondition has to be decided from the header,
//! which arrives first.

use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam::{Format, Header, Read, Reader, Writer};
use rust_htslib::htslib;

/// The conventional name for standard input or output.
pub const STREAM: &str = "-";

pub fn is_stream(p: &Path) -> bool {
    p.as_os_str() == STREAM
}

/// Open a BAM for reading; `-` reads standard input.
pub fn open_reader(path: &Path, threads: usize) -> Result<Reader> {
    let mut r = if is_stream(path) {
        Reader::from_stdin().context("reading BAM from standard input")?
    } else {
        Reader::from_path(path).with_context(|| format!("opening {}", path.display()))?
    };
    // Decompression threads help on a pipe as much as on a file.
    if threads > 1 {
        r.set_threads(threads)?;
    }
    Ok(r)
}

/// Tell a CRAM reader which parts of each record will be read, so it can skip
/// decoding the rest.
///
/// `required` is a union of htslib `SAM_*` bits; see
/// [`RecordField::sam_fields`](crate::record_field::RecordField::sam_fields).
/// A field left out comes back as a placeholder rather than an error -- an
/// empty name, a zero mapq -- which is why the bits come from one exhaustive
/// match rather than being assembled by hand at each call site.
///
/// `md_nm` says whether an `MD` or `NM` tag will be read. CRAM usually omits
/// both and htslib regenerates them on decode by comparing every aligned base
/// against the reference, which is wasted work for a run that reads neither.
///
/// Both settings are CRAM-only, and htslib ignores them on any other format,
/// so this is safe to call on whatever [`open_reader`] returned. Call it
/// before the first record is read.
pub fn limit_cram_decoding(r: &mut Reader, required: u32, md_nm: bool) -> Result<()> {
    r.set_cram_options(htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS, required)
        .context("restricting the CRAM fields to decode")?;
    if !md_nm {
        r.set_cram_options(htslib::hts_fmt_option_CRAM_OPT_DECODE_MD, 0)
            .context("turning off MD/NM regeneration")?;
    }
    Ok(())
}

/// Open a BAM for writing; `-` writes standard output.
///
/// Always BAM, never SAM, even to a pipe: the next stage is another tool, and
/// uncompressed text would be a strange default for a binary format. Use
/// `samtools view` if a human needs to look.
pub fn open_writer(path: &Path, header: &Header, threads: usize) -> Result<Writer> {
    let mut w = if is_stream(path) {
        Writer::from_stdout(header, Format::Bam).context("writing BAM to standard output")?
    } else {
        Writer::from_path(path, header, Format::Bam)
            .with_context(|| format!("creating {}", path.display()))?
    };
    if threads > 1 {
        w.set_threads(threads)?;
    }
    Ok(w)
}