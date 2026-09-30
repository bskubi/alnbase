//! Apply an overlap plan to real records.
//!
//! Everything that decides *what* to do is in [`crate::overlap`], which is pure
//! and tested. This file reads records, hands their shape to the planner, and
//! writes the result back as a valid SAM template.
//!
//! Two judgements are made here rather than there, because they are about the
//! file rather than the molecule.
//!
//! **What "drop a record" means.** A supplementary with nothing left is
//! not written. A primary cannot vanish: SAM requires exactly one primary line
//! per read. If another alignment of the same read survives, the planner marks
//! it `promote`, and its alignment is moved onto the primary's record, which
//! carries the whole read (the survivor's hard clips become soft clips); the
//! supplementary record itself is not written. If nothing of the read survives,
//! the primary is written the way SAM describes an unmapped read.
//!
//! **What a changed alignment invalidates.** Whenever a template has a record
//! clipped, dropped or promoted, its mate fields (RNEXT, PNEXT, TLEN, flags
//! 0x2/0x8/0x20, MC and any MQ) are rebuilt from the final primary records the
//! way `samtools fixmate` builds them, and the SA tags of a read whose
//! alignments changed are rebuilt from its surviving alignments.

use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{CompressionLevel, Header, HeaderView, Read, Record};

use crate::bam_io::{open_reader, open_writer};
use crate::cli_overlap::{OverlapArgs, StaleTags};
use crate::overlap::{self, Aln, CigarOp, RecordPlan, Unresolved, Verdict, QUAL_MISSING};
use crate::contig_map::ContigMap;
use crate::aref::Aref;
use crate::seq::Seq;

pub fn run(args: OverlapArgs) -> Result<()> {
    args.validate().map_err(|e| anyhow!(e))?;
    let opts = args.options();

    let mut reader = open_reader(&args.input, args.threads)?;
    // The header arrives first, which is the only chance a pipe gives us to
    // reject bad input before consuming it.
    let text = String::from_utf8_lossy(reader.header().as_bytes()).into_owned();
    overlap::check_grouping(&text).map_err(|e| anyhow!(e))?;

    let header = Header::from_template(reader.header());
    let view = reader.header().clone();
    let mut writer = open_writer(&args.output, &header, args.threads)?;
    // `validate` has already held the level to 0-9.
    writer
        .set_compression_level(CompressionLevel::Level(args.compression_level as u32))
        .context("setting the output compression level")?;

    let refr = match &args.refr {
        Some(p) => Some(Aref::open(p)?),
        None => None,
    };
    // Resolved once, by name, and length-checked: a BAM tid indexes the BAM's
    // @SQ list and a Aref tid indexes the FASTA's order, which are not the
    // same list. Doing it here also takes the per-record name decode and hash
    // out of `fix_tags`.
    let contigs = match &refr {
        Some(r) => {
            let m = ContigMap::build(&view, r)?;
            if let Some(w) = m.warning() {
                eprintln!("{w}");
            }
            Some(m)
        }
        None => None,
    };
    let mut group: Vec<Record> = Vec::new();
    let mut name: Vec<u8> = Vec::new();
    let mut stats = overlap::Stats::default();
    let ctx = Ctx { opts: &opts, args: &args, refr: refr.as_ref(), contigs: contigs.as_ref(), view: &view };

    for rec in reader.records() {
        let rec = rec?;
        if !rec.qname().eq(&name[..]) && !group.is_empty() {
            for r in resolve(std::mem::take(&mut group), &ctx, &mut stats)? {
                writer.write(&r)?;
            }
        }
        if group.is_empty() {
            name = rec.qname().to_vec();
        }
        group.push(rec);
    }
    if !group.is_empty() {
        for r in resolve(std::mem::take(&mut group), &ctx, &mut stats)? {
            writer.write(&r)?;
        }
    }

    let unresolved: u64 = stats.unresolved.iter().sum();
    if unresolved > 0 {
        eprintln!(
            "{unresolved} of {} templates had a fragment length that could not be \
             established, and were {}. See the `unresolved` section of the statistics \
             for the breakdown by cause.",
            stats.templates,
            match opts.on_unresolved {
                Unresolved::Pass => "passed through unmodified",
                Unresolved::Drop => "discarded",
            }
        );
    }
    if let Some(path) = &args.stats {
        write_stats(&stats, path)?;
    }
    Ok(())
}

/// Write the run statistics. `-` goes to standard error, not standard output,
/// which may be carrying the BAM.
fn write_stats(stats: &overlap::Stats, path: &std::path::Path) -> Result<()> {
    let tsv = stats.to_tsv();
    if crate::bam_io::is_stream(path) {
        eprint!("{tsv}");
        return Ok(());
    }
    std::fs::write(path, tsv).with_context(|| format!("writing {}", path.display()))
}

/// Everything a template is resolved against that does not change per template.
struct Ctx<'a> {
    opts: &'a overlap::Options,
    args: &'a OverlapArgs,
    refr: Option<&'a Aref>,
    contigs: Option<&'a ContigMap>,
    view: &'a HeaderView,
}

/// 0 for R1, 1 for R2, by the same rule the planner uses.
fn end_of(rec: &Record) -> usize {
    rec.is_last_in_template() as usize
}

/// What a record said before anything was edited, kept for the SA rebuild.
struct Orig {
    end: usize,
    ops: Vec<CigarOp>,
    md: Option<String>,
    nm: Option<i64>,
    had_sa: bool,
}

/// Plan one template and return the records to write, in input order.
fn resolve(group: Vec<Record>, ctx: &Ctx, stats: &mut overlap::Stats) -> Result<Vec<Record>> {
    let (opts, args) = (ctx.opts, ctx.args);
    // A secondary alignment is the same bases placed elsewhere -- an
    // alternative hypothesis, not a second observation -- so resolving it
    // against the primary is meaningless. Supplementaries do participate.
    let live: Vec<usize> = (0..group.len())
        .filter(|&i| !group[i].is_secondary() && !group[i].is_unmapped())
        .collect();

    let alns: Vec<Aln> = live.iter().map(|&i| to_aln(i, &group[i])).collect();
    let out = overlap::plan(&alns, opts);

    stats.add(&out);

    if args.report_unresolved {
        if let Verdict::Unresolved { reason, l } = &out.verdict {
            let q = String::from_utf8_lossy(group[0].qname()).into_owned();
            eprintln!("{q}: {} (best fragment length {l:?})", reason.as_str());
        }
    }

    // "Discard the template" means discard it: writing unmapped records for a
    // template nothing could be established about would leave the ambiguity in
    // the file under a different name.
    let discard_all = matches!(out.verdict, Verdict::Unresolved { .. })
        && opts.on_unresolved == Unresolved::Drop;
    if discard_all {
        return Ok(Vec::new());
    }

    let n = group.len();
    let mut plan_of: Vec<Option<&RecordPlan>> = vec![None; n];
    for p in &out.plans {
        plan_of[p.id] = Some(p);
    }

    // Which ends had an alignment move, and the records that anchor the
    // drop-and-promote case.
    let mut end_moved = [false; 2];
    let mut dropped_primary: [Option<usize>; 2] = [None; 2];
    let mut promoted: [Option<usize>; 2] = [None; 2];
    for i in 0..n {
        let Some(p) = plan_of[i] else { continue };
        let e = end_of(&group[i]);
        if p.clip_front > 0 || p.clip_back > 0 || p.drop || p.promote {
            end_moved[e] = true;
        }
        if p.drop && !group[i].is_supplementary() {
            dropped_primary[e] = Some(i);
        }
        if p.promote {
            promoted[e] = Some(i);
        }
    }
    let moved = end_moved[0] || end_moved[1];

    let orig: Vec<Orig> = if moved {
        group
            .iter()
            .map(|r| Orig {
                end: end_of(r),
                ops: to_ops(r),
                md: match r.aux(b"MD") {
                    Ok(Aux::String(s)) => Some(s.to_string()),
                    _ => None,
                },
                nm: aux_int(r, b"NM"),
                had_sa: r.aux(b"SA").is_ok(),
            })
            .collect()
    } else {
        Vec::new()
    };

    // The whole read, for an end whose primary has nothing left, gathered
    // before anything is clipped. The survivor's own copy of a base comes
    // first, since the consensus edited the aligned copy.
    let whole: Vec<Option<(Vec<u8>, Vec<u8>)>> = (0..2)
        .map(|e| {
            let d = dropped_primary[e]?;
            let first: Vec<usize> = promoted[e].into_iter().chain([d]).collect();
            let rest = (0..n)
                .filter(|&i| plan_of[i].is_some() && end_of(&group[i]) == e && !first.contains(&i));
            let order: Vec<usize> = first.iter().copied().chain(rest).collect();
            whole_read(order.iter().map(|&i| (&group[i], plan_of[i])))
        })
        .collect();

    let mut slots: Vec<Option<Record>> = Vec::with_capacity(n);
    let mut edits: Vec<Option<usize>> = vec![None; n];
    for (i, mut rec) in group.into_iter().enumerate() {
        if let Some(p) = plan_of[i] {
            let k = p.clip_front + p.clip_back + p.qual_set.len();
            if k > 0 || p.drop {
                edits[i] = Some(k);
            }
            if p.drop && rec.is_supplementary() {
                slots.push(None);
                continue;
            }
            // A dropped primary whose read has a survivor is rebuilt below
            // from the survivor's alignment, not unmapped.
            let replaced = p.drop && promoted[end_of(&rec)].is_some();
            if !p.is_empty() && !replaced {
                apply(&mut rec, p, args, ctx.refr, ctx.contigs)?;
            }
        }
        slots.push(Some(rec));
    }

    // Where each written alignment came from, for its NM in SA.
    let mut origin: Vec<usize> = (0..n).collect();
    for e in 0..2 {
        let (Some(d), Some(s)) = (dropped_primary[e], promoted[e]) else { continue };
        let (Some(prim), Some(sup)) = (slots[d].take(), slots[s].take()) else { continue };
        slots[d] = Some(graft(prim, sup, plan_of[s].expect("promoted records are planned"), whole[e].as_ref())?);
        edits[d] = Some(edits[d].unwrap_or(0) + edits[s].unwrap_or(0));
        edits[s] = None;
        origin[d] = s;
    }

    if moved {
        rebuild_sa(&mut slots, &orig, &origin, &plan_of, end_moved, ctx.view)?;
        sync_mates(&mut slots)?;
    }

    let mut written = Vec::with_capacity(n);
    for (i, slot) in slots.into_iter().enumerate() {
        if let Some(mut rec) = slot {
            tag(&mut rec, &out.verdict, edits[i], args)?;
            written.push(rec);
        }
    }
    Ok(written)
}

fn to_aln(id: usize, rec: &Record) -> Aln {
    let cigar = to_ops(rec);

    // htslib's 4-bit nibbles are bit-identical to `Seq`, so no decode step.
    let s = rec.seq();
    let seq = (0..s.len()).map(|i| Seq(s.encoded_base(i))).collect();

    Aln {
        id,
        is_r1: !rec.is_last_in_template(),
        is_reverse: rec.is_reverse(),
        is_supplementary: rec.is_supplementary(),
        mapq: rec.mapq(),
        tid: rec.tid(),
        pos: rec.pos(),
        cigar,
        seq,
        qual: rec.qual().to_vec(),
    }
}

/// SEQ (ASCII) and QUAL with the plan's base and quality edits applied.
///
/// Edits are absolute rather than relative, so applying them twice cannot
/// compound. Indices are in the record's own SEQ space, before any clip moves
/// anything, which is the space the planner spoke in. A record stored without
/// qualities (QUAL `*`) keeps none: one edited byte would make the string
/// neither `*` nor valid.
fn edited(rec: &Record, p: Option<&RecordPlan>) -> (Vec<u8>, Vec<u8>) {
    let mut qual = rec.qual().to_vec();
    let mut seq = rec.seq().as_bytes();
    if let Some(p) = p {
        if qual.first() != Some(&QUAL_MISSING) {
            for &(i, q) in &p.qual_set {
                if i < qual.len() {
                    qual[i] = q.min(93);
                }
            }
        }
        for &(i, b) in &p.base_set {
            if i < seq.len() {
                seq[i] = b.to_char() as u8;
            }
        }
    }
    (seq, qual)
}

fn apply(
    rec: &mut Record,
    p: &RecordPlan,
    args: &OverlapArgs,
    refr: Option<&Aref>,
    contigs: Option<&ContigMap>,
) -> Result<()> {
    let (mut seq, mut qual) = edited(rec, Some(p));
    let qname = rec.qname().to_vec();

    let clipped = p.clip_front > 0 || p.clip_back > 0;
    if p.drop {
        // Nothing of this read survives anywhere (a supplementary never
        // reaches here, and a primary with a survivor is grafted instead), so
        // it is written as SAM describes an unmapped read: back in the
        // orientation it was sequenced in, with no alignment, MAPQ 0 and no
        // proper-pair bit. Its placement beside its mate is set with the other
        // mate fields.
        if rec.is_reverse() {
            reverse_complement(&mut seq);
            qual.reverse();
            rec.unset_reverse();
        }
        rec.set(&qname, None, &seq, &qual);
        rec.set_unmapped();
        rec.unset_proper_pair();
        rec.set_mapq(0);
        strip_tags(rec);
    } else if clipped {
        let ops: Vec<CigarOp> = to_ops(rec);
        let (mut new_ops, new_pos) =
            overlap::reclip(&ops, rec.pos(), p.clip_front, p.clip_back);
        // A hard clip has to take SEQ and QUAL with it, or the record no longer
        // describes itself.
        if p.hard_front || p.hard_back {
            let (h, cut_front, cut_back) =
                overlap::harden(&new_ops, p.hard_front, p.hard_back);
            new_ops = h;
            let keep = cut_front..seq.len().saturating_sub(cut_back);
            if keep.start < keep.end {
                seq = seq[keep.clone()].to_vec();
                qual = qual[keep].to_vec();
            }
        }
        let cigar = CigarString(new_ops.iter().map(to_htslib).collect());
        rec.set(&qname, Some(&cigar), &seq, &qual);
        // Clipping from the front moves the alignment start; clipping from the
        // back does not. `reclip` returns the adjusted position.
        rec.set_pos(new_pos);
        fix_tags(rec, args, refr, contigs)?;
    } else {
        let cigar = rec.cigar().take();
        rec.set(&qname, Some(&cigar), &seq, &qual);
    }
    Ok(())
}

/// Put a promoted supplementary's alignment on its read's primary record.
///
/// The primary record is kept because it is the one that carries the whole
/// read and every read-level tag (read group, barcodes, base modifications);
/// only the alignment moves. The survivor's hard clips become soft clips over
/// bases taken from the whole read, except on a side this run hard-clipped
/// itself (`--clip-mode hard`), which stays hard. NM, MD and AS come from the
/// survivor, whose alignment they describe; SA and mate fields are rebuilt
/// afterwards. If the whole read is not available (its primary was itself
/// hard-clipped), the survivor becomes the primary as it stands.
fn graft(
    mut prim: Record,
    mut sup: Record,
    p: &RecordPlan,
    whole: Option<&(Vec<u8>, Vec<u8>)>,
) -> Result<Record> {
    let ops = to_ops(&sup);
    let len = overlap::read_length(&ops);
    let Some((wseq, wqual)) = whole.filter(|w| w.0.len() == len) else {
        sup.unset_supplementary();
        return Ok(sup);
    };
    let keep_front = p.hard_front && p.clip_front > 0;
    let keep_back = p.hard_back && p.clip_back > 0;
    let (lead, trail) = hard_clips(&ops);

    let (mut seq, mut qual) = (wseq.clone(), wqual.clone());
    if sup.is_reverse() {
        reverse_complement(&mut seq);
        qual.reverse();
    }
    let span = (if keep_front { lead } else { 0 })..(len - if keep_back { trail } else { 0 });
    let (seq, qual) = (seq[span.clone()].to_vec(), qual[span].to_vec());
    let new_ops = soften(&ops, !keep_front, !keep_back);
    let cigar = CigarString(new_ops.iter().map(to_htslib).collect());

    let qname = prim.qname().to_vec();
    prim.set(&qname, Some(&cigar), &seq, &qual);
    prim.set_tid(sup.tid());
    prim.set_pos(sup.pos());
    prim.set_mapq(sup.mapq());
    if sup.is_reverse() {
        prim.set_reverse();
    } else {
        prim.unset_reverse();
    }
    strip_tags(&mut prim);
    for t in [b"NM", b"MD", b"AS"] {
        if let Ok(v) = sup.aux(t) {
            prim.push_aux(t, v)?;
        }
    }
    Ok(prim)
}

/// A read as it came off the sequencer -- ASCII bases and qualities, 5' to 3'
/// -- assembled from the edited copies its records hold, earlier records
/// winning. `None` if a base is held by no record, or the records disagree
/// about the read's length.
fn whole_read<'a>(
    recs: impl Iterator<Item = (&'a Record, Option<&'a RecordPlan>)>,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let mut bases: Vec<Option<(u8, u8)>> = Vec::new();
    for (rec, p) in recs {
        let ops = to_ops(rec);
        let len = overlap::read_length(&ops);
        if bases.is_empty() {
            bases = vec![None; len];
        } else if bases.len() != len {
            return None;
        }
        let (lead, _) = hard_clips(&ops);
        let (seq, qual) = edited(rec, p);
        for (q, (&b, &bq)) in seq.iter().zip(&qual).enumerate() {
            let at = q + lead;
            if at >= len {
                break;
            }
            let (e, b) = if rec.is_reverse() { (len - 1 - at, complement(b)) } else { (at, b) };
            bases[e].get_or_insert((b, bq));
        }
    }
    if bases.is_empty() {
        return None;
    }
    let (seq, mut qual): (Vec<u8>, Vec<u8>) = bases.into_iter().collect::<Option<Vec<_>>>()?.into_iter().unzip();
    // Records that disagree about having qualities cannot make a valid QUAL
    // between them; `*` is the honest answer.
    if qual.contains(&QUAL_MISSING) {
        qual.iter_mut().for_each(|q| *q = QUAL_MISSING);
    }
    Some((seq, qual))
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other, // N, S, W, =
    }
}

fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    seq.iter_mut().for_each(|b| *b = complement(*b));
}

/// Total leading and trailing `H`.
fn hard_clips(ops: &[CigarOp]) -> (usize, usize) {
    let h = |o: &&CigarOp| matches!(o, CigarOp::H(_));
    let lead: usize = ops.iter().take_while(h).map(|o| o.len()).sum();
    let trail: usize = if ops.iter().all(|o| matches!(o, CigarOp::H(_))) {
        0
    } else {
        ops.iter().rev().take_while(h).map(|o| o.len()).sum()
    };
    (lead, trail)
}

/// Turn the leading and/or trailing hard clips into soft clips, merging with
/// any soft clip beside them.
fn soften(ops: &[CigarOp], front: bool, back: bool) -> Vec<CigarOp> {
    let n = ops.len();
    let lead = ops.iter().take_while(|o| matches!(o, CigarOp::H(_))).count();
    let trail = ops.iter().rev().take_while(|o| matches!(o, CigarOp::H(_))).count().min(n - lead);
    let mut out: Vec<CigarOp> = Vec::with_capacity(n);
    for (i, op) in ops.iter().enumerate() {
        let op = match *op {
            CigarOp::H(k) if (front && i < lead) || (back && i >= n - trail) => CigarOp::S(k),
            other => other,
        };
        match (out.last_mut(), op) {
            (Some(CigarOp::S(a)), CigarOp::S(b)) => *a += b,
            (_, op) => out.push(op),
        }
    }
    out
}

/// Rebuild SA on every surviving alignment of a read whose alignments moved.
///
/// SA lists the read's *other* alignments as `rname,pos,strand,CIGAR,mapQ,NM;`,
/// primary first. CIGARs are written with soft clips, as BWA and minimap2 write
/// them. A read none of whose records carried SA is left without one; a record
/// whose read has no other alignment left loses it.
fn rebuild_sa(
    slots: &mut [Option<Record>],
    orig: &[Orig],
    origin: &[usize],
    plan_of: &[Option<&RecordPlan>],
    end_moved: [bool; 2],
    view: &HeaderView,
) -> Result<()> {
    for e in 0..2 {
        if !end_moved[e] {
            continue;
        }
        if !orig.iter().any(|o| o.end == e && o.had_sa) {
            continue;
        }
        let mut members: Vec<usize> = (0..slots.len())
            .filter(|&i| {
                slots[i].as_ref().is_some_and(|r| {
                    end_of(r) == e && !r.is_secondary() && !r.is_unmapped()
                })
            })
            .collect();
        members.sort_by_key(|&i| slots[i].as_ref().unwrap().is_supplementary());
        let entries: Vec<String> = members
            .iter()
            .map(|&i| {
                let r = slots[i].as_ref().unwrap();
                let o = origin[i];
                let ops = soften(&to_ops(r), true, true);
                format!(
                    "{},{},{},{},{},{};",
                    String::from_utf8_lossy(view.tid2name(r.tid() as u32)),
                    r.pos() + 1,
                    if r.is_reverse() { '-' } else { '+' },
                    cigar_string(&ops),
                    r.mapq(),
                    sa_nm(r, &orig[o], plan_of[o]),
                )
            })
            .collect();
        for (k, &i) in members.iter().enumerate() {
            let r = slots[i].as_mut().unwrap();
            let _ = r.remove_aux(b"SA");
            let sa: String = entries
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != k)
                .map(|(_, s)| s.as_str())
                .collect();
            if !sa.is_empty() {
                r.push_aux(b"SA", Aux::String(&sa))?;
            }
        }
    }
    Ok(())
}

/// The NM to give an alignment in SA.
///
/// The record's own NM when it has one. Otherwise -- `--stale-tags strip`
/// removed it from a clipped record -- it is worked out from the input MD over
/// the bases that survived clipping, which is exact. Failing that, the input
/// NM, which clipping can only have lowered; failing everything, 0.
fn sa_nm(rec: &Record, o: &Orig, p: Option<&RecordPlan>) -> i64 {
    if let Some(nm) = aux_int(rec, b"NM") {
        return nm;
    }
    let (front, back) = p.map_or((0, 0), |p| (p.clip_front, p.clip_back));
    o.md
        .as_deref()
        .and_then(|md| nm_from_md(&o.ops, md, front, back))
        .or(o.nm)
        .unwrap_or(0)
}

/// Edit distance of an alignment after `front` and `back` aligned (`M`/`I`)
/// bases are clipped, from its unclipped CIGAR and MD.
///
/// MD says which `M` bases mismatch and which reference bases are deleted, so
/// the edits that survive can be counted without the reference: mismatches and
/// inserted bases that are kept, and deletions with a kept base on each side
/// (a deletion left at an edge is removed by `reclip`). `None` if MD does not
/// fit the CIGAR.
fn nm_from_md(ops: &[CigarOp], md: &str, front: usize, back: usize) -> Option<i64> {
    enum Tok {
        Run(usize),
        Mis,
        Del(usize),
    }
    let mut toks: Vec<Tok> = Vec::new();
    let b = md.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i].is_ascii_digit() {
            let s = i;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            let n: usize = md[s..i].parse().ok()?;
            if n > 0 {
                toks.push(Tok::Run(n));
            }
        } else if b[i] == b'^' {
            let s = i + 1;
            i = s;
            while i < b.len() && b[i].is_ascii_alphabetic() {
                i += 1;
            }
            toks.push(Tok::Del(i - s));
        } else if b[i].is_ascii_alphabetic() {
            toks.push(Tok::Mis);
            i += 1;
        } else {
            return None;
        }
    }

    let total: usize = ops.iter().filter(|o| o.is_aligned()).map(|o| o.len()).sum();
    let kept = |a: usize| a >= front && a + back < total;
    let (mut nm, mut a) = (0i64, 0usize);
    let mut toks = toks.into_iter();
    let mut run = 0usize;
    for op in ops {
        match *op {
            CigarOp::M(k) => {
                for _ in 0..k {
                    if run == 0 {
                        match toks.next()? {
                            Tok::Run(r) => run = r,
                            Tok::Mis => {
                                nm += kept(a) as i64;
                                a += 1;
                                continue;
                            }
                            Tok::Del(_) => return None,
                        }
                    }
                    run -= 1;
                    a += 1;
                }
            }
            CigarOp::I(k) => {
                for _ in 0..k {
                    nm += kept(a) as i64;
                    a += 1;
                }
            }
            CigarOp::D(k) => {
                if run != 0 {
                    return None;
                }
                match toks.next()? {
                    Tok::Del(d) if d == k as usize => {}
                    _ => return None,
                }
                if a > 0 && kept(a - 1) && a < total && kept(a) {
                    nm += k as i64;
                }
            }
            CigarOp::N(_) | CigarOp::S(_) | CigarOp::H(_) | CigarOp::P(_) => {}
        }
    }
    Some(nm)
}

fn aux_int(rec: &Record, tag: &[u8]) -> Option<i64> {
    match rec.aux(tag).ok()? {
        Aux::I8(v) => Some(v as i64),
        Aux::U8(v) => Some(v as i64),
        Aux::I16(v) => Some(v as i64),
        Aux::U16(v) => Some(v as i64),
        Aux::I32(v) => Some(v as i64),
        Aux::U32(v) => Some(v as i64),
        _ => None,
    }
}

fn cigar_string(ops: &[CigarOp]) -> String {
    ops.iter()
        .map(|o| {
            let c = match o {
                CigarOp::M(_) => 'M',
                CigarOp::I(_) => 'I',
                CigarOp::D(_) => 'D',
                CigarOp::N(_) => 'N',
                CigarOp::S(_) => 'S',
                CigarOp::H(_) => 'H',
                CigarOp::P(_) => 'P',
            };
            format!("{}{}", o.len(), c)
        })
        .collect()
}

/// One past the last reference base, as htslib's `bam_endpos` counts it: an
/// alignment consuming no reference is one base long.
fn end_pos(rec: &Record) -> i64 {
    let rlen: i64 = if rec.is_unmapped() {
        0
    } else {
        rec.cigar()
            .iter()
            .map(|c| match *c {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) | Cigar::RefSkip(n) => n as i64,
                _ => 0,
            })
            .sum()
    };
    rec.pos() + rlen.max(1)
}

/// The primary line's facts its mate's fields are built from.
struct MateInfo {
    tid: i32,
    pos: i64,
    reverse: bool,
    unmapped: bool,
    mapq: u8,
    cigar: String,
    /// 5' position as `samtools fixmate` takes it: POS forward, end reverse.
    five: i64,
}

impl MateInfo {
    fn of(r: &Record) -> Self {
        let cigar = if r.is_unmapped() || r.cigar_len() == 0 {
            "*".to_string()
        } else {
            r.cigar().iter().map(|c| format!("{c}")).collect()
        };
        MateInfo {
            tid: r.tid(),
            pos: r.pos(),
            reverse: r.is_reverse(),
            unmapped: r.is_unmapped(),
            mapq: r.mapq(),
            cigar,
            five: if r.is_reverse() { end_pos(r) } else { r.pos() },
        }
    }
}

/// Make every record's mate fields agree with the mate's final primary line.
///
/// The rules are `samtools fixmate`'s, so running it on this output changes
/// none of these fields:
///
/// - an unmapped primary takes its mapped mate's RNAME and POS (`*`/0 when
///   both are unmapped, as SAM recommends);
/// - RNEXT/PNEXT are the mate primary's RNAME/POS, 0x20 its 0x10, 0x8 its 0x4;
/// - TLEN is the distance between the two 5' ends (POS for a forward record,
///   the end for a reverse one), signed towards the mate, and 0 when either is
///   unmapped or they are on different references;
/// - MC is the mate primary's CIGAR (`*` if unmapped), and is removed when both
///   are unmapped; an existing MQ becomes the mate's MAPQ, or is removed when
///   the mate is unmapped;
/// - 0x2 is cleared on the whole template unless both primaries are mapped to
///   one reference in forward-reverse orientation. It is never set.
///
/// Supplementary and secondary lines get the same mate fields relative to the
/// mate's primary, as BWA writes them (`fixmate` itself leaves them alone).
/// Nothing is done unless each end has exactly one primary line and both are
/// flagged paired.
fn sync_mates(slots: &mut [Option<Record>]) -> Result<()> {
    let primaries = |e: usize, slots: &[Option<Record>]| -> Vec<usize> {
        (0..slots.len())
            .filter(|&i| {
                slots[i].as_ref().is_some_and(|r| {
                    end_of(r) == e && !r.is_secondary() && !r.is_supplementary()
                })
            })
            .collect()
    };
    let (p1, p2) = (primaries(0, slots), primaries(1, slots));
    if p1.len() != 1 || p2.len() != 1 {
        return Ok(());
    }
    let (a, b) = (p1[0], p2[0]);
    if !slots[a].as_ref().unwrap().is_paired() || !slots[b].as_ref().unwrap().is_paired() {
        return Ok(());
    }

    let (ua, ub) = (slots[a].as_ref().unwrap().is_unmapped(), slots[b].as_ref().unwrap().is_unmapped());
    let place = |slots: &mut [Option<Record>], to: usize, tid: i32, pos: i64| {
        let r = slots[to].as_mut().unwrap();
        r.set_tid(tid);
        r.set_pos(pos);
    };
    match (ua, ub) {
        (true, false) => {
            let r = slots[b].as_ref().unwrap();
            let (t, p) = (r.tid(), r.pos());
            place(slots, a, t, p);
        }
        (false, true) => {
            let r = slots[a].as_ref().unwrap();
            let (t, p) = (r.tid(), r.pos());
            place(slots, b, t, p);
        }
        (true, true) => {
            place(slots, a, -1, -1);
            place(slots, b, -1, -1);
        }
        (false, false) => {}
    }

    let mates = [MateInfo::of(slots[b].as_ref().unwrap()), MateInfo::of(slots[a].as_ref().unwrap())];
    let (m1, m2) = (&mates[1], &mates[0]);
    let proper = !m1.unmapped && !m2.unmapped && m1.tid == m2.tid && {
        let (first, second) = if m1.five > m2.five { (m2, m1) } else { (m1, m2) };
        !first.reverse && second.reverse
    };

    for rec in slots.iter_mut().flatten() {
        let mate = &mates[end_of(rec)];
        rec.set_mtid(mate.tid);
        rec.set_mpos(mate.pos);
        if mate.reverse {
            rec.set_mate_reverse();
        } else {
            rec.unset_mate_reverse();
        }
        if mate.unmapped {
            rec.set_mate_unmapped();
        } else {
            rec.unset_mate_unmapped();
        }
        let unmapped = rec.is_unmapped();
        let _ = rec.remove_aux(b"MC");
        if !mate.unmapped || !unmapped {
            rec.push_aux(b"MC", Aux::String(&mate.cigar))?;
        }
        if rec.aux(b"MQ").is_ok() {
            let _ = rec.remove_aux(b"MQ");
            if !mate.unmapped {
                rec.push_aux(b"MQ", Aux::I32(mate.mapq as i32))?;
            }
        }
        let tlen = if !unmapped && !mate.unmapped && rec.tid() == mate.tid {
            let five = if rec.is_reverse() { end_pos(rec) } else { rec.pos() };
            mate.five - five
        } else {
            0
        };
        rec.set_insert_size(tlen);
        if !proper {
            rec.unset_proper_pair();
        }
    }
    Ok(())
}

/// The tags this step writes about the template it just resolved.
///
/// Deliberately outside the `X` namespace: `XM`, `XR` and `XG` are Bismark's
/// methylation call string and its companions, and this can run upstream of a
/// methylation caller. Writing a mismatch count into `XM` would not be a
/// collision so much as a silent replacement of the field such a pipeline
/// exists to read.
fn tag(
    rec: &mut Record,
    verdict: &Verdict,
    edits: Option<usize>,
    args: &OverlapArgs,
) -> Result<()> {
    if args.no_tag {
        return Ok(());
    }
    for c in ['O', 'N', 'S', 'X', 'K', 'V', 'v'] {
        let _ = rec.remove_aux(&args.tag(c));
    }
    let _ = rec.remove_aux(&args.length_tag());
    match verdict {
        Verdict::NoOverlap { .. } => {
            rec.push_aux(&args.tag('O'), Aux::I32(0))?;
        }
        Verdict::Unresolved { reason, .. } => {
            rec.push_aux(&args.tag('V'), Aux::String(reason.as_str()))?;
        }
        Verdict::Resolved(r) => {
            rec.push_aux(&args.length_tag(), Aux::I32(r.l as i32))?;
            rec.push_aux(&args.tag('O'), Aux::I32(r.overlap as i32))?;
            rec.push_aux(&args.tag('N'), Aux::I32(r.support as i32))?;
            rec.push_aux(&args.tag('S'), Aux::Float(r.span_frac))?;
            rec.push_aux(&args.tag('X'), Aux::I32(r.mismatches as i32))?;
            rec.push_aux(&args.tag('K'), Aux::Char(if r.kept_r1 { b'1' } else { b'2' }))?;
        }
    }
    // An upstream step that silently edits quality strings is unauditable.
    // When a concordance check disagrees with another caller by a fraction of a
    // percent, this is how you find out whether this step caused it.
    if let Some(n) = edits {
        rec.push_aux(&args.tag('v'), Aux::I32(n as i32))?;
    }
    Ok(())
}

fn to_htslib(o: &CigarOp) -> Cigar {
    match *o {
        CigarOp::M(n) => Cigar::Match(n),
        CigarOp::I(n) => Cigar::Ins(n),
        CigarOp::D(n) => Cigar::Del(n),
        CigarOp::N(n) => Cigar::RefSkip(n),
        CigarOp::S(n) => Cigar::SoftClip(n),
        CigarOp::H(n) => Cigar::HardClip(n),
        CigarOp::P(n) => Cigar::Pad(n),
    }
}

fn from_htslib(c: &Cigar) -> CigarOp {
    match *c {
        Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => CigarOp::M(n),
        Cigar::Ins(n) => CigarOp::I(n),
        Cigar::Del(n) => CigarOp::D(n),
        Cigar::RefSkip(n) => CigarOp::N(n),
        Cigar::SoftClip(n) => CigarOp::S(n),
        Cigar::HardClip(n) => CigarOp::H(n),
        Cigar::Pad(n) => CigarOp::P(n),
    }
}

fn to_ops(rec: &Record) -> Vec<CigarOp> {
    rec.cigar().iter().map(from_htslib).collect()
}

fn strip_tags(rec: &mut Record) {
    let _ = rec.remove_aux(b"AS");
    let _ = rec.remove_aux(b"SA");
    let _ = rec.remove_aux(b"XA");
    let _ = rec.remove_aux(b"NM");
    let _ = rec.remove_aux(b"MD");
}

/// NM, MD and AS describe the alignment before it was clipped.
///
/// NM is the edit distance and MD encodes where the mismatches are; both are
/// defined against the reference and can be recomputed from it. AS is the
/// aligner's own score, with its gap penalties and heuristics, and cannot be
/// honestly regenerated -- so it is always removed rather than adjusted.
/// SA lists the other alignments of this read by CIGAR, which clipping has just
/// invalidated. Leaving any of them stale is worse than removing them, because
/// downstream trusts them silently.
fn fix_tags(
    rec: &mut Record,
    args: &OverlapArgs,
    refr: Option<&Aref>,
    contigs: Option<&ContigMap>,
) -> Result<()> {
    let _ = rec.remove_aux(b"AS");
    let _ = rec.remove_aux(b"SA");
    let _ = rec.remove_aux(b"XA");

    match (args.stale_tags, refr) {
        (StaleTags::Recompute, Some(r)) => {
            // By name, via the map built once in `run` -- never by BAM tid,
            // which indexes a different list. A record on a contig the
            // reference lacks keeps no NM/MD rather than getting wrong ones.
            let Some(tid) = contigs.and_then(|m| m.refr_tid(rec.tid())) else {
                let _ = rec.remove_aux(b"NM");
                let _ = rec.remove_aux(b"MD");
                return Ok(());
            };
            let (nm, md) = nm_and_md(rec, r, tid);
            let _ = rec.remove_aux(b"NM");
            let _ = rec.remove_aux(b"MD");
            rec.push_aux(b"NM", Aux::I32(nm))?;
            rec.push_aux(b"MD", Aux::String(&md))?;
        }
        _ => {
            let _ = rec.remove_aux(b"NM");
            let _ = rec.remove_aux(b"MD");
        }
    }
    Ok(())
}

/// Edit distance and the MD string for the surviving alignment.
///
/// MD is runs of matching bases as counts, mismatched *reference* bases as
/// letters, and deletions as `^` followed by the deleted reference bases. It
/// exists so a caller can reconstruct the reference without loading it.
///
/// A `CIGAR N` appears in neither: it is not an edit and its bases were never
/// in the fragment. A reader walks MD against the CIGAR, which already says
/// where the skip is, so the match run runs straight through a junction.
fn nm_and_md(rec: &Record, refr: &Aref, tid: usize) -> (i32, String) {
    let s = rec.seq();
    let mut nm = 0i32;
    let mut md = String::new();
    let mut run = 0usize;
    let mut r = rec.pos();
    let mut q = 0usize;

    for op in rec.cigar().iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for _ in 0..n {
                    let rb = refr.slice(tid, r, r + 1).and_then(|x| x.first().copied());
                    let qb = Seq(s.encoded_base(q));
                    match rb {
                        Some(rb) if rb == qb => run += 1,
                        Some(rb) => {
                            md.push_str(&run.to_string());
                            run = 0;
                            md.push(rb.to_char());
                            nm += 1;
                        }
                        None => run += 1,
                    }
                    r += 1;
                    q += 1;
                }
            }
            Cigar::Ins(n) => {
                nm += n as i32;
                q += n as usize;
            }
            Cigar::Del(n) => {
                md.push_str(&run.to_string());
                run = 0;
                md.push('^');
                if let Some(bases) = refr.slice(tid, r, r + n as i64) {
                    for b in bases {
                        md.push(b.to_char());
                    }
                }
                nm += n as i32;
                r += n as i64;
            }
            // An intron is not a deletion and does not belong in MD. MD spells
            // out deleted reference bases after a `^`, and a reader walks it
            // alongside the CIGAR -- which already says where the skip is -- so
            // the match run continues across the junction rather than breaking
            // on it. Nor does an intron count toward NM: nothing was edited,
            // the transcript never contained those bases.
            Cigar::RefSkip(n) => {
                r += n as i64;
            }
            Cigar::SoftClip(n) => q += n as usize,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    md.push_str(&run.to_string());
    (nm, md)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::{built, LAST_IN_TEMPLATE, REVERSE};
    use clap::Parser;
    use rust_htslib::bam::record::Cigar;


    /// `OverlapArgs` is a clap `Args`, so the only honest way to make one is to
    /// parse a command line. Building it field by field would let a test run
    /// against a combination the CLI would reject.
    #[derive(Parser, Debug)]
    struct Harness {
        #[command(flatten)]
        args: OverlapArgs,
    }

    fn args(extra: &[&str]) -> OverlapArgs {
        let mut argv = vec!["overlap"];
        argv.extend_from_slice(extra);
        argv.extend(["in.bam", "out.bam"]);
        Harness::try_parse_from(argv).expect("parse").args
    }

    fn cigar_text(rec: &Record) -> String {
        rec.cigar().iter().map(|c| format!("{c}")).collect()
    }

    /// The two CIGAR conversion tables are written out by hand in opposite
    /// directions, so nothing but a round trip stops them drifting apart.
    #[test]
    fn cigar_ops_survive_a_round_trip() {
        let all = [
            CigarOp::M(3),
            CigarOp::I(1),
            CigarOp::D(2),
            CigarOp::N(7),
            CigarOp::S(4),
            CigarOp::H(5),
            CigarOp::P(6),
        ];
        for op in all {
            assert_eq!(from_htslib(&to_htslib(&op)), op, "{op:?} did not survive");
        }
        // `=` and `X` fold into `M` on the way in, which is deliberate: the
        // planner only cares whether a base consumes reference, query or both.
        // Note this is the one place D and N stay interchangeable at all: they
        // both consume reference and nothing else, so the planner's arms that
        // only advance a coordinate keep them together on purpose. The places
        // where the two differ -- the MD string here, and the emitted columns
        // in `alignment` -- say so explicitly.
        assert_eq!(from_htslib(&Cigar::Equal(3)), CigarOp::M(3));
        assert_eq!(from_htslib(&Cigar::Diff(3)), CigarOp::M(3));
    }

    /// The planner speaks about fragments, so R1 is derived from the template
    /// flag rather than from the strand.
    #[test]
    fn to_aln_reads_the_record_into_the_planners_shape() {
        let cig = [Cigar::SoftClip(2), Cigar::Match(4), Cigar::Ins(1), Cigar::Match(3)];
        let rec = built(b"q", b"ACGTACGTAC", &cig, 10, REVERSE);
        let aln = to_aln(7, &rec);

        assert_eq!(aln.id, 7);
        assert!(aln.is_r1, "no LAST_IN_TEMPLATE bit means R1");
        assert!(aln.is_reverse);
        assert!(!aln.is_supplementary);
        assert_eq!(aln.pos, 10);
        assert_eq!(aln.tid, 0);
        assert_eq!(
            aln.cigar,
            vec![CigarOp::S(2), CigarOp::M(4), CigarOp::I(1), CigarOp::M(3)]
        );
        // SEQ is carried across as raw nibbles, one per base, in record order.
        assert_eq!(aln.seq.len(), 10);
        assert_eq!(aln.seq[0], Seq::A);
        assert_eq!(aln.seq[1], Seq::C);
        assert_eq!(aln.qual.len(), 10);

        let r2 = built(b"q", b"ACGT", &[Cigar::Match(4)], 0, LAST_IN_TEMPLATE);
        assert!(!to_aln(0, &r2).is_r1);
    }

    /// Base and quality edits are absolute SEQ-space assignments, so applying a
    /// plan twice lands in the same place as applying it once.
    #[test]
    fn base_and_qual_edits_are_idempotent() {
        let a = args(&[]);
        let plan = RecordPlan {
            id: 0,
            qual_set: vec![(1, 0), (3, 12)],
            base_set: vec![(2, Seq::N)],
            ..Default::default()
        };

        let mut rec = built(b"q", b"ACGTACGT", &[Cigar::Match(8)], 0, 0);
        apply(&mut rec, &plan, &a, None, None).unwrap();
        let once = (rec.seq().as_bytes(), rec.qual().to_vec());
        assert_eq!(once.0[2], b'N', "base 2 set to N");
        assert_eq!(once.1[1], 0);
        assert_eq!(once.1[3], 12);
        assert_eq!(once.1[0], 37, "untouched positions keep their quality");

        apply(&mut rec, &plan, &a, None, None).unwrap();
        assert_eq!((rec.seq().as_bytes(), rec.qual().to_vec()), once, "not idempotent");
    }

    /// Out-of-range indices are ignored rather than panicking: a hard clip can
    /// shorten SEQ after the planner chose an index against the original.
    #[test]
    fn edits_past_the_end_of_seq_are_ignored() {
        let a = args(&[]);
        let plan = RecordPlan {
            id: 0,
            qual_set: vec![(99, 0)],
            base_set: vec![(99, Seq::N)],
            ..Default::default()
        };
        let mut rec = built(b"q", b"ACGT", &[Cigar::Match(4)], 0, 0);
        apply(&mut rec, &plan, &a, None, None).unwrap();
        assert_eq!(rec.seq().as_bytes(), b"ACGT");
    }

    /// A primary with nothing left is written as SAM describes an unmapped read:
    /// 0x4 set, 0x2 and 0x10 cleared, MAPQ 0, no CIGAR, and SEQ/QUAL back in the
    /// orientation the read was sequenced in.
    #[test]
    fn a_dropped_primary_is_written_as_an_unmapped_read() {
        let a = args(&[]);
        let plan = RecordPlan { id: 0, drop: true, ..Default::default() };
        let mut rec = built(b"q", b"AACGTC", &[Cigar::Match(6)], 5, REVERSE | 0x1 | 0x2);
        rec.set(b"q", Some(&CigarString(vec![Cigar::Match(6)])), b"AACGTC", &[10, 11, 12, 13, 14, 15]);
        rec.push_aux(b"NM", Aux::I32(3)).unwrap();

        apply(&mut rec, &plan, &a, None, None).unwrap();
        assert!(rec.is_unmapped());
        assert!(!rec.is_reverse(), "stored as sequenced");
        assert!(!rec.is_proper_pair());
        assert_eq!(rec.mapq(), 0);
        assert_eq!(rec.cigar_len(), 0);
        assert_eq!(rec.seq().as_bytes(), b"GACGTT", "reverse-complemented, not lost");
        assert_eq!(rec.qual(), &[15, 14, 13, 12, 11, 10]);
        assert!(rec.aux(b"NM").is_err(), "alignment tags go with the alignment");
    }

    /// Clipping the front moves POS; clipping the back does not.
    #[test]
    fn clipping_the_front_moves_the_alignment_start() {
        let a = args(&[]);

        let mut front = built(b"q", b"ACGTACGT", &[Cigar::Match(8)], 20, 0);
        apply(
            &mut front,
            &RecordPlan { id: 0, clip_front: 3, ..Default::default() },
            &a, None, None,
        ).unwrap();
        assert_eq!(front.pos(), 23, "three aligned bases removed from the front");
        assert_eq!(cigar_text(&front), "3S5M");

        let mut back = built(b"q", b"ACGTACGT", &[Cigar::Match(8)], 20, 0);
        apply(
            &mut back,
            &RecordPlan { id: 0, clip_back: 3, ..Default::default() },
            &a, None, None,
        ).unwrap();
        assert_eq!(back.pos(), 20, "clipping the back leaves POS alone");
        assert_eq!(cigar_text(&back), "5M3S");
    }

    /// A hard clip removes the bases from SEQ too, or the record stops
    /// describing itself.
    #[test]
    fn a_hard_clip_takes_seq_and_qual_with_it() {
        let a = args(&[]);
        let mut rec = built(b"q", b"ACGTACGT", &[Cigar::Match(8)], 20, 0);
        apply(
            &mut rec,
            &RecordPlan { id: 0, clip_front: 3, hard_front: true, ..Default::default() },
            &a, None, None,
        ).unwrap();
        assert_eq!(rec.seq().as_bytes(), b"TACGT", "the clipped bases are gone");
        assert_eq!(rec.qual().len(), 5, "QUAL stays parallel to SEQ");
        assert_eq!(cigar_text(&rec), "3H5M");
        assert_eq!(rec.pos(), 23);
    }

    /// `--stale-tags strip` removes what clipping invalidated rather than
    /// recomputing it, and AS is removed either way because it cannot be.
    #[test]
    fn clipping_removes_tags_it_invalidated() {
        let a = args(&["--stale-tags", "strip"]);
        let mut rec = built(b"q", b"ACGTACGT", &[Cigar::Match(8)], 20, 0);
        rec.push_aux(b"NM", Aux::I32(2)).unwrap();
        rec.push_aux(b"MD", Aux::String("8")).unwrap();
        rec.push_aux(b"AS", Aux::I32(99)).unwrap();

        apply(
            &mut rec,
            &RecordPlan { id: 0, clip_back: 2, ..Default::default() },
            &a, None, None,
        ).unwrap();
        for t in [b"NM".as_slice(), b"MD".as_slice(), b"AS".as_slice()] {
            assert!(rec.aux(t).is_err(), "{} survived clipping", String::from_utf8_lossy(t));
        }
    }

    /// An empty plan leaves the record byte for byte as it was.
    #[test]
    fn an_empty_plan_changes_nothing() {
        let a = args(&[]);
        let cig = [Cigar::SoftClip(1), Cigar::Match(6), Cigar::SoftClip(1)];
        let mut rec = built(b"q", b"ACGTACGT", &cig, 20, 0);
        let before = (rec.seq().as_bytes(), rec.qual().to_vec(), rec.pos(), cigar_text(&rec));

        let plan = RecordPlan { id: 0, ..Default::default() };
        assert!(plan.is_empty());
        apply(&mut rec, &plan, &a, None, None).unwrap();
        assert_eq!(
            (rec.seq().as_bytes(), rec.qual().to_vec(), rec.pos(), cigar_text(&rec)),
            before
        );
    }

    // ---- whole templates ------------------------------------------------

    /// A deterministic reference, so reads cut from it agree with each other.
    fn genome(seed: u64, len: usize) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        let mut v = s.to_vec();
        reverse_complement(&mut v);
        v
    }

    /// A paired record: `ops` as text-free htslib CIGAR, placed on `tid`.
    #[allow(clippy::too_many_arguments)]
    fn rec(
        name: &[u8],
        tid: i32,
        pos: i64,
        mapq: u8,
        flags: u16,
        cigar: &[Cigar],
        seq: &[u8],
        mate: (i32, i64),
    ) -> Record {
        let mut r = built(name, seq, cigar, pos, flags | 0x1);
        r.set_tid(tid);
        r.set_mapq(mapq);
        r.set_mtid(mate.0);
        r.set_mpos(mate.1);
        r
    }

    fn resolve_with(extra: &[&str], recs: Vec<Record>) -> Vec<Record> {
        let a = args(extra);
        let opts = a.options();
        let view = crate::test_support::header_many(&[("chr1", 1000), ("chr2", 1000)]);
        let ctx = Ctx { opts: &opts, args: &a, refr: None, contigs: None, view: &view };
        let mut stats = overlap::Stats::default();
        resolve(recs, &ctx, &mut stats).unwrap()
    }

    fn string_tag(r: &Record, t: &[u8]) -> Option<String> {
        match r.aux(t) {
            Ok(Aux::String(s)) => Some(s.to_string()),
            _ => None,
        }
    }

    /// Clipping R2's 3' end moves its POS and changes its CIGAR, so R1's PNEXT
    /// and MC must follow. Neither 5' end moved, so TLEN (5' to 5', as
    /// `samtools fixmate` measures it) does not.
    #[test]
    fn mate_fields_follow_a_clipped_mate() {
        let g = genome(1, 200);
        let r1 = rec(b"p", 0, 40, 60, 0x40 | 0x20 | 0x2, &[Cigar::Match(40)], &g[40..80], (0, 60));
        let mut r2 = rec(b"p", 0, 60, 60, 0x80 | REVERSE | 0x2, &[Cigar::Match(40)], &g[60..100], (0, 40));
        r2.push_aux(b"MC", Aux::String("40M")).unwrap();
        let mut r1 = r1;
        r1.push_aux(b"MC", Aux::String("40M")).unwrap();
        r1.push_aux(b"MQ", Aux::I32(60)).unwrap();

        let out = resolve_with(&[], vec![r1, r2]);
        let (r1, r2) = (&out[0], &out[1]);
        assert_eq!(cigar_text(r2), "20S20M");
        assert_eq!(r2.pos(), 80);
        assert_eq!((r1.mtid(), r1.mpos()), (0, 80), "PNEXT is the mate's new POS");
        assert_eq!(string_tag(r1, b"MC").as_deref(), Some("20S20M"));
        assert_eq!((r1.insert_size(), r2.insert_size()), (60, -60));
        assert!(r1.is_proper_pair() && r2.is_proper_pair(), "still forward-reverse");
        assert!(r1.is_mate_reverse() && !r2.is_mate_reverse());
        assert_eq!(aux_int(r1, b"MQ"), Some(60));
    }

    /// Read-through: R1 is kept and R2 has nothing aligned left. R2 becomes an
    /// unmapped read placed at R1, and R1 learns that its mate is unmapped.
    #[test]
    fn an_unmapped_mate_is_placed_beside_its_mate_and_flagged_on_it() {
        let g = genome(2, 300);
        let adapter = b"AGATCGGAAG";
        let r1_seq: Vec<u8> = [&g[160..190], &adapter[..]].concat();
        let r2_read: Vec<u8> = [revcomp(&g[160..190]), adapter.to_vec()].concat();
        let r2_seq = revcomp(&r2_read);
        let r1 = rec(b"t", 0, 160, 60, 0x40 | 0x20 | 0x2, &[Cigar::Match(30), Cigar::SoftClip(10)], &r1_seq, (0, 160));
        let r2 = rec(b"t", 0, 160, 60, 0x80 | REVERSE | 0x2, &[Cigar::SoftClip(10), Cigar::Match(30)], &r2_seq, (0, 160));

        let out = resolve_with(&[], vec![r1, r2]);
        assert_eq!(out.len(), 2, "both ends are still in the file");
        let (r1, r2) = (&out[0], &out[1]);
        assert!(r2.is_unmapped() && !r2.is_reverse() && !r2.is_proper_pair());
        assert_eq!((r2.tid(), r2.pos(), r2.mapq(), r2.cigar_len()), (0, 160, 0, 0));
        assert_eq!(r2.seq().as_bytes(), r2_read, "SEQ as sequenced");
        assert_eq!(string_tag(r2, b"MC").as_deref(), Some("30M10S"));

        assert!(r1.is_mate_unmapped() && !r1.is_mate_reverse() && !r1.is_proper_pair());
        assert_eq!((r1.mtid(), r1.mpos()), (0, 160));
        assert_eq!(string_tag(r1, b"MC").as_deref(), Some("*"));
        assert_eq!((r1.insert_size(), r2.insert_size()), (0, 0));
    }

    /// The aligner made R1's chr1 segment primary, and it lies wholly inside
    /// the overlap R2 keeps. Its chr2 supplementary survives. The read must
    /// keep exactly one primary line, carrying the whole read: the
    /// supplementary's alignment moves onto the primary record with its hard
    /// clip softened, and no supplementary or unmapped R1 record is left.
    #[test]
    fn promotion_leaves_one_primary_carrying_the_whole_read() {
        let (c1, c2) = (genome(3, 300), genome(4, 300));
        let read1: Vec<u8> = [&c2[140..160], &c1[40..100]].concat();
        let mut prim = rec(b"h", 0, 40, 20, 0x40, &[Cigar::SoftClip(20), Cigar::Match(60)], &read1, (0, 40));
        prim.push_aux(b"SA", Aux::String("chr2,141,+,20M60S,60,0;")).unwrap();
        prim.push_aux(b"RG", Aux::String("lib1")).unwrap();
        let mut sup = rec(b"h", 1, 140, 60, 0x40 | 0x800, &[Cigar::Match(20), Cigar::HardClip(60)], &c2[140..160], (0, 40));
        sup.push_aux(b"SA", Aux::String("chr1,41,+,20S60M,20,0;")).unwrap();
        sup.push_aux(b"NM", Aux::I32(0)).unwrap();
        let r2 = rec(b"h", 0, 40, 60, 0x80 | REVERSE, &[Cigar::Match(80)], &c1[40..120], (0, 40));

        let out = resolve_with(&[], vec![sup, prim, r2]);
        let r1s: Vec<&Record> = out.iter().filter(|r| !r.is_last_in_template()).collect();
        assert_eq!(r1s.len(), 1, "one R1 line, not an unmapped primary plus a promoted one");
        let r1 = r1s[0];
        assert!(!r1.is_supplementary() && !r1.is_unmapped());
        assert_eq!((r1.tid(), r1.pos(), r1.mapq()), (1, 140, 60));
        assert_eq!(cigar_text(r1), "20M60S");
        assert_eq!(r1.seq().as_bytes(), read1, "the whole read, not the 20 supplementary bases");
        assert!(r1.aux(b"SA").is_err(), "no other alignment of R1 is left");
        assert_eq!(string_tag(r1, b"RG").as_deref(), Some("lib1"), "read-level tags stay");
        assert_eq!(aux_int(r1, b"NM"), Some(0), "NM describes the survivor's alignment");

        let r2 = out.iter().find(|r| r.is_last_in_template()).unwrap();
        assert_eq!((r2.mtid(), r2.mpos()), (1, 140));
        assert_eq!(string_tag(r2, b"MC").as_deref(), Some("20M60S"));
        assert!(!r2.is_proper_pair() && !r1.is_proper_pair(), "mates on different references");
        assert_eq!(r2.insert_size(), 0);
    }

    /// When R1's supplementary is clipped, the SA on R1's primary lists its
    /// new CIGAR and its NM over the bases that are left (worked out from MD,
    /// since `--stale-tags strip` removed the record's own NM).
    #[test]
    fn sa_lists_the_final_state_of_the_other_alignments() {
        let (c1, c2) = (genome(5, 300), genome(6, 300));
        let read1: Vec<u8> = [&c1[40..100], &c2[100..120]].concat();
        let mut prim = rec(b"s", 0, 40, 60, 0x40, &[Cigar::Match(60), Cigar::SoftClip(20)], &read1, (1, 105));
        prim.push_aux(b"NM", Aux::I32(0)).unwrap();
        prim.push_aux(b"SA", Aux::String("chr2,101,+,60S20M,40,1;")).unwrap();
        let mut sup = rec(b"s", 1, 100, 40, 0x40 | 0x800, &[Cigar::HardClip(60), Cigar::Match(20)], &c2[100..120], (1, 105));
        // One mismatch, at the supplementary's 18th base: inside the part R1
        // gives up, so the clipped alignment has none.
        sup.push_aux(b"NM", Aux::I32(1)).unwrap();
        sup.push_aux(b"MD", Aux::String("17A2")).unwrap();
        sup.push_aux(b"SA", Aux::String("chr1,41,+,60M20S,60,0;")).unwrap();
        let r2 = rec(b"s", 1, 105, 60, 0x80 | REVERSE, &[Cigar::Match(35)], &c2[105..140], (0, 40));

        let out = resolve_with(&[], vec![prim, sup, r2]);
        assert_eq!(cigar_text(&out[1]), "60H5M15S");
        assert_eq!(string_tag(&out[0], b"SA").as_deref(), Some("chr2,101,+,60S5M15S,40,0;"));
        assert_eq!(string_tag(&out[1], b"SA").as_deref(), Some("chr1,41,+,60M20S,60,0;"));
    }

    #[test]
    fn nm_after_clipping_is_counted_from_md() {
        use CigarOp::*;
        // aligned index: M 0-9, I 10-11, M 12-16, (3D), M 17-21; mismatch at 18.
        let ops = [M(10), I(2), M(5), D(3), M(5)];
        let md = "15^ACG1T3";
        assert_eq!(nm_from_md(&ops, md, 0, 0), Some(6));
        assert_eq!(nm_from_md(&ops, md, 0, 4), Some(5), "mismatch clipped");
        assert_eq!(nm_from_md(&ops, md, 0, 5), Some(2), "deletion left at the edge goes too");
        assert_eq!(nm_from_md(&ops, md, 11, 0), Some(5), "one inserted base clipped");
        assert_eq!(nm_from_md(&ops, "20", 0, 0), None, "MD that does not fit the CIGAR");
    }

    /// QUAL `*` has no qualities to combine. It stays `*`, and the mate that
    /// does have qualities is not edited against a missing value.
    #[test]
    fn a_missing_quality_string_stays_missing() {
        let g = genome(7, 200);
        let mut r1 = rec(b"q", 0, 40, 60, 0x40, &[Cigar::Match(40)], &g[40..80], (0, 60));
        let mut r2 = rec(b"q", 0, 60, 60, 0x80 | REVERSE, &[Cigar::Match(40)], &g[60..100], (0, 40));
        r2.set(b"q", Some(&CigarString(vec![Cigar::Match(40)])), &g[60..100], &[QUAL_MISSING; 40]);
        let out = resolve_with(&[], vec![r1.clone(), r2.clone()]);
        assert!(out[1].qual().iter().all(|&q| q == QUAL_MISSING), "R2 QUAL is still *");
        assert!(out[0].qual().iter().all(|&q| q == 37), "R1 not combined with a missing value");
        assert_eq!(cigar_text(&out[1]), "20S20M", "the overlap is still resolved");

        r1.set(b"q", Some(&CigarString(vec![Cigar::Match(40)])), &g[40..80], &[QUAL_MISSING; 40]);
        let out = resolve_with(&[], vec![r1, r2]);
        assert!(out.iter().all(|r| r.qual().iter().all(|&q| q == QUAL_MISSING)));
    }

    /// `--compression-level` reaches the writer: level 0 writes uncompressed
    /// blocks, so the same input comes out much larger than at level 9.
    #[test]
    fn the_compression_level_is_honoured() {
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::{Format, Writer};
        let dir = crate::test_support::temp_dir();
        let input = dir.join("alnbase_test_overlap_level_in.bam");
        let mut h = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"GO", "query");
        h.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        h.push_record(&sq);
        let g = genome(8, 1000);
        {
            let mut w = Writer::from_path(&input, &h, Format::Bam).unwrap();
            for t in 0..300 {
                let name = format!("t{t}");
                let at = (t * 7) % 800;
                w.write(&rec(name.as_bytes(), 0, at as i64, 60, 0x40, &[Cigar::Match(40)], &g[at..at + 40], (0, 0))).unwrap();
                w.write(&rec(name.as_bytes(), 0, at as i64 + 90, 60, 0x80 | REVERSE, &[Cigar::Match(40)], &g[at + 90..at + 130], (0, 0))).unwrap();
            }
        }
        let size = |level: &str| {
            let out = dir.join(format!("alnbase_test_overlap_level_{level}.bam"));
            let argv = ["overlap", "--compression-level", level, input.to_str().unwrap(), out.to_str().unwrap()];
            run(Harness::try_parse_from(argv).unwrap().args).unwrap();
            std::fs::metadata(&out).unwrap().len()
        };
        let (l0, l9) = (size("0"), size("9"));
        assert!(l0 > 2 * l9, "level 0 gave {l0} bytes and level 9 gave {l9}");
        assert!(args(&["--compression-level", "10"]).validate().is_err());
    }
}