use radsort::sort_by_key;

use crate::ruranges_structs::{GroupType, PositionType};

#[allow(clippy::too_many_arguments)]
/// A projected fragment, plus the annotation row that produced it.
///
/// Deliberately private: `StrandInterval` is public API, and widening it would
/// break any downstream struct literal or exhaustive destructuring.
#[derive(Debug, Clone)]
struct ProjectedFragment<T: PositionType> {
    start: T,
    end: T,
    idx: u32,
    fwd: bool,
    exon: u32,
}

/// Project local coordinates onto genomic ones, reporting only the query rows.
///
/// [`map_to_global_with_exons`] additionally reports which annotation row
/// produced each fragment, which callers need when they carry annotation
/// metadata onto the output.
#[allow(clippy::too_many_arguments)]
pub fn map_to_global<G: GroupType, T: PositionType>(
    /* left  table (exons / annotation) */
    ex_tx: &[G],
    ex_local_start: &[T],
    ex_local_end: &[T],

    /* right table (queries / local intervals) */
    q_tx: &[G],
    q_start: &[T],
    q_end: &[T],

    /* extra, still in strict order */
    ex_chr_code: &[G],
    ex_genome_start: &[T],
    ex_genome_end: &[T],
    ex_fwd: &[bool],
    q_fwd: &[bool],
    sort_output: bool,
) -> (Vec<u32>, Vec<T>, Vec<T>, Vec<bool>) {
    let (idx, start, end, fwd, _exon) = map_to_global_with_exons(
        ex_tx,
        ex_local_start,
        ex_local_end,
        q_tx,
        q_start,
        q_end,
        ex_chr_code,
        ex_genome_start,
        ex_genome_end,
        ex_fwd,
        q_fwd,
        sort_output,
    );
    (idx, start, end, fwd)
}

/// [`map_to_global`], additionally reporting the annotation row that produced
/// each fragment as a fifth output.
#[allow(clippy::too_many_arguments)]
pub fn map_to_global_with_exons<G: GroupType, T: PositionType>(
    /* left  table (exons / annotation) */
    ex_tx: &[G],
    ex_local_start: &[T],
    ex_local_end: &[T],

    /* right table (queries / local intervals) */
    q_tx: &[G],
    q_start: &[T],
    q_end: &[T],

    /* extra, still in strict order */
    ex_chr_code: &[G],
    ex_genome_start: &[T],
    ex_genome_end: &[T],
    ex_fwd: &[bool],
    q_fwd: &[bool],
    sort_output: bool,
) -> (Vec<u32>, Vec<T>, Vec<T>, Vec<bool>, Vec<u32>) {
    // ------------------- sanity checks (debug-only) ------------------------
    debug_assert_eq!(ex_tx.len(), ex_local_start.len());
    debug_assert_eq!(ex_tx.len(), ex_local_end.len());
    debug_assert_eq!(ex_tx.len(), ex_chr_code.len());
    debug_assert_eq!(ex_tx.len(), ex_genome_start.len());
    debug_assert_eq!(ex_tx.len(), ex_genome_end.len());
    debug_assert_eq!(ex_tx.len(), ex_fwd.len());

    debug_assert_eq!(q_tx.len(), q_start.len());
    debug_assert_eq!(q_tx.len(), q_end.len());
    debug_assert_eq!(q_tx.len(), q_fwd.len());

    // ------------------- output buffers -----------------------------------
    let mut results = Vec::new();

    // ------------------- two-pointer sweep ---------------------------------
    let mut ei = 0usize; // exon pointer
    let mut qi = 0usize; // query pointer
    let ex_n = ex_tx.len();
    let q_n = q_tx.len();

    while qi < q_n {
        let tx_code = q_tx[qi];

        // move exon pointer to this transcript (or beyond)
        while ei < ex_n && ex_tx[ei] < tx_code {
            ei += 1;
        }

        // if no exons for this transcript, skip its queries
        if ei >= ex_n || ex_tx[ei] != tx_code {
            while qi < q_n && q_tx[qi] == tx_code {
                qi += 1;
            }
            continue;
        }

        // ------------------------------------------------------------
        // process all queries with transcript == tx_code
        // ------------------------------------------------------------
        let mut ej = ei; // exon cursor inside tx

        while qi < q_n && q_tx[qi] == tx_code {
            let mut l = q_start[qi];
            let lend = q_end[qi];
            let idx = qi as u32; // row number into query table
            let local_f = q_fwd[qi];

            // advance exon cursor until its end is after l
            while ej < ex_n && ex_tx[ej] == tx_code && ex_local_end[ej] <= l {
                ej += 1;
            }

            let mut ek = ej;
            while l < lend && ek < ex_n && ex_tx[ek] == tx_code {
                let el_start = ex_local_start[ek];
                let el_end = ex_local_end[ek];

                if l >= el_end {
                    ek += 1;
                    continue;
                }

                // clip to current exon
                let seg_end_local = if lend < el_end { lend } else { el_end };

                // translate to genome
                let offset1 = l - el_start;
                let offset2 = seg_end_local - el_start;

                let (g_start, g_end) = if ex_fwd[ek] {
                    (ex_genome_start[ek] + offset1, ex_genome_start[ek] + offset2)
                } else {
                    (ex_genome_end[ek] - offset2, ex_genome_end[ek] - offset1)
                };

                // push result
                results.push(ProjectedFragment {
                    start: g_start,
                    end: g_end,
                    idx,
                    fwd: local_f == ex_fwd[ek],
                    exon: ek as u32,
                });

                // advance inside query
                l = seg_end_local;
                if l >= lend {
                    break;
                }
                ek += 1;
            }

            qi += 1; // next query row
        }

        // skip remaining exons of this transcript
        while ei < ex_n && ex_tx[ei] == tx_code {
            ei += 1;
        }
    }

    if sort_output {
        sort_by_key(&mut results, |i| i.idx);
    }

    let mut out_idxs = Vec::with_capacity(results.len());
    let mut out_starts = Vec::with_capacity(results.len());
    let mut out_ends = Vec::with_capacity(results.len());
    let mut out_strands = Vec::with_capacity(results.len());
    let mut out_exons = Vec::with_capacity(results.len());

    for rec in results {
        out_idxs.push(rec.idx);
        out_starts.push(rec.start);
        out_ends.push(rec.end);
        out_strands.push(rec.fwd);
        out_exons.push(rec.exon);
    }

    (out_idxs, out_starts, out_ends, out_strands, out_exons)
}

#[cfg(test)]
mod tests {
    use super::{map_to_global, map_to_global_with_exons};

    /// Two exons of one transcript; the query crosses their junction.
    fn junction_case() -> (
        Vec<u32>,
        Vec<i64>,
        Vec<i64>,
        Vec<u32>,
        Vec<i64>,
        Vec<i64>,
        Vec<u32>,
        Vec<i64>,
        Vec<i64>,
        Vec<bool>,
        Vec<bool>,
    ) {
        (
            vec![0, 0],
            vec![0, 100],
            vec![100, 200],
            vec![0],
            vec![95],
            vec![105],
            vec![0, 0],
            vec![1000, 3000],
            vec![1100, 3100],
            vec![true, true],
            vec![true],
        )
    }

    #[test]
    fn exon_output_names_the_row_that_produced_each_fragment() {
        let (ex_tx, ex_ls, ex_le, q_tx, q_s, q_e, ex_chr, ex_gs, ex_ge, ex_fwd, q_fwd) =
            junction_case();
        let (idx, starts, ends, fwd, exons) = map_to_global_with_exons(
            &ex_tx, &ex_ls, &ex_le, &q_tx, &q_s, &q_e, &ex_chr, &ex_gs, &ex_ge, &ex_fwd, &q_fwd,
            false,
        );
        assert_eq!(idx, vec![0, 0]);
        assert_eq!(starts, vec![1095, 3000]);
        assert_eq!(ends, vec![1100, 3005]);
        assert_eq!(fwd, vec![true, true]);
        // The first fragment comes from exon row 0, the second from row 1.
        assert_eq!(exons, vec![0, 1]);
    }

    #[test]
    fn the_four_output_form_is_unchanged() {
        let (ex_tx, ex_ls, ex_le, q_tx, q_s, q_e, ex_chr, ex_gs, ex_ge, ex_fwd, q_fwd) =
            junction_case();
        let short = map_to_global(
            &ex_tx, &ex_ls, &ex_le, &q_tx, &q_s, &q_e, &ex_chr, &ex_gs, &ex_ge, &ex_fwd, &q_fwd,
            false,
        );
        let (idx, starts, ends, fwd, _) = map_to_global_with_exons(
            &ex_tx, &ex_ls, &ex_le, &q_tx, &q_s, &q_e, &ex_chr, &ex_gs, &ex_ge, &ex_fwd, &q_fwd,
            false,
        );
        assert_eq!(short, (idx, starts, ends, fwd));
    }

    /// A reverse-strand exon projects the fragment onto the mirrored coordinates.
    #[test]
    fn reverse_strand_exons_report_their_own_row() {
        let (idx, starts, ends, fwd, exons) = map_to_global_with_exons(
            &[0_u32, 0],
            &[0_i64, 100],
            &[100_i64, 200],
            &[0_u32],
            &[95_i64],
            &[105_i64],
            &[0_u32, 0],
            &[1000_i64, 3000],
            &[1100_i64, 3100],
            &[false, false],
            &[true],
            false,
        );
        assert_eq!(idx, vec![0, 0]);
        assert_eq!(starts, vec![1000, 3095]);
        assert_eq!(ends, vec![1005, 3100]);
        assert_eq!(fwd, vec![false, false]);
        assert_eq!(exons, vec![0, 1]);
    }
}
