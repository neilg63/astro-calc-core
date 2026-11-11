//! Aspects between ecliptic longitudes, and the periods during which two charts are in
//! aspect.
//!
//! The crate defines no aspects of its own. A client supplies the angles and orbs it cares
//! about as [`AspectOrb`] values (or `(angle, orb)` pairs), so any harmonic can be used.
//!
//! Two levels are offered:
//! * instant tests on longitudes: [`separation`], [`match_aspects`], [`best_aspect`],
//!   [`aspects_between`];
//! * a search over time, [`Ephemeris::progressed_aspects`], which compares two people's
//!   natal or secondary-progressed charts and returns when each aspect entered its orb,
//!   became exact and left it.

use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::math::diff_deg;
use crate::position::{BodyPos, CoordSystem, Perspective};
use crate::spk::Ephemeris;
use crate::time::JdUt;

/// Mean tropical year in days, the default "day for a year" ratio of secondary progressions.
pub const TROPICAL_YEAR: f64 = 365.242_19;

/// Times found by the period search are refined to this many days (about 9 seconds).
const TIME_TOL: f64 = 1e-4;

/// An aspect to look for: a target separation and the allowed deviation either side.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AspectOrb {
    /// Target separation in degrees, 0–180.
    pub angle: f64,
    /// Orb in degrees: the aspect holds while the separation is within `angle ± orb`.
    pub orb: f64,
}

impl AspectOrb {
    /// Any angle is folded into 0–180 (so 240 becomes 120); the orb's sign is ignored.
    pub fn new(angle: f64, orb: f64) -> Self {
        Self { angle: diff_deg(0.0, angle).abs(), orb: orb.abs() }
    }
}

impl From<(f64, f64)> for AspectOrb {
    fn from((angle, orb): (f64, f64)) -> Self {
        Self::new(angle, orb)
    }
}

/// Converts `(angle, orb)` pairs into aspect definitions.
pub fn aspect_orbs(pairs: &[(f64, f64)]) -> Vec<AspectOrb> {
    pairs.iter().map(|&p| p.into()).collect()
}

/// Shortest angular distance between two longitudes, 0–180 degrees.
#[inline]
pub fn separation(lng_a: f64, lng_b: f64) -> f64 {
    diff_deg(lng_a, lng_b).abs()
}

/// One aspect found between two longitudes.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AspectMatch {
    /// The matched definition.
    pub aspect: AspectOrb,
    /// Its index in the list of definitions that was searched.
    pub index: usize,
    /// Actual separation, 0–180 degrees.
    pub separation: f64,
    /// `separation − aspect.angle`: negative when closer than the exact aspect.
    pub deviation: f64,
}

/// All aspects from `targets` that the two longitudes make, in the order of `targets`.
/// Several can match at once when orbs overlap. Does not allocate.
pub fn match_aspects(lng_a: f64, lng_b: f64, targets: &[AspectOrb]) -> impl Iterator<Item = AspectMatch> + '_ {
    let sep = separation(lng_a, lng_b);
    targets.iter().enumerate().filter_map(move |(index, &aspect)| {
        let deviation = sep - aspect.angle;
        (deviation.abs() <= aspect.orb).then_some(AspectMatch { aspect, index, separation: sep, deviation })
    })
}

/// The matching aspect with the smallest deviation, if any.
pub fn best_aspect(lng_a: f64, lng_b: f64, targets: &[AspectOrb]) -> Option<AspectMatch> {
    match_aspects(lng_a, lng_b, targets).min_by(|x, y| x.deviation.abs().total_cmp(&y.deviation.abs()))
}

/// An aspect between a body of one chart and a body of another at one moment.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BodyAspect {
    pub body_a: Body,
    pub body_b: Body,
    #[cfg_attr(feature = "serde", serde(flatten))]
    pub matched: AspectMatch,
}

/// Every aspect between each position in `a` and each in `b` (by `lng`). Pass the same slice
/// twice for the aspects within one chart; a body is then also compared with itself, so
/// leave out a zero-degree target or filter `body_a != body_b`.
pub fn aspects_between(a: &[BodyPos], b: &[BodyPos], targets: &[AspectOrb]) -> Vec<BodyAspect> {
    let mut out = Vec::new();
    for pa in a {
        for pb in b {
            out.extend(match_aspects(pa.lng, pb.lng, targets).map(|matched| BodyAspect { body_a: pa.body, body_b: pb.body, matched }));
        }
    }
    out
}

/// Which chart of a person takes part in a comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ChartKind {
    /// The birth chart: positions fixed at the moment of birth.
    Natal,
    /// Secondary progression: at age *Y* years the chart is the sky *Y* days after birth.
    #[default]
    Progressed,
}

/// A person in a chart comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChartSubject {
    pub birth: JdUt,
    /// Birthplace. It only affects topocentric results.
    pub geo: GeoPos,
    pub kind: ChartKind,
}

impl ChartSubject {
    pub fn natal(birth: impl Into<JdUt>, geo: GeoPos) -> Self {
        Self { birth: birth.into(), geo, kind: ChartKind::Natal }
    }

    pub fn progressed(birth: impl Into<JdUt>, geo: GeoPos) -> Self {
        Self { birth: birth.into(), geo, kind: ChartKind::Progressed }
    }

    /// The Julian Day whose sky represents this person on the real date `jd`:
    /// the birth moment for a natal chart, `birth + age / year_days` for a progressed one.
    pub fn chart_jd(&self, jd: f64, year_days: f64) -> f64 {
        match self.kind {
            ChartKind::Natal => self.birth.value(),
            ChartKind::Progressed => self.birth.value() + (jd - self.birth.value()) / year_days,
        }
    }
}

/// Settings for [`Ephemeris::progressed_aspects`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProgressionOptions {
    pub perspective: Perspective,
    /// Days of real life represented by one day of progression. Default [`TROPICAL_YEAR`].
    pub year_days: f64,
    /// Sampling step of the search in real days. Default 30. An aspect that stays in orb
    /// for less than a step can be missed if its closest approach is far from a sample;
    /// with progressions even the Moon moves only about 1° in 30 days.
    pub step_days: f64,
}

impl Default for ProgressionOptions {
    fn default() -> Self {
        Self { perspective: Perspective::Geocentric, year_days: TROPICAL_YEAR, step_days: 30.0 }
    }
}

/// A stretch of time during which a body of chart A is in aspect to a body of chart B.
/// Times are real (calendar) UT Julian Days, not progressed ones.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AspectPeriod {
    pub body_a: Body,
    pub body_b: Body,
    pub aspect: AspectOrb,
    /// When the orb was entered; `None` if already in orb at the start of the search.
    pub enter_jd: Option<f64>,
    /// When the aspect was exact. Empty if it never got there, several entries if a
    /// retrograde brought it back.
    pub exact_jds: Vec<f64>,
    /// When the orb was left; `None` if still in orb at the end of the search.
    pub leave_jd: Option<f64>,
    /// Smallest deviation from the exact aspect during the period, in degrees.
    pub min_orb: f64,
    /// When that closest approach occurred.
    pub min_orb_jd: f64,
}

type Eval<'a> = dyn Fn(f64) -> Result<f64, CalcError> + 'a;

/// Root of `f` between `a` and `b`, where `f(a)` (= `fa`) and `f(b)` differ in sign.
fn bisect(f: &Eval, mut a: f64, mut fa: f64, mut b: f64) -> Result<f64, CalcError> {
    while b - a > TIME_TOL {
        let m = 0.5 * (a + b);
        let fm = f(m)?;
        if (fm > 0.0) == (fa > 0.0) {
            a = m;
            fa = fm;
        } else {
            b = m;
        }
    }
    Ok(0.5 * (a + b))
}

/// Golden-section minimum of `f` on `[a, b]`: `(t, f(t))`.
fn golden_min(f: &Eval, mut a: f64, mut b: f64) -> Result<(f64, f64), CalcError> {
    const R: f64 = 0.618_033_988_749_894_9;
    let mut c = b - R * (b - a);
    let mut d = a + R * (b - a);
    let (mut fc, mut fd) = (f(c)?, f(d)?);
    while b - a > TIME_TOL {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - R * (b - a);
            fc = f(c)?;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + R * (b - a);
            fd = f(d)?;
        }
    }
    Ok(if fc < fd { (c, fc) } else { (d, fd) })
}

type OpenPeriod = (Option<f64>, Vec<(f64, f64)>);

struct RawPeriod {
    enter: Option<f64>,
    leave: Option<f64>,
    exact: Vec<f64>,
    min_orb: f64,
    min_jd: f64,
}

/// Exact moments: sign changes of the deviation between consecutive `(time, deviation)` nodes.
fn crossings(nodes: &[(f64, f64)], h: &Eval) -> Result<Vec<f64>, CalcError> {
    let mut exact = Vec::new();
    for w in nodes.windows(2) {
        let ((t0, h0), (t1, h1)) = (w[0], w[1]);
        if (h0 > 0.0) != (h1 > 0.0) && h0.abs() < 90.0 && h1.abs() < 90.0 {
            exact.push(bisect(h, t0, h0, t1)?);
        }
    }
    Ok(exact)
}

/// Completes a period from its nodes (boundaries and the samples inside, in time order).
fn finish(mut nodes: Vec<(f64, f64)>, enter: Option<f64>, leave: Option<f64>, h: &Eval) -> Result<RawPeriod, CalcError> {
    let mut exact = crossings(&nodes, h)?;
    if let Some(&t) = exact.first() {
        return Ok(RawPeriod { enter, leave, exact, min_orb: 0.0, min_jd: t });
    }
    // No crossing seen: locate the closest approach around the nearest node.
    let k = (0..nodes.len()).min_by(|&i, &j| nodes[i].1.abs().total_cmp(&nodes[j].1.abs())).unwrap_or(0);
    let (mut min_jd, mut min_orb) = (nodes[k].0, nodes[k].1.abs());
    let (lo, hi) = (nodes[k.saturating_sub(1)].0, nodes[(k + 1).min(nodes.len() - 1)].0);
    if hi > lo {
        let (t, m) = golden_min(&|t| h(t).map(f64::abs), lo, hi)?;
        if m < min_orb {
            (min_jd, min_orb) = (t, m);
            // The deviation may have crossed zero and come back between two samples.
            let at = nodes.partition_point(|n| n.0 < t);
            nodes.insert(at, (t, h(t)?));
            exact = crossings(&nodes, h)?;
            if let Some(&t) = exact.first() {
                (min_jd, min_orb) = (t, 0.0);
            }
        }
    }
    Ok(RawPeriod { enter, leave, exact, min_orb, min_jd })
}

/// Finds the in-orb periods of one signed deviation function `h`, given its samples `hs` at
/// `times`.
fn scan(times: &[f64], hs: &[f64], orb: f64, h: &Eval) -> Result<Vec<RawPeriod>, CalcError> {
    let g = |t: f64| h(t).map(|v| v.abs() - orb);
    let inside = |i: usize| hs[i].abs() <= orb;
    let mut out = Vec::new();
    // The period being followed: its entry time and its (time, deviation) nodes.
    let mut open: Option<OpenPeriod> = None;
    for i in 0..times.len() {
        if inside(i) {
            if open.is_none() {
                if i == 0 {
                    open = Some((None, Vec::new()));
                } else {
                    let t = bisect(&g, times[i - 1], hs[i - 1].abs() - orb, times[i])?;
                    open = Some((Some(t), vec![(t, h(t)?)]));
                }
            }
            if let Some((_, nodes)) = open.as_mut() {
                nodes.push((times[i], hs[i]));
            }
        } else if let Some((enter, mut nodes)) = open.take() {
            let t = bisect(&g, times[i - 1], hs[i - 1].abs() - orb, times[i])?;
            nodes.push((t, h(t)?));
            out.push(finish(nodes, enter, Some(t), h)?);
        } else if i > 0 && i + 1 < times.len() && !inside(i - 1) && !inside(i + 1) {
            // A brief aspect can fall between samples: look where the deviation bottoms out
            // close enough to the orb to have dipped inside it.
            let (a, m, b) = (hs[i - 1].abs(), hs[i].abs(), hs[i + 1].abs());
            if m <= a && m <= b && m - orb <= (a - m).max(b - m) {
                let (t, v) = golden_min(&|t| h(t).map(f64::abs), times[i - 1], times[i + 1])?;
                if v <= orb {
                    let enter = bisect(&g, times[i - 1], a - orb, t)?;
                    let leave = bisect(&g, t, v - orb, times[i + 1])?;
                    let nodes = vec![(enter, h(enter)?), (t, h(t)?), (leave, h(leave)?)];
                    out.push(finish(nodes, Some(enter), Some(leave), h)?);
                }
            }
        }
    }
    if let Some((enter, nodes)) = open {
        out.push(finish(nodes, enter, None, h)?);
    }
    Ok(out)
}

impl Ephemeris {
    /// Finds when the bodies of two people's charts were in aspect between the real dates
    /// `from` and `to`.
    ///
    /// On each date both charts are taken as they stand for that person on that date: a
    /// progressed chart is advanced by the person's own age (one day per `year_days`), a
    /// natal chart stays fixed. Every body in `bodies_a` is compared with every body in
    /// `bodies_b` by ecliptic longitude against every definition in `targets`, and each
    /// stretch within orb is returned as one [`AspectPeriod`], ordered by its start.
    ///
    /// Dates before a progressed subject's birth progress backwards (converse); start the
    /// search at the later birth to avoid them.
    #[allow(clippy::too_many_arguments)]
    pub fn progressed_aspects(
        &self,
        a: &ChartSubject,
        bodies_a: &[Body],
        b: &ChartSubject,
        bodies_b: &[Body],
        targets: &[AspectOrb],
        from: impl Into<JdUt>,
        to: impl Into<JdUt>,
        opts: &ProgressionOptions,
    ) -> Result<Vec<AspectPeriod>, CalcError> {
        let (from, to) = (from.into().value(), to.into().value());
        let (span, step) = (to - from, opts.step_days);
        // Written so that NaN fails too.
        if !(span >= 0.0 && step > 0.0 && opts.year_days > 0.0 && span.is_finite()) {
            return Err(CalcError::InvalidInterval { days: span, interval: step });
        }
        let n = (span / step).ceil() as usize;
        let times: Vec<f64> = (0..=n).map(|i| (from + i as f64 * step).min(to)).collect();

        let lng = |s: &ChartSubject, body: Body, t: f64| -> Result<f64, CalcError> {
            let jd = s.chart_jd(t, opts.year_days);
            Ok(self.calc_position(body, jd, s.geo, CoordSystem::Ecliptic, opts.perspective)?.lng)
        };
        // Longitudes on the sampling grid, [time][body]. A natal chart is computed once.
        let grid = |s: &ChartSubject, bodies: &[Body]| -> Result<Vec<Vec<f64>>, CalcError> {
            let at = |t: f64| -> Result<Vec<f64>, CalcError> {
                let jd = s.chart_jd(t, opts.year_days);
                let pos = self.calc_positions(bodies, jd, s.geo, CoordSystem::Ecliptic, opts.perspective)?;
                Ok(pos.iter().map(|p| p.lng).collect())
            };
            match s.kind {
                ChartKind::Natal => Ok(vec![at(from)?; times.len()]),
                ChartKind::Progressed => times.iter().map(|&t| at(t)).collect(),
            }
        };
        let (grid_a, grid_b) = (grid(a, bodies_a)?, grid(b, bodies_b)?);

        let mut out = Vec::new();
        let mut hs = vec![0.0; times.len()];
        for (ia, &body_a) in bodies_a.iter().enumerate() {
            for (ib, &body_b) in bodies_b.iter().enumerate() {
                for &aspect in targets {
                    // B may stand `angle` ahead of A or behind it: two separate alignments,
                    // which coincide for the conjunction and the opposition.
                    let one_sided = aspect.angle < 1e-9 || aspect.angle > 180.0 - 1e-9;
                    let offsets = [aspect.angle, -aspect.angle];
                    for &offset in &offsets[..if one_sided { 1 } else { 2 }] {
                        for (i, v) in hs.iter_mut().enumerate() {
                            *v = diff_deg(offset, grid_a[i][ia] - grid_b[i][ib]);
                        }
                        let h = |t: f64| Ok(diff_deg(offset, lng(a, body_a, t)? - lng(b, body_b, t)?));
                        for p in scan(&times, &hs, aspect.orb, &h)? {
                            out.push(AspectPeriod {
                                body_a,
                                body_b,
                                aspect,
                                enter_jd: p.enter,
                                exact_jds: p.exact,
                                leave_jd: p.leave,
                                min_orb: p.min_orb,
                                min_orb_jd: p.min_jd,
                            });
                        }
                    }
                }
            }
        }
        out.sort_by(|x, y| x.enter_jd.unwrap_or(from).total_cmp(&y.enter_jd.unwrap_or(from)));
        Ok(out)
    }
}
