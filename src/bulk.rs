//! Bulk calculation over many records: charts for a list of births, and aspect searches
//! between charts.
//!
//! Every record carries a caller-chosen key (a name, a database id, a reference to the
//! caller's own record) which is never inspected and is handed back with the result, so a
//! caller can filter its records first and still tell which ones came back.
//!
//! Two kinds of function are offered:
//! * [`Ephemeris::bulk_charts`] calculates a chart (chosen bodies, the ascendant among them
//!   if listed) for each record;
//! * the aspect searches compare charts and return only the pairs that make one of the
//!   wanted aspects within its orb: [`Ephemeris::aspects_in_pairs`] within given pairs,
//!   [`Ephemeris::aspects_one_to_many`] between one chart and a list, and
//!   [`Ephemeris::aspects_many_to_many`] between every two records of one list.
//!
//! Records may arrive in any order. They are worked through in date order, so that
//! neighbouring records read neighbouring parts of the ephemeris files, and results are
//! returned in the order of the input. Aspects are matched by the `aspect360` crate, as
//! everywhere in this crate. With the `parallel` feature the work is spread over
//! the threads of the `rayon` pool. A record that cannot be calculated (a date outside the
//! available range, say) never fails the batch: it is reported under `errors` with its key.

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use aspect360::{Aspect360, ToRing360};

use crate::aspects::{AspectOrb, ChartKind, TROPICAL_YEAR};
use crate::ayanamsha::Ayanamsha;
use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::position::{BodyPos, CoordSystem, Perspective};
use crate::spk::{CalcEngine, Ephemeris};
use crate::time::JdUt;

/// One record of a bulk calculation: a moment, a place and the caller's key for it.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChartInput<K> {
    /// Handed back untouched with the result.
    pub key: K,
    pub jd: JdUt,
    pub geo: GeoPos,
}

impl<K> ChartInput<K> {
    pub fn new(key: K, jd: impl Into<JdUt>, geo: GeoPos) -> Self {
        Self { key, jd: jd.into(), geo }
    }
}

/// What [`Ephemeris::bulk_charts`] calculates for each record.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChartOptions {
    /// The bodies wanted, in the order they are returned. List [`Body::Ascendant`] for the
    /// ascendant.
    pub bodies: Vec<Body>,
    pub perspective: Perspective,
    /// If set, longitudes are sidereal.
    pub ayanamsha: Option<Ayanamsha>,
}

impl Default for ChartOptions {
    /// [`Body::DEFAULT_BODIES`], geocentric and tropical.
    fn default() -> Self {
        Self { bodies: Body::DEFAULT_BODIES.to_vec(), perspective: Perspective::Geocentric, ayanamsha: None }
    }
}

/// The chart of one record.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KeyedChart<K> {
    pub key: K,
    pub jd: JdUt,
    /// Ecliptic positions in the order of [`ChartOptions::bodies`].
    pub positions: Vec<BodyPos>,
    /// How the positions were calculated.
    pub engine: CalcEngine,
}

/// A record that could not be calculated.
#[derive(Debug, Clone, PartialEq)]
pub struct KeyedError<K> {
    pub key: K,
    pub error: CalcError,
}

/// Result of [`Ephemeris::bulk_charts`]. Both lists keep the order of the input.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkCharts<K> {
    pub charts: Vec<KeyedChart<K>>,
    pub errors: Vec<KeyedError<K>>,
}

/// Which aspects to look for between two charts, A and B.
// Not `Serialize`: `aspect360::AspectOrb` (0.1.x) does not derive it.
#[derive(Debug, Clone)]
pub struct AspectCriteria {
    /// Bodies taken from chart A.
    pub bodies_a: Vec<Body>,
    /// Bodies taken from chart B. Each body of A is compared with each of these.
    pub bodies_b: Vec<Body>,
    /// The aspects wanted, each with its own orb; see [`aspect_orbs`](crate::aspect_orbs).
    /// An aspect counts on either side: 150° matches B standing 150° ahead of A or behind.
    pub aspects: Vec<AspectOrb>,
    pub perspective: Perspective,
    /// If set, each chart's longitudes are sidereal for its own date.
    pub ayanamsha: Option<Ayanamsha>,
    /// Chart of A to use; only relevant when [`at`](Self::at) is set.
    pub kind_a: ChartKind,
    /// Chart of B to use; only relevant when [`at`](Self::at) is set.
    pub kind_b: ChartKind,
    /// The real date for which progressed charts are taken: a progressed chart is the sky
    /// `age in years` days after birth. `None` compares birth charts whatever the kinds.
    pub at: Option<JdUt>,
    /// Days of real life represented by one day of progression. Default [`TROPICAL_YEAR`].
    pub year_days: f64,
}

impl AspectCriteria {
    /// Birth charts compared geocentrically in the tropical zodiac.
    pub fn new(bodies_a: &[Body], bodies_b: &[Body], aspects: &[AspectOrb]) -> Self {
        Self {
            bodies_a: bodies_a.to_vec(),
            bodies_b: bodies_b.to_vec(),
            aspects: aspects.to_vec(),
            perspective: Perspective::Geocentric,
            ayanamsha: None,
            kind_a: ChartKind::Natal,
            kind_b: ChartKind::Natal,
            at: None,
            year_days: TROPICAL_YEAR,
        }
    }

    /// Uses charts of the given kinds as they stand on the real date `at`.
    pub fn progressed_to(mut self, at: impl Into<JdUt>, kind_a: ChartKind, kind_b: ChartKind) -> Self {
        (self.at, self.kind_a, self.kind_b) = (Some(at.into()), kind_a, kind_b);
        self
    }

    /// The Julian Day whose sky is the chart of a person born at `birth`.
    fn chart_jd(&self, birth: f64, kind: ChartKind) -> f64 {
        match (self.at, kind) {
            (Some(at), ChartKind::Progressed) => birth + (at.value() - birth) / self.year_days,
            _ => birth,
        }
    }

    /// The aspects made between longitudes of chart A and of chart B. Matching is done by
    /// `aspect360`, one target at a time so that nothing is allocated for a pair that fails.
    fn hits(&self, lngs_a: &[f64], lngs_b: &[f64]) -> Vec<AspectHit> {
        let mut out = Vec::new();
        for (&body_a, &lng_a) in self.bodies_a.iter().zip(lngs_a) {
            for (&body_b, &lng_b) in self.bodies_b.iter().zip(lngs_b) {
                let ring_a = lng_a.to_360();
                for aspect in &self.aspects {
                    let found = ring_a.calc_aspect_f64(lng_b, aspect.target(), aspect.orb());
                    if found.matched() {
                        out.push(AspectHit { body_a, body_b, lng_a, lng_b, angle: found.target(), orb: found.orb(), deviation: found.distance() });
                    }
                }
            }
        }
        out
    }
}

/// One aspect found between a body of chart A and a body of chart B.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AspectHit {
    pub body_a: Body,
    pub body_b: Body,
    /// Longitude of the body of A, in degrees.
    pub lng_a: f64,
    /// Longitude of the body of B, in degrees.
    pub lng_b: f64,
    /// The aspect matched, as given in the criteria.
    pub angle: f64,
    /// The orb allowed for it.
    pub orb: f64,
    /// How far the separation is from the exact aspect, signed as `aspect360` reports it
    /// (`AspectResult::distance`). Its absolute value never exceeds `orb`.
    pub deviation: f64,
}

/// Two charts in aspect, with every aspect found between them.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AspectMatch<K> {
    /// Key of the record taken as chart A.
    pub a: K,
    /// Key of the record taken as chart B.
    pub b: K,
    /// Never empty.
    pub hits: Vec<AspectHit>,
}

impl<K> AspectMatch<K> {
    /// The smallest absolute deviation among the hits: a measure for ranking matches.
    pub fn closest(&self) -> f64 {
        self.hits.iter().map(|h| h.deviation.abs()).fold(f64::INFINITY, f64::min)
    }
}

/// Result of an aspect search: the matching pairs only, and the records that failed.
#[derive(Debug, Clone, PartialEq)]
pub struct AspectMatches<K> {
    pub matches: Vec<AspectMatch<K>>,
    pub errors: Vec<KeyedError<K>>,
}

/// Calls `f` for every index in ascending order of `jds`, on several threads with the
/// `parallel` feature, and returns the results in index order.
fn in_date_order<R: Send>(jds: &[f64], f: impl Fn(usize) -> R + Sync) -> Vec<R> {
    let mut order: Vec<usize> = (0..jds.len()).collect();
    order.sort_unstable_by(|&a, &b| jds[a].total_cmp(&jds[b]));
    #[cfg(feature = "parallel")]
    let done: Vec<(usize, R)> = order.par_iter().map(|&i| (i, f(i))).collect();
    #[cfg(not(feature = "parallel"))]
    let done: Vec<(usize, R)> = order.iter().map(|&i| (i, f(i))).collect();
    let mut out: Vec<Option<R>> = jds.iter().map(|_| None).collect();
    for (i, r) in done {
        out[i] = Some(r);
    }
    out.into_iter().flatten().collect()
}

/// Calls `f` for every index below `n` and joins the results in index order.
fn each_index<R: Send>(n: usize, f: impl Fn(usize) -> Vec<R> + Sync) -> Vec<R> {
    #[cfg(feature = "parallel")]
    return (0..n).into_par_iter().flat_map_iter(&f).collect();
    #[cfg(not(feature = "parallel"))]
    (0..n).flat_map(f).collect()
}

type Lngs = Result<Vec<f64>, CalcError>;

impl Ephemeris {
    /// Longitudes of `bodies` for each `(jd, place)`, in the order given.
    fn bulk_longitudes(&self, points: &[(f64, GeoPos)], bodies: &[Body], c: &AspectCriteria) -> Vec<Lngs> {
        let jds: Vec<f64> = points.iter().map(|p| p.0).collect();
        in_date_order(&jds, |i| {
            let (jd, geo) = points[i];
            Ok(self.calc_longitudes(bodies, jd, geo, c.perspective, c.ayanamsha)?.iter().map(|l| l.lng).collect())
        })
    }

    /// Calculates a chart for every record: the ecliptic positions of the bodies in `opts`.
    ///
    /// ```no_run
    /// # use astro_calc_core::*;
    /// # let eph = Ephemeris::open_dir("data").unwrap();
    /// let records = [ChartInput::new("A", 2443121.2359, GeoPos::new(44.9, 12.4, 0.0).unwrap())];
    /// let found = eph.bulk_charts(&records, &ChartOptions::default());
    /// assert_eq!(found.charts[0].key, "A");
    /// ```
    pub fn bulk_charts<K: Clone>(&self, items: &[ChartInput<K>], opts: &ChartOptions) -> BulkCharts<K> {
        let jds: Vec<f64> = items.iter().map(|i| i.jd.value()).collect();
        let places: Vec<GeoPos> = items.iter().map(|i| i.geo).collect();
        let done = in_date_order(&jds, |i| -> Result<Vec<BodyPos>, CalcError> {
            let mut positions = self.calc_positions(&opts.bodies, jds[i], places[i], CoordSystem::Ecliptic, opts.perspective)?;
            if let Some(aya) = opts.ayanamsha {
                let offset = self.ayanamsha(jds[i], aya)?;
                positions.iter_mut().for_each(|p| p.lng = (p.lng - offset).rem_euclid(360.0));
            }
            Ok(positions)
        });
        let mut out = BulkCharts { charts: Vec::with_capacity(items.len()), errors: Vec::new() };
        for (item, result) in items.iter().zip(done) {
            match result {
                Ok(positions) => {
                    let engine = positions.first().map(|p| p.engine).unwrap_or_else(|| self.engine_used(item.jd.value()));
                    out.charts.push(KeyedChart { key: item.key.clone(), jd: item.jd, positions, engine });
                }
                Err(error) => out.errors.push(KeyedError { key: item.key.clone(), error }),
            }
        }
        out
    }

    /// Looks for aspects within each given pair, the first record as chart A and the second
    /// as chart B, and returns the pairs that make at least one. A pair of which either
    /// record fails is left out and the failing record reported under `errors`.
    pub fn aspects_in_pairs<K: Clone>(&self, pairs: &[(ChartInput<K>, ChartInput<K>)], criteria: &AspectCriteria) -> AspectMatches<K> {
        let side = |first: bool, kind, bodies: &[Body]| {
            let points: Vec<_> = pairs
                .iter()
                .map(|p| if first { &p.0 } else { &p.1 })
                .map(|item| (criteria.chart_jd(item.jd.value(), kind), item.geo))
                .collect();
            self.bulk_longitudes(&points, bodies, criteria)
        };
        let lngs_a = side(true, criteria.kind_a, &criteria.bodies_a);
        let lngs_b = side(false, criteria.kind_b, &criteria.bodies_b);
        let mut out = AspectMatches { matches: Vec::new(), errors: Vec::new() };
        for ((pair, a), b) in pairs.iter().zip(lngs_a).zip(lngs_b) {
            match (a, b) {
                (Ok(a), Ok(b)) => {
                    let hits = criteria.hits(&a, &b);
                    if !hits.is_empty() {
                        out.matches.push(AspectMatch { a: pair.0.key.clone(), b: pair.1.key.clone(), hits });
                    }
                }
                (a, b) => {
                    out.errors.extend(a.err().map(|error| KeyedError { key: pair.0.key.clone(), error }));
                    out.errors.extend(b.err().map(|error| KeyedError { key: pair.1.key.clone(), error }));
                }
            }
        }
        out
    }

    /// Looks for aspects between one chart (A) and each of `targets` (B), and returns the
    /// targets that make at least one. Fails only if the source chart cannot be calculated.
    pub fn aspects_one_to_many<K: Clone>(
        &self,
        source: &ChartInput<K>,
        targets: &[ChartInput<K>],
        criteria: &AspectCriteria,
    ) -> Result<AspectMatches<K>, CalcError> {
        let source_jd = criteria.chart_jd(source.jd.value(), criteria.kind_a);
        let source_lngs: Vec<f64> = self
            .calc_longitudes(&criteria.bodies_a, source_jd, source.geo, criteria.perspective, criteria.ayanamsha)?
            .iter()
            .map(|l| l.lng)
            .collect();
        let points: Vec<_> = targets.iter().map(|t| (criteria.chart_jd(t.jd.value(), criteria.kind_b), t.geo)).collect();
        let jds: Vec<f64> = points.iter().map(|p| p.0).collect();
        // Longitudes are dropped as soon as they have been tested.
        let done = in_date_order(&jds, |i| -> Result<Vec<AspectHit>, CalcError> {
            let (jd, geo) = points[i];
            let lngs: Vec<f64> =
                self.calc_longitudes(&criteria.bodies_b, jd, geo, criteria.perspective, criteria.ayanamsha)?.iter().map(|l| l.lng).collect();
            Ok(criteria.hits(&source_lngs, &lngs))
        });
        let mut out = AspectMatches { matches: Vec::new(), errors: Vec::new() };
        for (target, result) in targets.iter().zip(done) {
            match result {
                Ok(hits) if hits.is_empty() => {}
                Ok(hits) => out.matches.push(AspectMatch { a: source.key.clone(), b: target.key.clone(), hits }),
                Err(error) => out.errors.push(KeyedError { key: target.key.clone(), error }),
            }
        }
        Ok(out)
    }

    /// Looks for aspects between every two records of one list and returns the pairs that
    /// make at least one.
    ///
    /// Each two records are tried both ways round, since a body of A to a body of B is not
    /// the same test as the reverse; the record whose bodies were taken from
    /// [`bodies_a`](AspectCriteria::bodies_a) is `a` in the match. When both body lists and
    /// both chart kinds are the same the two ways are equivalent and each pair is reported
    /// once, the earlier record of the list as `a`.
    ///
    /// Positions are calculated once per record, but the number of comparisons grows with
    /// the square of the list: 10,000 records make 50 million pairs.
    pub fn aspects_many_to_many<K: Clone>(&self, items: &[ChartInput<K>], criteria: &AspectCriteria) -> AspectMatches<K> {
        let symmetric = criteria.bodies_a == criteria.bodies_b && (criteria.at.is_none() || criteria.kind_a == criteria.kind_b);
        let side = |kind, bodies: &[Body]| {
            let points: Vec<_> = items.iter().map(|i| (criteria.chart_jd(i.jd.value(), kind), i.geo)).collect();
            self.bulk_longitudes(&points, bodies, criteria)
        };
        let lngs_a = side(criteria.kind_a, &criteria.bodies_a);
        let lngs_b = if symmetric { Vec::new() } else { side(criteria.kind_b, &criteria.bodies_b) };
        let lngs_b = if symmetric { &lngs_a } else { &lngs_b };

        let n = items.len();
        let found: Vec<(usize, usize, Vec<AspectHit>)> = each_index(n, |i| {
            let Ok(a) = &lngs_a[i] else { return Vec::new() };
            let start = if symmetric { i + 1 } else { 0 };
            (start..n)
                .filter(|&j| j != i)
                .filter_map(|j| {
                    let hits = criteria.hits(a, lngs_b[j].as_ref().ok()?);
                    (!hits.is_empty()).then_some((i, j, hits))
                })
                .collect()
        });

        let matches = found.into_iter().map(|(i, j, hits)| AspectMatch { a: items[i].key.clone(), b: items[j].key.clone(), hits }).collect();
        let mut errors = Vec::new();
        for (i, item) in items.iter().enumerate() {
            let failed = lngs_a[i].as_ref().err().or_else(|| lngs_b[i].as_ref().err());
            errors.extend(failed.map(|error| KeyedError { key: item.key.clone(), error: error.clone() }));
        }
        AspectMatches { matches, errors }
    }
}
