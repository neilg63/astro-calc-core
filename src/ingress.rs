//! Ingresses: the moments a body crosses from one equal division of the zodiac into the
//! next — a sign (30°), a nakshatra (13°20′), a pada (3°20′) or any other span that divides
//! the circle — on the tropical zodiac or a sidereal one.
//!
//! Unlike the Moon − Sun elongation behind the tithis, one body's longitude need not keep
//! growing: a planet may enter a division, turn retrograde and leave by the way it came. The
//! search therefore walks along the body's path in steps too short to reach a boundary, and
//! closes in once a step lands in another division.

use crate::ayanamsha::Ayanamsha;
use crate::body::Body;
use crate::error::CalcError;
use crate::geo::GeoPos;
use crate::position::{CoordSystem, Perspective};
use crate::spk::{CalcEngine, Ephemeris};
use crate::time::JdUt;

/// A body's passage from one division of the zodiac into a neighbouring one.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LongitudeCrossing {
    pub jd: JdUt,
    /// The division entered, counted from 0 at 0° of the zodiac.
    pub into: u32,
    /// The longitude of the boundary crossed, in degrees.
    pub lng: f64,
    /// Whether the body was moving forwards through the zodiac. A retrograde body enters a
    /// division at its end.
    pub direct: bool,
    /// The engine that produced the time.
    #[cfg_attr(feature = "serde", serde(default))]
    pub engine: CalcEngine,
}

/// A body's stay in the division of the zodiac it is in at some moment.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LongitudeStay {
    /// The division, counted from 0 at 0° of the zodiac.
    pub division: u32,
    /// How the body last came into it. `None` if that is further back than the search went.
    pub entered: Option<LongitudeCrossing>,
    /// How it next leaves; `into` is the division it moves on to. `None` if that is further
    /// ahead than the search went.
    pub left: Option<LongitudeCrossing>,
}

/// The shortest step of the walk, in days: under ten seconds.
const MIN_STEP: f64 = 1e-4;
/// How closely a crossing is timed, in days: about a millisecond.
const TIME_TOLERANCE: f64 = 1e-8;

/// An upper bound on how fast a body's geocentric speed in longitude can change, in degrees
/// per day per day. Generous: it only shortens the steps taken where a body is slow.
fn max_acceleration(body: Body) -> f64 {
    match body {
        Body::Jupiter | Body::Saturn | Body::Uranus | Body::Neptune | Body::Pluto => 0.01,
        Body::Sun | Body::Venus | Body::Mars => 0.06,
        Body::Mercury => 0.25,
        // The Moon, its nodes and apogees, the small bodies and anything unknown.
        _ => 1.0,
    }
}

/// The equal divisions of the circle that a span in degrees makes.
#[derive(Clone, Copy)]
struct Divisions {
    span: f64,
    count: u32,
}

impl Divisions {
    fn new(span: f64) -> Self {
        Self { span, count: ((360.0 / span).round() as u32).max(1) }
    }

    fn of(self, lng: f64) -> u32 {
        ((lng.rem_euclid(360.0) / self.span) as u32).min(self.count - 1)
    }

    /// Degrees from a longitude to the nearer boundary of a division it lies in.
    fn room(self, lng: f64, division: u32) -> f64 {
        let within = lng.rem_euclid(360.0) - division as f64 * self.span;
        within.min(self.span - within).max(0.0)
    }
}

impl Ephemeris {
    /// Geocentric longitude on the zodiac asked for, and its speed. The speed is the tropical
    /// one: an ayanamsha drifts by under 0.00004° a day, which no step taken here notices.
    fn zodiac_lng_speed(&self, body: Body, t: f64, ayanamsha: Option<Ayanamsha>) -> Result<(f64, f64), CalcError> {
        let p = self.calc_position(body, t, GeoPos::default(), CoordSystem::Ecliptic, Perspective::Geocentric)?;
        let offset = match ayanamsha {
            Some(aya) => self.ayanamsha(t, aya)?,
            None => 0.0,
        };
        Ok(((p.lng - offset).rem_euclid(360.0), p.lng_speed))
    }

    /// Walks from `start` forwards (`forwards`) or backwards in time until the body is in
    /// another division, for at most `max_days`, and times the crossing. `None` if the body
    /// stays put that long, or the ephemeris ends first.
    fn longitude_crossing_from(
        &self,
        body: Body,
        start: f64,
        forwards: bool,
        divisions: Divisions,
        ayanamsha: Option<Ayanamsha>,
        max_days: f64,
    ) -> Result<Option<LongitudeCrossing>, CalcError> {
        let dir = if forwards { 1.0 } else { -1.0 };
        let accel = max_acceleration(body);
        let (mut lng, mut speed) = self.zodiac_lng_speed(body, start, ayanamsha)?;
        let home = divisions.of(lng);
        let mut t = start;
        loop {
            // Neither the present speed nor the greatest change in it can carry the body more
            // than half the way to a boundary each, so a step never passes over a division.
            let room = divisions.room(lng, home);
            let step = (0.5 * room / speed.abs()).min((room / accel).sqrt()).max(MIN_STEP);
            let next = t + dir * step;
            if (next - start).abs() > max_days {
                return Ok(None);
            }
            let (next_lng, next_speed) = match self.zodiac_lng_speed(body, next, ayanamsha) {
                Err(e) if e.is_out_of_range() => return Ok(None),
                found => found?,
            };
            if divisions.of(next_lng) != home {
                // `inside` is in the home division and `outside` is not: halve the gap.
                let (mut inside, mut outside) = (t, next);
                while (outside - inside).abs() > TIME_TOLERANCE {
                    let mid = 0.5 * (inside + outside);
                    if divisions.of(self.zodiac_lng_speed(body, mid, ayanamsha)?.0) == home {
                        inside = mid;
                    } else {
                        outside = mid;
                    }
                }
                // The crossing is the first moment in the division entered: the far side of
                // the gap going forwards, the near side going back.
                let jd = if forwards { outside } else { inside };
                let (lng, speed) = self.zodiac_lng_speed(body, jd, ayanamsha)?;
                let into = divisions.of(lng);
                let direct = speed >= 0.0;
                let boundary = if direct { into } else { (into + 1) % divisions.count };
                return Ok(Some(LongitudeCrossing { jd: JdUt(jd), into, lng: boundary as f64 * divisions.span, direct, engine: self.engine_used(jd) }));
            }
            (t, lng, speed) = (next, next_lng, next_speed);
        }
    }

    /// The stay of a body in the division of the zodiac it is in at `jd`: when it last came
    /// in and when it next leaves, each looked for over at most `max_days`.
    ///
    /// `span` is the size of a division in degrees and must divide 360: 30 for the signs,
    /// 360/27 for the nakshatras. With an ayanamsha the zodiac is sidereal. Geocentric.
    ///
    /// A body that turns retrograde may leave by the boundary it came in at; see
    /// [`LongitudeCrossing::direct`].
    pub fn longitude_stay(&self, body: Body, jd: impl Into<JdUt>, span: f64, ayanamsha: Option<Ayanamsha>, max_days: f64) -> Result<LongitudeStay, CalcError> {
        let jd = jd.into().value();
        let divisions = Divisions::new(span);
        Ok(LongitudeStay {
            division: divisions.of(self.zodiac_lng_speed(body, jd, ayanamsha)?.0),
            entered: self.longitude_crossing_from(body, jd, false, divisions, ayanamsha, max_days)?,
            left: self.longitude_crossing_from(body, jd, true, divisions, ayanamsha, max_days)?,
        })
    }

    /// The first time after `jd` that a body moves into another division of the zodiac, if
    /// it does within `max_days`. See [`longitude_stay`](Self::longitude_stay) for `span`.
    pub fn next_longitude_crossing(
        &self,
        body: Body,
        jd: impl Into<JdUt>,
        span: f64,
        ayanamsha: Option<Ayanamsha>,
        max_days: f64,
    ) -> Result<Option<LongitudeCrossing>, CalcError> {
        self.longitude_crossing_from(body, jd.into().value(), true, Divisions::new(span), ayanamsha, max_days)
    }

    /// The next `count` moves of a body into another division of the zodiac after `jd`, in
    /// order: with a span of 360/27 and the Moon, the beginnings of its next `count`
    /// nakshatras. Fewer are returned if the rest lie more than `max_days` after `jd`, or
    /// beyond the end of the ephemeris. A retrograde body's re-entries count among them.
    /// See [`longitude_stay`](Self::longitude_stay) for `span`.
    pub fn next_longitude_crossings(
        &self,
        body: Body,
        jd: impl Into<JdUt>,
        count: usize,
        span: f64,
        ayanamsha: Option<Ayanamsha>,
        max_days: f64,
    ) -> Result<Vec<LongitudeCrossing>, CalcError> {
        let start = jd.into().value();
        let divisions = Divisions::new(span);
        let mut t = start;
        let mut out = Vec::with_capacity(count);
        while out.len() < count {
            let Some(crossing) = self.longitude_crossing_from(body, t, true, divisions, ayanamsha, max_days - (t - start))? else { break };
            t = crossing.jd.value();
            out.push(crossing);
        }
        Ok(out)
    }

    /// Every move of a body into another division of the zodiac between two times, in order:
    /// with a span of 360/27 and the Moon, the beginnings of its nakshatras. A retrograde
    /// body's re-entries are among them. See [`longitude_stay`](Self::longitude_stay) for `span`.
    pub fn longitude_crossings(
        &self,
        body: Body,
        from: impl Into<JdUt>,
        to: impl Into<JdUt>,
        span: f64,
        ayanamsha: Option<Ayanamsha>,
    ) -> Result<Vec<LongitudeCrossing>, CalcError> {
        let (mut t, to) = (from.into().value(), to.into().value());
        let divisions = Divisions::new(span);
        let mut out = Vec::new();
        while let Some(crossing) = self.longitude_crossing_from(body, t, true, divisions, ayanamsha, to - t)? {
            t = crossing.jd.value();
            out.push(crossing);
        }
        Ok(out)
    }
}
