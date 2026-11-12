//! Chart angles and house cusps.
//!
//! All angles derive from the ARMC (apparent sidereal time at the observer's meridian), the
//! true obliquity and the geographic latitude, using the classical spherical-astronomy
//! constructions (e.g. Holden, "The Elements of House Division"; Meeus ch. 13–14):
//!
//! * `cusp(θ, P)` is the ecliptic point on the eastern side of the great circle through the north
//!   and south points of a horizon with pole height `P`, at oblique ascension `θ + 90°`.
//!   The ascendant is `cusp(ARMC, φ)`; most quadrant systems differ only in how they choose
//!   `θ` and `P` for each intermediate cusp.
//! * Placidus and Koch fail where the ecliptic becomes circumpolar (|φ| ≥ 90° − ε); like the
//!   reference implementation they then fall back to Porphyry and report it.

use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt;
use std::str::FromStr;

use crate::geo::GeoPos;
use crate::math::norm_deg;
use crate::models::{DeltaTModel, EarthOrientation, ut_to_tt};
use crate::time::JdUt;

/// Chart angles and auxiliary points, all ecliptic longitudes in degrees except `armc`
/// (right ascension of the meridian, degrees).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Angles {
    pub ascendant: f64,
    pub mc: f64,
    pub armc: f64,
    pub vertex: f64,
    /// Equatorial ascendant (east point).
    pub equasc: f64,
    /// Co-ascendant (W. Koch).
    pub coasc1: f64,
    /// Co-ascendant (M. Munkasey).
    pub coasc2: f64,
    /// Polar ascendant (M. Munkasey).
    pub polasc: f64,
}

impl Angles {
    pub fn descendant(&self) -> f64 {
        norm_deg(self.ascendant + 180.0)
    }

    pub fn ic(&self) -> f64 {
        norm_deg(self.mc + 180.0)
    }
}

/// House division methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[non_exhaustive]
pub enum HouseSystem {
    /// `W`: each sign is a house, starting with the ascendant's sign.
    WholeSign,
    /// `E` (also `A`): 30° houses from the ascendant.
    Equal,
    /// `V`: equal houses with the ascendant in the middle of the 1st.
    Vehlow,
    /// `O`: quadrants trisected in ecliptic longitude.
    Porphyry,
    /// `P`: trisection of semi-arcs in time.
    Placidus,
    /// `K`: birthplace system, trisection of the MC degree's semi-arc.
    Koch,
    /// `B`: trisection of the ascendant's semi-arc in right ascension.
    Alcabitius,
    /// `C`: prime vertical divided into 30° arcs.
    Campanus,
    /// `R`: celestial equator divided into 30° arcs.
    Regiomontanus,
    /// `T`: Polich/Page "topocentric" system.
    Topocentric,
    /// `M`: equator divided from the ARMC, projected via the ecliptic poles.
    Morinus,
    /// `X`: axial rotation / meridian houses.
    Meridian,
    /// `H`: horizon (azimuthal) houses.
    Horizontal,
    /// `U`: Krusinski-Pisa-Goelzer; the great circle through the ascendant and the zenith
    /// divided into 30° arcs, projected onto the ecliptic along meridians.
    Krusinski,
}

impl HouseSystem {
    pub const ALL: [HouseSystem; 14] = [
        HouseSystem::WholeSign,
        HouseSystem::Equal,
        HouseSystem::Vehlow,
        HouseSystem::Porphyry,
        HouseSystem::Placidus,
        HouseSystem::Koch,
        HouseSystem::Alcabitius,
        HouseSystem::Campanus,
        HouseSystem::Regiomontanus,
        HouseSystem::Topocentric,
        HouseSystem::Morinus,
        HouseSystem::Meridian,
        HouseSystem::Horizontal,
        HouseSystem::Krusinski,
    ];

    /// Conventional one-letter code.
    pub const fn code(self) -> char {
        match self {
            HouseSystem::WholeSign => 'W',
            HouseSystem::Equal => 'E',
            HouseSystem::Vehlow => 'V',
            HouseSystem::Porphyry => 'O',
            HouseSystem::Placidus => 'P',
            HouseSystem::Koch => 'K',
            HouseSystem::Alcabitius => 'B',
            HouseSystem::Campanus => 'C',
            HouseSystem::Regiomontanus => 'R',
            HouseSystem::Topocentric => 'T',
            HouseSystem::Morinus => 'M',
            HouseSystem::Meridian => 'X',
            HouseSystem::Horizontal => 'H',
            HouseSystem::Krusinski => 'U',
        }
    }

    pub fn from_code(c: char) -> Option<Self> {
        Some(match c.to_ascii_uppercase() {
            'W' => HouseSystem::WholeSign,
            'E' | 'A' => HouseSystem::Equal,
            'V' => HouseSystem::Vehlow,
            'O' => HouseSystem::Porphyry,
            'P' => HouseSystem::Placidus,
            'K' => HouseSystem::Koch,
            'B' => HouseSystem::Alcabitius,
            'C' => HouseSystem::Campanus,
            'R' => HouseSystem::Regiomontanus,
            'T' => HouseSystem::Topocentric,
            'M' => HouseSystem::Morinus,
            'X' => HouseSystem::Meridian,
            'H' => HouseSystem::Horizontal,
            'U' => HouseSystem::Krusinski,
            _ => return None,
        })
    }

    pub const fn name(self) -> &'static str {
        match self {
            HouseSystem::WholeSign => "Whole Sign",
            HouseSystem::Equal => "Equal",
            HouseSystem::Vehlow => "Vehlow Equal",
            HouseSystem::Porphyry => "Porphyry",
            HouseSystem::Placidus => "Placidus",
            HouseSystem::Koch => "Koch",
            HouseSystem::Alcabitius => "Alcabitius",
            HouseSystem::Campanus => "Campanus",
            HouseSystem::Regiomontanus => "Regiomontanus",
            HouseSystem::Topocentric => "Polich/Page (Topocentric)",
            HouseSystem::Morinus => "Morinus",
            HouseSystem::Meridian => "Meridian (Axial Rotation)",
            HouseSystem::Horizontal => "Horizontal (Azimuthal)",
            HouseSystem::Krusinski => "Krusinski-Pisa-Goelzer",
        }
    }
}

impl fmt::Display for HouseSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.code())
    }
}

impl FromStr for HouseSystem {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut chars = s.trim().chars();
        match (chars.next(), chars.next()) {
            (Some(c), None) => Self::from_code(c).ok_or_else(|| format!("unknown house system: {s}")),
            _ => Err(format!("unknown house system: {s}")),
        }
    }
}

/// House cusps 1–12 (ecliptic longitudes, degrees).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Houses {
    pub system: HouseSystem,
    pub cusps: [f64; 12],
    /// Set when the requested system is undefined at this latitude and another was used.
    pub fallback: Option<HouseSystem>,
}

impl Houses {
    /// House number (1–12) containing an ecliptic longitude.
    pub fn house_of(&self, lng: f64) -> u8 {
        for i in 0..12 {
            let start = self.cusps[i];
            let end = self.cusps[(i + 1) % 12];
            let span = (end - start).rem_euclid(360.0);
            if (lng - start).rem_euclid(360.0) < span {
                return i as u8 + 1;
            }
        }
        1
    }
}

/// The frame shared by all angle and house calculations at one moment and place.
#[derive(Debug, Clone, Copy)]
struct Sky {
    /// ARMC, radians.
    ramc: f64,
    /// True obliquity, radians.
    eps: f64,
    /// Geographic latitude, radians.
    phi: f64,
}

impl Sky {
    fn new(jd: f64, geo: &GeoPos) -> Self {
        let eo = EarthOrientation::new(jd, ut_to_tt(DeltaTModel::Observed, jd));
        Self {
            ramc: (eo.gast + geo.lng.to_radians()).rem_euclid(2.0 * PI),
            eps: eo.eps_true,
            phi: geo.lat.to_radians(),
        }
    }

    /// Ecliptic point at oblique ascension θ + 90° for pole height `pole` (radians in, degrees out).
    fn cusp(&self, theta: f64, pole: f64) -> f64 {
        let (st, ct) = theta.sin_cos();
        let (se, ce) = self.eps.sin_cos();
        norm_deg(ct.atan2(-(st * ce + pole.tan() * se)).to_degrees())
    }

    /// Ecliptic longitude of the point with right ascension `ra` (radians in, degrees out).
    fn lng_from_ra(&self, ra: f64) -> f64 {
        norm_deg(ra.sin().atan2(ra.cos() * self.eps.cos()).to_degrees())
    }

    /// Declination of the ecliptic point at `lng` degrees (radians).
    fn decl(&self, lng: f64) -> f64 {
        (self.eps.sin() * lng.to_radians().sin()).asin()
    }

    fn mc(&self) -> f64 {
        self.lng_from_ra(self.ramc)
    }

    /// Ascendant. Beyond the polar circles the formula can return the descending point of
    /// the ecliptic; the ascendant is kept in the eastern half, i.e. 0–180° after the MC.
    fn asc(&self) -> f64 {
        self.asc_at(self.ramc)
    }

    /// Ascendant for another ARMC, kept east of that ARMC's MC.
    fn asc_at(&self, ramc: f64) -> f64 {
        let asc = self.cusp(ramc, self.phi);
        if (asc - self.lng_from_ra(ramc)).rem_euclid(360.0) > 180.0 { norm_deg(asc + 180.0) } else { asc }
    }

    /// True when the raw ascendant formula returned the western point and was corrected.
    fn asc_flipped(&self) -> bool {
        (self.cusp(self.ramc, self.phi) - self.mc()).rem_euclid(360.0) > 180.0
    }

    /// Vertex: the western intersection of ecliptic and prime vertical. Undefined on the
    /// equator, where by convention it is 0° or 180° depending on the ARMC's half-circle.
    fn vertex(&self) -> f64 {
        if self.phi == 0.0 {
            return if self.ramc < PI { 0.0 } else { 180.0 };
        }
        let vx = self.cusp(self.ramc + PI, FRAC_PI_2 - self.phi);
        // keep it in the western half (180–360° after the MC)
        if (vx - self.mc()).rem_euclid(360.0) < 180.0 { norm_deg(vx + 180.0) } else { vx }
    }

    fn angles(&self) -> Angles {
        let r = self.ramc;
        Angles {
            ascendant: self.asc(),
            mc: self.mc(),
            armc: r.to_degrees(),
            vertex: self.vertex(),
            equasc: self.cusp(r, 0.0),
            coasc1: self.cusp(r, -self.phi),
            coasc2: self.cusp(r, FRAC_PI_2 - self.phi),
            polasc: norm_deg(self.cusp(r, -self.phi) + 180.0),
        }
    }

    /// Ascendant's ecliptic semi-arc is undefined where the ecliptic can be circumpolar.
    fn is_polar(&self) -> bool {
        self.phi.abs() >= FRAC_PI_2 - self.eps
    }

    fn houses(&self, system: HouseSystem) -> Houses {
        let asc = self.asc();
        let mc = self.mc();
        let r = self.ramc;
        let d30 = PI / 6.0;
        let mut fallback = None;
        let quad = |c11: f64, c12: f64, c2: f64, c3: f64| quadrant_cusps(asc, mc, c11, c12, c2, c3);
        let cusps = match system {
            HouseSystem::WholeSign => equal_from((asc / 30.0).floor() * 30.0),
            HouseSystem::Equal => equal_from(asc),
            HouseSystem::Vehlow => equal_from(asc - 15.0),
            HouseSystem::Porphyry => porphyry(asc, mc),
            HouseSystem::Regiomontanus => {
                let pole = |h: f64| (self.phi.tan() * h.sin()).atan();
                quad(
                    self.cusp(r - 2.0 * d30, pole(d30)),
                    self.cusp(r - d30, pole(2.0 * d30)),
                    self.cusp(r + d30, pole(4.0 * d30)),
                    self.cusp(r + 2.0 * d30, pole(5.0 * d30)),
                )
            }
            HouseSystem::Campanus => {
                // h measured on the prime vertical from the meridian.
                let c = |h: f64| {
                    let a = (self.phi.cos() * h.tan()).atan();
                    let a = if h > FRAC_PI_2 { a + PI } else { a };
                    let pole = (self.phi.sin() * h.sin()).asin();
                    self.cusp(r + a - FRAC_PI_2, pole)
                };
                quad(c(d30), c(2.0 * d30), c(4.0 * d30), c(5.0 * d30))
            }
            HouseSystem::Topocentric => {
                let pole = |k: f64| (k / 3.0 * self.phi.tan()).atan();
                quad(
                    self.cusp(r - 2.0 * d30, pole(1.0)),
                    self.cusp(r - d30, pole(2.0)),
                    self.cusp(r + d30, pole(2.0)),
                    self.cusp(r + 2.0 * d30, pole(1.0)),
                )
            }
            HouseSystem::Alcabitius => {
                let asc_ra = self.ra_of(asc);
                let sa = (asc_ra - r).rem_euclid(2.0 * PI); // diurnal semi-arc of the ascendant
                let na = PI - sa;
                quad(
                    self.lng_from_ra(r + sa / 3.0),
                    self.lng_from_ra(r + 2.0 * sa / 3.0),
                    self.lng_from_ra(r + PI - 2.0 * na / 3.0),
                    self.lng_from_ra(r + PI - na / 3.0),
                )
            }
            HouseSystem::Koch | HouseSystem::Placidus if self.is_polar() => {
                fallback = Some(HouseSystem::Porphyry);
                porphyry(asc, mc)
            }
            HouseSystem::Koch => {
                let ad = (self.phi.tan() * self.decl(mc).tan()).asin();
                // Diurnal semi-arc of the MC degree, trisected on both sides of the ascendant.
                let dsa = FRAC_PI_2 + ad;
                quad(
                    self.asc_at(r - 2.0 * dsa / 3.0),
                    self.asc_at(r - dsa / 3.0),
                    self.asc_at(r + dsa / 3.0),
                    self.asc_at(r + 2.0 * dsa / 3.0),
                )
            }
            HouseSystem::Placidus => quad(
                self.placidus(1.0 / 3.0, true),
                self.placidus(2.0 / 3.0, true),
                self.placidus(2.0 / 3.0, false),
                self.placidus(1.0 / 3.0, false),
            ),
            HouseSystem::Morinus => {
                let mut c = [0.0; 12];
                for (i, v) in c.iter_mut().enumerate() {
                    let ra = r + FRAC_PI_2 + i as f64 * d30;
                    *v = norm_deg((ra.sin() * self.eps.cos()).atan2(ra.cos()).to_degrees());
                }
                c
            }
            HouseSystem::Meridian => {
                let mut c = [0.0; 12];
                for (i, v) in c.iter_mut().enumerate() {
                    *v = self.lng_from_ra(r + FRAC_PI_2 + i as f64 * d30);
                }
                c
            }
            HouseSystem::Horizontal => self.horizontal(),
            HouseSystem::Krusinski => self.krusinski(asc),
        };
        let mut cusps = cusps;
        // Beyond the polar circles, when the ascendant had to be moved to the eastern side, the
        // great-circle systems keep their raw intermediate cusps, which then lie opposite.
        if self.asc_flipped()
            && matches!(system, HouseSystem::Regiomontanus | HouseSystem::Campanus | HouseSystem::Topocentric)
        {
            for (i, c) in cusps.iter_mut().enumerate() {
                if i != 0 && i != 6 {
                    *c = norm_deg(*c + 180.0);
                }
            }
        }
        Houses { system, cusps, fallback }
    }

    /// Horizon (azimuthal) houses: cusp k is the ecliptic point on the vertical circle at
    /// azimuth A₁₀ ± 30°·(10 − k), on the half-circle through that azimuth's horizon point.
    /// A₁₀ is the meridian side on which the MC culminates; the sequence runs clockwise
    /// (decreasing azimuth) in the northern hemisphere and anticlockwise in the southern
    /// (the equator follows the southern convention, matching the reference).
    fn horizontal(&self) -> [f64; 12] {
        let north_hemi = self.phi > 0.0;
        let mc_dec = self.decl(self.mc());
        let mc_south = mc_dec < self.phi;
        let a10 = if mc_south { PI } else { 0.0 };
        let dir = if north_hemi { 1.0 } else { -1.0 };
        let (sl, cl) = self.ramc.sin_cos();
        let (sp, cp) = self.phi.sin_cos();
        let up = [cp * cl, cp * sl, sp];
        let north = [-sp * cl, -sp * sl, cp];
        let east = [-sl, cl, 0.0];
        let (se, ce) = self.eps.sin_cos();
        let ecl_pole = [0.0, -se, ce];
        std::array::from_fn(|i| {
            let k = (i + 1) as f64;
            let az = a10 + dir * (10.0 - k) * PI / 6.0;
            let (sa, ca) = az.sin_cos();
            let h = [0, 1, 2].map(|j| ca * north[j] + sa * east[j]);
            let n = crate::math::cross(up, h);
            let mut d = crate::math::cross(n, ecl_pole);
            if crate::math::dot(d, h) < 0.0 {
                d = crate::math::scale(d, -1.0);
            }
            norm_deg((d[1] * ce + d[2] * se).atan2(d[0]).to_degrees())
        })
    }

    /// Krusinski-Pisa-Goelzer houses (B. Krusinski, 1995): the great circle through the
    /// ascendant and the zenith is divided into twelve 30° arcs starting at the ascendant,
    /// and each division point is carried to the ecliptic along its meridian (the great
    /// circle through the celestial poles). The zenith and nadir fall on the 10th and 4th
    /// cusps, which are therefore the MC and IC.
    fn krusinski(&self, asc: f64) -> [f64; 12] {
        let (sl, cl) = self.ramc.sin_cos();
        let (sp, cp) = self.phi.sin_cos();
        let up = [cp * cl, cp * sl, sp];
        let (se, ce) = self.eps.sin_cos();
        let (sa, ca) = asc.to_radians().sin_cos();
        let east = [ca, sa * ce, sa * se];
        std::array::from_fn(|i| {
            // Houses are numbered against the diurnal motion: downwards from the ascendant.
            let (st, ct) = (-(i as f64) * PI / 6.0).sin_cos();
            let p = [0, 1, 2].map(|j| ct * east[j] + st * up[j]);
            self.lng_from_ra(p[1].atan2(p[0]))
        })
    }

    fn ra_of(&self, lng: f64) -> f64 {
        let l = lng.to_radians();
        (l.sin() * self.eps.cos()).atan2(l.cos())
    }

    /// Placidus cusp at semi-arc fraction `f` (1/3 or 2/3) above (`upper`) or below the horizon,
    /// on the eastern side. Iterates on the cusp's declination.
    fn placidus(&self, f: f64, upper: bool) -> f64 {
        let r = self.ramc;
        let ra_for = |lng: f64| {
            let ad = (self.phi.tan() * self.decl(lng).tan()).clamp(-1.0, 1.0).asin();
            if upper {
                r + f * (FRAC_PI_2 + ad)
            } else {
                r + PI - f * (FRAC_PI_2 - ad)
            }
        };
        // Start from the equal-division right ascension.
        let mut lng = self.lng_from_ra(if upper { r + f * FRAC_PI_2 } else { r + PI - f * FRAC_PI_2 });
        for _ in 0..50 {
            let next = self.lng_from_ra(ra_for(lng));
            let delta = crate::math::diff_deg(lng, next).abs();
            lng = next;
            if delta < 1e-10 {
                break;
            }
        }
        lng
    }
}

fn equal_from(start: f64) -> [f64; 12] {
    std::array::from_fn(|i| norm_deg(start + 30.0 * i as f64))
}

fn porphyry(asc: f64, mc: f64) -> [f64; 12] {
    let east = (asc - mc).rem_euclid(360.0) / 3.0; // MC → ASC quadrant
    let west = 60.0 - east; // ASC → IC is 180 − 3·east, per third
    quadrant_cusps(
        asc,
        mc,
        norm_deg(mc + east),
        norm_deg(mc + 2.0 * east),
        norm_deg(asc + west),
        norm_deg(asc + 2.0 * west),
    )
}

/// Assembles 12 cusps from the angles and the four eastern intermediate cusps (11, 12, 2, 3);
/// the western ones are their opposites.
fn quadrant_cusps(asc: f64, mc: f64, c11: f64, c12: f64, c2: f64, c3: f64) -> [f64; 12] {
    let opp = |x: f64| norm_deg(x + 180.0);
    [asc, c2, c3, opp(mc), opp(c11), opp(c12), opp(asc), opp(c2), opp(c3), mc, c11, c12]
}

/// Ascendant (ecliptic longitude, degrees) for a UT Julian Day and place.
pub fn ascendant(jd: impl Into<JdUt>, geo: GeoPos) -> f64 {
    Sky::new(jd.into().value(), &geo).asc()
}

/// Ascendant, MC, ARMC, vertex and auxiliary points.
pub fn angles(jd: impl Into<JdUt>, geo: GeoPos) -> Angles {
    Sky::new(jd.into().value(), &geo).angles()
}

/// House cusps for one system.
pub fn houses(jd: impl Into<JdUt>, geo: GeoPos, system: HouseSystem) -> Houses {
    Sky::new(jd.into().value(), &geo).houses(system)
}

/// Angles plus cusps for several systems, sharing one sidereal-time calculation.
pub fn house_sets(jd: impl Into<JdUt>, geo: GeoPos, systems: &[HouseSystem]) -> (Angles, Vec<Houses>) {
    let sky = Sky::new(jd.into().value(), &geo);
    (sky.angles(), systems.iter().map(|&s| sky.houses(s)).collect())
}
