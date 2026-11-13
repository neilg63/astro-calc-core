//! Celestial body identifiers.
//!
//! [`Body`] is a small `Copy + Eq + Hash` value (16 bytes) so it can be used as a map key,
//! passed by value through bulk calculations and shared across threads freely.
//!
//! * Well-known bodies and calculated points are plain variants with a two-letter key.
//! * Lesser-known minor bodies and stars use [`Body::Catalog`] with a structured catalogue id.
//! * User-defined points use [`Body::Derived`], a handle into a
//!   [`DerivedRegistry`](crate::derived::DerivedRegistry) that holds the recipe.
//!
//! The string form round-trips: `Display` writes the two-letter key (`"ju"`), `prefix:id`
//! for catalogue objects (`"hip:32349"`) or `derived:n`, and `FromStr` parses all of them
//! plus full-name aliases, case-insensitively.

use std::fmt;
use std::str::FromStr;

use crate::derived::DerivedId;
use crate::error::BodyParseError;

/// Source catalogue for [`Body::Catalog`] objects.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Catalog {
    /// JPL Horizons / NAIF SPK id, e.g. `20002060` for (2060) Chiron.
    JplHorizons,
    /// Minor Planet Center number, e.g. `90377` for Sedna.
    MinorPlanetCenter,
    /// Gaia DR3 `source_id`.
    Gaia,
    /// Henry Draper catalogue number.
    HenryDraper,
    /// Hipparcos catalogue number.
    Hipparcos,
    /// Yale Bright Star (HR) number.
    BrightStar,
    /// Index into a file of Kepler orbital elements for hypothetical bodies.
    Hypothetical,
}

impl Catalog {
    pub const ALL: [Catalog; 7] = [
        Catalog::JplHorizons,
        Catalog::MinorPlanetCenter,
        Catalog::Gaia,
        Catalog::HenryDraper,
        Catalog::Hipparcos,
        Catalog::BrightStar,
        Catalog::Hypothetical,
    ];

    /// Canonical string prefix used in `prefix:id` keys.
    pub const fn prefix(self) -> &'static str {
        match self {
            Catalog::JplHorizons => "jpl",
            Catalog::MinorPlanetCenter => "mpc",
            Catalog::Gaia => "gaia",
            Catalog::HenryDraper => "hd",
            Catalog::Hipparcos => "hip",
            Catalog::BrightStar => "hr",
            Catalog::Hypothetical => "hyp",
        }
    }

    /// Parses a prefix (case-insensitive), accepting a few common aliases.
    pub fn from_prefix(prefix: &str) -> Option<Self> {
        Some(match prefix.trim().to_ascii_lowercase().as_str() {
            "jpl" | "horizons" | "naif" => Catalog::JplHorizons,
            "mpc" | "mp" => Catalog::MinorPlanetCenter,
            "gaia" | "dr3" => Catalog::Gaia,
            "hd" => Catalog::HenryDraper,
            "hip" | "hipparcos" => Catalog::Hipparcos,
            "hr" | "bsc" => Catalog::BrightStar,
            "hyp" | "fict" => Catalog::Hypothetical,
            _ => return None,
        })
    }

    /// True for star catalogues (positions come from a star catalogue, not an ephemeris file).
    pub const fn is_stellar(self) -> bool {
        matches!(
            self,
            Catalog::Gaia | Catalog::HenryDraper | Catalog::Hipparcos | Catalog::BrightStar
        )
    }
}

impl fmt::Display for Catalog {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.prefix())
    }
}

/// Broad classification, used to dispatch to the right data source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BodyKind {
    /// Sun and Moon (DE ephemeris).
    Luminary,
    /// Planets, Earth and Pluto (DE ephemeris).
    Planet,
    /// Lunar nodes and apsides (calculated, no data file).
    LunarPoint,
    /// Asteroids, centaurs and other minor bodies (small-body SPK files).
    MinorBody,
    /// Moons of other planets (JPL satellite SPK files).
    Satellite,
    /// Uranian / hypothetical bodies (Kepler orbital elements).
    Hypothetical,
    /// Fixed stars (star catalogue).
    Star,
    /// User-defined point (derived registry).
    Derived,
    /// Chart angles, which depend on the time and the place (the ascendant).
    Angle,
}

/// A celestial body or calculated point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[non_exhaustive]
pub enum Body {
    // Sun, Moon, planets and Pluto
    Sun,
    Moon,
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    // Lunar nodes and apsides
    MeanNode,
    TrueNode,
    MeanSouthNode,
    SouthNode,
    MeanApogee,
    OscuApogee,
    IntpApogee,
    IntpPerigee,
    // Main asteroids and centaurs
    Chiron,
    Pholus,
    Ceres,
    Pallas,
    Juno,
    Vesta,
    // Major moons of other planets
    Io,
    Europa,
    Ganymede,
    Callisto,
    Titan,
    Triton,
    Charon,
    // Uranian / Transpluto
    Cupido,
    Hades,
    Zeus,
    Kronos,
    Apollon,
    Admetos,
    Vulkanus,
    Poseidon,
    Isis,
    /// The ascendant: the degree of the ecliptic rising in the east at the given place.
    /// Not a body, but listed alongside them in a chart. Needs a location.
    Ascendant,
    /// A minor body or star identified by catalogue number.
    Catalog { catalog: Catalog, id: u64 },
    /// A user-defined point; its recipe lives in a `DerivedRegistry`.
    Derived(DerivedId),
}

impl Body {
    /// Every variant that has a two-letter key, in declaration order.
    pub const NAMED: [Body; 42] = [
        Body::Sun,
        Body::Moon,
        Body::Mercury,
        Body::Venus,
        Body::Earth,
        Body::Mars,
        Body::Jupiter,
        Body::Saturn,
        Body::Uranus,
        Body::Neptune,
        Body::Pluto,
        Body::MeanNode,
        Body::TrueNode,
        Body::MeanSouthNode,
        Body::SouthNode,
        Body::MeanApogee,
        Body::OscuApogee,
        Body::IntpApogee,
        Body::IntpPerigee,
        Body::Chiron,
        Body::Pholus,
        Body::Ceres,
        Body::Pallas,
        Body::Juno,
        Body::Vesta,
        Body::Io,
        Body::Europa,
        Body::Ganymede,
        Body::Callisto,
        Body::Titan,
        Body::Triton,
        Body::Charon,
        Body::Cupido,
        Body::Hades,
        Body::Zeus,
        Body::Kronos,
        Body::Apollon,
        Body::Admetos,
        Body::Vulkanus,
        Body::Poseidon,
        Body::Isis,
        Body::Ascendant,
    ];

    /// The two luminaries: the Sun and the Moon. Neither is a planet.
    pub const LUMINARIES: [Body; 2] = [Body::Sun, Body::Moon];

    /// The planets as seen from the Earth, Mercury to Pluto. Pluto is listed with them by
    /// astrological convention, although astronomers class it as a dwarf planet. The Earth
    /// is not included; see [`Body::set_for`] for views from elsewhere.
    pub const PLANETS: [Body; 8] = [
        Body::Mercury,
        Body::Venus,
        Body::Mars,
        Body::Jupiter,
        Body::Saturn,
        Body::Uranus,
        Body::Neptune,
        Body::Pluto,
    ];

    /// The luminaries followed by the planets: Sun, Moon and Mercury to Pluto.
    pub const LUMINARIES_AND_PLANETS: [Body; 10] = [
        Body::Sun,
        Body::Moon,
        Body::Mercury,
        Body::Venus,
        Body::Mars,
        Body::Jupiter,
        Body::Saturn,
        Body::Uranus,
        Body::Neptune,
        Body::Pluto,
    ];

    /// The default set of bodies, known by the key `default`: ascendant, Sun, Moon, Mercury to
    /// Pluto, Ketu and Rahu (the true lunar nodes). See [`Body::set_for`] for views from elsewhere than the Earth.
    pub const DEFAULT_BODIES: [Body; 13] = [
        Body::Ascendant,
        Body::Sun,
        Body::Moon,
        Body::Mercury,
        Body::Venus,
        Body::Mars,
        Body::Jupiter,
        Body::Saturn,
        Body::Uranus,
        Body::Neptune,
        Body::Pluto,
        Body::SouthNode,
        Body::TrueNode,
    ];

    /// The `basic` set: the default bodies without Uranus, Neptune and Pluto, as used in Indian astrology.
    pub const BASIC: [Body; 10] = [
        Body::Ascendant,
        Body::Sun,
        Body::Moon,
        Body::Mercury,
        Body::Venus,
        Body::Mars,
        Body::Jupiter,
        Body::Saturn,
        Body::SouthNode,
        Body::TrueNode,
    ];

    /// Adapts a set of bodies to a perspective. Seen from the Earth the set is unchanged.
    /// From anywhere else the Earth is added (after Venus), and what cannot be seen from
    /// there is dropped: the body at the centre, the ascendant and the Earth's lunar points.
    pub fn set_for(set: &[Body], perspective: crate::position::Perspective) -> Vec<Body> {
        if perspective.is_earth_based() {
            return set.to_vec();
        }
        let mut list: Vec<Body> = set
            .iter()
            .copied()
            .filter(|b| Some(*b) != perspective.center() && !matches!(b.kind(), BodyKind::Angle | BodyKind::LunarPoint))
            .collect();
        if !list.contains(&Body::Earth) && perspective.center() != Some(Body::Earth) {
            let at = list.iter().position(|b| !matches!(b, Body::Sun | Body::Moon | Body::Mercury | Body::Venus)).unwrap_or(list.len());
            list.insert(at, Body::Earth);
        }
        list
    }

    /// Shorthand for a catalogue object.
    pub const fn catalog(catalog: Catalog, id: u64) -> Self {
        Body::Catalog { catalog, id }
    }

    /// Two-letter key for named variants; `None` for catalogue and derived objects.
    pub const fn key(&self) -> Option<&'static str> {
        use Body::*;
        Some(match self {
            Sun => "su",
            Moon => "mo",
            Mercury => "me",
            Venus => "ve",
            Earth => "ea",
            Mars => "ma",
            Jupiter => "ju",
            Saturn => "sa",
            Uranus => "ur",
            Neptune => "ne",
            Pluto => "pl",
            MeanNode => "mn",
            TrueNode => "ra",
            MeanSouthNode => "ms",
            SouthNode => "ke",
            MeanApogee => "ml",
            OscuApogee => "ol",
            IntpApogee => "ia",
            IntpPerigee => "ip",
            Chiron => "ch",
            Pholus => "ph",
            Ceres => "ce",
            Pallas => "pa",
            Juno => "jn",
            Vesta => "vs",
            Io => "io",
            Europa => "eu",
            Ganymede => "ga",
            Callisto => "cl",
            Titan => "ti",
            Triton => "tr",
            Charon => "cn",
            Cupido => "cu",
            Hades => "ha",
            Zeus => "ze",
            Kronos => "kr",
            Apollon => "ap",
            Admetos => "ad",
            Vulkanus => "vu",
            Poseidon => "po",
            Isis => "is",
            Ascendant => "as",
            Catalog { .. } | Derived(_) => return None,
        })
    }

    /// English display name for named variants.
    pub const fn name(&self) -> Option<&'static str> {
        use Body::*;
        Some(match self {
            Sun => "Sun",
            Moon => "Moon",
            Mercury => "Mercury",
            Venus => "Venus",
            Earth => "Earth",
            Mars => "Mars",
            Jupiter => "Jupiter",
            Saturn => "Saturn",
            Uranus => "Uranus",
            Neptune => "Neptune",
            Pluto => "Pluto",
            MeanNode => "Mean Node",
            TrueNode => "True Node",
            MeanSouthNode => "Mean South Node",
            SouthNode => "South Node",
            MeanApogee => "Mean Apogee",
            OscuApogee => "Osculating Apogee",
            IntpApogee => "Interpolated Apogee",
            IntpPerigee => "Interpolated Perigee",
            Chiron => "Chiron",
            Pholus => "Pholus",
            Ceres => "Ceres",
            Pallas => "Pallas",
            Juno => "Juno",
            Vesta => "Vesta",
            Io => "Io",
            Europa => "Europa",
            Ganymede => "Ganymede",
            Callisto => "Callisto",
            Titan => "Titan",
            Triton => "Triton",
            Charon => "Charon",
            Cupido => "Cupido",
            Hades => "Hades",
            Zeus => "Zeus",
            Kronos => "Kronos",
            Apollon => "Apollon",
            Admetos => "Admetos",
            Vulkanus => "Vulkanus",
            Poseidon => "Poseidon",
            Isis => "Isis",
            Ascendant => "Ascendant",
            Catalog { .. } | Derived(_) => return None,
        })
    }

    /// Classification used to choose a data source.
    pub const fn kind(&self) -> BodyKind {
        use Body::*;
        match self {
            Sun | Moon => BodyKind::Luminary,
            Mercury | Venus | Earth | Mars | Jupiter | Saturn | Uranus | Neptune | Pluto => {
                BodyKind::Planet
            }
            MeanNode | TrueNode | MeanSouthNode | SouthNode | MeanApogee | OscuApogee
            | IntpApogee | IntpPerigee => BodyKind::LunarPoint,
            Ascendant => BodyKind::Angle,
            Chiron | Pholus | Ceres | Pallas | Juno | Vesta => BodyKind::MinorBody,
            Io | Europa | Ganymede | Callisto | Titan | Triton | Charon => BodyKind::Satellite,
            Cupido | Hades | Zeus | Kronos | Apollon | Admetos | Vulkanus | Poseidon | Isis => {
                BodyKind::Hypothetical
            }
            Catalog { catalog, .. } => match catalog {
                self::Catalog::JplHorizons | self::Catalog::MinorPlanetCenter => {
                    BodyKind::MinorBody
                }
                self::Catalog::Hypothetical => BodyKind::Hypothetical,
                _ => BodyKind::Star,
            },
            Derived(_) => BodyKind::Derived,
        }
    }

    /// The point diametrically opposite, for the node pairs; `None` otherwise.
    pub const fn opposite_node(&self) -> Option<Body> {
        match self {
            Body::MeanNode => Some(Body::MeanSouthNode),
            Body::MeanSouthNode => Some(Body::MeanNode),
            Body::TrueNode => Some(Body::SouthNode),
            Body::SouthNode => Some(Body::TrueNode),
            _ => None,
        }
    }
}

impl fmt::Display for Body {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Body::Catalog { catalog, id } => write!(f, "{}:{}", catalog.prefix(), id),
            Body::Derived(id) => write!(f, "derived:{}", id.0),
            named => f.write_str(named.key().expect("named variants have keys")),
        }
    }
}

impl FromStr for Body {
    type Err = BodyParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        if s.is_empty() {
            return Err(BodyParseError::Empty);
        }
        if let Some((prefix, id)) = s.split_once(':') {
            let id_num: u64 = id
                .trim()
                .parse()
                .map_err(|_| BodyParseError::InvalidId(s.to_string()))?;
            if prefix.trim().eq_ignore_ascii_case("derived") {
                let n = u32::try_from(id_num)
                    .map_err(|_| BodyParseError::InvalidId(s.to_string()))?;
                return Ok(Body::Derived(DerivedId(n)));
            }
            let catalog = self::Catalog::from_prefix(prefix)
                .ok_or_else(|| BodyParseError::UnknownCatalog(prefix.to_string()))?;
            return Ok(Body::Catalog { catalog, id: id_num });
        }
        use Body::*;
        Ok(match s.to_ascii_lowercase().as_str() {
            "su" | "sun" => Sun,
            "mo" | "moon" => Moon,
            "me" | "mercury" => Mercury,
            "ve" | "venus" => Venus,
            "ea" | "earth" => Earth,
            "as" | "asc" | "ascendant" => Ascendant,
            "ma" | "mars" => Mars,
            "ju" | "jupiter" => Jupiter,
            "sa" | "saturn" => Saturn,
            "ur" | "uranus" => Uranus,
            "ne" | "neptune" => Neptune,
            "pl" | "pluto" => Pluto,
            "mn" | "meannode" | "mean_node" => MeanNode,
            "ra" | "rahu" | "truenode" | "true_node" => TrueNode,
            "ms" | "meansouthnode" | "mean_south_node" => MeanSouthNode,
            "ke" | "ketu" | "sn" | "southnode" | "south_node" => SouthNode,
            "ml" | "lilith" | "meanapogee" | "mean_apogee" | "meanapog" => MeanApogee,
            "ol" | "osculatingapogee" | "oscuapogee" | "oscu_apogee" | "oscuapog" => OscuApogee,
            "ia" | "interpolatedapogee" | "intpapogee" | "intp_apogee" | "intpapog" => IntpApogee,
            "ip" | "interpolatedperigee" | "intpperigee" | "intp_perigee" | "intpperg" => IntpPerigee,
            "ch" | "chiron" => Chiron,
            "ph" | "pholus" => Pholus,
            "ce" | "ceres" => Ceres,
            "pa" | "pallas" => Pallas,
            "jn" | "juno" => Juno,
            "vs" | "vesta" => Vesta,
            "io" => Io,
            "eu" | "europa" => Europa,
            "ga" | "ganymede" => Ganymede,
            "cl" | "callisto" => Callisto,
            "ti" | "titan" => Titan,
            "tr" | "triton" => Triton,
            "cn" | "charon" => Charon,
            "cu" | "cupido" => Cupido,
            "ha" | "hades" => Hades,
            "ze" | "zeus" => Zeus,
            "kr" | "kronos" => Kronos,
            "ap" | "apollon" => Apollon,
            "ad" | "admetos" => Admetos,
            "vu" | "vulkanus" => Vulkanus,
            "po" | "poseidon" => Poseidon,
            "is" | "isis" | "transpluto" => Isis,
            _ => return Err(BodyParseError::Unknown(s.to_string())),
        })
    }
}

/// Parses a comma-separated list such as `"su,mo,ma"`. The set names `default`, `basic` and
/// `all` may stand alone or among the keys (`"basic,ch"`); a body named twice is kept once.
pub fn parse_body_list(s: &str) -> Result<Vec<Body>, BodyParseError> {
    parse_body_list_for(s, crate::position::Perspective::Geocentric)
}

/// As [`parse_body_list`], with `default` and `basic` adapted to the perspective by
/// [`Body::set_for`]: away from the Earth they gain the Earth and lose what cannot be seen.
pub fn parse_body_list_for(s: &str, perspective: crate::position::Perspective) -> Result<Vec<Body>, BodyParseError> {
    let mut list: Vec<Body> = Vec::new();
    let mut add = |bodies: Vec<Body>| {
        for b in bodies {
            if !list.contains(&b) {
                list.push(b);
            }
        }
    };
    for token in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
        match token.to_ascii_lowercase().as_str() {
            // `core` is the former name of the default set, still accepted for now.
            "default" | "core" => add(Body::set_for(&Body::DEFAULT_BODIES, perspective)),
            "basic" => add(Body::set_for(&Body::BASIC, perspective)),
            "all" => add(Body::NAMED.to_vec()),
            _ => add(vec![token.parse()?]),
        }
    }
    Ok(list)
}

#[cfg(feature = "serde")]
mod serde_impl {
    use super::Body;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    impl Serialize for Body {
        fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
            s.collect_str(self)
        }
    }

    impl<'de> Deserialize<'de> for Body {
        fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
            let s = <std::borrow::Cow<'de, str>>::deserialize(d)?;
            s.parse().map_err(serde::de::Error::custom)
        }
    }
}
