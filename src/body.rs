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
    /// A minor body or star identified by catalogue number.
    Catalog { catalog: Catalog, id: u64 },
    /// A user-defined point; its recipe lives in a `DerivedRegistry`.
    Derived(DerivedId),
}

impl Body {
    /// Every variant that has a two-letter key, in declaration order.
    pub const NAMED: [Body; 41] = [
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
    ];

    /// Sun, Moon, Mercury–Pluto and the true node: the API's `core` set.
    pub const CORE: [Body; 11] = [
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
        Body::TrueNode,
    ];

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

/// Parses a comma-separated list such as `"su,mo,ma"`, or the specials `"core"` and `"all"`.
pub fn parse_body_list(s: &str) -> Result<Vec<Body>, BodyParseError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "core" => Ok(Body::CORE.to_vec()),
        "all" => Ok(Body::NAMED.to_vec()),
        _ => s
            .split(',')
            .filter(|p| !p.trim().is_empty())
            .map(str::parse)
            .collect(),
    }
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
