//! User-defined (derived) points.
//!
//! A derived point is a recipe for computing an ecliptic position from other bodies, angles
//! or fixed times, such as the Part of Fortune, a midpoint or a natal position frozen in time.
//! [`Body::Derived`] holds only a [`DerivedId`], so `Body` stays `Copy + Eq + Hash`; the
//! recipe lives in a [`DerivedRegistry`], normally shared as an `Arc` by all threads.
//!
//! A definition may only reference points already registered, so the registry is acyclic by
//! construction and evaluation always terminates.

use std::collections::HashMap;
use std::sync::Arc;

use crate::body::Body;
use crate::error::DerivedError;
use crate::geo::GeoPos;
use crate::time::JdUt;

/// Handle to a point in a [`DerivedRegistry`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize), serde(transparent))]
pub struct DerivedId(pub u32);

/// Chart angles usable as inputs to derived points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Angle {
    Ascendant,
    Descendant,
    Mc,
    Ic,
    Vertex,
}

/// An input to a derived point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Ref {
    /// A body at the evaluation time.
    Body(Body),
    /// A chart angle at the evaluation time and place.
    Angle(Angle),
    /// A body at a fixed time, optionally topocentric from a fixed place.
    At(Body, JdUt, Option<GeoPos>),
}

/// Recipe for a derived point. All results are ecliptic longitudes (and latitudes where
/// meaningful); speeds are obtained by numeric differentiation.
#[derive(Debug, Clone, PartialEq)]
pub enum Derived {
    /// A fixed ecliptic position, e.g. a natal degree for transposed rise times.
    Fixed { lng: f64, lat: f64 },
    /// Midpoint of two inputs along the shorter arc.
    Midpoint { a: Ref, b: Ref },
    /// Arabic part: `base + plus − minus`, with `plus` and `minus` swapped at night if asked
    /// (e.g. Fortuna = Asc + Moon − Sun).
    Arabic {
        base: Ref,
        plus: Ref,
        minus: Ref,
        reverse_at_night: bool,
    },
    /// Another input shifted by a fixed number of degrees.
    Offset { of: Ref, degrees: f64 },
    /// Position of `body` at the midpoint in time and place of two charts (Davison).
    Davison {
        body: Body,
        a: (JdUt, GeoPos),
        b: (JdUt, GeoPos),
    },
    /// Secondary-style progression: `days_per_year` days after `natal` per year elapsed.
    Progressed {
        body: Body,
        natal: JdUt,
        days_per_year: f64,
    },
    /// Another derived point, for building recipes from recipes.
    Of(DerivedId),
}

impl Derived {
    /// Derived ids this recipe depends on directly.
    fn dependencies(&self) -> impl Iterator<Item = DerivedId> + '_ {
        let refs: [Option<&Ref>; 3] = match self {
            Derived::Midpoint { a, b } => [Some(a), Some(b), None],
            Derived::Arabic {
                base, plus, minus, ..
            } => [Some(base), Some(plus), Some(minus)],
            Derived::Offset { of, .. } => [Some(of), None, None],
            _ => [None, None, None],
        };
        let own = match self {
            Derived::Of(id) => Some(*id),
            Derived::Davison { body, .. } | Derived::Progressed { body, .. } => derived_id(body),
            _ => None,
        };
        refs.into_iter()
            .flatten()
            .filter_map(|r| match r {
                Ref::Body(b) | Ref::At(b, ..) => derived_id(b),
                Ref::Angle(_) => None,
            })
            .chain(own)
    }
}

fn derived_id(body: &Body) -> Option<DerivedId> {
    match body {
        Body::Derived(id) => Some(*id),
        _ => None,
    }
}

/// Named store of derived-point recipes.
#[derive(Debug, Clone, Default)]
pub struct DerivedRegistry {
    defs: Vec<(Arc<str>, Derived)>,
    by_name: HashMap<Arc<str>, DerivedId>,
}

impl DerivedRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Registers a recipe under a unique, case-insensitive name and returns its body handle.
    pub fn register(&mut self, name: &str, def: Derived) -> Result<Body, DerivedError> {
        let key: Arc<str> = Arc::from(name.trim().to_ascii_lowercase());
        if key.parse::<Body>().is_ok() {
            return Err(DerivedError::ReservedName(name.to_string()));
        }
        if self.by_name.contains_key(&key) {
            return Err(DerivedError::DuplicateName(name.to_string()));
        }
        if let Some(missing) = def.dependencies().find(|d| self.get(*d).is_none()) {
            return Err(DerivedError::UnknownReference(missing.0));
        }
        let id = DerivedId(self.defs.len() as u32);
        self.defs.push((key.clone(), def));
        self.by_name.insert(key, id);
        Ok(Body::Derived(id))
    }

    pub fn get(&self, id: DerivedId) -> Option<&Derived> {
        self.defs.get(id.0 as usize).map(|(_, d)| d)
    }

    pub fn name(&self, id: DerivedId) -> Option<&str> {
        self.defs.get(id.0 as usize).map(|(n, _)| &**n)
    }

    pub fn len(&self) -> usize {
        self.defs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.defs.is_empty()
    }

    /// Resolves a registered name first, then any built-in key or `prefix:id`.
    pub fn resolve(&self, key: &str) -> Option<Body> {
        let k = key.trim().to_ascii_lowercase();
        if let Some(id) = self.by_name.get(k.as_str()) {
            return Some(Body::Derived(*id));
        }
        match k.parse::<Body>() {
            Ok(Body::Derived(id)) => self.get(id).map(|_| Body::Derived(id)),
            Ok(b) => Some(b),
            Err(_) => None,
        }
    }
}
