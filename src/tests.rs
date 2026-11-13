use std::collections::HashSet;

use crate::*;

fn assert_send_sync_copy<T: Send + Sync + Copy + 'static>() {}

#[test]
fn body_is_small_copy_and_thread_safe() {
    assert_send_sync_copy::<Body>();
    assert_send_sync_copy::<JdUt>();
    assert_send_sync_copy::<GeoPos>();
    assert!(std::mem::size_of::<Body>() <= 16);
}

#[test]
fn registry_is_thread_safe() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<DerivedRegistry>();
}

#[test]
fn named_keys_round_trip_and_are_unique() {
    let mut seen = HashSet::new();
    for b in Body::NAMED {
        let key = b.key().expect("named body has a key");
        assert_eq!(key.len(), 2, "{b:?} key {key} is not two letters");
        assert!(seen.insert(key), "duplicate key {key}");
        assert_eq!(key.parse::<Body>().unwrap(), b);
        assert_eq!(b.to_string(), key);
        let name = b.name().unwrap().replace(' ', "");
        assert_eq!(name.parse::<Body>().unwrap(), b, "full name {name}");
    }
}

#[test]
fn parsing_is_case_insensitive_and_trims() {
    assert_eq!(" JU ".parse::<Body>().unwrap(), Body::Jupiter);
    assert_eq!("Jupiter".parse::<Body>().unwrap(), Body::Jupiter);
    assert_eq!("LILITH".parse::<Body>().unwrap(), Body::MeanApogee);
}

#[test]
fn node_aliases_match_api_semantics() {
    assert_eq!("ra".parse::<Body>().unwrap(), Body::TrueNode);
    assert_eq!("ke".parse::<Body>().unwrap(), Body::SouthNode);
    assert_eq!("sn".parse::<Body>().unwrap(), Body::SouthNode);
    assert_eq!("mn".parse::<Body>().unwrap(), Body::MeanNode);
    assert_eq!(Body::TrueNode.opposite_node(), Some(Body::SouthNode));
    assert_eq!(Body::MeanSouthNode.opposite_node(), Some(Body::MeanNode));
}

#[test]
fn unknown_keys_are_errors_not_earth() {
    assert_eq!(
        "xx".parse::<Body>(),
        Err(BodyParseError::Unknown("xx".into()))
    );
    assert_eq!("".parse::<Body>(), Err(BodyParseError::Empty));
}

#[test]
fn catalog_objects_round_trip() {
    for (s, cat, id) in [
        ("hip:32349", Catalog::Hipparcos, 32349),
        ("mpc:90377", Catalog::MinorPlanetCenter, 90377),
        ("jpl:20002060", Catalog::JplHorizons, 20002060),
        ("gaia:4295806720", Catalog::Gaia, 4295806720),
        ("hd:48915", Catalog::HenryDraper, 48915),
        ("hr:2491", Catalog::BrightStar, 2491),
        ("hyp:3", Catalog::Hypothetical, 3),
    ] {
        let b: Body = s.parse().unwrap();
        assert_eq!(b, Body::catalog(cat, id));
        assert_eq!(b.to_string(), s);
        assert_eq!(Catalog::from_prefix(cat.prefix()), Some(cat));
    }
    assert_eq!(
        "HIP:32349".parse::<Body>().unwrap(),
        Body::catalog(Catalog::Hipparcos, 32349)
    );
    assert!(matches!(
        "zz:1".parse::<Body>(),
        Err(BodyParseError::UnknownCatalog(_))
    ));
    assert!(matches!(
        "hip:abc".parse::<Body>(),
        Err(BodyParseError::InvalidId(_))
    ));
}

#[test]
fn kinds_route_to_data_sources() {
    assert_eq!(Body::Sun.kind(), BodyKind::Luminary);
    assert_eq!(Body::Pluto.kind(), BodyKind::Planet);
    assert_eq!(Body::OscuApogee.kind(), BodyKind::LunarPoint);
    assert_eq!(Body::Vesta.kind(), BodyKind::MinorBody);
    assert_eq!(Body::Titan.kind(), BodyKind::Satellite);
    assert_eq!(Body::Isis.kind(), BodyKind::Hypothetical);
    assert_eq!(
        Body::catalog(Catalog::Hipparcos, 1).kind(),
        BodyKind::Star
    );
    assert_eq!(
        Body::catalog(Catalog::MinorPlanetCenter, 1).kind(),
        BodyKind::MinorBody
    );
    assert_eq!(
        Body::catalog(Catalog::Hypothetical, 1).kind(),
        BodyKind::Hypothetical
    );
}

#[test]
fn body_lists() {
    assert_eq!(
        parse_body_list("su, mo,ma,").unwrap(),
        vec![Body::Sun, Body::Moon, Body::Mars]
    );
    assert_eq!(parse_body_list("default").unwrap(), Body::DEFAULT_BODIES.to_vec());
    assert!(parse_body_list("su,xx").is_err());
}

#[test]
fn derived_registry_names_and_references() {
    let mut reg = DerivedRegistry::new();
    let fortuna = reg
        .register(
            "Fortuna",
            Derived::Arabic {
                base: Ref::Angle(Angle::Ascendant),
                plus: Ref::Body(Body::Moon),
                minus: Ref::Body(Body::Sun),
                reverse_at_night: true,
            },
        )
        .unwrap();
    assert_eq!(reg.resolve("fortuna"), Some(fortuna));
    assert_eq!(reg.resolve("ju"), Some(Body::Jupiter));
    assert_eq!(fortuna.to_string(), "derived:0");
    assert_eq!(reg.resolve("derived:0"), Some(fortuna));
    assert_eq!(reg.resolve("derived:9"), None);

    // later recipes may build on earlier ones
    let mid = reg
        .register(
            "fortuna-sun",
            Derived::Midpoint {
                a: Ref::Body(fortuna),
                b: Ref::Body(Body::Sun),
            },
        )
        .unwrap();
    assert_eq!(mid, Body::Derived(DerivedId(1)));

    // forward references are rejected, so cycles are impossible
    assert_eq!(
        reg.register("bad", Derived::Of(DerivedId(7))),
        Err(DerivedError::UnknownReference(7))
    );
    assert_eq!(
        reg.register("FORTUNA", Derived::Fixed { lng: 0.0, lat: 0.0 }),
        Err(DerivedError::DuplicateName("FORTUNA".into()))
    );
    assert_eq!(
        reg.register("mars", Derived::Fixed { lng: 0.0, lat: 0.0 }),
        Err(DerivedError::ReservedName("mars".into()))
    );
}

#[test]
fn geo_parses_api_loc_format() {
    let g: GeoPos = "45.336,13.278".parse().unwrap();
    assert_eq!((g.lat, g.lng, g.alt), (45.336, 13.278, 0.0));
    let g: GeoPos = "-33.865, 151.209, 58".parse().unwrap();
    assert_eq!(g.alt, 58.0);
    assert!("91,0".parse::<GeoPos>().is_err());
    assert!("1,2,3,4".parse::<GeoPos>().is_err());
    assert_eq!(GeoPos::new(0.0, -3.0, 0.0).unwrap().solar_offset_secs(), -720);
}

#[test]
fn julian_day_conversions_match_api_examples() {
    // From the API docs: 1996-01-12T11:00:00 UTC -> jd 2450094.9583333335, unix 821444400
    let jd = JdUt::from_unix(821_444_400);
    assert!((jd.value() - 2_450_094.958_333_333_5).abs() < 1e-9);
    let dt = chrono::NaiveDateTime::parse_from_str("1996-01-12T11:00:00", "%Y-%m-%dT%H:%M:%S")
        .unwrap();
    assert!((JdUt::from(dt).value() - jd.value()).abs() < 1e-9);
    assert_eq!(jd.to_naive().unwrap(), dt);
    // millisecond precision survives the round trip
    let ms = 821_444_400_123;
    assert_eq!(JdUt::from_unix_millis(ms).to_unix_millis(), ms);
}

#[cfg(feature = "serde")]
#[test]
fn serde_uses_string_keys() {
    let v = vec![Body::Jupiter, Body::catalog(Catalog::Hipparcos, 32349)];
    let json = serde_json::to_string(&v).unwrap();
    assert_eq!(json, r#"["ju","hip:32349"]"#);
    let back: Vec<Body> = serde_json::from_str(&json).unwrap();
    assert_eq!(back, v);
}

#[test]
fn iso_strings_are_completed_and_converted_once() {
    let jd = JdUt::from_iso("1990-08-23T12").unwrap();
    assert_eq!(jd.to_iso().unwrap(), "1990-08-23T12:00:00");
    assert_eq!(JdUt::from_iso("1996-01-12T11:25:00").unwrap().value(), 2450094.9756944445);
    assert_eq!("1996-01-12T11:25".parse::<JdUt>().unwrap().value(), 2450094.9756944445);
    assert_eq!(JdUt::from_iso("2024-01").unwrap().to_iso().unwrap(), "2024-01-01T00:00:00");
    assert_eq!(JdUt::from_iso("2024").unwrap().to_iso().unwrap(), "2024-01-01T00:00:00");
    assert_eq!(JdUt::from_iso("2450094.5").unwrap().value(), 2450094.5);
    assert!(matches!(JdUt::from_iso("not a date"), Err(TimeError::Parse(_))));
}
